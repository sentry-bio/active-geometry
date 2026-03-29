#!/usr/bin/env python3
"""Fine κ scan [0.70, 1.32] over valid pre-cliff range for fungal telescope."""

import json, numpy as np, torch, sys, math, os, random, re
import pandas as pd
from scipy.stats import spearmanr, pearsonr
import dendropy

sys.path.insert(0, '/home/rohit/biosphere_inference')
sys.path.insert(0, '/home/rohit')
from model_v15_5 import V15Model


def poincare_dist_mat(Z, c, eps=1e-7):
    B = Z.shape[0]
    u = Z.unsqueeze(1).expand(B, B, -1)
    v = Z.unsqueeze(0).expand(B, B, -1)
    diff_sq = ((u - v) ** 2).sum(-1)
    u_sq = (u ** 2).sum(-1)
    v_sq = (v ** 2).sum(-1)
    denom = ((1 - c * u_sq) * (1 - c * v_sq)).clamp(min=eps)
    arg = (1 + 2 * c * diff_sq / denom).clamp(min=1 + eps)
    return torch.acosh(arg) / torch.sqrt(c + eps)


print('Loading checkpoint...')
ckpt = torch.load('/fast/sentrybio/checkpoints/v15_5_compact2_diagnostic/best.pt',
                   map_location='cpu', weights_only=False)
state = ckpt.get('model_state_dict', ckpt)
latent_dim = state['ode_flow.field.0.weight_orig'].shape[1] - 1
ode_hidden  = state['ode_flow.field.0.bias'].shape[0]
counts = {k: state[f'{k}.prototypes'].shape[0]
          for k in ['bact_fam','arch_fam','euk_fam','bact_gen','arch_gen','euk_gen']}
vocab_size = state['encoder.encoder.embed.weight'].shape[0]
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = V15Model(vocab_size=vocab_size, latent_dim=latent_dim,
                 counts=counts, ode_hidden=ode_hidden).to(device)
state_f = {k: v for k, v in state.items() if 'curvature_history' not in k}
model.load_state_dict(state_f, strict=False)
model.eval()
print(f'live_kappa={model.live_kappa:.6f}')

df = pd.read_csv('/home/rohit/eggnog/manifest_go_kegg_v2.csv',
                  usecols=['accession','domain','species','phylum','tokenized_path'],
                  low_memory=False)
fungi = df[
    (df['domain'] == 'Eukaryota') &
    (df['phylum'].isin(['Ascomycota','Basidiomycota','Mucoromycota','Microsporidia']))
].copy()
fungi = fungi[fungi['tokenized_path'].apply(
    lambda p: bool(p) and isinstance(p, str) and os.path.exists(p)
)].reset_index(drop=True)

print('Loading tree...')
tree = dendropy.Tree.get(path='/home/rohit/li2021_fungal_tree.treefile',
                          schema='newick', preserve_underscores=True)
taxon_labels = {t.label: t for t in tree.taxon_namespace}
taxon_norm = {lbl.replace('_', ' ').lower(): lbl for lbl in taxon_labels}

fungi['species_norm'] = fungi['species'].fillna('').str.lower()
matched_df = fungi[fungi['species_norm'].isin(taxon_norm)].copy().reset_index(drop=True)
n = len(matched_df)
print(f'Matched: {n} genomes')
print('Phyla:', matched_df['phylum'].value_counts().to_dict())

print('Building patristic matrix...')
pdm = tree.phylogenetic_distance_matrix()
D = np.zeros((n, n), dtype=np.float32)
for i in range(n):
    lbl_i = taxon_norm[matched_df.loc[i, 'species_norm']]
    tax_i = taxon_labels[lbl_i]
    for j in range(i + 1, n):
        lbl_j = taxon_norm[matched_df.loc[j, 'species_norm']]
        tax_j = taxon_labels[lbl_j]
        d = pdm.distance(tax_i, tax_j)
        D[i, j] = d
        D[j, i] = d
    if (i + 1) % 100 == 0:
        print(f'  PDM row {i+1}/{n}')
D_max = D.max()
D /= D_max
print(f'Patristic done. raw_max={D_max:.4f}')

print('Extracting embeddings...')
random.seed(42)
np.random.seed(42)
embs = []
with torch.no_grad():
    for i in range(n):
        tokens = np.load(matched_df.loc[i, 'tokenized_path']).astype(np.int64)
        if len(tokens) > 8192:
            s = random.randint(0, len(tokens) - 8192)
            tokens = tokens[s:s + 8192]
        t = torch.from_numpy(tokens).unsqueeze(0).to(device)
        z = model.encode_angular_only(t)
        embs.append(z.squeeze(0).cpu())
        if (i + 1) % 200 == 0:
            print(f'  {i+1}/{n} embedded')
Z = torch.stack(embs)

radii = Z.norm(dim=1).numpy()
training_kappa = model.live_kappa
ball_radius = 1 / training_kappa**0.5
print(f'Radii: mean={radii.mean():.4f} max={radii.max():.4f}')
print(f'Ball radius at κ={training_kappa:.4f}: {ball_radius:.4f}')
pct_outside = (radii > ball_radius).mean() * 100
print(f'Fraction outside ball: {pct_outside:.1f}%')

# Fine scan over valid range [0.70, 1.32] step=0.02
kappas = list(np.arange(0.70, 1.33, 0.02))
print(f'\nFine scan κ=[{kappas[0]:.2f}, {kappas[-1]:.2f}] step=0.02:')

Z_dev = Z.to(device)
D_dev = torch.tensor(D, dtype=torch.float32, device=device)
mask = torch.triu(torch.ones(n, n, dtype=torch.bool, device=device), diagonal=1)
d_ref = D_dev[mask].cpu().numpy()

rows = []
for kappa in kappas:
    c = torch.tensor(kappa, dtype=torch.float32, device=device)
    with torch.no_grad():
        d_hyp = poincare_dist_mat(Z_dev, c)[mask].cpu().numpy()
    # Check fraction of pairs with valid denominators
    rho, _ = spearmanr(d_hyp, d_ref)
    pr, _ = pearsonr(np.log(d_hyp + 1e-8), np.log(d_ref + 1e-8))
    rows.append({'kappa': float(kappa), 'spearman': float(rho), 'pearson_log': float(pr)})
    print(f'  κ={kappa:.2f}  ρ={rho:.4f}  pearson_log={pr:.4f}')

best_s = max(rows, key=lambda x: x['spearman'])
best_p = max(rows, key=lambda x: x['pearson_log'])
print(f'\nSpearman peak:    κ={best_s["kappa"]:.2f}, ρ={best_s["spearman"]:.4f}')
print(f'Pearson-log peak: κ={best_p["kappa"]:.2f}, r={best_p["pearson_log"]:.4f}')

result = {
    'n_genomes': n,
    'phyla': matched_df['phylum'].value_counts().to_dict(),
    'training_kappa': float(training_kappa),
    'embedding_radius_mean': float(radii.mean()),
    'embedding_radius_max': float(radii.max()),
    'ball_radius_training_kappa': float(ball_radius),
    'pct_outside_ball': float(pct_outside),
    'patristic_raw_max': float(D_max),
    'fine_valid_scan': rows,
    'spearman_peak_kappa': best_s['kappa'],
    'spearman_peak_rho': best_s['spearman'],
    'pearson_log_peak_kappa': best_p['kappa'],
    'pearson_log_peak': best_p['pearson_log'],
}
with open('/home/rohit/fungal_telescope_fine.json', 'w') as f:
    json.dump(result, f, indent=2)
print('Saved -> /home/rohit/fungal_telescope_fine.json')
