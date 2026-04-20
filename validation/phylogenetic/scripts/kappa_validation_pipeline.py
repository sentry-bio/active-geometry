# === One-Cell, Lightweight κ-Validation Pipeline ===
# What it does:
#   - Installs deps
#   - Uploads phylogenetic trees (Newick/Nexus)
#   - Auto-detects schema per file
#   - Profiles κ around 1.25 for d in {2,3,4}
#   - Plots stress vs κ and prints Δstress at κ=1.25
#   - Optional mini-bootstrap (B=30) for quick CI (set MINI_BOOTSTRAP=True)

# --------- Install deps (Colab-safe) ----------
import sys, subprocess, pkgutil
def _pip_install(pkgs):
    subprocess.check_call([sys.executable, "-m", "pip", "-q", "install", *pkgs])
need = []
for p in ["numpy","scipy","matplotlib","biopython","dendropy","scikit-learn","torch","tqdm"]:
    if pkgutil.find_loader(p) is None:
        need.append(p)
if need:
    _pip_install(need)

# --------- Imports & config ----------
import numpy as np, scipy as sp
import scipy.optimize
import matplotlib.pyplot as plt
from tqdm import tqdm
import torch
from sklearn.manifold import MDS
import dendropy
import argparse
import os

np.random.seed(1337)
torch.manual_seed(1337)

parser = argparse.ArgumentParser(description='Lightweight κ-validation pipeline')
parser.add_argument('--files', nargs='*', help='Paths to tree files (.nwk/.tre/.treorg/.nex)')
parser.add_argument('--no_bootstrap', action='store_true', help='Disable mini-bootstrap')
parser.add_argument('--max_taxa', type=int, default=400, help='Subsample cap for large trees (None to disable)')
parser.add_argument('--dims', type=int, nargs='*', default=[2,3,4], help='Embedding dimensions to test')
parser.add_argument('--kmin', type=float, default=0.8, help='Min κ in grid')
parser.add_argument('--kmax', type=float, default=1.7, help='Max κ in grid')
parser.add_argument('--ksteps', type=int, default=91, help='Number of κ grid points')
parser.add_argument('--no_refine', action='store_true', help='Skip local κ refinement and use grid minimum')
parser.add_argument('--scale', type=str, default='median', choices=['none','median','mean','max'], help='Normalize patristic distances')
args, _ = parser.parse_known_args()

TARGET_K = 1.25
K_GRID   = np.linspace(args.kmin, args.kmax, args.ksteps)
DIMS     = args.dims
MINI_BOOTSTRAP = not args.no_bootstrap
MAX_TAXA = args.max_taxa

print('Versions → numpy', np.__version__, '| torch', torch.__version__)

# --------- Upload trees ----------
if args.files:
    file_paths = args.files
else:
    try:
        from google.colab import files
        print('Upload one or more phylogenetic trees: .nwk/.tre/.treorg (Newick) or .nex/.nexorg (Nexus)')
        uploaded = files.upload()
        file_paths = list(uploaded.keys())
    except Exception as e:
        print("Colab 'files.upload()' not available. Using local file paths.")
        file_paths = [
            "solanaceae_tree_fixed.nwk",
            "bayesian_tree_fixed.nwk"
        ]

# --------- I/O: Read tree → patristic distance matrix ----------
def load_tree_any(path):
    """
    Try Newick first, then Nexus. Works for .nwk/.tre/.treorg and .nex/.nexorg
    """
    tried = []
    
    # If looks like NEXUS but missing leading '#', normalize
    try:
        with open(path, 'r') as f:
            head = f.read(16)
        if head.startswith('NEXUS'):
            with open(path, 'r') as f:
                content = f.read()
            if not content.lstrip().startswith('#NEXUS'):
                content = '#NEXUS\n' + content
            try:
                tree_collection = dendropy.TreeList.get(data=content, schema="nexus", preserve_underscores=True)
                if len(tree_collection) > 0:
                    return tree_collection[0]
            except Exception:
                pass
    except Exception:
        pass

    # First try as pure Newick
    try:
        t = dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True)
        return t
    except Exception as ex:
        tried.append(("newick", str(ex)))
    
    # Try as NEXUS format
    try:
        # For NEXUS files, try to extract trees
        tree_collection = dendropy.TreeList.get(path=path, schema="nexus", preserve_underscores=True)
        if len(tree_collection) > 0:
            return tree_collection[0]  # Return first tree
        else:
            raise RuntimeError("No trees found in NEXUS file")
    except Exception as ex:
        tried.append(("nexus", str(ex)))
    
    # Try reading as NEXUS but extracting tree manually
    try:
        with open(path, 'r') as f:
            content = f.read()
        
        # Look for tree definition
        import re
        tree_match = re.search(r'TREE\s+\w+\s*=\s*(.+?);', content, re.IGNORECASE | re.DOTALL)
        if tree_match:
            newick_str = tree_match.group(1).strip()
            # Clean up any extra annotations
            newick_str = re.sub(r'\[%[^\]]*\]', '', newick_str)
            if not newick_str.endswith(';'):
                newick_str = newick_str + ';'
            
            # Parse as Newick string
            t = dendropy.Tree.get(data=newick_str, schema="newick", preserve_underscores=True)
            return t
        else:
            raise RuntimeError("No TREE definition found in file")
    except Exception as ex:
        tried.append(("manual_nexus", str(ex)))
    
    raise RuntimeError(f"Failed to parse {path}. Tried: {tried}")

def patristic_dm_from_tree(path):
    tree = load_tree_any(path)
    leaves = list(tree.leaf_node_iter())
    n = len(leaves)
    pdm = tree.phylogenetic_distance_matrix()
    D = np.zeros((n,n), dtype=float)
    for i, ti in enumerate(leaves):
        for j, tj in enumerate(leaves):
            # Use Taxon objects (not Node) for distance lookup
            D[i,j] = pdm.patristic_distance(ti.taxon, tj.taxon)
    # Prefer taxon labels if available
    labels = [
        str(t.taxon.label) if (t.taxon is not None and t.taxon.label is not None) else str(t.label)
        for t in leaves
    ]
    return labels, D

trees = []
for p in file_paths:
    try:
        labels, D = patristic_dm_from_tree(p)
        if D.shape[0] >= 8:   # ignore tiny trees
            trees.append({'name': p, 'labels': labels, 'D': D})
            print(f"Loaded: {p} | n={D.shape[0]}")
        else:
            print(f"Skipped small tree: {p} (n={D.shape[0]})")
    except Exception as e:
        print(f"[WARN] Could not parse {p}: {e}")

if not trees:
    raise SystemExit("No valid trees loaded. Please re-run and upload .nwk/.tre/.nex files with branch lengths.")

# --------- Hyperbolic geometry (Poincaré ball, curvature −κ²) ----------
def _arcosh(x):
    x = np.maximum(x, 1.0 + 1e-12)
    return np.log(x + np.sqrt(x*x - 1.0))

def poincare_distance_matrix(X, kappa):
    # X in ball; we keep points near origin, so this is numerically stable
    X2 = np.sum(X*X, axis=1, keepdims=True)             # (n,1)
    UV = np.sum((X[:,None,:]-X[None,:,:])**2, axis=2)   # (n,n) squared Euclid
    denom = (1 - (kappa**2)*X2)
    denom_mat = denom @ denom.T
    z = 1 + 2*(kappa**2)*UV/denom_mat
    z = np.maximum(z, 1.0 + 1e-12)
    return (1.0/kappa) * _arcosh(z)

def stress_loss(X, D_target, kappa):
    D_hat = poincare_distance_matrix(X, kappa)
    return float(np.mean((D_hat - D_target)**2))

def init_from_mds(D, d, kappa, seed=0):
    mds = MDS(n_components=d, dissimilarity='precomputed', random_state=seed, normalized_stress='auto')
    X_euc = mds.fit_transform(D)
    # scale down safely into the ball
    X = X_euc / (10.0 * (np.max(np.linalg.norm(X_euc, axis=1)) + 1e-9))
    return X.astype(np.float64)

# --------- Optimizer (LBFGS) for given (κ,d) ----------
def optimize_embedding(D, d, kappa, n_restarts=2, max_iter=120, seed=123):
    best_loss, best_X = np.inf, None
    for r in range(n_restarts):
        X0 = init_from_mds(D, d, kappa, seed=seed+r)
        X  = torch.tensor(X0, dtype=torch.double, requires_grad=True)
        D_t= torch.tensor(D, dtype=torch.double)
        k  = torch.tensor(float(kappa), dtype=torch.double)
        opt= torch.optim.LBFGS([X], lr=1.0, max_iter=20, history_size=50, line_search_fn='strong_wolfe')

        def arcosh_t(z):
            z = torch.clamp(z, min=1.0 + 1e-12)
            return torch.log(z + torch.sqrt(z*z - 1.0))

        def closure():
            opt.zero_grad()
            X2 = (X*X).sum(dim=1, keepdim=True)
            diff = X.unsqueeze(1) - X.unsqueeze(0)
            UV = (diff*diff).sum(dim=2)
            denom = (1 - (k*k)*X2)
            denom_mat = denom @ denom.T
            z = 1 + 2*(k*k)*UV/denom_mat
            D_hat = (1.0/k) * arcosh_t(z)
            loss = torch.mean((D_hat - D_t)**2)
            if loss.requires_grad and torch.is_grad_enabled():
                loss.backward()
            return loss

        opt.step(closure)
        with torch.no_grad():
            loss_t = closure()
            if not torch.isfinite(loss_t):
                loss = np.inf
            else:
                loss = float(loss_t.item())
            if loss < best_loss:
                best_loss = loss
                best_X = X.detach().numpy()
    return best_X, best_loss

# --------- κ profiling (fast) ----------
def profile_kappa_fast(D, k_grid=K_GRID, dims=DIMS):
    records = []
    for kappa in tqdm(k_grid, desc='κ grid'):
        best_loss, best_d = np.inf, None
        for d in dims:
            _, loss = optimize_embedding(D, d, kappa, n_restarts=2, max_iter=120)
            if np.isfinite(loss) and loss < best_loss:
                best_loss, best_d = loss, d
        if best_d is not None and np.isfinite(best_loss):
            records.append((float(kappa), int(best_d), float(best_loss)))
        else:
            # skip this κ if no valid fit
            continue
    if len(records) == 0:
        raise RuntimeError("No valid κ records; check data or reduce MAX_TAXA/dims/grid.")
    arr = np.array(records, dtype=object)
    i0 = int(np.argmin(arr[:,2].astype(float)))
    k0, d0, L0 = float(arr[i0,0]), int(arr[i0,1]), float(arr[i0,2])
    if getattr(args, 'no_refine', False):
        return records, (k0, d0, L0)
    # local refine in κ with d fixed (fast)
    def f(k):
        _, L = optimize_embedding(D, d0, float(k), n_restarts=2, max_iter=120)
        return L
    a, b = max(K_GRID[0], k0-0.05), min(K_GRID[-1], k0+0.05)
    res = sp.optimize.minimize_scalar(f, bounds=(a,b), method='bounded', options={'xatol':1e-3})
    k_hat, L_hat = float(res.x), float(res.fun)
    return records, (k_hat, d0, L_hat)

# --------- Optional mini-bootstrap for quick CI ----------
def bootstrap_distance_matrix(D, rng, sigma=0.05):
    n = D.shape[0]
    noise = rng.lognormal(mean=0.0, sigma=sigma, size=(n,n))
    noise = 0.5*(noise + noise.T)
    np.fill_diagonal(noise, 1.0)
    return D*noise

def bootstrap_kappa_fast(D, B=30):
    rng = np.random.RandomState(2024)
    ks = []
    for _ in tqdm(range(B), desc='mini-bootstrap'):
        Db = bootstrap_distance_matrix(D, rng)
        _, (k_hat, _, _) = profile_kappa_fast(Db, k_grid=np.linspace(0.9,1.6,49), dims=DIMS)
        ks.append(k_hat)
    ks = np.array(ks)
    lo, hi = np.percentile(ks, [2.5, 97.5])
    se = ks.std(ddof=1)
    return ks, (lo, hi, se)

# --------- Run per-tree analysis ----------
for t in trees:
    name, D = t['name'], t['D']
    print(f"\n=== {name} ===")
    # Subsample very large trees for speed and stability
    do_bootstrap = MINI_BOOTSTRAP
    if MAX_TAXA is not None and D.shape[0] > MAX_TAXA:
        rng = np.random.RandomState(123)
        idx = np.sort(rng.choice(D.shape[0], size=MAX_TAXA, replace=False))
        D = D[np.ix_(idx, idx)]
        t['labels'] = [t['labels'][i] for i in idx]
        print(f"Subsampled tree from n={t['D'].shape[0]} to n={D.shape[0]} for speed.")
        do_bootstrap = False
    # Distance normalization
    if args.scale and args.scale != 'none':
        iu = np.triu_indices(D.shape[0], 1)
        if args.scale == 'median':
            s = np.median(D[iu])
        elif args.scale == 'mean':
            s = np.mean(D[iu])
        elif args.scale == 'max':
            s = np.max(D[iu])
        else:
            s = 1.0
        s = float(s) if s > 0 else 1.0
        D = D / s
    prof, (k_hat, d_hat, L_hat) = profile_kappa_fast(D)
    t['profile'] = prof; t['k_hat'] = k_hat; t['d_hat'] = d_hat; t['L_hat'] = L_hat

    # stress at target κ (use best d for speed)
    _, L_target = optimize_embedding(D, d_hat, TARGET_K, n_restarts=3, max_iter=150)
    delta = L_target - L_hat

    # Plot stress vs κ
    arr = np.array(prof, dtype=object)
    kappas = arr[:,0].astype(float)
    losses = arr[:,2].astype(float)
    plt.figure()
    plt.plot(kappas, losses)
    plt.axvline(TARGET_K, linestyle='--')
    plt.title(f'Stress vs κ — {name}')
    plt.xlabel('κ'); plt.ylabel('Mean squared stress')
    base = os.path.basename(name)
    safe = base.replace('/', '_').replace('\\', '_').replace('.', '_')
    out_png = f"stress_vs_kappa_{safe}.png"
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved plot → {out_png}")

    print(f"κ̂ = {k_hat:.4f} (best d={d_hat}) | Stress@κ̂ = {L_hat:.3e}")
    print(f"Stress@κ=1.25 = {L_target:.3e} | Δstress(1.25 − κ̂) = {delta:.3e}")

    if do_bootstrap:
        print(f"[Bootstrap] {name} — running quick B=30 (≈ few minutes)...")
        ks, (lo, hi, se) = bootstrap_kappa_fast(D, B=30)
        print(f"κ̂ bootstrap: mean={ks.mean():.4f}, SE={se:.4f}, 95% CI=[{lo:.4f}, {hi:.4f}]")

print("\nDone. If most trees have Δstress near 0 and κ̂ within ~±0.02 of 1.25 (with CI overlap), you have a solid lightweight invariance signal.")