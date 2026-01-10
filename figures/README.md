# Figure Generation

Reproducible scripts for all publication figures.

## Figures

| Script | Figure | Description |
|--------|--------|-------------|
| `fig1_state_equation.py` | Fig 1 | Geometric state equation diagram |
| `fig2_convergence.py` | Fig 2 | 5-seed coordinate convergence |
| `fig3_viral_depth.py` | Fig 3 | κ vs phylogenetic depth |
| `fig4_curvature_entropy.py` | Fig 4 | Curvature-entropy law validation |

## Usage

### Generate All Figures

```bash
make figures
```

### Generate Single Figure

```bash
python figures/fig4_curvature_entropy.py
```

## Output

Generated figures are saved to `outputs/`:
- PNG (300 DPI) for manuscripts
- PDF (vector) for publication

## Dependencies

Figures read from:
1. `../constants.yaml` - Canonical values
2. `../validation/*/results/*.yaml` - Experimental results

This ensures figures always reflect the latest validated results.

## Style

All figures use:
- Okabe-Ito colorblind-safe palette
- Sans-serif fonts (Arial/Helvetica)
- Nature-style minimalist design
