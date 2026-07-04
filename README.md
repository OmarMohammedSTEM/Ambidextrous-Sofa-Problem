# Ambidextrous Moving Sofa: Path-Space Optimization

Code and data for *"An Independent Path-Space Optimization Attack on the
Ambidextrous Moving Sofa Problem: Numerical evidence for the local optimality
of Romik's sofa"* (Omar Mohammed, 2026).

The ambidextrous moving sofa problem asks for the planar shape of maximal area
that can negotiate 90° corridor turns of **both** handedness (corridor width 1).
The best known shape is Romik's ambidextrous sofa (2017), area
≈ 1.644955218, whose optimality — even *local* optimality — is an open problem.

This project reformulates the problem exactly as an optimization over two
independent rotation-center paths, attacks the record with a free-form Fourier
parametrization that releases Romik's symmetry assumptions, and finds that the
optimizer reconverges to Romik's sofa from below (best area
1.6449449 ± 3.4×10⁻⁶) with a numerically negative-definite Hessian on the full
36-parameter section — evidence that Romik's sofa is a strict local maximum.

## Requirements

Python 3.11+ with `numpy`, `scipy`, `shapely`, `matplotlib`
(tested: NumPy 2.3, SciPy 1.16, Shapely 2.1).

```
pip install numpy scipy shapely matplotlib
```

## Scripts (run in this order to regenerate everything)

| script | what it does | output |
|---|---|---|
| `src/ambi.py` | geometry kernel (polygon-intersection area evaluator, Richardson extrapolation) + staged optimization: double-Hammersley baseline, symmetric Fourier ladder K=1..3, fully asymmetric K=3 run, directional probes | `data/ambi_results.json` |
| `src/polish.py` | K=4 refinement + shape figure | `data/ambi_results.json` (key `B_K4`), figure |
| `src/k5.py` | K=5 continuation | key `B_K5` |
| `src/newton.py` | curvature-corrected Newton polish → **headline configuration** | key `B_K5_newton` |
| `src/hess.py` | full 36×36 finite-difference Hessian eigendecomposition at h=0.015 | `data/hess.json` |
| `src/hess2.py <h>` | Hessian at another step size (run with `0.0075` and `0.03`); saves full matrix | `data/hess_h<h>.json` |
| `src/linescan.py` | direct area scans along the softest Hessian eigenvectors: machine-zero gauge check + effective-λ fits (requires `data/hess_h0.03.json`) | `data/linescan.json` |

Total runtime ≈ 2.5 h on a desktop, dominated by the three Hessians
(~2600 area evaluations each). Run all scripts from the repository root
(e.g. `python src/ambi.py`). All results are already checkpointed in `data/`,
so you can inspect every number in the paper without rerunning anything.

## Key numbers

| quantity | value |
|---|---|
| Romik's exact area | 1.644955218425440 |
| best reconstruction (K=5 + Newton) | 1.6449449 ± 3.4×10⁻⁶ |
| Hessian spectrum (36-dim section) | all eigenvalues in [−170.5, −0.0107] |
| translation-gauge line scan | machine zero (exactness check) |
| soft-mode effective λ (direct scans) | ≈ −0.018, −0.022 (strictly negative) |

## Validation

The evaluator is calibrated on the single-turn Hammersley sofa (exact area
π/2 + 2/π): discretization bias is 0.64/N to within 2% over four doublings of
N, and one Richardson step 2·A(2N) − A(N) leaves an O(N⁻²) residual;
successive extrapolations agree to 3.4×10⁻⁶, which is the error bar used
throughout.

## Files

- `paper_preprint.pdf` — the paper
- `figures/` — shape and convergence figures
- `data/` — checkpointed results: coefficients, Hessians, scans (JSON)

## License

MIT
