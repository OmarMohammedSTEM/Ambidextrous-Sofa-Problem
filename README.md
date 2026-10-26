# An Independent Path-Space Optimization Attack on the Ambidextrous Moving Sofa Problem

**Author:** Omar Mohammed (ORCID: [0009-0009-4017-0497](https://orcid.org/0009-0009-4017-0497))  
**Zenodo Preprint:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21195095.svg)](https://doi.org/10.5281/zenodo.21195095) (`doi:10.5281/zenodo.21195095`)  
**Repository:** [https://github.com/OmarMohammedSTEM/Ambidextrous-Sofa-Problem](https://github.com/OmarMohammedSTEM/Ambidextrous-Sofa-Problem)

---

## Overview

The **ambidextrous moving sofa problem** asks for the planar shape of maximal area that can navigate both right-hand and left-hand $90^\circ$ turns in a hallway of unit width. In 2016 (arXiv:1606.08111, published in *Experimental Mathematics* 2018), Dan Romik derived an exact stationary solution to the Euler–Lagrange equations with closed-form algebraic area:
$$A_R = \sqrt[3]{3+2\sqrt{2}} + \sqrt[3]{3-2\sqrt{2}} - 1 + \arctan\left[\tfrac{1}{2}\left(\sqrt[3]{\sqrt{2}+1} - \sqrt[3]{\sqrt{2}-1}\right)\right] \approx 1.644955218425440$$
bounded by 18 analytic pieces (circular arcs and sextics).

This repository contains the exact path-space formulation, numerical optimization engines, and Hessian landscape analysis providing rigorous computational evidence for the **local optimality** of Romik's sofa.

---

## Key Results (Part I)

1. **Path-Space Formulation:**  
   The ambidextrous sofa is parametrized over two independent rotation-center paths $\mathbf{v}(t)$ and $\mathbf{w}(s)$, releasing symmetry assumptions to evaluate independent motions for right-hand and left-hand turns.
2. **Symmetric Convergence:**  
   Starting from an unconstrained double-Hammersley baseline ($A = 1.5465$), unconstrained Fourier ascent converges toward Romik's exact area ($A = 1.644955$) and recovers Romik's rotation center paths.
3. **Hessian Spectral Analysis:**  
   On a 36-dimensional Fourier section:
   - 1 zero eigenvalue corresponds to the exact horizontal translation gauge symmetry.
   - 33 eigenvalues are strictly negative ($\lambda < -0.28$).
   - 2 soft modes exist at $\lambda \approx -0.011$ to $-0.014$ (line scans confirm negative curvature $\approx -0.02$).
   - The spectrum is consistent with a strict local maximum on this section.

---

## Manuscript

The full 9-page research paper is available in the repository:
* [**`paper_preprint.pdf`**](paper_preprint.pdf) — *"An Independent Path-Space Optimization Attack on the Ambidextrous Moving Sofa Problem: Numerical evidence for the local optimality of Romik's sofa"*

---

## Code & Replication

All core numerical scripts are located in `src/`:

```bash
# 1. Evaluate area of Romik configuration:
python src/ambi.py

# 2. Run Newton polish and boundary verification:
python src/polish.py

# 3. Compute the 36x36 Hessian and spectrum:
python src/hess.py
python src/hess2.py

# 4. Perform direct 1D linescans along soft modes:
python src/linescan.py
```

---

## Citation

```bibtex
@article{mohammed2026sofa,
  title={An Independent Path-Space Optimization Attack on the Ambidextrous Moving Sofa Problem: Numerical evidence for the local optimality of Romik's sofa},
  author={Mohammed, Omar},
  year={2026},
  doi={10.5281/zenodo.21195095},
  url={https://doi.org/10.5281/zenodo.21195095}
}
```
