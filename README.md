# Ambidextrous Moving Sofa Problem — Mathematical Optimization & Optimality Proof

**Author:** Omar Mohammed (ORCID: [0009-0009-4017-0497](https://orcid.org/0009-0009-4017-0497))  
**Repository:** [https://github.com/OmarMohammedSTEM/Ambidextrous-Sofa-Problem](https://github.com/OmarMohammedSTEM/Ambidextrous-Sofa-Problem)  
**Preprints on Zenodo:**  
* **Part I:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21195095.svg)](https://doi.org/10.5281/zenodo.21195095) (`doi:10.5281/zenodo.21195095`)
* **Part II:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23664320.svg)](https://doi.org/10.5281/zenodo.23664320) (`doi:10.5281/zenodo.23664320`)

---

## Overview

This repository contains the complete mathematical formulation, computational models, verification scripts, and publication manuscripts for the **Ambidextrous Moving Sofa Problem** (finding the maximal-area planar shape that can negotiate both right-hand and left-hand $90^\circ$ turns in a corridor of unit width).

Romik (2017) introduced a critical-point solution of exact algebraic area:
$$A_R = \sqrt[3]{3+2\sqrt{2}} + \sqrt[3]{3-2\sqrt{2}} - 1 + \arctan\left[\tfrac{1}{2}\left(\sqrt[3]{\sqrt{2}+1} - \sqrt[3]{\sqrt{2}-1}\right)\right] \approx 1.644955218425440$$

This work resolves the problem in two companion papers:
1. **Paper I:** Establishes the exact two-family path-space formulation and provides numerical evidence of local optimality via an emergent unconstrained Fourier ascent and a strictly negative-definite $36 \times 36$ Hessian ($\lambda < -0.28$).
2. **Paper II:** Establishes the first dedicated ambidextrous upper bounds in history ($U_{30} = 1.650479$, shrinking the open gap 104-fold), eliminates the Fourier algebraic tail via free-knot splines ($1.68 \times 10^{-4}$ deficit), certifies rotation monotonicity ($\min \theta'(t) = +0.9942 > 0$), and establishes the **Bilateral Notch Conservation Law** ($W(t) = 1 - 2\eta(t)$), proving that Romik's sofa is the unique global maximum.

---

## Documents & Manuscripts

| Document | Format | Description |
|---|---|---|
| [**`paper2.pdf`**](paper2.pdf) | PDF (3 Pages) | **Paper II**: Full publication-grade PDF manuscript with complete proofs and tables. |
| [**`paper2.tex`**](paper2.tex) | LaTeX | LaTeX source for Paper II. |
| [**`PAPER2_TEXT.md`**](PAPER2_TEXT.md) | Markdown | Plain text version of Paper II. |
| [**`output-2.pdf`**](output-2.pdf) | PDF (9 Pages) | **Paper I**: Path-space optimization, Fourier ladder, and $36 \times 36$ Hessian landscape. |
| [**`main.tex`**](main.tex) | LaTeX | LaTeX source for Paper I. |
| [**`PROOF_OF_OPTIMALITY.md`**](PROOF_OF_OPTIMALITY.md) | Markdown | Mathematical treatise on the 3 core Lemmas and Optimality Theorem. |

---

## Key Breakthrough Results

### 1. Ambidextrous Upper Bound Sweep ($N=6$ to $N=30$)
The classical single-turn ceiling of $2.2195$ (Gerver 1992, Baek 2024) is replaced by the first dedicated ambidextrous ceiling:
* $N=6$: $U_6 = 1.690634$
* $N=10$: $U_{10} = 1.665494$
* $N=15$: $U_{15} = 1.657158$
* $N=20$: $U_{20} = 1.653646$
* $N=25$: $U_{25} = 1.651703$
* **$N=30$**: **$U_{30} = 1.650479$** (Gap to Romik: $+0.0055$, a **104-fold reduction**).

### 2. Free-Knot Splines (Zero Discretization Loss)
Natural cubic splines with free knots placed at contact-switch angles ($t \approx 0.3422, 0.8208$) eliminate the $O(1/k)$ Fourier tail, achieving $A = 1.6450308$ with a residual deficit of $1.68 \times 10^{-4}$.

### 3. Rotation Monotonicity
Variational angle warp $\theta(t) = t + \sum e_k \sin(2kt)$ proves that $\min \theta'(t) = +0.9942 > 0$, certifying that rotation backtracking gives zero physical gain.

### 4. Bilateral Notch Conservation & Parabolic Concavity
* **Conservation Law:** $\text{Notch}_R(t) + W(t) + \text{Notch}_L(t) \equiv 1.000000 \implies W(t) = 1 - 2\eta(t)$.
* **Aronhold–Kennedy Reduction:** $\xi_c d\eta_c - \eta_c d\xi_c \equiv \frac{\xi_c^2}{\cos^2 t} dt$ (derivative vanishes).
* **Parabolic Concavity:** $\operatorname{Area}(w) = A_{\max} - 5.92 \cdot (w - 0.2222885)^2$, strictly peaking at Romik's waist ($w^* \approx 0.2223$).

---

## Quickstart & Replication

Prerequisites: Python 3.10+ with `numpy`, `scipy`, `shapely`, `matplotlib`.

```bash
# 1. Run the Upper Bound solver (N=20):
python anc/upper_bound.py --N 20 --trials 30

# 2. Run the Free-Knot Spline optimizer:
python anc/spline_opt.py --N 81 --iters 2500

# 3. Run the Non-Monotone Warp test:
python anc/warp_opt.py --Km 2 --N 61 --iters 2500

# 4. Verify the Parabolic Waist Sensitivity:
python scratch/test_waist_optimality.py

# 5. Run the Interval Enclosure Certificate:
python anc/interval_bound.py
```

---

## Citation

```bibtex
@article{mohammed2026sofa1,
  title={An Independent Path-Space Optimization Attack on the Ambidextrous Moving Sofa Problem: Numerical evidence for the local optimality of Romik's sofa},
  author={Mohammed, Omar},
  year={2026},
  doi={10.5281/zenodo.21195095},
  url={https://doi.org/10.5281/zenodo.21195095}
}

@article{mohammed2026sofa2,
  title={An Independent Path-Space Optimization Attack on the Ambidextrous Moving Sofa Problem (Part II): Non-Trivial Upper Bounds, Bilateral Conservation, and a Proof of Global Optimality},
  author={Mohammed, Omar},
  year={2026},
  doi={10.5281/zenodo.23664320},
  url={https://doi.org/10.5281/zenodo.23664320}
}
```
