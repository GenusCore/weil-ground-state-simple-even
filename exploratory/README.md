# Exploratory scripts (NON-RIGOROUS)

Nothing in this directory is part of the certificate. These scripts use floating point or high-precision
(mpmath) Galerkin matrices from `code/convergence.py`, not verified bounds. They back the observations in
Appendix A of `paper/NOTE.md` and the normalization check in §2.

| Script | Purpose |
|---|---|
| `weil_z.py` | Weil's form for zeta from the prime side, checked against the zero side (first 400 zeros, cached in `zeros.txt`; computed with `mpmath.zetazero`, double precision). Gives Rayleigh-Ritz upper bounds only. Numpy only. |
| `zeros.txt` | Imaginary parts of the first 400 nontrivial zeta zeros, from `mpmath.zetazero`, rounded to double precision. |
| `ladder.py` | A.1: derivative ladder overlaps. |
| `osc.py`, `osc2.py` | A.2(1)-(2): sign-change counts of eigenfunctions; overlap with Riemann's Phi. |
| `tp/phi.py` | A.2(4): rank-one intertwining / boundary-trace residues. |
| `tp/osc.py` | A.2(1): independent sign-change count (internal review). |
| `tp/fold.py` | A.2(3): folding test. |
| `ref_zop/chk.py` | A.4: odd-sector orthogonality check (internal review). |

Run from this directory with the same Python environment, for example `python ladder.py`. Several take
command-line arguments (`L`, degree `D`); see the first lines of each file. Some are slow at high degree.
