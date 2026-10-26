# Simple, even ground state of the localized Weil form for L in (0, 1.19]

This is a computer-assisted proof, using Arb/FLINT ball arithmetic, of the following statement. For every window
half-width `L` in `(0, 1.19]`, the lowest eigenvalue of Weil's quadratic form for the Riemann zeta function on
`L^2(-L, L)` is simple, and its eigenfunction is even.

**Please read this first.**

- The work is **unrefereed**. Nobody outside the project has checked the proofs or the code.
- It was produced **largely by an AI assistant** (Anthropic's Claude, used through Claude Code). The author directed
  the work but is not a professional mathematician. See [AI assistance](#ai-assistance).
- It is computer-assisted. Its correctness depends on the code in `code/`, on Arb/FLINT, and on the lemmas in the
  paper.
- It is **not** a result about the Riemann hypothesis and brings no proof of RH closer in any quantified sense. In the
  Connes–van Suijlekom / Connes–Consani–Moscovici programme, "the lowest eigenvalue of `Q_L` is simple with an even
  eigenfunction" is one hypothesis. Assumed for a given `L`, it implies that the Fourier transform of that ground state
  has only real zeros. Reaching RH would also need that hypothesis for every `L` and control of the limit `L -> inf`
  (that those zeros converge to the zeta zeros). Neither is addressed here. We verify the hypothesis only on the
  bounded range `(0, 1.19]`.

The full write-up is [`paper/NOTE.md`](paper/NOTE.md).

## What is proved

Take real `f` supported in `[-L, L]` (log scale) and form `Q(f) = W(f * f~)`. Here `W` is Weil's explicit-formula
functional. It has the pole terms, the archimedean term `Re psi(1/4 + it/2) - log pi`, and every prime power `n` with
`log n < 2L` (paper §2). Let `A_L` be the self-adjoint operator of the closed form on `L^2(-L, L)`.

**Theorem.** For every `L` in `(0, 1.19]`, the lowest eigenvalue of `A_L` is simple, and its eigenfunction is even.

The proof has two pieces.

- **Piece I, `(0, 0.05]`.** A Lipschitz argument for a rescaled form, anchored at the certified window `0.05`. The
  gaps are `>= 0.4172` (second even eigenvalue) and `>= 0.6624` (odd sector).
- **Piece II, `[0.05, 1.19]`.** 31 grid windows. At each one we certify a two-sided bracket for `lambda_0`, a lower
  bound on the second even eigenvalue and a lower bound on the odd sector. Adjacent windows are glued by min-max
  monotonicity in `L`. Every one of the 30 steps passes. The smallest step ratio is `1.288`, and the smallest past
  `L = 1` is `2.148`.

The prime powers 2, 3, 4, 5, 7, 8 and 9 all enter inside this range, at `L = (1/2) log n`. The next one, `n = 11`,
enters at `1.1989`, just past the range.

The certified numbers are in `data/theorem_bounds.json` (Table 2 of the paper). The code's grid ends at `1.19`;
windows past `1.19` are **not** certified and are not claimed.

## What is new, and what is not

The convention throughout is that `L` is the half-width and `p^k` enters when `log p^k < 2L`.

- **Believed new: `(0.8, 1.19]`.** Our search, last refreshed 2026-10-25, found no public proof or certificate of
  simple-even for any continuum window `L > 0.8`. The only exception is the non-explicit neighbourhood of `0.8` in
  mathia PL-287. Andrade (Zenodo 20887028) and Andrews Jr. (Zenodo 22715669) certify simple-even for *finite* Galerkin matrices at
  `L ≈ 1.28, 1.68, 1.99`; those are not statements about the continuum form. Much of this literature is unrefereed
  Zenodo/GitHub work and appears quickly, so the search may have missed items.
- **Not new: `(0, 0.8]`.** Kramarenko-Byrd (Zenodo 23007100, v1 of 2026-09-27, unrefereed) claims a rigorous proof of
  simple-even for all `0 < L <= 4/5`, in the same convention. We have not checked that proof. Our certificate on
  `(0, 0.8]` is an independent computer-assisted confirmation of that range (independent of Kramarenko-Byrd; not
  methodologically independent of Zhu, whose reduction it uses). Earlier drafts of this work claimed `[0.1406, 1]` as new; we found
  Kramarenko-Byrd afterwards and **withdraw that claim**.
- **Not new: small `L`.**
  - `L <= 0.1406`: an explicit Beurling–Deny / positivity-improving argument (mathia PL-294; cf. Andrade, Zenodo
    20682834).
  - Non-explicit small `L`: Suzuki, arXiv:2606.09096.
  - `L <= 9/50`: an analytic argument is reported by grandchallenge / MATHSOLVE. We have not read it in full or
    verified its window convention.
- **Not new: `L = 0.8`.** Zhu (arXiv:2608.24827) certified simple-even at this single window. Zhu also certified upper
  bounds on `lambda_0` for `L` in `[0.5, 2]`.
- **Not new: the reduction.** The whole frequency-cutoff ("one-stroke") reduction is Zhu's (Thm 1.1, Lemma 3.1). That
  includes the head/tail block estimate, the envelope bound and the normalization. Every window certificate here rests
  on it.
- **Not new: positivity.** Positivity, i.e. `lambda_0 > 0`, is certified or claimed elsewhere on this whole range,
  including past `L = 1`. Sources include jose-paul-b, Zhu, Liu (`L = 1`, `17/16`), Mori (up to `(log 15)/2`), Senke
  (`77/64`) and Elston (per-sector bounds at `1.19`–`1.495`). Apart from Zhu's single window `L = 0.8`, none of these
  positivity results addresses simplicity or parity.

Our own contributions are:

- the certified covering of a continuum of windows by monotonicity;
- the Lipschitz piece at small `L`;
- the numerics that reach `L = 1.19`: a hybrid rigorous spherical Bessel evaluation (classical
  recurrence / continued-fraction scheme, cf. Gautschi 1967), and a sharper cutoff certified by
  interval subdivision (`psi_floor`, paper §4.6);
- write-ups of standard estimates and re-derivations of Zhu's lemmas.

See paper §1.3–1.4 for the full list with citations.

## Caveats

These are taken from paper §7.

1. **Normalization.** Our `Q` follows Zhu §6 and was checked against the zero side of the explicit formula. It is
   matched to the Connes–van Suijlekom / Connes–Consani–Moscovici form `QW_lambda` only up to a positive affine
   rescaling `aQ + b||.||^2`, which does not affect simple-even. It has not been checked line by line. A mismatch in
   the relative weights of the pole, archimedean or prime terms would not be covered.
2. **Self-adjoint realization.** The operator is the one of the closure of `Q` on `C_c^inf(-L, L)`. The trial
   functions are polynomials and are argued to lie in the form domain (paper §4.4).
3. **Python `min` on balls.** The gluing and a few steps in `certify.py` use Python `min` on Arb balls. That is not
   rigorous when balls overlap. In every current use they are separated by many orders of magnitude.
4. **Piece I precision anomaly.** The certified upper end of one integral grows with working precision (paper §7.5).
   We use the largest value. The cause is not explained.
5. **Past `L = 1`.** The step ratios there are smaller (10.00, 6.548, 3.774, 3.933, 2.857, 2.148, 264.7), and the
   sharp-cutoff path is used only there. Zhu retracted an earlier positivity claim at `L = 1.19`, which is exactly our
   last window: our bracket there (`lambda_0 >= 7.823e-48`) implies positivity at that `L`, rests on the unreviewed
   sharp-cutoff and triangular code, and is the window most in need of independent checking. One error constant (Gauss–Legendre, off by one node)
   has already been found and corrected in this project; see the erratum in the paper. The tail bounds, the drop
   bound, the transport precision, the Bessel seeds and `psi_floor` deserve an independent audit. The sharp-cutoff
   path and the windows 1.03–1.19, which carry the novelty claim, have not been through any review round, internal or
   external. Nor has the upper-triangle panel summation in `certify.py` / `node_partial.py` that produced the windows
   1.12–1.19 (byte-identical to the previous full-product code on spot checks at `L = 0.3`, `0.8`, `1.0`).
6. **Sources not fully read.** Liu, the MATHSOLVE argument and the Kramarenko-Byrd proof were not read or checked in
   full. Zhu was read from a text extraction of v2.

## Reproducing

### Requirements

- Python 3.12 with the packages in `code/requirements.txt` (`python-flint==0.9.0`, `mpmath==1.4.1`; `numpy` is
  needed only by `convergence.py`, `gl_correction.py` and `exploratory/`).
- FLINT 3 development headers, MPFR and GMP for the C tools. The certificate used the system FLINT 3.4.
  python-flint 0.9.0 bundles its own FLINT 3.6.0.

```sh
python3.12 -m venv .venv && . .venv/bin/activate
pip install -r code/requirements.txt
cd code && make            # builds vchol and veig
```

### Check the glued theorem from the stored bounds (a few seconds)

```sh
cd code && python theorem.py --check
```

This re-runs Piece I, including the Lipschitz constants at 768 bits, and every Piece II step on
`data/theorem_bounds.json`. It uses the longest prefix of the grid present in the file. It should end with
`THEOREM CERTIFIED for all L in (0, 1.19]` (compare `data/logs/check.log`). It does **not** recompute any window.

### Recompute windows

On a single machine:

```sh
cd code
python certify.py 0.1 0.2                 # single windows (seconds to minutes each at L <= 0.8)
UPTO=0.8 CERT_PREC=768 python theorem.py  # windows 0.05 ... 0.8 into ./theorem_bounds.json, then glue
```

The rerun took 18 minutes for stage 1 on a 20-thread laptop. A run of `certify.py 0.1` takes
about 2 minutes and 0.3 GB of RAM. `theorem.py` (and `certify.py` when run as a script)
use a fixed pool of 20 processes; edit `Pool(20)` there on smaller machines. Windows past `0.8` are far more
expensive. Locally they will take days, so use a cluster.

On a cluster (optional): copy `code/` to each host as `~/weil`, or set `CERT_REMOTE_DIR`. Put a `.venv` inside it with
the same requirements. Make the hosts reachable by key-based ssh (`BatchMode=yes`). The coordinator refuses to use a
host whose `certify.py`, `bessel_rec.py` or `node_partial.py` differs by SHA-256.

```sh
export CLUSTER=user@host1,user@host2,user@host3
UPTO=0.97 CERT_PREC=1536 python theorem.py                                     # 0.84 ... 0.97
CERT_BESSEL_PREC=6144 CERT_MATMUL_PREC=320 CERT_VCHOL_PREC=384 python theorem.py  # 0.98 onward
python validate10.py   # optional: sharp-path recertification of L = 1.0, compared with ./theorem_bounds.json
```

We used a 3-node cluster with 64 cores per node. Observed times on it were:

| Windows | Time |
|---|---|
| `0.84`–`0.97` | about 10–14 min per window |
| `0.98`, `0.99` | about 85 min per window |
| `1.03` | 48 min |
| `1.06`, `1.09` | 3 h 04 min and 3 h 16 min |
| `1.12` | 6 h 53 min (cluster assembly, then about 44 min local certification) |
| `1.15`, `1.18`, `1.19` | pipelined: cluster assembly (`prefetch.py`) about 6 h 33 min, 6 h 58 min and 6 h 59 min; local certification 41, 64 and 79 min, overlapped with assembly of the next window |

Memory per worker grows like `N^2`: the cap in `node_partial.py` budgets about `100 N^2` bytes per worker (the
upper-triangle sums, measured 0.87 GB at `N = 4305`, plus the ~0.75 GB shared-memory file) and limits the number of
workers by available memory (`CERT_NODE_PROCS` overrides it). The cached head matrices in `code/cache/` also grow like
`N^2`: about 160 MB each at `L = 1.0`, 1.06 GB each at `L = 1.19`.

Other environment variables:

- `CERT_MARGIN` (default 1.05) and `CERT_DDIV` (default 100) set the slack of the verified Cholesky.
- The rerun needed `CERT_MARGIN=1.2 CERT_DDIV=8` at `L = 1.0` (paper §7.2).
- `CERT_DIGITS` sets the transport precision.

With the released code, `L = 1.0` is computed by the sharp path, so a fresh run reproduces `data/validate_1.0.json`
at that window rather than the rerun's row.

## File map

```
README.md               this file
LICENSE                 MIT (code)
CITATION.cff            citation metadata (author placeholder)
paper/NOTE.md           the paper (CC BY 4.0, see paper/LICENSE)
code/certify.py         per-window certificate: lower bounds (frequency side, Zhu's reduction), upper bound
                        (position-side Rayleigh quotient), sharp-cutoff path (SHARP_TS, psi_floor), cluster dispatch
code/theorem.py         grid + Piece I + Piece II gluing; --check mode
code/lipschitz.py       Piece I constants
code/bessel_rec.py      rigorous hybrid spherical Bessel evaluation (used for L >= 0.98)
code/node_partial.py    cluster worker (upper-triangle panel sums), memory-capped
code/prefetch.py        assembles head matrices on the cluster ahead of the certifying run (pipelining)
code/vchol.c            verified Cholesky (C, FLINT)
code/veig.c             inverse-iteration estimates (non-rigorous, only to choose what to certify)
code/validate10.py      sharp-path validation at L = 1.0
code/gl_correction.py   historical: Gauss-Legendre erratum correction of the original runs
code/convergence.py     non-rigorous Galerkin matrices (trial vectors in gl_correction.py; exploratory scripts)
code/Makefile, code/requirements.txt
data/                   certified bounds, history, scrubbed logs (see data/README.md)
exploratory/            NON-RIGOROUS scripts behind paper Appendix A and the normalization check
```

## AI assistance

The code, the computations, the paper drafts and the internal "reviews" in the paper's review log were produced
largely by an AI assistant: Anthropic's Claude, used through Claude Code, in sessions directed by the author, who
decides what to publish. The author is not a professional mathematician. The mathematical arguments have been
checked only by AI-run internal reviews, not by an independent human expert.

Errors were found during the work. Known ones are documented in the paper's erratum and §7:

- the Gauss–Legendre off-by-one;
- a thread-unsafe precision change;
- an out-of-memory failure.

There may be others. Independent checking of the code and of the lemmas is welcome. That includes checking against
the cited prior work, especially Zhu and Kramarenko-Byrd.

## Credits

- **A. Weil.** The explicit formula and the quadratic form.
- **A. Connes and W. van Suijlekom** (arXiv:2511.23257) and **A. Connes, C. Consani and H. Moscovici**
  (arXiv:2511.22755). They posed simple-even as the missing step. Also A. Connes, survey arXiv:2602.04022.
- **M. Suzuki** (arXiv:2606.09096). The small-`L` theorem.
- **Zhu** (arXiv:2608.24827). The reduction strategy this proof rests on (one-stroke frequency cutoff, head/tail block
  estimate, envelope Lemma 3.1, normalization). Also the certificate at `L = 0.8` and certified upper bounds on
  `lambda_0` for `L` in `[0.5, 2]`.
- **Kramarenko-Byrd** (Zenodo 23007100). Claimed simple-even on `(0, 0.8]`.
- **Other prior work.**
  - mathia project notes PL-282, PL-287, PL-294 (GitHub `murillo128/mathia`).
  - Andrade (Zenodo 20682834, 20887028).
  - grandchallenge / MATHSOLVE.
  - jose-paul-b (`weil-positivity-certificates`).
  - Liu.
  - Mori (Zenodo 23281236).
  - Senke (Zenodo 22871384), positivity at `L = 77/64`.
  - Andrews Jr. (Zenodo 22715669), finite-matrix simple-even at `L ≈ 1.282`.
  - Elston (`weil-windows`).
  - Elder (Zenodo 23113186).
  - Eberwein (Zenodo 23177438).
  - "Eureka" (arXiv:2608.19047).

  See paper §1.3 and §9 for what each one did.
- **Software.**
  - FLINT/Arb (F. Johansson and the FLINT developers).
  - python-flint (its contributors).
  - GMP and MPFR (used by the C tools).
  - mpmath (F. Johansson et al.). The zeta zeros in `exploratory/zeros.txt` were computed with `mpmath.zetazero`.
  - NumPy.
- **W. Gautschi** (SIAM Review 9, 1967). The continued-fraction scheme behind the Bessel evaluation.
- **A. Connes and C. Consani** (arXiv:2106.01715, Lemma 2.7). Pole splitting (paper Appendix A).
- **Quadrature bound.** L. N. Trefethen, *Approximation Theory and Approximation Practice*, Thm 19.3.

## Author and contact

Johann Alexander Botha, GenusCore (not a professional mathematician). Contact: info@genuscore.com.
Corrections and error reports are very welcome.

## License

- **Code:** MIT; see `LICENSE`.
- **Paper and data:** CC BY 4.0; see `paper/LICENSE`.

## How to cite

See `CITATION.cff`:

> Johann Alexander Botha, *A certified simple, even ground state of the localized Weil quadratic form for all windows
> L in (0, 1.19]*, computer-assisted note and code, version 2026-10-26, unrefereed, produced with AI assistance.

Please also cite Zhu (arXiv:2608.24827), whose reduction the proof uses, and the prior work listed above where
relevant.
