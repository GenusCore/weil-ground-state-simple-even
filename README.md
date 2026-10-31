# Simple, even ground state of the localized Weil form for L in (0, 1.4035]

This is a computer-assisted proof, using Arb/FLINT ball arithmetic, of the following statement. For every window
half-width `L` in `(0, 1.4035]`, the lowest eigenvalue of Weil's quadratic form for the Riemann zeta function on
`L^2(-L, L)` is simple, and its eigenfunction is even.

**Please read this first.**

- The work is **unrefereed**. Nobody outside the project has checked the proofs or the code.
- It was produced **largely by an AI assistant** (Anthropic's Claude, used through Claude Code). The author directed
  the work but is not a professional mathematician. See [AI assistance](#ai-assistance).
- It is computer-assisted. Its correctness depends on the code in `code/` (including Zac Elston's certifier in
  `code/gap/elston/`), on Arb/FLINT, and on the lemmas in the paper.
- **Release 1.1.0 (2026-10-31) is an interim release.** The goal is `(0, 1.5]`; the computations for
  `(1.4035, 1.5]` are still running, and nothing is claimed there. The new part `(1.19, 1.4035]` has not been
  audited (see [Caveats](#caveats)).
- It is **not** a result about the Riemann hypothesis and brings no proof of RH closer in any quantified sense. In the
  Connes–van Suijlekom / Connes–Consani–Moscovici programme, "the lowest eigenvalue of `Q_L` is simple with an even
  eigenfunction" is one hypothesis. Assumed for a given `L`, it implies that the Fourier transform of that ground state
  has only real zeros. Reaching RH would also need that hypothesis for every `L` and control of the limit `L -> inf`
  (that those zeros converge to the zeta zeros). Neither is addressed here. We verify the hypothesis only on the
  bounded range `(0, 1.4035]`.

The paper is [`paper/main.pdf`](paper/main.pdf) (source `paper/main.tex`, `paper/refs.bib`). The note of the
previous release, which covers `(0, 1.19]`, is kept as [`paper/NOTE.md`](paper/NOTE.md).

## What is new in 1.1.0

- The certified range grows from `(0, 1.19]` to `(0, 1.4035]` (**Piece III**, paper §5.4 and Table 4).
- Lower bounds come from Elston's window-compressed prime-norm certifier (`weil-windows`, MIT, upstream commit
  `3cc4cb7b`), run unchanged except for one added option, `--deflate`. It adds `v v^T` (`v` an approximate
  ground vector) to the certified matrix, so the certificate bounds the **second** even eigenvalue (paper
  Lemma 5.14, Corollary 5.15).
- Upper bounds `U(L)` come from Legendre–Galerkin Rayleigh quotients (`K = 640`, 3000 bits) with a rigorous
  Gauss–Legendre bound (paper §7). `U(1.19) <= 8.712e-48` is sharper than Piece II's `1.312e-47`, and the first step
  needs it.
- 8 grid steps `1.19, 1.2205, 1.25, 1.2785, 1.3055, 1.3315, 1.3565, 1.3805, 1.4035` all pass. The smallest step
  ratio is `1.131`, on `[1.2205, 1.25]`, and the odd sector binds everywhere.
- Data point, **not** part of the theorem: the step `[1.4255, 1.447]` also passes (ratio `1.379`). The step
  `[1.4035, 1.4255]` before it waits for the even run at `1.4255`.
- Still running: the even sector at `1.4255`, and both sectors at `1.468`, `1.4885` and `1.5`.
- New files: `code/gap/` (certifier, step test, parameter rule), `code/upper/` (upper-bound code),
  `data/piece3/` (19 certificate files, the upper bounds, the check log and SHA-256 sums), the LaTeX paper.

## What is proved

Take real `f` supported in `[-L, L]` (log scale) and form `Q(f) = W(f * f~)`. Here `W` is Weil's explicit-formula
functional. It has the pole terms, the archimedean term `Re psi(1/4 + it/2) - log pi`, and every prime power `n` with
`log n < 2L` (paper §2). Let `A_L` be the self-adjoint operator of the closed form on `L^2(-L, L)`.

**Theorem.** For every `L` in `(0, 1.4035]`, the lowest eigenvalue of `A_L` is simple, and its eigenfunction is
even.

The proof has three pieces.

- **Piece I, `(0, 0.05]`.** A Lipschitz argument for a rescaled form, anchored at the certified window `0.05`. The
  gaps are `>= 0.4172` (second even eigenvalue) and `>= 0.6624` (odd sector).
- **Piece II, `[0.05, 1.19]`.** 31 grid windows. At each one we certify a two-sided bracket for `lambda_0`, a lower
  bound on the second even eigenvalue and a lower bound on the odd sector. Adjacent windows are glued by min-max
  monotonicity in `L`. Every one of the 30 steps passes. The smallest step ratio is `1.288`, and the smallest past
  `L = 1` is `2.148`.

- **Piece III, `[1.19, 1.4035]`.** 8 grid steps. At each window, certified lower bounds on the second even
  eigenvalue and on the odd sector (Elston's certifier with `--deflate`) and a rigorous upper bound on the lowest
  eigenvalue (Galerkin Rayleigh quotient). Every step passes; the smallest step ratio is `1.131`.

The prime powers 2, 3, 4, 5, 7, 8 and 9 enter inside Piece II, at `L = (1/2) log n`; 11, 13 and 16 enter inside
Piece III (at `1.1989`, `1.2825`, `1.3863`).

The certified numbers are in `data/theorem_bounds.json` (Pieces I–II, Table 2 of the paper) and `data/piece3/`
(Piece III, Table 4). Windows past `1.4035` are **not** part of the theorem.

## What is new, and what is not

The convention throughout is that `L` is the half-width and `p^k` enters when `log p^k < 2L`.

- **Believed new: `(0.8, 1.4035]`.** Our search (arXiv, Zenodo, GitHub and rStation; refreshed 2026-10-31) found no
  public proof or certificate of simple-even for any continuum window `L > 0.8`. The only exception is the non-explicit neighbourhood of `0.8` in
  mathia PL-287. Andrade (Zenodo 20887028) and Andrews Jr. (Zenodo 22715669) certify simple-even for *finite* Galerkin matrices at
  `L ≈ 1.28, 1.68, 1.99`; those are not statements about the continuum form. Much of this literature is unrefereed
  Zenodo/GitHub work and appears quickly, so the search may have missed items. On rStation, lihongmin
  (riemann_hypothesis paper 45 v1, published 2026-10-23) and bramble_trace (paper 60 v1, published 2026-10-26)
  certify positivity of both parity sectors on `(0, (1/2) log 11]` and `(0, 2]` respectively; both use our form, and
  their support radius `a` is our `L`. Neither identifies the ground state or compares the sectors (paper 45 §1: "our
  lower certificate neither identifies the ground state nor estimates the sharp minimum"). What we claim as new on
  `(0.8, 1.4035]` is therefore simple-even itself (the ground state is simple and even) and the certified lower
  bounds on the second even eigenvalue `lambda_2^even` that prove it (Piece II, and the rank-one deflation of
  Piece III), not positivity.
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
  (`77/64`), Elston (per-sector bounds at `1.19`–`1.495`), lihongmin (rStation riemann_hypothesis paper 45 v1, both
  sectors, `(0, (1/2) log 11]`) and bramble_trace (rStation paper 60 v1, both sectors, `(0, 2]`, full-source Schur
  certificate). Apart from Zhu's single window `L = 0.8`, none of these positivity results addresses simplicity or
  parity.

Our own contributions are:

- the certified covering of a continuum of windows by monotonicity;
- the Lipschitz piece at small `L`;
- the rank-one deflation that turns Elston's lowest-eigenvalue certifier into a second-eigenvalue certifier, and
  the Piece III covering built on it;
- the numerics that reach `L = 1.19`: a hybrid rigorous spherical Bessel evaluation (classical
  recurrence / continued-fraction scheme, cf. Gautschi 1967), and a sharper cutoff certified by
  interval subdivision (`psi_floor`, paper §5.2);
- write-ups of standard estimates and re-derivations of Zhu's lemmas.

See paper §1.3–1.4 for the full list with citations.

## Caveats

These are taken from paper §10.

1. **Normalization.** Our `Q` follows Zhu §6 and was checked against the zero side of the explicit formula. It is
   matched to the Connes–van Suijlekom / Connes–Consani–Moscovici form `QW_lambda` only up to a positive affine
   rescaling `aQ + b||.||^2`, which does not affect simple-even. It has not been checked line by line. A mismatch in
   the relative weights of the pole, archimedean or prime terms would not be covered.
2. **Self-adjoint realization.** The operator is the one of the closure of `Q` on `C_c^inf(-L, L)`. The trial
   functions are polynomials and are argued to lie in the form domain (paper §5.2).
3. **Python `min` on balls.** The gluing and a few steps in `certify.py` use Python `min` on Arb balls. That is not
   rigorous when balls overlap. In every current use they are separated by many orders of magnitude.
4. **Piece I precision anomaly.** The certified upper end of one integral grows with working precision (paper Remark 5.2).
   We use the largest value. The cause is not explained.
5. **Past `L = 1`.** The step ratios there are smaller (10.00, 6.548, 3.774, 3.933, 2.857, 2.148, 264.7), and the
   sharp-cutoff path is used only there. Zhu retracted an earlier positivity claim at `L = 1.19`, which is the last
   window of Piece II: our bracket there (`lambda_0 >= 7.823e-48`) implies positivity at that `L`, rests on the unreviewed
   sharp-cutoff and triangular code, and is the window most in need of independent checking. Positivity at `L = 1.19`
   is independently certified by rStation papers 45 and 60; what rests on our unreviewed code there is the
   simple-even bracket, not positivity. One error constant (Gauss–Legendre, off by one node)
   has already been found and corrected in this project; see Remark 5.9 of the paper. The tail bounds, the drop
   bound, the transport precision, the Bessel seeds and `psi_floor` deserve an independent audit. The sharp-cutoff
   path and the windows 1.03–1.19, which carry the novelty claim, have not been through any review round, internal or
   external. Nor has the upper-triangle panel summation in `certify.py` / `node_partial.py` that produced the windows
   1.12–1.19 (byte-identical to the previous full-product code on spot checks at `L = 0.3`, `0.8`, `1.0`).
6. **Piece III is unaudited beyond `L = 1.19`.** Elston's certifier was audited independently only at `L = 1.19`.
   The runs at larger `L`, the per-window parameters (`beta_K`, `T#`, `N`, Gauss points) and our `--deflate`
   addition have not been reviewed by anyone. The step ratios are small (`1.131`–`1.552`), so an error of about 13%
   in one bound could break a step. The lower bounds (Elston's code) and the upper bounds (our code) must evaluate
   the same form with the same norm (paper §5.4, "One form, two codes"). The revision of
   `code/upper/fast_recurrence.py` used for the `U` runs was not recorded; re-running `U(1.19)` with the shipped code
   reproduces the stored line exactly (`data/piece3/rerun_U1.19.jsonl`).
7. **Sources not fully read.** Liu, the MATHSOLVE argument and the Kramarenko-Byrd proof were not read or checked in
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

### Check Piece III from the stored certificates (seconds)

```sh
(cd data/piece3 && sha256sum -c SHA256SUMS)
python code/gap/check_steps.py
```

This compares, in Arb, `min(lambda_2^even, lambda_1^odd)(L_{i+1})` with `U(L_i)` for every step. It reads
`data/piece3/` and should end with `CERTIFIED: ... simple-even for all L in (0, 1.4035]` (compare
`data/piece3/check.log`). Together with `theorem.py --check` above this is the whole theorem.

### Recompute Piece III

Each certificate file stores its options in the field `args`. In `code/gap/elston/`:

```sh
C="--panels-per-unit 4 --schur-k 4000 --schur-iters 1500"
W="--a 1.2205 --delta 0.3 --beta-k 64 --tsharp 850 --nmodes 849 --gl 56 --prec 400"   # from g_1.2205_*.json
python grid_certificate_rigorous.py $W $C --sector odd --workers 28 --json g_1.2205_odd.json
python grid_certificate_rigorous.py $W $C --sector even --deflate --workers 28 --json g_1.2205_l2.json
```

The windows took 12 minutes to 7 hours each with 15–52 workers on 64-vCPU machines (paper Table 3). Memory
grows like `workers x N^2`: about 2 GB per worker at `N = 1757`. Upper bounds, in `code/upper/` (about 4–14 minutes
per window):

```sh
THREADS=8 python fr_sweep.py upper.jsonl "640" "1.19,1.2205,1.25" 3000
```

### Recompute Pieces I–II windows

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
- The rerun needed `CERT_MARGIN=1.2 CERT_DDIV=8` at `L = 1.0` (`paper/NOTE.md` §7.2).
- `CERT_DIGITS` sets the transport precision.

With the released code, `L = 1.0` is computed by the sharp path, so a fresh run reproduces `data/validate_1.0.json`
at that window rather than the rerun's row.

## File map

```
README.md               this file
REVIEWING.md            how to review: weakest links and checks by effort
LICENSE                 MIT (code)
CITATION.cff            citation metadata
paper/main.tex, main.pdf the paper, release 1.1.0 (CC BY 4.0, see paper/LICENSE); refs.bib
paper/NOTE.md           the note of the (0, 1.19] release of 2026-10-26
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
code/gap/elston/        Elston's certifier (MIT, own LICENSE, UPSTREAM_COMMIT) with the added --deflate option
code/gap/check_steps.py Piece III step test (gluing)
code/gap/params.py      design rule for the Piece III window parameters (not used by the proof)
code/upper/             rigorous Legendre-Galerkin upper bounds (fast_recurrence.py, hp.py, fr_sweep.py)
data/                   certified bounds, history, scrubbed logs (see data/README.md)
data/piece3/            Piece III certificates, upper bounds, check log, SHA256SUMS
exploratory/            NON-RIGOROUS scripts behind `paper/NOTE.md` Appendix A and the normalization check
```

## AI assistance

The code, the computations, the paper drafts and the internal "reviews" in the paper's review log were produced
largely by an AI assistant: Anthropic's Claude, used through Claude Code, in sessions directed by the author, who
decides what to publish. The author is not a professional mathematician. The mathematical arguments have been
checked only by AI-run internal reviews, not by an independent human expert.

Errors were found during the work. Known ones are documented in the paper (Remark 5.9, §10) and in `paper/NOTE.md` (erratum, §7):

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
  - Elston (`weil-windows`, Zenodo 10.5281/zenodo.23174403). His certifier produces every Piece III lower bound.
  - Elder (Zenodo 23113186).
  - Eberwein (Zenodo 23177438).
  - "Eureka" (arXiv:2608.19047).
  - lihongmin, "Cutoff-free positivity of localized Weil forms through the (log 11)/2 threshold", rStation
    riemann_hypothesis paper 45, v1, published 2026-10-23T09:48:37Z,
    https://rstation.ai/problems/riemann_hypothesis/papers/45/v1. Positivity of both parity sectors on
    `(0, (1/2) log 11]`.
  - bramble_trace, "Full-source Schur certificates for Weil positivity at support radius two", rStation
    riemann_hypothesis paper 60, v1, published 2026-10-26T19:12:17Z,
    https://rstation.ai/problems/riemann_hypothesis/papers/60/v1. Positivity of both parity sectors on `(0, 2]`.

  See paper §1.3 for what each one did.
- **Software.**
  - FLINT/Arb (F. Johansson and the FLINT developers).
  - python-flint (its contributors).
  - GMP and MPFR (used by the C tools).
  - mpmath (F. Johansson et al.). The zeta zeros in `exploratory/zeros.txt` were computed with `mpmath.zetazero`.
  - NumPy.
- **W. Gautschi** (SIAM Review 9, 1967). The continued-fraction scheme behind the Bessel evaluation.
- **A. Connes and C. Consani** (arXiv:2106.01715, Lemma 2.7). Pole splitting (`paper/NOTE.md` Appendix A).
- **Quadrature bound.** L. N. Trefethen, *Approximation Theory and Approximation Practice*, Thm 19.3.

## Author and contact

Johann Alexander Botha, GenusCore (not a professional mathematician). Contact: info@genuscore.com.
Corrections and error reports are very welcome.

## License

- **Code:** MIT; see `LICENSE`. `code/gap/elston/` is Zac Elston's code under his MIT licence (`code/gap/elston/LICENSE`).
- **Paper and data:** CC BY 4.0; see `paper/LICENSE`.

## How to cite

See `CITATION.cff`:

> Johann Alexander Botha, *A certified simple, even ground state of the localized Weil quadratic form for all windows
> L in (0, 1.4035]*, computer-assisted paper and code, version 1.1.0 (2026-10-31), unrefereed, produced with AI
> assistance.

Please also cite Zhu (arXiv:2608.24827), whose reduction the proof uses, Elston (`weil-windows`), whose certifier
Piece III uses, and the prior work listed above where
relevant.
