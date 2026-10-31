# Reviewing and auditing this work

This is a computer-assisted proof that has **not been peer-reviewed**. All checking so far was done by the
author with AI assistance: independent recomputation, two separate code paths, and AI-run review rounds
(see the paper, §"Verification and audits"). No independent human expert has checked the proofs or the code.

**Independent review is welcome from anyone.** You do not need to check everything: one verified step,
one reproduced window, or one well-founded objection is a real contribution and will be credited
(by name or anonymously, as you prefer) in the paper and in this file.

## How to report

- Open a GitHub issue, with the label `review` if you can. Please say which section, file or table it
  concerns, what you checked, and what you found, including "checked, no problem found".
- Or write to info@genuscore.com.
- Errors are the most valuable kind of feedback. Confirmed errors are fixed in public, with credit,
  and listed in the paper's erratum.

## Where to look first: the weakest links

These are the parts we consider least checked, roughly in order of how much the result depends on them:

1. **Piece III (L in (1.19, 1.4035], new in release 1.1.0).** It uses Elston's certifier
   (`code/gap/elston/`, MIT), which was independently audited only at L = 1.19, and our `--deflate` addition
   (about 15 lines, option `--deflate` in `grid_certificate_rigorous.py`). Neither the larger-L runs, nor the
   per-window parameter changes (beta_K, T#, N, Gauss points), nor the deflation have been audited. The step
   ratios are small (1.131 to 1.552), so a 13% error in one bound could break a step.
2. **One form, two codes.** The step test compares Elston's lower bounds with our upper bounds. Both must
   evaluate exactly the same quadratic form with the same norm (paper, §5.4, "One form, two codes").
3. **The sharp-cutoff path and the windows 1.03–1.19** of Piece II (paper §5.2 and the README), which have
   not been through any review round.
4. **The rigorous Gauss–Legendre bound** in `code/upper/fast_recurrence.py`, `quad_bounds` (Bernstein
   ellipse, n-point rule exponent ρ^(−2n+2)), which all upper bounds U(L) rely on. The revision of that file
   used for the U runs was not recorded (paper, Remark 7.2).

## Checks by effort

### Minutes (no heavy computation)

- **Re-check the gluing.** `python code/theorem.py --check` (Pieces I–II, `data/theorem_bounds.json`) and
  `python code/gap/check_steps.py` (Piece III, `data/piece3/`); confirm every step ratio
  min(λ₂^even, λ₁^odd)(L_{i+1}) / U(L_i) exceeds 1, and compare with `data/logs/check.log` and
  `data/piece3/check.log`. `sha256sum -c data/piece3/SHA256SUMS` (run in `data/piece3/`) checks the files.
- **Read the deflation lemma and its corollary** (paper, Lemma 5.14 and Corollary 5.15, half a page).
- **Read the monotonicity lemma and the step criterion** (paper §3).
- **Spot-check certificate files** against the tables in the paper (every number in a table is traceable to a
  file listed in the reproducibility table).

### Hours (one machine)

- **Re-run one Piece III window**, e.g. L = 1.2205 (a few core-hours), with and without `--deflate`, and
  compare with the shipped JSON.
- **Re-run one upper bound U(L)** with `code/upper/fr_sweep.py` (K = 640, 3000 bits, about 12 minutes) and
  compare `lam0_upper` with `data/piece3/upper640.jsonl`.
- **Re-run one Piece II window** with `code/certify.py`.

### Days (cluster or patience)

- Full re-run of Piece III up to 1.4035 (about 800 worker-process hours; parameters in each JSON's `args`).
- A line-by-line audit of Elston's certifier at larger L, or of our upper-bound code.
- An independent implementation of any one component (the best check of all).

## What a review is not expected to do

You do not need to assess the Riemann hypothesis, the Connes–van Suijlekom programme, or the novelty
claims; those are stated with sources in the paper's introduction. The question is only:
**are the stated bounds correct and does the argument prove the theorem as stated?**
