# Lean pilot: the matrix step of two certified windows

This directory is a pilot. It checks one step of the proof in
[`paper/NOTE.md`](../../paper/NOTE.md) in Lean 4 (with Mathlib) on the real data: the linear
algebra behind the eigenvalue lower bounds at the two smallest windows, `L = 0.05` and
`L = 0.1`. Main project:
[GenusCore/weil-ground-state-simple-even](https://github.com/GenusCore/weil-ground-state-simple-even).

Everything else in the proof is outside this pilot. The analytic inputs that the checked step
uses appear as explicit hypotheses (`HeadHyp`), not as axioms.

## What Lean checks

For each window there are three `97 × 97` certificates: the even head matrix `A`, the
deflated `A + vvᵀ` (for `λ₁^even`), and the odd head matrix `A'`. `export_window.py` (untrusted)
rounds the cached Arb balls of the head matrices outward to boxes
`|A_ij − C_ij/2⁶⁴| ≤ rad/2⁶⁴` with integer `C_ij`, and supplies an approximate Cholesky factor
`R`, a shift `M`, a residual budget `K` and a lower bound `μ`.

`check` ([`Checker.lean`](SpectralCertificate/Checker.lean)) computes the residual
`C + vvᵀ − M·I − RRᵀ` exactly in `Int` and requires `|residual_ij| + rad ≤ K` and
`μ + n·K ≤ M`. The kernel runs it via `decide +kernel` (no `native_decide`, no compiled
code). Its soundness is proved once, for every `n` and every real matrix in the box:

```lean
theorem check_sound (c : Cert) (h : check c = true) (A : Fin c.n → Fin c.n → ℝ) (hA : InBox c.C c.rad A)
    (x : Fin c.n → ℝ) :
    (c.mu : ℝ) / 2 ^ 64 * ∑ i, x i ^ 2 ≤ quad A x + (∑ i, vec c.n c.v i * x i) ^ 2
```

With the deflation vector, `check_sound_perp` gives `A ≥ μ` on `v⊥`.
[`Window.lean`](SpectralCertificate/Window.lean) turns the three certificates of a window
into eigenvalue lower bounds, assuming `HeadHyp`:

```lean
def HeadHyp (w : Window) (lam0 lam1 lamOdd : ℝ) : Prop :=
  ∃ (A A' : Fin w.n → Fin w.n → ℝ) (β eD eB eD' eB' : ℝ),
    InBox w.Ce w.rade A ∧ InBox w.Co w.rado A' ∧
    (w.betaLo : ℝ) ≤ β ∧ eD ≤ w.eDe ∧ eB ≤ w.eBe ∧ eD' ≤ w.eDo ∧ eB' ≤ w.eBo ∧
    (∀ μ : ℝ, (∀ x, μ * ∑ i, x i ^ 2 ≤ quad A x) → min μ (β - eD) - eB ≤ lam0) ∧
    (∀ μ : ℝ, (∀ x, ∑ i, vec w.n w.v i * x i = 0 → μ * ∑ i, x i ^ 2 ≤ quad A x) →
      min μ (β - eD) - eB ≤ lam1) ∧
    (∀ μ : ℝ, (∀ x, μ * ∑ i, x i ^ 2 ≤ quad A' x) → min μ (β - eD') - eB' ≤ lamOdd)

theorem window_sound (w : Window) (he : check w.even = true) (hd : check w.deflated = true)
    (ho : check w.odd = true) {lam0 lam1 lamOdd : ℝ} (h : HeadHyp w lam0 lam1 lamOdd) :
    bound w.fe.mu w.betaLo w.eDe w.eBe ≤ lam0 ∧ bound w.fd.mu w.betaLo w.eDe w.eBe ≤ lam1 ∧
      bound w.fo.mu w.betaLo w.eDo w.eBo ≤ lamOdd
```

`bound mu β e_D e_B = min(mu/2⁶⁴, β − e_D) − e_B`. The six `*_ok` theorems in
`Data/W005.lean` and `Data/W01.lean` are the kernel checks of the real certificates.
[`Rows.lean`](SpectralCertificate/Rows.lean) proves (`W005_rows`, `W01_rows`) that the
resulting bounds imply the stored lower ends of `λ₁^even` and `λ₀^odd` in
`data/theorem_bounds.json` at these two windows; `test_certificate.py` checks that the
constants in `Rows.lean` equal the JSON values exactly. Each `μ` sits about 4.8% below the
eigenvalue estimate, and the Lean bounds exceed the stored lower ends by `1e-20` to `1.5e-19`:

| L | certificate | μ (Lean) | estimate | stored lower bound |
|---|---|---|---|---|
| 0.05 | even `A` (λ₀) | 0.9351626876 | 0.981921 | 0.93516268756623234487 ± 3.43e-21 |
| 0.05 | deflated (λ₁^even) | 1.8875580286 | 1.981936 | 1.8875580286496939558 ± 2.46e-20 |
| 0.05 | odd `A'` (λ₀^odd) | 1.9269729395 | 2.023322 | 1.9269729394615598321 ± 2.89e-21 |
| 0.1 | even `A` | 0.4348359561 | 0.456578 | 0.43483595610756670056 ± 4.83e-21 |
| 0.1 | deflated | 1.3872336837 | 1.456595 | 1.3872336837137494079 ± 4.89e-20 |
| 0.1 | odd `A'` | 1.2664541279 | 1.329777 | 1.2664541279327078659 ± 4.40e-20 |

Each window is certified at the exact dyadic point `L̃ = arb(L).upper()` (768 bits,
`L ≤ L̃ < L + 2⁻⁷⁰⁰`), stored in `certificates.json`.

## What is trusted

The hypotheses of `HeadHyp`, and the tools:

- that the Arb balls in the head-matrix cache enclose the true head matrices at `L̃`
  (`InBox`): the enclosures of `certify.py`, Lemmas (B), (M), (J), (G), Arb and FLINT.
  SHA-256 sums of the cache files are in `certificates.json`;
- the constants `β`, `e_D`, `e_B` (rational bounds in each data file), from Lemmas (E), (R),
  (S) and the tail estimates;
- Lemmas (R), (S), (D) as operator statements: the three implications in `HeadHyp`
  (min-max, the head-tail reduction, rank-one interlacing);
- Lean's kernel and Mathlib. The axiom audit allows only `propext`, `Classical.choice` and
  `Quot.sound`.

Lean does not define Weil's form; `lam0`, `lam1`, `lamOdd` are arbitrary reals, and Lean
does not prove that the hypotheses are satisfiable.

## Timings

Local run, details in [`VALIDATION.txt`](VALIDATION.txt): clean `lake build` 92 s with
0 warnings (about 38 s per data file, so about 10 s per certificate; peak about 7.4 GB RSS per
Lean process, of which about 5.1 GB is memory-mapped Mathlib); `check_lean.py` 63 s;
29 Python tests in under 1 s. The GitHub workflow runs the same steps; it has not run yet.

## Run it

Lean `v4.19.0` and Mathlib `v4.19.0` are pinned; `lake-manifest.json` is committed.

```sh
python3 checker.py            # exact Python preflight of the data files (stdlib only)
python3 -m unittest -v        # adversarial and synchronisation tests
lake exe cache get            # Mathlib build cache
lake build                    # kernel-checks the six certificates
python3 check_lean.py         # build, axiom audit, boundary acceptance, 4 tampered inputs rejected
```

`check_lean.py` rejects forbidden keywords (`sorry`, `native_decide`, `axiom`, `extern`, …),
prints the axioms of 11 theorems, asks Lean to accept `μ = M − n·K` exactly, and asks it to
reject four tampered inputs (a shifted matrix entry, a radius one unit `2⁻⁶⁴` past the budget,
a halved factor entry, `μ` one unit past `M − n·K`). On success it writes
`validation-local.json`.

To regenerate the data files (needs the Weil code with its head cache and `python-flint`):

```sh
WEIL_RERUN=/path/to/weil-mechanism/rerun /path/to/venv/bin/python export_window.py
```

## Files

| File | Role |
|---|---|
| `SpectralCertificate/Checker.lean` | Integer checker and its soundness proof (`check_sound`, `check_sound_perp`) |
| `SpectralCertificate/Window.lean` | Packed-data decoders, `HeadHyp`, `window_sound` |
| `SpectralCertificate/Data/W005.lean`, `W01.lean` | Generated: boxes, factors, constants; three `decide +kernel` checks each |
| `SpectralCertificate/Rows.lean` | `W005_rows`, `W01_rows`: link to the stored table |
| `export_window.py` | Untrusted producer (needs Arb and the Weil cache) |
| `checker.py` | Exact Python mirror of `check`, data emission and parsing |
| `check_lean.py`, `test_certificate.py` | Lean audit and negative tests; Python tests |
| `certificates.json` | Provenance: cache hashes, constants, `μ`, dyadic window points |

## Scope / next steps

This pilot covers 2 of the 27 certified windows (the certified range currently reaches
`L = 1.09`; the run is extending toward `1.19`), and only the matrix step. Planned future work:

- **All windows.** Near `L = 1` the head matrices have `N` up to about 4,500 per sector,
  dense factors, and margins of about `10⁻³⁰`, so the residual needs integers of 100+ bits and
  about `N³/6` multiply-adds per certificate. That needs new performance work beyond kernel
  `decide` (efficient verified arithmetic, structured or blocked certificates).
- **The analytic lemmas.** Formalizing (R), (S), (D), (E) and the tail bounds, so that the
  hypotheses above become theorems.
- **The full theorem.** Formalizing how the windows combine into a machine-checked
  conditional theorem over the whole certified range.
- **A reusable certificate library** for PSD and eigenvalue certificates of this kind.

Code in this directory uses the repository's MIT licence. It was prepared with AI assistance
and has not received independent expert review.
