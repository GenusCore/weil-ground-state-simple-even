# A certified simple, even ground state of the localized Weil quadratic form for all windows $L \in (0, 1.19]$

**Johann Alexander Botha** (GenusCore; contact: info@genuscore.com)

*Computer-assisted note. Not refereed. Release version of 2026-10-26 (internal reviews of 2026-10-23/24, extension past $L = 1$ of 2026-10-24/25, extension to $L = 1.19$ of 2026-10-25/26). The sharp-cutoff path (§4.6), the windows 1.03–1.19 and the triangular assembly code used for 1.12–1.19 have not been through any review round, internal or external.*
*Code, data and logs: this repository (`code/`, `data/`, `data/logs/`). The certificate of record is `data/theorem_bounds.json`; `python code/theorem.py --check` re-runs Piece I and Piece II on it without recomputing any window (§6).*
*File names.* Sections that describe the history of the certificate name development files that are not part of this release, for example `run_theorem*.log`, `run_step3*.log`, `cache/`, `NOVELTY.md` (our prior-art survey), `REDUCTION.md` and the directories `weil-mechanism/`, `rerun/`, `extend/`. Released counterparts: the `rerun/` and `extend/` code is `code/` (the `extend/` revision); `rerun/run_rerun.log`, `extend/run_extend.log` and `run_gl_correction.log` are in `data/logs/`; the `theorem_bounds.json` of the original runs is `data/history/theorem_bounds_original.json`, and `theorem_bounds_corrected.json` is in `data/history/`; `weil_z.py` and the Appendix A scripts are in `exploratory/`.
*AI-assistance disclosure: this note, the certification code, the correction script, and every internal review of them (including the reviews in the Review log) were produced largely by an AI assistant (Anthropic's Claude, used through Claude Code), directed by the author, who is not a professional mathematician. No human referee has checked any of it.*

---

## Abstract

We consider Weil's quadratic form $Q$ for the Riemann zeta function, restricted to real test functions supported in a window $[-L, L]$ on the logarithmic scale. The form keeps the pole terms, the archimedean term, and every prime power $n$ with $\log n < 2L$. We prove that for **every** $L \in (0, 1.19]$ the lowest eigenvalue of the associated self-adjoint operator is simple and its eigenfunction is even. The proof is computer-assisted and has two parts.

- **Piece I, $(0, 0.05]$.** A Lipschitz estimate for a rescaled form. It uses explicit derivative bounds and one certified integral, and it is anchored at the certified window $L = 0.05$.
- **Piece II, $[0.05, 1.19]$.** A covering by 31 grid windows, glued together by min-max monotonicity in $L$.

At each grid window we certify three things in Arb ball arithmetic:

- a two-sided bracket for $\lambda_0$;
- a lower bound for the second even eigenvalue $\lambda_1^{\mathrm{even}}$;
- a lower bound for the bottom of the odd sector, $\lambda_0^{\mathrm{odd}}$.

**Small windows.** For $L \le (1/2)\log\rho_{\mathrm{plastic}} \approx 0.1406$ the statement was already available. The unrefereed notes [PL294] give an explicit Beurling–Deny (positivity-improving) argument. The unrefereed [An1] states the same threshold for the full form in a remark; its theorem covers only the pole-removed form. [Su] covers a non-explicit range of small $L$.

**What is new, and what is not.** Zhu certified simple-even at the single window $L = 0.8$ [Zhu, Thm 6.2] (§1.3). Kramarenko-Byrd [KB] (Zenodo deposit of 2026-09-27; unrefereed) claims a rigorous proof of simple-even for **all** $0 < L \le 4/5$, in the same half-width convention. Earlier drafts of this note, written before we found [KB], claimed the interval $[0.1406, 1]$ as new; **that claim is withdrawn.** The part we believe to be new is $(0.8, 1.19]$: our search (§1.3; last refreshed 2026-10-25) found no public proof or certificate of simple-even for any window $L > 0.8$ of the continuum form, other than the non-explicit neighbourhood of $0.8$ in [PL287]. In this range the prime powers $5$, $7$, $8$ and $9$ enter (at $L = 0.8047$, $0.9730$, $1.0397$, $1.0986$). On $(0, 0.8]$ our certificate is an independent computer-assisted confirmation of the range claimed in [KB] (independent of [KB]; not methodologically independent of [Zhu], whose reduction it uses). Every window certificate rests on Zhu's frequency-cutoff ("one-stroke") reduction, including its head/tail block estimate [Zhu, Thm 1.1].

**Reproducibility status.**

- **Certificate of record.** `data/theorem_bounds.json` holds the 31 windows. The 24 windows $L \le 1$ come from a clean end-to-end rerun with a single code revision (2026-10-24): empty cache, empty bounds file, `vchol`/`veig` rebuilt from source, SHA-256 check of the code on every cluster host (log `data/logs/run_rerun.log`). The windows $1.03$, $1.06$, $1.09$ come from the extension run (2026-10-24/25), and $1.12$, $1.15$, $1.18$, $1.19$ from its continuation (2026-10-25/26) with the triangular head assembly (Erratum, further issues, item 4); both use the sharp-cutoff code path of §4.6 and share the log `data/logs/run_extend.log`. `python code/theorem.py --check` re-runs Piece I and all 30 Piece II steps on this file and ends with "THEOREM CERTIFIED for all L in (0, 1.19]" (`data/logs/check.log`). Table 2 (§5) lists the numbers.
- **History.** The windows $L \le 1$ were first certified in several runs by successive revisions of the code, with a Gauss–Legendre error constant that was off by one node. That certificate was corrected by a Weyl-shift argument (see Erratum and correction). Table 1 and many numbers quoted in §§4–6 are from the corrected original certificate (`data/history/theorem_bounds_corrected.json`). The rerun supersedes it and does not depend on it; the two agree in the lower bounds to about three significant digits, except at $L = 1$ (§7.2).

A further caveat: our $Q$ is identified with the Connes–van Suijlekom / CCM form only up to positive affine rescaling (§7.3).

---

## Erratum and correction: Gauss–Legendre error constant

**The bug.** The original runs of `certify.py` applied the bound $\frac{64}{15}M\rho^{-2n}/(\rho^2-1)$ [Tr, Thm 19.3] to the **$n$-point** Gauss–Legendre rule. [Tr] states that bound for the **$(n+1)$-point** rule. The correct bound for $n$ points is $\frac{64}{15}M\rho^{-2n+2}/(\rho^2-1)$. The bug affected two places:

- **Head panels.** These use $NG = 64$, with $\rho \approx 2.0806$. The quadrature ball was too small by the factor $\rho^2 \approx 4.33$.
- **Rayleigh outer rule.** This uses $n_{\mathrm{out}} = 450$ with $\rho = 3$, so the error term was too small by a factor of 9.

**The fix in the code.** `certify.py` now uses `rho**(-2*(NG-1))` and `rho**(-2*(n_out-1))`.

**How the certificate was corrected** (`gl_correction.py`, `run_gl_correction.log`, about 76 s wall time). No head matrix was rebuilt.

1. **Lower bounds** ($\lambda_0$, $\lambda_1^{\mathrm{even}}$, $\lambda_0^{\mathrm{odd}}$).
   - Let $q_{\mathrm{old}}$ be the per-entry quadrature radius the original runs used. The true head entries lie within $\rho^2 q_{\mathrm{old}}$ of the midpoints.
   - So the true head matrix lies within a further $(\rho^2-1)q_{\mathrm{old}}$ per entry of the certified ball. That is a perturbation of spectral norm at most $\delta = N(\rho^2-1)q_{\mathrm{old}}$.
   - By Weyl, every certified head eigenvalue bound drops by at most $\delta$. So does the final bound $\min(\mu, \beta - e_D) - e_B$. This also covers the deflated head $A + vv^T$, which carries the same quadrature ball.
   - We recomputed $\delta$ per window and per parity from the code's formula and subtracted it from every lower bound. The window parameters ($T^\#$, $N$, $\beta$, $M$) entering $q_{\mathrm{old}}$ are the ones the **final** code computes. We assume they equal those of the original runs (§7.2).
   - The per-window values are in Table 1. They are identical for the two parities to the digits shown. Their range is:
     - $\le 6.4\cdot10^{-35}$ for $L \le 0.8$;
     - $\le 1.6\cdot10^{-33}$ for $0.84 \le L \le 0.97$;
     - $\le 3.2\cdot10^{-32}$ for $0.98 \le L \le 1$.
2. **Upper bounds $U$.** These were recomputed from scratch with the fixed `rayleigh_upper` at every window.
   - **Trial vectors.** Each window uses a fresh trial vector: the even ground vector of a non-rigorous Galerkin matrix (`convergence.galerkin`). Its degree is $D_G = \min(2\lceil eL\cdot2\pi e^{2L}/8\rceil + 60,\ 298)$, which on the grid runs from $62$ ($L = 0.05$) to $92$ ($L = 1$). The vector therefore has $D_G/2 \le 46$ even coefficients, and the trial polynomial has degree $\le 90$. The code's caps of degree 298 and 150 coefficients never bind. Any trial vector gives a rigorous upper bound, so the choice affects only sharpness.
   - **Comparison with the old $U$.** The new $U$ is used at all windows. It is smaller than the old one at every window except $L = 0.99$, where it is slightly larger: $3.640\cdot10^{-29}$ against the old $3.551\cdot10^{-29}$ (log midpoint). Hence every bracket became narrower or stayed the same, except at $L = 0.99$.
3. **Re-check of the theorem.**
   - Piece I was re-run with freshly computed Lipschitz constants at 768 bits: $C_{\mathrm{even}} \le 4.876362$ and $C_{\mathrm{odd}} \le 0.761666$.
   - All 23 Piece II steps were re-checked on the corrected numbers.

**Result.**

- The log ends with "THEOREM CERTIFIED (corrected) for all L in (0, 1.0]".
- Piece I gaps: $\ge 0.4173$ (even) and $\ge 0.6625$ (odd).
- Step ratios: $1.288$ to $380.0$.
- **Status.** The correction is rigorous given the parameter assumption above. That assumption is checked for $L \ge 0.98$ only as far as $N$ goes ($T^\#$ enters through $N$): the cached head files `cache/head_<L>_<parity>_<N>.txt` from the original runs carry $N = 2447, 2471, 2496$, and these are the values the final code computes. It is not checked for $\beta$, $M$, or for $L \le 0.97$. Under this assumption the theorem holds with the corrected numbers.
- Numbers labelled "corrected" in this note are from this corrected certificate. The certificate of record is the rerun plus the extension (Table 2).

**Further issues found during the extension past $L = 1$.** None of them affects a certified bound.

1. **Process-global precision and threads.** python-flint's working precision `ctx.prec` is global to the process. When `CLUSTER` is set, `window_bounds` assembles the even and odd heads concurrently in two threads. A development version of `psi_floor` (§4.6) used the context manager `ctx.workprec`; its save and restore of the global value, interleaved across the two threads, left the coordinator at the lower precision. The released `psi_floor` does not change the precision, and `_head` asserts that `ctx.prec` is unchanged after it. The remaining `ctx.workprec` uses run in worker processes (`_partial`) or in the main thread after both head threads have joined (`lowest`, `rayleigh_upper`). Ball arithmetic at a lower precision still returns valid, only wider, enclosures, so this was a robustness issue, not a rigor issue.
2. **Out of memory at $L = 1.06$.** At $N = 2962$ each worker process holds about 1.75 GB of ball sums (measured) plus a shared-memory file. With one worker per two cores and both sectors running at once, the cluster nodes ran out of memory and the run stopped with `BrokenProcessPool` (`data/logs/run_extend.log`). No bounds were written. The fix in `node_partial.py` caps the workers per sector at $0.4 \cdot \mathrm{MemAvailable}/(240 N^2)$ (override: `CERT_NODE_PROCS`) and streams the pairwise merges line by line. The fix changes how panels are split among processes and in which order the ball sums are added, not what is summed. The window $1.03$ was certified before the fix, $1.06$ and $1.09$ after it. The SHA-256 check covers `node_partial.py`, so in each run the hosts used the coordinator's version.
3. **A decimal constant.** The floor $\beta = 3/16$ of the sharp path is an exact dyadic number. A first version used `arb("0.2")`, whose parsing radius entered every head entry. This affected sharpness only.
4. **Memory past $L = 1.09$: triangular assembly.** At $N \approx 4300$–$4600$ ($L = 1.12$–$1.19$) the full $N \times N$ ball products of item 2 no longer fit. Before the windows $1.12$–$1.19$ were run (`certify.py` edited 2026-10-25 15:08, run restarted 15:12), `_partial` was changed to accumulate only the **upper triangle** of each head block, in row blocks of 256 rows (block $b$ = rows $[256b, 256b+256)$ × columns $[256b, N)$), and `node_partial.py` now streams these lines to its shared-memory file one at a time. Each entry is the same sum over the same nodes as before; only about $N^2/2$ entries are formed, with half the multiplications. The worker cap now budgets about $100 N^2$ bytes per worker (the upper-triangle sums, measured $0.87$ GB at $N = 4305$, plus the ~$0.75$ GB shared-memory file) instead of $240 N^2$, and the merge stage uses a fresh process pool so that worker memory is returned. **Equivalence check** (one-off, not in the repository): the old and new `_partial` were run on $(L = 0.3$, even, panels 0–2$)$, $(L = 0.8$, odd, panels 5–6$)$ and $(L = 1.0$, even, panel 40, $N = 1745)$; all $4{,}851$, $11{,}476$ and $1{,}523{,}385$ entries came out as byte-identical strings. The windows $L \le 1.09$ in the certificate of record were certified before this change. A new helper `prefetch.py` assembles head matrices into `cache/` ahead of the certifying run; it was used to pipeline the runs (§4.6).

---

## 1. Introduction

### 1.1 The Weil form

For a real function $f \in L^2(\mathbb{R})$ with support in $[-L, L]$, put $\tilde f(x) = f(-x)$, $g = f * \tilde f$ and $F(z) = \hat f(z) = \int f(x) e^{izx}\,dx$. Weil's explicit-formula functional, evaluated on $g$, gives the quadratic form

$$
Q(f) = W(f * \tilde f), \qquad
W(g) = \hat g(\tfrac i2) + \hat g(-\tfrac i2) - g(0)\log\pi + \frac{1}{2\pi}\int_{\mathbb{R}} \hat g(r)\,\mathrm{Re}\,\psi\!\left(\tfrac14 + \tfrac{ir}{2}\right) dr - \sum_{n \ge 2} \frac{\Lambda(n)}{\sqrt n}\,\big[g(\log n) + g(-\log n)\big].
$$

Weil's criterion says that RH is equivalent to $Q(f) \ge 0$ for all such $f$ and all $L$. Since $g$ is supported in $[-2L, 2L]$, only prime powers with $\log n < 2L$ contribute. For $L \le 1.19$ these are $n = 2, 3, 4, 5, 7, 8, 9$. Each $n$ enters at $L = \tfrac12\log n$:

| $n$ | 2 | 3 | 4 | 5 | 7 | 8 | 9 | 11 |
|---|---|---|---|---|---|---|---|---|
| enters at $L = \tfrac12 \log n$ | 0.3466 | 0.5493 | 0.6931 | 0.8047 | 0.9730 | 1.0397 | 1.0986 | 1.1989 (outside range) |

### 1.2 The simple-even hypothesis

Connes and van Suijlekom [CvS] and Connes–Consani–Moscovici [CCM, §8] consider the operator $A_L$ attached to this form on $L^2(-L, L)$. They show that if its lowest eigenvalue is **simple** with an **even** eigenfunction $\xi_L$, then $\hat\xi_L$ has only real zeros. The simple-even property is assumed there and is listed as a missing step. Even for all $L$, this would not by itself prove RH: the programme also needs the zeros of $\hat\xi_L$ to converge to those of $\zeta$ as $L\to\infty$. This note addresses neither the all-$L$ statement nor the limit; it verifies the hypothesis on $(0, 1.19]$ only.

- Elder [El] (Zenodo deposit; unrefereed) also treats it as open ("Step A").
- Numerically it is consistently observed, for example in [Zhu, Table 1] and in Eberwein [Eb] (Zenodo deposit; unrefereed).

### 1.3 Prior work

The prior-art survey is recorded in our internal survey `NOVELTY.md` (not released). The convention there matches ours: $L$ is the half-width, and $p^k$ enters when $\log p^k < 2L$. For each source we give its status in brackets, using the labels of §9.

- **Suzuki** [Su, Thm 1.4] [arXiv preprint] proves, for the full form, that the lowest eigenvalue is positive and simple with an even eigenfunction for $L \in (0, a_0)$. The proof uses Kato perturbation theory, and $a_0$ is not explicit.
- **mathia PL-294** [PL294] [GitHub notes; unrefereed] shows that the full form satisfies the Beurling–Deny criterion if and only if $L \le a_{BD} = \tfrac12\log\rho_{\mathrm{plastic}} \approx 0.1406$. In that case the semigroup is positivity improving, so the ground state is simple and of fixed sign, hence even. The notes leave evenness implicit.
- **Andrade** [An1] [Zenodo deposit; unrefereed] proves a Perron-type theorem for the *pole-removed* form for all windows. A remark there gives the same threshold for the full kernel ($2a \le t^* = 0.28120$). In that remark the full-form statement is reduced to an equivalence condition that is checked numerically.
- **Zhu** [Zhu] [arXiv preprint] certifies simple-even at the single window $L = 0.8$, in Arb [Zhu, Thm 6.2]. In Zhu's indexing, $\lambda_1^{\mathrm{even}}, \lambda_2^{\mathrm{even}}, \lambda_1^{\mathrm{odd}}$ are our $\lambda_0, \lambda_1^{\mathrm{even}}, \lambda_0^{\mathrm{odd}}$. The theorem states
  $$8.9\cdot10^{-18} \le \lambda_0 \le 2.523\cdot10^{-16}, \qquad \lambda_1^{\mathrm{even}} \ge 2.085\cdot10^{-12}, \qquad 8.206\cdot10^{-15} \le \lambda_0^{\mathrm{odd}} \le 2.347\cdot10^{-14}.$$
  - **Second even bound.** It is rigorous: it comes from a verified Cholesky on a codimension-one restriction at $s_2 = 2.0857\cdot10^{-12}$ (Courant–Fischer), lifted to the full form by the block estimate of Thm 1.1 [Zhu, §6.2].
  - **Positivity.** Zhu certifies $Q \ge 8.9\cdot10^{-18}\|f\|^2$ on even $f$ at $L = 0.8$ [Zhu, Thm 1.2], and for all $f$ [Zhu, Cor 6.3].
  - **Certified upper bounds.** Zhu certifies unconditional variational upper bounds on $\lambda_0$ at $L = 0.5, 0.6, \dots, 1.2, 1.4, \dots, 2.0$ [Zhu, §1.3, §11.4, Table 3], for example $\lambda_0(0.8) \le 2.27\cdot10^{-17}$ and $\lambda_0(1.0) \le 8.08\cdot10^{-30}$. Combined with Thm 1.2 this gives a certified two-sided enclosure $8.9\cdot10^{-18} \le \lambda_0(0.8) \le 2.27\cdot10^{-17}$ [Zhu, §1.4, §8].
  - **The reduction.** Zhu's frequency-cutoff ("one-stroke") reduction with its head/tail block estimate [Zhu, Thm 1.1, proved in §4] is the basis of our Lemmas (R) and (S). The envelope bound behind it is [Zhu, Lemma 3.1], our Lemma (E).
  - Zhu's Table 1 is numerical.
  - Zhu §7 retracts an earlier certified-positivity claim at $L = 1.19$, which is our last window (see §7, item 8).
- **Positivity certificates.** Positivity, or lower bounds on $\lambda_0$, on parts of $(0, 0.973]$ have been certified or claimed by the sources below. None of them certifies simplicity, the parity of the ground state, or an upper bound on $\lambda_0$.
  - jose-paul-b [JPB] [GitHub repository; "AI-audited" per `NOVELTY.md`; unrefereed]: all $L \le (\log 7)/2$.
  - Eureka [Eu, §5.3] [arXiv preprint]: a self-described interval-Cholesky positivity "candidate" on $(0, 0.345]$, not replayed by us.
  - Liu [Liu] [alphaXiv preprint; not located or read]: $L = 1$ and $17/16$.
  - Eriksson/AIRR, via [PL282] [GitHub notes; unrefereed].
  - Elston [Els] [GitHub repository; unrefereed] and Elder [El] [Zenodo deposit; unrefereed]: per-sector lower bounds at $L = 0.8$ (a Zhu re-run), with even-sector upper bounds numerical only.
- **Kramarenko-Byrd** [KB] [Zenodo deposit 23007100, v1 of 2026-09-27; unrefereed] claims a rigorous proof that the ground state is simple and even, and that its Fourier transform has only real zeros, for **all** $0 < L \le 4/5$, in our half-width convention. Its value $\lambda_0(0.8) > 10^{-17}$ is consistent with ours. We have not checked its proof. Our earlier prior-art searches missed it. If it is correct, our result on $(0, 0.8]$ is an independent confirmation, not a new result.
- **grandchallenge / MATHSOLVE** [MS] [unrefereed; not read in full] reportedly gives an analytic argument for windows $L \le 9/50 = 0.18$. We have not verified its content or its window convention.
- **Positivity past $L = 1$.** We found only positivity (or per-sector lower-bound) results past $L = 1$, none of them on simplicity or parity of the ground state: Liu [Liu] ($L = 1$ and $17/16$; not read); Mori [Mo] (Zenodo 23281236 and related deposits; up to $L = (\log 15)/2$); Senke [Se] ($L = 77/64$; per-sector exact Gram certificates; not read in full); Elston [Els] (per-sector lower bounds at $L = 1.19$, $1.3$, $1.495$).
- **Andrade** [An2] [Zenodo deposit; unrefereed] certifies simple-even for three *finite* Galerkin matrices $QW^N$ (with pole) at $L \approx 1.28, 1.68, 1.99$. These are not statements about the continuum form. Andrews Jr. [AnJ] does the same for finite matrices at $L \approx 1.282$.
- **Perturbative simplicity near $L = 0.8$** [PL287] [GitHub notes; unrefereed] holds on a neighbourhood whose radius is not explicit.

### 1.4 What is new

> For $L \le 0.1406$ the ground state is simple and even by the Beurling–Deny argument (explicit threshold in the unrefereed [PL294], cf. the remark in the unrefereed [An1], which builds on [Su §5]; a non-explicit small-$L$ range in [Su, Thm 1.4]; an analytic argument up to $L \le 0.18$ is reported in [MS]). Zhu certified simple-even at the single window $L = 0.8$ [Zhu, Thm 6.2], and Kramarenko-Byrd claims it for every $L \le 0.8$ [KB]. **We give a rigorous, computer-assisted (Arb) proof that the ground state of the full localized Weil form is simple with even eigenfunction for every $L \in (0, 1.19]$. The part we believe to be new is $(0.8, 1.19]$**, where the prime powers $5$, $7$, $8$ and $9$ enter (at $L = 0.8047$, $0.9730$, $1.0397$, $1.0986$). On $(0, 0.8]$ our certificate is an independent computer-assisted confirmation of the range claimed in [KB] (independent of [KB]; not methodologically independent of [Zhu], whose reduction it uses), and at $L = 0.8$ a re-certification of [Zhu, Thm 6.2] with independent code. We have not compared methods with [KB].

**Credit to Zhu.** The whole reduction strategy is Zhu's "one-stroke reduction" [Zhu, Thm 1.1, proved in §4]:

- **Frequency cutoff.** Cut the frequency integral at $T^\#$ and replace $\Psi_L$ by the constant $\beta^*$ above it (our Lemma (R)). The envelope behind $\beta^*$ is [Zhu, Lemma 3.1] (our Lemma (E)).
- **Head and tail.** Bound the infinite Legendre tail by the head/tail block estimate $\min(\lambda_0, \beta^* - \varepsilon_D) - \varepsilon_B$ (our Lemma (S)). Zhu applies the same block estimate to the second even eigenvalue [Zhu, §6.2], and we use it in that form for every $k$.

The scale $2\pi e^{A_L}$ in our cutoff rule is Zhu's $T_1$. The per-window certificate follows the pattern of [Zhu, Thm 6.2]: certified lower bounds on $\lambda_1^{\mathrm{even}}$ and $\lambda_0^{\mathrm{odd}}$ compared against an upper bound on $\lambda_0$.

**Our own contributions:**

- self-contained proofs of the standard estimates (B), (M), (J), (D) and a Binet-based proof of Zhu's envelope (E) (these are classical facts, not new results), and our own re-derivation and implementation of (R) and (S);
- a certified covering of a continuum of windows by monotonicity in $L$;
- a Lipschitz argument at small $L$;
- the numerics (hybrid Bessel evaluation, and the sharp-cutoff path of §4.6) that reach $L = 1.19$, past the entry of the prime powers 7, 8 and 9.

**By-products:**

- certified two-sided brackets for $\lambda_0(L)$ at all 31 grid windows. Certified upper bounds on $\lambda_0$ are not new for $L \ge 0.5$ [Zhu, Table 3]. What is new is only the brackets on this particular grid, including the windows below $0.5$. Where the grids overlap ($0.5, 0.6, 0.8, 0.9, 1.0$), the $U$ of the corrected original certificate (Table 1) is smaller than Zhu's Table 3 values. In the certificate of record (Table 2) this holds at $0.5, 0.6, 0.9$ and $1.0$ but not at $0.8$ ($3.681\cdot10^{-17}$ against $2.27\cdot10^{-17}$);
- a reproduction of [Zhu, Thm 1.2]. It uses independent code but not an independent method, since it has the same normalization and cutoff strategy.

**Not claimed as new:**

- the statement on $(0, 0.8]$, claimed by [KB]; on $(0, 0.1406]$ by [PL294] (cf. [An1], [Su]) and, reportedly, on $(0, 0.18]$ by [MS];
- the window $L = 0.8$ by itself;
- simplicity on a non-explicit neighbourhood of $L = 0.8$ [PL287];
- positivity, i.e. $\lambda_0 > 0$, on this range;
- certified upper bounds on $\lambda_0$ for $L \in [0.5, 2]$ [Zhu, Table 3], and the two-sided enclosure at $L = 0.8$ [Zhu, §1.4];
- the frequency-cutoff reduction and its head/tail block estimate [Zhu, Thm 1.1];
- the expected ordering of the spectrum, which was anticipated numerically.

---

## 2. Setup and normalization

### 2.1 Frequency side

For real $f$ supported in $[-L, L]$ we have $\hat g(t) = F(t)F(-t)$, which equals $|F(t)|^2$ for real $t$. By Plancherel, $g(0) = \|f\|^2 = \frac1\pi \int_0^\infty |F|^2$, so

$$
Q(f) = [\mathrm{pole}] + \frac1\pi \int_0^\infty \Psi_L(t)\,|F(t)|^2\,dt, \qquad
\Psi_L(t) = \mathrm{Re}\,\psi\!\left(\tfrac14 + \tfrac{it}2\right) - \log\pi - \sum_{\log n < 2L} \frac{2\Lambda(n)}{\sqrt n}\cos(t\log n),
$$

with $[\mathrm{pole}] = \hat g(\tfrac i2) + \hat g(-\tfrac i2) = 2F(\tfrac i2)F(-\tfrac i2)$. Here $F(\tfrac i2) = \int f e^{-x/2}$ is real.

- For even $f$: $[\mathrm{pole}] = +2F(\tfrac i2)^2$.
- For odd $f$: $[\mathrm{pole}] = -2F(\tfrac i2)^2$.

This normalization is taken from Zhu [Zhu, §6]. It was checked numerically against the zero side of the explicit formula (first 400 zeros) in `weil_z.py`.

### 2.2 Position side

Let $\gamma$ be Euler's constant. For $g = f * \tilde f$ (even):

$$
Q(f) = 4\int_0^{2L} g(x)\cosh\tfrac x2\,dx - g(0)\big(\log\pi + \gamma + \log(1 - e^{-4L})\big) + \int_0^{2L} \frac{2\,[e^{-2x}g(0) - e^{-x/2}g(x)]}{1 - e^{-2x}}\,dx - 2\sum_{\log n < 2L} \frac{\Lambda(n)}{\sqrt n}\,g(\log n). \tag{P}
$$

*Derivation.* Take $f$ smooth, so that $g$ is smooth.

1. Insert $\psi(s) = -\gamma + \int_0^\infty \frac{e^{-t} - e^{-st}}{1 - e^{-t}}\,dt$ and apply Fourier inversion. Fubini applies because $g$ is smooth.
2. For $x > 2L$ only the $g(0)$ part of the integrand survives, and it integrates in closed form. This produces the $\log(1 - e^{-4L})$ term.
3. The pole term gives $4\int_0^{2L} g\cosh(x/2)$.

**Extension to the trial functions.** We need (P) for $f$ a polynomial on $[-L, L]$. Then $g$ is a polynomial on $[0, 2L]$ and Lipschitz on $\mathbb{R}$.

- Mollify, $f_\varepsilon = f * \phi_\varepsilon$, and apply (P) on the slightly larger window $L + \varepsilon$.
- The formula is consistent across windows: no $\log n$ lies at $2L$, and the $\log(1 - e^{-4L})$ term absorbs the change of upper limit.

Now let $\varepsilon \to 0$:

- The frequency side converges by dominated convergence, since $\Psi_L|F|^2 = O(t^{-2}\log t)$.
- On the position side, $g_\varepsilon \to g$ uniformly with a uniform Lipschitz constant. This dominates the integrand near $x = 0$, where it behaves like $(g(0) - g(x))/(2x)$.

Formula (P) was checked against `weil_z.py`.

### 2.3 The operator, parity, discrete spectrum

Let $A_L$ be the self-adjoint operator associated with the closure of $Q$ on $C_c^\infty(-L, L)$ in $H_L = L^2(-L, L)$.

**The spectrum is discrete.**

- By Lemma (E) below, $\Psi_L(t) \ge \log(t/2\pi) - 1/t - A_L \to \infty$.
- The pole term is bounded by $C\|f\|^2$.
- Hence the form norm controls $\int\log(2+t)|F|^2$.

For functions supported in $[-L, L]$, the set where this quantity and $\|f\|$ are bounded is compact in $L^2$:

- high frequencies are uniformly small, because the weight tends to $\infty$;
- on compact $t$-sets the $F$ are equicontinuous entire functions (Paley–Wiener).

So $A_L$ has compact resolvent. This agrees with [CCM, Thm 3.6].

Write $\lambda_0(L) \le \lambda_1(L) \le \dots$ for the min-max values, which here are eigenvalues. The form has a real kernel, so real and imaginary parts decouple, and we work with real $f$ throughout.

**Lemma 2.1 (parity splitting).** Let $(Rf)(x) = f(-x)$. Then $Q(Rf) = Q(f)$, and $Q(f_e + f_o) = Q(f_e) + Q(f_o)$ for even $f_e$ and odd $f_o$. Hence $A_L = A_L^{\mathrm{even}} \oplus A_L^{\mathrm{odd}}$, and the spectrum of $A_L$ is the union of the two sector spectra.

*Proof.* Replacing $f$ by $Rf$ replaces $F(t)$ by $F(-t)$. This leaves $|F(t)|$ ($t$ real) and the product $F(\tfrac i2)F(-\tfrac i2)$ unchanged. Let $B$ be the bilinear form of $Q$. Then $B(f_e, f_o) = B(Rf_e, Rf_o) = B(f_e, -f_o) = -B(f_e, f_o)$, so $B(f_e, f_o) = 0$. ∎

**Corollary 2.2.** The ground state of $A_L$ is simple and even if and only if
$$\lambda_0^{\mathrm{even}}(L) < \min\big(\lambda_1^{\mathrm{even}}(L),\ \lambda_0^{\mathrm{odd}}(L)\big).$$

The simple-even property is invariant under $Q \mapsto aQ + b\|\cdot\|^2$ with $a > 0$. So conventions that differ from ours by a positive scale or a constant shift do not affect the theorem. Examples are where $\log\pi$ is placed and a global factor of $2\pi$.

### 2.4 Monotonicity in $L$

**Lemma 2.3.** For $0 < L < L'$ and every $k$:
$$\lambda_k^{\mathrm{even}}(L') \le \lambda_k^{\mathrm{even}}(L), \qquad \lambda_k^{\mathrm{odd}}(L') \le \lambda_k^{\mathrm{odd}}(L).$$

*Proof.* Extension by zero, $\iota : L^2(-L, L) \to L^2(-L', L')$, is an isometry. It preserves parity and maps the form domain into the form domain. The form itself does not change: $Q_{L'}(\iota f) = Q_L(f)$.

- The pole and archimedean terms depend only on $F = \hat f$.
- A prime power with $2L \le \log n < 2L'$ contributes $g(\log n)$, which is $0$: $g = f * \tilde f$ is continuous and supported in $[-2L, 2L]$, and $g(\pm 2L) = 0$.

By min-max, $\lambda_k = \inf_{\dim V = k+1}\sup_{f \in V} Q(f)/\|f\|^2$. Taking the infimum over a larger family of subspaces (those of $H_{L'}$, which include the images $\iota V$) can only decrease it. ∎

**Consequence.** On a step $[L_i, L_{i+1}]$:
$$\lambda_0(L) \le \lambda_0(L_i) \le U(L_i), \qquad \lambda_1^{\mathrm{even}}(L) \ge \lambda_1^{\mathrm{even}}(L_{i+1}), \qquad \lambda_0^{\mathrm{odd}}(L) \ge \lambda_0^{\mathrm{odd}}(L_{i+1}).$$

---

## 3. Main theorem

**Theorem.** For every $L \in (0, 1.19]$, the lowest eigenvalue of Weil's quadratic form for $\zeta$ on $L^2(-L, L)$ is simple and its eigenfunction is even. Here the form is taken on the log scale, with pole terms and all prime powers with $\log n < 2L$ (§2). Quantitatively:

- On $(0, 0.05]$: $\lambda_1^{\mathrm{even}} - \lambda_0 \ge 0.4172$ and $\lambda_0^{\mathrm{odd}} - \lambda_0 \ge 0.6624$.
- On each grid step $[L_i, L_{i+1}] \subset [0.05, 1.19]$: $\min(\lambda_1^{\mathrm{even}}, \lambda_0^{\mathrm{odd}})/\lambda_0 \ge \rho_i$, with $\rho_i \ge 1.288$ for every step (Table 2).

These numbers are from the certificate of record (Table 2), outward-rounded lower ends. (The corrected original certificate gives $0.4173$, $0.6625$ and $1.288$ on $(0, 1]$; Table 1.)

The proof is computer-assisted and uses Arb ball arithmetic in two places:

- python-flint 0.9.0, which bundles FLINT 3.6.0;
- the C tools `vchol` and `veig`, linked against the system FLINT 3.4 (FLINT on the cluster hosts not recorded).

**Reproducibility status of the certificate.** The certificate of record is `data/theorem_bounds.json`: the clean rerun for $L \le 1$ (one code revision, corrected Gauss–Legendre radius folded into every entry, remote code checked by SHA-256) and the extension run for $L = 1.03, 1.06, 1.09$ and its continuation for $L = 1.12, 1.15, 1.18, 1.19$ (triangular head assembly, Erratum, further issues, item 4). The windows $L > 1$ were certified with the sharp-cutoff path (§4.6); its validation at $L = 1$ is in §4.6. The original runs and their correction (Erratum) are history; the theorem does not depend on them.

The theorem is also subject to the normalization caveat: the form is matched to CvS/CCM only up to positive affine rescaling (§7.3).

---

## 4. Proof

### 4.1 Piece I: small windows $L \in (0, 0.05]$ (`lipschitz.py`)

Here $L \le 0.05 < \tfrac12\log 2$, so no prime powers occur.

**Rescaling.** Put $f(x) = L^{-1/2} f_1(x/L)$ with $\|f_1\|_{L^2(-1,1)} = 1$, and $F_1 = \hat f_1$. Substituting $t = \tau/L$ gives

$$
Q(f) = \log\tfrac1L + R_L(f_1), \qquad
R_L(f_1) = P_L(f_1) + \frac1\pi \int_0^\infty m_L(\tau)\,|F_1(\tau)|^2\,d\tau, \qquad
m_L(\tau) = \mathrm{Re}\,\psi\!\left(\tfrac14 + \tfrac{i\tau}{2L}\right) - \log\pi + \log L,
$$

with

$$
P_L(f_1) = +2L\Big(\int_{-1}^1 f_1(s)e^{-Ls/2}ds\Big)^2 \ \ (f_1 \text{ even}), \qquad
P_L(f_1) = -2L\Big(\int_{-1}^1 f_1(s)\sinh(Ls/2)ds\Big)^2 \ \ (f_1 \text{ odd}).
$$

The shift $\log(1/L)$ is the same for all $f_1$. So the **gaps** of $Q_L$ equal the gaps of $R_L$, sector by sector, and $R_L$ acts on the fixed space $L^2(-1, 1)$.

**Derivative of the multiplier.** Put $u = \tau/(2L)$. Then
$$\frac{\partial m_L}{\partial L} = \frac1L\,h(u), \qquad h(u) = 1 - \mathrm{Re}\big(iu\,\psi'(\tfrac14 + iu)\big).$$

**Lemma 4.1 (bound on $h$).** For $u > 0$, $|h(u)| \le 0.8125/u^2$.

*Proof.* Expand $\psi'(w) = \sum_{n \ge 0}(n + w)^{-2}$. This gives $h(u) = 1 - \sum_{n \ge 0}\varphi(n + \tfrac14)$, with $\varphi(a) = 2au^2/(a^2 + u^2)^2$ and $\int_0^\infty \varphi = 1$. The sum is the midpoint rule for $\int_{-1/4}^\infty \varphi$ on cells of width 1.

- The midpoint Peano kernel is bounded by $1/8$, and $\int_{\mathbb{R}}|\varphi''| = 6/u^2$. So the quadrature error is at most $\tfrac34 u^{-2}$.
- The piece of the integral below $0$ satisfies $|\int_{-1/4}^0 \varphi| \le \int_{-1/4}^0 2|a|u^{-2}\,da = \tfrac1{16}u^{-2}$.

Adding the two gives $(\tfrac34 + \tfrac1{16})u^{-2}$. ∎

(Numerically $u^2 h(u) \to -1/48$, so the bound has a large margin.)

**Bound on $\int|h|$.**

- On $[0, 40]$ we bound $\int|h|$ by an upper Riemann sum: cells of width $1/200$, with $\sup|h|$ on each cell enclosed by evaluating $h$ on the ball covering the cell.
- On $[40, \infty)$, Lemma 4.1 contributes $0.8125/40$.

The resulting enclosure depends on the working precision. All entries below are upper ends of the Arb enclosures, rounded up.

| Precision | How it is run | $\int\|h\| \le$ | $C_{\mathrm{even}} \le$ | $C_{\mathrm{odd}} \le$ |
|---|---|---|---|---|
| 768 bits | `theorem.py` and `gl_correction.py`, which import `certify` and so run at `CERT_PREC` | 0.593979 | 4.876362 | 0.761666 |
| 53 bits | standalone `lipschitz.py` | 0.592773 | 4.87484 | 0.760130 |

The corrected certificate recomputed the 768-bit constants (`run_gl_correction.log`: $C_{\mathrm{even}} \in [4.87636 \pm 1.7\cdot10^{-6}]$, $C_{\mathrm{odd}} \in [0.761665 \pm 4.2\cdot10^{-7}]$), and we quote those.

Both sets of constants are upper bounds under Arb's containment guarantee. It is unexplained why the higher-precision enclosure is the larger one; see §7.5.

**Lipschitz constants.** Since $|F_1|^2 \le (\int|f_1|)^2 \le 2$, the multiplier part of $\frac{d}{dL}R_L(f_1)$ is bounded by
$$\frac1\pi \int_0^\infty \frac1L\,|h(\tfrac\tau{2L})|\cdot 2\,d\tau = \frac4\pi\int_0^\infty |h|.$$

For the even pole term, write $J = \int f_1 e^{-Ls/2}$. Then $\frac{d}{dL}(2LJ^2) = 2J^2 + 4LJJ'$, and Cauchy–Schwarz gives
$$J^2 \le \frac{2\sinh L}{L}, \qquad |J'| \le \tfrac12\Big(\tfrac23 e^L\Big)^{1/2}.$$

For the odd pole term, write $K = \int f_1\sinh(Ls/2)$. Then
$$K^2 \le 2\sinh^2(L/2), \qquad |K'| \le \tfrac12\cosh(L/2)\sqrt{2/3}.$$

All of these bounds increase in $L$. Evaluating them at $L_0 = 0.05$ gives, uniformly for $L \le 0.05$,
$$C_{\mathrm{even}} = \tfrac4\pi\textstyle\int|h| + 2J^2 + 4L_0|J||J'| \le 4.876362, \qquad C_{\mathrm{odd}} \le 0.761666.$$

The forms $R_L$ share the core $C_c^\infty(-1, 1)$ and satisfy $|R_L(f_1) - R_{L'}(f_1)| \le C_{\mathrm{sector}}|L - L'|\,\|f_1\|^2$. By min-max, every sector eigenvalue of $R_L$ is therefore $C_{\mathrm{sector}}$-Lipschitz in $L$.

**Conclusion of Piece I.** Let $U(L_0)$ be the certified upper bound on $\lambda_0(L_0)$, and let $\ell_1, \ell_{\mathrm{odd}}$ be the certified lower bounds at $L_0 = 0.05$ (Table 1, corrected values). For $L \le L_0$:

$$
\lambda_1^{\mathrm{even}}(L) - \lambda_0(L) \ge \ell_1 - U - 2C_{\mathrm{even}}L_0 \ge 0.4173, \qquad
\lambda_0^{\mathrm{odd}}(L) - \lambda_0(L) \ge \ell_{\mathrm{odd}} - U - (C_{\mathrm{even}} + C_{\mathrm{odd}})L_0 \ge 0.6625.
$$

The exact lower ends are $0.417374$ and $0.662537$, recomputed from `theorem_bounds_corrected.json` with the constants above. Both are positive, so Piece I holds by Corollary 2.2. With the certificate of record (Table 2) the same computation gives $\ge 0.4172$ and $\ge 0.6624$ (§3); this is what `theorem.py --check` verifies. ∎

### 4.2 Piece II: the grid (`theorem.py`, re-checked by `gl_correction.py`)

The grid is
$$0.05,\ 0.1,\ 0.2,\ 0.3,\ 0.35,\ 0.4,\ 0.45,\ 0.5,\ 0.55,\ 0.6,\ 0.64,\ 0.68,\ 0.72,\ 0.76,\ 0.8,\ 0.84,\ 0.87,\ 0.9,\ 0.93,\ 0.96,\ 0.97,\ 0.98,\ 0.99,\ 1.0,\ 1.03,\ 1.06,\ 1.09,\ 1.12,\ 1.15,\ 1.18,\ 1.19.$$

For each step $[L_i, L_{i+1}]$ we check
$$\min\big(\ell_1(L_{i+1}),\ \ell_{\mathrm{odd}}(L_{i+1})\big) > U(L_i).$$

By Lemma 2.3 this gives $\lambda_0(L) < \min(\lambda_1^{\mathrm{even}}(L), \lambda_0^{\mathrm{odd}}(L))$ on the whole step. Prime-power entries inside a step are harmless, because Lemma 2.3 holds across them.

All 23 steps on $[0.05, 1]$ pass with the corrected bounds (ratios in Table 1), and all 30 steps on $[0.05, 1.19]$ pass with the certificate of record (ratios in Table 2).

- The smallest margins are at the low end: $1.288$ on $[0.05, 0.1]$ and $1.320$ on $[0.1, 0.2]$.
- From $[0.96, 0.97]$ to $[0.99, 1]$ the ratios are $332.7$–$380.0$ (Table 1; in Table 2: $255.3$–$354.9$, and past $L = 1$: $10.00$, $6.548$, $3.774$, $3.933$, $2.857$, $2.148$, $264.7$).

**How the minimum is taken.** Both `theorem.py` and `gl_correction.py` use Python's built-in `min` on Arb balls. That is not a rigorous minimum when the balls overlap. Here it is harmless: in `theorem_bounds_corrected.json`, $\ell_1$ and $\ell_{\mathrm{odd}}$ are disjoint at every window.

- **Relative separation.** The larger of the two exceeds the smaller by a factor of at least $1.021$ (at $L = 0.05$), and by at least $1.09$ from $L = 0.1$ on.
- **Radii.** Every radius is below $5\cdot10^{-20}$ in absolute terms, and at every window the absolute separation exceeds the radii by a factor of more than $10^{17}$.

The ratios in Table 1 were recomputed with explicit endpoint comparisons.

**Remark (dyadic windows).** The code certifies each window at the dyadic number $\tilde L_i$ = `arb(L_i).upper()`, which satisfies $L_i \le \tilde L_i < L_i + 2^{-p}$ at working precision $p \ge 768$.

- All four bounds of a window, including the trial function in §4.4, refer to $[-\tilde L_i, \tilde L_i]$.
- The covering is therefore $(0, \tilde L_1] \cup \bigcup_i [\tilde L_i, \tilde L_{i+1}]$. Since $\tilde L_{31} = 1.19$ up to $2^{-p}$ and $\tilde L_{31} \ge 1.19$, this covers $(0, 1.19]$.
- In Piece I, $L_0$ enters as the ball `arb("0.05")`, which contains $\tilde L_1$, and all bounds there are monotone in $L$.

### 4.3 Per-window lower bounds (`certify.py`)

Fix $L$ and a parity. We use the orthonormal Legendre basis $e_n(x) = c_n P_n(x/L)$ with $c_n = \sqrt{(2n+1)/(2L)}$, taking orders $n = 2m$ (even sector) or $n = 2m+1$ (odd sector). Up to a unit factor common to the sector,
$$\hat e_n(t) = c_n\,2L\,(-1)^m j_n(Lt),$$
where $j_n$ is the spherical Bessel function. The pole vector is
$$p_m = \sqrt2\,c_n\,2L\,i_n(L/2),$$
so that $[\mathrm{pole}] = \pm(p\cdot a)^2$ for $f = \sum a_m e_{n(m)}$.

**Lemma (B), Binet bound.** For $\mathrm{Re}\,w > 0$,
$$\big|\psi(w) - \log w + \tfrac1{2w}\big| \le \frac1{12|w|\,\mathrm{Re}\,w}.$$

*Proof.* Binet's formula reads $\psi(w) = \log w - \tfrac1{2w} - \int_0^\infty b(x)e^{-wx}\,dx$ with $b(x) = \tfrac1{e^x - 1} - \tfrac1x + \tfrac12$. We have $b(0) = 0$ and $0 \le b' \le \tfrac1{12}$.

For the upper bound on $b'$, write $y = x/2$. Then $b'(x) = \tfrac14(y^{-2} - \sinh^{-2}y)$.

- For $y < \sqrt3$, the inequality $b' \le 1/12$ reduces to $\sinh^2 y \le 3y^2/(3 - y^2)$. This holds coefficient by coefficient: $2^{2k-1}/(2k)! \le 3^{1-k}$, with equality at $k = 1, 2$, and the term ratio is $4/((2k+1)(2k+2)) \le 1/3$.
- For $y \ge \sqrt3$ it is trivial.

Integrating by parts, $\int b\,e^{-wx} = w^{-1}\int b'e^{-wx}$, whose modulus is at most $1/(12|w|\,\mathrm{Re}\,w)$. ∎

**Lemma (E), envelope** (the envelope of [Zhu, Lemma 3.1]; we prove it via (B)). For $w = \tfrac14 + \tfrac{it}2$ and $t \ge 3/2$,
$$\mathrm{Re}\,\psi(w) \ge \log\tfrac t2 - \tfrac2{3t} - \tfrac1{2t^2} \ge \log\tfrac t2 - \tfrac1t.$$

*Proof.* Apply (B) with $|w| \ge t/2$ and $\mathrm{Re}(1/(2w)) \le 1/(2t^2)$. The second inequality is equivalent to $t \ge 3/2$. ∎

Put $A_L = \sum_{\log n < 2L} 2\Lambda(n)/\sqrt n$. The prime comb is $\ge -A_L$, and the envelope increases in $t$. Hence for $t \ge T^\#$,
$$\Psi_L(t) \ge \beta := \log\frac{T^\#}{2\pi} - \frac1{T^\#} - A_L > 0.$$
This is Zhu's $\beta^*$ [Zhu, Thm 1.1].

**Lemma (R), frequency reduction (Zhu's one-stroke reduction [Zhu, Thm 1.1]).** Let
$$\mathcal{R}(f) = [\mathrm{pole}] + \frac1\pi\int_0^{T^\#}(\Psi_L - \beta)|F|^2\,dt + \beta\|f\|^2.$$
Then $Q \ge \mathcal{R}$ as forms, and therefore $\lambda_k(Q) \ge \lambda_k(\mathcal{R})$ in each sector.

*Proof (as in [Zhu, §4]).* Write $\Psi_L = (\Psi_L - \beta) + \beta$ and use $\frac1\pi\int_0^\infty|F|^2 = \|f\|^2$. Drop the nonnegative integral over $[T^\#, \infty)$, then apply min-max. ∎

Zhu states the reduction for the even sector and notes that it holds verbatim in the odd sector [Zhu, §6.1].

**Cutoff rule.** The code sets
$$T^\# = \max\Big(200,\ 5\big\lceil 80/(5L)\big\rceil,\ 5\big\lceil 2\pi e^{A_L + 0.5}/5\big\rceil\Big).$$

- The floor value $200$ is Zhu's choice at $L = 0.8$.
- The second term makes $LT^\# \ge 80$, which keeps $\beta$ well above the eigenvalues at small $L$.
- The third term is Zhu's scale $T_1 = 2\pi e^{A_L}$, inflated so that $\beta \gtrsim 0.5$ at large $L$.

The head has $N = \lceil eLT^\#/4\rceil + 42$ modes per sector. Values per window are listed in Table 1; for example, $T^\# = 3610$ and $N = 2496$ at $L = 1$.

**Lemma (M), monotonicity of $\mathrm{Re}\,\psi$.**
$$\mathrm{Re}\,\psi(s + iy) = \psi(s) + \sum_{n \ge 0}\frac{y^2}{(n+s)\big((n+s)^2 + y^2\big)}$$
is increasing in $|y|$ for $s > 0$.

*Proof.* Take the real part of $\psi(w) = -\gamma + \sum_n\big(\tfrac1{n+1} - \tfrac1{n+w}\big)$. ∎

So on $[0, T^\#]$,
$$|\Psi_L - \beta| \le M_\psi := \max\big(|\mathrm{Re}\,\psi(\tfrac14)|,\ |\mathrm{Re}\,\psi(\tfrac14 + \tfrac{iT^\#}2)|\big) + \log\pi + A_L + \beta.$$

**Lemma (J), Bessel bounds.**
$$|j_n(x)| \le \frac{x^n}{(2n+1)!!}, \qquad 0 \le i_n(x) \le \frac{x^n}{(2n+1)!!}\,e^{x^2/2}.$$

*Proof.* Use the Poisson integral $j_n(x) = \frac{x^n}{2^{n+1}n!}\int_{-1}^1\cos(xs)(1 - s^2)^n\,ds$ together with $\int_{-1}^1(1 - s^2)^n\,ds = 2^{2n+1}(n!)^2/(2n+1)!$. For $i_n$, replace $\cos$ by $\cosh$ and use $\cosh(xs) \le e^{x^2/2}$ on $[-1, 1]$. ∎

**Lemma (G), Gauss–Legendre and Bernstein [Tr, Thm 19.3].** Suppose $g$ is analytic in the Bernstein ellipse $E_\rho$ and $|g| \le M$ there. Then the **$n$-point** Gauss–Legendre rule on $[-1, 1]$ has error at most
$$\tfrac{64}{15}M\rho^{-2n+2}/(\rho^2 - 1).$$
[Tr] states the bound with $\rho^{-2n}$ for the $(n+1)$-point rule. Also (Bernstein), a polynomial of degree $d$ bounded by $M$ on $[-1, 1]$ is bounded by $M\rho^d$ on $E_\rho$.

The current code uses this bound. The original window runs used $\rho^{-2n}$; see Erratum and correction.

**Head matrix and quadrature.** We compute
$$C_{jk} = \frac1\pi\int_0^{T^\#}(\Psi_L - \beta)\,\hat e_j\hat e_k\,dt$$
by 64-point Gauss–Legendre on each unit panel $[k, k+1]$, $k < T^\#$, in Arb.

- **Ellipse.** The ellipse has half-height $B = 2/5$ in $t$, so $\rho = 2B + \sqrt{4B^2 + 1}$.
- **Continuation of the integrand.** For complex $t$, the analytic continuation of $\mathrm{Re}\,\psi(\tfrac14 + \tfrac{it}2)$ is $\tfrac12\big(\psi(w_+) + \psi(w_-)\big)$ with $w_\pm = \tfrac14 \pm \tfrac{it}2$.
  - On the strip both arguments satisfy $\mathrm{Re}\,w_\pm \ge \tfrac14 - \tfrac B2 = 0.05 > 0$ and $|w_\pm| \le |w|_{\max}$.
  - So $\psi(w_\pm)$ is analytic there, and (B) bounds it by $\max(\log|w|_{\max}, -\log\mathrm{Re}\,w) + \pi/2 + 1/(2\mathrm{Re}\,w) + 1/(12(\mathrm{Re}\,w)^2)$.
- **Basis functions.** On the strip, $|\hat e_j| \le \sqrt{2L}\,e^{LB}$.
- **Prime comb.** It is bounded by $\sum 2\Lambda(n)n^{-1/2}\cosh(B\log n)$.

The per-entry error is $T^\#\cdot\frac{32}{15}M\rho^{-126}/(\rho^2 - 1)$, where the factor $\tfrac12$ is the Jacobian of a unit panel. (The original runs used $\rho^{-128}$; the corrected certificate accounts for the difference by the shift $\delta$.) The error is added as a ball to every entry. The head is
$$A = C \pm pp^{T} + \beta I + \text{quadrature ball},$$
with the sign $+$ in the even sector and $-$ in the odd sector.

**Tails.** For tail orders $n \ge n_0 = n(N)$ on $[0, T^\#]$, put $X = LT^\#$. Then
$$|\hat e_n| \le \sigma_n := c_n\,2L\,X^n/(2n+1)!!,$$
and $\sigma_{n+2}/\sigma_n \le r$ with
$$r = \sqrt{\tfrac{2n_0+5}{2n_0+1}}\cdot\frac{X^2}{(2n_0+3)(2n_0+5)} < 1$$
(asserted in the code). Hence $\sum\sigma^2 \le \sigma_{n_0}^2/(1 - r^2)$.

The tail of the pole vector satisfies $\|p_{\mathrm{tail}}\| \le p_{\mathrm{tail}} := \sqrt2\,c_{n_0}2L\,y^{n_0}/(2n_0+1)!!\cdot e^{y^2/2}/\sqrt{1 - y^2}$, with $y = L/2$. Set

$$
e_D = \frac{M_\psi}{\pi}T^\#\sum\sigma^2 + [\text{odd}]\,p_{\mathrm{tail}}^2, \qquad
e_B = \frac{M_\psi}{\pi}T^\#\sqrt{2L}\sqrt N\Big(\sum\sigma^2\Big)^{1/2} + 2\|p_{\mathrm{head}}\|\,p_{\mathrm{tail}}.
$$

Then:

- the tail–tail block satisfies $D \ge (\beta - e_D)I$; in the even sector the tail pole term $+pp^T \ge 0$ is simply dropped;
- the head–tail block satisfies $\|X\| \le e_B$, by Hilbert–Schmidt bounds.

**Lemma (S), block bound.** This is the head/tail block estimate of [Zhu, Thm 1.1], stated for the $k$-th eigenvalue as Zhu uses it for the second even eigenvalue in [Zhu, §6.2]. Suppose
$$\mathcal{R} = \begin{pmatrix} A & X \\ X^{T} & D\end{pmatrix}, \qquad D \ge (\beta - e_D)I, \qquad \|X\| \le e_B.$$
Then $\lambda_k(\mathcal{R}) \ge \min(\lambda_k(A), \beta - e_D) - e_B$.

*Proof (our re-derivation).* The off-diagonal part satisfies $\big\|\begin{pmatrix}0 & X\\ X^T & 0\end{pmatrix}\big\| = \|X\|$, so $\mathcal{R} \ge \mathrm{diag}(A, D) - e_B$. The $k$-th eigenvalue of $\mathrm{diag}(A, D)$ is at least $\min(\lambda_k(A), \lambda_0(D))$. Conclude with Weyl. ∎

**Lemma (D), deflation.** If $A + vv^T - \mu I \ge 0$, then $A$ has at most one eigenvalue below $\mu$, i.e. $\lambda_1(A) \ge \mu$.

*Proof.* Rank-one interlacing: $\lambda_1(A) \ge \lambda_0(A + vv^T)$. ∎

**Verified Cholesky** (`certified_pd_py`, and `vchol.c` for $N > 120$). To show $A - \mu I \ge 0$ for a ball matrix $A$, set $d = \mu/100$.

1. Compute an approximate Cholesky factor $R$ of $\mathrm{mid}(A) - (\mu + d)I$ on midpoints only, rounding every entry to its midpoint. We avoid the ball Cholesky `arb_mat_cho` because it blows up at $n \approx 2500$.
2. Compute $E = A - (\mu + d)I - RR^T$ rigorously in ball arithmetic.
3. If $\|E\|_F < d$, then $A - \mu I = RR^T + E + dI \ge 0$ for every matrix in the ball.

**Assembling the lower bounds.** Estimates come from inverse iteration on midpoints (`lowest`; `veig.c` for $N > 300$). These estimates are used **only** to choose $\mu = \mathrm{est}/1.05$ and the deflation vector $v$; they carry no rigor. With $\lambda_k(A) \ge \mu$ certified as above, Lemma (S) gives:

- $\lambda_0 \ge \min(\mu_0, \beta - e_D) - e_B$, from the even head;
- $\lambda_1^{\mathrm{even}} \ge \min(\mu_1, \beta - e_D) - e_B$, from the deflated head $A + vv^T$ via Lemma (D);
- $\lambda_0^{\mathrm{odd}} \ge \min(\mu^{\mathrm{odd}}, \beta - e_D^{\mathrm{odd}}) - e_B^{\mathrm{odd}}$.

Each of these bounds $\lambda_k(Q)$ from below by Lemma (R). In the corrected certificate each is further lowered by the per-window $\delta$ (Erratum and correction).

The code computes $\min(\mu, \beta - e_D)$ with Python's `min` on balls. That is harmless here because $\beta - e_D \ge 0.5$ well exceeds every $\mu$:

- at small $L$, where $\beta \approx 5.5$;
- at large $L$, where $\mu \ll 1$.

### 4.4 Per-window upper bound on $\lambda_0$ (`rayleigh_upper`)

**Trial function.** In the corrected certificate, every $U$ comes from `gl_correction.py`.

- The trial vector $v$ is the even ground vector of a non-rigorous Galerkin matrix (`convergence.galerkin`) of degree $D_G = \min(2\lceil eL\cdot 2\pi e^{2L}/8\rceil + 60,\ 298)$. On the grid this ranges from $62$ at $L = 0.05$ to $92$ at $L = 1$ (for example 66 at 0.5, 78 at 0.8, 84 at 0.9).
- $v$ has $m^* = D_G/2 \le 46$ even coefficients. The code truncates to the first 150, which is a no-op on this grid.
- The trial function is $f = \sum_{m < m^*} v_m e_{2m}$, an even polynomial in $s = x/L$ of degree $D = 2m^* - 2 \le 90$.
- (The original runs used the head's approximate ground vector with $m^* = \min(N, 150)$. Those $U$ values are superseded.)

All computations below are at 1536 bits.

1. **The autocorrelation.** $g(Ls) = L\int_{-1}^{1-s} f(y)f(y+s)\,dy$ is computed exactly by $(D+1)$-point Gauss–Legendre, which is exact for the degree-$2D$ integrand.
2. **The integral in (P).** It is computed by 450-point Gauss–Legendre in $s \in [0, 2]$.
3. **Quadrature error.** We bound it by (G) with the correct exponent $\rho^{-898}$, on the ellipse with foci $0, 2$ and $\rho = 3$.
   - On $[0, 2]$, $|g| \le \|f\|^2$ by Cauchy–Schwarz. By Bernstein, $|g| \le \rho^{2D+1}\|f\|^2$ on the ellipse, since $g$ has degree $2D+1$ in $s$.
   - The integrand's maximum is enclosed on the boundary by 2048 ball arcs, with an asserted lower bound $|1 - e^{-2Lz}| > 0$ on each arc.
   - The point $z = 0$ is a removable singularity. The other zeros, $i\pi k/L$, lie outside the ellipse, whose semi-minor axis is $4/3 < \pi$.
4. **The remaining terms of (P)** are evaluated exactly in Arb.

This gives $U \ge Q(f)/\|f\|^2 \ge \lambda_0$, with $\|f\|^2 = \sum v_m^2$. `gl_correction.py` asserts $\lambda_0^{\mathrm{lo}} \le U$ at every window.

**Form domain.** $\Psi_L(t) \sim \log t$, so the form domain is $\{f : \int\log(2+t)|\hat f|^2 < \infty\}$.

- A polynomial restricted to $[-L, L]$ has $|\hat f|^2 = O(t^{-2})$, so it lies in this domain.
- It also lies in the closure of $C_c^\infty(-L, L)$, because the weight is weaker than any $H^s$ with $s < 1/2$.

See §7.4.

### 4.5 Numerics

**Precision.** At $x = LT^\#$, Arb's Bessel evaluation amplifies input radii by roughly $e^x$, so the precision must exceed about $1.45x$ bits. The settings used were:

| Windows | Settings |
|---|---|
| $L \le 0.8$ | `CERT_PREC=768` |
| $0.84 \le L \le 0.97$ (prime 5 present, $T^\# = 830$) | `CERT_PREC=1536` |
| $0.98 \le L \le 1.0$ (prime 7, $T^\# = 3610$) | `CERT_BESSEL_PREC=6144`, `CERT_MATMUL_PREC=320`, `CERT_VCHOL_PREC=384` |

For the last group:

- nodes and Bessel values are computed at 6144 bits;
- the $N \times N$ panel products are formed at 320 bits, with rounding absorbed into the radii;
- the verified Cholesky runs at 384 bits.

The base `CERT_PREC` of that run is not recorded in the log.

**Hybrid spherical Bessel** (`bessel_rec.sph_bessel_hybrid`, used only for $L \ge 0.98$, under `CERT_BESSEL_PREC`). This computes $j_0, \dots, j_{n_{\max}}(x)$ in $O(n_{\max})$ ball operations.

- **Orders $n \le n_1 = \lfloor x\rfloor$.** Arb's own $j_s, j_{s+1}$ are used as seeds every 400 orders, followed by upward recurrence $j_{n+1} = \frac{2n+1}x j_n - j_{n-1}$ in ball arithmetic. Radii grow by at most about $2.5\times$ per step, i.e. at most about $10^{200}$ between seeds. This is harmless at 6144 bits and is tracked in the balls anyway.
- **Orders $n > n_1$.** The ratios $r_k = j_k/j_{k-1}$ come from the backward continued fraction $r_k = x/(2k + 1 - x\,r_{k+1})$. It is started from $r_{M+1} \in [0, 1]$ with $M = \max(n_{\max}, \lfloor 2.72x/2\rfloor + 60)$ and anchored at $j_{n_1}$.

**Lemma (enclosure).** For $k \ge n_1 + 1$ every true ratio $r_k$ lies in $(0, 1]$, so the interval iteration started from $[0, 1]$ encloses it.

*Proof.* If $2k + 1 \ge 2x$, the map $\tau \mapsto x/(2k + 1 - x\tau)$ sends $[0, 1]$ into $[x/(2k+1),\ x/(2k+1-x)] \subset [0, 1]$.

Now argue by downward induction from a large $K$, where $j_K(x) > 0$ and $r_K \in [0, 1]$. Suppose $2k + 1 \ge 2x$ and $r_{k+1} \in [0, 1]$. The recurrence gives
$$j_{k-1}/j_k = (2k+1)/x - r_{k+1} \ge 1 > 0,$$
so $r_k \in (0, 1]$ and $j_{k-1} > 0$. The condition $2k + 1 \ge 2x$ holds for all $k \ge n_1 + 1$, so the induction runs down to $k = n_1 + 1$. ∎

This also covers the anchor ratio $r_{n_1+1} = j_{n_1+1}/j_{n_1}$. The enclosures are rigorous.

They are not necessarily sharp at the top orders.

- For $L \ge 0.98$, $M$ equals $n_{\max}$, so there is no contraction slack above $n_{\max}$.
- In a reviewer's unrecorded spot check (2048 bits, $x = 1000.7$, $n_{\max} = 2100$), relative radii reached $0.16$ at the top orders.
- Pincherle's theorem [Ga] guarantees convergence, and hence sharpness, only when $M \gg n_{\max}$.

**Tail truncation with a drop bound.** At each node, values with $|\hat e_m| < 10^{-100}$ are set to zero, because mixed magnitudes slow FLINT's product down by about $30\times$. Each entry then changes by at most $2|s|\,F_{\max}d_{\max}$, where $s$ is the node weight times $(\Psi - \beta)/\pi$. The sum of these over all nodes is added to every entry's radius.

**Transport.** Partial sums are serialized as ball strings: 240 digits by default and 60 digits between cluster hosts. Rounding is absorbed into the radii, and parsing is outward.

**Remote code check (added after the window runs).** Before dispatching any panel group, the current `certify.py` computes SHA-256 hashes of `certify.py`, `bessel_rec.py` and `node_partial.py`. Each cluster host must pass `sha256sum -c` on byte-identical copies before it runs. This was validated end-to-end at $L = 0.4$; no log of that validation run is in the repository. The original cluster runs ($0.84 \le L \le 1$) predate this check (§7.2).

### 4.6 Sharp cutoff for $L \ge 1$ (`SHARP_TS`, `psi_floor`; extension windows)

With the cutoff rule of §4.3, $T^\#$ grows like $2\pi e^{A_L + 1/2}$: $3610$ at $L = 1$, $5890$ at $L = 1.06$, where $n = 8$ enters, and $12250$ at $L = 1.12$, after $n = 9$ has entered. The head size $N \approx eLT^\#/4$ and the Bessel precision grow with it. The prime comb reaches $-A_L$ only near isolated large $t$, so $\Psi_L$ is usually well above the envelope bound. For the windows listed in `certify.SHARP_TS` the code therefore uses a smaller cutoff:

- **Cutoff.** $T^\# = 2450$ for $L = 1.0, 1.03$, $T^\# = 4000$ for $L = 1.06, 1.09$ and $T^\# = 5550$ for $L = 1.12, 1.15, 1.18, 1.19$. The crude cutoff of §4.3 is kept as $T_2$ ($3610$, $5890$, resp. $12250$).
- **Floor $\beta$.** $\beta = 3/16$, an exact dyadic number. The code asserts $\log(T_2/2\pi) - 1/T_2 - A_L \ge 3/16$, so by Lemma (E) $\Psi_L \ge \beta$ on $[T_2, \infty)$.
- **The gap $[T^\#, T_2]$** (`psi_floor`). $\Psi_L > \beta$ is certified by Arb enclosures of $\Psi_L$ on cells of width $0.01$; a cell whose enclosure is not above $\beta$ is bisected, and the code stops with an error if a cell would become narrower than $10^{-9}$.
- Together, $\Psi_L \ge \beta$ on $[T^\#, \infty)$. That is all Lemma (R) needs, so (R), (S) and the tail bounds of §4.3 apply unchanged with this $T^\#$ and $\beta$.
- **Quadrature and head.** $NG = 112$ Gauss–Legendre nodes per unit panel (instead of 64) and $N = \lceil eLT^\#/4\rceil + 80$ head modes (instead of $+42$), so that $N \cdot$(quadrature error) and $e_B$ stay far below $\lambda_0$. The Bessel precision is raised to at least $512\lceil 1.75 LT^\#/512\rceil$ bits.

| $L$ | $n$ | $T^\#$ | $T_2$ | $N$ | Bessel bits | wall time |
|---|---|---|---|---|---|---|
| 1.0 (validation) | 2,3,4,5,7 | 2450 | 3610 | 1745 | 6144 | 45 min 46 s |
| 1.03 | 2,3,4,5,7 | 2450 | 3610 | 1795 | 6144 | 2859 s |
| 1.06 | 2,3,4,5,7,8 | 4000 | 5890 | 2962 | 7680 | 11022 s |
| 1.09 | 2,3,4,5,7,8 | 4000 | 5890 | 3043 | 7680 | 11760 s |
| 1.12 | 2,3,4,5,7,8,9 | 5550 | 12250 | 4305 | 11264 | 24794 s |
| 1.15 | 2,3,4,5,7,8,9 | 5550 | 12250 | 4418 | 11264 | ≈ 6 h 33 min + 2471 s |
| 1.18 | 2,3,4,5,7,8,9 | 5550 | 12250 | 4531 | 11776 | ≈ 6 h 58 min + 3817 s |
| 1.19 | 2,3,4,5,7,8,9 | 5550 | 12250 | 4569 | 11776 | ≈ 6 h 59 min + 4769 s |

The wall times are from `data/logs/run_extend.log` (for 1.06 and 1.09, the differences of the elapsed times printed after the resume). All windows ran with `CERT_BESSEL_PREC=6144 CERT_MATMUL_PREC=320 CERT_VCHOL_PREC=384` and the default `CERT_MARGIN`, `CERT_DDIV`; 1.12–1.19 also with the default `CERT_PREC` (768) and the triangular head assembly (Erratum, further issues, item 4). The windows up to 1.12 ran one at a time: cluster assembly of both heads, then certification on the coordinator. At 1.12 the 24794 s (25 Oct 15:12–22:06) split into about 6 h 10 min of assembly and about 44 min of local certification. For 1.15, 1.18 and 1.19 the runs were **pipelined**: `prefetch.py`, coordinated from one cluster host, assembled the heads of window $k+1$ on the cluster while the laptop certified window $k$ from the cached heads. Their wall-time entries are therefore two numbers that overlap in time with neighbouring windows and should not be added across rows: the cluster assembly time (approximate, from log timestamps) plus the elapsed time printed by the local certification. The whole continuation $1.12$–$1.19$ took from 25 Oct 15:12 to 26 Oct 19:41 (SAST), about 28.5 h.

**Validation at $L = 1$** (`validate10.py`, `data/validate_1.0.json`). The sharp path recertifies $L = 1$ with $T^\# = 2450$, $N = 1745$:
$$4.978\cdot10^{-30} \le \lambda_0 \le 9.406\cdot10^{-30}, \qquad \lambda_1^{\mathrm{even}} \ge 1.808\cdot10^{-23}, \qquad \lambda_0^{\mathrm{odd}} \ge 1.255\cdot10^{-26},$$
against the rerun's $4.531\cdot10^{-30} \le \lambda_0 \le 6.601\cdot10^{-30}$, $\lambda_1^{\mathrm{even}} \ge 1.649\cdot10^{-23}$, $\lambda_0^{\mathrm{odd}} \ge 1.141\cdot10^{-26}$ (crude $T^\# = 3610$, $N = 2496$). The brackets overlap, as they must, and the sharp lower bounds are about 10% higher. The sharp $U$ is larger; in both runs the trial vector is the approximate ground vector of the even head, truncated to 150 coefficients, and the two heads differ. The certificate of record keeps the rerun row at $L = 1$; the step $[1, 1.03]$ uses its $U$.

**Credit.** This is a refinement inside Zhu's reduction [Zhu, Thm 1.1]: only the way $\Psi_L \ge \beta$ above the cutoff is established changes. We do not know whether Zhu or others use the same device and make no novelty claim for it.

---

## 5. Results

**Table 1.** Corrected certified bounds per window, from `theorem_bounds_corrected.json`. Rounding is outward: lower bounds down, upper bounds up.

- **$\delta$** is the Weyl shift subtracted from every lower bound at that window, rounded up. It is the same for both parities to the digits shown.
- **Other columns.** The table lists the prime powers $n$ with $\log n < 2L$, the cutoff $T^\#$, and the head size $N$ per sector. "Binding" says which of $\lambda_1^{\mathrm{even}}$, $\lambda_0^{\mathrm{odd}}$ has the smaller lower bound.
- **Step ratio.** The step ratio in row $L_{i+1}$ is $\min(\lambda_1^{\mathrm{even}}, \lambda_0^{\mathrm{odd}})_{\mathrm{lo}}(L_{i+1}) / U(L_i)$. It was computed as a true lower end from the JSON enclosures and floored to 4 significant digits. Each value is at or below the lower end printed in `run_gl_correction.log`.

| $L$ | $n$ | $T^\#$ | $N$ | $\delta \le$ | $\lambda_0 \ge$ | $\lambda_0 \le U$ | $\lambda_1^{\mathrm{even}} \ge$ | $\lambda_0^{\mathrm{odd}} \ge$ | binding | step ratio $\ge$ |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.05 | – | 1600 | 97 | 1.3e-35 | 0.9351 | 0.9826 | 1.887 | 1.926 | even | (Piece I) |
| 0.1 | – | 800 | 97 | 1.3e-35 | 0.4348 | 0.4571 | 1.387 | 1.266 | odd | 1.288 |
| 0.2 | – | 400 | 97 | 1.4e-35 | 0.09014 | 0.09493 | 1.042 | 0.6037 | odd | 1.320 |
| 0.3 | – | 270 | 98 | 1.5e-35 | 7.143e-3 | 7.572e-3 | 0.7437 | 0.2110 | odd | 2.223 |
| 0.35 | 2 | 230 | 97 | 1.5e-35 | 1.110e-3 | 1.194e-3 | 0.6153 | 0.06086 | odd | 8.038 |
| 0.4 | 2 | 200 | 97 | 1.6e-35 | 1.662e-4 | 1.814e-4 | 0.2866 | 0.01361 | odd | 11.40 |
| 0.45 | 2 | 200 | 104 | 1.9e-35 | 1.465e-5 | 1.618e-5 | 0.08888 | 2.205e-3 | odd | 12.16 |
| 0.5 | 2 | 200 | 110 | 2.4e-35 | 8.482e-7 | 9.351e-7 | 0.01650 | 1.752e-4 | odd | 10.83 |
| 0.55 | 2,3 | 200 | 117 | 2.9e-35 | 4.424e-8 | 5.320e-8 | 1.234e-3 | 1.213e-5 | odd | 12.98 |
| 0.6 | 2,3 | 200 | 124 | 3.4e-35 | 1.321e-9 | 1.612e-9 | 8.661e-5 | 4.921e-7 | odd | 9.252 |
| 0.64 | 2,3 | 200 | 129 | 3.9e-35 | 7.330e-11 | 9.103e-11 | 5.914e-6 | 3.008e-8 | odd | 18.67 |
| 0.68 | 2,3 | 200 | 135 | 4.5e-35 | 1.984e-12 | 2.481e-12 | 3.870e-7 | 1.180e-9 | odd | 12.96 |
| 0.72 | 2,3,4 | 200 | 140 | 5.0e-35 | 5.422e-14 | 7.754e-14 | 1.013e-8 | 3.167e-11 | odd | 12.77 |
| 0.76 | 2,3,4 | 200 | 146 | 5.7e-35 | 7.533e-16 | 1.170e-15 | 3.574e-10 | 7.456e-13 | odd | 9.617 |
| 0.8 | 2,3,4 | 200 | 151 | 6.4e-35 | 9.787e-18 | 1.679e-17 | 4.895e-12 | 1.086e-14 | odd | 9.295 |
| 0.84 | 2,3,4,5 | 830 | 516 | 1.1e-33 | 1.134e-19 | 1.384e-19 | 8.754e-14 | 1.322e-16 | odd | 7.878 |
| 0.87 | 2,3,4,5 | 830 | 533 | 1.2e-33 | 2.636e-21 | 3.400e-21 | 2.568e-15 | 3.892e-18 | odd | 28.14 |
| 0.9 | 2,3,4,5 | 830 | 550 | 1.3e-33 | 3.401e-23 | 4.299e-23 | 5.866e-17 | 5.819e-20 | odd | 17.12 |
| 0.93 | 2,3,4,5 | 830 | 567 | 1.4e-33 | 4.303e-25 | 5.716e-25 | 8.697e-19 | 9.536e-22 | odd | 22.18 |
| 0.96 | 2,3,4,5 | 830 | 584 | 1.6e-33 | 4.510e-27 | 5.799e-27 | 9.996e-21 | 8.439e-24 | odd | 14.76 |
| 0.97 | 2,3,4,5 | 830 | 590 | 1.6e-33 | 8.936e-28 | 1.177e-27 | 2.109e-21 | 1.929e-24 | odd | 332.7 |
| 0.98 | 2,3,4,5,7 | 3610 | 2447 | 3.0e-32 | 1.681e-28 | 2.130e-28 | 5.262e-22 | 4.471e-25 | odd | 380.0 |
| 0.99 | 2,3,4,5,7 | 3610 | 2471 | 3.1e-32 | 2.779e-29 | 3.640e-29 | 1.101e-22 | 7.684e-26 | odd | 360.7 |
| 1.0 | 2,3,4,5,7 | 3610 | 2496 | 3.2e-32 | 5.147e-30 | 6.497e-30 | 1.885e-23 | 1.304e-26 | odd | 358.5 |

The gaps for Piece I on $(0, 0.05]$ are $\lambda_1^{\mathrm{even}} - \lambda_0 \ge 0.4173$ and $\lambda_0^{\mathrm{odd}} - \lambda_0 \ge 0.6625$.

**Observations** (all from the corrected table).

1. **The odd sector binds.** For $L \ge 0.1$ the binding constraint is always the **odd** sector. The certified lower bound $\ell_{\mathrm{odd}}/U$ on the ratio $\lambda_0^{\mathrm{odd}}/\lambda_0$ is about $2.77$ at $L = 0.1$ and at least $2.0\cdot10^3$ at $L = 1$. It grows overall, though not monotonically.
2. **Decay of $\lambda_0$.** $\lambda_0$ falls from $\approx 0.96$ at $L = 0.05$ to $\approx 6\cdot10^{-30}$ at $L = 1$, roughly like $10^{-cL^2}$.
3. **Bracket widths.** The relative width $U/\lambda_0^{\mathrm{lo}}$ is at most $1.72\times$.
   - The widest brackets are $1.715\times$ at $L = 0.8$, $1.552\times$ at $L = 0.76$ and $1.430\times$ at $L = 0.72$. All other windows are at most $1.33\times$, and those with $L \le 0.5$ are at most $1.11\times$.
   - The widest brackets all occur at windows with $T^\# = 200$ and $L \ge 0.72$.
   - **Where the slack lies at $L = 0.8$.** Zhu reports a non-certified window floor $\lambda_{\min}(Q|_{\mathrm{window}}) = 1.656\cdot10^{-17}$, from exact time-domain assembly [Zhu, §5.5(b)]. Our $U = 1.679\cdot10^{-17}$ is within about $1.4\%$ of it. If that floor is accurate, almost all of the factor $1.715$ lies on the lower-bound side.
   - **A plausible but unverified cause.** The reduction (R) loses more when $T^\# = 200$ is close to Zhu's threshold $T_1(0.8) \approx 119$. Zhu's own reduced forms show the same loss: $1.2\cdot10^{-18}$ at $T^\# = 150$ and $9\cdot10^{-18}$ at $T^\# = 200$ [Zhu, §5.5(b)].
   - Our $U$ equals the Galerkin value at its own degree by construction, so that agreement is not evidence of sharpness.
   - With the old, superseded $U$ the widest bracket was $3.76\times$, also at $L = 0.8$.
4. **Where the covering has room.** The steps have the most room near $L = 1$ (ratios $332.7$–$380.0$ on the 0.01-wide steps from $0.96$ on), where the grid could be coarser. The tightest steps are at small $L$ ($1.288$, $1.320$, $2.223$).

**Table 2. Certificate of record** (`data/theorem_bounds.json`: clean rerun for $L \le 1$, extension run for $L > 1$). Rounding is outward: lower bounds down, upper bounds up, to 4 significant digits. "Width" is $U/\lambda_0^{\mathrm{lo}}$, rounded up. The step ratio in row $L_{i+1}$ is $\min(\lambda_1^{\mathrm{even}}, \lambda_0^{\mathrm{odd}})_{\mathrm{lo}}(L_{i+1})/U(L_i)$, lower end floored. Values for $L \le 1$ differ from Table 1: the lower bounds slightly (the corrected quadrature radius is folded into the entries instead of a Weyl shift; at $L = 1$ also the `CERT_MARGIN`/`CERT_DDIV` change of §7.2), and $U$ more, because the rerun uses the approximate head ground vector as trial vector, not the Galerkin vector of `gl_correction.py`. $T^\#$ and $N$ for $L \le 1$ are as in Table 1; for $L > 1$ see §4.6.

| $L$ | $\lambda_0 \ge$ | $\lambda_0 \le U$ | $\lambda_1^{\mathrm{even}} \ge$ | $\lambda_0^{\mathrm{odd}} \ge$ | binding | width $\le$ | step ratio $\ge$ |
|---|---|---|---|---|---|---|---|
| 0.05 | 0.9351 | 0.9827 | 1.887 | 1.926 | even | 1.051 | (Piece I) |
| 0.1 | 0.4348 | 0.4572 | 1.387 | 1.266 | odd | 1.052 | 1.288 |
| 0.2 | 0.09014 | 0.09499 | 1.046 | 0.6037 | odd | 1.054 | 1.320 |
| 0.3 | 7.143e-3 | 7.587e-3 | 0.7437 | 0.2110 | odd | 1.063 | 2.222 |
| 0.35 | 1.110e-3 | 1.204e-3 | 0.6153 | 0.06086 | odd | 1.084 | 8.022 |
| 0.4 | 1.662e-4 | 1.841e-4 | 0.2866 | 0.01361 | odd | 1.107 | 11.31 |
| 0.45 | 1.465e-5 | 1.647e-5 | 0.08888 | 2.205e-3 | odd | 1.125 | 11.98 |
| 0.5 | 8.482e-7 | 9.506e-7 | 0.01650 | 1.752e-4 | odd | 1.121 | 10.63 |
| 0.55 | 4.424e-8 | 6.066e-8 | 1.234e-3 | 1.213e-5 | odd | 1.371 | 12.76 |
| 0.6 | 1.321e-9 | 1.813e-9 | 8.661e-5 | 4.921e-7 | odd | 1.372 | 8.114 |
| 0.64 | 7.330e-11 | 1.036e-10 | 5.914e-6 | 3.008e-8 | odd | 1.414 | 16.59 |
| 0.68 | 1.984e-12 | 2.832e-12 | 3.870e-7 | 1.180e-9 | odd | 1.427 | 11.39 |
| 0.72 | 5.422e-14 | 1.207e-13 | 1.013e-8 | 3.167e-11 | odd | 2.225 | 11.18 |
| 0.76 | 7.533e-16 | 1.802e-15 | 3.574e-10 | 7.456e-13 | odd | 2.392 | 6.181 |
| 0.8 | 9.787e-18 | 3.681e-17 | 4.895e-12 | 1.086e-14 | odd | 3.761 | 6.031 |
| 0.84 | 1.134e-19 | 1.757e-19 | 8.754e-14 | 1.322e-16 | odd | 1.550 | 3.593 |
| 0.87 | 2.636e-21 | 4.329e-21 | 2.568e-15 | 3.892e-18 | odd | 1.643 | 22.15 |
| 0.9 | 3.401e-23 | 5.545e-23 | 5.866e-17 | 5.819e-20 | odd | 1.630 | 13.44 |
| 0.93 | 4.303e-25 | 7.594e-25 | 8.697e-19 | 9.536e-22 | odd | 1.765 | 17.20 |
| 0.96 | 4.510e-27 | 7.555e-27 | 9.996e-21 | 8.439e-24 | odd | 1.675 | 11.11 |
| 0.97 | 8.936e-28 | 1.613e-27 | 2.109e-21 | 1.929e-24 | odd | 1.805 | 255.3 |
| 0.98 | 1.682e-28 | 2.165e-28 | 5.262e-22 | 4.471e-25 | odd | 1.287 | 277.2 |
| 0.99 | 2.782e-29 | 3.552e-29 | 1.101e-22 | 7.684e-26 | odd | 1.277 | 354.9 |
| 1.0 | 4.531e-30 | 6.601e-30 | 1.649e-23 | 1.141e-26 | odd | 1.457 | 321.5 |
| 1.03 | 2.096e-32 | 4.110e-32 | 1.055e-25 | 6.604e-29 | odd | 1.961 | 10.00 |
| 1.06 | 7.275e-35 | 1.249e-34 | 4.435e-28 | 2.691e-31 | odd | 1.716 | 6.548 |
| 1.09 | 1.306e-37 | 2.162e-37 | 1.034e-30 | 4.711e-34 | odd | 1.655 | 3.774 |
| 1.12 | 1.903e-40 | 3.001e-40 | 2.049e-33 | 8.502e-37 | odd | 1.577 | 3.933 |
| 1.15 | 1.752e-43 | 2.799e-43 | 2.501e-36 | 8.572e-40 | odd | 1.597 | 2.857 |
| 1.18 | 9.792e-47 | 1.636e-46 | 2.088e-39 | 6.011e-43 | odd | 1.670 | 2.148 |
| 1.19 | 7.823e-48 | 1.312e-47 | 1.439e-40 | 4.329e-44 | odd | 1.677 | 264.7 |

Piece I on the certificate of record: $\lambda_1^{\mathrm{even}} - \lambda_0 \ge 0.4172$ and $\lambda_0^{\mathrm{odd}} - \lambda_0 \ge 0.6624$ on $(0, 0.05]$, with $C_{\mathrm{even}} \le 4.876362$, $C_{\mathrm{odd}} \le 0.761666$. The smallest step ratios are $1.288$ ($[0.05, 0.1]$), $1.320$ ($[0.1, 0.2]$) and, past $L = 1$, $2.148$ ($[1.15, 1.18]$); past $L = 1$ the margins shrink step by step from $[1.09, 1.12]$ on ($3.933$, $2.857$, $2.148$; the short step $[1.18, 1.19]$ gives $264.7$). The log prints midpoints (e.g. $3.934$ at $[1.09, 1.12]$); the table floors lower ends. The odd sector binds at every window from $L = 0.1$ on.

---

## 6. Reproducibility

### Software

- **Python side.** Python 3.12.13 with python-flint 0.9.0. The Python-side Arb is the FLINT 3.6.0 bundled in the wheel (`flint.__FLINT_VERSION__` reports 3.6.0).
- **C tools.** `vchol` and `veig` are built with gcc 15.2.0 against the system FLINT 3.4 (`libflint.so.22`), using `gcc -O2 -o vchol vchol.c -lflint -lmpfr -lgmp` and likewise for `veig`.
- **Cluster.** A 3-node cluster with 64 cores per node, reached over ssh. Its FLINT and python-flint versions were not recorded; the SHA-256 check covers the code files, not the libraries.

### What the certificate consists of

1. **Certificate of record.** `data/theorem_bounds.json` (31 windows) and the gluing in `code/theorem.py` (Piece I and Piece II). The rows $L \le 1$ are copied verbatim from the clean rerun, the rows $1.03$–$1.09$ verbatim from the extension run (pre-triangular code), and the rows $1.12$–$1.19$ verbatim from its continuation (triangular code). The continuation's own bounds file repeats the older rows; they differ from the certificate of record only in the last digit of some radii (the certificate of record's are the same or slightly wider), and the step $[1.09, 1.12]$ was glued with the certificate-of-record row at $1.09$. Provenance is in `data/README.md`. `python code/theorem.py --check` re-checks the gluing (`data/logs/check.log`).
2. **History.** The original window runs (`run_theorem*.log`, not released) produced `data/history/theorem_bounds_original.json` with successive revisions of the code, all with the Gauss–Legendre off-by-one (§7.2). `code/gl_correction.py` corrected them (`data/history/theorem_bounds_corrected.json`, `data/logs/run_gl_correction.log`). It ran with the code revision of the clean rerun, which differs from `code/` in the sharp-cutoff path (§4.6) and in `node_partial.py` (§Erratum, further issues); re-running it with `code/` would use the sharp parameters at $L = 1$.

### Commands

All commands run in `code/` after `make` (which builds `vchol` and `veig`), with the Python environment of `requirements.txt`.

| Command | What it does |
|---|---|
| `python theorem.py --check` | Piece I and all Piece II steps on `../data/theorem_bounds.json`, over the longest certified prefix of the grid. A few seconds; no window is recomputed. |
| `python lipschitz.py` | Piece I constants at the default 53-bit precision. The log prints the enclosures $[0.592772 \pm 2.65\cdot10^{-7}]$, $[4.87483 \pm 4.86\cdot10^{-6}]$ and $[0.760129 \pm 1.10\cdot10^{-7}]$, i.e. $\int\|h\| \le 0.592773$, $C_{\mathrm{even}} \le 4.87484$, $C_{\mathrm{odd}} \le 0.760130$. At 768 bits (inside `theorem.py` / `gl_correction.py`) the values are 0.593979, 4.876362 and 0.761666 (§4.1). |
| `python certify.py [L ...]` | Single windows. The default $L = 0.4, 0.5, 0.6, 0.7, 0.8$ includes the Zhu Thm 1.2 check. |
| `UPTO=0.8 CERT_PREC=768 python theorem.py` | Stage 1 of the rerun: windows $0.05$–$0.8$ (local). Windows already in `./theorem_bounds.json` are reused. |
| `UPTO=0.97 CLUSTER=user@host1,user@host2,user@host3 CERT_PREC=1536 python theorem.py` | Stage 2: windows $0.84$–$0.97$. Panels are dispatched over ssh to `node_partial.py` on each host, in the directory `CERT_REMOTE_DIR` (default `~/weil`), after a SHA-256 check of the three code files. |
| `CLUSTER=... CERT_BESSEL_PREC=6144 CERT_MATMUL_PREC=320 CERT_VCHOL_PREC=384 python theorem.py` | Stage 3 and the extension: windows $0.98$ onward. In the rerun, $L = 1$ needed `CERT_MARGIN=1.2 CERT_DDIV=8` (§7.2). With the released code, $L = 1$ is computed by the sharp path (§4.6). |
| `CLUSTER=... CERT_BESSEL_PREC=6144 CERT_MATMUL_PREC=320 CERT_VCHOL_PREC=384 python prefetch.py 1.15 1.18` | Assembles the heads of the listed windows into `cache/` on the cluster (both sectors concurrently), with the same environment as the certifying run; a later `theorem.py` run then certifies those windows from the cache. Used to pipeline $1.15$–$1.19$ (§4.6). |
| `python validate10.py` | The sharp-path validation at $L = 1$ (§4.6). Reads `./theorem_bounds.json` for the comparison. |
| `python gl_correction.py` | Historical: the Gauss–Legendre correction and re-check of the original runs (reads `./theorem_bounds.json`, writes `theorem_bounds_corrected.json`). To re-run it, first copy `../data/history/theorem_bounds_original.json` to `./theorem_bounds.json`. It now imports the released `certify.py`, not the revision it originally ran with (§6, item 2), so byte-exact reproduction of `theorem_bounds_corrected.json` is not guaranteed. |

### Hardware and run times

The local machine (a 20-thread laptop) uses a pool of 20 processes (12 in `gl_correction.py`). The cluster has 3 nodes with 64 cores each, reached over ssh. The table lists the original runs; the clean rerun and the extension follow below it.

| Run | Windows | Where | Wall time |
|---|---|---|---|
| `run_certify.log` | 0.4, 0.5, 0.6, 0.7, 0.8 + Zhu Thm 1.2 | local | 10 min 41 s |
| `run_theorem.log` | 0.05 – 0.8 (15) | local | 32 min 40 s |
| `run_theorem_09.log` | 0.84, 0.87, 0.9 | cluster | 92 min 48 s (≈ 29–33 min per window) |
| `run_theorem_097.log` | 0.93, 0.96, 0.97 | cluster | 87 min 12 s (≈ 27–30 min per window) |
| `run_theorem_10.log` | 0.98, 0.99, 1.0 | cluster | 266 min 27 s (4917 s, 5054 s, 6008 s) |
| `run_gl_correction.log` | all 24 (correction, new $U$, re-check) | local | 1 min 16 s |

Clean rerun (`data/logs/run_rerun.log`): stage 1 local, 18 min 20 s; stage 2 on the cluster, 70 min 48 s; stage 3 on the cluster, 247 min 20 s until $L = 1$ failed with the default margins, then 10 min 34 s to finish $L = 1$ (from the cached heads) with `CERT_MARGIN=1.2 CERT_DDIV=8` and the gluing. Extension (`data/logs/run_extend.log`): validation at $L = 1$ 45 min 46 s; $L = 1.03$ 2859 s; the first attempt at $1.06$ ran out of memory after 230 min 41 s in total; after the fix, $1.06$ and $1.09$ took 11022 s and 11760 s. Continuation with the triangular code: $1.12$ took 24794 s; $1.15$, $1.18$, $1.19$ were pipelined (cluster assembly about 6.5–7 h per window, local certification 2471 s, 3817 s, 4769 s; §4.6); the continuation as a whole ran from 25 Oct 15:12 to 26 Oct 19:41. The grid ends at $1.19$; nothing past it is certified.

In the cluster runs the `user` times count only the coordinator. For $L \ge 0.98$ the head matrices are cached as ball strings in `cache/head_<L>_<parity>_<N>.txt`; their size grows like $N^2$, about 160 MB each at $L = 1.0$ and 1.06 GB each at $L = 1.19$.

### Cross-checks

1. **Zhu Thm 1.2 reproduced** (independent code, same normalization and cutoff strategy).
   - At $L = 0.8$ we certify $\lambda_0 \ge 9.787\cdot10^{-18} \ge 8.9\cdot10^{-18}$ on the dyadic window $\tilde L \ge 0.8$. By Lemma 2.3 this gives $Q \ge 8.9\cdot10^{-18}\|f\|^2$ for even $f$ supported in $[-0.8, 0.8]$.
   - The code's own assertion is weaker: `certified_pd(A_even, 8.9e-18)` and $\min(8.9\cdot10^{-18}, \beta - e_D) - e_B > 0$. That proves only $Q \ge (8.9\cdot 10^{-18} - e_B)\|f\|^2$.
2. **Consistency with Zhu at $L = 0.8$** (Thm 6.2 and Table 3). The second even value is a rigorous lower bound in both columns. Zhu's own non-rigorous reference value for his reduced head is $2.3175\cdot10^{-12}$. That is below our certified $4.895\cdot10^{-12}$, but there is no contradiction: an eigenvalue of a reduced form only bounds the true value from below, and Zhu's reduced matrices use different cutoffs and heads ($T^\# = 150$ or $200$, 200 modes).

   | Quantity | Zhu | Ours | Agreement |
   |---|---|---|---|
   | $\lambda_0$ | $[8.9\cdot 10^{-18},\ 2.523\cdot 10^{-16}]$ [Thm 6.2]; $\le 2.27\cdot10^{-17}$ [Table 3] | $[9.787\cdot 10^{-18},\ 1.679\cdot 10^{-17}]$ | Table 1 bracket lies inside Zhu's sharper enclosure $[8.9\cdot10^{-18}, 2.27\cdot10^{-17}]$; the Table 2 upper bound $3.681\cdot10^{-17}$ is weaker than Zhu's |
   | $\lambda_1^{\mathrm{even}}$ | $\ge 2.085\cdot 10^{-12}$ [Thm 6.2] | $\ge 4.895\cdot 10^{-12}$ | consistent |
   | $\lambda_0^{\mathrm{odd}}$ | $[8.206\cdot 10^{-15},\ 2.347\cdot 10^{-14}]$ [Thm 6.2] | $\ge 1.086\cdot 10^{-14}$ | consistent |

   (Table 1 values.) At the other shared windows the corrected original $U$ also lies below Zhu's certified Table 3 bounds: $9.351\cdot10^{-7}$ against $1.05\cdot10^{-6}$ at 0.5; $1.612\cdot10^{-9}$ against $1.82\cdot10^{-9}$ at 0.6; $4.299\cdot10^{-23}$ against $5.91\cdot10^{-23}$ at 0.9; and $6.497\cdot10^{-30}$ against $8.08\cdot10^{-30}$ at 1.0. With the certificate of record (Table 2) the values are $9.506\cdot10^{-7}$, $1.813\cdot10^{-9}$, $5.545\cdot10^{-23}$ and $6.601\cdot10^{-30}$ (below Zhu), and $3.681\cdot10^{-17}$ at 0.8 (above Zhu).
3. **Non-rigorous Galerkin values** (see `REDUCTION.md`):
   - The values are $\lambda_0 \approx 9.350\cdot 10^{-7}$ at $0.5$ ($D_G = 66$, `run_step3.log`), $1.679\cdot 10^{-17}$ at $0.8$ ($D_G = 78$, `run_step3.log`), $4.140\cdot 10^{-23}$ at $0.9$ ($D_G = 124$, `run_step3_09.log`), and $5.885\cdot 10^{-30}$ at $1.0$ ($D_G = 132$, `run_step3b.log`).
   - Each lies inside our corrected bracket.
   - At $0.5, 0.6, 0.8, 0.9, 1.0$ the new $U$ equals the `run_step3.log` Galerkin value at the degree used by `gl_correction.py` ($D_G = 66, 70, 78, 84, 92$): $9.350\cdot10^{-7}$, $1.611\cdot10^{-9}$, $1.679\cdot10^{-17}$, $4.299\cdot10^{-23}$ and $6.497\cdot10^{-30}$. This holds by construction, since $U$ is the Rayleigh quotient of that Galerkin vector.
   - The Galerkin values are not converged in $D_G$. At $L = 1$, going from $D_G = 92$ to $132$ lowers the value by about $9\%$.
   - `REDUCTION.md` records that these values and $\lambda_1/\lambda_0$ agree with Zhu's Table 1. We did not re-check that comparison for this note.
4. **Explicit formula test.** `weil_z.py` checks the prime side of $Q$ against the zero side (first 400 zeros) for the Legendre basis. The position-side formula (P) was checked against `weil_z.py`.
5. **Internal assertions.** The code asserts:
   - $\beta > 0$;
   - $r < 1$ for the tail sums;
   - no $\log n$ overlaps $2L$;
   - a positive denominator on every arc of the Rayleigh ellipse;
   - $\lambda_0^{\mathrm{lo}} \le U$ at every window (also asserted on the corrected numbers in `gl_correction.py`);
   - $h(0) = 1$.

   Across the grid, $\lambda_0^{\mathrm{lo}}(L_{i+1}) \le U(L_i)$ always holds, as monotonicity requires.
6. **Bessel validation.** The `__main__` test of `bessel_rec.py` checks `sph_bessel_all` against Arb's $j_n$ at $x = 0.37, 12.5, 160.25, 3000.5$.
   - No recorded test exists for `sph_bessel_hybrid`, the routine used for $L \ge 0.98$.
   - A reviewer's unrecorded spot check at $x = 200.3$ and $x = 1000.7$ found it consistent with Arb.
   - Its rigor rests on ball arithmetic and the enclosure lemma (§4.5).
7. **Repeated windows agree.** Windows computed in several runs (0.4, 0.5, 0.6, 0.8) give identical bounds in `run_certify.log` and `run_theorem.log`.
8. **Window parameters (partial).** For $L \ge 0.98$, the head sizes $N = 2447, 2471, 2496$ recorded in the cache filenames of the original run equal the $N$ computed by the final code. The final code also gives $T^\# = 830$, $N = 516$–$590$ for $0.84 \le L \le 0.97$, matching Table 1, but no record of the original runs' values exists for those windows.

---

## 7. Limitations and caveats

1. **Computer-assisted, unrefereed, AI-assisted.** Correctness depends on:
   - the code: about 460 lines in `certify.py`, plus `lipschitz.py`, `theorem.py`, `gl_correction.py`, `convergence.py`, `bessel_rec.py` and `vchol.c`;
   - Arb/FLINT;
   - the lemmas in §4.

   The note, the code and all reviews so far were produced with AI assistance. No independent human party has reviewed the code or the proofs.
2. **Code revisions between the original runs; clean end-to-end rerun.**
   - **Successive revisions.** The window runs finished on 2026-10-22 ($L \le 0.97$; the cluster runs ended at 19:05 and 20:56) and on 2026-10-23 at 09:10 ($L \ge 0.98$). File timestamps show these changes in between:
     - **2026-10-22 19:30 and 19:31:** `bessel_rec.py`, then `veig.c`/`veig`. This is about half a minute after `run_theorem_097` started (about 19:29). `bessel_rec` is imported only under `CERT_BESSEL_PREC`, i.e. only by the $L \ge 0.98$ run. `veig` supplies only non-rigorous estimates.
     - **2026-10-22 20:57:** `theorem.py`, after `run_theorem_097`. It does gluing only, and that gluing is re-done by `gl_correction.py`.
     - **2026-10-23:** `node_partial.py` (03:04), `vchol.c`/`vchol` (04:38) and `certify.py` (04:43, and again 11:07). The 11:07 change is the Gauss–Legendre fix plus the SHA-256 check.

     So no window was certified with the final `certify.py`. The windows $0.84$–$0.97$ used earlier revisions of `certify.py`, `node_partial.py` and `vchol`, the last because $N = 516$–$590 > 120$.
   - **How the correction bridges this.** `gl_correction.py` uses the stored bounds and recomputes only $\delta$, $U$, Piece I and Piece II with the final code. Its $\delta$ uses the window parameters ($T^\#$, $N$, $\beta$, $M$) as the final code computes them. We assume these equal those of the original runs, since the parameter code path is meant to be unchanged. This is checked only for $N$ (and hence $T^\#$) at $L \ge 0.98$ (§6, cross-check 8). The correction of the lower bounds is rigorous only under this assumption. One could remove the assumption by reading the entry radii directly from the cached heads ($L \ge 0.98$), but we have not done so.
   - **Remote code in the original runs.** In the cluster runs $0.84 \le L \le 1$, the coordinator did not check which code the remote hosts ran. The only check on remote output was a line count. The current code refuses to run on a host whose `certify.py`, `bessel_rec.py` or `node_partial.py` differs by SHA-256 (§4.5), but that check post-dates those runs.
   - **Clean end-to-end rerun (done, 2026-10-24).** Directory `rerun/`, script `run.sh`, log `run_rerun.log`. Fresh copies of the final code, with no cache and no `theorem_bounds.json` carried over; `vchol`/`veig` were rebuilt from source; the SHA-256 hashes of `certify.py`, `bessel_rec.py` and `node_partial.py` matched on all three hosts. The run has three stages with the precisions listed in §6.
     - All 24 windows certified. Every Piece II step holds; the smallest ratio is $\ge 1.288$ at $[0.05, 0.1]$ and the step $[0.99, 1]$ gives $\ge 321.5$ (floored lower ends; the log shows midpoints). Piece I holds, with gaps $\ge 0.4172$ (even) and $\ge 0.6624$ (odd).
     - Lower bounds agree with `theorem_bounds_corrected.json` to three significant digits at every window except $L = 1$. Upper bounds $U$ differ because a different, equally valid, trial vector is used.
     - **One parameter change at $L = 1$.** With the corrected Gauss–Legendre bound, entry radii are $\rho^2 \approx 4.33$ times larger than in the original runs (about $1.6\cdot10^{-35}$ instead of $3.7\cdot10^{-36}$). The Frobenius residual test in `vchol` then exceeded its slack $d = \mu/100$, so $L = 1$ failed at first ("could not certify"). The fix: the slack is now a parameter (`CERT_DDIV`, default 100), and so is the margin below the estimate (`CERT_MARGIN`, default 1.05). $L = 1$ was certified with `CERT_MARGIN=1.2 CERT_DDIV=8`, i.e. $\mu = \hat\lambda_0/1.2$ and $d = \mu/8$. The test is unchanged: $\|E\|_F < d$ still implies $A - \mu I \succeq 0$. The cost is a weaker bound at $L = 1$ only: $\lambda_0 \ge 4.532\cdot10^{-30}$, $\lambda_1^{\mathrm{even}} \ge 1.650\cdot10^{-23}$ and $\lambda_0^{\mathrm{odd}} \ge 1.142\cdot10^{-26}$, about 12% below the corrected original values. All other windows used the defaults. The certificate change is local (`vchol` runs on the coordinator), so the remote-hashed assembly code is unaffected.
3. **Normalization.** Our $Q$ follows Zhu [Zhu, §6] and was validated against the zero side of the explicit formula. It has **not** been checked line by line against the $QW_\lambda$ of Connes–van Suijlekom / CCM ($\lambda = e^L$). The conclusion is invariant under $Q \mapsto aQ + b\|\cdot\|^2$ with $a > 0$. A discrepancy in the relative weights of the pole, archimedean or prime terms would not be covered by that invariance.
4. **Form domain.** The operator is taken to be the one associated with the closure of $Q$ on $C_c^\infty(-L, L)$. Discreteness is proved in §2.3. The Rayleigh trial functions are polynomials that do not vanish at $\pm L$. We argue in §4.4 that they lie in the form domain, since the log-weighted space is weaker than $H^s$ for every $s < 1/2$. If a different self-adjoint realization is meant, this point should be re-examined.
5. **Lipschitz step.** Piece I uses the min-max Lipschitz property for the family $R_L$ on the common core. The derivative bounds are crude but explicit, and the argument assumes the form domain does not depend on $L$, which holds because the weights are comparable.

   The certified upper end of $\int|h|$ grows with working precision: $\le 0.592773$ at 53 bits, $\le 0.593583$ at 128 bits, $\le 0.593979$ at 768 bits. All of these are meant to be upper bounds, and we use the largest. Higher precision should not normally give a larger upper bound, and this behaviour is not yet explained. It should be investigated by evaluating $h$ on subdivided cells.
6. **Latent comparison fragility.** Python `min`/`max` on Arb balls are not rigorous when the balls overlap. They occur in the step checks of `theorem.py` and `gl_correction.py`, and in `certify.py` in the $\min(\mu, \beta - e_D)$ steps and the Bessel-continuation bound. Every current use is separated by many orders of magnitude more than the ball radii (§4.2, §4.3). The code should use explicit endpoint comparisons.
7. **Sources not read.**
   - Liu's alphaXiv paper was not located; its contents are reported as positivity only.
   - Connes' CIRM 2026 slides could not be fetched.
   - Eberwein's ledger proofs were not inspected.
   - Zhu's Thm 1.1, Lemma 3.1, Thm 6.2, §5.5, §6 and Table 3 were read from a text extraction of arXiv:2608.24827v2.
   - [PL294] and Andrade's remark are unrefereed, so the overlap on $(0, 0.1406]$ should be read with that in mind.
   - Whether the mathia notes are AI-generated is not documented in `NOVELTY.md`, and we have not checked it.
8. **Reason for extra care near and past $L \approx 1$.** Zhu §7 retracted an earlier certified-positivity claim at $L = 1.19$. That is exactly the last window of our grid: our bracket there ($\lambda_0(1.19) \ge 7.823\cdot10^{-48}$) implies positivity at that $L$, rests on the unreviewed sharp-cutoff and triangular code, and is the window most in need of independent checking. One constant error has already been found and corrected here (Gauss–Legendre off-by-one). The following bounds near $L \approx 1$ deserve an independent audit:
   - $e_D$ and $e_B$;
   - the $10^{-100}$ drop bound;
   - the 60-digit transport;
   - the hybrid Bessel seeds;
   - `psi_floor` and the sharp-cutoff parameters (§4.6), used only for $L > 1$ in the certificate of record. This path and the windows $1.03$–$1.19$ have not been through any review round, internal or external (the review log, §10, predates them);
   - the triangular head assembly in `_partial`/`node_partial.py` (Erratum, further issues, item 4), used only for $1.12$–$1.19$. Its only check is the one-off byte-identical comparison with the old code at three small cases.

   The step ratios near $L = 1$ ($255.3$–$354.9$ in Table 2) leave a large margin, but past $L = 1$ they are smaller, and from $[1.09, 1.12]$ on shrinking ($10.00$, $6.548$, $3.774$, $3.933$, $2.857$, $2.148$), and margin does not make up for a wrong bound.
9. **Stale docstrings.** Where they differ from this note, this note supersedes them:
   - the `head()` docstring, which gives the cache filename without `_<N>`;
   - the `bessel_rec.py` module docstring, "STATUS: prototype". It refers to `sph_bessel_all`; the re-seeded hybrid is the fix;
   - the `vchol.c` header, which says $d = \lambda/100$; the divisor is now the optional fourth argument (`CERT_DDIV`, default 100).

   Development notes not in this release (`README.md`, `REDUCTION.md`, `NOVELTY.md`) predate this note. `REDUCTION.md`'s statement "Certified for every $L \in (0, 0.8]$" and its reliance on an unverified survey "Fact 6.4" are not part of this theorem.
10. **Scope.** The theorem says nothing for $L > 1.19$. The prime $11$ enters at $L = \tfrac12\log 11 \approx 1.1989$, just past the last window, and raises $A_L$ and hence the cutoff; going past it needs new parameters and new runs. The step margins past $L = 1.09$ are shrinking ($3.933$, $2.857$, $2.148$). Each further window costs about 7 h of cluster assembly plus about an hour of local certification, and memory per worker still grows like $N^2$ (Erratum, further issues).
11. **Future work, not part of this theorem.** `REDUCTION.md` reduces a convergence statement towards RH to the simple-even hypothesis plus a gap-ratio condition along some sequence $L_j \to \infty$. That reduction is conditional on unverified normalization and on the survey's Fact 6.4. This note supplies the hypothesis only on $(0, 1.19]$ and makes no claim about that programme.

---

## 8. Conclusion

We certify, by Arb ball arithmetic, that the full localized Weil form for $\zeta$ on $L^2(-L, L)$ has a simple ground state with even eigenfunction for every $L \in (0, 1.19]$. The part we believe to be new is $(0.8, 1.19]$; simple-even on $(0, 0.8]$ is claimed in [KB].

**Beyond $L = 1$ (non-rigorous, Appendix A).** The staircase states are, to $1 - \cos \le 10^{-4}$, the Gram–Schmidt derivative ladder $k_L', k_L'', \dots$ of the ground state, which explains the parity alternation. Oscillation (Sturm) structure fails from the 5th state. Every all-$L$ reformulation we tried hides a scalar margin that shrinks about as fast as $\lambda_0$ (e.g. $1 - F \approx 5\cdot10^{-25}$ at $L = 1$, first-order estimate, A.3), which is why soft arguments cannot extend the theorem to all $L$.

**Margins** (certificate of record, Table 2):

- Piece I gaps: $\ge 0.4172$ (even) and $\ge 0.6624$ (odd).
- Piece II step ratios: $\ge 1.288$; past $L = 1$, $\ge 2.148$.
- Two-sided $\lambda_0$ brackets at 31 windows, each of relative width at most $3.761\times$ (at most $1.72\times$ at the 24 windows of the corrected original certificate, Table 1).

**Erratum.** A Gauss–Legendre error constant was off by one node in the original runs. It is fixed in the code. The certificate was corrected by a Weyl shift of all lower bounds by $\delta \le 3.2\cdot10^{-32}$ and freshly certified $U$ values. Under the assumption that the window parameters are unchanged between code revisions (§7.2), the theorem holds with the corrected numbers. The clean rerun (§7.2) confirms it without that assumption.

**Conditional ingredients and open items:**

- A clean end-to-end rerun with the final code certifies the theorem on its own (§7.2). It removes the earlier conditions: the $\delta$ parameter assumption and the unverified remote code apply only to the original runs.
- The identification with the CvS/CCM form is checked only up to positive affine rescaling (§7.3).
- The choice of self-adjoint realization (§7.4).
- An unexplained precision dependence of the Piece I integral bound (§7.5).
- Correctness of Arb/FLINT and of the code.

**Not claimed as new:**

- the range $(0, 0.8]$ [KB], and in particular $(0, 0.1406]$ [PL294];
- the window $L = 0.8$ by itself [Zhu, Thm 6.2];
- certified upper bounds on $\lambda_0$ for $L \in [0.5, 2]$ and the two-sided enclosure at $L = 0.8$ [Zhu, §1.4, Table 3];
- the frequency-cutoff reduction and its head/tail block estimate [Zhu, Thm 1.1];
- the perturbative neighbourhood of $0.8$ [PL287];
- positivity on this range;
- the expected spectral ordering.

---

## 9. References

**Status labels.**

- [preprint]: arXiv or alphaXiv, not refereed.
- [deposit]: Zenodo, not refereed.
- [repo]: GitHub repository, not refereed.
- [GitHub notes]: project notes in a GitHub repository, not refereed. `NOVELTY.md` describes the mathia notes as "analytic, unrefereed notes". We do not know whether they are AI-generated.
- "AI-audited" and "candidate" are the sources' own descriptions as recorded in `NOVELTY.md`.

We have not verified exact titles or years for most entries. Where a title was not recorded, we give a description in plain text, not as a title.

- **[An1]** Andrade, Zenodo 20682834 [deposit; unrefereed]. Perron-type theorem for the pole-removed Weil form; remark with the Beurling–Deny threshold $t^* = 0.28120$ for the full kernel.
- **[An2]** Andrade, Zenodo 20887028 [deposit; unrefereed]. Certified simple-even for finite Galerkin matrices $QW^N$.
- **[KB]** Kramarenko-Byrd, Zenodo 23007100, v1 of 2026-09-27 [deposit; unrefereed]. Claims rigorous simple-even (and real Fourier zeros of the ground state) for all $0 < L \le 4/5$.
- **[CCM]** A. Connes, C. Consani, H. Moscovici, arXiv:2511.22755 [preprint]. Thm 3.6 (discrete spectrum); §8 (simple-even as a missing step).
- **[CvS]** A. Connes, W. van Suijlekom, arXiv:2511.23257 [preprint].
- **[Co]** A. Connes, arXiv:2602.04022 [preprint]. Prolate construction, eq. (17); cited for context and for the non-rigorous `code/convergence.py`; not used in the proof.
- **[Eb]** Eberwein, Zenodo 23177438 [deposit; unrefereed]. Numerical simple-even observations, $L \approx 0.55$–$1.15$.
- **[El]** Elder, Zenodo 23113186 [deposit; unrefereed]. Per-sector lower bounds; treats simple-even ("Step A") as open.
- **[Els]** Elston, *weil-windows* [repo; unrefereed]. Per-sector lower bounds at $L = 0.8$, $1.19$, $1.3$, $1.495$.
- **[Eu]** "Eureka", arXiv:2608.19047, §5.3 [preprint; author as listed in `NOVELTY.md`]. Self-described interval-Cholesky positivity "candidate".
- **[Ga]** W. Gautschi, *Computational aspects of three-term recurrence relations*, SIAM Review 9 (1967).
- **[JPB]** jose-paul-b, *weil-positivity-certificates* [repo; "AI-audited"; unrefereed].
- **[Mo]** Mori, Zenodo 23281236 and related deposits [deposits; unrefereed]. Positivity at $L = (\log 10)/2$, $(\log 12)/2$, $(\log 15)/2$.
- **[MS]** grandchallenge / MATHSOLVE [unrefereed; not read in full]. Reported analytic argument for $L \le 9/50$.
- **[Liu]** Liu, positivity/coercivity of the Weil form at $L = 1, 17/16$ [alphaXiv preprint; not located, not read].
- **[PL282], [PL287], [PL294]** mathia project notes PL-282, PL-287, PL-294, GitHub `murillo128/mathia` [GitHub notes; unrefereed].
- **[Su]** M. Suzuki, arXiv:2606.09096 [preprint]. Thm 1.4 and §5.
- **[Tr]** L. N. Trefethen, *Approximation Theory and Approximation Practice*, SIAM. Thm 19.3, stated for the $(n+1)$-point rule.
- **[Wa]** G. N. Watson, *A Treatise on the Theory of Bessel Functions*, §15.3. No longer needed for the enclosure lemma.
- **[Zhu]** Zhu, arXiv:2608.24827 v2 [preprint]. Thm 1.1 (one-stroke reduction with head/tail block estimate, proved in §4), Lemma 3.1 (envelope), Thm 1.2, §1.3–1.4 and §8 (certified upper bounds; two-sided enclosure at $L = 0.8$), §5.5 (window floor), Thm 6.2, Cor 6.3, §6 (normalization, odd sector, certificates at $L = 0.8$), §11.4 and Table 3 (certified upper bounds for $L \in [0.5, 2.0]$), Table 1, §7 (retraction at $L = 1.19$).
- Also listed in `NOVELTY.md`, with no overlap on simplicity: Groskin, arXiv:2605.20224 and 2607.02828 [preprints]; Mehennaoui, Zenodo 23442026 [deposit; unrefereed]; Caceres, Zenodo 23070557 [deposit; unrefereed].
- **[We]** A. Weil, *Sur les « formules explicites » de la théorie des nombres premiers*, Comm. Sém. Math. Univ. Lund, tome supplémentaire (1952), 252–265. The explicit formula and the quadratic form.
- **[Se]** A. Senke, *Exact Finite Gram Closure for Weil Positivity on the Fixed Window (−77/64, 77/64)*, Zenodo 22871384, v0.3 (2026-09-21); unrefereed, restricted-access package. Positivity only.
- **[AnJ]** R. Andrews Jr., *Independent Reproduction and Convergence Analysis of the CCM Zeta Spectral Triple*, Zenodo 22715669, v2.5 (2026-09-11); unrefereed. Certifies simple-even ground states of *finite* matrices at $\lambda^2 = 13$ ($L \approx 1.282$); not a continuum statement.
- **Software:** FLINT/Arb (F. Johansson and the FLINT developers; F. Johansson, "Arb: efficient arbitrary-precision midpoint-radius interval arithmetic", IEEE Trans. Comput. 66(8) (2017), 1281–1292); python-flint 0.9.0 (bundling FLINT 3.6.0); system FLINT 3.4 for the C tools; mpmath (F. Johansson et al.) for the non-rigorous computations and for the zeta zeros (`mpmath.zetazero`) used in `weil_z.py`; NumPy.

---

## Appendix A. Structural observations and obstructions (non-rigorous)

Everything in this appendix is **non-rigorous**. It uses high-precision Galerkin computations (mpmath, 120–150 digits, Legendre basis on $[-L, L]$, `convergence.galerkin`), not Arb certificates. Its aim is to record what we learned when we looked for a mechanism that would give simple-even for **all** $L$ at once, including what failed. Scripts are listed in each item.

### A.1 The derivative ladder (positive observation)

Let $k_L$ be the even ground state and $\psi_1, \psi_2, \dots$ the next states in order of eigenvalue. Gram–Schmidt the derivatives $k_L', k_L'', \dots$ within each parity sector against the lower ladder elements, giving $w_1, w_2, \dots$. Then $w_j \approx \psi_j$ (`ladder.py`). The baseline ladder $x^j k_L$ does much worse.

| $L$ | $D$ | $1 - \lvert\langle w_j, \psi_j\rangle\rvert$ for $j = 1, \dots, 5$ (derivatives) | same, baseline $x^j k_L$ |
|---|---|---|---|
| 0.6 | 70 | 2.4e-5, 9.2e-4, 0.77, 0.91, 0.94 | 7.2e-3, 7.0e-2, 0.35, 0.94, 0.29 |
| 0.8 | 78 | 2.5e-6, 2.2e-5, 1.2e-4, 3.7e-2, 0.96 | 4.3e-3, 3.9e-2, 0.18, 0.57, 0.79 |
| 1.0 | 92 | 3.5e-7, 2.6e-6, 1.3e-5, 6.2e-5, 2.5e-4 | 3.1e-3, 2.8e-2, 0.12, 0.40, 0.95 |
| 1.2 | 118 | 9.7e-7, 3.6e-6, 2.4e-5, 6.7e-5, 2.2e-4 | 2.4e-3, 2.2e-2, 9.6e-2, 0.31, 0.79 |
| 1.4 | 160 | 1.7e-7, 6.3e-7, 5.5e-6, 1.5e-5, 8.3e-5 | 2.1e-3, 1.9e-2, 8.1e-2, 0.26, 0.68 |

- The ladder holds for as many steps as there are "tiny" eigenvalues (the staircase): 2 steps at $L = 0.6$, 3–4 at $0.8$, and all 5 tested from $L = 1$ on. It improves as $L$ grows.
- **It explains the strict even/odd alternation of the staircase** (fact 3 in §5): $d/dx$ flips parity.
- **Heuristic reason.** On the Fourier side, $\mathrm{span}\{k_L, k_L', \dots, k_L^{(j)}\} = \{p(t)\,\hat k_L(t) : \deg p \le j\}$. Every such function vanishes wherever $\hat k_L$ does, i.e. (by CCM and our §5 observation) very nearly at the zeta zeros. By the explicit formula, $Q(f)$ is small when $\hat f$ nearly vanishes at the zeros. The polynomial factor $t^j$ amplifies the residual at high zeros, which is why the ladder eigenvalues increase with $j$.
- **Relation to prior work.** CCM use the commutator of $Q_L$ with $D = d/dx$ (2511.22755, §5, Lemma 5.2). The ladder is the spectral shadow of that structure. Neither CCM nor the other sources we checked report the numerical ladder overlaps, so we record them. A related sub-claim from internal review: $\psi_1 \to k_L'/\lVert k_L'\rVert$ in $L^2$ as $L \to \infty$.
- **Caveat from refereeing (2026-10-24).** The ladder is an $L^2$ statement only. Its energies are not converged: at $L = 0.6$, $Q(w_2) = 5.5\cdot10^{-4}, 1.6\cdot10^{-2}, 1.29$ for $D = 40, 80, 160$. The edge values $k_L(\pm L)$ also drift down with $D$. Both fit a log-type boundary layer, $k_L \sim (\log 1/d)^{-1/2}$ at distance $d$ from $\pm L$, as for log-Laplacian forms, so that $k_L'$ is probably not in the form domain. Three proof routes failed on this point: the commutator, explicit-formula Rayleigh–Ritz and scaling-limit approaches (`structure/ladder3/`). A theorem would first need a regularised ladder.
- **What it does not give.** It does not prove the ordering $\lambda_0 < \lambda_0^{\mathrm{odd}}$. That would need the Rayleigh quotients of the ladder to be ordered, and to precision $\lambda_0$ (A.3).

### A.2 Negative results

1. **No Sturm / Gantmacher–Krein oscillation.** For an oscillation kernel the $k$-th eigenfunction has exactly $k$ sign changes. That would give simplicity and alternation for all $L$. Here it fails from $k = 5$ (`osc.py`, `osc2.py`). At $L = 1$ the sign-change counts for $k = 0, \dots, 12$ are $0, 1, 2, 3, 4, 7, 8, 9, 10, 15, 16, 17, 22$. The extra zeros lie in the bulk at amplitude about $10^{-3}$ of the maximum, and are identical at $D = 92$ and $D = 140$. A second, independent computation (an internal-review agent, $D = 60, 80, 90$) found the same count at $k = 5$ for $L = 0.8$ and $1.0$. Only parity survives: count $\equiv k \pmod 2$. Counts for $k \le 4$ are exactly $k$ at every $L$ tested.
2. **The ground state does not converge to Riemann's $\Phi$.** $\Phi$ is the classical kernel with $\Xi(t) = \int \Phi(u) e^{itu}\,du$. The $L^2$ overlap between $k_L$ and $\Phi$ on $[-L, L]$ is $0.958, 0.922, 0.910, 0.903, 0.897$ at $L = 0.3, 0.6, 0.8, 1.0, 1.2$. It decreases with $L$, and the peaks differ ($\approx 1.60$ against $2.23$). $\hat k_L$ shares the zeros of $\Xi$ but not its envelope.
3. **Folding.** The Beurling–Deny route is closed for $L > 0.1406$ [PL294]. The cheap substitute, comparing $Q(\lvert\phi\rvert)$ with $Q(\phi)$ for the odd ground state $\phi$, has the wrong sign at $L = 0.3$ and $0.7$.
4. **Boundary-trace criteria degenerate.** A rank-one intertwining $Q_{\mathrm{odd}}\,\partial - \partial\, Q_{\mathrm{even}}$ gives a Weyl function whose residues $c_k = \psi_k(L)\int\psi_k$ were positive across the staircase. Their sign also changed exactly where alternation breaks ($L = 0.8$, $\lambda \approx 0.85$). But the $c_k$ drift like $(\log D)^{-1/2}$ as $D \to \infty$, because point evaluation at $\pm L$ is not continuous for this log-Laplacian-type form. In the continuum, the criterion needs a renormalised trace that has not been constructed.
5. **Pole splitting is known.** The pole kernel $2\cosh((x - y)/2)$ is rank two: $+$ on the even sector, $-$ on the odd. The pole-free form has a simple, positive, even ground state for all $L$ by Perron–Frobenius. Both facts are in [An1] (with Connes–Consani, arXiv 2106.01715, Lemma 2.7). We do not claim them.

### A.3 The common obstruction: the knife-edge margin

Every all-$L$ reformulation we tried converts the gap into a scalar inequality. These were: Krein / Aronszajn secular equations in the even or odd sector, boundary-trace Weyl functions, Hadamard-type eigenvalue flow in $L$, and de Branges / explicit-formula comparisons. The gap is large in ratio: $\min(\lambda_1^{\mathrm{even}}, \lambda_0^{\mathrm{odd}})/\lambda_0 \approx 2000$ at $L = 1$. But the scalar inequality's margin shrinks as fast as $\lambda_0$ itself.

- **Example: the odd-sector secular function** $F = 2\langle s, (M_o - \lambda_0)^{-1} s\rangle$ with $s = \sinh(x/2)$. Here $M_o$ is the odd block of the pole-free form, and "$\lambda_0$ is below the odd spectrum" is equivalent to $F < 1$. We find $F = 0.0354$ at $L = 0.3$ and $0.99997511$ at $L = 0.6$. At $L = 1$ it is $1.0000000$ to 8 digits, and a first-order estimate gives $1 - F \approx 5\cdot10^{-25}$.
- **Example: the eigenvalue flow.** Hadamard-type formulas tie $d\lambda/dL$ to an edge coefficient $A$ with $2A^2 \approx r\lambda_0$, where $r \approx 15$–$150$. A virial identity writes this as an $O(1)$ archimedean term cancelling an $O(1)$ prime term to precision $\lambda_0$.

**Conclusion of the search.** Ruling out a level crossing at window $L$ requires controlling the prime terms to relative precision about $\lambda_0(L)$, which falls doubly exponentially. So soft structure cannot decide parity for all $L$: positivity, oscillation, interlacing, monotone flows. This suggests that an all-$L$ proof of simple-even would need arithmetic input comparable to Weil positivity itself. That is a heuristic observation, not a theorem.

### A.4 Open, testable sub-claims

- **Positivity.** $k_L > 0$ on $(-L, L)$ for all $L$. True numerically at every $L$ we computed; listed as open in [An1]. On its own it does not order the odd sector, because $Q$ is not positivity-preserving.
- **Ladder limit.** The overlaps in A.1 tend to 1 as $L \to \infty$ for each fixed $j$.
- **Odd-sector orthogonality.** The odd near-kernel modes of the pole-free form are nearly orthogonal to $\sinh(x/2)$. At $L = 1$, $\lvert\langle s, q_0\rangle\rvert \approx 3\cdot10^{-20}$.

Scripts (in `exploratory/`; non-rigorous): `osc.py`, `osc2.py`, `ladder.py`; review-agent scripts `tp/phi.py`, `tp/osc.py`, `tp/fold.py`, `ref_zop/chk.py`. The three proof attempts cited in A.1 (`structure/ladder3/`) are internal notes and are not included.

---

## 10. Review log

All reviewers below were AI agents (Claude) run by the author; none is an independent human referee. All rounds predate the sharp-cutoff path (§4.6), the windows $1.03$–$1.19$ and the triangular head assembly used for $1.12$–$1.19$ (Erratum, further issues, item 4); none of these has been through any review round.

**Verification method, round 1.** Each issue was checked against the code, logs and JSON before it was applied. We also ran three light computations:

- `lipschitz.py` at 53, 128 and 768 bits;
- step ratios and bracket widths recomputed from `theorem_bounds.json`;
- the Gauss–Legendre correction $(\rho^2-1)Nq$ per window, from the code's own formula (no matrices rebuilt).

**Reviewer 1**

| # | Issue | Handling |
|---|---|---|
| R1.1 | Gauss–Legendre constant off by one node | **Fixed, confirmed real bug.** Verified against [Tr]. In round 1 this was documented as a Known issue with a Weyl-shift repair. Superseded in round 2 by the code fix and `gl_correction.py` (see C.1). |
| R1.2 | Piece I constants do not match `lipschitz.py` | **Fixed (partly valid).** Both reviewers were right in their own setting: 53-bit standalone against 768-bit inside `theorem.py`. The note quotes the 768-bit constants and records the unexplained precision dependence (§4.1, §7.5). |
| R1.3 | Observation 3 false | **Fixed in round 1 for the old $U$:** worst ratio 3.761 at $L = 0.8$. Observations recomputed again from the corrected table in round 2 (C.3), and the argument revised in round 3 (D.4). |
| R1.4 | Bounds rounded up instead of down | **Fixed.** Step ratios and gaps are computed as true lower ends and floored (again in round 2). |
| R1.5 | Bessel enclosure gap and sharpness claim | **Fixed.** Downward-induction proof; top-order looseness stated when $M = n_{\max}$. |
| R1.6 | Continuation of $\mathrm{Re}\,\psi$ off the real axis | **Fixed.** Definition via $w_\pm$ and the bound on both arguments. |
| R1.7 | Python `min` on balls | **Fixed (documented)** in §4.2, §4.3 and §7.6; also applies to `gl_correction.py`. |
| R1.8 | FLINT versions; stale docstring; base `CERT_PREC` | **Fixed.** FLINT 3.6.0 (Python) and 3.4 (C tools); stale docstrings in §7.9; unrecorded `CERT_PREC` noted in §4.5. |
| R1.9 | Fubini for (P) with Lipschitz $g$ | **Fixed.** Mollification argument (§2.2). |
| R1.10 | Self-contained discreteness proof | **Fixed.** Compact-resolvent argument (§2.3). |

**Reviewer 2**

| # | Issue | Handling |
|---|---|---|
| R2.1 | Gauss–Legendre off-by-one | **Fixed** as R1.1. The claim that $L = 1$ "may no longer pass" was rejected; confirmed in round 2: corrected $\lambda_0 \ge 5.147\cdot10^{-30}$ and all steps pass. |
| R2.2 | Remote code unchecked; files modified after runs | **Fixed in code (round 2):** SHA-256 check on every host, validated at $L = 0.4$. It post-dates the original cluster runs; §7.2, the abstract and §3 say so (D.2). |
| R2.3 | Zhu Thm 1.2 assertion proves only $8.9\cdot10^{-18} - e_B$ | **Fixed.** The cross-check cites $\lambda_0^{\mathrm{lo}}(0.8) \ge 9.787\cdot10^{-18}$. |
| R2.4 | Observation 3 false | **Fixed** as R1.3 / C.3 / D.4. |
| R2.5 | Trial degree misstated | **Fixed;** §4.4 now gives the trial degrees actually used (round 3, D.9). |
| R2.6 | Python `min`/`max` on balls | **Fixed** as R1.7. |
| R2.7 | `sph_bessel_hybrid` untested; "prototype" docstring | **Fixed (documented)** in §6 and §7.9. |
| R2.8 | Stale docstrings | **Fixed (documented)** in §7.9. |
| R2.9 | Observation 1 ratio at $L = 1$ | **Fixed (round 3, D.10).** Now quoted as the certified lower bound $\ell_{\mathrm{odd}}/U \ge 2.0\cdot10^3$, not as a ratio of two lower bounds. |

**Reviewer 3**

| # | Issue | Handling |
|---|---|---|
| R3.1 | Observation 3 false | **Fixed** as R1.3 / C.3 / D.4. |
| R3.2 | Frequency-cutoff strategy not credited to Zhu | **Fixed in rounds 2 and 3 (C.4, D.7).** |
| R3.3 | Inconsistent about what Zhu certified | **Fixed in round 2 (C.5)** from the source text, and completed in round 3 with Zhu's certified upper bounds (D.1). The round-1 weakening ("rigor unconfirmed") is withdrawn. |
| R3.4 | Unrefereed and AI-generated sources not labelled; no AI disclosure | **Fixed;** labels made consistent in rounds 2 and 3 (C.6, D.5). |
| R3.5 | Abstract and theorem state the result without the caveats | **Fixed;** reproducibility status made explicit in round 2 (C.7) and extended in round 3 (D.2). |
| R3.6 | Rounding not outward | **Fixed** as R1.4; 53-bit Lipschitz values corrected in round 3 (D.3). |
| R3.7 | "Analytic" Piece I overstated | **Fixed.** |
| R3.8 | "Perron–Frobenius route fails" overgeneralizes | **Fixed.** Now refers to the Beurling–Deny criterion [PL294]. |
| R3.9 | PL-287 missing from "not claimed" | **Fixed.** |
| R3.10 | Elston and Elder missing | **Fixed.** |
| R3.11 | "Independent reproduction" too strong | **Fixed.** |
| R3.12 | Hybrid Bessel rigor vs "prototype" status | **Fixed (documented).** |
| R3.13 | Stale README and REDUCTION | **Not fixed in repo.** Documented in §7.9 and §7.11; the note supersedes them. |
| R3.14 | Bibliography looks invented | **Fixed partly.** Status labels added; exact titles marked unverified rather than invented. |
| R3.15 | No conclusion section | **Fixed.** §8. |

**Round 2 (correction and final revision)**

*Verification method, round 2.* We read `certify.py` (fix at the panel and Rayleigh rules; SHA-256 check in the cluster dispatcher), `gl_correction.py` and `run_gl_correction.log`. We recomputed every Table 1 entry, step ratio and Piece I gap from `theorem_bounds_corrected.json` with explicit endpoint comparisons and outward rounding, and checked Zhu's Thm 1.1 and Thm 6.2 against a text extraction of the source.

| # | Item | Handling |
|---|---|---|
| C.1 | Known issue → Erratum and correction | **Done.** The bug is fixed in `certify.py` (`rho**(-2*(NG-1))`, `rho**(-2*(n_out-1))`). `gl_correction.py` subtracts a per-window Weyl shift $\delta$ from every lower bound ($\delta \le 6.4\cdot10^{-35}$ / $1.6\cdot10^{-33}$ / $3.2\cdot10^{-32}$ by group; per-window values in Table 1). It also recomputes every $U$ with the fixed `rayleigh_upper` and fresh Galerkin trial vectors. Result: "THEOREM CERTIFIED (corrected)". The round-1 † annotations are removed, since all values are now corrected. The rigor qualification was added in round 3 (D.2). |
| C.2 | Table 1 and Piece I | **Replaced** by the corrected values. Piece I: $C_{\mathrm{even}} \le 4.876362$, $C_{\mathrm{odd}} \le 0.761666$; gaps $\ge 0.4173$ / $0.6625$ (exact lower ends 0.417374 / 0.662537). The log's display "0.4174" is a midpoint; we floor the lower end. Step ratios 1.288 … 380.0, floored lower ends. These are slightly below the log's midpoints (e.g. 1.289, 380.1). |
| C.3 | Observation 3 (false with old $U$) | **Recomputed from the corrected table.** Widest bracket $1.715\times$ at $L = 0.8$, then $1.552\times$ (0.76) and $1.430\times$ (0.72); all others $\le 1.33\times$. The supporting argument was replaced in round 3 (D.4). Observations 1, 2 and 4 re-checked against the corrected numbers. |
| C.4 | Credit Zhu for the frequency-cutoff reduction | **Done.** Lemma (R) is stated as Zhu's one-stroke reduction [Zhu, Thm 1.1], and $\beta$ and the $2\pi e^{A_L}$ cutoff scale are attributed. Also in the abstract, §1.3, §1.4 and the "not claimed" lists. Extended in round 3 (D.7). |
| C.5 | What Zhu certified at $L = 0.8$ | **Made consistent with [Zhu, Thm 6.2]:** $8.9\cdot10^{-18} \le \lambda_0 \le 2.523\cdot10^{-16}$, $\lambda_1^{\mathrm{even}} \ge 2.085\cdot10^{-12}$ (rigorous, Courant–Fischer Cholesky at $2.0857\cdot10^{-12}$), $8.206\cdot10^{-15} \le \lambda_0^{\mathrm{odd}} \le 2.347\cdot10^{-14}$. The index shift relative to our notation is stated. "Reported" was replaced by "certified" throughout. The comparison table in §6 was updated, with a note on Zhu's non-rigorous $2.3175\cdot10^{-12}$ reference value. |
| C.6 | Consistent labels for unrefereed / AI-generated sources | **Done** in §1.2, §1.3 and §9: [JPB] "AI-audited"; [Eu] self-described "candidate"; all Zenodo deposits unrefereed. The round-2 label "AI research log" for mathia was withdrawn in round 3 (D.5). |
| C.7 | Honest reproducibility status in abstract and theorem | **Done.** The abstract, §3, §6 and §8 state that the windows were certified by successive code revisions, that the certificate is the original runs plus `gl_correction.py` with the fixed code, and that a full end-to-end rerun with the final code has not been done and is recommended. Extended in round 3 (D.2). |
| C.8 | SHA-256 check of remote code | **Documented** in §4.5, §6 and §7.2. It covers future cluster runs and was validated at $L = 0.4$ (no log in the repository). It does not retroactively cover the original runs for $0.84 \le L \le 1$. |
| C.9 | Side observations | The new $U$ at $L = 0.99$ is slightly larger than the old one; the certificate uses the new value uniformly. `gl_correction.py` assumes the final code reproduces the original window parameters; this is stated in §7.2. |

**Round 3 (checker issues)**

*Verification method, round 3.* We read `gl_correction.py`, `NOVELTY.md`, `run_step3*.log`, the `lipschitz.py` output, the file timestamps and the `cache/` filenames. We recomputed separations, radii, $\ell_{\mathrm{odd}}/U$ and $U/\lambda_0^{\mathrm{lo}}$ from `theorem_bounds_corrected.json` at 768 bits, evaluated the final code's window parameters for $0.84 \le L \le 1$, and checked Zhu §1.3–1.4, Lemma 3.1, Thm 1.1, §4, §5.5, §6.1–6.2 and Table 3 against the text extraction.

| # | Item | Handling |
|---|---|---|
| D.1 | Zhu's certified upper bounds omitted | **Valid; fixed.** [Zhu, §1.3, §11.4, Table 3] certifies upper bounds on $\lambda_0$ for $L \in [0.5, 2]$, and [Zhu, §1.4] gives the enclosure $8.9\cdot10^{-18} \le \lambda_0(0.8) \le 2.27\cdot10^{-17}$. These are now cited in §1.3. The by-product claim is requalified in §1.4, the bounds are added to the "not claimed" lists in §1.4 and §8, and a Table 3 row and comparison are added to §6 cross-check 2. `NOVELTY.md`'s contrary statement is flagged in §7.9. |
| D.2 | Rigor caveats missing from abstract and theorem | **Valid; fixed.** The abstract, §3, the Erratum and §8 now state that the cluster windows' remote code was unverified and that the $\delta$ correction assumes unchanged window parameters. "Rigorously corrected" is replaced by "corrected, under the stated assumption". Partial check: the final code's $N$ at $L \ge 0.98$ matches the cache filenames of the original run (§6, cross-check 8). |
| D.3 | 53-bit / 128-bit Lipschitz values were midpoints | **Valid; fixed.** These are now upper ends rounded up: 0.592773, 4.87484 and 0.760130 (53 bits), and 0.593583 (128 bits). |
| D.4 | Observation 3 argument circular | **Valid; fixed.** $U$ equals the Galerkin value by construction. We now use Zhu's non-certified window floor $1.656\cdot10^{-17}$ [Zhu, §5.5(b)], which puts $U$ within about 1.4% of it. The cause remains labelled unverified. |
| D.5 | "AI research log" label for mathia unsupported | **Valid; fixed.** `NOVELTY.md` says "analytic, unrefereed notes", and no file in the repository says the notes are AI-generated. The label is now "[GitHub notes; unrefereed]", and §7.7 says their AI status is unknown. The round-2 review instruction that attributed the label to `NOVELTY.md` misreported it. |
| D.6 | Abstract overstated [An1] | **Valid; fixed.** [An1]'s full-form statement is now described as a remark; its theorem covers only the pole-removed form. |
| D.7 | Lemma (S) presented as going beyond Zhu | **Valid; fixed.** (S) is now stated as the block estimate of [Zhu, Thm 1.1], used for the $k$-th eigenvalue as in [Zhu, §6.2]. (E) is credited to [Zhu, Lemma 3.1]. Only our re-derivations and the implementation are listed as ours. |
| D.8 | Revision timeline incomplete | **Valid; fixed.** §7.2 now includes `bessel_rec.py` (19:30, used only by the $L \ge 0.98$ run), `veig.c`/`veig` (19:31, estimates only) and `theorem.py` (20:57, gluing only). |
| D.9 | Trial-vector caps never bind | **Valid; fixed.** $D_G = 62$–$92$, at most 46 coefficients, trial degree $\le 90$; the truncation is a no-op. |
| D.10 | Observation 1 used a ratio of lower bounds | **Valid; fixed.** It now quotes the certified $\ell_{\mathrm{odd}}/U$: about 2.77 at $L = 0.1$ and at least $2.0\cdot10^3$ at $L = 1$. |
| D.11 | §4.2 separation stated as absolute 0.039 | **Valid; fixed.** It is now a relative separation of at least 1.021 ($\ge 1.09$ from $L = 0.1$), radii below $5\cdot10^{-20}$, and absolute separation more than $10^{17}$ times the radii. |
| D.12 | Galerkin cross-check sources unspecified | **Valid; fixed.** The degree and log are named for each value. The new $U$ is stated to equal the Galerkin value at its own degree at 0.5, 0.6, 0.8, 0.9 and 1.0 by construction, and the values are noted as not converged. |
| D.13 | Erratum step 2 bullets redundant; old $U$ at 0.99 | **Fixed.** The bullets are merged, and the old $U$ is quoted as $3.551\cdot10^{-29}$ (log midpoint). |