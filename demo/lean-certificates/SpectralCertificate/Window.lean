import SpectralCertificate.Checker

/-!
# One window: three head-matrix certificates and the trusted constants

For a window `[-L, L]` the paper (`paper/NOTE.md`, §4.3) bounds three eigenvalues of Weil's form `Q` from
below through the `N × N` head matrices `A` (even sector) and `A'` (odd sector):

* `λ₀(Q) ≥ min(μ₀, β - e_D) - e_B` if `A ≥ μ₀`,                        Lemmas (R), (S);
* `λ₁^even(Q) ≥ min(μ₁, β - e_D) - e_B` if `A ≥ μ₁` on `v⊥`,             Lemmas (R), (S), (D);
* `λ₀^odd(Q) ≥ min(μ₂, β - e_D') - e_B'` if `A' ≥ μ₂`,                   Lemmas (R), (S).

`HeadHyp` collects exactly what is NOT checked in Lean: that the true head matrices lie in the exported boxes
(Arb enclosures), the bounds on `β`, `e_D`, `e_B`, and the three implications above. `window_sound` is the
Lean-checked part: the three matrix inequalities (from `check`) and the arithmetic.
-/

namespace SpectralCertificate

/-- Decode `k` integers packed into the natural number `x`, base `B`, offset `off`. -/
def unpack (B : ℕ) (off : ℤ) : ℕ → ℕ → List ℤ
  | 0, _ => []
  | k + 1, x => (((x % B : ℕ) : ℤ) - off) :: unpack B off k (x / B)

/-- Decode `f` lower-triangular rows (row `i` has `i + 1` entries). -/
def unpackTri (B : ℕ) (off : ℤ) : ℕ → ℕ → ℕ → List (List ℤ)
  | 0, _, _ => []
  | f + 1, i, x => unpack B off (i + 1) x :: unpackTri B off f (i + 1) (x / B ^ (i + 1))

/-- Decode factor rows `(s, vals)`, row `i` holding columns `s..i`. -/
def unpackRows (B : ℕ) (off : ℤ) : List ℕ → ℕ → ℕ → List (ℕ × List ℤ)
  | [], _, _ => []
  | s :: ss, i, x => (s, unpack B off (i + 1 - s) x) :: unpackRows B off ss (i + 1) (x / B ^ (i + 1 - s))

structure Factor where
  Rc : List (ℕ × List ℤ)
  M : ℤ
  K : ℕ
  mu : ℤ

structure Window where
  n : ℕ
  Ce : List (List ℤ)
  Co : List (List ℤ)
  rade : ℕ
  rado : ℕ
  v : List ℤ
  fe : Factor
  fd : Factor
  fo : Factor
  betaLo : ℚ
  eDe : ℚ
  eBe : ℚ
  eDo : ℚ
  eBo : ℚ

def Window.even (w : Window) : Cert := ⟨w.n, w.Ce, w.fe.Rc, [], w.rade, w.fe.M, w.fe.K, w.fe.mu⟩
def Window.deflated (w : Window) : Cert := ⟨w.n, w.Ce, w.fd.Rc, w.v, w.rade, w.fd.M, w.fd.K, w.fd.mu⟩
def Window.odd (w : Window) : Cert := ⟨w.n, w.Co, w.fo.Rc, [], w.rado, w.fo.M, w.fo.K, w.fo.mu⟩

/-- The trusted inputs of one window (see the module docstring). `lam0`, `lam1`, `lamOdd` stand for
`λ₀(Q)`, `λ₁^even(Q)`, `λ₀^odd(Q)` on the dyadic window `[-L̃, L̃]`. The witnesses are existential, but the
hypothesis is meant to be discharged by the true objects at `L̃`: `A`, `A'` the true head matrices,
`β, e_D, e_B` the true constants of Lemmas (E), (R), (S). -/
def HeadHyp (w : Window) (lam0 lam1 lamOdd : ℝ) : Prop :=
  ∃ (A A' : Fin w.n → Fin w.n → ℝ) (β eD eB eD' eB' : ℝ),
    InBox w.Ce w.rade A ∧ InBox w.Co w.rado A' ∧
    (w.betaLo : ℝ) ≤ β ∧ eD ≤ w.eDe ∧ eB ≤ w.eBe ∧ eD' ≤ w.eDo ∧ eB' ≤ w.eBo ∧
    (∀ μ : ℝ, (∀ x, μ * ∑ i, x i ^ 2 ≤ quad A x) → min μ (β - eD) - eB ≤ lam0) ∧
    (∀ μ : ℝ, (∀ x, ∑ i, vec w.n w.v i * x i = 0 → μ * ∑ i, x i ^ 2 ≤ quad A x) →
      min μ (β - eD) - eB ≤ lam1) ∧
    (∀ μ : ℝ, (∀ x, μ * ∑ i, x i ^ 2 ≤ quad A' x) → min μ (β - eD') - eB' ≤ lamOdd)

/-- The lower bound produced by a certified `μ` and the trusted constants. -/
noncomputable def bound (mu : ℤ) (betaLo eD eB : ℚ) : ℝ := min ((mu : ℝ) / 2 ^ 64) ((betaLo : ℝ) - eD) - (eB : ℝ)

theorem window_sound (w : Window) (he : check w.even = true) (hd : check w.deflated = true)
    (ho : check w.odd = true) {lam0 lam1 lamOdd : ℝ} (h : HeadHyp w lam0 lam1 lamOdd) :
    bound w.fe.mu w.betaLo w.eDe w.eBe ≤ lam0 ∧ bound w.fd.mu w.betaLo w.eDe w.eBe ≤ lam1 ∧
      bound w.fo.mu w.betaLo w.eDo w.eBo ≤ lamOdd := by
  obtain ⟨A, A', β, eD, eB, eD', eB', hA, hA', hβ, heD, heB, heD', heB', h0, h1, h2⟩ := h
  have mono : ∀ (μ β eD eB : ℝ) (b d e : ℚ), (b : ℝ) ≤ β → eD ≤ d → eB ≤ e →
      min μ (b - d) - e ≤ min μ (β - eD) - eB := by
    intro μ β eD eB b d e hb hd he
    have := min_le_min_left μ (show (b : ℝ) - d ≤ β - eD by linarith)
    linarith
  simp only [bound]
  refine ⟨le_trans (mono _ _ _ _ _ _ _ hβ heD heB) (h0 _ fun x => ?_),
    le_trans (mono _ _ _ _ _ _ _ hβ heD heB) (h1 _ fun x hx => ?_),
    le_trans (mono _ _ _ _ _ _ _ hβ heD' heB') (h2 _ fun x => ?_)⟩
  · simpa [Window.even, vec] using check_sound w.even he A hA x
  · exact check_sound_perp w.deflated hd A hA x hx
  · simpa [Window.odd, vec] using check_sound w.odd ho A' hA' x

end SpectralCertificate
