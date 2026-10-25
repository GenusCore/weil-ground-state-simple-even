import SpectralCertificate.Data.W005
import SpectralCertificate.Data.W01

/-!
# Link to the stored table

`data/theorem_bounds.json` stores, per window, balls for `λ₀`, `λ₁^even` and `λ₀^odd`. Below are the lower ends
of the `λ₁^even` and `λ₀^odd` balls at `L = 0.05` and `0.1`, as exact integers over `S = 10^59`
(`test_certificate.py` checks them against the JSON file). The two theorems show that the bounds Lean derives
from the certificates (`window_sound`) imply these stored lower ends.
-/

namespace SpectralCertificate

def S : ℕ := 10 ^ 59

def l1_005 : ℤ := 188755802864969395577540000000000000000000000000000000000000
def od_005 : ℤ := 192697293946155983209711000000000000000000000000000000000000
def l1_01 : ℤ := 138723368371374940785110000000000000000000000000000000000000
def od_01 : ℤ := 126645412793270786585600000000000000000000000000000000000000

/-- The Lean-checked window `L̃ ≥ 0.05` implies the stored lower ends of `λ₁^even, λ₀^odd`. -/
theorem W005_rows {lam0 lam1 lamOdd : ℝ} (h : HeadHyp W005.window lam0 lam1 lamOdd) :
    (l1_005 : ℝ) / S ≤ lam1 ∧ (od_005 : ℝ) / S ≤ lamOdd := by
  obtain ⟨-, h1, h2⟩ := window_sound _ W005.even_ok W005.deflated_ok W005.odd_ok h
  refine ⟨le_trans ?_ h1, le_trans ?_ h2⟩ <;>
    simp only [bound, W005.window, l1_005, od_005, S, le_sub_iff_add_le, le_min_iff] <;> norm_num

/-- Same at `L̃ ≥ 0.1`. -/
theorem W01_rows {lam0 lam1 lamOdd : ℝ} (h : HeadHyp W01.window lam0 lam1 lamOdd) :
    (l1_01 : ℝ) / S ≤ lam1 ∧ (od_01 : ℝ) / S ≤ lamOdd := by
  obtain ⟨-, h1, h2⟩ := window_sound _ W01.even_ok W01.deflated_ok W01.odd_ok h
  refine ⟨le_trans ?_ h1, le_trans ?_ h2⟩ <;>
    simp only [bound, W01.window, l1_01, od_01, S, le_sub_iff_add_le, le_min_iff] <;> norm_num

end SpectralCertificate
