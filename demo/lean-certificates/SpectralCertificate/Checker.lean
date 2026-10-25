import Mathlib

/-!
# A positive-semidefiniteness checker for interval matrices, proved sound once for all `n`

Data (all integers, so the checker is pure `Int`/`Nat` arithmetic that the Lean kernel can evaluate):

* `C` : lower-triangular rows of the center matrix, entry `(i, j)` (`j ≤ i`) means `C i j / 2^64`;
  the upper triangle is read by symmetry (`entry`).
* `rad` : one radius, `rad / 2^64`, for every entry.
* `Rc` : rows of a factor `R`, entries `/ 2^32` (an untrusted, approximate Cholesky factor). Row `i` is
  stored as `(s, vals)`, meaning `s` zeros followed by `vals` (`expand`).
* `v` : optional deflation vector, entries `/ 2^32` (`[]` for none).
* `M`, `K`, `mu` : shift, entrywise residual bound, and the certified lower bound, all `/ 2^64`.

`check` accepts iff `mu + n K ≤ M`, `C` has `n` rows, every row of `R` has length `≤ n`, and for every
`j ≤ i < n`: `|C i j + v i v j - M δ_ij - (R Rᵀ) i j| + rad ≤ K`.
`check_sound`: then `x ᵀ A x + (v ⬝ x)² ≥ mu ‖x‖²` for every real `A` in the box.
-/

namespace SpectralCertificate

/-- Dot product of two integer lists (stops at the shorter one). -/
def dot : List ℤ → List ℤ → ℤ
  | a :: as, b :: bs => a * b + dot as bs
  | [], _ => 0
  | _ :: _, [] => 0

/-- Symmetric lookup in a lower-triangular list matrix. -/
def entry (C : List (List ℤ)) (i j : ℕ) : ℤ := (C.getD (max i j) []).getD (min i j) 0

/-- A row stored as `(s, vals)`: `s` zeros, then `vals`. -/
def expand (r : ℕ × List ℤ) : List ℤ := List.replicate r.1 0 ++ r.2

/-- `dot` of two expanded rows, computed without the leading zeros. -/
def dotC (a b : ℕ × List ℤ) : ℤ :=
  if a.1 ≤ b.1 then dot (a.2.drop (b.1 - a.1)) b.2 else dot a.2 (b.2.drop (a.1 - b.1))

structure Cert where
  n : ℕ
  C : List (List ℤ)
  Rc : List (ℕ × List ℤ)
  v : List ℤ
  rad : ℕ
  M : ℤ
  K : ℕ
  mu : ℤ

/-- The factor as a dense list matrix. -/
def Cert.R (c : Cert) : List (List ℤ) := c.Rc.map expand

/-- Integer residual `2^64 (C + v vᵀ - M I - R Rᵀ)_{ij}` (the specification). -/
def resid (c : Cert) (i j : ℕ) : ℤ :=
  entry c.C i j + c.v.getD i 0 * c.v.getD j 0 - (if i = j then c.M else 0) -
    dot (c.R.getD i []) (c.R.getD j [])

/-- Row `i` against the previous rows `j < i` (structural recursion: no indexing in the kernel). -/
def rowCheck (rad K : ℕ) (M : ℤ) (Ri : ℕ × List ℤ) (vi : ℤ) : List ℤ → List (ℕ × List ℤ) → List ℤ → Bool
  | [], _, _ => false
  | [c], [], _ => decide ((c + vi * vi - M - dotC Ri Ri).natAbs + rad ≤ K)
  | _ :: _ :: _, [], _ => false
  | c :: cs, Rj :: Rs, vs =>
      decide ((c + vi * vs.headD 0 - dotC Ri Rj).natAbs + rad ≤ K) && rowCheck rad K M Ri vi cs Rs vs.tail

def rowsCheck (rad K : ℕ) (M : ℤ) (v : List ℤ) :
    List (List ℤ) → List (ℕ × List ℤ) → List (ℕ × List ℤ) → List ℤ → Bool
  | Ci :: Cs, Ri :: Rs, prev, vs =>
      rowCheck rad K M Ri (vs.headD 0) Ci prev v && rowsCheck rad K M v Cs Rs (prev ++ [Ri]) vs.tail
  | [], [], _, _ => true
  | [], _ :: _, _, _ => false
  | _ :: _, [], _, _ => false

/-- The checker. Evaluated by the kernel (`decide +kernel`), never by compiled code. -/
def check (c : Cert) : Bool :=
  decide (c.mu + c.n * c.K ≤ c.M) && decide (c.C.length = c.n) &&
    c.Rc.all (fun r => decide (r.1 + r.2.length ≤ c.n)) && rowsCheck c.rad c.K c.M c.v c.C c.Rc [] c.v

/-- A real matrix lies in the box of center `C` and radius `rad`. -/
def InBox {n : ℕ} (C : List (List ℤ)) (rad : ℕ) (A : Fin n → Fin n → ℝ) : Prop :=
  ∀ i j, |A i j - (entry C i j : ℝ) / 2 ^ 64| ≤ (rad : ℝ) / 2 ^ 64

/-- The quadratic form `x ᵀ A x`. -/
def quad {n : ℕ} (A : Fin n → Fin n → ℝ) (x : Fin n → ℝ) : ℝ := ∑ i, ∑ j, x i * A i j * x j

/-- The real vector encoded by an integer list at scale `2^32`. -/
noncomputable def vec (n : ℕ) (v : List ℤ) : Fin n → ℝ := fun i => (v.getD i 0 : ℝ) / 2 ^ 32

/-! ### Pure real algebra -/

theorem quad_lower {n : ℕ} (A R : Fin n → Fin n → ℝ) (v x : Fin n → ℝ) (m k : ℝ) (hk : 0 ≤ k)
    (hE : ∀ i j, |A i j + v i * v j - (if i = j then m else 0) - ∑ l, R i l * R j l| ≤ k) :
    (m - n * k) * ∑ i, x i ^ 2 ≤ quad A x + (∑ i, v i * x i) ^ 2 := by
  set D : Fin n → Fin n → ℝ :=
    fun i j => A i j + v i * v j - (if i = j then m else 0) - ∑ l, R i l * R j l with hD
  have expand : quad A x + (∑ i, v i * x i) ^ 2 =
      m * ∑ i, x i ^ 2 + ∑ l, (∑ i, R i l * x i) ^ 2 + ∑ i, ∑ j, x i * D i j * x j := by
    have h1 : (∑ i, v i * x i) ^ 2 = ∑ i, ∑ j, x i * (v i * v j) * x j := by
      rw [sq, Finset.sum_mul_sum]; refine Finset.sum_congr rfl fun i _ => ?_
      refine Finset.sum_congr rfl fun j _ => ?_; ring
    have h2 : m * ∑ i, x i ^ 2 = ∑ i, ∑ j, x i * (if i = j then m else 0) * x j := by
      rw [Finset.mul_sum]; refine Finset.sum_congr rfl fun i _ => ?_
      simp only [mul_ite, mul_zero, ite_mul, zero_mul, Finset.sum_ite_eq, Finset.mem_univ, if_true]
      ring
    have h3 : ∑ l, (∑ i, R i l * x i) ^ 2 = ∑ i, ∑ j, x i * (∑ l, R i l * R j l) * x j := by
      simp only [sq, Finset.sum_mul_sum, Finset.mul_sum, Finset.sum_mul]
      rw [Finset.sum_comm]; refine Finset.sum_congr rfl fun i _ => ?_
      rw [Finset.sum_comm]; refine Finset.sum_congr rfl fun j _ => ?_
      refine Finset.sum_congr rfl fun l _ => ?_; ring
    rw [h1, h2, h3, quad, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun j _ => ?_
    simp only [hD]; ring
  have hcross : -(n * k * ∑ i, x i ^ 2) ≤ ∑ i, ∑ j, x i * D i j * x j := by
    have habs : -(k * (∑ i, |x i|) ^ 2) ≤ ∑ i, ∑ j, x i * D i j * x j := by
      rw [sq, Finset.sum_mul_sum, Finset.mul_sum, ← Finset.sum_neg_distrib]
      refine Finset.sum_le_sum fun i _ => ?_
      rw [Finset.mul_sum, ← Finset.sum_neg_distrib]
      refine Finset.sum_le_sum fun j _ => ?_
      have := hE i j
      have hb : |x i * D i j * x j| ≤ k * (|x i| * |x j|) := by
        rw [abs_mul, abs_mul]
        calc |x i| * |D i j| * |x j| = |D i j| * (|x i| * |x j|) := by ring
          _ ≤ k * (|x i| * |x j|) := mul_le_mul_of_nonneg_right this (by positivity)
      linarith [neg_abs_le (x i * D i j * x j)]
    have hcs : (∑ i, |x i|) ^ 2 ≤ n * ∑ i, x i ^ 2 := by
      have := sq_sum_le_card_mul_sum_sq (s := Finset.univ) (f := fun i : Fin n => |x i|)
      simpa [sq_abs] using this
    nlinarith [mul_le_mul_of_nonneg_left hcs hk]
  have hR : 0 ≤ ∑ l, (∑ i, R i l * x i) ^ 2 := Finset.sum_nonneg fun l _ => sq_nonneg _
  rw [expand]; nlinarith

/-! ### From the integer checker to the real hypothesis -/

theorem dot_comm : ∀ a b : List ℤ, dot a b = dot b a
  | _ :: _, _ :: _ => by simp only [dot]; rw [dot_comm, mul_comm]
  | [], [] | [], _ :: _ | _ :: _, [] => by simp [dot]

theorem dot_eq_sum : ∀ (a b : List ℤ) (N : ℕ), a.length ≤ N →
    dot a b = ∑ k ∈ Finset.range N, a.getD k 0 * b.getD k 0
  | [], b, N, _ => by simp [dot]
  | _ :: _, [], N, _ => by simp [dot]
  | x :: a, y :: b, N, h => by
    obtain ⟨N, rfl⟩ : ∃ N', N = N' + 1 := ⟨N - 1, by simp at h; omega⟩
    rw [Finset.sum_range_succ', dot, dot_eq_sum a b N (by simpa using h)]
    simp [add_comm]

theorem dot_rep_right : ∀ (d : ℕ) (x y : List ℤ), dot x (List.replicate d 0 ++ y) = dot (x.drop d) y
  | 0, x, y => by simp
  | _ + 1, [], y => by simp [dot]
  | d + 1, a :: x, y => by simp [List.replicate_succ, dot, dot_rep_right d x y]

theorem dot_rep_both : ∀ (k : ℕ) (x y : List ℤ),
    dot (List.replicate k 0 ++ x) (List.replicate k 0 ++ y) = dot x y
  | 0, x, y => by simp
  | k + 1, x, y => by simp [List.replicate_succ, dot, dot_rep_both k x y]

theorem dotC_eq (a b : ℕ × List ℤ) : dotC a b = dot (expand a) (expand b) := by
  unfold dotC expand
  split_ifs with h
  · obtain ⟨d, hd⟩ := Nat.exists_eq_add_of_le h
    rw [hd, List.replicate_add, List.append_assoc, dot_rep_both, dot_rep_right]; congr 2; omega
  · obtain ⟨d, hd⟩ := Nat.exists_eq_add_of_le (le_of_lt (not_le.mp h))
    rw [hd, List.replicate_add, List.append_assoc, dot_rep_both, dot_comm (List.replicate d 0 ++ a.2),
      dot_rep_right, dot_comm]
    congr 2; omega

theorem getD_tail (vs : List ℤ) (j : ℕ) : vs.tail.getD j 0 = vs.getD (j + 1) 0 := by
  cases vs <;> simp

theorem headD_eq (vs : List ℤ) : vs.headD 0 = vs.getD 0 0 := by
  cases vs <;> simp

theorem getD_take {α : Type*} (d : α) : ∀ (l : List α) (i j : ℕ), j < i → (l.take i).getD j d = l.getD j d
  | [], _, _, _ => by simp
  | _ :: _, 0, _, h => absurd h (Nat.not_lt_zero _)
  | _ :: _, _ + 1, 0, _ => by simp
  | _ :: l, i + 1, j + 1, h => by simpa using getD_take d l i j (by omega)

theorem rowCheck_spec (rad K : ℕ) (M : ℤ) (Ri : ℕ × List ℤ) (vi : ℤ) :
    ∀ (prev : List (ℕ × List ℤ)) (Ci vs : List ℤ), rowCheck rad K M Ri vi Ci prev vs = true →
      Ci.length = prev.length + 1 ∧
      (∀ j < prev.length, (Ci.getD j 0 + vi * vs.getD j 0 - dotC Ri (prev.getD j (0, []))).natAbs + rad ≤ K) ∧
      (Ci.getD prev.length 0 + vi * vi - M - dotC Ri Ri).natAbs + rad ≤ K
  | [], [], _, h => by simp [rowCheck] at h
  | [], [c], _, h => by simpa [rowCheck] using h
  | [], _ :: _ :: _, _, h => by simp [rowCheck] at h
  | _ :: _, [], _, h => by simp [rowCheck] at h
  | Rj :: Rs, c :: cs, vs, h => by
    simp only [rowCheck, Bool.and_eq_true, decide_eq_true_eq] at h
    obtain ⟨hlen, hoff, hdiag⟩ := rowCheck_spec rad K M Ri vi Rs cs vs.tail h.2
    refine ⟨by simp [hlen], fun j hj => ?_, by simpa using hdiag⟩
    rcases j with _ | j
    · have h1 := h.1; rw [headD_eq] at h1; simpa using h1
    · have := hoff j (by simpa using hj)
      simpa [getD_tail] using this

theorem rowsCheck_spec (rad K : ℕ) (M : ℤ) (v : List ℤ) :
    ∀ (Cs : List (List ℤ)) (Rs prev : List (ℕ × List ℤ)) (vs : List ℤ),
      rowsCheck rad K M v Cs Rs prev vs = true →
      Cs.length = Rs.length ∧ ∀ i < Cs.length,
        rowCheck rad K M (Rs.getD i (0, [])) (vs.getD i 0) (Cs.getD i []) (prev ++ Rs.take i) v = true
  | [], [], _, _, _ => by simp
  | [], _ :: _, _, _, h => by simp [rowsCheck] at h
  | _ :: _, [], _, _, h => by simp [rowsCheck] at h
  | Ci :: Cs, Ri :: Rs, prev, vs, h => by
    simp only [rowsCheck, Bool.and_eq_true] at h
    obtain ⟨hlen, hrow⟩ := rowsCheck_spec rad K M v Cs Rs (prev ++ [Ri]) vs.tail h.2
    refine ⟨by simp [hlen], fun i hi => ?_⟩
    rcases i with _ | i
    · have h1 := h.1; rw [headD_eq] at h1; simpa using h1
    · have := hrow i (by simpa using hi)
      simpa [getD_tail] using this

theorem R_getD (c : Cert) (i : ℕ) : c.R.getD i [] = expand (c.Rc.getD i (0, [])) := by
  unfold Cert.R
  induction c.Rc generalizing i with
  | nil => simp [expand]
  | cons r rs ih => cases i <;> simp only [List.map_cons, List.getD_cons_zero, List.getD_cons_succ, ih]

theorem resid_symm (c : Cert) (i j : ℕ) : resid c i j = resid c j i := by
  simp only [resid, entry, max_comm i j, min_comm i j, mul_comm (c.v.getD i 0), dot_comm (c.R.getD i []),
    eq_comm (a := i)]

theorem check_entry (c : Cert) (h : check c = true) (i j : ℕ) (hi : i < c.n) (hj : j < c.n) :
    (resid c i j).natAbs + c.rad ≤ c.K := by
  simp only [check, Bool.and_eq_true, decide_eq_true_eq] at h
  obtain ⟨⟨⟨-, hC⟩, -⟩, hrows⟩ := h
  obtain ⟨hlenR, hrow⟩ := rowsCheck_spec _ _ _ _ _ _ _ _ hrows
  -- the lower triangle `j ≤ i`
  have low : ∀ i j, i < c.n → j ≤ i → (resid c i j).natAbs + c.rad ≤ c.K := by
    intro i j hi hji
    have hiC : i < c.C.length := by omega
    have hspec := rowCheck_spec _ _ _ _ _ _ _ _ (hrow i hiC)
    simp only [List.nil_append, List.length_take, show min i c.Rc.length = i by omega] at hspec
    obtain ⟨-, hoff, hdiag⟩ := hspec
    rcases Nat.lt_or_eq_of_le hji with hlt | rfl
    · have := hoff j hlt
      rw [getD_take _ _ _ _ hlt, dotC_eq, ← R_getD, ← R_getD] at this
      simpa [resid, entry, max_eq_left hji, min_eq_right hji, Nat.ne_of_gt hlt] using this
    · rw [dotC_eq, ← R_getD] at hdiag
      simpa [resid, entry] using hdiag
  rcases le_total j i with hji | hij
  · exact low i j hi hji
  · rw [resid_symm]; exact low j i hj hij

theorem check_sound (c : Cert) (h : check c = true) (A : Fin c.n → Fin c.n → ℝ) (hA : InBox c.C c.rad A)
    (x : Fin c.n → ℝ) :
    (c.mu : ℝ) / 2 ^ 64 * ∑ i, x i ^ 2 ≤ quad A x + (∑ i, vec c.n c.v i * x i) ^ 2 := by
  have hcheck := h
  simp only [check, Bool.and_eq_true, List.all_eq_true, decide_eq_true_eq] at hcheck
  obtain ⟨⟨⟨hmu, -⟩, hlen⟩, -⟩ := hcheck
  set R : Fin c.n → Fin c.n → ℝ := fun i l => ((c.R.getD i []).getD l 0 : ℝ) / 2 ^ 32
  have hdot : ∀ i j : Fin c.n, (dot (c.R.getD i []) (c.R.getD j []) : ℝ) / 2 ^ 64 = ∑ l, R i l * R j l := by
    intro i j
    have hl : (c.R.getD i []).length ≤ c.n := by
      rw [R_getD]
      rcases Nat.lt_or_ge (i : ℕ) c.Rc.length with hi | hi
      · rw [List.getD_eq_getElem _ _ hi]; simpa [expand] using hlen _ (List.getElem_mem hi)
      · rw [List.getD_eq_default _ _ hi]; simp [expand]
    rw [dot_eq_sum _ _ c.n hl, ← Fin.sum_univ_eq_sum_range]
    push_cast; rw [Finset.sum_div]
    refine Finset.sum_congr rfl fun l _ => ?_
    simp only [R]; ring
  have hk : (0 : ℝ) ≤ (c.K : ℝ) / 2 ^ 64 := by positivity
  have key := quad_lower A R (vec c.n c.v) x ((c.M : ℝ) / 2 ^ 64) ((c.K : ℝ) / 2 ^ 64) hk (by
    intro i j
    have he := check_entry c h i j i.2 j.2
    have hres : ((resid c i j : ℤ) : ℝ) / 2 ^ 64 = (entry c.C i j : ℝ) / 2 ^ 64 + vec c.n c.v i * vec c.n c.v j -
        (if i = j then (c.M : ℝ) / 2 ^ 64 else 0) - ∑ l, R i l * R j l := by
      rw [← hdot]; simp only [resid, vec, Fin.val_eq_val, apply_ite (fun z : ℤ => (z : ℝ))]
      split_ifs <;> push_cast <;> ring
    have habs : |((resid c i j : ℤ) : ℝ)| + c.rad ≤ c.K := by
      rw [← Int.cast_abs, ← Int.cast_natAbs]; exact_mod_cast he
    have hbox := hA i j
    have hsplit : A i j + vec c.n c.v i * vec c.n c.v j - (if i = j then (c.M : ℝ) / 2 ^ 64 else 0) -
        ∑ l, R i l * R j l = (A i j - (entry c.C i j : ℝ) / 2 ^ 64) + ((resid c i j : ℤ) : ℝ) / 2 ^ 64 := by
      rw [hres]; ring
    rw [hsplit]
    calc _ ≤ |A i j - (entry c.C i j : ℝ) / 2 ^ 64| + |((resid c i j : ℤ) : ℝ) / 2 ^ 64| := abs_add _ _
      _ ≤ (c.rad : ℝ) / 2 ^ 64 + |((resid c i j : ℤ) : ℝ)| / 2 ^ 64 := by
          rw [abs_div, abs_of_pos (by positivity : (0 : ℝ) < 2 ^ 64)]; linarith
      _ ≤ (c.K : ℝ) / 2 ^ 64 := by rw [div_add_div_same, div_le_div_iff_of_pos_right (by positivity)]; linarith)
  have hm : (c.mu : ℝ) / 2 ^ 64 ≤ (c.M : ℝ) / 2 ^ 64 - c.n * ((c.K : ℝ) / 2 ^ 64) := by
    have : (c.mu : ℝ) + c.n * c.K ≤ c.M := by exact_mod_cast hmu
    rw [mul_div_assoc', div_sub_div_same, div_le_div_iff_of_pos_right (by positivity)]; linarith
  have hx : 0 ≤ ∑ i, x i ^ 2 := Finset.sum_nonneg fun i _ => sq_nonneg _
  nlinarith [mul_le_mul_of_nonneg_right hm hx]

/-- Deflation form of the certificate: on the hyperplane `v ⬝ x = 0` the form is at least `mu`. -/
theorem check_sound_perp (c : Cert) (h : check c = true) (A : Fin c.n → Fin c.n → ℝ) (hA : InBox c.C c.rad A)
    (x : Fin c.n → ℝ) (hx : ∑ i, vec c.n c.v i * x i = 0) :
    (c.mu : ℝ) / 2 ^ 64 * ∑ i, x i ^ 2 ≤ quad A x := by
  simpa [hx] using check_sound c h A hA x

end SpectralCertificate
