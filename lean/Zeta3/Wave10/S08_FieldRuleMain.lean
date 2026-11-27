import Mathlib
import Zeta3.Wave10.S03_FieldRule
/-!
# S08 -- continuation of S03: the field rule, nested-sum and `𝒢_B` statements

Target id: S08, the continuation of `Zeta3.Wave10.S03_FieldRule` (planner split, 2026-11-26).
Paper label: `conj:field` (paper/ieee/main.tex, Theorem "Field rule", lines 1486-1494; the
"More precisely" sentence, line 1493).
Supplement: paper/ieee/proofs/sec_g6_field.tex (Lemma `g6:lem:A` 73-102; Lemma `g6:prop:R` 133-137;
proof of the last sentence 184-188). Statements `stmt_4, 8, 9, 11` keep their S03 numbers; the notation
(`S03.InDB`, `S03.InG`, `S03.ST`, …) is that of S03.
-/
-- BEGIN STATEMENTS
namespace S08

open Real MeasureTheory Set Filter

/-- Lemma `g6:lem:A`: for `k ≥ 0` and a monomial `𝔷^s` in `𝔷_2, …, 𝔷_k` of weight `W ≤ k` there is
`T ∈ 𝔇^{(k+1-W)}` with `[𝔷^s] P_{k,n} = ⟨T⟩_n = 𝒮_T(2n)` for every `n ≥ 1`. -/
def stmt_4 : Prop :=
  ∀ k : ℕ, ∀ s : ℕ →₀ ℕ, (∀ w ∈ s.support, 2 ≤ w ∧ w ≤ k) → S03.mwt s ≤ k →
    ∃ T : List S03.DTerm, S03.InDB (k + 1 - S03.mwt s) T ∧
      ∀ n : ℕ, 1 ≤ n → (Q04.Pf k n).coeff s = S03.ST T (2 * n)

/-- Lemma `g6:prop:R`: for `T ∈ 𝔇^{(B)}` there is `F ∈ 𝒢_B` with `[x^N]F = 𝒮_T(N)` for every even
`N ≥ 2`: `F(x) = q(1/x) + Σ_{N≥0} a_N x^N` on `0 < |x| < 1` (a Laurent series with finitely many
negative powers) and `a_N = 𝒮_T(N)` for even `N ≥ 2`. -/
def stmt_8 : Prop :=
  ∀ B : ℕ, ∀ T : List S03.DTerm, S03.InDB B T →
    ∃ F : ℝ → ℝ, S03.InG B F ∧ ∃ q : Polynomial ℚ, ∃ a : ℕ → ℚ,
      (∀ x : ℝ, 0 < |x| → |x| < 1 →
        HasSum (fun N : ℕ => (a N : ℝ) * x ^ N) (F x - Polynomial.aeval (1 / x) q)) ∧
      ∀ N : ℕ, 2 ≤ N → Even N → a N = S03.ST T N

/-- Proof of the theorem, first paragraph: for every `k` and every monomial `𝔷^s` in
`𝔷_2, …, 𝔷_k` of weight `W ≤ k` there is `F_s ∈ 𝒢_{k+1-W}` with
`F_s(x) = Σ_{n≥1} [𝔷^s]P_{k,n} x^{2n}` for `0 < |x| < 1`. -/
def stmt_9 : Prop :=
  ∀ k : ℕ, ∀ s : ℕ →₀ ℕ, (∀ w ∈ s.support, 2 ≤ w ∧ w ≤ k) → S03.mwt s ≤ k →
    ∃ F : ℝ → ℝ, S03.InG (k + 1 - S03.mwt s) F ∧
      ∀ x : ℝ, 0 < |x| → |x| < 1 →
        HasSum (fun n : ℕ => (((Q04.Pf k (n + 1)).coeff s : ℚ) : ℝ) * x ^ (2 * (n + 1))) (F x)

/-- The "More precisely" sentence of `conj:field`: for every `k` there is a finite list of terms
`(s_j, p_j, a_j, b_j, d_j, w_j)`, independent of the angle, with `𝔷^{s_j}` a monomial in
`𝔷_2, …, 𝔷_k`, `p_j ∈ ℚ[x]`, `w_j` a regular word over `{0,1,-1}` and
`wt(𝔷^{s_j}) + |w_j| ≤ k + 1`, such that for every `0 < α < π`, with `c = cos(α/2)`,
`M_k(α) = Σ_j ζ^{s_j} · p_j(c) / (c^{a_j} (1-c)^{b_j} (1+c)^{d_j}) · H_{w_j}(c)`. -/
def stmt_11 : Prop :=
  ∀ k : ℕ, ∃ L : List S03.FTerm,
    (∀ g ∈ L, (∀ w ∈ g.1.support, 2 ≤ w ∧ w ≤ k) ∧ S03.IsLetters g.2.2.2.2.2 ∧
      S03.Regular g.2.2.2.2.2 ∧ S03.mwt g.1 + g.2.2.2.2.2.length ≤ k + 1) ∧
    ∀ α : ℝ, 0 < α → α < π →
      V05.Mk k α = (L.map fun g => S03.zmon g.1 * S03.kval g.2.1 g.2.2.1 g.2.2.2.1 g.2.2.2.2.1 (D04.c α)
        * S03.H g.2.2.2.2.2 (D04.c α)).sum

end S08
-- END STATEMENTS

-- ===== merged from work/S08/PartA.lean =====
namespace S08.PA

open S03

/-- `f n = 𝒮_T(2n)` for `n ≥ 1` with `T ∈ 𝔇^{(B)}`. -/
def Rep (B : ℕ) (f : ℕ → ℚ) : Prop :=
  ∃ T : List DTerm, InDB B T ∧ ∀ n : ℕ, 1 ≤ n → f n = S03.ST T (2 * n)

lemma ST_nil (N : ℕ) : S03.ST [] N = 0 := by simp [S03.ST]

lemma ST_cons (τ : DTerm) (T : List DTerm) (N : ℕ) :
    S03.ST (τ :: T) N = ratEval τ.1 τ.2.1 (N : ℚ) * Z τ.2.2 N + S03.ST T N := by
  simp [S03.ST]

lemma ST_append (T T' : List DTerm) (N : ℕ) : S03.ST (T ++ T') N = S03.ST T N + S03.ST T' N := by
  simp [S03.ST, List.map_append, List.sum_append]

lemma InDB_nil (B : ℕ) : InDB B [] := by intro τ h; simp at h

lemma InDB_append {B : ℕ} {T T' : List DTerm} (h : InDB B T) (h' : InDB B T') :
    InDB B (T ++ T') := by
  intro τ hτ
  rcases List.mem_append.mp hτ with h1 | h1
  · exact h τ h1
  · exact h' τ h1

lemma InDB_mono {B B' : ℕ} {T : List DTerm} (h : InDB B T) (hB : B ≤ B') : InDB B' T := by
  intro τ hτ
  obtain ⟨h1, h2, h3⟩ := h τ hτ
  exact ⟨h1, by omega, fun p hp => ⟨(h3 p hp).1, by have := (h3 p hp).2; omega⟩⟩

lemma Rep_congr {B : ℕ} {f g : ℕ → ℚ} (h : Rep B f) (hfg : ∀ n, 1 ≤ n → g n = f n) :
    Rep B g := by
  obtain ⟨T, hT, hf⟩ := h
  exact ⟨T, hT, fun n hn => by rw [hfg n hn, hf n hn]⟩

lemma Rep_mono {B B' : ℕ} {f : ℕ → ℚ} (h : Rep B f) (hB : B ≤ B') : Rep B' f := by
  obtain ⟨T, hT, hf⟩ := h
  exact ⟨T, InDB_mono hT hB, hf⟩

lemma Rep_zero (B : ℕ) : Rep B (fun _ => 0) := ⟨[], InDB_nil B, fun n _ => by simp [S03.ST]⟩

lemma Rep_zero' (B : ℕ) {f : ℕ → ℚ} (h : ∀ n, 1 ≤ n → f n = 0) : Rep B f :=
  Rep_congr (Rep_zero B) h

lemma Rep_add {B : ℕ} {f g : ℕ → ℚ} (hf : Rep B f) (hg : Rep B g) :
    Rep B (fun n => f n + g n) := by
  obtain ⟨T, hT, hf⟩ := hf
  obtain ⟨T', hT', hg⟩ := hg
  exact ⟨T ++ T', InDB_append hT hT', fun n hn => by show f n + g n = _; rw [ST_append, hf n hn, hg n hn]⟩

lemma Rep_sum {ι : Type*} (S : Finset ι) {B : ℕ} (f : ι → ℕ → ℚ)
    (h : ∀ i ∈ S, Rep B (f i)) : Rep B (fun n => ∑ i ∈ S, f i n) := by
  classical
  induction S using Finset.induction_on with
  | empty => simpa using Rep_zero B
  | insert a S ha ih =>
    have h1 := h a (Finset.mem_insert_self a S)
    have h2 := ih (fun i hi => h i (Finset.mem_insert_of_mem hi))
    refine Rep_congr (Rep_add h1 h2) (fun n _ => ?_)
    rw [Finset.sum_insert ha]

/-- A rational function with poles at integers `≤ 1` of order `≤ B`. -/
lemma Rep_rat {B : ℕ} {f : ℕ → ℚ} (Q : Polynomial ℚ) (e : ℤ →₀ ℕ)
    (he : ∀ p ∈ e.support, p ≤ 1 ∧ e p ≤ B)
    (hf : ∀ n, 1 ≤ n → f n = ratEval Q e ((2 * n : ℕ) : ℚ)) : Rep B f := by
  refine ⟨[(Q, e, [])], ?_, fun n hn => ?_⟩
  · intro τ hτ
    simp only [List.mem_singleton] at hτ
    subst hτ
    refine ⟨fun a ha => by simp at ha, by simp [wt], fun p hp => ⟨(he p hp).1, ?_⟩⟩
    simp only [wt, List.map_nil, List.sum_nil, zero_add]
    exact (he p hp).2
  · rw [hf n hn]; simp [S03.ST, Z]

lemma Rep_const (c : ℚ) : Rep 0 (fun _ => c) :=
  Rep_rat (Polynomial.C c) 0 (fun p hp => by simp at hp) (fun n _ => by simp [ratEval])

lemma Rep_Z (v : List ℤ) (hv : IsIndex v) : Rep (wt v) (fun n => Z v (2 * n)) := by
  refine ⟨[(1, 0, v)], ?_, fun n _ => ?_⟩
  · intro τ hτ
    simp only [List.mem_singleton] at hτ
    subst hτ
    exact ⟨hv, le_refl _, fun p hp => by simp at hp⟩
  · simp [S03.ST, ratEval]

lemma ratEval_mul (Q1 Q2 : Polynomial ℚ) (e1 e2 : ℤ →₀ ℕ) (N : ℚ) :
    ratEval (Q1 * Q2) (e1 + e2) N = ratEval Q1 e1 N * ratEval Q2 e2 N := by
  unfold ratEval
  rw [Finsupp.prod_add_index' (fun _ => by simp) (fun _ _ _ => by rw [pow_add]),
    Polynomial.eval_mul, div_mul_div_comm]

lemma term_mul (B1 B2 : ℕ) (τ1 τ2 : DTerm) (h1 : InDB B1 [τ1]) (h2 : InDB B2 [τ2]) :
    ∃ T : List DTerm, InDB (B1 + B2) T ∧ ∀ N : ℕ, S03.ST T N =
      (ratEval τ1.1 τ1.2.1 (N : ℚ) * Z τ1.2.2 N) * (ratEval τ2.1 τ2.2.1 (N : ℚ) * Z τ2.2.2 N) := by
  obtain ⟨i1, w1b, p1⟩ := h1 τ1 (List.mem_singleton_self _)
  obtain ⟨i2, w2b, p2⟩ := h2 τ2 (List.mem_singleton_self _)
  obtain ⟨L, hL, hZ⟩ := stuffle_list τ1.2.2 τ2.2.2 i1 i2
  refine ⟨L.map (fun w => (τ1.1 * τ2.1, τ1.2.1 + τ2.2.1, w)), ?_, fun N => ?_⟩
  · intro τ hτ
    obtain ⟨w, hw, rfl⟩ := List.mem_map.mp hτ
    obtain ⟨hwi, hww⟩ := hL w hw
    refine ⟨hwi, by simp only; omega, fun p hp => ?_⟩
    simp only at hp ⊢
    have hp' := Finsupp.support_add hp
    rw [Finset.mem_union] at hp'
    rw [Finsupp.add_apply, hww]
    by_cases ha : p ∈ τ1.2.1.support
    · by_cases hb : p ∈ τ2.2.1.support
      · exact ⟨(p1 p ha).1, by have := (p1 p ha).2; have := (p2 p hb).2; omega⟩
      · rw [Finsupp.notMem_support_iff.mp hb]
        exact ⟨(p1 p ha).1, by have := (p1 p ha).2; omega⟩
    · rw [Finsupp.notMem_support_iff.mp ha]
      have hb : p ∈ τ2.2.1.support := by tauto
      exact ⟨(p2 p hb).1, by have := (p2 p hb).2; omega⟩
  · have : S03.ST (L.map (fun w => (τ1.1 * τ2.1, τ1.2.1 + τ2.2.1, w))) N =
        ratEval (τ1.1 * τ2.1) (τ1.2.1 + τ2.2.1) (N : ℚ) * (L.map fun w => Z w N).sum := by
      unfold S03.ST
      rw [List.map_map, ← List.sum_map_mul_left]
      rfl
    rw [this, ← hZ N, ratEval_mul]
    ring

lemma InDB_cons {B : ℕ} {τ : DTerm} {T : List DTerm} (h : InDB B (τ :: T)) :
    InDB B [τ] ∧ InDB B T :=
  ⟨fun σ hσ => h σ (by simp at hσ; simp [hσ]), fun σ hσ => h σ (List.mem_cons_of_mem _ hσ)⟩

lemma list_mul_single (B1 B2 : ℕ) (τ1 : DTerm) (h1 : InDB B1 [τ1]) :
    ∀ T2 : List DTerm, InDB B2 T2 → ∃ T : List DTerm, InDB (B1 + B2) T ∧ ∀ N : ℕ,
      S03.ST T N = (ratEval τ1.1 τ1.2.1 (N : ℚ) * Z τ1.2.2 N) * S03.ST T2 N := by
  intro T2
  induction T2 with
  | nil => intro _; exact ⟨[], InDB_nil _, fun N => by simp [S03.ST]⟩
  | cons τ2 T2 ih =>
    intro h2
    obtain ⟨ha, hb⟩ := InDB_cons h2
    obtain ⟨Ta, hTa, hSa⟩ := term_mul B1 B2 τ1 τ2 h1 ha
    obtain ⟨Tb, hTb, hSb⟩ := ih hb
    exact ⟨Ta ++ Tb, InDB_append hTa hTb, fun N => by
      rw [ST_append, hSa, hSb, ST_cons]; ring⟩

lemma list_mul (B1 B2 : ℕ) : ∀ T1 : List DTerm, InDB B1 T1 → ∀ T2 : List DTerm, InDB B2 T2 →
    ∃ T : List DTerm, InDB (B1 + B2) T ∧ ∀ N : ℕ, S03.ST T N = S03.ST T1 N * S03.ST T2 N := by
  intro T1
  induction T1 with
  | nil => intro _ T2 _; exact ⟨[], InDB_nil _, fun N => by simp [S03.ST]⟩
  | cons τ1 T1 ih =>
    intro h1 T2 h2
    obtain ⟨ha, hb⟩ := InDB_cons h1
    obtain ⟨Ta, hTa, hSa⟩ := list_mul_single B1 B2 τ1 ha T2 h2
    obtain ⟨Tb, hTb, hSb⟩ := ih hb T2 h2
    exact ⟨Ta ++ Tb, InDB_append hTa hTb, fun N => by
      rw [ST_append, hSa, hSb, ST_cons]; ring⟩

lemma Rep_mul {B1 B2 : ℕ} {f g : ℕ → ℚ} (hf : Rep B1 f) (hg : Rep B2 g) :
    Rep (B1 + B2) (fun n => f n * g n) := by
  obtain ⟨T1, h1, hf⟩ := hf
  obtain ⟨T2, h2, hg⟩ := hg
  obtain ⟨T, hT, hS⟩ := list_mul B1 B2 T1 h1 T2 h2
  exact ⟨T, hT, fun n hn => by show f n * g n = _; rw [hS, hf n hn, hg n hn]⟩

lemma Rep_cmul {B : ℕ} (c : ℚ) {f : ℕ → ℚ} (hf : Rep B f) : Rep B (fun n => c * f n) := by
  have := Rep_mul (Rep_const c) hf
  rw [zero_add] at this
  exact this

lemma Rep_neg {B : ℕ} {f : ℕ → ℚ} (hf : Rep B f) : Rep B (fun n => -f n) :=
  Rep_congr (Rep_cmul (-1) hf) (fun n _ => by ring)

/-! ## Harmonic numbers -/

lemma HrQ_shift (r M : ℕ) : ∀ j : ℕ,
    V05.HrQ r (M + j) = V05.HrQ r M + ∑ l ∈ Finset.range j, 1 / ((M : ℚ) + l + 1) ^ r := by
  intro j
  induction j with
  | zero => simp
  | succ j ih =>
    rw [← add_assoc, V05.HrQ, Finset.sum_range_succ, ← V05.HrQ, ih, Finset.sum_range_succ]
    push_cast; ring

lemma hn_eq_HrQ (r : ℕ) : ∀ N : ℕ, hn r N = V05.HrQ r N := by
  intro N
  induction N with
  | zero => simp [hn, V05.HrQ]
  | succ N ih => rw [hn_succ, ih, V05.HrQ, V05.HrQ, Finset.sum_range_succ]

lemma Rep_inv_pow (r l : ℕ) : Rep r (fun n => 1 / ((n : ℚ) + l) ^ r) := by
  refine Rep_rat (Polynomial.C ((2 : ℚ) ^ r)) (Finsupp.single (-(2 * (l : ℤ))) r) ?_ ?_
  · intro p hp
    have := Finsupp.support_single_subset hp
    rw [Finset.mem_singleton] at this
    subst this
    refine ⟨by omega, by simp⟩
  · intro n hn
    unfold ratEval
    rw [Finsupp.prod_single_index (by simp)]
    have h1 : ((n : ℚ) + l) ≠ 0 := by positivity
    simp only [Polynomial.eval_C]
    push_cast
    rw [show (2 * (n : ℚ) - -(2 * (l : ℚ))) = 2 * ((n : ℚ) + l) by ring, mul_pow]
    field_simp

lemma Rep_inv_pow2 (r l : ℕ) : Rep r (fun n => 1 / (2 * (n : ℚ) + l) ^ r) := by
  refine Rep_rat 1 (Finsupp.single (-(l : ℤ)) r) ?_ ?_
  · intro p hp
    have := Finsupp.support_single_subset hp
    rw [Finset.mem_singleton] at this
    subst this
    refine ⟨by omega, by simp⟩
  · intro n hn
    unfold ratEval
    rw [Finsupp.prod_single_index (by simp)]
    simp only [Polynomial.eval_one]
    push_cast
    rw [show (2 * (n : ℚ) - -(l : ℚ)) = 2 * (n : ℚ) + l by ring]

lemma wt_single (a : ℤ) : wt [a] = a.natAbs := by simp [wt]

lemma Rep_H2 (r j : ℕ) (hr : 1 ≤ r) : Rep r (fun n => V05.HrQ r (n + j - 1)) := by
  have hZ1 := Rep_cmul ((2 : ℚ) ^ (r - 1)) (Rep_Z [(r : ℤ)] (by intro a ha; simp at ha; omega))
  have hZ2 := Rep_cmul ((2 : ℚ) ^ (r - 1)) (Rep_Z [-(r : ℤ)] (by intro a ha; simp at ha; omega))
  rw [wt_single, Int.natAbs_natCast] at hZ1
  rw [wt_single, Int.natAbs_neg, Int.natAbs_natCast] at hZ2
  have hI := Rep_neg (Rep_inv_pow r 0)
  have hS := Rep_sum (Finset.range j) (B := r) (fun l n => 1 / ((n : ℚ) + l) ^ r)
    (fun l _ => Rep_inv_pow r l)
  refine Rep_congr (Rep_add (Rep_add (Rep_add hZ1 hZ2) hI) hS) (fun n hn => ?_)
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  rw [show m + 1 + j - 1 = m + j by omega, HrQ_shift r m j, ← hn_eq_HrQ]
  have h3 := (proof3 r (m + 1) hr (by omega)).1
  rw [show m + 1 - 1 = m by omega] at h3
  rw [← h3]
  simp only [Nat.cast_add, Nat.cast_one, Nat.cast_zero, add_zero]
  have : ∀ l ∈ Finset.range j, 1 / ((m : ℚ) + l + 1) ^ r = 1 / ((m : ℚ) + 1 + l) ^ r :=
    fun l _ => by ring_nf
  rw [Finset.sum_congr rfl this]
  ring

lemma Rep_H1 (r : ℕ) (hr : 1 ≤ r) : Rep r (fun n => V05.HrQ r (n - 1)) :=
  Rep_congr (Rep_H2 r 0 hr) (fun n _ => by simp)

lemma Rep_H3 (r j : ℕ) (hr : 1 ≤ r) : Rep r (fun n => V05.HrQ r (n + j + n - 1)) := by
  have hZ1 := Rep_Z [(r : ℤ)] (by intro a ha; simp at ha; omega)
  rw [wt_single, Int.natAbs_natCast] at hZ1
  have hI := Rep_neg (Rep_inv_pow2 r 0)
  have hS := Rep_sum (Finset.range j) (B := r) (fun l n => 1 / (2 * (n : ℚ) + l) ^ r)
    (fun l _ => Rep_inv_pow2 r l)
  refine Rep_congr (Rep_add (Rep_add hZ1 hI) hS) (fun n hn => ?_)
  rw [show n + j + n - 1 = (2 * n - 1) + j by omega, HrQ_shift r (2 * n - 1) j, ← hn_eq_HrQ]
  have h3 := (proof3 r n hr hn).2
  rw [← h3]
  have h2 : ((2 * n - 1 : ℕ) : ℚ) + 1 = 2 * (n : ℚ) := by
    rw [Nat.cast_sub (by omega)]; push_cast; ring
  have : ∀ l ∈ Finset.range j, 1 / (((2 * n - 1 : ℕ) : ℚ) + l + 1) ^ r = 1 / (2 * (n : ℚ) + l) ^ r :=
    fun l _ => by rw [add_right_comm, h2]
  rw [Finset.sum_congr rfl this]
  simp only [Nat.cast_zero, add_zero]
  ring

/-! ## The rational prefactor -/

lemma esum_apply (S : Finset ℤ) (p : ℤ) :
    (∑ q ∈ S, Finsupp.single q (1 : ℕ)) p = if p ∈ S then 1 else 0 := by
  rw [Finsupp.finsetSum_apply]
  rw [Finset.sum_congr rfl (fun q _ => show Finsupp.single q (1 : ℕ) p =
    if q = p then 1 else 0 by rw [Finsupp.single_apply])]
  rw [Finset.sum_ite_eq']

lemma esum_prod (S : Finset ℤ) (N : ℚ) :
    (∑ q ∈ S, Finsupp.single q (1 : ℕ)).prod (fun p k => (N - (p : ℚ)) ^ k) =
      ∏ q ∈ S, (N - (q : ℚ)) := by
  rw [← Finsupp.prod_finsetSum_index (h := fun (p : ℤ) (k : ℕ) => (N - (p : ℚ)) ^ k)
    (fun _ => by simp) (fun _ _ _ => by rw [pow_add])]
  refine Finset.prod_congr rfl (fun q _ => ?_)
  rw [Finsupp.prod_single_index (by simp), pow_one]

lemma Rep_rat1 {f : ℕ → ℚ} (S : Finset ℤ) (hS : ∀ p ∈ S, p ≤ 1) (P : Polynomial ℚ)
    (hf : ∀ n, 1 ≤ n → f n * ∏ q ∈ S, ((2 * n : ℕ) - (q : ℚ)) = P.eval (n : ℚ)) : Rep 1 f := by
  refine Rep_rat (P.comp (Polynomial.C (1 / 2 : ℚ) * Polynomial.X))
    (∑ q ∈ S, Finsupp.single q (1 : ℕ)) ?_ ?_
  · intro p hp
    rw [Finsupp.mem_support_iff, esum_apply] at hp
    rw [esum_apply]
    split_ifs at hp ⊢ with h
    · exact ⟨hS p h, le_refl _⟩
    · exact absurd rfl hp
  · intro n hn
    unfold ratEval
    rw [esum_prod, Polynomial.eval_comp]
    simp only [Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_X]
    rw [show (1 / 2 : ℚ) * ((2 * n : ℕ) : ℚ) = n by push_cast; ring, ← hf n hn]
    have hne : ∏ q ∈ S, ((2 * n : ℕ) - (q : ℚ)) ≠ 0 := by
      rw [Finset.prod_ne_zero_iff]
      intro q hq
      have h1 : (q : ℚ) ≤ 1 := by exact_mod_cast hS q hq
      have h2 : (2 : ℚ) ≤ ((2 * n : ℕ) : ℚ) := by
        push_cast; have : (1 : ℚ) ≤ n := by exact_mod_cast hn
        linarith
      linarith
    field_simp

/-- The pole set `{1, 0, -1, …, -j}`. -/
def Spol (j : ℕ) : Finset ℤ := insert 1 ((Finset.range (j + 1)).image (fun t : ℕ => -(t : ℤ)))

lemma Spol_le (j : ℕ) : ∀ p ∈ Spol j, p ≤ 1 := by
  intro p hp
  simp only [Spol, Finset.mem_insert, Finset.mem_image] at hp
  rcases hp with rfl | ⟨t, _, rfl⟩
  · exact le_refl _
  · omega

lemma Spol_prod (j : ℕ) (N : ℚ) :
    ∏ q ∈ Spol j, (N - (q : ℚ)) = (N - 1) * ∏ t ∈ Finset.range (j + 1), (N + t) := by
  unfold Spol
  rw [Finset.prod_insert (by simp; omega), Finset.prod_image (by intro a _ b _ h; simp at h; exact h)]
  push_cast
  congr 1
  refine Finset.prod_congr rfl (fun t _ => by ring)

lemma Rep_phi (k j i : ℕ) : Rep 1 (fun n =>
    -(1 / ((n : ℚ) * (2 * (n : ℚ) - 1))) * (Q04.eF k n j i (k - i) * V05.rhoQ j n)) := by
  obtain ⟨Pe, hPe⟩ := Q04.eF_poly k j i (k - i)
  rcases j with _ | j
  · refine Rep_rat1 (Spol 0) (Spol_le 0) (Polynomial.C (-2 : ℚ) * Pe) (fun n hn => ?_)
    rw [Spol_prod, hPe n]
    have h1 : (n : ℚ) ≠ 0 := by positivity
    have h2 : (2 * (n : ℚ) - 1) ≠ 0 := by
      have : (1 : ℚ) ≤ n := by exact_mod_cast hn
      linarith
    simp only [V05.rhoQ, Finset.range_zero, Finset.prod_empty, Finset.prod_range_succ,
      Polynomial.eval_mul, Polynomial.eval_C]
    push_cast
    have h2' : (n : ℚ) * 2 - 1 ≠ 0 := by rw [mul_comm]; exact h2
    rw [show ((2 : ℚ) * n - 1) * (1 * (2 * n + 0)) = 2 * (n * (2 * n - 1)) by ring]
    field_simp
  · refine Rep_rat1 (Spol j) (Spol_le j)
      (Polynomial.C (-1 : ℚ) * Pe * ∏ t ∈ Finset.range j, (Polynomial.X + Polynomial.C ((t : ℚ) + 1)))
      (fun n hn => ?_)
    rw [Spol_prod, hPe n]
    have h1 : (n : ℚ) ≠ 0 := by positivity
    have h2 : (2 * (n : ℚ) - 1) ≠ 0 := by
      have : (1 : ℚ) ≤ n := by exact_mod_cast hn
      linarith
    have h3 : ∏ t ∈ Finset.range (j + 1), (2 * (n : ℚ) + t) ≠ 0 := by
      rw [Finset.prod_ne_zero_iff]; intro t _; positivity
    unfold V05.rhoQ
    rw [Finset.prod_div_distrib, Finset.prod_range_succ' (fun t => (n : ℚ) + t)]
    simp only [Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_prod, Polynomial.eval_add,
      Polynomial.eval_X]
    push_cast
    have h2' : (n : ℚ) * 2 - 1 ≠ 0 := by rw [mul_comm]; exact h2
    field_simp
    ring

/-! ## Families in `ℚ[𝔷][x, x']` -/

/-- `f` is represented in `𝔇^{(d - w)}` when `w ≤ d`, and vanishes when `d < w`. -/
def Rp (d w : ℕ) (f : ℕ → ℚ) : Prop :=
  (w ≤ d → Rep (d - w) f) ∧ (d < w → ∀ n, 1 ≤ n → f n = 0)

lemma Rp_zero' (d w : ℕ) {f : ℕ → ℚ} (h : ∀ n, 1 ≤ n → f n = 0) : Rp d w f :=
  ⟨fun _ => Rep_zero' _ h, fun _ => h⟩

lemma Rp_congr {d w : ℕ} {f g : ℕ → ℚ} (h : Rp d w f) (hfg : ∀ n, 1 ≤ n → g n = f n) :
    Rp d w g :=
  ⟨fun hw => Rep_congr (h.1 hw) hfg, fun hw n hn => by rw [hfg n hn, h.2 hw n hn]⟩

lemma Rp_add {d w : ℕ} {f g : ℕ → ℚ} (hf : Rp d w f) (hg : Rp d w g) :
    Rp d w (fun n => f n + g n) :=
  ⟨fun hw => Rep_add (hf.1 hw) (hg.1 hw), fun hw n hn => by
    show f n + g n = 0; rw [hf.2 hw n hn, hg.2 hw n hn, add_zero]⟩

lemma Rp_sum {ι : Type*} (S : Finset ι) {d w : ℕ} (f : ι → ℕ → ℚ)
    (h : ∀ i ∈ S, Rp d w (f i)) : Rp d w (fun n => ∑ i ∈ S, f i n) :=
  ⟨fun hw => Rep_sum S f (fun i hi => (h i hi).1 hw), fun hw n hn =>
    Finset.sum_eq_zero (fun i hi => (h i hi).2 hw n hn)⟩

lemma Rp_cmul {d w : ℕ} (c : ℚ) {f : ℕ → ℚ} (hf : Rp d w f) : Rp d w (fun n => c * f n) :=
  ⟨fun hw => Rep_cmul c (hf.1 hw), fun hw n hn => by show c * f n = 0; rw [hf.2 hw n hn, mul_zero]⟩

lemma Rp_mul {d1 w1 d2 w2 : ℕ} {f g : ℕ → ℚ} (hf : Rp d1 w1 f) (hg : Rp d2 w2 g) :
    Rp (d1 + d2) (w1 + w2) (fun n => f n * g n) := by
  by_cases h1 : w1 ≤ d1
  · by_cases h2 : w2 ≤ d2
    · refine ⟨fun _ => ?_, fun h => by omega⟩
      have := Rep_mul (hf.1 h1) (hg.1 h2)
      rwa [show d1 - w1 + (d2 - w2) = d1 + d2 - (w1 + w2) by omega] at this
    · exact Rp_zero' _ _ (fun n hn => by
        show f n * g n = 0; rw [hg.2 (by omega) n hn, mul_zero])
  · exact Rp_zero' _ _ (fun n hn => by
      show f n * g n = 0; rw [hf.2 (by omega) n hn, zero_mul])

lemma Rp_const (d w : ℕ) (c : ℚ) (h : c ≠ 0 → w ≤ d) : Rp d w (fun _ => c) := by
  by_cases hc : c = 0
  · exact Rp_zero' _ _ (fun _ _ => hc)
  · exact ⟨fun _ => Rep_mono (Rep_const c) (Nat.zero_le _), fun hw => absurd (h hc) (by omega)⟩

lemma mwt_add (s1 s2 : ℕ →₀ ℕ) : mwt (s1 + s2) = mwt s1 + mwt s2 := by
  unfold mwt
  rw [Finsupp.sum_add_index' (fun _ => by simp) (fun _ _ _ => by ring)]

/-- Every coefficient `[x^{t_0} x'^{t_1} 𝔷^s]` is represented in `𝔇^{(t_0 + t_1 - wt s)}`. -/
def GF (F : ℕ → Q04.S2) : Prop :=
  ∀ (t : Fin 2 →₀ ℕ) (s : ℕ →₀ ℕ), Rp (t 0 + t 1) (mwt s) (fun n => ((F n).coeff t).coeff s)

lemma GF_add {F G : ℕ → Q04.S2} (hF : GF F) (hG : GF G) : GF (fun n => F n + G n) :=
  fun t s => Rp_congr (Rp_add (hF t s) (hG t s)) (fun n _ => by
    rw [MvPolynomial.coeff_add, MvPolynomial.coeff_add])

lemma GF_cmul (c : ℚ) {F : ℕ → Q04.S2} (hF : GF F) : GF (fun n => c • F n) :=
  fun t s => Rp_congr (Rp_cmul c (hF t s)) (fun n _ => by
    simp only [MvPolynomial.coeff_smul, smul_eq_mul])

lemma GF_sub {F G : ℕ → Q04.S2} (hF : GF F) (hG : GF G) : GF (fun n => F n - G n) := by
  have := GF_add hF (GF_cmul (-1) hG)
  refine fun t s => Rp_congr (this t s) (fun n _ => ?_)
  show ((F n - G n).coeff t).coeff s = ((F n + (-1 : ℚ) • G n).coeff t).coeff s
  rw [sub_eq_add_neg, ← neg_one_smul ℚ (G n)]

lemma GF_sum {ι : Type*} (S : Finset ι) (F : ι → ℕ → Q04.S2) (h : ∀ i ∈ S, GF (F i)) :
    GF (fun n => ∑ i ∈ S, F i n) :=
  fun t s => Rp_congr (Rp_sum S (fun i n => ((F i n).coeff t).coeff s) (fun i hi => h i hi t s))
    (fun n _ => by simp only [MvPolynomial.coeff_sum])

lemma GF_const (G : Q04.S2) (hG : Q04.GoodD G) : GF (fun _ => G) := by
  intro t s
  refine Rp_const _ _ _ (fun hc => ?_)
  exact ((hG t) s (MvPolynomial.mem_support_iff.mpr hc)).2

lemma GF_mul {F G : ℕ → Q04.S2} (hF : GF F) (hG : GF G) : GF (fun n => F n * G n) := by
  classical
  intro t s
  have key : ∀ n, ((F n * G n).coeff t).coeff s = ∑ x ∈ Finset.HasAntidiagonal.antidiagonal t,
      ∑ y ∈ Finset.HasAntidiagonal.antidiagonal s, ((F n).coeff x.1).coeff y.1 * ((G n).coeff x.2).coeff y.2 := by
    intro n
    rw [MvPolynomial.coeff_mul, MvPolynomial.coeff_sum]
    refine Finset.sum_congr rfl (fun x _ => ?_)
    rw [MvPolynomial.coeff_mul]
  refine Rp_congr (Rp_sum _ _ (fun x hx => Rp_sum _ _ (fun y hy => ?_))) (fun n _ => key n)
  rw [Finset.HasAntidiagonal.mem_antidiagonal] at hx hy
  have := Rp_mul (hF x.1 y.1) (hG x.2 y.2)
  rw [← mwt_add, hy] at this
  have ht : x.1 0 + x.1 1 + (x.2 0 + x.2 1) = t 0 + t 1 := by
    rw [← hx]; simp only [Finsupp.add_apply]; ring
  rw [ht] at this
  exact this

lemma GF_one : GF (fun _ => (1 : Q04.S2)) := GF_const 1 Q04.GoodD_one

lemma GF_pow {F : ℕ → Q04.S2} (hF : GF F) (r : ℕ) : GF (fun n => F n ^ r) := by
  induction r with
  | zero => simpa using GF_one
  | succ r ih =>
    have := GF_mul ih hF
    simpa [pow_succ] using this

lemma GF_expT {F : ℕ → Q04.S2} (hF : GF F) (K : ℕ) : GF (fun n => Q04.expT K (F n)) := by
  unfold Q04.expT
  exact GF_sum _ _ (fun r _ => GF_cmul _ (GF_pow hF r))

lemma mwt_eq_zero_of_WL0 (q : MvPolynomial ℕ ℚ) (hq : Q04.WeightLE 0 q) (s : ℕ →₀ ℕ)
    (hs : q.coeff s ≠ 0) : mwt s = 0 := by
  have := (hq s (MvPolynomial.mem_support_iff.mpr hs)).2
  unfold mwt; omega

lemma GF_cG {r : ℕ} {c : ℕ → ℚ} (hc : Rep r c) (G : Q04.S2) (hG : Q04.Hm r G) :
    GF (fun n => c n • G) := by
  intro t s
  have hco : ∀ n, ((c n • G).coeff t).coeff s = (G.coeff t).coeff s * c n := by
    intro n; simp only [MvPolynomial.coeff_smul, smul_eq_mul]; ring
  by_cases h0 : (G.coeff t).coeff s = 0
  · exact Rp_zero' _ _ (fun n _ => by rw [hco, h0, zero_mul])
  · have hw := mwt_eq_zero_of_WL0 _ (hG t).1 s h0
    have ht := (hG t).2 (fun h => h0 (by rw [h]; simp))
    rw [hw, ht]
    refine ⟨fun _ => ?_, fun h => by omega⟩
    exact Rep_congr (Rep_cmul ((G.coeff t).coeff s) hc) (fun n _ => hco n)

lemma GF_Delta (K j : ℕ) : GF (fun n => Q04.DeltaT K (n + j) n) := by
  unfold Q04.DeltaT
  dsimp only
  refine GF_add (GF_sum _ _ fun r hr => GF_cmul _ ?_) (GF_const _ ?_)
  · have hr1 : 1 ≤ r := (Finset.mem_Icc.mp hr).1
    have hx := Q04.Hm_X 0
    have hx' := Q04.Hm_X 1
    refine GF_sub (GF_add ?_ ?_) ?_
    · have := GF_cG (Rep_H2 r j hr1) _ (by simpa using Q04.Hm_pow hx r)
      simpa using this
    · have := GF_cG (Rep_H1 r hr1) _ (by simpa using Q04.Hm_pow hx' r)
      simpa using this
    · have := GF_cG (Rep_H3 r j hr1) _ (by simpa using Q04.Hm_pow (Q04.Hm_add hx hx') r)
      simpa using this
  · refine Q04.GoodD_sum _ _ fun w hw => Q04.GoodD_smul _ ?_
    have h := Q04.Hm_sub (Q04.Hm_sub (Q04.Hm_pow (Q04.Hm_add (Q04.Hm_X 0) (Q04.Hm_X 1)) w)
      (Q04.Hm_pow (Q04.Hm_X 0) w)) (Q04.Hm_pow (Q04.Hm_X 1) w)
    simp only [mul_one] at h
    exact Q04.GoodD_CX (Finset.mem_Icc.mp hw).1 h

end S08.PA

namespace S08.PA

open S03

lemma Rp_Dsharp (i m j : ℕ) (s : ℕ →₀ ℕ) :
    Rp (i + m) (mwt s) (fun n => (Q04.DsharpF i m (n + j) n).coeff s) := by
  have h := GF_expT (GF_Delta (i + m) j) (i + m) (Finsupp.single 0 i + Finsupp.single 1 m) s
  simp only [Finsupp.add_apply, Finsupp.single_apply] at h
  norm_num at h
  refine Rp_congr (Rp_cmul (((i.factorial * m.factorial : ℕ) : ℚ)) h) (fun n _ => ?_)
  unfold Q04.DsharpF
  rw [MvPolynomial.coeff_smul, smul_eq_mul]

lemma Pf_coeff (k n : ℕ) (s : ℕ →₀ ℕ) : (Q04.Pf k n).coeff s =
    ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
      (-(1 / ((n : ℚ) * (2 * (n : ℚ) - 1))) * (Q04.eF k n j i (k - i) * V05.rhoQ j n)) *
        (Q04.DsharpF i (k - i) (n + j) n).coeff s := by
  unfold Q04.Pf
  rw [MvPolynomial.coeff_smul, MvPolynomial.coeff_sum, smul_eq_mul, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  rw [MvPolynomial.coeff_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [MvPolynomial.coeff_smul, smul_eq_mul]
  ring

theorem main : S08.stmt_4 := by
  intro k s _ hW
  have h : Rep (k + 1 - mwt s) (fun n => (Q04.Pf k n).coeff s) := by
    refine Rep_congr (Rep_sum _ _ (fun j _ => Rep_sum _ _ (fun i hi => ?_)))
      (fun n _ => Pf_coeff k n s)
    have hik : i + (k - i) = k := Nat.add_sub_cancel' (Nat.lt_succ_iff.mp (Finset.mem_range.mp hi))
    have hD := Rp_Dsharp i (k - i) j s
    rw [hik] at hD
    have := Rep_mul (Rep_phi k j i) (hD.1 hW)
    exact Rep_mono this (by omega)
  obtain ⟨T, hT, hf⟩ := h
  exact ⟨T, hT, hf⟩

end S08.PA

-- ===== merged from work/S08/PartB.lean =====
/-!
# S08, Part B: Lemma `g6:prop:R` (`S08.stmt_8`)

`Good B c`: some `F ∈ 𝒢_B` has a Laurent expansion `F(x) = q(1/x) + Σ a_N x^N` on `0 < |x| < 1`
with `a_N = c N` for even `N ≥ 2`. `stmt_8` is `∀ B T, InDB B T → Good B (ST T)`.
-/

namespace S08.PB

open Real MeasureTheory Set Filter Polynomial S03

/-- `PS B a`: some `F ∈ 𝒢_B` equals `q(1/x) + Σ a_N x^N` on `0 < |x| < 1`. -/
def PS (B : ℕ) (a : ℕ → ℚ) : Prop :=
  ∃ F : ℝ → ℝ, InG B F ∧ ∃ q : Polynomial ℚ,
    ∀ x : ℝ, 0 < |x| → |x| < 1 →
      HasSum (fun N : ℕ => (a N : ℝ) * x ^ N) (F x - Polynomial.aeval (1 / x) q)

/-- `Good B c`: the coefficients of some `F ∈ 𝒢_B` agree with `c` at even `N ≥ 2`. -/
def Good (B : ℕ) (c : ℕ → ℚ) : Prop :=
  ∃ a : ℕ → ℚ, PS B a ∧ ∀ N : ℕ, 2 ≤ N → Even N → a N = c N

lemma kval_C_mul (r : ℚ) (p : Polynomial ℚ) (a b d : ℕ) (x : ℝ) :
    kval (C r * p) a b d x = (r : ℝ) * kval p a b d x := by
  unfold kval
  simp only [map_mul, aeval_C, eq_ratCast]
  ring

lemma InG_zero (B : ℕ) : InG B (fun _ => 0) :=
  ⟨[], by simp, fun x _ _ => by simp⟩

lemma InG_add (B : ℕ) (F G : ℝ → ℝ) (hF : InG B F) (hG : InG B G) :
    InG B (fun x => F x + G x) := by
  obtain ⟨L1, h1, e1⟩ := hF
  obtain ⟨L2, h2, e2⟩ := hG
  refine ⟨L1 ++ L2, fun g hg => ?_, fun x hx0 hx1 => ?_⟩
  · rcases List.mem_append.mp hg with hg | hg
    · exact h1 g hg
    · exact h2 g hg
  · show F x + G x = _
    rw [e1 x hx0 hx1, e2 x hx0 hx1, List.map_append, List.sum_append]

lemma InG_smul (B : ℕ) (r : ℚ) (F : ℝ → ℝ) (hF : InG B F) :
    InG B (fun x => (r : ℝ) * F x) := by
  obtain ⟨L, h, e⟩ := hF
  refine ⟨L.map (fun g => (C r * g.1, g.2)), fun g hg => ?_, fun x hx0 hx1 => ?_⟩
  · rw [List.mem_map] at hg
    obtain ⟨g0, hg0, rfl⟩ := hg
    exact h g0 hg0
  · show (r : ℝ) * F x = _
    rw [e x hx0 hx1, List.map_map, ← List.sum_map_mul_left]
    congr 1
    refine List.map_congr_left (fun g _ => ?_)
    simp only [Function.comp, kval_C_mul]
    ring

lemma InG_mono (B B' : ℕ) (h : B ≤ B') (F : ℝ → ℝ) (hF : InG B F) : InG B' F := by
  obtain ⟨L, hL, e⟩ := hF
  exact ⟨L, fun g hg => ⟨(hL g hg).1, (hL g hg).2.1, (hL g hg).2.2.trans h⟩, e⟩

lemma PS_add (B : ℕ) (a b : ℕ → ℚ) (ha : PS B a) (hb : PS B b) :
    PS B (fun N => a N + b N) := by
  obtain ⟨F, hF, q, hq⟩ := ha
  obtain ⟨G, hG, r, hr⟩ := hb
  refine ⟨_, InG_add B F G hF hG, q + r, fun x hx0 hx1 => ?_⟩
  have := (hq x hx0 hx1).add (hr x hx0 hx1)
  convert this using 1
  · ext N; push_cast; ring
  · rw [map_add]; ring

lemma PS_smul (B : ℕ) (r : ℚ) (a : ℕ → ℚ) (ha : PS B a) :
    PS B (fun N => r * a N) := by
  obtain ⟨F, hF, q, hq⟩ := ha
  refine ⟨_, InG_smul B r F hF, C r * q, fun x hx0 hx1 => ?_⟩
  have := (hq x hx0 hx1).mul_left (r : ℝ)
  convert this using 1
  · ext N; push_cast; ring
  · rw [map_mul, aeval_C, eq_ratCast]; ring

lemma PS_mono (B B' : ℕ) (h : B ≤ B') (a : ℕ → ℚ) (ha : PS B a) : PS B' a := by
  obtain ⟨F, hF, q, hq⟩ := ha
  exact ⟨F, InG_mono B B' h F hF, q, hq⟩

lemma Good_zero (B : ℕ) : Good B (fun _ => 0) :=
  ⟨fun _ => 0, ⟨fun _ => 0, InG_zero B, 0, fun x _ _ => by simp⟩,
    fun _ _ _ => rfl⟩

lemma Good_add (B : ℕ) (c d : ℕ → ℚ) (hc : Good B c) (hd : Good B d) :
    Good B (fun N => c N + d N) := by
  obtain ⟨a, ha, ea⟩ := hc
  obtain ⟨b, hb, eb⟩ := hd
  exact ⟨_, PS_add B a b ha hb, fun N h2 he => by simp only [ea N h2 he, eb N h2 he]⟩

lemma Good_smul (B : ℕ) (r : ℚ) (c : ℕ → ℚ) (hc : Good B c) :
    Good B (fun N => r * c N) := by
  obtain ⟨a, ha, ea⟩ := hc
  exact ⟨_, PS_smul B r a ha, fun N h2 he => by simp only [ea N h2 he]⟩

lemma Good_congr (B : ℕ) (c d : ℕ → ℚ) (hc : Good B c)
    (h : ∀ N : ℕ, 2 ≤ N → Even N → c N = d N) : Good B d := by
  obtain ⟨a, ha, ea⟩ := hc
  exact ⟨a, ha, fun N h2 he => (ea N h2 he).trans (h N h2 he)⟩

lemma Good_sub (B : ℕ) (c d : ℕ → ℚ) (hc : Good B c) (hd : Good B d) :
    Good B (fun N => c N - d N) :=
  Good_congr B _ _ (Good_add B c _ hc (Good_smul B (-1) d hd)) (fun N _ _ => by ring)

lemma Good_mono (B B' : ℕ) (h : B ≤ B') (c : ℕ → ℚ) (hc : Good B c) : Good B' c := by
  obtain ⟨a, ha, ea⟩ := hc
  exact ⟨a, PS_mono B B' h a ha, ea⟩

lemma Good_finsum {ι : Type*} (B : ℕ) (s : Finset ι) (c : ι → ℕ → ℚ)
    (h : ∀ i ∈ s, Good B (c i)) : Good B (fun N => ∑ i ∈ s, c i N) := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using Good_zero B
  | insert i s hi ih =>
    refine Good_congr B _ _ (Good_add B (c i) _ (h i (by simp))
      (ih (fun j hj => h j (by simp [hj])))) (fun N _ _ => ?_)
    rw [Finset.sum_insert hi]

lemma Good_list {α : Type*} (B : ℕ) (T : List α) (c : α → ℕ → ℚ)
    (h : ∀ τ ∈ T, Good B (c τ)) : Good B (fun N => (T.map fun τ => c τ N).sum) := by
  induction T with
  | nil => simpa using Good_zero B
  | cons τ T ih =>
    refine Good_congr B _ _ (Good_add B (c τ) _ (h τ (by simp))
      (ih (fun j hj => h j (by simp [hj])))) (fun N _ _ => ?_)
    simp

lemma kval_shift (p : Polynomial ℚ) (a b d m : ℕ) (x : ℝ) :
    kval p (a + m) b d x = kval p a b d x / x ^ m := by
  unfold kval
  rw [pow_add, div_div]
  congr 1
  ring

lemma list_sum_div {α : Type*} (L : List α) (f : α → ℝ) (c : ℝ) :
    (L.map f).sum / c = (L.map fun g => f g / c).sum := by
  induction L with
  | nil => simp
  | cons g L ih => simp only [List.map_cons, List.sum_cons, add_div, ih]

lemma InG_shift (B m : ℕ) (F : ℝ → ℝ) (hF : InG B F) : InG B (fun x => F x / x ^ m) := by
  obtain ⟨L, hL, e⟩ := hF
  refine ⟨L.map (fun g => (g.1, g.2.1 + m, g.2.2)), fun g hg => ?_, fun x hx0 hx1 => ?_⟩
  · rw [List.mem_map] at hg
    obtain ⟨g0, hg0, rfl⟩ := hg
    exact hL g0 hg0
  · show F x / x ^ m = _
    rw [e x hx0 hx1, list_sum_div, List.map_map]
    congr 1
    refine List.map_congr_left (fun g _ => ?_)
    simp only [Function.comp, kval_shift]
    ring

lemma pow_div_eq (x : ℝ) (hx : x ≠ 0) (i m : ℕ) (h : i < m) :
    x ^ i / x ^ m = (1 / x) ^ (m - i) := by
  rw [div_pow, one_pow, show x ^ m = x ^ i * x ^ (m - i) by rw [← pow_add]; congr 1; omega]
  have : x ^ i ≠ 0 := pow_ne_zero _ hx
  have : x ^ (m - i) ≠ 0 := pow_ne_zero _ hx
  field_simp

lemma PS_shift (B m : ℕ) (a : ℕ → ℚ) (ha : PS B a) : PS B (fun N => a (N + m)) := by
  obtain ⟨F, hF, q, hq⟩ := ha
  refine ⟨_, InG_shift B m F hF, X ^ m * q + ∑ i ∈ Finset.range m, C (a i) * X ^ (m - i),
    fun x hx0 hx1 => ?_⟩
  have hx : x ≠ 0 := abs_pos.mp hx0
  have h2 := ((hasSum_nat_add_iff' m).mpr (hq x hx0 hx1)).div_const (x ^ m)
  convert h2 using 1
  · ext N
    rw [pow_add]
    have : x ^ m ≠ 0 := pow_ne_zero _ hx
    field_simp
  · simp only [map_add, map_mul, map_pow, aeval_X, map_sum, aeval_C, eq_ratCast]
    have hs : (∑ i ∈ Finset.range m, (a i : ℝ) * x ^ i) / x ^ m =
        ∑ i ∈ Finset.range m, (a i : ℝ) * (1 / x) ^ (m - i) := by
      rw [Finset.sum_div]
      exact Finset.sum_congr rfl (fun i hi => by
        rw [mul_div_assoc, pow_div_eq x hx i m (Finset.mem_range.mp hi)])
    rw [sub_div, sub_div, hs]
    have : x ^ m ≠ 0 := pow_ne_zero _ hx
    rw [div_pow, one_pow]
    field_simp
    ring

lemma PS_of_series (B k : ℕ) (b : ℕ → ℚ) (F : ℝ → ℝ) (hF : InG B F)
    (h : ∀ x : ℝ, 0 < |x| → |x| < 1 → HasSum (fun n => (b n : ℝ) * x ^ (n + k)) (F x)) :
    PS B (fun N => if N < k then 0 else b (N - k)) := by
  refine ⟨F, hF, 0, fun x hx0 hx1 => ?_⟩
  rw [map_zero, sub_zero]
  refine (hasSum_nat_add_iff' k).mp ?_
  have h0 : ∑ i ∈ Finset.range k, (((if i < k then 0 else b (i - k) : ℚ)) : ℝ) * x ^ i = 0 := by
    refine Finset.sum_eq_zero (fun i hi => ?_)
    simp [Finset.mem_range.mp hi]
  rw [h0, sub_zero]
  refine (h x hx0 hx1).congr_fun (fun n => ?_)
  simp [show ¬ (n + k < k) by omega]

lemma open_dom (x : ℝ) (hx0 : 0 < |x|) (hx1 : |x| < 1) :
    ∀ᶠ y in nhds x, 0 < |y| ∧ |y| < 1 := by
  have hopen : IsOpen (Ioo (-1 : ℝ) 1 \ {0}) := isOpen_Ioo.sdiff isClosed_singleton
  have hxs : x ∈ Ioo (-1 : ℝ) 1 \ {0} := ⟨abs_lt.mp hx1, by simpa [abs_pos] using hx0⟩
  filter_upwards [hopen.mem_nhds hxs] with y hy
  obtain ⟨hy1, hy0⟩ := hy
  have hy0' : y ≠ 0 := hy0
  exact ⟨abs_pos.mpr hy0', abs_lt.mpr hy1⟩

lemma series_deriv (c : ℕ → ℝ) (D : ℕ) (hc : ∀ n, |c n| ≤ ((n : ℝ) + 1) ^ D) (F : ℝ → ℝ)
    (F' : ℝ) (x : ℝ) (hx0 : 0 < |x|) (hx1 : |x| < 1)
    (hF : ∀ y : ℝ, 0 < |y| → |y| < 1 → HasSum (fun n => c n * y ^ (n + 1)) (F y))
    (hd : HasDerivAt F F' x) :
    HasSum (fun n => c n * ((n : ℝ) + 1) * x ^ (n + 1)) (x * F') := by
  obtain ⟨_, hsum, hder⟩ := psD c D hc x hx1
  have hev : F =ᶠ[nhds x] (fun y => ∑' n, c n * y ^ (n + 1)) := by
    filter_upwards [open_dom x hx0 hx1] with y hy
    exact (hF y hy.1 hy.2).tsum_eq.symm
  have h2 : HasDerivAt F (∑' n, c n * ((n : ℝ) + 1) * x ^ n) x := hder.congr_of_eventuallyEq hev
  rw [hd.unique h2, ← tsum_mul_left]
  have := hsum.hasSum.mul_left x
  rw [this.tsum_eq]
  refine this.congr_fun (fun n => ?_)
  ring

lemma theta_step (B D : ℕ) (c : ℕ → ℝ) (hc : ∀ n, |c n| ≤ ((n : ℝ) + 1) ^ D) (F : ℝ → ℝ)
    (hF : InG B F)
    (hs : ∀ y : ℝ, 0 < |y| → |y| < 1 → HasSum (fun n => c n * y ^ (n + 1)) (F y)) :
    ∃ G : ℝ → ℝ, InG B G ∧ ∀ x : ℝ, 0 < |x| → |x| < 1 →
      HasSum (fun n => c n * ((n : ℝ) + 1) * x ^ (n + 1)) (G x) := by
  obtain ⟨G, hG, hd⟩ := InG_theta B F hF
  refine ⟨G, hG, fun x hx0 hx1 => ?_⟩
  have := series_deriv c D hc F _ x hx0 hx1 hs (hd x hx0 hx1)
  have hx : x ≠ 0 := abs_pos.mp hx0
  rwa [mul_div_cancel₀ _ hx] at this

lemma theta_iter (B D : ℕ) (c : ℕ → ℝ) (hc : ∀ n, |c n| ≤ ((n : ℝ) + 1) ^ D) (F : ℝ → ℝ)
    (hF : InG B F)
    (hs : ∀ y : ℝ, 0 < |y| → |y| < 1 → HasSum (fun n => c n * y ^ (n + 1)) (F y)) :
    ∀ j : ℕ, ∃ G : ℝ → ℝ, InG B G ∧ ∀ x : ℝ, 0 < |x| → |x| < 1 →
      HasSum (fun n => c n * ((n : ℝ) + 1) ^ j * x ^ (n + 1)) (G x) := by
  intro j
  induction j with
  | zero => exact ⟨F, hF, fun x hx0 hx1 => by simpa using hs x hx0 hx1⟩
  | succ j ih =>
    obtain ⟨G, hG, hGs⟩ := ih
    have hb : ∀ n, |c n * ((n : ℝ) + 1) ^ j| ≤ ((n : ℝ) + 1) ^ (D + j) := by
      intro n
      rw [abs_mul, abs_pow, abs_of_pos (by positivity : (0 : ℝ) < n + 1), pow_add]
      exact mul_le_mul_of_nonneg_right (hc n) (by positivity)
    obtain ⟨G2, hG2, hG2s⟩ := theta_step B (D + j) _ hb G hG hGs
    refine ⟨G2, hG2, fun x hx0 hx1 => (hG2s x hx0 hx1).congr_fun (fun n => ?_)⟩
    ring

lemma C_sq (C : ℤ) (h : C = 1 ∨ C = -1) : (C : ℚ) * C = 1 := by
  rcases h with rfl | rfl <;> norm_num

/-- The blocks `β_{p,κ,v} = x^p H_{w[κ;v]}`, `p = 1 - m`: coefficient `(N-p)^{-κ} Z_v(N-p-1)`. -/
lemma block_beta (v : List ℤ) (hv : IsIndex v) (κ m : ℕ) (hκ : 1 ≤ κ) :
    Good (κ + wt v) (fun N => ratEval 1 (Finsupp.single (1 - (m : ℤ)) κ) (N : ℚ) *
      Z v (N - 2 + m)) := by
  obtain ⟨h1, h2, h3, h4, h5, h6⟩ := build_props v hv 1 κ (Or.inl rfl) hκ
  obtain ⟨⟨_, _, hC, _, hlen⟩, hS⟩ := dict (build 1 κ v) h1 h2 h3
  rw [h5, h6] at hlen
  have hPS : PS (κ + wt v) (fun N => if N < 2 then 0 else
      ((dd (build 1 κ v)).2.2.1 : ℚ) * (1 / (((N - 2 : ℕ) : ℚ) + 1) ^ κ) * Z v (N - 2)) := by
    refine PS_of_series _ 2 (fun n => ((dd (build 1 κ v)).2.2.1 : ℚ) * (1 / ((n : ℚ) + 1) ^ κ) *
      Z v n) (fun x => x * H (build 1 κ v) x) ⟨[(X, 0, 0, 0, build 1 κ v)], ?_, ?_⟩ ?_
    · intro g hg
      simp only [List.mem_singleton] at hg
      subst hg
      exact ⟨h2, h3, by simp [hlen]⟩
    · intro x _ _
      simp [kval]
    · intro x _ hx1
      refine ((hS x hx1).mul_left x).congr_fun (fun n => ?_)
      unfold cf
      rw [h4, h5, h6]
      push_cast
      ring
  have := PS_smul _ ((dd (build 1 κ v)).2.2.1 : ℚ) _ (PS_shift _ m _ hPS)
  refine ⟨_, this, fun N hN _ => ?_⟩
  simp only [show ¬ (N + m < 2) by omega, ite_false]
  rw [← mul_assoc, ← mul_assoc, C_sq _ hC, one_mul]
  unfold ratEval
  rw [Finsupp.prod_single_index (by simp), eval_one, show N - 2 + m = N + m - 2 by omega]
  congr 3
  rw [Nat.cast_sub (by omega)]
  push_cast
  ring

lemma build_one (v : List ℤ) : ∃ w', build 1 1 v = 1 :: w' := by
  cases v <;> simp [build]

/-- `θ H_{(1,w')} = x/(1-x) H_{w'}`, with its series. -/
lemma gamma_first (v : List ℤ) (hv : IsIndex v) :
    ∃ G : ℝ → ℝ, InG (wt v) G ∧ ∀ x : ℝ, 0 < |x| → |x| < 1 →
      HasSum (fun n => cf (build 1 1 v) n * ((n : ℝ) + 1) * x ^ (n + 1)) (G x) := by
  obtain ⟨h1, h2, h3, _, h5, h6⟩ := build_props v hv 1 1 (Or.inl rfl) le_rfl
  obtain ⟨hF, hS⟩ := dict (build 1 1 v) h1 h2 h3
  have hlen := hF.2.2.2.2
  rw [h5, h6] at hlen
  obtain ⟨w', hw'⟩ := build_one v
  rw [hw'] at h2 h3 hlen hS hF
  have hl' : IsLetters w' := fun b hb => h2 b (by simp [hb])
  have hr' : Regular w' := reg_tail 1 w' h3
  refine ⟨fun x => x * (f 1 x * H w' x), ⟨[(X, 0, 1, 0, w')], ?_, ?_⟩, ?_⟩
  · intro g hg
    simp only [List.mem_singleton] at hg
    subst hg
    refine ⟨hl', hr', ?_⟩
    simp only [List.length_cons] at hlen
    simp only; omega
  · intro x _ _
    simp [kval, f]
    ring
  · intro x hx0 hx1
    rw [hw']
    exact series_deriv _ _ (cf_bound _ hF) _ _ x hx0 hx1 (fun y _ hy1 => hS y hy1)
      (H_deriv 1 w' h2 h3 x hx0 hx1)

/-- The blocks `γ_{j,v} = θ^{j+1} H_{w[1;v]}`: coefficient `N^j Z_v(N-1)`. -/
lemma block_gamma (v : List ℤ) (hv : IsIndex v) (j : ℕ) :
    Good (wt v) (fun N => ratEval (X ^ j) 0 (N : ℚ) * Z v (N - 2 + 1)) := by
  obtain ⟨h1, h2, h3, h4, h5, h6⟩ := build_props v hv 1 1 (Or.inl rfl) le_rfl
  obtain ⟨hF, _⟩ := dict (build 1 1 v) h1 h2 h3
  have hC := hF.2.2.1
  obtain ⟨G0, hG0, hG0s⟩ := gamma_first v hv
  have hb : ∀ n, |cf (build 1 1 v) n * ((n : ℝ) + 1)| ≤
      ((n : ℝ) + 1) ^ ((dd (build 1 1 v)).2.2.2.length + 1) := by
    intro n
    rw [abs_mul, abs_of_pos (by positivity : (0 : ℝ) < n + 1), pow_succ]
    exact mul_le_mul_of_nonneg_right (cf_bound _ hF n) (by positivity)
  obtain ⟨G, hG, hGs⟩ := theta_iter _ _ _ hb G0 hG0 hG0s j
  have hPS := PS_of_series (wt v) 1 (fun n => ((dd (build 1 1 v)).2.2.1 : ℚ) *
      ((n : ℚ) + 1) ^ j * Z v n) G hG (fun x hx0 hx1 => (hGs x hx0 hx1).congr_fun (fun n => by
        unfold cf
        rw [h4, h5, h6]
        push_cast
        have : ((n : ℝ) + 1) ≠ 0 := by positivity
        field_simp
        simp))
  refine ⟨_, PS_smul _ ((dd (build 1 1 v)).2.2.1 : ℚ) _ hPS, fun N hN _ => ?_⟩
  simp only [show ¬ (N < 1) by omega, ite_false]
  rw [← mul_assoc, ← mul_assoc, C_sq _ hC, one_mul]
  unfold ratEval
  rw [Finsupp.prod_zero_index, eval_pow, eval_X, div_one, show N - 2 + 1 = N - 1 by omega]
  congr 2
  rw [Nat.cast_sub (by omega)]
  push_cast
  ring

lemma ratEval_mul (Q Q' : Polynomial ℚ) (e e' : ℤ →₀ ℕ) (N : ℚ) :
    ratEval Q e N * ratEval Q' e' N = ratEval (Q * Q') (e + e') N := by
  unfold ratEval
  rw [Finsupp.prod_add_index' (fun _ => pow_zero _) (fun _ _ _ => pow_add _ _ _), eval_mul,
    div_mul_div_comm]

lemma wt_take_drop (v : List ℤ) (t : ℕ) : wt (v.take t) + wt (v.drop t) = wt v := by
  conv_rhs => rw [← List.take_append_drop t v]
  simp [wt]

/-- The bound `T ∈ 𝔇^{(B)}` for one term `R_v [v]`. -/
def Bd (B : ℕ) (v : List ℤ) (e : ℤ →₀ ℕ) : Prop :=
  wt v ≤ B ∧ ∀ p ∈ e.support, p ≤ 1 ∧ wt v + e p ≤ B

lemma Bd_all (B : ℕ) (v : List ℤ) (e : ℤ →₀ ℕ) (h : Bd B v e) (p : ℤ) : wt v + e p ≤ B := by
  by_cases hp : p ∈ e.support
  · exact (h.2 p hp).2
  · rw [Finsupp.notMem_support_iff.mp hp]; simpa using h.1

/-- The statement for one index word, in base `N - 2`. -/
def Claim (v : List ℤ) : Prop :=
  ∀ B : ℕ, ∀ Q : Polynomial ℚ, ∀ e : ℤ →₀ ℕ, Bd B v e →
    Good B (fun N => ratEval Q e (N : ℚ) * Z v (N - 2))

lemma lower (v : List ℤ) (hv : IsIndex v) (IH : ∀ t < v.length, Claim (v.drop (t + 1)))
    (B τ : ℕ) (Q : Polynomial ℚ) (e : ℤ →₀ ℕ) (hb : Bd B v e) :
    Good B (fun N => ratEval Q e (N : ℚ) * Z v (N - 2 + τ) - ratEval Q e (N : ℚ) * Z v (N - 2)) := by
  choose Qt et het hWt using fun t : ℕ => proof7b v hv τ (t + 1)
  refine Good_congr B _ _ (Good_finsum B (Finset.range v.length)
    (fun t N => ratEval (Q * Qt t) (e + et t) (N : ℚ) * Z (v.drop (t + 1)) (N - 2))
    (fun t ht => IH t (Finset.mem_range.mp ht) B _ _ ⟨?_, ?_⟩)) (fun N hN he => ?_)
  · have := wt_take_drop v (t + 1)
    have := hb.1
    omega
  · intro p hp
    have hsupp := Finsupp.support_add hp
    have hp1 : p ≤ 1 := by
      rcases Finset.mem_union.mp hsupp with h | h
      · exact (hb.2 p h).1
      · exact (het t p h).1
    refine ⟨hp1, ?_⟩
    have h1 := Bd_all B v e hb p
    have h2 : et t p ≤ wt (v.take (t + 1)) := by
      by_cases hq : p ∈ (et t).support
      · exact (het t p hq).2
      · rw [Finsupp.notMem_support_iff.mp hq]; exact Nat.zero_le _
    have := wt_take_drop v (t + 1)
    rw [Finsupp.add_apply]
    omega
  · rw [proof7a v hv τ N hN he, Finset.sum_range_succ', List.take_zero, List.drop_zero]
    rw [show Wc [] ((τ : ℤ) - 1) (N : ℚ) = 1 from rfl, one_mul, mul_add, add_sub_cancel_right,
      Finset.mul_sum]
    refine Finset.sum_congr rfl (fun t _ => ?_)
    rw [hWt t N (by exact_mod_cast hN), ← mul_assoc, ratEval_mul]

lemma shift_down (v : List ℤ) (hv : IsIndex v) (IH : ∀ t < v.length, Claim (v.drop (t + 1)))
    (B τ : ℕ) (Q : Polynomial ℚ) (e : ℤ →₀ ℕ) (hb : Bd B v e)
    (h : Good B (fun N => ratEval Q e (N : ℚ) * Z v (N - 2 + τ))) :
    Good B (fun N => ratEval Q e (N : ℚ) * Z v (N - 2)) :=
  Good_congr B _ _ (Good_sub B _ _ h (lower v hv IH B τ Q e hb)) (fun N _ _ => by ring)

lemma shift_up (v : List ℤ) (hv : IsIndex v) (IH : ∀ t < v.length, Claim (v.drop (t + 1)))
    (B τ : ℕ) (Q : Polynomial ℚ) (e : ℤ →₀ ℕ) (hb : Bd B v e)
    (h : Good B (fun N => ratEval Q e (N : ℚ) * Z v (N - 2))) :
    Good B (fun N => ratEval Q e (N : ℚ) * Z v (N - 2 + τ)) :=
  Good_congr B _ _ (Good_add B _ _ h (lower v hv IH B τ Q e hb)) (fun N _ _ => by ring)

lemma prod_split (a : ℤ) (j : ℕ) (f : ℤ →₀ ℕ) (N : ℚ) :
    (Finsupp.single a j + f).prod (fun p k => (N - (p : ℚ)) ^ k) =
      (N - (a : ℚ)) ^ j * f.prod (fun p k => (N - (p : ℚ)) ^ k) := by
  rw [Finsupp.prod_add_index' (fun _ => pow_zero _) (fun _ _ _ => pow_add _ _ _)]
  congr 1
  exact Finsupp.prod_single_index (pow_zero _)

/-- One step of the partial fraction decomposition. -/
lemma pf_split (a : ℤ) (k : ℕ) (f : ℤ →₀ ℕ) (haf : a ∉ f.support)
    (hf1 : ∀ p ∈ f.support, p ≤ 1) (ha1 : a ≤ 1) (Q : Polynomial ℚ) :
    ∃ c0 : ℚ, ∃ Q2 : Polynomial ℚ, ∀ N : ℚ, 2 ≤ N →
      ratEval Q (Finsupp.single a (k + 1) + f) N =
        c0 * ratEval 1 (Finsupp.single a (k + 1)) N + ratEval Q2 (Finsupp.single a k + f) N := by
  set E : Polynomial ℚ := ∏ p ∈ f.support, (X - C (p : ℚ)) ^ (f p) with hE
  have hEN : ∀ N : ℚ, E.eval N = f.prod (fun p k => (N - (p : ℚ)) ^ k) := by
    intro N
    simp [hE, Finsupp.prod, eval_prod]
  have hEa : E.eval (a : ℚ) ≠ 0 := by
    rw [hEN, Finsupp.prod, Finset.prod_ne_zero_iff]
    intro p hp
    apply pow_ne_zero
    have : a ≠ p := fun h => haf (h ▸ hp)
    rw [sub_ne_zero]
    exact_mod_cast this
  set c0 := Q.eval (a : ℚ) / E.eval (a : ℚ) with hc0
  have hroot : (Q - C c0 * E).IsRoot (a : ℚ) := by
    simp only [IsRoot, eval_sub, eval_mul, eval_C, hc0]
    field_simp
    ring
  set Q2 := (Q - C c0 * E) /ₘ (X - C (a : ℚ))
  have hQ : (X - C (a : ℚ)) * Q2 = Q - C c0 * E := mul_divByMonic_eq_iff_isRoot.mpr hroot
  refine ⟨c0, Q2, fun N hN => ?_⟩
  have hQN := congrArg (eval N) hQ
  simp only [eval_mul, eval_sub, eval_X, eval_C] at hQN
  have hEN0 : E.eval N ≠ 0 := by
    rw [hEN, Finsupp.prod, Finset.prod_ne_zero_iff]
    intro p hp
    apply pow_ne_zero
    have : (p : ℚ) ≤ 1 := by exact_mod_cast hf1 p hp
    linarith
  have hNa : N - (a : ℚ) ≠ 0 := by
    have : (a : ℚ) ≤ 1 := by exact_mod_cast ha1
    linarith
  unfold ratEval
  rw [prod_split, prod_split, show (Finsupp.single a (k + 1)).prod (fun p k => (N - (p : ℚ)) ^ k)
    = (N - (a : ℚ)) ^ (k + 1) from Finsupp.prod_single_index (pow_zero _), eval_one, ← hEN]
  have hq : Q.eval N = c0 * E.eval N + (N - (a : ℚ)) * Q2.eval N := by linarith
  rw [hq]
  field_simp
  ring

lemma stepA (v : List ℤ) (hv : IsIndex v) (IH : ∀ t < v.length, Claim (v.drop (t + 1)))
    (B : ℕ) (Q : Polynomial ℚ) (hB : wt v ≤ B) :
    Good B (fun N => ratEval Q 0 (N : ℚ) * Z v (N - 2)) := by
  have hi : ∀ i : ℕ, Good B (fun N => ratEval (X ^ i) 0 (N : ℚ) * Z v (N - 2)) := fun i =>
    shift_down v hv IH B 1 (X ^ i) 0 ⟨hB, by simp⟩ (Good_mono _ _ hB _ (block_gamma v hv i))
  refine Good_congr B _ _ (Good_finsum B (Finset.range (Q.natDegree + 1))
    (fun i N => Q.coeff i * (ratEval (X ^ i) 0 (N : ℚ) * Z v (N - 2)))
    (fun i _ => Good_smul _ _ _ (hi i))) (fun N _ _ => ?_)
  unfold ratEval
  simp only [Finsupp.prod_zero_index, div_one, eval_pow, eval_X]
  rw [eval_eq_sum_range, Finset.sum_mul]
  refine Finset.sum_congr rfl (fun i _ => by ring)

lemma stepB (v : List ℤ) (hv : IsIndex v) (IH : ∀ t < v.length, Claim (v.drop (t + 1)))
    (B : ℕ) (p : ℤ) (κ : ℕ) (hp : p ≤ 1) (hκ : 1 ≤ κ) (hb : wt v + κ ≤ B) :
    Good B (fun N => ratEval 1 (Finsupp.single p κ) (N : ℚ) * Z v (N - 2)) := by
  obtain ⟨m, rfl⟩ : ∃ m : ℕ, p = 1 - m := ⟨(1 - p).toNat, by omega⟩
  refine shift_down v hv IH B m 1 _ ⟨by omega, fun q hq => ?_⟩
    (Good_mono _ _ (by omega) _ (block_beta v hv κ m hκ))
  have := Finsupp.support_single_subset hq
  rw [Finset.mem_singleton] at this
  subst this
  rw [Finsupp.single_eq_same]
  exact ⟨hp, by omega⟩

lemma stepC (v : List ℤ) (hv : IsIndex v) (IH : ∀ t < v.length, Claim (v.drop (t + 1))) :
    Claim v := by
  intro B Q e
  revert Q
  refine Finsupp.induction e ?_ ?_
  · intro Q hb
    exact stepA v hv IH B Q hb.1
  · intro a b f haf hb0 ihf Q hb
    have hfa : f a = 0 := Finsupp.notMem_support_iff.mp haf
    have hval : (Finsupp.single a b + f) a = b := by
      rw [Finsupp.add_apply, Finsupp.single_eq_same, hfa, add_zero]
    have hamem : a ∈ (Finsupp.single a b + f).support := by
      rw [Finsupp.mem_support_iff, hval]; exact hb0
    obtain ⟨ha1, hab⟩ := hb.2 a hamem
    rw [hval] at hab
    have hBf : Bd B v f := by
      refine ⟨hb.1, fun p hp => ?_⟩
      have hpa : p ≠ a := fun h => haf (h ▸ hp)
      have hv' : (Finsupp.single a b + f) p = f p := by
        rw [Finsupp.add_apply, Finsupp.single_eq_of_ne hpa, zero_add]
      have hmem : p ∈ (Finsupp.single a b + f).support := by
        rw [Finsupp.mem_support_iff, hv']; exact Finsupp.mem_support_iff.mp hp
      have := hb.2 p hmem
      rw [hv'] at this
      exact this
    have hf1 : ∀ p ∈ f.support, p ≤ 1 := fun p hp => (hBf.2 p hp).1
    suffices H : ∀ b' : ℕ, b' ≤ b → ∀ Q : Polynomial ℚ,
        Good B (fun N => ratEval Q (Finsupp.single a b' + f) (N : ℚ) * Z v (N - 2)) from
      H b le_rfl Q
    intro b'
    induction b' with
    | zero =>
      intro _ Q
      rw [Finsupp.single_zero, zero_add]
      exact ihf Q hBf
    | succ b' ihb =>
      intro hle Q
      obtain ⟨c0, Q2, hsplit⟩ := pf_split a b' f haf hf1 ha1 Q
      refine Good_congr B _ _ (Good_add B _ _
        (Good_smul B c0 _ (stepB v hv IH B a (b' + 1) ha1 (by omega) (by omega)))
        (ihb (by omega) Q2)) (fun N hN _ => ?_)
      rw [hsplit N (by exact_mod_cast hN)]
      ring

lemma claim_all : ∀ n : ℕ, ∀ v : List ℤ, v.length = n → IsIndex v → Claim v := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    intro v hlen hv
    refine stepC v hv (fun t ht => ?_)
    exact ih _ (by rw [List.length_drop]; omega) _ rfl
      (fun b hb => hv b (List.mem_of_mem_drop hb))

lemma term_good (B : ℕ) (τ : DTerm) (h : IsIndex τ.2.2 ∧ wt τ.2.2 ≤ B ∧
    ∀ p ∈ τ.2.1.support, p ≤ 1 ∧ wt τ.2.2 + τ.2.1 p ≤ B) :
    Good B (fun N => ratEval τ.1 τ.2.1 (N : ℚ) * Z τ.2.2 N) := by
  have hc := claim_all _ τ.2.2 rfl h.1
  have hb : Bd B τ.2.2 τ.2.1 := ⟨h.2.1, h.2.2⟩
  have := shift_up τ.2.2 h.1 (fun t _ => claim_all _ _ rfl
    (fun b hb => h.1 b (List.mem_of_mem_drop hb))) B 2 τ.1 τ.2.1 hb (hc B τ.1 τ.2.1 hb)
  refine Good_congr B _ _ this (fun N hN _ => ?_)
  rw [Nat.sub_add_cancel hN]

theorem main : S08.stmt_8 := by
  intro B T hT
  obtain ⟨a, ⟨F, hF, q, hq⟩, ha⟩ := Good_list B T
    (fun τ N => ratEval τ.1 τ.2.1 (N : ℚ) * Z τ.2.2 N) (fun τ hτ => term_good B τ (hT τ hτ))
  exact ⟨F, hF, q, a, hq, ha⟩

end S08.PB

-- ===== merged from work/S08/PartC.lean =====
namespace S08.PC

open Real Polynomial

lemma InG_add (B : ℕ) (F G : ℝ → ℝ) (hF : S03.InG B F) (hG : S03.InG B G) :
    S03.InG B (fun x => F x + G x) := by
  obtain ⟨L1, h1, e1⟩ := hF
  obtain ⟨L2, h2, e2⟩ := hG
  refine ⟨L1 ++ L2, ?_, ?_⟩
  · intro g hg
    rcases List.mem_append.mp hg with h | h
    · exact h1 g h
    · exact h2 g h
  · intro x hx0 hx1
    show F x + G x = _
    rw [e1 x hx0 hx1, e2 x hx0 hx1, List.map_append, List.sum_append]

lemma kval_C_mul (c : ℚ) (p : Polynomial ℚ) (a b d : ℕ) (x : ℝ) :
    S03.kval (C c * p) a b d x = (c : ℝ) * S03.kval p a b d x := by
  simp only [S03.kval, map_mul, aeval_C, eq_ratCast]
  ring

lemma InG_smul (B : ℕ) (c : ℚ) (F : ℝ → ℝ) (hF : S03.InG B F) :
    S03.InG B (fun x => (c : ℝ) * F x) := by
  obtain ⟨L, hL, e⟩ := hF
  refine ⟨L.map (fun g => (C c * g.1, g.2)), ?_, ?_⟩
  · intro g hg
    rw [List.mem_map] at hg
    obtain ⟨g0, hg0, rfl⟩ := hg
    exact hL g0 hg0
  · intro x hx0 hx1
    show (c : ℝ) * F x = _
    rw [e x hx0 hx1, List.map_map, ← List.sum_map_mul_left]
    congr 1
    refine List.map_congr_left (fun g _ => ?_)
    simp only [Function.comp, kval_C_mul]
    ring

lemma InG_mono (B : ℕ) (c : ℚ) (i : ℕ) :
    S03.InG B (fun x => Polynomial.aeval (1 / x) (C c * X ^ i)) := by
  refine ⟨[(C c, i, 0, 0, [])], ?_, ?_⟩
  · intro g hg
    simp only [List.mem_singleton] at hg
    subst hg
    refine ⟨?_, ?_, ?_⟩
    · intro a ha; simp at ha
    · intro a ha; simp at ha
    · simp
  · intro x _ _
    simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, S03.kval, S03.H,
      map_mul, aeval_C, map_pow, aeval_X, eq_ratCast, pow_zero, mul_one, add_zero, div_pow]
    ring

lemma InG_poly (B : ℕ) (q : Polynomial ℚ) :
    S03.InG B (fun x => Polynomial.aeval (1 / x) q) := by
  induction q using Polynomial.induction_on' with
  | add p r hp hr =>
    have := InG_add B _ _ hp hr
    simpa only [map_add] using this
  | monomial i c =>
    have := InG_mono B c i
    simpa only [← Polynomial.C_mul_X_pow_eq_monomial] using this

/-- The even part with the Laurent tail removed. -/
noncomputable def Gs (F : ℝ → ℝ) (q : Polynomial ℚ) (a0 : ℚ) (x : ℝ) : ℝ :=
  ((1 / 2 : ℚ) : ℝ) * F x + ((1 / 2 : ℚ) : ℝ) * F (-x) +
    ((-1 / 2 : ℚ) : ℝ) * Polynomial.aeval (1 / x) q +
    ((-1 / 2 : ℚ) : ℝ) * Polynomial.aeval (1 / -x) q + ((-a0 : ℚ) : ℝ)

lemma InG_Gs (B : ℕ) (F : ℝ → ℝ) (hF : S03.InG B F) (q : Polynomial ℚ) (a0 : ℚ) :
    S03.InG B (Gs F q a0) := by
  have h1 := InG_smul B (1 / 2) F hF
  have h2 := InG_smul B (1 / 2) _ (S03.InG_neg B F hF)
  have h3 := InG_smul B (-1 / 2) _ (InG_poly B q)
  have h4 := InG_smul B (-1 / 2) _ (S03.InG_neg B _ (InG_poly B q))
  have h5 : S03.InG B (fun x => Polynomial.aeval (1 / x) (C (-a0))) := InG_poly B _
  simp only [aeval_C, eq_ratCast] at h5
  exact InG_add B _ _ (InG_add B _ _ (InG_add B _ _ (InG_add B _ _ h1 h2) h3) h4) h5

lemma even_part_coef (a : ℕ → ℚ) (x : ℝ) (N : ℕ) :
    ((a N : ℝ) * x ^ N + (a N : ℝ) * (-x) ^ N) / 2 =
      if Even N then (a N : ℝ) * x ^ N else 0 := by
  split_ifs with h
  · rw [h.neg_pow]; ring
  · rw [Nat.not_even_iff_odd] at h
    rw [h.neg_pow]; ring

lemma series_Gs (F : ℝ → ℝ) (q : Polynomial ℚ) (a : ℕ → ℚ)
    (hser : ∀ x : ℝ, 0 < |x| → |x| < 1 →
      HasSum (fun N : ℕ => (a N : ℝ) * x ^ N) (F x - Polynomial.aeval (1 / x) q))
    (x : ℝ) (hx0 : 0 < |x|) (hx1 : |x| < 1) :
    HasSum (fun n : ℕ => (a (2 * (n + 1)) : ℝ) * x ^ (2 * (n + 1))) (Gs F q (a 0) x) := by
  have hA := hser x hx0 hx1
  have hB := hser (-x) (by rwa [abs_neg]) (by rwa [abs_neg])
  have hC := (hA.add hB).div_const 2
  set f : ℕ → ℝ := fun N => if Even N then (a N : ℝ) * x ^ N else 0 with hf
  have hC' : HasSum f ((F x - Polynomial.aeval (1 / x) q +
      (F (-x) - Polynomial.aeval (1 / -x) q)) / 2) := by
    refine hC.congr_fun ?_
    intro N
    simp only [hf]
    exact (even_part_coef a x N).symm
  have hD := hC'.update 0 0
  have hinj : Function.Injective (fun n : ℕ => 2 * (n + 1)) := by
    intro m n h; simp only at h; omega
  have hE := (hinj.hasSum_iff (f := Function.update f 0 0) ?_).mpr hD
  · convert hE using 1
    · funext n
      simp only [Function.comp, Function.update, hf]
      rw [dif_neg (by omega), if_pos (by exact ⟨n + 1, by ring⟩)]
    · simp only [hf, Gs, pow_zero, mul_one]
      rw [if_pos (by exact ⟨0, rfl⟩)]
      push_cast
      ring
  · intro N hN
    simp only [Set.mem_range, not_exists] at hN
    by_cases h0 : N = 0
    · subst h0; simp
    · rw [Function.update_of_ne h0]
      simp only [hf]
      rw [if_neg]
      rintro ⟨m, hm⟩
      exact hN (m - 1) (by omega)

theorem main : S08.stmt_4 → S08.stmt_8 → S08.stmt_9 := by
  intro h4 h8 k s hs hw
  obtain ⟨T, hT, hTP⟩ := h4 k s hs hw
  obtain ⟨F, hF, q, a, hser, ha⟩ := h8 _ T hT
  refine ⟨Gs F q (a 0), InG_Gs _ F hF q (a 0), ?_⟩
  intro x hx0 hx1
  refine (series_Gs F q a hser x hx0 hx1).congr_fun ?_
  intro n
  rw [hTP (n + 1) (by omega), ha (2 * (n + 1)) (by omega) ⟨n + 1, by ring⟩]

end S08.PC

-- ===== merged from work/S08/PartD.lean =====
namespace S08.PD

open Real

/-- The generating series value `Σ_{n≥1} [𝔷^s]P_{k,n} x^{2n}` as a `tsum`. -/
noncomputable def G (k : ℕ) (s : ℕ →₀ ℕ) (x : ℝ) : ℝ :=
  ∑' n : ℕ, (((Q04.Pf k (n + 1)).coeff s : ℚ) : ℝ) * x ^ (2 * (n + 1))

/-- The admissibility condition on a monomial. -/
def Cond (k : ℕ) (s : ℕ →₀ ℕ) : Prop :=
  (∀ w ∈ s.support, 2 ≤ w ∧ w ≤ k) ∧ S03.mwt s ≤ k

/-- The value of a final term. -/
noncomputable def fv (x : ℝ) (g : S03.FTerm) : ℝ :=
  S03.zmon g.1 * S03.kval g.2.1 g.2.2.1 g.2.2.2.1 g.2.2.2.2.1 x * S03.H g.2.2.2.2.2 x

/-- The properties of a final term. -/
def FProp (k : ℕ) (g : S03.FTerm) : Prop :=
  (∀ w ∈ g.1.support, 2 ≤ w ∧ w ≤ k) ∧ S03.IsLetters g.2.2.2.2.2 ∧
    S03.Regular g.2.2.2.2.2 ∧ S03.mwt g.1 + g.2.2.2.2.2.length ≤ k + 1

lemma sum_map_pair (s : ℕ →₀ ℕ) (x : ℝ) (Ls : List S03.GTerm) :
    ((Ls.map fun g => ((s, g) : S03.FTerm)).map (fv x)).sum =
      S03.zmon s * (Ls.map fun g => S03.kval g.1 g.2.1 g.2.2.1 g.2.2.2.1 x * S03.H g.2.2.2.2 x).sum := by
  induction Ls with
  | nil => simp
  | cons g Ls ih =>
    simp only [List.map_cons, List.sum_cons] at ih ⊢
    rw [ih]; unfold fv; ring

lemma collect (h9 : S08.stmt_9) (k : ℕ) : ∀ l : List (ℕ →₀ ℕ), (∀ s ∈ l, Cond k s) →
    ∃ L : List S03.FTerm, (∀ g ∈ L, FProp k g) ∧
      ∀ x : ℝ, 0 < |x| → |x| < 1 →
        (l.map fun s => S03.zmon s * G k s x).sum = (L.map (fv x)).sum := by
  intro l
  induction l with
  | nil => intro _; exact ⟨[], by simp, by simp⟩
  | cons s l ih =>
    intro hl
    obtain ⟨L', hL'p, hL'v⟩ := ih (fun t ht => hl t (List.mem_cons_of_mem _ ht))
    have hs := hl s List.mem_cons_self
    obtain ⟨F, ⟨Ls, hLsp, hLsv⟩, hF⟩ := h9 k s hs.1 hs.2
    refine ⟨(Ls.map fun g => ((s, g) : S03.FTerm)) ++ L', ?_, ?_⟩
    · intro g hg
      rcases List.mem_append.mp hg with hg | hg
      · obtain ⟨g', hg', rfl⟩ := List.mem_map.mp hg
        obtain ⟨h1, h2, h3⟩ := hLsp g' hg'
        refine ⟨hs.1, h1, h2, ?_⟩
        show S03.mwt s + g'.2.2.2.2.length ≤ k + 1
        have := hs.2
        omega
      · exact hL'p g hg
    · intro x hx0 hx1
      rw [List.map_append, List.sum_append, ← hL'v x hx0 hx1, sum_map_pair,
        List.map_cons, List.sum_cons]
      congr 1
      rw [← hLsv x hx0 hx1]
      unfold G
      rw [(hF x hx0 hx1).tsum_eq]

/-- The finite set of admissible monomials. -/
noncomputable def Sk (k : ℕ) : Finset (ℕ →₀ ℕ) :=
  ((Finset.Icc 2 k).finsupp fun _ => Finset.range (k + 1)).filter fun s => S03.mwt s ≤ k

lemma Sk_cond (k : ℕ) (s : ℕ →₀ ℕ) (hs : s ∈ Sk k) : Cond k s := by
  unfold Sk at hs
  rw [Finset.mem_filter, Finset.mem_finsupp_iff] at hs
  refine ⟨fun w hw => ?_, hs.2⟩
  have := hs.1.1 hw
  simpa using this

lemma single_le_mwt (s : ℕ →₀ ℕ) (w : ℕ) (hw : w ∈ s.support) : s w * w ≤ S03.mwt s := by
  unfold S03.mwt Finsupp.sum
  exact Finset.single_le_sum (f := fun w => s w * w) (fun _ _ => Nat.zero_le _) hw

lemma supp_sub (k n : ℕ) (s : ℕ →₀ ℕ) (hs : s ∈ (Q04.Pf k n).support) : s ∈ Sk k := by
  obtain ⟨h2, hw⟩ := Q04.proof_5.1 k n s hs
  have hw' : S03.mwt s ≤ k := hw
  unfold Sk
  rw [Finset.mem_filter, Finset.mem_finsupp_iff]
  refine ⟨⟨fun w hws => ?_, fun i hi => ?_⟩, hw'⟩
  · have := single_le_mwt s w hws
    have h1 : 1 ≤ s w := Nat.one_le_iff_ne_zero.mpr (Finsupp.mem_support_iff.mp hws)
    have := h2 w hws
    rw [Finset.mem_Icc]
    constructor
    · exact this
    · nlinarith
  · rw [Finset.mem_range]
    by_cases h0 : s i = 0
    · omega
    · have hws : i ∈ s.support := Finsupp.mem_support_iff.mpr h0
      have := single_le_mwt s i hws
      have := h2 i hws
      nlinarith

lemma aeval_eq (k n : ℕ) :
    MvPolynomial.aeval X05.z (Q04.Pf k n) =
      ∑ s ∈ Sk k, (((Q04.Pf k n).coeff s : ℚ) : ℝ) * S03.zmon s := by
  rw [MvPolynomial.aeval_def, MvPolynomial.eval₂_eq]
  apply Finset.sum_subset (fun s hs => supp_sub k n s hs)
  · intro s _ hs
    rw [MvPolynomial.notMem_support_iff.mp hs]
    simp

lemma c_pos (α : ℝ) (h0 : 0 < α) (h1 : α < π) : 0 < D04.c α ∧ D04.c α < 1 := by
  unfold D04.c
  constructor
  · apply Real.cos_pos_of_mem_Ioo; constructor <;> linarith
  · have : Real.cos (α / 2) < Real.cos 0 := by
      apply Real.cos_lt_cos_of_nonneg_of_le_pi_div_two (le_refl 0) (by linarith) (by linarith)
    simpa using this

theorem main : S08.stmt_9 → S08.stmt_11 := by
  intro h9 k
  obtain ⟨L, hLp, hLv⟩ := collect h9 k (Sk k).toList
    (fun s hs => Sk_cond k s (Finset.mem_toList.mp hs))
  refine ⟨L, hLp, ?_⟩
  intro α h0 h1
  obtain ⟨hc0, hc1⟩ := c_pos α h0 h1
  have hx0 : 0 < |D04.c α| := by rw [abs_of_pos hc0]; exact hc0
  have hx1 : |D04.c α| < 1 := by rw [abs_of_pos hc0]; exact hc1
  have key := (Q04.proof_5.2 α h0 h1 k).2
  have hsum : HasSum (fun n : ℕ => MvPolynomial.aeval X05.z (Q04.Pf k (n + 1)) * D04.c α ^ (2 * (n + 1)))
      (∑ s ∈ Sk k, S03.zmon s * G k s (D04.c α)) := by
    have : (fun n : ℕ => MvPolynomial.aeval X05.z (Q04.Pf k (n + 1)) * D04.c α ^ (2 * (n + 1))) =
        fun n => ∑ s ∈ Sk k, S03.zmon s *
          ((((Q04.Pf k (n + 1)).coeff s : ℚ) : ℝ) * D04.c α ^ (2 * (n + 1))) := by
      funext n
      rw [aeval_eq, Finset.sum_mul]
      apply Finset.sum_congr rfl
      intro s _; ring
    rw [this]
    apply hasSum_sum
    intro s hs
    have hc := Sk_cond k s hs
    obtain ⟨F, _, hF⟩ := h9 k s hc.1 hc.2
    have hF' := hF (D04.c α) hx0 hx1
    unfold G
    rw [hF'.tsum_eq]
    exact hF'.mul_left _
  rw [key.unique hsum, ← Finset.sum_map_toList, hLv _ hx0 hx1]
  rfl

end S08.PD

theorem S08.proof_4 : S08.stmt_4 := S08.PA.main

theorem S08.proof_8 : S08.stmt_8 := S08.PB.main

theorem S08.proof_9 : S08.stmt_9 := S08.PC.main S08.PA.main S08.PB.main

theorem S08.proof_11 : S08.stmt_11 := S08.PD.main (S08.PC.main S08.PA.main S08.PB.main)
