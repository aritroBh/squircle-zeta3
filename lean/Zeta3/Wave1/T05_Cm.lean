import Mathlib
/-!
# T05 — the constant `C_m` of Theorem `thm:Cm` (regular `2m`-gons)

Target id: T05.
Paper label: `thm:Cm`, paper/ieee/main.tex lines 1657–1663
(proof supplement: paper/ieee/proofs/sec_g4_Cm.tex, label `g4:thm:Cm`).

Claim quoted verbatim (main.tex lines 1657–1663):

  \begin{theorem}[Regular $2m$-gons$^\dagger$]\label{thm:Cm}
  For the smoothed regular $2m$-gon, the $\zeta(3)$ part of what the mixed terms
  contribute to $P'''(1)$ is exactly $C_m\,\zeta(3)$, where
  \[
  C_m=48m\sin^4\frac{\pi}{2m}\,\cos\frac{\pi}{2m}\Bigl(3\sin^2\frac{\pi}{2m}-1\Bigr).
  \]
  So $C_2=6\sqrt2$, $C_3=-\frac{9\sqrt3}{8}$, and $C_m<0$ for every $m\ge3$.
  \end{theorem}

and the linking sentence after Proposition `prop:edgelaw` (main.tex line 1686):
"For the regular $2m$-gon $K_1$ of Theorem~\ref{thm:Cm}, $\sum_eZ_e=C_m$", where
(main.tex, Proposition `prop:edgelaw`)
  $Z_e=6\,\xi_0\xi_1\frac{h^2}{L}\Bigl(3\,\frac{\xi_0L_0+\xi_1L_1}{L}-2\Bigr)$,
with $L$ the edge length, $h$ its distance from the center, $L_0$, $L_1=L-L_0$ the
signed lengths from $V_0$, $V_1$ to the foot of the perpendicular from the center,
$\xi_0=|\langle v_{j_0},V_1\rangle|$, $\xi_1=|\langle v_{j_1},V_0\rangle|$, and
$\langle v_{j_0},x\rangle=0$, $\langle v_{j_1},x\rangle=0$ the kink lines through
$V_0$, $V_1$. The kink normals are (main.tex line 1646)
$n_i=(\cos\frac{i\pi}m,\sin\frac{i\pi}m)$, $i=1,\dots,m$, and
$K_1=\{\sum_{i=1}^m|\langle n_i,x\rangle|\le1\}$.

Formalization choices:
* `Cm m` is the closed formula, a real number for `m : ℕ` (cast to `ℝ`).
  `m ≥ 2` is the paper's range (`2m`-gon); for `m = 0` the formula is junk but no
  statement uses it.
* The corners of `K_1` are taken as in the proof supplement (sec_g4_Cm.tex):
  `C_k = ρ₀ · (-sin(kπ/m), cos(kπ/m))`, `ρ₀ = tan(π/(2m))`
  (so `C_0 = ρ₀(0,1)` and `C_1 = ρ₀(-sin(π/m), cos(π/m))`). `stmt_5` checks that
  these corners lie on `∂K_1` (`Σ_i |⟨n_i, C_k⟩| = 1`).
* `Z_e` for the edge `C_k C_{k+1}` is computed from the primitive geometry:
  `L = |C_{k+1} - C_k|`, `h = |det(C_k, C_{k+1})| / L` (distance from the origin to
  the line of the edge), `L_0 = ⟨-C_k, C_{k+1} - C_k⟩ / L` (signed), `L_1 = L - L_0`.
  The kink-line indices `j_0, j_1 ∈ {1,…,m}` are quantified with the hypotheses
  `⟨n_{j_0}, C_k⟩ = 0` and `⟨n_{j_1}, C_{k+1}⟩ = 0` (for each corner exactly one such
  index exists, so the hypotheses are not vacuous). `stmt_4` says every edge has
  `2m · Z_e = C_m`; since there are `2m` edges, this is `Σ_e Z_e = C_m`.
* Divisions by `L` are harmless: `L > 0` for `m ≥ 2`.

NOT COVERED:
* The analytic content of `thm:Cm`: that `C_m ζ(3)` is the `ζ(3)` part of the
  contribution of the mixed local terms (`ℓ_0^2 ℓ_1`, `ℓ_0 ℓ_1^2`) of the `d^3`
  term of the length element of `∂K_p` to `P'''(1)`. This needs the perimeter
  `P(p)` of `K_p`, its third derivative at `p = 1`, the edge-wise decomposition into
  local / mixed / far terms, and a notion of "the `ζ(3)` part" of a number, none of
  which can be stated faithfully in Mathlib today. Only the closed form, the special
  values, the sign, and the identity `Σ_e Z_e = C_m` (Proposition `prop:edgelaw`'s
  per-edge coefficient summed over the regular polygon) are formalized.
* That the points `C_k` are exactly the vertices of `K_1` (only that they lie on
  `∂K_1` is stated), and Proposition `prop:edgelaw` itself.
-/
-- BEGIN STATEMENTS
namespace T05

/-- `s = sin(π/(2m))`. -/
noncomputable def s (m : ℕ) : ℝ := Real.sin (Real.pi / (2 * (m : ℝ)))

/-- `c = cos(π/(2m))`. -/
noncomputable def c (m : ℕ) : ℝ := Real.cos (Real.pi / (2 * (m : ℝ)))

/-- The paper's constant
`C_m = 48 m sin^4(π/(2m)) cos(π/(2m)) (3 sin^2(π/(2m)) - 1)`. -/
noncomputable def Cm (m : ℕ) : ℝ :=
  48 * (m : ℝ) * s m ^ 4 * c m * (3 * s m ^ 2 - 1)

/-- Euclidean inner product on `ℝ × ℝ`. -/
def ip (a b : ℝ × ℝ) : ℝ := a.1 * b.1 + a.2 * b.2

/-- Kink normal `n_i = (cos(iπ/m), sin(iπ/m))`. -/
noncomputable def n (m i : ℕ) : ℝ × ℝ :=
  (Real.cos ((i : ℝ) * Real.pi / m), Real.sin ((i : ℝ) * Real.pi / m))

/-- `ρ₀ = tan(π/(2m))`, the circumradius of `K_1`. -/
noncomputable def rho0 (m : ℕ) : ℝ := Real.tan (Real.pi / (2 * (m : ℝ)))

/-- Corner `C_k = ρ₀ (-sin(kπ/m), cos(kπ/m))` of the regular `2m`-gon `K_1`. -/
noncomputable def corner (m k : ℕ) : ℝ × ℝ :=
  (rho0 m * (-Real.sin ((k : ℝ) * Real.pi / m)), rho0 m * Real.cos ((k : ℝ) * Real.pi / m))

/-- Edge length `L = |V₁ - V₀|`. -/
noncomputable def edgeL (V0 V1 : ℝ × ℝ) : ℝ :=
  Real.sqrt ((V1.1 - V0.1) ^ 2 + (V1.2 - V0.2) ^ 2)

/-- Distance `h` from the center to the line through `V₀, V₁`. -/
noncomputable def edgeH (V0 V1 : ℝ × ℝ) : ℝ :=
  |V0.1 * V1.2 - V0.2 * V1.1| / edgeL V0 V1

/-- Signed length `L₀` from `V₀` to the foot of the perpendicular from the center. -/
noncomputable def edgeL0 (V0 V1 : ℝ × ℝ) : ℝ :=
  ip (-V0) (V1 - V0) / edgeL V0 V1

/-- Proposition `prop:edgelaw`'s coefficient
`Z_e = 6 ξ₀ ξ₁ (h²/L) (3 (ξ₀ L₀ + ξ₁ L₁)/L - 2)` for the edge `V₀V₁` with kink
normals `v₀` (kink line through `V₀`) and `v₁` (kink line through `V₁`). -/
noncomputable def Ze (V0 V1 v0 v1 : ℝ × ℝ) : ℝ :=
  let L := edgeL V0 V1
  let h := edgeH V0 V1
  let L0 := edgeL0 V0 V1
  let L1 := L - L0
  let xi0 := |ip v0 V1|
  let xi1 := |ip v1 V0|
  6 * xi0 * xi1 * (h ^ 2 / L) * (3 * (xi0 * L0 + xi1 * L1) / L - 2)

/-- "So $C_2=6\sqrt2$". -/
def stmt_1 : Prop := Cm 2 = 6 * Real.sqrt 2

/-- "$C_3=-\frac{9\sqrt3}{8}$". -/
def stmt_2 : Prop := Cm 3 = -(9 * Real.sqrt 3 / 8)

/-- "and $C_m<0$ for every $m\ge3$". -/
def stmt_3 : Prop := ∀ m : ℕ, 3 ≤ m → Cm m < 0

/-- Line 1686: "For the regular $2m$-gon $K_1$ of Theorem~\ref{thm:Cm},
$\sum_eZ_e=C_m$" — every one of the `2m` edges `C_k C_{k+1}` has `2m · Z_e = C_m`,
with `Z_e` computed from the geometry and the kink lines through its corners. -/
def stmt_4 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∀ k j0 j1 : ℕ,
    j0 ∈ Finset.Icc 1 m → j1 ∈ Finset.Icc 1 m →
    ip (n m j0) (corner m k) = 0 → ip (n m j1) (corner m (k + 1)) = 0 →
    2 * (m : ℝ) * Ze (corner m k) (corner m (k + 1)) (n m j0) (n m j1) = Cm m

/-- Setup check (proof supplement, "$\mathcal N=1$ at all vertices"): the corners
`C_k` lie on `∂K_1 = {Σ_{i=1}^m |⟨n_i,x⟩| = 1}`. -/
def stmt_5 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∀ k : ℕ,
    ∑ i ∈ Finset.Icc 1 m, |ip (n m i) (corner m k)| = 1

end T05
-- END STATEMENTS
namespace T05

lemma theta_pos {m : ℕ} (hm : 2 ≤ m) : 0 < Real.pi / (2 * (m : ℝ)) := by
  have : (2:ℝ) ≤ m := by exact_mod_cast hm
  have := Real.pi_pos
  positivity

lemma theta_le {m : ℕ} (hm : 2 ≤ m) : Real.pi / (2 * (m : ℝ)) ≤ Real.pi / 4 := by
  have h2 : (2:ℝ) ≤ m := by exact_mod_cast hm
  have := Real.pi_pos
  rw [div_le_div_iff₀ (by positivity) (by norm_num)]
  nlinarith

lemma s_pos {m : ℕ} (hm : 2 ≤ m) : 0 < s m := by
  have := theta_le hm
  have := Real.pi_pos
  exact Real.sin_pos_of_pos_of_lt_pi (theta_pos hm) (by linarith)

lemma c_pos {m : ℕ} (hm : 2 ≤ m) : 0 < c m := by
  have := theta_le hm
  have := theta_pos hm
  have := Real.pi_pos
  exact Real.cos_pos_of_mem_Ioo ⟨by linarith, by linarith⟩

lemma s_le_half {m : ℕ} (hm : 3 ≤ m) : s m ≤ 1 / 2 := by
  have h3 : (3:ℝ) ≤ m := by exact_mod_cast hm
  have hp := Real.pi_pos
  have h2 : 2 ≤ m := by omega
  have hle : Real.pi / (2 * (m : ℝ)) ≤ Real.pi / 6 := by
    rw [div_le_div_iff₀ (by positivity) (by norm_num)]
    nlinarith
  have := theta_pos h2
  have := Real.sin_le_sin_of_le_of_le_pi_div_two (by linarith) (by linarith) hle
  unfold s
  rwa [Real.sin_pi_div_six] at this

end T05

namespace T05

lemma pi_div_m {m : ℕ} (hm : 2 ≤ m) : Real.pi / (m : ℝ) = 2 * (Real.pi / (2 * (m : ℝ))) := by
  have : (m:ℝ) ≠ 0 := by
    have : (2:ℝ) ≤ m := by exact_mod_cast hm
    linarith
  field_simp

lemma rho0_eq {m : ℕ} : rho0 m = s m / c m := by
  unfold rho0 s c; exact Real.tan_eq_sin_div_cos _

lemma ip_corner (m j k : ℕ) :
    ip (n m j) (corner m k) = rho0 m * Real.sin (((j:ℝ) - k) * Real.pi / m) := by
  unfold ip n corner
  simp only
  have : ((j:ℝ) - k) * Real.pi / m = j * Real.pi / m - k * Real.pi / m := by ring
  rw [this, Real.sin_sub]; ring

lemma sum_Icc_range (f : ℕ → ℝ) (m : ℕ) :
    ∑ i ∈ Finset.Icc 1 m, f i = ∑ i ∈ Finset.range m, f (i + 1) := by
  induction m with
  | zero => simp
  | succ m ih =>
    rw [Finset.sum_Icc_succ_top (by omega), Finset.sum_range_succ, ih]

lemma sum_shift (m : ℕ) (G : ℤ → ℝ) (hG : ∀ x : ℤ, G (x + m) = G x) (k : ℕ) :
    ∑ i ∈ Finset.range m, G ((i : ℤ) + 1 - k) = ∑ i ∈ Finset.range m, G ((i : ℤ) + 1) := by
  induction k with
  | zero => simp
  | succ k ih =>
    rw [← ih]
    have h1 := Finset.sum_range_succ' (fun i : ℕ => G ((i : ℤ) - k)) m
    have h2 := Finset.sum_range_succ (fun i : ℕ => G ((i : ℤ) - k)) m
    have h3 := hG (-(k:ℤ))
    simp only [Nat.cast_zero, zero_sub] at h1
    have h4 : G ((m:ℤ) - k) = G (-(k:ℤ)) := by
      have : (m:ℤ) - k = -(k:ℤ) + m := by ring
      rw [this, h3]
    have h5 : ∀ i : ℕ, G (((i:ℕ):ℤ) + 1 - ((k+1 : ℕ) : ℤ)) = G ((i:ℤ) - k) := by
      intro i; congr 1; push_cast; ring
    have h6 : ∀ i : ℕ, G (((i + 1 : ℕ):ℤ) - k) = G ((i:ℤ) + 1 - k) := by
      intro i; congr 1
    simp only [h5]
    simp only [h6] at h1
    linarith

lemma sin_sum (θ : ℝ) (n : ℕ) :
    2 * Real.sin θ * ∑ i ∈ Finset.range n, Real.sin (2 * ((i:ℝ) + 1) * θ)
      = Real.cos θ - Real.cos ((2 * (n:ℝ) + 1) * θ) := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [Finset.sum_range_succ, mul_add, ih]
    have e1 : (2 * ((n+1:ℕ):ℝ) + 1) * θ = 2 * ((n:ℝ) + 1) * θ + θ := by push_cast; ring
    have e2 : (2 * (n:ℝ) + 1) * θ = 2 * ((n:ℝ) + 1) * θ - θ := by ring
    rw [e1, e2, Real.cos_add, Real.cos_sub]
    ring

lemma sin_total {m : ℕ} (hm : 2 ≤ m) :
    ∑ i ∈ Finset.range m, Real.sin (((i:ℝ) + 1) * Real.pi / m) = c m / s m := by
  have hs := s_pos hm
  have h := sin_sum (Real.pi / (2 * (m:ℝ))) m
  have hmne : (m:ℝ) ≠ 0 := by
    have : (2:ℝ) ≤ m := by exact_mod_cast hm
    linarith
  have e : (2 * (m:ℝ) + 1) * (Real.pi / (2 * (m:ℝ))) = Real.pi + Real.pi / (2 * (m:ℝ)) := by
    field_simp
  rw [e, add_comm, Real.cos_add_pi] at h
  have e2 : ∀ i : ℕ, Real.sin (((i:ℝ) + 1) * Real.pi / m)
      = Real.sin (2 * ((i:ℝ) + 1) * (Real.pi / (2 * (m:ℝ)))) := by
    intro i; congr 1; field_simp
  simp only [e2]
  unfold s at hs ⊢
  unfold c
  rw [eq_div_iff hs.ne']
  nlinarith [h]

end T05

namespace T05

lemma mpos {m : ℕ} (hm : 2 ≤ m) : (0:ℝ) < m := by
  have : (2:ℝ) ≤ m := by exact_mod_cast hm
  linarith

lemma sc_sq (m : ℕ) : s m ^ 2 + c m ^ 2 = 1 := by
  unfold s c; exact Real.sin_sq_add_cos_sq _

lemma sin2 (m : ℕ) : Real.sin (2 * (Real.pi / (2 * (m:ℝ)))) = 2 * s m * c m := by
  rw [Real.sin_two_mul]; rfl

lemma cos2 (m : ℕ) : Real.cos (2 * (Real.pi / (2 * (m:ℝ)))) = 1 - 2 * s m ^ 2 := by
  rw [Real.cos_two_mul]
  have := sc_sq m
  have e : Real.cos (Real.pi / (2 * (m:ℝ))) = c m := rfl
  rw [e]; nlinarith

lemma rho0_pos {m : ℕ} (hm : 2 ≤ m) : 0 < rho0 m := by
  rw [rho0_eq]; have := s_pos hm; have := c_pos hm; positivity

lemma rho0_c {m : ℕ} (hm : 2 ≤ m) : rho0 m * c m = s m := by
  rw [rho0_eq]; have := c_pos hm; field_simp

lemma corner_succ {m : ℕ} (hm : 2 ≤ m) (k : ℕ) :
    corner m (k + 1) = (rho0 m * (-Real.sin ((k:ℝ) * Real.pi / m + 2 * (Real.pi / (2 * (m:ℝ))))),
      rho0 m * Real.cos ((k:ℝ) * Real.pi / m + 2 * (Real.pi / (2 * (m:ℝ))))) := by
  have e : (((k + 1 : ℕ) : ℝ)) * Real.pi / m = k * Real.pi / m + 2 * (Real.pi / (2 * (m:ℝ))) := by
    rw [← pi_div_m hm]; push_cast; ring
  unfold corner
  rw [e]

lemma abs_sin_sub (x y : ℝ) (hx : Real.sin x = 0) : |Real.sin (x - y)| = |Real.sin y| := by
  have h := Real.sin_sq_add_cos_sq x
  have hc : |Real.cos x| = 1 := by
    have h1 : Real.cos x ^ 2 = 1 := by rw [hx] at h; linarith
    nlinarith [abs_nonneg (Real.cos x), sq_abs (Real.cos x)]
  rw [Real.sin_sub, hx, zero_mul, zero_sub, abs_neg, abs_mul, hc, one_mul]

lemma abs_sin_add (x y : ℝ) (hx : Real.sin x = 0) : |Real.sin (x + y)| = |Real.sin y| := by
  have h := Real.sin_sq_add_cos_sq x
  have hc : |Real.cos x| = 1 := by
    have h1 : Real.cos x ^ 2 = 1 := by rw [hx] at h; linarith
    nlinarith [abs_nonneg (Real.cos x), sq_abs (Real.cos x)]
  rw [Real.sin_add, hx, zero_mul, zero_add, abs_mul, hc, one_mul]

lemma sin_pi_m {m : ℕ} (hm : 2 ≤ m) : 0 < Real.sin (Real.pi / (m:ℝ)) := by
  rw [pi_div_m hm, sin2]
  have := s_pos hm; have := c_pos hm; positivity

lemma xi0 {m : ℕ} (hm : 2 ≤ m) (k j : ℕ) (h : ip (n m j) (corner m k) = 0) :
    |ip (n m j) (corner m (k + 1))| = 2 * s m ^ 2 := by
  rw [ip_corner] at h ⊢
  have ht := rho0_pos hm
  have hsin : Real.sin (((j:ℝ) - k) * Real.pi / m) = 0 := by
    rcases mul_eq_zero.mp h with h | h
    · exact absurd h ht.ne'
    · exact h
  have e : (((j:ℝ) - ((k + 1 : ℕ) : ℝ)) * Real.pi / m)
      = ((j:ℝ) - k) * Real.pi / m - Real.pi / m := by push_cast; ring
  rw [e, abs_mul, abs_of_pos ht, abs_sin_sub _ _ hsin, abs_of_pos (sin_pi_m hm), pi_div_m hm,
    sin2]
  have := rho0_c hm
  linear_combination (2 * s m) * this

lemma xi1 {m : ℕ} (hm : 2 ≤ m) (k j : ℕ) (h : ip (n m j) (corner m (k + 1)) = 0) :
    |ip (n m j) (corner m k)| = 2 * s m ^ 2 := by
  rw [ip_corner] at h ⊢
  have ht := rho0_pos hm
  have hsin : Real.sin (((j:ℝ) - ((k + 1 : ℕ) : ℝ)) * Real.pi / m) = 0 := by
    rcases mul_eq_zero.mp h with h | h
    · exact absurd h ht.ne'
    · exact h
  have e : (((j:ℝ) - k) * Real.pi / m)
      = (((j:ℝ) - ((k + 1 : ℕ) : ℝ)) * Real.pi / m) + Real.pi / m := by push_cast; ring
  rw [e, abs_mul, abs_of_pos ht, abs_sin_add _ _ hsin, abs_of_pos (sin_pi_m hm), pi_div_m hm,
    sin2]
  have := rho0_c hm
  linear_combination (2 * s m) * this

lemma geomL {m : ℕ} (hm : 2 ≤ m) (k : ℕ) :
    edgeL (corner m k) (corner m (k + 1)) = 2 * rho0 m * s m := by
  have ht := rho0_pos hm
  have hs := s_pos hm
  have hA := Real.sin_sq_add_cos_sq ((k:ℝ) * Real.pi / m)
  have hS := sc_sq m
  rw [corner_succ hm k]
  unfold edgeL corner
  simp only [Real.sin_add, Real.cos_add, sin2, cos2]
  rw [← Real.sqrt_sq (show 0 ≤ 2 * rho0 m * s m by positivity)]
  congr 1
  linear_combination (rho0 m ^ 2 * ((2 * s m ^ 2) ^ 2 + (2 * s m * c m) ^ 2)) * hA
    + (4 * rho0 m ^ 2 * s m ^ 2) * hS

lemma geomDet {m : ℕ} (hm : 2 ≤ m) (k : ℕ) :
    (corner m k).1 * (corner m (k + 1)).2 - (corner m k).2 * (corner m (k + 1)).1
      = 2 * rho0 m ^ 2 * s m * c m := by
  have hA := Real.sin_sq_add_cos_sq ((k:ℝ) * Real.pi / m)
  rw [corner_succ hm k]
  unfold corner
  simp only [Real.sin_add, Real.cos_add, sin2, cos2]
  linear_combination (2 * rho0 m ^ 2 * s m * c m) * hA

lemma geomL0 {m : ℕ} (hm : 2 ≤ m) (k : ℕ) :
    ip (-(corner m k)) (corner m (k + 1) - corner m k) = 2 * rho0 m ^ 2 * s m ^ 2 := by
  have hA := Real.sin_sq_add_cos_sq ((k:ℝ) * Real.pi / m)
  rw [corner_succ hm k]
  unfold corner
  simp only [ip, Prod.fst_neg, Prod.snd_neg, Prod.fst_sub, Prod.snd_sub, Real.sin_add,
    Real.cos_add, sin2, cos2]
  linear_combination (2 * rho0 m ^ 2 * s m ^ 2) * hA

end T05

theorem T05.proof_1 : T05.stmt_1 := by
  unfold T05.stmt_1 T05.Cm T05.s T05.c
  have e : Real.pi / (2 * ((2:ℕ):ℝ)) = Real.pi / 4 := by push_cast; ring
  rw [e, Real.sin_pi_div_four, Real.cos_pi_div_four]
  have hs : (Real.sqrt 2 / 2) ^ 2 = 1 / 2 := by
    rw [div_pow, Real.sq_sqrt (by norm_num)]; norm_num
  have h4 : (Real.sqrt 2 / 2) ^ 4 = 1 / 4 := by
    rw [show (Real.sqrt 2 / 2) ^ 4 = ((Real.sqrt 2 / 2) ^ 2) ^ 2 by ring, hs]; norm_num
  rw [hs, h4]
  push_cast
  ring

theorem T05.proof_2 : T05.stmt_2 := by
  unfold T05.stmt_2 T05.Cm T05.s T05.c
  have e : Real.pi / (2 * ((3:ℕ):ℝ)) = Real.pi / 6 := by push_cast; ring
  rw [e, Real.sin_pi_div_six, Real.cos_pi_div_six]
  push_cast
  ring

theorem T05.proof_3 : T05.stmt_3 := by
  intro m hm
  have h2 : 2 ≤ m := by omega
  have hs := T05.s_pos h2
  have hc := T05.c_pos h2
  have hh := T05.s_le_half hm
  have hmpos : (0:ℝ) < m := by
    have : (2:ℝ) ≤ m := by exact_mod_cast h2
    linarith
  unfold T05.Cm
  have hneg : 3 * T05.s m ^ 2 - 1 < 0 := by nlinarith
  have hpos : 0 < 48 * (m:ℝ) * T05.s m ^ 4 * T05.c m := by positivity
  exact mul_neg_of_pos_of_neg hpos hneg

theorem T05.proof_4 : T05.stmt_4 := by
  intro m hm k j0 j1 _ _ h0 h1
  have hs := T05.s_pos hm
  have hc := T05.c_pos hm
  have ht := T05.rho0_pos hm
  have hmp := T05.mpos hm
  have htc := T05.rho0_c hm
  have hL := T05.geomL hm k
  have hD := T05.geomDet hm k
  have hL0 := T05.geomL0 hm k
  have hx0 := T05.xi0 hm k j0 h0
  have hx1 := T05.xi1 hm k j1 h1
  have hLpos : 0 < 2 * T05.rho0 m * T05.s m := by positivity
  have hH : T05.edgeH (T05.corner m k) (T05.corner m (k + 1)) = T05.s m := by
    unfold T05.edgeH
    rw [hL, hD, abs_of_pos (by positivity), div_eq_iff hLpos.ne']
    linear_combination (2 * T05.rho0 m * T05.s m) * htc
  have hL0' : T05.edgeL0 (T05.corner m k) (T05.corner m (k + 1)) = T05.rho0 m * T05.s m := by
    unfold T05.edgeL0
    rw [hL, hL0, div_eq_iff hLpos.ne']
    ring
  unfold T05.Ze
  simp only [hL, hH, hL0', hx0, hx1]
  unfold T05.Cm
  rw [T05.rho0_eq]
  field_simp
  ring
theorem T05.proof_5 : T05.stmt_5 := by
  intro m hm k
  have hs := T05.s_pos hm
  have hc := T05.c_pos hm
  have hmne : (m:ℝ) ≠ 0 := by
    have : (2:ℝ) ≤ m := by exact_mod_cast hm
    linarith
  set G : ℤ → ℝ := fun x => |Real.sin ((x:ℝ) * Real.pi / m)| with hGdef
  have hG : ∀ x : ℤ, G (x + m) = G x := by
    intro x
    simp only [hGdef]
    have : (((x + (m:ℤ) : ℤ) : ℝ)) * Real.pi / m = (x:ℝ) * Real.pi / m + Real.pi := by
      push_cast; field_simp
    rw [this, Real.sin_add_pi, abs_neg]
  have hterm : ∀ i : ℕ, |T05.ip (T05.n m i) (T05.corner m k)|
      = T05.rho0 m * G ((i:ℤ) - k) := by
    intro i
    rw [T05.ip_corner, abs_mul, abs_of_pos (by rw [T05.rho0_eq]; positivity)]
    simp only [hGdef]
    push_cast; rfl
  simp only [hterm]
  rw [← Finset.mul_sum, T05.sum_Icc_range (fun i => G ((i:ℤ) - k))]
  have := T05.sum_shift m G hG k
  push_cast at this ⊢
  rw [this]
  have h2 : ∑ i ∈ Finset.range m, G ((i:ℤ) + 1)
      = ∑ i ∈ Finset.range m, Real.sin (((i:ℝ) + 1) * Real.pi / m) := by
    apply Finset.sum_congr rfl
    intro i hi
    simp only [hGdef]
    push_cast
    apply abs_of_nonneg
    apply Real.sin_nonneg_of_nonneg_of_le_pi
    · positivity
    · have : (i:ℝ) + 1 ≤ m := by
        have : i < m := Finset.mem_range.mp hi
        exact_mod_cast this
      rw [div_le_iff₀ (by positivity)]
      nlinarith [Real.pi_pos]
  rw [h2, T05.sin_total hm, T05.rho0_eq]
  field_simp
