import Mathlib
/-!
# T08 — edge sums `Σ_e Z_e`, `Σ_e Z_e^{(5)}`, `Σ_e Z_e^{(2,3)}` over the regular `2m`-gon

Target id: T08.
Paper labels: `prop:edgelaw` (the coefficient `Z_e`), `thm:Cm` (the constant `C_m`),
`prop:order5` (`Z_e^{(5)}`, `Z_e^{(2,3)}`, `C_m^{(5)}`, `C_m^{(2,3)}`).
Paper lines: main.tex 1657–1663 (`thm:Cm`), 1671–1684 (`prop:edgelaw`), 1686 (the
sentence `\sum_eZ_e=C_m`), 1753 (definitions of `\bar\xi`, `\xi_+`, `\xi_\times`),
1755–1770 (`prop:order5`). The regular polygon and its corners are those of the proof
supplement, paper/ieee/proofs/sec_g4_poly.tex (setting) and sec_g4_Cm.tex (proof of
`thm:Cm`); the per-edge reduction for order five is in sec_g9_order5.tex lines 128–140.

Claims quoted verbatim.

main.tex 1657–1663:
  \begin{theorem}[Regular $2m$-gons$^\dagger$]\label{thm:Cm}
  For the smoothed regular $2m$-gon, the $\zeta(3)$ part of what the mixed terms
  contribute to $P'''(1)$ is exactly $C_m\,\zeta(3)$, where
  \[
  C_m=48m\sin^4\frac{\pi}{2m}\,\cos\frac{\pi}{2m}\Bigl(3\sin^2\frac{\pi}{2m}-1\Bigr).
  \]

main.tex 1671–1678 (Proposition `prop:edgelaw`, the part used here):
  Let $e$ be an edge of a centrally symmetric convex polygon, from the corner $V_0$ to
  the corner $V_1$, of length $L$ and at distance $h$ from the center. The foot of the
  perpendicular from the center cuts the line of $e$ into signed lengths $L_0$ and
  $L_1=L-L_0$, measured from $V_0$ and $V_1$. If $\langle v_{j_0},x\rangle=0$ is the
  kink line through $V_0$ and $\langle v_{j_1},x\rangle=0$ the one through $V_1$, put
  $\xi_0=|\langle v_{j_0},V_1\rangle|$ and $\xi_1=|\langle v_{j_1},V_0\rangle|$. [...]
  \[
  Z_e=6\,\xi_0\xi_1\frac{h^2}{L}\Bigl(3\,\frac{\xi_0L_0+\xi_1L_1}{L}-2\Bigr).
  \]

main.tex 1686:
  For the regular $2m$-gon $K_1$ of Theorem~\ref{thm:Cm}, $\sum_eZ_e=C_m$

main.tex 1753 (one word, the noun after "Keep the", elided as [...] so that this file
contains no command keyword):
  Keep the [...] of Proposition~\ref{prop:edgelaw} and write
  $\bar\xi=\frac{\xi_0L_0+\xi_1L_1}L$, $\xi_+=\xi_0+\xi_1$, $\xi_\times=\xi_0\xi_1$.

main.tex 1757–1770 (Proposition `prop:order5`, the formulas):
  \begin{align*}
  Z_e^{(5)}&=360\,\xi_\times\frac{h^2}L\Bigl[20\bar\xi^3-32\bar\xi^2+15\bar\xi-2\\
  &\quad-4(5\bar\xi-4)\,\xi_\times\frac{L_0L_1}{L^2}+\frac{h^2}{L^2}\Bigl(5\xi_\times\xi_++8\xi_+^2\\
  &\quad-28\xi_\times-15(\xi_+^2-3\xi_\times)\bar\xi\Bigr)\Bigr],\\
  Z_e^{(2,3)}&=120\,\xi_\times\frac{h^2}L\Bigl[12\bar\xi^2-12\bar\xi+2\\
  &\quad-12(5\bar\xi-4)\,\xi_\times\frac{L_0L_1}{L^2}+\frac{h^2}{L^2}\Bigl(15\xi_\times\xi_+\\
  &\quad-3\xi_+^2+24\xi_\times-45\xi_\times\bar\xi\Bigr)\Bigr].
  \end{align*}
  [...] For the regular $2m$-gon $K_1$ of Theorem~\ref{thm:Cm}, with
  $s=\sin\frac\pi{2m}$ and $c=\cos\frac\pi{2m}$,
  \begin{align*}
  \sum_eZ_e^{(5)}&=C_m^{(5)}=2880\,m\,s^4c\,(65s^6-63s^4+17s^2-1),\\
  \sum_eZ_e^{(2,3)}&=C_m^{(2,3)}=-960\,m\,c\,s^4(45s^6-27s^4+6s^2-1).
  \end{align*}

The regular `2m`-gon (sec_g4_poly.tex, Setting of the paper):
  $n_i=(\cos\frac{i\pi}m,\sin\frac{i\pi}m)$, $i=1,\dots,m$, and
  $K_p=\{x\in\RR^2:\sum_{i=1}^m|\ip{n_i}{x}|^p\le1\}$; so $K_1=\{\sum_i|\ip{n_i}{x}|\le1\}$
  and the kink lines are $\ip{n_i}{x}=0$ (main.tex: "its kink lines are the lines
  $\langle v_j,x\rangle=0$").
Corners (sec_g4_Cm.tex): $C_k=\varrho^kC_0$ with $C_0=\rho_0(0,1)$,
  $\rho_0=\tan\frac\pi{2m}$, i.e. `C_k = ρ₀ (-sin(kπ/m), cos(kπ/m))`; the edge `e_k`
  runs from `C_k` to `C_{k+1}`, `k = 0, …, 2m-1` (and `C_{2m} = C_0`).

Formalization choices:
* The edge data are NOT inserted by hand into the sums. For the edge `C_k C_{k+1}`,
  `L = |C_{k+1}-C_k|`, `h = |det(C_k,C_{k+1})|/L` (distance from the origin to the line
  of the edge), `L_0 = ⟨-C_k, C_{k+1}-C_k⟩/L` (signed distance from `C_k` to the foot of
  the perpendicular), `L_1 = L - L_0`, and `ξ_0 = |⟨n_{j_0},C_{k+1}⟩|`,
  `ξ_1 = |⟨n_{j_1},C_k⟩|`, where `n_{j_0}`, `n_{j_1}` are kink normals with
  `⟨n_{j_0},C_k⟩ = 0`, `⟨n_{j_1},C_{k+1}⟩ = 0`, `j_0, j_1 ∈ {1,…,m}`. In the sums the
  kink index of each corner is an arbitrary function `j` satisfying these hypotheses
  (`stmt_2` says such an index exists and is unique, so the hypotheses are not vacuous
  and determine the data).
* `stmt_1` and `stmt_2` tie the corners to the polygon: every segment `C_k C_{k+1}` lies
  on `∂K_1 = {Σ_i |⟨n_i,x⟩| = 1}`, and each `C_k` lies on exactly one kink line.
* `stmt_3` records the edge data the paper derives (sec_g4_Cm.tex):
  `L = 2s²/c`, `h = s`, `L_0 = L_1 = L/2`, `ξ_0 = ξ_1 = 2s²`.
* `stmt_4`–`stmt_6` are the three sums over the `2m` edges (the paper sums over all
  `2m` edges: "$\sum_eZ_e=2m\,Z_e$", sec_g4_Cm.tex), each against the paper's closed
  form in `s = sin(π/(2m))`, `c = cos(π/(2m))`.
* `stmt_7`–`stmt_9` are the same three identities as polynomial identities in real
  parameters `s, c` with `s² + c² = 1`, `s > 0`, `c > 0`, after inserting the edge data
  of `stmt_3` into the paper's per-edge formulas (`s, c > 0` is the paper's range; it
  only rules out the Lean junk value of division by `L = 2s²/c = 0`).
* `m ≥ 2` throughout (the paper's range: `2m`-gon with `m ≥ 2`).
* `ζ(3)`, `ζ(5)`, `ζ(2)ζ(3)` themselves do not appear: the claims formalized are the
  identities between the edge coefficients and the closed forms.

NOT COVERED:
* The analytic content of `thm:Cm`, `prop:edgelaw` and `prop:order5`: that `Z_e`
  (resp. `Z_e^{(5)}`, `Z_e^{(2,3)}`) is the `ζ(3)` (resp. `ζ(5)`, `ζ(2)ζ(3)`) part of
  what the local terms of `e` contribute to `P'''(1)` (resp. `P^{(5)}(1)`), and that
  `C_m ζ(3)` is the `ζ(3)` part of the mixed terms. This needs the perimeter `P(p)` of
  `K_p`, its derivatives at `p = 1`, the decomposition of the length element into
  local/mixed/far terms and a notion of "the `ζ(k)` part" of a number (formal
  coefficients after Beta-function evaluation); none of this can be stated faithfully
  in Mathlib today.
* That the `C_k` are exactly the vertices of `K_1` (only that consecutive `C_k` span
  boundary segments, and that each lies on a kink line, is stated).
* The parallelogram and affine-regular hexagon laws and the special values
  (`C_2`, `C_3`, sign of `C_m`: those are target T05).
-/
-- BEGIN STATEMENTS
namespace T08

/-- `s = sin(π/(2m))`. -/
noncomputable def s (m : ℕ) : ℝ := Real.sin (Real.pi / (2 * (m : ℝ)))

/-- `c = cos(π/(2m))`. -/
noncomputable def c (m : ℕ) : ℝ := Real.cos (Real.pi / (2 * (m : ℝ)))

/-- Euclidean inner product on `ℝ × ℝ`. -/
def ip (a b : ℝ × ℝ) : ℝ := a.1 * b.1 + a.2 * b.2

/-- Kink normal `n_i = (cos(iπ/m), sin(iπ/m))`. -/
noncomputable def n (m i : ℕ) : ℝ × ℝ :=
  (Real.cos ((i : ℝ) * Real.pi / m), Real.sin ((i : ℝ) * Real.pi / m))

/-- `ρ₀ = tan(π/(2m))`. -/
noncomputable def rho0 (m : ℕ) : ℝ := Real.tan (Real.pi / (2 * (m : ℝ)))

/-- Corner `C_k = ρ₀ (-sin(kπ/m), cos(kπ/m))` of the regular `2m`-gon `K_1`. -/
noncomputable def corner (m k : ℕ) : ℝ × ℝ :=
  (rho0 m * (-Real.sin ((k : ℝ) * Real.pi / m)), rho0 m * Real.cos ((k : ℝ) * Real.pi / m))

/-- The point `V₀ + t (V₁ - V₀)` of the segment `V₀V₁`. -/
def seg (V0 V1 : ℝ × ℝ) (t : ℝ) : ℝ × ℝ :=
  (V0.1 + t * (V1.1 - V0.1), V0.2 + t * (V1.2 - V0.2))

/-- Edge length `L = |V₁ - V₀|`. -/
noncomputable def edgeL (V0 V1 : ℝ × ℝ) : ℝ :=
  Real.sqrt ((V1.1 - V0.1) ^ 2 + (V1.2 - V0.2) ^ 2)

/-- Distance `h` from the center to the line through `V₀, V₁`. -/
noncomputable def edgeH (V0 V1 : ℝ × ℝ) : ℝ :=
  |V0.1 * V1.2 - V0.2 * V1.1| / edgeL V0 V1

/-- Signed length `L₀` from `V₀` to the foot of the perpendicular from the center. -/
noncomputable def edgeL0 (V0 V1 : ℝ × ℝ) : ℝ :=
  ip (-V0) (V1 - V0) / edgeL V0 V1

/-- `Z_e` of `prop:edgelaw` in the edge data `ξ₀, ξ₁, h, L, L₀` (`L₁ = L - L₀`):
`6 ξ₀ ξ₁ (h²/L) (3 (ξ₀L₀ + ξ₁L₁)/L - 2)`. -/
noncomputable def ZData (xi0 xi1 h L L0 : ℝ) : ℝ :=
  let L1 := L - L0
  6 * xi0 * xi1 * (h ^ 2 / L) * (3 * ((xi0 * L0 + xi1 * L1) / L) - 2)

/-- `Z_e^{(5)}` of `prop:order5` in the edge data, with
`ξ̄ = (ξ₀L₀ + ξ₁L₁)/L`, `ξ₊ = ξ₀ + ξ₁`, `ξₓ = ξ₀ξ₁`. -/
noncomputable def Z5Data (xi0 xi1 h L L0 : ℝ) : ℝ :=
  let L1 := L - L0
  let xb := (xi0 * L0 + xi1 * L1) / L
  let xp := xi0 + xi1
  let xx := xi0 * xi1
  360 * xx * (h ^ 2 / L) *
    (20 * xb ^ 3 - 32 * xb ^ 2 + 15 * xb - 2
      - 4 * (5 * xb - 4) * xx * (L0 * L1 / L ^ 2)
      + (h ^ 2 / L ^ 2) * (5 * xx * xp + 8 * xp ^ 2 - 28 * xx - 15 * (xp ^ 2 - 3 * xx) * xb))

/-- `Z_e^{(2,3)}` of `prop:order5` in the edge data. -/
noncomputable def Z23Data (xi0 xi1 h L L0 : ℝ) : ℝ :=
  let L1 := L - L0
  let xb := (xi0 * L0 + xi1 * L1) / L
  let xp := xi0 + xi1
  let xx := xi0 * xi1
  120 * xx * (h ^ 2 / L) *
    (12 * xb ^ 2 - 12 * xb + 2
      - 12 * (5 * xb - 4) * xx * (L0 * L1 / L ^ 2)
      + (h ^ 2 / L ^ 2) * (15 * xx * xp - 3 * xp ^ 2 + 24 * xx - 45 * xx * xb))

/-- `Z_e` for the edge `V₀V₁` with kink normal `v₀` through `V₀` and `v₁` through
`V₁`, all edge data computed from the geometry
(`ξ₀ = |⟨v₀,V₁⟩|`, `ξ₁ = |⟨v₁,V₀⟩|`). -/
noncomputable def Ze (V0 V1 v0 v1 : ℝ × ℝ) : ℝ :=
  ZData |ip v0 V1| |ip v1 V0| (edgeH V0 V1) (edgeL V0 V1) (edgeL0 V0 V1)

/-- `Z_e^{(5)}` for the edge `V₀V₁`, from the geometry. -/
noncomputable def Ze5 (V0 V1 v0 v1 : ℝ × ℝ) : ℝ :=
  Z5Data |ip v0 V1| |ip v1 V0| (edgeH V0 V1) (edgeL V0 V1) (edgeL0 V0 V1)

/-- `Z_e^{(2,3)}` for the edge `V₀V₁`, from the geometry. -/
noncomputable def Ze23 (V0 V1 v0 v1 : ℝ × ℝ) : ℝ :=
  Z23Data |ip v0 V1| |ip v1 V0| (edgeH V0 V1) (edgeL V0 V1) (edgeL0 V0 V1)

/-- `48 m s⁴ c (3s² - 1)` as a polynomial in `m, s, c`. -/
def CmPoly (m s c : ℝ) : ℝ := 48 * m * s ^ 4 * c * (3 * s ^ 2 - 1)

/-- `2880 m s⁴ c (65s⁶ - 63s⁴ + 17s² - 1)`. -/
def Cm5Poly (m s c : ℝ) : ℝ := 2880 * m * s ^ 4 * c * (65 * s ^ 6 - 63 * s ^ 4 + 17 * s ^ 2 - 1)

/-- `-960 m c s⁴ (45s⁶ - 27s⁴ + 6s² - 1)`. -/
def Cm23Poly (m s c : ℝ) : ℝ := -960 * m * c * s ^ 4 * (45 * s ^ 6 - 27 * s ^ 4 + 6 * s ^ 2 - 1)

/-- `C_m` of `thm:Cm`. -/
noncomputable def Cm (m : ℕ) : ℝ := CmPoly m (s m) (c m)

/-- `C_m^{(5)}` of `prop:order5`. -/
noncomputable def Cm5 (m : ℕ) : ℝ := Cm5Poly m (s m) (c m)

/-- `C_m^{(2,3)}` of `prop:order5`. -/
noncomputable def Cm23 (m : ℕ) : ℝ := Cm23Poly m (s m) (c m)

/-- Setup ("For the regular $2m$-gon $K_1$"): every segment `C_k C_{k+1}` lies on
`∂K_1 = {Σ_{i=1}^m |⟨n_i,x⟩| = 1}`, so the `C_k C_{k+1}` are edges of `K_1`. -/
def stmt_1 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∀ k : ℕ, ∀ t : ℝ, 0 ≤ t → t ≤ 1 →
    ∑ i ∈ Finset.Icc 1 m, |ip (n m i) (seg (corner m k) (corner m (k + 1)) t)| = 1

/-- Setup ("the kink line through $V_0$"): each corner `C_k` lies on exactly one kink
line `⟨n_j,x⟩ = 0`, `j ∈ {1,…,m}`. -/
def stmt_2 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∀ k : ℕ, ∃! j : ℕ, j ∈ Finset.Icc 1 m ∧ ip (n m j) (corner m k) = 0

/-- Edge data of the regular `2m`-gon (sec_g4_Cm.tex, proof of `thm:Cm`):
`L = 2s²/c`, `h = s`, `L₀ = L/2`, `ξ₀ = ξ₁ = 2s²`, derived from the corners and the
kink lines through them. -/
def stmt_3 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∀ k j0 j1 : ℕ,
    j0 ∈ Finset.Icc 1 m → j1 ∈ Finset.Icc 1 m →
    ip (n m j0) (corner m k) = 0 → ip (n m j1) (corner m (k + 1)) = 0 →
    edgeL (corner m k) (corner m (k + 1)) = 2 * s m ^ 2 / c m ∧
    edgeH (corner m k) (corner m (k + 1)) = s m ∧
    edgeL0 (corner m k) (corner m (k + 1)) = edgeL (corner m k) (corner m (k + 1)) / 2 ∧
    |ip (n m j0) (corner m (k + 1))| = 2 * s m ^ 2 ∧
    |ip (n m j1) (corner m k)| = 2 * s m ^ 2

/-- main.tex 1686: "For the regular $2m$-gon $K_1$ of Theorem~\ref{thm:Cm},
$\sum_eZ_e=C_m$", the sum over the `2m` edges `C_k C_{k+1}`, `k = 0,…,2m-1`, with
`j k` the kink line through `C_k`. -/
def stmt_4 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∀ j : ℕ → ℕ,
    (∀ k, j k ∈ Finset.Icc 1 m ∧ ip (n m (j k)) (corner m k) = 0) →
    ∑ k ∈ Finset.range (2 * m),
      Ze (corner m k) (corner m (k + 1)) (n m (j k)) (n m (j (k + 1))) = Cm m

/-- `prop:order5`: "$\sum_eZ_e^{(5)}=C_m^{(5)}=2880\,m\,s^4c\,(65s^6-63s^4+17s^2-1)$". -/
def stmt_5 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∀ j : ℕ → ℕ,
    (∀ k, j k ∈ Finset.Icc 1 m ∧ ip (n m (j k)) (corner m k) = 0) →
    ∑ k ∈ Finset.range (2 * m),
      Ze5 (corner m k) (corner m (k + 1)) (n m (j k)) (n m (j (k + 1))) = Cm5 m

/-- `prop:order5`: "$\sum_eZ_e^{(2,3)}=C_m^{(2,3)}=-960\,m\,c\,s^4(45s^6-27s^4+6s^2-1)$". -/
def stmt_6 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∀ j : ℕ → ℕ,
    (∀ k, j k ∈ Finset.Icc 1 m ∧ ip (n m (j k)) (corner m k) = 0) →
    ∑ k ∈ Finset.range (2 * m),
      Ze23 (corner m k) (corner m (k + 1)) (n m (j k)) (n m (j (k + 1))) = Cm23 m

/-- Algebraic core of `\sum_eZ_e=C_m`: with the edge data `ξ₀ = ξ₁ = 2s²`, `h = s`,
`L = 2s²/c`, `L₀ = L/2 = s²/c` and `s² + c² = 1`, `2m · Z_e = C_m`. -/
def stmt_7 : Prop :=
  ∀ m s c : ℝ, 0 < s → 0 < c → s ^ 2 + c ^ 2 = 1 →
    2 * m * ZData (2 * s ^ 2) (2 * s ^ 2) s (2 * s ^ 2 / c) (s ^ 2 / c) = CmPoly m s c

/-- Algebraic core of `\sum_eZ_e^{(5)}=C_m^{(5)}`. -/
def stmt_8 : Prop :=
  ∀ m s c : ℝ, 0 < s → 0 < c → s ^ 2 + c ^ 2 = 1 →
    2 * m * Z5Data (2 * s ^ 2) (2 * s ^ 2) s (2 * s ^ 2 / c) (s ^ 2 / c) = Cm5Poly m s c

/-- Algebraic core of `\sum_eZ_e^{(2,3)}=C_m^{(2,3)}`. -/
def stmt_9 : Prop :=
  ∀ m s c : ℝ, 0 < s → 0 < c → s ^ 2 + c ^ 2 = 1 →
    2 * m * Z23Data (2 * s ^ 2) (2 * s ^ 2) s (2 * s ^ 2 / c) (s ^ 2 / c) = Cm23Poly m s c

end T08
-- END STATEMENTS
namespace T08

lemma alpha_bounds (m : ℕ) (hm : 2 ≤ m) :
    0 < Real.pi / (2 * (m : ℝ)) ∧ Real.pi / (2 * (m : ℝ)) ≤ Real.pi / 4 := by
  have hm' : (2 : ℝ) ≤ m := by exact_mod_cast hm
  have hpi := Real.pi_pos
  constructor
  · exact div_pos hpi (by linarith)
  · rw [div_le_div_iff₀ (by linarith) (by norm_num)]
    nlinarith

lemma s_pos (m : ℕ) (hm : 2 ≤ m) : 0 < s m := by
  obtain ⟨h1, h2⟩ := alpha_bounds m hm
  exact Real.sin_pos_of_pos_of_lt_pi h1 (by linarith [Real.pi_pos])

lemma c_pos (m : ℕ) (hm : 2 ≤ m) : 0 < c m := by
  obtain ⟨h1, h2⟩ := alpha_bounds m hm
  exact Real.cos_pos_of_mem_Ioo ⟨by linarith [Real.pi_pos], by linarith [Real.pi_pos]⟩

lemma sc_sq (m : ℕ) : s m ^ 2 + c m ^ 2 = 1 := Real.sin_sq_add_cos_sq _

lemma rho_eq (m : ℕ) : rho0 m = s m / c m := Real.tan_eq_sin_div_cos _

lemma theta_eq (m : ℕ) (hm : 2 ≤ m) : Real.pi / (m : ℝ) = 2 * (Real.pi / (2 * (m : ℝ))) := by
  have : (m : ℝ) ≠ 0 := by
    have : (2 : ℝ) ≤ m := by exact_mod_cast hm
    linarith
  field_simp

lemma cos_theta (m : ℕ) (hm : 2 ≤ m) : Real.cos (Real.pi / (m : ℝ)) = 1 - 2 * s m ^ 2 := by
  rw [theta_eq m hm, Real.cos_two_mul]
  have := sc_sq m
  simp only [s, c] at this ⊢
  nlinarith

lemma sin_theta (m : ℕ) (hm : 2 ≤ m) : Real.sin (Real.pi / (m : ℝ)) = 2 * s m * c m := by
  rw [theta_eq m hm, Real.sin_two_mul]
  rfl

lemma ip_corner (m i k : ℕ) :
    ip (n m i) (corner m k) = rho0 m * Real.sin (((i : ℝ) - k) * Real.pi / m) := by
  unfold ip n corner
  simp only
  rw [show ((i : ℝ) - k) * Real.pi / m = (i : ℝ) * Real.pi / m - (k : ℝ) * Real.pi / m by ring,
    Real.sin_sub]
  ring

lemma corner_succ (m k : ℕ) :
    corner m (k + 1) = (rho0 m * (-Real.sin ((k : ℝ) * Real.pi / m + Real.pi / m)),
      rho0 m * Real.cos ((k : ℝ) * Real.pi / m + Real.pi / m)) := by
  have : (((k + 1 : ℕ) : ℝ) * Real.pi / m) = (k : ℝ) * Real.pi / m + Real.pi / m := by
    push_cast; ring
  unfold corner
  rw [this]

lemma sq_len (m k : ℕ) (hm : 2 ≤ m) :
    ((corner m (k + 1)).1 - (corner m k).1) ^ 2 + ((corner m (k + 1)).2 - (corner m k).2) ^ 2
      = 4 * s m ^ 2 * rho0 m ^ 2 := by
  rw [corner_succ]
  unfold corner
  simp only
  have e1 := Real.sin_sq_add_cos_sq ((k : ℝ) * Real.pi / m)
  have e2 := Real.sin_sq_add_cos_sq ((k : ℝ) * Real.pi / m + Real.pi / m)
  have hc := Real.cos_sub ((k : ℝ) * Real.pi / m + Real.pi / m) ((k : ℝ) * Real.pi / m)
  rw [show (k : ℝ) * Real.pi / m + Real.pi / m - (k : ℝ) * Real.pi / m = Real.pi / m by ring,
    cos_theta m hm] at hc
  linear_combination (rho0 m ^ 2) * (e1 + e2 + 2 * hc)

lemma det_edge (m k : ℕ) (hm : 2 ≤ m) :
    (corner m k).1 * (corner m (k + 1)).2 - (corner m k).2 * (corner m (k + 1)).1
      = rho0 m ^ 2 * (2 * s m * c m) := by
  rw [corner_succ, ← sin_theta m hm]
  unfold corner
  simp only
  have hc := Real.sin_sub ((k : ℝ) * Real.pi / m + Real.pi / m) ((k : ℝ) * Real.pi / m)
  rw [show (k : ℝ) * Real.pi / m + Real.pi / m - (k : ℝ) * Real.pi / m = Real.pi / m by ring] at hc
  rw [hc]
  ring

lemma dot_edge (m k : ℕ) (hm : 2 ≤ m) :
    ip (-(corner m k)) (corner m (k + 1) - corner m k) = rho0 m ^ 2 * (2 * s m ^ 2) := by
  have e1 := Real.sin_sq_add_cos_sq ((k : ℝ) * Real.pi / m)
  have hc := Real.cos_sub ((k : ℝ) * Real.pi / m + Real.pi / m) ((k : ℝ) * Real.pi / m)
  rw [show (k : ℝ) * Real.pi / m + Real.pi / m - (k : ℝ) * Real.pi / m = Real.pi / m by ring,
    cos_theta m hm] at hc
  rw [corner_succ]
  unfold ip corner
  simp only [Prod.fst_neg, Prod.snd_neg, Prod.fst_sub, Prod.snd_sub]
  linear_combination (rho0 m ^ 2) * (e1 + hc)

lemma edgeL_eq (m k : ℕ) (hm : 2 ≤ m) :
    edgeL (corner m k) (corner m (k + 1)) = 2 * s m ^ 2 / c m := by
  have hs := s_pos m hm
  have hc := c_pos m hm
  unfold edgeL
  rw [sq_len m k hm, rho_eq, show 4 * s m ^ 2 * (s m / c m) ^ 2 = (2 * s m ^ 2 / c m) ^ 2 by
    field_simp; ring]
  exact Real.sqrt_sq (by positivity)

lemma edgeH_eq (m k : ℕ) (hm : 2 ≤ m) :
    edgeH (corner m k) (corner m (k + 1)) = s m := by
  have hs := s_pos m hm
  have hc := c_pos m hm
  unfold edgeH
  rw [edgeL_eq m k hm, det_edge m k hm, rho_eq, abs_of_pos (by positivity)]
  field_simp

lemma edgeL0_eq (m k : ℕ) (hm : 2 ≤ m) :
    edgeL0 (corner m k) (corner m (k + 1)) = s m ^ 2 / c m := by
  have hs := s_pos m hm
  have hc := c_pos m hm
  unfold edgeL0
  rw [dot_edge m k hm, edgeL_eq m k hm, rho_eq]
  field_simp

lemma xi_aux (m : ℕ) (hm : 2 ≤ m) (x : ℝ) (hx : Real.sin x = 0) :
    |rho0 m * Real.sin (x - Real.pi / m)| = 2 * s m ^ 2 ∧
    |rho0 m * Real.sin (x + Real.pi / m)| = 2 * s m ^ 2 := by
  have hs := s_pos m hm
  have hc := c_pos m hm
  have hcx : |Real.cos x| = 1 := by
    have h2 : Real.cos x ^ 2 = 1 := by nlinarith [Real.sin_sq_add_cos_sq x]
    rw [← Real.sqrt_sq_eq_abs, h2, Real.sqrt_one]
  have hpos : 0 < rho0 m * (2 * s m * c m) := by rw [rho_eq]; positivity
  have hval : rho0 m * (2 * s m * c m) = 2 * s m ^ 2 := by
    rw [rho_eq]; field_simp
  constructor
  · rw [Real.sin_sub, hx, sin_theta m hm,
      show rho0 m * (0 * Real.cos (Real.pi / m) - Real.cos x * (2 * s m * c m))
        = -(rho0 m * (2 * s m * c m)) * Real.cos x by ring,
      abs_mul, abs_neg, abs_of_pos hpos, hcx, hval, mul_one]
  · rw [Real.sin_add, hx, sin_theta m hm,
      show rho0 m * (0 * Real.cos (Real.pi / m) + Real.cos x * (2 * s m * c m))
        = (rho0 m * (2 * s m * c m)) * Real.cos x by ring,
      abs_mul, abs_of_pos hpos, hcx, hval, mul_one]

lemma xi0_eq (m k j0 : ℕ) (hm : 2 ≤ m) (h0 : ip (n m j0) (corner m k) = 0) :
    |ip (n m j0) (corner m (k + 1))| = 2 * s m ^ 2 := by
  have hr : rho0 m ≠ 0 := by
    have := s_pos m hm; have := c_pos m hm
    rw [rho_eq]; positivity
  rw [ip_corner] at h0
  have hx : Real.sin (((j0 : ℝ) - k) * Real.pi / m) = 0 := by
    rcases mul_eq_zero.mp h0 with h | h
    · exact absurd h hr
    · exact h
  rw [ip_corner, show (((j0 : ℝ) - ((k + 1 : ℕ) : ℝ)) * Real.pi / m)
      = ((j0 : ℝ) - k) * Real.pi / m - Real.pi / m by push_cast; ring]
  exact (xi_aux m hm _ hx).1

lemma xi1_eq (m k j1 : ℕ) (hm : 2 ≤ m) (h1 : ip (n m j1) (corner m (k + 1)) = 0) :
    |ip (n m j1) (corner m k)| = 2 * s m ^ 2 := by
  have hr : rho0 m ≠ 0 := by
    have := s_pos m hm; have := c_pos m hm
    rw [rho_eq]; positivity
  rw [ip_corner] at h1
  have hx : Real.sin (((j1 : ℝ) - ((k + 1 : ℕ) : ℝ)) * Real.pi / m) = 0 := by
    rcases mul_eq_zero.mp h1 with h | h
    · exact absurd h hr
    · exact h
  rw [ip_corner, show (((j1 : ℝ) - (k : ℝ)) * Real.pi / m)
      = ((j1 : ℝ) - ((k + 1 : ℕ) : ℝ)) * Real.pi / m + Real.pi / m by push_cast; ring]
  exact (xi_aux m hm _ hx).2

end T08

theorem T08.geo_3 : T08.stmt_3 := by
  intro m hm k j0 j1 _ _ h0 h1
  refine ⟨T08.edgeL_eq m k hm, T08.edgeH_eq m k hm, ?_, T08.xi0_eq m k j0 hm h0,
    T08.xi1_eq m k j1 hm h1⟩
  rw [T08.edgeL0_eq m k hm, T08.edgeL_eq m k hm]
  ring

namespace T08

lemma telesc (a : ℝ) (n : ℕ) :
    2 * Real.sin a * ∑ i ∈ Finset.range n, Real.sin (2 * i * a)
      = Real.cos a - Real.cos ((2 * n - 1) * a) := by
  induction n with
  | zero =>
    have : ((2 * ((0 : ℕ) : ℝ) - 1) * a) = -a := by push_cast; ring
    rw [this, Real.cos_neg]; simp
  | succ n ih =>
    rw [Finset.sum_range_succ, mul_add, ih]
    have h1 : (2 * ((n + 1 : ℕ) : ℝ) - 1) * a = 2 * n * a + a := by push_cast; ring
    have h2 : (2 * (n : ℝ) - 1) * a = 2 * n * a - a := by ring
    rw [h1, h2, Real.cos_add, Real.cos_sub]; ring

lemma sum_sin (m : ℕ) (hm : 2 ≤ m) :
    ∑ i ∈ Finset.range m, |Real.sin ((i : ℝ) * Real.pi / m)| = c m / s m := by
  have hs := s_pos m hm
  have hmR : (0 : ℝ) < m := by
    have : (2 : ℝ) ≤ m := by exact_mod_cast hm
    linarith
  have hterm : ∀ i ∈ Finset.range m,
      |Real.sin ((i : ℝ) * Real.pi / m)| = Real.sin (2 * i * (Real.pi / (2 * m))) := by
    intro i hi
    have hi' : (i : ℝ) + 1 ≤ m := by
      have := Finset.mem_range.mp hi
      exact_mod_cast this
    have e : (i : ℝ) * Real.pi / m = 2 * i * (Real.pi / (2 * m)) := by field_simp
    rw [e]
    apply abs_of_nonneg
    apply Real.sin_nonneg_of_nonneg_of_le_pi
    · have := Real.pi_pos
      have : (0 : ℝ) ≤ i := Nat.cast_nonneg i
      positivity
    · rw [← e]
      rw [div_le_iff₀ hmR]
      nlinarith [Real.pi_pos]
  rw [Finset.sum_congr rfl hterm]
  have h := telesc (Real.pi / (2 * m)) m
  have e2 : (2 * (m : ℝ) - 1) * (Real.pi / (2 * m)) = Real.pi - Real.pi / (2 * m) := by
    field_simp
  rw [e2, Real.cos_pi_sub] at h
  rw [eq_div_iff hs.ne']
  have : s m = Real.sin (Real.pi / (2 * m)) := rfl
  have h' : c m = Real.cos (Real.pi / (2 * m)) := rfl
  rw [this, h']
  linarith

lemma sin_sign (m : ℕ) (hm : 0 < m) (j : ℤ) :
    0 ≤ Real.sin ((j : ℝ) * Real.pi / m) * Real.sin (((j : ℝ) + 1) * Real.pi / m) := by
  have hmR : (0 : ℝ) < m := by exact_mod_cast hm
  have hmZ : (0 : ℤ) < m := by exact_mod_cast hm
  have hjr : (m : ℤ) * (j / (m : ℤ)) + j % (m : ℤ) = j := Int.mul_ediv_add_emod j m
  have hr0 : 0 ≤ j % (m : ℤ) := Int.emod_nonneg _ hmZ.ne'
  have hr1 : j % (m : ℤ) < m := Int.emod_lt_of_pos _ hmZ
  generalize j / (m : ℤ) = q at hjr
  generalize j % (m : ℤ) = r at hjr hr0 hr1
  have hjR : (j : ℝ) = m * q + r := by rw [← hjr]; push_cast; ring
  have hr0R : (0 : ℝ) ≤ r := by exact_mod_cast hr0
  have hr1R : (r : ℝ) + 1 ≤ m := by exact_mod_cast hr1
  have ha : (j : ℝ) * Real.pi / m = r * Real.pi / m + q * Real.pi := by
    rw [hjR]; field_simp; ring
  have hb : ((j : ℝ) + 1) * Real.pi / m = (r + 1) * Real.pi / m + q * Real.pi := by
    rw [hjR]; field_simp; ring
  rw [ha, hb, Real.sin_add_int_mul_pi, Real.sin_add_int_mul_pi]
  have h1 : 0 ≤ Real.sin ((r : ℝ) * Real.pi / m) := by
    apply Real.sin_nonneg_of_nonneg_of_le_pi
    · exact div_nonneg (mul_nonneg hr0R Real.pi_pos.le) hmR.le
    · rw [div_le_iff₀ hmR]; nlinarith [Real.pi_pos]
  have h2 : 0 ≤ Real.sin (((r : ℝ) + 1) * Real.pi / m) := by
    apply Real.sin_nonneg_of_nonneg_of_le_pi
    · exact div_nonneg (mul_nonneg (by linarith) Real.pi_pos.le) hmR.le
    · rw [div_le_iff₀ hmR]; nlinarith [Real.pi_pos]
  rcases Int.even_or_odd q with h | h
  · rw [h.neg_one_zpow]; nlinarith [mul_nonneg h1 h2]
  · rw [h.neg_one_zpow]; nlinarith [mul_nonneg h1 h2]

lemma same_sign (a b : ℝ) (h : 0 ≤ a * b) : (0 ≤ a ∧ 0 ≤ b) ∨ (a ≤ 0 ∧ b ≤ 0) := by
  rcases le_total 0 a with ha | ha <;> rcases le_total 0 b with hb | hb
  · exact Or.inl ⟨ha, hb⟩
  · have : a * b = 0 := le_antisymm (mul_nonpos_of_nonneg_of_nonpos ha hb) h
    rcases mul_eq_zero.mp this with h0 | h0
    · exact Or.inr ⟨h0.le, hb⟩
    · exact Or.inl ⟨ha, h0.ge⟩
  · have : a * b = 0 := le_antisymm (mul_nonpos_of_nonpos_of_nonneg ha hb) h
    rcases mul_eq_zero.mp this with h0 | h0
    · exact Or.inl ⟨h0.ge, hb⟩
    · exact Or.inr ⟨ha, h0.le⟩
  · exact Or.inr ⟨ha, hb⟩

lemma abs_comb (a b t : ℝ) (h : 0 ≤ a * b) (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    |(1 - t) * a + t * b| = (1 - t) * |a| + t * |b| := by
  have h1 : 0 ≤ 1 - t := by linarith
  rcases same_sign a b h with ⟨ha, hb⟩ | ⟨ha, hb⟩
  · rw [abs_of_nonneg ha, abs_of_nonneg hb, abs_of_nonneg (by positivity)]
  · rw [abs_of_nonpos ha, abs_of_nonpos hb, abs_of_nonpos]
    · ring
    · nlinarith [mul_nonneg h1 (neg_nonneg.mpr ha), mul_nonneg ht0 (neg_nonneg.mpr hb)]

lemma sum_shift (m : ℕ) (hm : 2 ≤ m) (k : ℕ) :
    ∑ i ∈ Finset.range m, |Real.sin ((((i + 1 : ℕ) : ℝ) - k) * Real.pi / m)| = c m / s m := by
  have hmR : (0 : ℝ) < m := by
    have : (2 : ℝ) ≤ m := by exact_mod_cast hm
    linarith
  induction k with
  | zero =>
    have e1 := Finset.sum_range_succ' (fun i : ℕ => |Real.sin ((i : ℝ) * Real.pi / m)|) m
    have e2 := Finset.sum_range_succ (fun i : ℕ => |Real.sin ((i : ℝ) * Real.pi / m)|) m
    rw [← sum_sin m hm]
    simp only [Nat.cast_zero, sub_zero, zero_mul, zero_div, Real.sin_zero, abs_zero, add_zero] at e1 ⊢
    have : (m : ℝ) * Real.pi / m = Real.pi := by field_simp
    simp only [this, Real.sin_pi, abs_zero, add_zero] at e2
    rw [← e1, e2]
  | succ k ih =>
    rw [← ih]
    set g : ℕ → ℝ := fun i => |Real.sin (((i : ℝ) - k) * Real.pi / m)| with hg
    have e1 := Finset.sum_range_succ' g m
    have e2 := Finset.sum_range_succ g m
    have g0 : g 0 = g m := by
      simp only [hg, Nat.cast_zero]
      have : ((m : ℝ) - k) * Real.pi / m = Real.pi - k * Real.pi / m := by field_simp
      rw [this, Real.sin_pi_sub, show ((0 : ℝ) - k) * Real.pi / m = -(k * Real.pi / m) by ring,
        Real.sin_neg, abs_neg]
    have l : ∑ i ∈ Finset.range m, |Real.sin ((((i + 1 : ℕ) : ℝ) - k) * Real.pi / m)|
        = ∑ i ∈ Finset.range m, g (i + 1) := rfl
    have r : ∑ i ∈ Finset.range m, |Real.sin ((((i + 1 : ℕ) : ℝ) - ((k + 1 : ℕ) : ℝ)) * Real.pi / m)|
        = ∑ i ∈ Finset.range m, g i := by
      apply Finset.sum_congr rfl
      intro i _
      simp only [hg]
      congr 3
      push_cast; ring
    rw [r, ← l]
    linarith

end T08

namespace T08

lemma term_eq (m : ℕ) (hm : 2 ≤ m) (k i : ℕ) (t : ℝ) (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    |ip (n m i) (seg (corner m k) (corner m (k + 1)) t)|
      = rho0 m * ((1 - t) * |Real.sin (((i : ℝ) - k) * Real.pi / m)|
          + t * |Real.sin (((i : ℝ) - ((k + 1 : ℕ) : ℝ)) * Real.pi / m)|) := by
  have hs := s_pos m hm
  have hc := c_pos m hm
  have hρ : 0 < rho0 m := by rw [rho_eq]; positivity
  have hip : ip (n m i) (seg (corner m k) (corner m (k + 1)) t)
      = (1 - t) * ip (n m i) (corner m k) + t * ip (n m i) (corner m (k + 1)) := by
    unfold ip seg; ring
  rw [hip, ip_corner, ip_corner,
    show (1 - t) * (rho0 m * Real.sin (((i : ℝ) - k) * Real.pi / m))
        + t * (rho0 m * Real.sin (((i : ℝ) - ((k + 1 : ℕ) : ℝ)) * Real.pi / m))
      = rho0 m * ((1 - t) * Real.sin (((i : ℝ) - k) * Real.pi / m)
        + t * Real.sin (((i : ℝ) - ((k + 1 : ℕ) : ℝ)) * Real.pi / m)) by ring,
    abs_mul, abs_of_pos hρ]
  congr 1
  apply abs_comb _ _ _ _ ht0 ht1
  have := sin_sign m (by omega) ((i : ℤ) - k - 1)
  have e1 : (((i : ℤ) - k - 1 : ℤ) : ℝ) * Real.pi / m
      = (((i : ℝ) - ((k + 1 : ℕ) : ℝ)) * Real.pi / m) := by push_cast; ring
  have e2 : ((((i : ℤ) - k - 1 : ℤ) : ℝ) + 1) * Real.pi / m = ((i : ℝ) - k) * Real.pi / m := by
    push_cast; ring
  rw [e1, e2] at this
  linarith

lemma kink_iff (m : ℕ) (hm : 2 ≤ m) (i k : ℕ) :
    ip (n m i) (corner m k) = 0 ↔ (m : ℤ) ∣ (i : ℤ) - k := by
  have hs := s_pos m hm
  have hc := c_pos m hm
  have hmR : (0 : ℝ) < m := by
    have : (2 : ℝ) ≤ m := by exact_mod_cast hm
    linarith
  have hr : rho0 m ≠ 0 := by rw [rho_eq]; positivity
  rw [ip_corner, mul_eq_zero, or_iff_right hr, Real.sin_eq_zero_iff]
  constructor
  · rintro ⟨q, hq⟩
    have h1 : (q : ℝ) * m = (i : ℝ) - k := by
      field_simp at hq
      have := Real.pi_pos
      nlinarith [hq]
    have h2 : q * (m : ℤ) = (i : ℤ) - k := by exact_mod_cast h1
    exact ⟨q, by linarith⟩
  · rintro ⟨q, hq⟩
    have h1 : (i : ℝ) - k = m * q := by exact_mod_cast hq
    refine ⟨q, ?_⟩
    rw [h1]; field_simp

end T08

theorem T08.geo_1 : T08.stmt_1 := by
  intro m hm k t ht0 ht1
  have hs := T08.s_pos m hm
  have hc := T08.c_pos m hm
  have hρ : T08.rho0 m = T08.s m / T08.c m := T08.rho_eq m
  rw [Finset.sum_congr rfl (fun i _ => T08.term_eq m hm k i t ht0 ht1)]
  have hIcc : Finset.Icc 1 m = Finset.Ico 1 (m + 1) := by
    ext x; simp only [Finset.mem_Icc, Finset.mem_Ico]; omega
  rw [hIcc, Finset.sum_Ico_eq_sum_range]
  simp only [Nat.succ_sub_one]
  rw [← Finset.mul_sum, Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum]
  have h1 := T08.sum_shift m hm k
  have h2 := T08.sum_shift m hm (k + 1)
  simp only [add_comm 1] at *
  rw [h1, h2, hρ]
  field_simp
  ring

theorem T08.geo_2 : T08.stmt_2 := by
  intro m hm k
  have hm0 : 0 < m := by omega
  have h := Nat.mod_add_div (k + (m - 1)) m
  have hlt := Nat.mod_lt (k + (m - 1)) hm0
  generalize (k + (m - 1)) % m = x at h hlt ⊢
  generalize (k + (m - 1)) / m = q at h
  have h' : (x : ℤ) + (m : ℤ) * (q : ℤ) = (k : ℤ) + m - 1 := by
    have h1 := congrArg (Nat.cast : ℕ → ℤ) h
    push_cast at h1
    have h2 : ((m - 1 : ℕ) : ℤ) = (m : ℤ) - 1 := by omega
    rw [h2] at h1
    linarith
  have hj : (m : ℤ) ∣ ((x + 1 : ℕ) : ℤ) - k := by
    refine ⟨1 - (q : ℤ), ?_⟩
    push_cast
    linear_combination h'
  refine ⟨x + 1, ⟨?_, ?_⟩, ?_⟩
  · simp only [Finset.mem_Icc]; omega
  · rw [T08.kink_iff m hm]; exact hj
  · rintro y ⟨hy, hyk⟩
    rw [T08.kink_iff m hm] at hyk
    simp only [Finset.mem_Icc] at hy
    have hd : (m : ℤ) ∣ (y : ℤ) - ((x + 1 : ℕ) : ℤ) := by
      have := dvd_sub hyk hj
      have e : (y : ℤ) - k - (((x + 1 : ℕ) : ℤ) - k) = (y : ℤ) - ((x + 1 : ℕ) : ℤ) := by ring
      rwa [e] at this
    have := Int.eq_zero_of_abs_lt_dvd hd (by rw [abs_lt]; push_cast; omega)
    omega

theorem T08.alg_7 : T08.stmt_7 := by
  intro m s c hs hc h
  have hc2 : c ^ 2 = 1 - s ^ 2 := by linarith
  simp only [T08.ZData, T08.CmPoly]
  field_simp
  ring

theorem T08.alg_8 : T08.stmt_8 := by
  intro m s c hs hc h
  have hc2 : c ^ 2 = 1 - s ^ 2 := by linarith
  simp only [T08.Z5Data, T08.Cm5Poly]
  field_simp
  linear_combination (5760 * m * s ^ 2 - 14400 * m * s ^ 4) * hc2

theorem T08.alg_9 : T08.stmt_9 := by
  intro m s c hs hc h
  have hc2 : c ^ 2 = 1 - s ^ 2 := by linarith
  simp only [T08.Z23Data, T08.Cm23Poly]
  field_simp
  linear_combination (5760 * m * s ^ 2 - 14400 * m * s ^ 4) * hc2

namespace T08

lemma kink_facts (m : ℕ) (hm : 2 ≤ m) (j : ℕ → ℕ)
    (hj : ∀ k, j k ∈ Finset.Icc 1 m ∧ ip (n m (j k)) (corner m k) = 0) (k : ℕ) :
    |ip (n m (j k)) (corner m (k + 1))| = 2 * s m ^ 2 ∧
    |ip (n m (j (k + 1))) (corner m k)| = 2 * s m ^ 2 :=
  ⟨xi0_eq m k (j k) hm (hj k).2, xi1_eq m k (j (k + 1)) hm (hj (k + 1)).2⟩

end T08

theorem T08.proof_1 : T08.stmt_1 := T08.geo_1
theorem T08.proof_2 : T08.stmt_2 := T08.geo_2
theorem T08.proof_3 : T08.stmt_3 := T08.geo_3

theorem T08.proof_4 : T08.stmt_4 := by
  intro m hm j hj
  have hs := T08.s_pos m hm
  have hc := T08.c_pos m hm
  have hterm : ∀ k ∈ Finset.range (2 * m),
      T08.Ze (T08.corner m k) (T08.corner m (k + 1)) (T08.n m (j k)) (T08.n m (j (k + 1)))
        = T08.ZData (2 * T08.s m ^ 2) (2 * T08.s m ^ 2) (T08.s m)
          (2 * T08.s m ^ 2 / T08.c m) (T08.s m ^ 2 / T08.c m) := by
    intro k _
    obtain ⟨f0, f1⟩ := T08.kink_facts m hm j hj k
    unfold T08.Ze
    rw [f0, f1, T08.edgeH_eq m k hm, T08.edgeL_eq m k hm, T08.edgeL0_eq m k hm]
  rw [Finset.sum_congr rfl hterm, Finset.sum_const, Finset.card_range, nsmul_eq_mul]
  have := T08.alg_7 (m : ℝ) (T08.s m) (T08.c m) hs hc (T08.sc_sq m)
  push_cast
  rw [this]
  rfl

theorem T08.proof_5 : T08.stmt_5 := by
  intro m hm j hj
  have hs := T08.s_pos m hm
  have hc := T08.c_pos m hm
  have hterm : ∀ k ∈ Finset.range (2 * m),
      T08.Ze5 (T08.corner m k) (T08.corner m (k + 1)) (T08.n m (j k)) (T08.n m (j (k + 1)))
        = T08.Z5Data (2 * T08.s m ^ 2) (2 * T08.s m ^ 2) (T08.s m)
          (2 * T08.s m ^ 2 / T08.c m) (T08.s m ^ 2 / T08.c m) := by
    intro k _
    obtain ⟨f0, f1⟩ := T08.kink_facts m hm j hj k
    unfold T08.Ze5
    rw [f0, f1, T08.edgeH_eq m k hm, T08.edgeL_eq m k hm, T08.edgeL0_eq m k hm]
  rw [Finset.sum_congr rfl hterm, Finset.sum_const, Finset.card_range, nsmul_eq_mul]
  have := T08.alg_8 (m : ℝ) (T08.s m) (T08.c m) hs hc (T08.sc_sq m)
  push_cast
  rw [this]
  rfl

theorem T08.proof_6 : T08.stmt_6 := by
  intro m hm j hj
  have hs := T08.s_pos m hm
  have hc := T08.c_pos m hm
  have hterm : ∀ k ∈ Finset.range (2 * m),
      T08.Ze23 (T08.corner m k) (T08.corner m (k + 1)) (T08.n m (j k)) (T08.n m (j (k + 1)))
        = T08.Z23Data (2 * T08.s m ^ 2) (2 * T08.s m ^ 2) (T08.s m)
          (2 * T08.s m ^ 2 / T08.c m) (T08.s m ^ 2 / T08.c m) := by
    intro k _
    obtain ⟨f0, f1⟩ := T08.kink_facts m hm j hj k
    unfold T08.Ze23
    rw [f0, f1, T08.edgeH_eq m k hm, T08.edgeL_eq m k hm, T08.edgeL0_eq m k hm]
  rw [Finset.sum_congr rfl hterm, Finset.sum_const, Finset.card_range, nsmul_eq_mul]
  have := T08.alg_9 (m : ℝ) (T08.s m) (T08.c m) hs hc (T08.sc_sq m)
  push_cast
  rw [this]
  rfl

theorem T08.proof_7 : T08.stmt_7 := T08.alg_7
theorem T08.proof_8 : T08.stmt_8 := T08.alg_8
theorem T08.proof_9 : T08.stmt_9 := T08.alg_9
