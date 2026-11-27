import Mathlib
import Zeta3.Wave6.Z04_SuperellipsoidArea
/-!
# R04 — the superellipsoid area to every order (paper label `thm:S`, "Moreover")

Target id: R04.
Paper labels: `thm:S` (paper/ieee/main.tex lines 1921–1943), the sentence "Moreover ..." and the
parenthetical "($Q$ is smooth on $[0,\frac14]$, one-sided at $0$)"; detailed proof
paper/ieee/proofs/sec_g3_S.tex, subsection "No logarithms, and $s_1,s_2,s_3$" (lines 89–121:
Lemma `g3:lem:smooth`, Lemma `g3:lem:coeffs`), and line 36 ($Q_k=Q^{(k)}(0)/k!$).
Earlier Lean: Z04 (the exact formula, `S → 24`, the case `K = 1`, `g3:lem:smooth` at `j = 0`),
W02 (`P(p)/8 - 1 = Σ M_k/p^{k+1}`, `M_k = W02.Mk k`), W03 (`M₀`, `M₁`), Y01 (`M₂`).

Quoted verbatim.

main.tex 1921–1943 (`thm:S`, the parts formalized here):

  Let $\varepsilon=1/p$ and $G(\varepsilon)=\Gamma(1+\varepsilon)^2/\Gamma(1+2\varepsilon)$.
  ...
  Moreover $S(p)=\sum_{k\le K}s_kp^{-k}+O(p^{-K-1})$ for every $K$, with no $\log p$ terms,
  $s_0=24$ and
  \begin{gather*}
  s_1=48M_0,\qquad s_2=48M_1-4\pi^2-8Q_0,\\
  s_3=48M_2-8\pi^2M_0+48\zeta(3)-8Q_1,
  \end{gather*}
  where $Q_0=Q(0)$ and $Q_1=Q'(0)$ ($Q$ is smooth on $[0,\frac14]$, one-sided at $0$) and

main.tex 1944:

  The $48\zeta(3)$ in $s_3$ comes from the face factor,
  $G(\varepsilon)=1-\zeta(2)\varepsilon^2+2\zeta(3)\varepsilon^3+\cdots$

sec_g3_S.tex 36:

  The smoothness of $Q$ on $[0,\frac14]$ stated in the theorem is Lemma~\ref{g3:lem:smooth},
  applied to $H_3$ ... $Q_k=Q^{(k)}(0)/k!$.

sec_g3_S.tex, Lemma `g3:lem:smooth`:

  For $m=1,2,3$, $H_m$ is finite and $C^\infty$ on $[0,\frac14]$ (one-sided derivatives at $0$,
  where $\mathrm h$ is defined by the same formula with $\varepsilon=0$).

  ... (in the proof)
  |\partial_\varepsilon^j\mathrm h|\le C_j\,(1+|\log y|)^{j+1}\ (y\le1),\qquad
  |\partial_\varepsilon^j\mathrm h|\le C_jy^{-2}\ (y\ge1),
  uniformly in $\varepsilon\in[0,\frac14]$.

sec_g3_S.tex, Lemma `g3:lem:coeffs`:

  $S(p)=\sum_{k\le K}s_kp^{-k}+O(p^{-K-1})$ for every $K$, with $s_0=24$, $s_1=48M_0$,
  $s_2=48M_1-4\pi^2-8Q_0$, $s_3=48M_2-8\pi^2M_0+48\zeta(3)-8Q_1$.

  ... (in the proof) Taylor's formula with the Lagrange remainder gives the expansion for every
  $K$, with no logarithms, and $s_k$ is $\frac1{k!}$ times the $k$-th derivative at $0$. To
  compute $s_0,\dots,s_3$ use $S=24G+6G(P-8)-8\varepsilon^2Q$ ...
  so $G=1-\zeta(2)\varepsilon^2+2\zeta(3)\varepsilon^3+O(\varepsilon^4)$.

Formalization choices.
* `S = Y05.S`, `P = W02.Pint`, `G = Z04.G`, `h_ε(y) = Z04.h ε y`, `H_m = Z04.H m`,
  `Q = Z04.Q` (all imported, all inside audited statement blocks). `Z04.Q` is defined by the
  theorem's formula for every real `ε`; at `ε = 0` it is "the same formula with $\varepsilon=0$"
  (`Γ(1) = 1`, `x^{ε-1} = x^{-1}`), so no separate extension is needed: the smoothness claims are
  `ContDiffOn ℝ _ Z04.Q [0, 1/4]`, which in Mathlib means one-sided derivatives at the endpoints
  (the form the supplement proves: `C^∞` on `[0,1/4]`).
* `ζ(3) = (riemannZeta 3).re`; `ζ(2)` is written `π²/6` (`riemannZeta_two`).
* `M_k = (W02.Mk k).re`: `W02.Mk` is complex, and W02.stmt_4 shows `F` is real on real `ε < 1`,
  so `M_k` is real; with these `M_k`, W02.stmt_6 gives `P(p)/8 - 1 = Σ_k M_k/p^{k+1}`.
* `Q_j = Q^{(j)}(0)/j!` is the one-sided `iteratedDerivWithin j Z04.Q [0,1/4] 0 / j!`
  (sec_g3_S.tex line 36). `g_k = G^{(k)}(0)/k!` (G is smooth near 0, two-sided).
* The sequence `s_k` (`R04.s`) is defined by the supplement's formula: the `ε^k` coefficient
  of `S = 24G + 6G(P-8) - 8ε²Q` with `G = Σ g_k ε^k`, `P - 8 = 8 Σ_{k≥0} M_k ε^{k+1}`,
  `Q = Σ Q_j ε^j`, i.e.
  `s_k = 24 g_k + 48 Σ_{i<k} g_i M_{k-1-i} - 8 Q_{k-2}` (the last term only for `k ≥ 2`).
* "$O(p^{-K-1})$ ... with no $\log p$ terms" is stated as: for every `K` there are `C`, `p₀`
  with `|S(p) - Σ_{k≤K} s_k p^{-k}| ≤ C p^{-K-1}` for all `p ≥ p₀`, for the one sequence `s`.
* Numerical check before freezing (mpmath, 20 digits; `h`, `Q` by quadrature, `Q₁` by a
  5-point difference, `M₁ = 2.02095276965635835695/8`, `M₂ = -0.02216191108726939942`,
  `M₀ = log 2 - √2 log(1+√2)`): `Q₀ = -2.83293568349`, `Q₁ = 6.50432557372` (paper:
  `6.5043255737`), `s₁ = -26.5585583866`, `s₂ = -4.68921551846`, `s₃ = 48.2874324862` (paper:
  `-26.5585583865`, `-4.6892155184`, `48.2874324860`). The direct face integral of `Y05.S`
  (2D quadrature, independent of the exact formula) gives `S(2) = 12.5663706144 = 4π` and
  `(S(p) - Σ_{k≤3} s_k p^{-k}) p⁴ = -69.84, -73.75, -75.85` at `p = 10, 20, 40`, consistent with
  `s₄ + s₅/p` (`s₄ ≈ -78`, `s₅ ≈ 90.7`): the `p^{-3}` truncation error is `O(p^{-4})`.

Statements (ordered so that partial progress is recorded).
* stmt_1: face factor, `G(ε) = 1 - ζ(2)ε² + 2ζ(3)ε³ + O(ε⁴)` near `0` (explicit constant).
* stmt_2: `G` is `C^∞` on `(-1/2, 1)` and its Taylor coefficients are `g₀ = 1`, `g₁ = 0`,
  `g₂ = -ζ(2)`, `g₃ = 2ζ(3)`.
* stmt_3: `Q` is `C¹` on `[0, 1/4]` (one-sided at `0`).
* stmt_4: `Q` is `C²` on `[0, 1/4]`.
* stmt_5: proof of `g3:lem:smooth`, every `j`: for each `y > 0`, `ε ↦ h_ε(y)` is `C^∞` on
  `[0,1/4]`, and `|∂_ε^j h_ε(y)| ≤ C_j (1+|log y|)^{j+1}` (`y ≤ 1`), `≤ C_j y^{-2}` (`y ≥ 1`),
  uniformly in `ε ∈ [0,1/4]`.
* stmt_6: Lemma `g3:lem:smooth`: for `m = 1, 2, 3`, `h_ε^m` is integrable on `(0,∞)` for every
  `ε ∈ [0,1/4]` ("finite") and `H_m` is `C^∞` on `[0,1/4]`.
* stmt_7: the parenthetical of `thm:S`: `Q` is `C^∞` on `[0, 1/4]`, one-sided at `0`.
* stmt_8: Lemma `g3:lem:coeffs`, the expansion: for every `K`,
  `|S(p) - Σ_{k≤K} s_k p^{-k}| ≤ C p^{-K-1}` for `p ≥ p₀`, with the explicit `s_k` above.
* stmt_9: Lemma `g3:lem:coeffs`, the four numbers: `s₀ = 24`, `s₁ = 48M₀`,
  `s₂ = 48M₁ - 4π² - 8Q₀`, `s₃ = 48M₂ - 8π²M₀ + 48ζ(3) - 8Q₁`, `Q₀ = Q(0)`, `Q₁` the
  one-sided `Q'(0)`.
* stmt_10: `thm:S`, "Moreover", in full (THE TARGET), independent of the definition of `s`:
  `Q` is `C^∞` on `[0,1/4]`, has a right derivative `Q₁` at `0`, and there is a sequence `s`
  with the expansion to every order `K` and `s₀ = 24`, `s₁ = 48M₀`, `s₂ = 48M₁ - 4π² - 8Q(0)`,
  `s₃ = 48M₂ - 8π²M₀ + 48ζ(3) - 8Q₁`.

NOT COVERED:
* The closed form of the vertex constant `Q₀` (the rest of `thm:S`, Section `g3:sec:S:Q0`).
* The numerical values of `s₁, s₂, s₃, Q₁` printed after the theorem, and the sentences about
  where the `ζ(3)` of `s₃` comes from beyond `G` ("$M_2$ contains more, and $Q_1$ may too").
* That `Y05.S` is the Hausdorff area of the surface (as in Z04).
* `S_n(p)` and Proposition `prop:cuberadius`.
-/
-- BEGIN STATEMENTS
namespace R04

open MeasureTheory Real Set Filter Topology
open scoped ContDiff

/-- `ζ(3)` as a real number: `(riemannZeta 3).re`. -/
noncomputable def ζ3 : ℝ := (riemannZeta 3).re

/-- `M_k`, the coefficient of `p^{-(k+1)}` in `P(p)/8 - 1` (W02, eq:expansion), as a real. -/
noncomputable def M (k : ℕ) : ℝ := (W02.Mk k).re

/-- The interval `[0, 1/4]` of Lemma `g3:lem:smooth`. -/
noncomputable def I0 : Set ℝ := Set.Icc 0 (1 / 4)

/-- `g_k = G^{(k)}(0)/k!`, the Taylor coefficients of the face factor `G` at `0`. -/
noncomputable def g (k : ℕ) : ℝ := iteratedDeriv k Z04.G 0 / (Nat.factorial k : ℝ)

/-- `Q_j = Q^{(j)}(0)/j!` (sec_g3_S.tex line 36), one-sided at `0` (derivatives within
`[0,1/4]`). -/
noncomputable def Qc (j : ℕ) : ℝ := iteratedDerivWithin j Z04.Q I0 0 / (Nat.factorial j : ℝ)

/-- The coefficient `s_k` of `ε^k` in `S = 24G + 6G(P-8) - 8ε²Q` (proof of `g3:lem:coeffs`):
`s_k = 24 g_k + 48 Σ_{i<k} g_i M_{k-1-i} - 8 Q_{k-2}` (the last term for `k ≥ 2` only). -/
noncomputable def s (k : ℕ) : ℝ :=
  24 * g k + 48 * ∑ i ∈ Finset.range k, g i * M (k - 1 - i)
    - 8 * (if 2 ≤ k then Qc (k - 2) else 0)

/-- main.tex 1944 / proof of `g3:lem:coeffs`:
"$G=1-\zeta(2)\varepsilon^2+2\zeta(3)\varepsilon^3+O(\varepsilon^4)$", with `ζ(2) = π²/6`. -/
def stmt_1 : Prop :=
  ∃ C δ : ℝ, 0 < δ ∧ ∀ ε : ℝ, |ε| ≤ δ →
    |Z04.G ε - (1 - π ^ 2 / 6 * ε ^ 2 + 2 * ζ3 * ε ^ 3)| ≤ C * ε ^ 4

/-- Proof of `g3:lem:coeffs`: `G` is smooth ("$\Gamma(1+\varepsilon)$, $1/\Gamma(1+2\varepsilon)$
are analytic near $[0,\frac14]$"), and its Taylor coefficients at `0` are `1, 0, -ζ(2), 2ζ(3)`. -/
def stmt_2 : Prop :=
  ContDiffOn ℝ ∞ Z04.G (Set.Ioo (-1 / 2) 1) ∧
    g 0 = 1 ∧ g 1 = 0 ∧ g 2 = -(π ^ 2 / 6) ∧ g 3 = 2 * ζ3

/-- `thm:S` parenthetical / `g3:lem:smooth`, first derivative: `Q` is `C¹` on `[0,1/4]`
(one-sided at `0`). -/
def stmt_3 : Prop :=
  ContDiffOn ℝ 1 Z04.Q I0

/-- `thm:S` parenthetical / `g3:lem:smooth`, second derivative: `Q` is `C²` on `[0,1/4]`. -/
def stmt_4 : Prop :=
  ContDiffOn ℝ 2 Z04.Q I0

/-- Proof of `g3:lem:smooth`: for each `y > 0`, `ε ↦ h_ε(y)` is smooth on `[0,1/4]`, and
"$|\partial_\varepsilon^j\mathrm h|\le C_j\,(1+|\log y|)^{j+1}\ (y\le1)$,
$|\partial_\varepsilon^j\mathrm h|\le C_jy^{-2}\ (y\ge1)$, uniformly in
$\varepsilon\in[0,\frac14]$", for every `j` (one-sided derivatives at `0`). -/
def stmt_5 : Prop :=
  (∀ y : ℝ, 0 < y → ContDiffOn ℝ ∞ (fun e => Z04.h e y) I0) ∧
    ∀ j : ℕ, ∃ C : ℝ, ∀ y : ℝ, 0 < y → ∀ ε ∈ I0,
      (y ≤ 1 → |iteratedDerivWithin j (fun e => Z04.h e y) I0 ε| ≤ C * (1 + |Real.log y|) ^ (j + 1))
        ∧ (1 ≤ y → |iteratedDerivWithin j (fun e => Z04.h e y) I0 ε| ≤ C / y ^ 2)

/-- Lemma `g3:lem:smooth`: "For $m=1,2,3$, $H_m$ is finite and $C^\infty$ on $[0,\frac14]$
(one-sided derivatives at $0$ ...)". "Finite" is stated as integrability of `h_ε^m` on
`(0,∞)` for each `ε ∈ [0,1/4]`. -/
def stmt_6 : Prop :=
  ∀ m : ℕ, 1 ≤ m → m ≤ 3 →
    (∀ ε ∈ I0, IntegrableOn (fun y => Z04.h ε y ^ m) (Set.Ioi 0)) ∧
      ContDiffOn ℝ ∞ (Z04.H m) I0

/-- `thm:S`: "($Q$ is smooth on $[0,\frac14]$, one-sided at $0$)". -/
def stmt_7 : Prop :=
  ContDiffOn ℝ ∞ Z04.Q I0

/-- Lemma `g3:lem:coeffs`: "$S(p)=\sum_{k\le K}s_kp^{-k}+O(p^{-K-1})$ for every $K$" (no
`log p` terms), with the explicit `s_k` of the supplement's formula. -/
def stmt_8 : Prop :=
  ∀ K : ℕ, ∃ C p₀ : ℝ, ∀ p : ℝ, p₀ ≤ p →
    |Y05.S p - ∑ k ∈ Finset.range (K + 1), s k / p ^ k| ≤ C / p ^ (K + 1)

/-- Lemma `g3:lem:coeffs`: "$s_0=24$, $s_1=48M_0$, $s_2=48M_1-4\pi^2-8Q_0$,
$s_3=48M_2-8\pi^2M_0+48\zeta(3)-8Q_1$", with `Q₀ = Q(0)` and `Q₁` the one-sided `Q'(0)`. -/
def stmt_9 : Prop :=
  s 0 = 24 ∧ s 1 = 48 * M 0 ∧
    s 2 = 48 * M 1 - 4 * π ^ 2 - 8 * Z04.Q 0 ∧
    s 3 = 48 * M 2 - 8 * π ^ 2 * M 0 + 48 * ζ3 - 8 * derivWithin Z04.Q I0 0

/-- `thm:S`, "Moreover", in full: `Q` is smooth on `[0,1/4]` with right derivative `Q₁` at
`0`, and "$S(p)=\sum_{k\le K}s_kp^{-k}+O(p^{-K-1})$ for every $K$, with no $\log p$ terms,
$s_0=24$ and $s_1=48M_0$, $s_2=48M_1-4\pi^2-8Q_0$, $s_3=48M_2-8\pi^2M_0+48\zeta(3)-8Q_1$,
where $Q_0=Q(0)$ and $Q_1=Q'(0)$". -/
def stmt_10 : Prop :=
  ContDiffOn ℝ ∞ Z04.Q I0 ∧
    ∃ Q₁ : ℝ, HasDerivWithinAt Z04.Q Q₁ (Set.Ici 0) 0 ∧
      ∃ σ : ℕ → ℝ,
        (∀ K : ℕ, ∃ C p₀ : ℝ, ∀ p : ℝ, p₀ ≤ p →
          |Y05.S p - ∑ k ∈ Finset.range (K + 1), σ k / p ^ k| ≤ C / p ^ (K + 1)) ∧
        σ 0 = 24 ∧ σ 1 = 48 * M 0 ∧
        σ 2 = 48 * M 1 - 4 * π ^ 2 - 8 * Z04.Q 0 ∧
        σ 3 = 48 * M 2 - 8 * π ^ 2 * M 0 + 48 * ζ3 - 8 * Q₁

end R04
-- END STATEMENTS

open MeasureTheory Real Set Filter Topology
open scoped ContDiff

namespace R04

lemma Gamma_an {x : ℝ} (hx : 0 < x) : AnalyticAt ℝ Real.Gamma x := by
  have hopen : IsOpen {s : ℂ | 0 < s.re} := isOpen_lt continuous_const Complex.continuous_re
  have hd : DifferentiableOn ℂ Complex.Gamma {s : ℂ | 0 < s.re} := fun s hs =>
    (Complex.differentiableAt_Gamma s (fun m hm => by
      rw [hm] at hs; simp at hs; linarith [(m.cast_nonneg : (0:ℝ) ≤ m)])).differentiableWithinAt
  have hA : AnalyticAt ℂ Complex.Gamma (x : ℂ) :=
    hd.analyticAt (hopen.mem_nhds (by simpa using hx))
  have h1 : AnalyticAt ℝ (fun t : ℝ => Complex.Gamma (t : ℂ)) x :=
    hA.restrictScalars.comp (Complex.ofRealCLM.analyticAt x)
  have h2 := (Complex.reCLM.analyticAt _).comp h1
  convert h2 using 1
  funext t; simp [Complex.Gamma_ofReal]

lemma Gamma_cd {x : ℝ} (hx : 0 < x) {n : WithTop ℕ∞} : ContDiffAt ℝ n Real.Gamma x :=
  (Gamma_an hx).contDiffAt

lemma Gamma_comp_cd {a b x : ℝ} (hx : 0 < a + b * x) {n : WithTop ℕ∞} :
    ContDiffAt ℝ n (fun e => Real.Gamma (a + b * e)) x :=
  ContDiffAt.comp (g := Real.Gamma) x (Gamma_cd hx)
    (contDiffAt_const.add (contDiffAt_const.mul contDiffAt_id))

lemma G_cdAt {ε : ℝ} (h1 : 0 < 1 + ε) (h2 : 0 < 1 + 2 * ε) {n : WithTop ℕ∞} :
    ContDiffAt ℝ n Z04.G ε := by
  have a1 : ContDiffAt ℝ n (fun e => Real.Gamma (1 + 1 * e)) ε := Gamma_comp_cd (by linarith)
  have a2 : ContDiffAt ℝ n (fun e => Real.Gamma (1 + 2 * e)) ε := Gamma_comp_cd h2
  have := (a1.pow 2).div a2 (Real.Gamma_pos_of_pos h2).ne'
  convert this using 1
  funext e; simp [Z04.G]

lemma G_eq_A (ε : ℝ) (hε : -1 / 2 < ε) : Z04.G ε = X05.A (2 * ε) 0 := by
  rw [X05.proof_2 (2 * ε) 0 (by linarith) (by simp; linarith)]
  unfold Z04.G
  ring_nf

end R04

theorem R04.proof_1 : R04.stmt_1 := by
  obtain ⟨C, δ, hδ, hA⟩ := X05.proof_6
  refine ⟨16 * |(π ^ 2 / 6) ^ 2 / 32 - 7 * (π ^ 4 / 90) / 32| +
      32 * |-(π ^ 2 / 6) * X05.z 3 / 16 + 3 * X05.z 5 / 16| + 64 * |C|, min (δ / 2) (1 / 4),
      by positivity, fun ε hε => ?_⟩
  have hε1 : |ε| ≤ 1 / 4 := le_trans hε (min_le_right _ _)
  have hε2 : |2 * ε| ≤ δ := by
    rw [abs_mul, abs_two]; linarith [min_le_left (δ / 2) (1 / 4)]
  have h := (hA (2 * ε) hε2).1
  rw [← R04.G_eq_A ε (by cases abs_le.mp hε1; linarith)] at h
  have e3 : X05.z 3 = R04.ζ3 := by simp [X05.z, R04.ζ3]
  rw [e3] at h
  set c4 := (π ^ 2 / 6) ^ 2 / 32 - 7 * (π ^ 4 / 90) / 32
  set c5 := -(π ^ 2 / 6) * R04.ζ3 / 16 + 3 * X05.z 5 / 16
  have hc5 : -(π ^ 2 / 6) * X05.z 3 / 16 + 3 * X05.z 5 / 16 = c5 := by rw [e3]
  rw [hc5]
  have ha : |ε| ^ 5 ≤ ε ^ 4 := by
    rw [show ε ^ 4 = |ε| ^ 4 by rw [pow_abs]; exact (abs_of_nonneg (by positivity)).symm]
    calc |ε| ^ 5 = |ε| ^ 4 * |ε| := by ring
      _ ≤ |ε| ^ 4 * 1 := by gcongr; linarith
      _ = _ := by ring
  have hb : |ε| ^ 6 ≤ ε ^ 4 := by
    rw [show ε ^ 4 = |ε| ^ 4 by rw [pow_abs]; exact (abs_of_nonneg (by positivity)).symm]
    calc |ε| ^ 6 = |ε| ^ 4 * |ε| ^ 2 := by ring
      _ ≤ |ε| ^ 4 * 1 := by gcongr; nlinarith [abs_nonneg ε]
      _ = _ := by ring
  have e : Z04.G ε - (1 - π ^ 2 / 6 * ε ^ 2 + 2 * R04.ζ3 * ε ^ 3) =
      (Z04.G ε - (1 - (π ^ 2 / 6) / 4 * (2 * ε) ^ 2 + R04.ζ3 / 4 * (2 * ε) ^ 3
        + c4 * (2 * ε) ^ 4 + c5 * (2 * ε) ^ 5)) + 16 * c4 * ε ^ 4 + 32 * c5 * ε ^ 5 := by ring
  rw [e]
  have h6 : |2 * ε| ^ 6 = 64 * |ε| ^ 6 := by rw [abs_mul, abs_two]; ring
  rw [h6] at h
  have t1 : |16 * c4 * ε ^ 4| ≤ 16 * |c4| * ε ^ 4 := by
    rw [abs_mul, abs_mul, abs_of_pos (by norm_num : (0:ℝ) < 16),
      abs_of_nonneg (by positivity : (0:ℝ) ≤ ε ^ 4)]
  have t2 : |32 * c5 * ε ^ 5| ≤ 32 * |c5| * ε ^ 4 := by
    rw [abs_mul, abs_mul, abs_of_pos (by norm_num : (0:ℝ) < 32), abs_pow]
    gcongr
  have t3 : C * (64 * |ε| ^ 6) ≤ 64 * |C| * ε ^ 4 := by
    have := le_abs_self C
    calc C * (64 * |ε| ^ 6) ≤ |C| * (64 * |ε| ^ 6) := by gcongr
      _ ≤ |C| * (64 * ε ^ 4) := by gcongr
      _ = _ := by ring
  calc _ ≤ |Z04.G ε - (1 - (π ^ 2 / 6) / 4 * (2 * ε) ^ 2 + R04.ζ3 / 4 * (2 * ε) ^ 3
        + c4 * (2 * ε) ^ 4 + c5 * (2 * ε) ^ 5)| + |16 * c4 * ε ^ 4| + |32 * c5 * ε ^ 5| :=
        abs_add_three _ _ _
    _ ≤ _ := by nlinarith

namespace R04

/-- asymptotic expansion to every order on `(0, δ]` -/
def AE (f : ℝ → ℝ) (a : ℕ → ℝ) : Prop :=
  ∀ K : ℕ, ∃ C δ : ℝ, 0 < δ ∧ ∀ ε : ℝ, 0 < ε → ε ≤ δ →
    |f ε - ∑ k ∈ Finset.range (K + 1), a k * ε ^ k| ≤ C * ε ^ (K + 1)

lemma taylor_bd {f : ℝ → ℝ} {U : Set ℝ} (hU : IsOpen U) {δ : ℝ} (hδ : 0 < δ)
    (hsub : Icc 0 δ ⊆ U) (hf : ContDiffOn ℝ ∞ f U) (K : ℕ) :
    ∃ C : ℝ, ∀ ε ∈ Icc 0 δ,
      |f ε - ∑ k ∈ Finset.range (K + 1), iteratedDeriv k f 0 / (k.factorial : ℝ) * ε ^ k|
        ≤ C * ε ^ (K + 1) := by
  have hI : ContDiffOn ℝ ((K + 1 : ℕ) : WithTop ℕ∞) f (Icc 0 δ) :=
    (hf.of_le (by exact_mod_cast le_top)).mono hsub
  have hc : ContinuousOn (iteratedDerivWithin (K + 1) f (Icc 0 δ)) (Icc 0 δ) :=
    hI.continuousOn_iteratedDerivWithin le_rfl (uniqueDiffOn_Icc hδ)
  obtain ⟨C, hC⟩ := isCompact_Icc.exists_bound_of_continuousOn hc
  refine ⟨C / K.factorial, fun ε hε => ?_⟩
  have := taylor_mean_remainder_bound (n := K) (by simpa using hI) hε (fun y hy => hC y hy)
  rw [taylor_within_apply] at this
  have h0 : (0 : ℝ) ∈ Icc 0 δ := ⟨le_rfl, hδ.le⟩
  have e : ∑ k ∈ Finset.range (K + 1),
      ((k.factorial : ℝ)⁻¹ * (ε - 0) ^ k) • iteratedDerivWithin k f (Icc 0 δ) 0 =
      ∑ k ∈ Finset.range (K + 1), iteratedDeriv k f 0 / (k.factorial : ℝ) * ε ^ k := by
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_Icc hδ)
      ((hf.contDiffAt (hU.mem_nhds (hsub h0))).of_le (by exact_mod_cast le_top)) h0]
    simp only [sub_zero, smul_eq_mul]; ring
  rw [e, Real.norm_eq_abs, sub_zero] at this
  calc _ ≤ C * ε ^ (K + 1) / K.factorial := this
    _ = _ := by ring

lemma d0_zero {K : ℕ} {d : ℕ → ℝ} {C δ : ℝ} (hδ : 0 < δ)
    (h : ∀ ε, 0 < ε → ε ≤ δ → |∑ k ∈ Finset.range (K + 1), d k * ε ^ k| ≤ C * ε ^ (K + 1)) :
    d 0 = 0 := by
  have t1 : Tendsto (fun ε : ℝ => |∑ k ∈ Finset.range (K + 1), d k * ε ^ k|) (𝓝[>] 0)
      (𝓝 |d 0|) := by
    have hc : Continuous (fun ε : ℝ => |∑ k ∈ Finset.range (K + 1), d k * ε ^ k|) := by
      fun_prop
    have h0 : Tendsto (fun ε : ℝ => |∑ k ∈ Finset.range (K + 1), d k * ε ^ k|) (𝓝 0)
        (𝓝 |d 0|) := by
      have := hc.tendsto 0
      convert this using 2
      simp [Finset.sum_range_succ']
    exact h0.mono_left nhdsWithin_le_nhds
  have t2 : Tendsto (fun ε : ℝ => C * ε ^ (K + 1)) (𝓝[>] 0) (𝓝 0) := by
    have hc : Continuous (fun ε : ℝ => C * ε ^ (K + 1)) := by fun_prop
    simpa using (hc.tendsto 0).mono_left nhdsWithin_le_nhds
  have hle : |d 0| ≤ 0 := le_of_tendsto_of_tendsto t1 t2 (by
    filter_upwards [Ioo_mem_nhdsGT hδ] with ε hε
    exact h ε hε.1 hε.2.le)
  exact abs_nonpos_iff.mp hle

lemma coeff_unique : ∀ (K : ℕ) (d : ℕ → ℝ) (C δ : ℝ), 0 < δ →
    (∀ ε, 0 < ε → ε ≤ δ → |∑ k ∈ Finset.range (K + 1), d k * ε ^ k| ≤ C * ε ^ (K + 1)) →
    ∀ k ≤ K, d k = 0 := by
  intro K
  induction K with
  | zero =>
    intro d C δ hδ h k hk
    obtain rfl : k = 0 := by omega
    exact d0_zero hδ h
  | succ K ih =>
    intro d C δ hδ h k hk
    have h0 : d 0 = 0 := d0_zero hδ h
    have h' : ∀ ε, 0 < ε → ε ≤ δ →
        |∑ k ∈ Finset.range (K + 1), d (k + 1) * ε ^ k| ≤ C * ε ^ (K + 1) := by
      intro ε hε hεδ
      have := h ε hε hεδ
      rw [Finset.sum_range_succ', h0, zero_mul, add_zero] at this
      have e : ∑ k ∈ Finset.range (K + 1), d (k + 1) * ε ^ (k + 1) =
          ε * ∑ k ∈ Finset.range (K + 1), d (k + 1) * ε ^ k := by
        rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun k _ => by ring
      rw [e, abs_mul, abs_of_pos hε,
        show C * ε ^ (K + 1 + 1) = ε * (C * ε ^ (K + 1)) by ring] at this
      exact le_of_mul_le_mul_left this hε
    rcases k with _ | k
    · exact h0
    · exact ih (fun k => d (k + 1)) C δ hδ h' k (by omega)

lemma sum_abs_le {a : ℕ → ℝ} {K : ℕ} {ε : ℝ} (h0 : 0 ≤ ε) (h1 : ε ≤ 1) :
    |∑ k ∈ Finset.range (K + 1), a k * ε ^ k| ≤ ∑ k ∈ Finset.range (K + 1), |a k| := by
  refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun k _ => ?_)
  rw [abs_mul, abs_of_nonneg (pow_nonneg h0 k)]
  exact mul_le_of_le_one_right (abs_nonneg _) (pow_le_one₀ h0 h1)

lemma AE.add {f g : ℝ → ℝ} {a b : ℕ → ℝ} (hf : AE f a) (hg : AE g b) :
    AE (fun e => f e + g e) (fun k => a k + b k) := by
  intro K
  obtain ⟨C1, δ1, h1, hf⟩ := hf K
  obtain ⟨C2, δ2, h2, hg⟩ := hg K
  refine ⟨C1 + C2, min δ1 δ2, lt_min h1 h2, fun ε hε hεδ => ?_⟩
  have e : f ε + g ε - ∑ k ∈ Finset.range (K + 1), (a k + b k) * ε ^ k =
      (f ε - ∑ k ∈ Finset.range (K + 1), a k * ε ^ k) +
        (g ε - ∑ k ∈ Finset.range (K + 1), b k * ε ^ k) := by
    simp only [add_mul, Finset.sum_add_distrib]; ring
  rw [e]
  calc _ ≤ _ := abs_add_le _ _
    _ ≤ C1 * ε ^ (K + 1) + C2 * ε ^ (K + 1) :=
        add_le_add (hf ε hε (hεδ.trans (min_le_left _ _)))
          (hg ε hε (hεδ.trans (min_le_right _ _)))
    _ = _ := by ring

lemma AE.const_mul {f : ℝ → ℝ} {a : ℕ → ℝ} (c : ℝ) (hf : AE f a) :
    AE (fun e => c * f e) (fun k => c * a k) := by
  intro K
  obtain ⟨C1, δ1, h1, hf⟩ := hf K
  refine ⟨|c| * C1, δ1, h1, fun ε hε hεδ => ?_⟩
  have e : c * f ε - ∑ k ∈ Finset.range (K + 1), c * a k * ε ^ k =
      c * (f ε - ∑ k ∈ Finset.range (K + 1), a k * ε ^ k) := by
    rw [mul_sub, Finset.mul_sum]; congr 1; exact Finset.sum_congr rfl fun k _ => by ring
  rw [e, abs_mul, mul_assoc]
  exact mul_le_mul_of_nonneg_left (hf ε hε hεδ) (abs_nonneg c)

lemma prod_trunc (a b : ℕ → ℝ) (K : ℕ) {ε : ℝ} (h0 : 0 ≤ ε) (h1 : ε ≤ 1) :
    |(∑ k ∈ Finset.range (K + 1), a k * ε ^ k) * (∑ k ∈ Finset.range (K + 1), b k * ε ^ k)
      - ∑ k ∈ Finset.range (K + 1), (∑ i ∈ Finset.range (k + 1), a i * b (k - i)) * ε ^ k|
      ≤ (∑ k ∈ Finset.range (K + 1), |a k|) * (∑ k ∈ Finset.range (K + 1), |b k|)
        * ε ^ (K + 1) := by
  set F : ℕ → ℕ → ℝ := fun i j => a i * b j * ε ^ i * ε ^ j with hF
  have eT : ∑ k ∈ Finset.range (K + 1), (∑ i ∈ Finset.range (k + 1), a i * b (k - i)) * ε ^ k
      = ∑ m ∈ Finset.range (K + 1), ∑ j ∈ Finset.range (K + 1 - m), F m j := by
    rw [← Finset.sum_range_diag_flip]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl fun i hi => ?_
    have hik : i ≤ k := Nat.lt_succ_iff.mp (Finset.mem_range.mp hi)
    simp only [hF]
    rw [show ε ^ k = ε ^ i * ε ^ (k - i) by rw [← pow_add, Nat.add_sub_cancel' hik]]
    ring
  have eP : (∑ k ∈ Finset.range (K + 1), a k * ε ^ k) * (∑ k ∈ Finset.range (K + 1), b k * ε ^ k)
      = ∑ m ∈ Finset.range (K + 1), ∑ j ∈ Finset.range (K + 1), F m j := by
    rw [Finset.sum_mul_sum]
    refine Finset.sum_congr rfl fun m _ => Finset.sum_congr rfl fun j _ => ?_
    simp only [hF]; ring
  have eD : ∑ m ∈ Finset.range (K + 1), ∑ j ∈ Finset.range (K + 1), F m j
      - ∑ m ∈ Finset.range (K + 1), ∑ j ∈ Finset.range (K + 1 - m), F m j
      = ∑ m ∈ Finset.range (K + 1), ∑ j ∈ Finset.Ico (K + 1 - m) (K + 1), F m j := by
    rw [← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun m _ => ?_
    rw [← Finset.sum_range_add_sum_Ico _ (Nat.sub_le (K + 1) m)]
    ring
  rw [eT, eP, eD, Finset.sum_mul_sum, Finset.sum_mul]
  refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun m hm => ?_)
  rw [Finset.sum_mul]
  refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
  refine le_trans (Finset.sum_le_sum fun j hj => (?_ : |F m j| ≤ |a m| * |b j| * ε ^ (K + 1)))
    (Finset.sum_le_sum_of_subset_of_nonneg (fun j hj => ?_) (fun j _ _ => by positivity))
  · simp only [hF]
    rw [abs_mul, abs_mul, abs_mul, abs_of_nonneg (pow_nonneg h0 _),
      abs_of_nonneg (pow_nonneg h0 _), mul_assoc, ← pow_add]
    have : K + 1 ≤ m + j := by
      have := (Finset.mem_Ico.mp hj).1; omega
    exact mul_le_mul_of_nonneg_left (pow_le_pow_of_le_one h0 h1 this) (by positivity)
  · exact Finset.mem_range.mpr (Finset.mem_Ico.mp hj).2

lemma AE.mul {f g : ℝ → ℝ} {a b : ℕ → ℝ} (hf : AE f a) (hg : AE g b) :
    AE (fun e => f e * g e) (fun k => ∑ i ∈ Finset.range (k + 1), a i * b (k - i)) := by
  intro K
  obtain ⟨C1, δ1, h1, hf⟩ := hf K
  obtain ⟨C2, δ2, h2, hg⟩ := hg K
  set SA := ∑ k ∈ Finset.range (K + 1), |a k|
  set SB := ∑ k ∈ Finset.range (K + 1), |b k|
  refine ⟨|C1| * (SB + |C2|) + SA * |C2| + SA * SB, min (min δ1 δ2) 1,
    lt_min (lt_min h1 h2) one_pos, fun ε hε hεδ => ?_⟩
  have hε1 : ε ≤ δ1 := hεδ.trans ((min_le_left _ _).trans (min_le_left _ _))
  have hε2 : ε ≤ δ2 := hεδ.trans ((min_le_left _ _).trans (min_le_right _ _))
  have hε3 : ε ≤ 1 := hεδ.trans (min_le_right _ _)
  set A := ∑ k ∈ Finset.range (K + 1), a k * ε ^ k
  set B := ∑ k ∈ Finset.range (K + 1), b k * ε ^ k
  set T := ∑ k ∈ Finset.range (K + 1), (∑ i ∈ Finset.range (k + 1), a i * b (k - i)) * ε ^ k
  have hA : |A| ≤ SA := sum_abs_le hε.le hε3
  have hB : |B| ≤ SB := sum_abs_le hε.le hε3
  have hAB := prod_trunc a b K hε.le hε3
  have r1 := hf ε hε hε1
  have r2 := hg ε hε hε2
  have hp : 0 < ε ^ (K + 1) := pow_pos hε _
  have hp1 : ε ^ (K + 1) ≤ 1 := pow_le_one₀ hε.le hε3
  have r1' : |f ε - A| ≤ |C1| * ε ^ (K + 1) :=
    r1.trans (mul_le_mul_of_nonneg_right (le_abs_self _) hp.le)
  have r2' : |g ε - B| ≤ |C2| * ε ^ (K + 1) :=
    r2.trans (mul_le_mul_of_nonneg_right (le_abs_self _) hp.le)
  have hg' : |g ε| ≤ SB + |C2| := by
    have : |g ε| ≤ |B| + |g ε - B| := by
      calc |g ε| = |B + (g ε - B)| := by ring_nf
        _ ≤ _ := abs_add_le _ _
    nlinarith [abs_nonneg C2]
  have e : f ε * g ε - T = (f ε - A) * g ε + A * (g ε - B) + (A * B - T) := by ring
  rw [e]
  have q1 : |(f ε - A) * g ε| ≤ |C1| * ε ^ (K + 1) * (SB + |C2|) := by
    rw [abs_mul]; exact mul_le_mul r1' hg' (abs_nonneg _) (by positivity)
  have q2 : |A * (g ε - B)| ≤ SA * (|C2| * ε ^ (K + 1)) := by
    rw [abs_mul]; exact mul_le_mul hA r2' (abs_nonneg _) ((abs_nonneg _).trans hA)
  calc _ ≤ |(f ε - A) * g ε| + |A * (g ε - B)| + |A * B - T| := abs_add_three _ _ _
    _ ≤ |C1| * ε ^ (K + 1) * (SB + |C2|) + SA * (|C2| * ε ^ (K + 1)) + SA * SB * ε ^ (K + 1) :=
        add_le_add (add_le_add q1 q2) hAB
    _ = _ := by ring

lemma AE_sq : AE (fun e => e ^ 2) (fun k => if k = 2 then 1 else 0) := by
  intro K
  refine ⟨1, 1, one_pos, fun ε hε hε1 => ?_⟩
  have e : ∑ k ∈ Finset.range (K + 1), (if k = 2 then (1 : ℝ) else 0) * ε ^ k =
      if 2 < K + 1 then ε ^ 2 else 0 := by
    simp only [ite_mul, one_mul, zero_mul]
    rw [Finset.sum_ite_eq' (Finset.range (K + 1)) 2 (fun k => ε ^ k)]
    simp [Finset.mem_range]
  rw [e]
  rcases K with _ | _ | K
  · simp; nlinarith
  · simp
  · simp only [show 2 < K + 2 + 1 by omega, ite_true, sub_self, abs_zero, one_mul]
    positivity

lemma chain_smooth {U : Set ℝ} (hU : IsOpen U) (F : ℕ → ℝ → ℝ)
    (h : ∀ n, ∀ ε ∈ U, HasDerivAt (F n) (F (n + 1) ε) ε) :
    ContDiffOn ℝ ∞ (F 0) U ∧ ∀ n, Set.EqOn (iteratedDeriv n (F 0)) (F n) U := by
  have key : ∀ N : ℕ, ∀ k, ContDiffOn ℝ N (F k) U := by
    intro N
    induction N with
    | zero =>
      intro k
      exact contDiffOn_zero.mpr (fun ε hε => (h k ε hε).continuousAt.continuousWithinAt)
    | succ N ih =>
      intro k
      rw [show ((N + 1 : ℕ) : WithTop ℕ∞) = (N : WithTop ℕ∞) + 1 by push_cast; rfl,
        contDiffOn_succ_iff_deriv_of_isOpen hU]
      refine ⟨fun ε hε => (h k ε hε).differentiableAt.differentiableWithinAt,
        fun hN => by simp at hN, ?_⟩
      exact (ih (k + 1)).congr (fun ε hε => (h k ε hε).deriv)
  refine ⟨contDiffOn_infty.mpr fun N => key N 0, ?_⟩
  intro n
  induction n with
  | zero => intro ε _; rfl
  | succ n ih =>
    intro ε hε
    rw [iteratedDeriv_succ]
    have hev : iteratedDeriv n (F 0) =ᶠ[𝓝 ε] F n := Filter.eventuallyEq_of_mem (hU.mem_nhds hε) ih
    rw [hev.deriv_eq]; exact (h n ε hε).deriv

lemma integral_chain {α : Type*} [MeasurableSpace α] {μ : Measure α} {U : Set ℝ} (hU : IsOpen U)
    (D : ℕ → ℝ → α → ℝ)
    (hm : ∀ ε ∈ U, AEStronglyMeasurable (D 0 ε) μ)
    (hd : ∀ᵐ x ∂μ, ∀ n, ∀ ε ∈ U, HasDerivAt (fun e => D n e x) (D (n + 1) ε x) ε)
    (hb : ∀ n, ∃ B : α → ℝ, Integrable B μ ∧ ∀ᵐ x ∂μ, ∀ ε ∈ U, |D n ε x| ≤ B x) :
    ∀ n, ∀ ε ∈ U, HasDerivAt (fun e => ∫ x, D n e x ∂μ) (∫ x, D (n + 1) ε x ∂μ) ε := by
  classical
  have hmeas : ∀ n, ∀ ε ∈ U, AEStronglyMeasurable (D n ε) μ := by
    intro n
    induction n with
    | zero => exact hm
    | succ n ih =>
      intro ε hε
      let f : ℝ → α → ℝ := fun t x =>
        if ε + t ∈ U then t⁻¹ * (D n (ε + t) x - D n ε x) else 0
      have hf : ∀ t, AEStronglyMeasurable (f t) μ := by
        intro t
        by_cases ht : ε + t ∈ U
        · simp only [f, ht, ite_true]; exact ((ih _ ht).sub (ih ε hε)).const_mul _
        · simp only [f, ht, ite_false]; exact aestronglyMeasurable_const
      refine aestronglyMeasurable_of_tendsto_ae (𝓝[≠] (0:ℝ)) hf ?_
      filter_upwards [hd] with x hx
      have h1 := hasDerivAt_iff_tendsto_slope_zero.mp (hx n ε hε)
      refine h1.congr' ?_
      have hev : ∀ᶠ t in 𝓝[≠] (0:ℝ), ε + t ∈ U := by
        have ht : Tendsto (fun t => ε + t) (𝓝 0) (𝓝 ε) := by
          have : Tendsto (fun t => ε + t) (𝓝 0) (𝓝 (ε + 0)) := tendsto_const_nhds.add tendsto_id
          simpa using this
        exact (ht.eventually (hU.mem_nhds hε)).filter_mono nhdsWithin_le_nhds
      filter_upwards [hev] with t ht
      simp [f, ht, smul_eq_mul]
  intro n ε hε
  obtain ⟨B, hBi, hB⟩ := hb (n + 1)
  obtain ⟨B0, hB0i, hB0⟩ := hb n
  obtain ⟨r, hr, hball⟩ := Metric.isOpen_iff.mp hU ε hε
  have := hasDerivAt_integral_of_dominated_loc_of_deriv_le (μ := μ)
    (F := fun e x => D n e x) (F' := fun e x => D (n + 1) e x) (x₀ := ε) (bound := B)
    (s := Metric.ball ε r)
    (Metric.ball_mem_nhds ε hr)
    (Filter.eventually_of_mem (hU.mem_nhds hε) (fun e he => hmeas n e he))
    (hB0i.mono' (hmeas n ε hε) (by
      filter_upwards [hB0] with x hx; simpa [Real.norm_eq_abs] using hx ε hε))
    (hmeas (n + 1) ε hε)
    (by filter_upwards [hB] with x hx e he; simpa [Real.norm_eq_abs] using hx e (hball he))
    hBi
    (by filter_upwards [hd] with x hx e he; exact hx n e (hball he))
  exact this.2

end R04

theorem R04.proof_2 : R04.stmt_2 := by
  have hU : IsOpen (Ioo (-1 / 2 : ℝ) 1) := isOpen_Ioo
  have hcd : ContDiffOn ℝ ∞ Z04.G (Ioo (-1 / 2) 1) := fun ε hε =>
    (R04.G_cdAt (by linarith [hε.1]) (by linarith [hε.1])).contDiffWithinAt
  refine ⟨hcd, ?_⟩
  have g0 : R04.g 0 = 1 := by simp [R04.g, Z04.G]
  obtain ⟨C, hC⟩ := R04.taylor_bd hU (δ := 1 / 4) (by norm_num)
    (fun x hx => ⟨by linarith [hx.1], by linarith [hx.2]⟩) hcd 3
  obtain ⟨C', δ', hδ', h1⟩ := R04.proof_1
  set c : ℕ → ℝ := fun k => if k = 0 then 1 else if k = 2 then -(π ^ 2 / 6) else
    if k = 3 then 2 * R04.ζ3 else 0 with hc
  have key := R04.coeff_unique 3 (fun k => R04.g k - c k) (C + C') (min δ' (1 / 4))
    (lt_min hδ' (by norm_num)) (by
      intro ε hε hεδ
      have a1 := hC ε ⟨hε.le, hεδ.trans (min_le_right _ _)⟩
      have a2 := h1 ε (by rw [abs_of_pos hε]; exact hεδ.trans (min_le_left _ _))
      have e : ∑ k ∈ Finset.range (3 + 1), (R04.g k - c k) * ε ^ k =
          -(Z04.G ε - ∑ k ∈ Finset.range (3 + 1),
            iteratedDeriv k Z04.G 0 / (k.factorial : ℝ) * ε ^ k)
          + (Z04.G ε - (1 - π ^ 2 / 6 * ε ^ 2 + 2 * R04.ζ3 * ε ^ 3)) := by
        simp only [Finset.sum_range_succ, Finset.sum_range_zero, hc, R04.g]
        norm_num
        ring
      rw [e]
      calc _ ≤ _ := abs_add_le _ _
        _ ≤ C * ε ^ (3 + 1) + C' * ε ^ 4 := by rw [abs_neg]; exact add_le_add a1 a2
        _ = _ := by ring)
  have k1 := key 1 (by norm_num)
  have k2 := key 2 (by norm_num)
  have k3 := key 3 (by norm_num)
  simp only [hc] at k1 k2 k3
  norm_num at k1 k2 k3
  exact ⟨g0, by linarith, by linarith, by linarith⟩

namespace R04

/-- the polynomials `R_n` with `∂_ε^n h-integrand = log^n x · x^{ε-1} e^{-x} Ψ_n(u)` -/
noncomputable def Rp : ℕ → Polynomial ℝ
  | 0 => 0
  | n + 1 => -Rp n - 2 + 2 * Polynomial.X * Rp n - 2 * Polynomial.X * Polynomial.derivative (Rp n)

/-- `Ψ_n(u) = 1 - e^{-u} + u e^{-u} R_n(u)` -/
noncomputable def Ψ (n : ℕ) (u : ℝ) : ℝ :=
  1 - Real.exp (-u) + u * Real.exp (-u) * (Rp n).eval u

/-- `Ψ_n'` -/
noncomputable def Ψd (n : ℕ) (u : ℝ) : ℝ :=
  Real.exp (-u) + Real.exp (-u) * (Rp n).eval u - u * Real.exp (-u) * (Rp n).eval u
    + u * Real.exp (-u) * (Polynomial.derivative (Rp n)).eval u

/-- the `n`-th `ε`-derivative of the integrand of `h_ε(y)` -/
noncomputable def Dh (n : ℕ) (ε y x : ℝ) : ℝ :=
  Real.log x ^ n * x ^ (ε - 1) * Real.exp (-x) * Ψ n (x ^ (2 - 2 * ε) / y ^ 2)

lemma Ψ_hasDeriv (n : ℕ) (u : ℝ) : HasDerivAt (Ψ n) (Ψd n u) u := by
  have h1 : HasDerivAt (fun u : ℝ => Real.exp (-u)) (-Real.exp (-u)) u := by
    simpa using (hasDerivAt_neg u).exp
  have h2 := ((hasDerivAt_id' u).mul h1).mul (Polynomial.hasDerivAt (Rp n) u)
  have h3 := (h1.const_sub 1).add h2
  unfold Ψ Ψd
  exact h3.congr_deriv (by simp only [Pi.mul_apply]; ring)

lemma Ψ_succ (n : ℕ) (u : ℝ) : Ψ (n + 1) u = Ψ n u - 2 * u * Ψd n u := by
  unfold Ψ Ψd
  simp only [Rp, Polynomial.eval_sub, Polynomial.eval_add, Polynomial.eval_neg,
    Polynomial.eval_mul, Polynomial.eval_X, Polynomial.eval_ofNat]
  ring

lemma Dh_zero (ε y : ℝ) : (fun x => Dh 0 ε y x) = Z04.hIntegrand ε y := by
  funext x
  simp [Dh, Ψ, Rp, Z04.hIntegrand]
  ring

lemma Dh_hasDeriv (n : ℕ) {x : ℝ} (hx : 0 < x) (y ε : ℝ) :
    HasDerivAt (fun e => Dh n e y x) (Dh (n + 1) ε y x) ε := by
  have a : HasDerivAt (fun e : ℝ => x ^ (e - 1)) (Real.log x * 1 * x ^ (ε - 1)) ε :=
    ((hasDerivAt_id ε).sub_const 1).const_rpow hx
  have b : HasDerivAt (fun e : ℝ => x ^ (2 - 2 * e)) (Real.log x * (-(2 * 1)) * x ^ (2 - 2 * ε)) ε :=
    (((hasDerivAt_id ε).const_mul 2).const_sub 2).const_rpow hx
  have bu := b.div_const (y ^ 2)
  have c := (Ψ_hasDeriv n (x ^ (2 - 2 * ε) / y ^ 2)).comp ε bu
  have := ((a.const_mul (Real.log x ^ n)).mul_const (Real.exp (-x))).mul c
  have e : Dh (n + 1) ε y x = Real.log x ^ n * (Real.log x * 1 * x ^ (ε - 1)) * Real.exp (-x) *
      Ψ n (x ^ (2 - 2 * ε) / y ^ 2) + Real.log x ^ n * x ^ (ε - 1) * Real.exp (-x) *
      (Ψd n (x ^ (2 - 2 * ε) / y ^ 2) *
        (Real.log x * (-(2 * 1)) * x ^ (2 - 2 * ε) / y ^ 2)) := by
    unfold Dh; rw [Ψ_succ]; ring
  rw [e]; exact this

lemma exp_poly_bd (P : Polynomial ℝ) :
    ∃ c : ℝ, 0 ≤ c ∧ ∀ u : ℝ, 0 ≤ u → |Real.exp (-u) * P.eval u| ≤ c := by
  refine ⟨∑ i ∈ Finset.range (P.natDegree + 1), |P.coeff i| * (i.factorial : ℝ), by positivity,
    fun u hu => ?_⟩
  rw [Polynomial.eval_eq_sum_range, Finset.mul_sum]
  refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun i _ => ?_)
  have h1 := Real.pow_div_factorial_le_exp u hu i
  have hf : (0 : ℝ) < i.factorial := by exact_mod_cast i.factorial_pos
  have h2 : u ^ i * Real.exp (-u) ≤ i.factorial := by
    rw [div_le_iff₀ hf] at h1
    calc u ^ i * Real.exp (-u) ≤ (Real.exp u * i.factorial) * Real.exp (-u) := by gcongr
      _ = i.factorial := by
        rw [mul_comm (Real.exp u), mul_assoc, ← Real.exp_add]; simp
  rw [abs_mul, abs_mul, abs_of_pos (Real.exp_pos _), abs_of_nonneg (pow_nonneg hu i)]
  calc Real.exp (-u) * (|P.coeff i| * u ^ i) = |P.coeff i| * (u ^ i * Real.exp (-u)) := by ring
    _ ≤ |P.coeff i| * i.factorial := by gcongr

lemma Ψ_bd (n : ℕ) : ∃ c : ℝ, 0 ≤ c ∧ ∀ u : ℝ, 0 ≤ u → ∀ θ : ℝ, 0 ≤ θ → θ ≤ 1 →
    |Ψ n u| ≤ c * u ^ θ := by
  obtain ⟨c1, hc1, h1⟩ := exp_poly_bd (Rp n)
  obtain ⟨c2, hc2, h2⟩ := exp_poly_bd (Polynomial.X * Rp n)
  refine ⟨1 + c1 + c2, by positivity, fun u hu θ hθ0 hθ1 => ?_⟩
  have e1 : 0 ≤ 1 - Real.exp (-u) := Z04.one_sub_exp_neg_nn hu
  have e2 : 1 - Real.exp (-u) ≤ u := Z04.one_sub_exp_neg_le
  have e3 : 1 - Real.exp (-u) ≤ 1 := Z04.one_sub_exp_neg_le_one u
  have eA : Ψ n u = (1 - Real.exp (-u)) + u * (Real.exp (-u) * (Rp n).eval u) := by
    unfold Ψ; ring
  have eB : Ψ n u = (1 - Real.exp (-u)) + Real.exp (-u) * (Polynomial.X * Rp n).eval u := by
    unfold Ψ; rw [Polynomial.eval_mul, Polynomial.eval_X]; ring
  have hA : |Ψ n u| ≤ (1 + c1) * u := by
    have h1' := h1 u hu
    have t : |u * (Real.exp (-u) * (Rp n).eval u)| ≤ c1 * u := by
      rw [abs_mul, abs_of_nonneg hu]
      calc u * |Real.exp (-u) * (Rp n).eval u| ≤ u * c1 := by gcongr
        _ = c1 * u := by ring
    have t2 := abs_add_le (1 - Real.exp (-u)) (u * (Real.exp (-u) * (Rp n).eval u))
    rw [← eA, abs_of_nonneg e1] at t2
    linarith
  have hB : |Ψ n u| ≤ 1 + c2 := by
    have t2 := abs_add_le (1 - Real.exp (-u)) (Real.exp (-u) * (Polynomial.X * Rp n).eval u)
    rw [← eB, abs_of_nonneg e1] at t2
    linarith [h2 u hu]
  rcases le_total u 1 with hu1 | hu1
  · have t : u ≤ u ^ θ := by
      have := Real.rpow_le_rpow_of_exponent_ge' hu hu1 hθ0 hθ1
      simpa using this
    calc |Ψ n u| ≤ (1 + c1) * u := hA
      _ ≤ (1 + c1 + c2) * u ^ θ := by
        apply mul_le_mul (by linarith) t hu (by positivity)
  · have t : 1 ≤ u ^ θ := Real.one_le_rpow hu1 hθ0
    calc |Ψ n u| ≤ 1 + c2 := hB
      _ ≤ (1 + c1 + c2) * 1 := by linarith
      _ ≤ (1 + c1 + c2) * u ^ θ := by gcongr

lemma pow_alg {x y : ℝ} (hx : 0 < x) (hy : 0 < y) (ε θ : ℝ) :
    x ^ (ε - 1) * (x ^ (2 - 2 * ε) / y ^ 2) ^ θ = x ^ (θ * (2 - 2 * ε) + ε - 1) / y ^ (2 * θ) := by
  rw [Real.div_rpow (by positivity) (by positivity), ← Real.rpow_mul hx.le,
    show (y ^ 2 : ℝ) = y ^ (2 : ℝ) by norm_cast, ← Real.rpow_mul hy.le, ← mul_div_assoc,
    ← Real.rpow_add hx]
  congr 2; ring

lemma Dh_dom (n : ℕ) : ∃ c : ℝ, 0 ≤ c ∧ ∀ ε y x : ℝ, 0 < y → 0 < x → ∀ θ : ℝ, 0 ≤ θ → θ ≤ 1 →
    |Dh n ε y x| ≤ c * (|Real.log x| ^ n * x ^ (θ * (2 - 2 * ε) + ε - 1) * Real.exp (-x))
      / y ^ (2 * θ) := by
  obtain ⟨c, hc, hΨ⟩ := Ψ_bd n
  refine ⟨c, hc, fun ε y x hy hx θ hθ0 hθ1 => ?_⟩
  have hu : 0 ≤ x ^ (2 - 2 * ε) / y ^ 2 := by positivity
  have h := hΨ _ hu θ hθ0 hθ1
  unfold Dh
  rw [abs_mul, abs_mul, abs_mul, abs_pow, abs_of_pos (Real.rpow_pos_of_pos hx _),
    abs_of_pos (Real.exp_pos _)]
  have e := pow_alg hx hy ε θ
  calc |Real.log x| ^ n * x ^ (ε - 1) * Real.exp (-x) * |Ψ n (x ^ (2 - 2 * ε) / y ^ 2)|
      ≤ |Real.log x| ^ n * x ^ (ε - 1) * Real.exp (-x) * (c * (x ^ (2 - 2 * ε) / y ^ 2) ^ θ) := by
        gcongr
    _ = c * (|Real.log x| ^ n * Real.exp (-x)) * (x ^ (ε - 1) * (x ^ (2 - 2 * ε) / y ^ 2) ^ θ) := by
        ring
    _ = _ := by rw [e]; ring

lemma logpow_le {x r : ℝ} (hx : 0 < x) (hr : 0 < r) (n : ℕ) :
    |Real.log x| ^ n ≤ (n.factorial : ℝ) / r ^ n * (x ^ (-r) + x ^ r) := by
  have hf : (0 : ℝ) < n.factorial := by exact_mod_cast n.factorial_pos
  have key : ∀ t : ℝ, 0 ≤ t → t ^ n ≤ (n.factorial : ℝ) / r ^ n * Real.exp (r * t) := by
    intro t ht
    have := Real.pow_div_factorial_le_exp (r * t) (by positivity) n
    rw [div_le_iff₀ hf, mul_pow] at this
    rw [div_mul_eq_mul_div, le_div_iff₀ (by positivity)]
    linarith
  rcases le_total x 1 with hx1 | hx1
  · have hl : Real.log x ≤ 0 := Real.log_nonpos hx.le hx1
    rw [abs_of_nonpos hl]
    have h := key (-Real.log x) (by linarith)
    have e : Real.exp (r * -Real.log x) = x ^ (-r) := by
      rw [Real.rpow_def_of_pos hx]; ring_nf
    rw [e] at h
    have : 0 ≤ (n.factorial : ℝ) / r ^ n * x ^ r := by positivity
    nlinarith
  · have hl : 0 ≤ Real.log x := Real.log_nonneg hx1
    rw [abs_of_nonneg hl]
    have h := key (Real.log x) hl
    have e : Real.exp (r * Real.log x) = x ^ r := by
      rw [Real.rpow_def_of_pos hx]; ring_nf
    rw [e] at h
    have : 0 ≤ (n.factorial : ℝ) / r ^ n * x ^ (-r) := by positivity
    nlinarith

lemma J_meas (n : ℕ) (s : ℝ) :
    Measurable (fun x : ℝ => |Real.log x| ^ n * x ^ (s - 1) * Real.exp (-x)) := by
  fun_prop

lemma J_int (n : ℕ) {s : ℝ} (hs : 0 < s) :
    IntegrableOn (fun x => |Real.log x| ^ n * x ^ (s - 1) * Real.exp (-x)) (Ioi 0) ∧
      ∫ x in Ioi 0, |Real.log x| ^ n * x ^ (s - 1) * Real.exp (-x) ≤
        (n.factorial : ℝ) / (s / 2) ^ n * (Real.Gamma (s / 2) + Real.Gamma (3 * s / 2)) := by
  have g1 := Real.GammaIntegral_convergent (s := s / 2) (by positivity)
  have g2 := Real.GammaIntegral_convergent (s := 3 * s / 2) (by positivity)
  have hg : IntegrableOn (fun x => (n.factorial : ℝ) / (s / 2) ^ n *
      (Real.exp (-x) * x ^ (s / 2 - 1) + Real.exp (-x) * x ^ (3 * s / 2 - 1))) (Ioi 0) :=
    (g1.add g2).const_mul ((n.factorial : ℝ) / (s / 2) ^ n)
  have hpt : ∀ x ∈ Ioi (0 : ℝ), |Real.log x| ^ n * x ^ (s - 1) * Real.exp (-x) ≤
      (n.factorial : ℝ) / (s / 2) ^ n *
        (Real.exp (-x) * x ^ (s / 2 - 1) + Real.exp (-x) * x ^ (3 * s / 2 - 1)) := by
    intro x hx
    have hx0 : (0 : ℝ) < x := hx
    have h := logpow_le hx0 (show 0 < s / 2 by positivity) n
    have e1 : x ^ (-(s / 2)) * x ^ (s - 1) = x ^ (s / 2 - 1) := by
      rw [← Real.rpow_add hx0]; ring_nf
    have e2 : x ^ (s / 2) * x ^ (s - 1) = x ^ (3 * s / 2 - 1) := by
      rw [← Real.rpow_add hx0]; ring_nf
    calc |Real.log x| ^ n * x ^ (s - 1) * Real.exp (-x)
        ≤ (n.factorial : ℝ) / (s / 2) ^ n * (x ^ (-(s / 2)) + x ^ (s / 2)) * x ^ (s - 1)
            * Real.exp (-x) := by gcongr
      _ = (n.factorial : ℝ) / (s / 2) ^ n * (Real.exp (-x) * (x ^ (-(s / 2)) * x ^ (s - 1))
            + Real.exp (-x) * (x ^ (s / 2) * x ^ (s - 1))) := by ring
      _ = _ := by rw [e1, e2]
  have hint : IntegrableOn (fun x => |Real.log x| ^ n * x ^ (s - 1) * Real.exp (-x)) (Ioi 0) := by
    refine Integrable.mono' hg (J_meas n s).aestronglyMeasurable ?_
    refine (ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun x hx => ?_)
    have hx0 : (0 : ℝ) < x := hx
    have := Real.rpow_pos_of_pos hx0 (s - 1)
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
    exact hpt x hx
  refine ⟨hint, ?_⟩
  calc _ ≤ ∫ x in Ioi 0, (n.factorial : ℝ) / (s / 2) ^ n *
        (Real.exp (-x) * x ^ (s / 2 - 1) + Real.exp (-x) * x ^ (3 * s / 2 - 1)) :=
        setIntegral_mono_on hint hg measurableSet_Ioi hpt
    _ = _ := by
        rw [integral_const_mul, integral_add g1 g2, ← Real.Gamma_eq_integral (by positivity),
          ← Real.Gamma_eq_integral (by positivity)]

lemma rpow_le_sum {x a b c : ℝ} (hx : 0 < x) (hba : b ≤ a) (hac : a ≤ c) :
    x ^ a ≤ x ^ b + x ^ c := by
  rcases le_total x 1 with hx1 | hx1
  · have := Real.rpow_le_rpow_of_exponent_ge hx hx1 hba
    have : 0 ≤ x ^ c := by positivity
    linarith
  · have := Real.rpow_le_rpow_of_exponent_le hx1 hac
    have : 0 ≤ x ^ b := by positivity
    linarith

end R04

namespace R04

/-- an open interval around `[0, 1/4]` -/
noncomputable def U0 : Set ℝ := Ioo (-1 / 20) (3 / 10)

lemma U0_open : IsOpen U0 := isOpen_Ioo

lemma I0_sub : I0 ⊆ U0 := fun x hx => ⟨by linarith [hx.1], by linarith [hx.2]⟩

/-- `∂_ε^n h_ε(y)` as an integral -/
noncomputable def hd (n : ℕ) (ε y : ℝ) : ℝ := ∫ x in Ioi (0 : ℝ), Dh n ε y x

lemma h_eq (ε y : ℝ) : Z04.h ε y = hd 0 ε y := by
  unfold Z04.h hd; rw [Dh_zero]

lemma hd_chain {y : ℝ} (hy : 0 < y) :
    ∀ n, ∀ ε ∈ U0, HasDerivAt (fun e => hd n e y) (hd (n + 1) ε y) ε := by
  refine integral_chain (μ := volume.restrict (Ioi 0)) U0_open (fun n ε x => Dh n ε y x)
    (fun ε _ => ?_) ?_ (fun n => ?_)
  · show AEStronglyMeasurable (fun x => Dh 0 ε y x) _
    rw [Dh_zero]; exact (Z04.hI_meas ε y).aestronglyMeasurable
  · exact (ae_restrict_iff' measurableSet_Ioi).mpr
      (Eventually.of_forall fun x hx n ε _ => Dh_hasDeriv n hx y ε)
  · obtain ⟨c, hc, hD⟩ := Dh_dom n
    refine ⟨fun x => c * (|Real.log x| ^ n * x ^ ((1 : ℝ) - 1) * Real.exp (-x)) /
      y ^ (2 * (1 / 2 : ℝ)), ((J_int n one_pos).1.const_mul c).div_const _, ?_⟩
    refine (ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun x hx ε _ => ?_)
    have := hD ε y x hy hx (1 / 2) (by norm_num) (by norm_num)
    rwa [show (1 / 2 : ℝ) * (2 - 2 * ε) + ε - 1 = (1 : ℝ) - 1 by ring] at this

lemma h_smooth {y : ℝ} (hy : 0 < y) :
    ContDiffOn ℝ ∞ (fun e => Z04.h e y) U0 ∧
      ∀ n, Set.EqOn (iteratedDeriv n (fun e => Z04.h e y)) (fun e => hd n e y) U0 := by
  have := chain_smooth U0_open (fun n e => hd n e y) (hd_chain hy)
  have e : (fun e => hd 0 e y) = (fun e => Z04.h e y) := by funext e; rw [h_eq]
  simp only [e] at this
  exact this

lemma hd_bd_gen (n : ℕ) : ∃ c : ℝ, 0 ≤ c ∧ ∀ ε y : ℝ, 0 < y → ∀ θ : ℝ, 0 ≤ θ → θ ≤ 1 →
    0 < θ * (2 - 2 * ε) + ε →
    |hd n ε y| ≤ c / y ^ (2 * θ) *
      ∫ x in Ioi 0, |Real.log x| ^ n * x ^ (θ * (2 - 2 * ε) + ε - 1) * Real.exp (-x) := by
  obtain ⟨c, hc, hD⟩ := Dh_dom n
  refine ⟨c, hc, fun ε y hy θ hθ0 hθ1 ha => ?_⟩
  have hJ := (J_int n ha).1
  have hJ' : ∀ x, |Real.log x| ^ n * x ^ (θ * (2 - 2 * ε) + ε - 1) * Real.exp (-x) =
      |Real.log x| ^ n * x ^ ((θ * (2 - 2 * ε) + ε) - 1) * Real.exp (-x) := fun x => rfl
  have := norm_integral_le_of_norm_le (μ := volume.restrict (Ioi 0)) (f := fun x => Dh n ε y x)
    ((hJ.const_mul c).div_const (y ^ (2 * θ)))
    ((ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun x hx => by
      rw [Real.norm_eq_abs]; exact hD ε y x hy hx θ hθ0 hθ1))
  rw [Real.norm_eq_abs, integral_div, integral_const_mul] at this
  unfold hd
  calc _ ≤ _ := this
    _ = _ := by ring

lemma J_mono (n : ℕ) {a b c : ℝ} (hb : 0 < b) (hba : b ≤ a) (hac : a ≤ c) :
    ∫ x in Ioi 0, |Real.log x| ^ n * x ^ (a - 1) * Real.exp (-x) ≤
      (∫ x in Ioi 0, |Real.log x| ^ n * x ^ (b - 1) * Real.exp (-x)) +
        ∫ x in Ioi 0, |Real.log x| ^ n * x ^ (c - 1) * Real.exp (-x) := by
  have i1 := (J_int n hb).1
  have i2 := (J_int n (hb.trans_le (hba.trans hac))).1
  rw [← integral_add i1 i2]
  refine setIntegral_mono_on (J_int n (hb.trans_le hba)).1 (i1.add i2) measurableSet_Ioi
    fun x hx => ?_
  have hx0 : (0 : ℝ) < x := hx
  have := rpow_le_sum hx0 (show b - 1 ≤ a - 1 by linarith) (show a - 1 ≤ c - 1 by linarith)
  have h0 : 0 ≤ |Real.log x| ^ n * Real.exp (-x) := by positivity
  calc |Real.log x| ^ n * x ^ (a - 1) * Real.exp (-x)
      = (|Real.log x| ^ n * Real.exp (-x)) * x ^ (a - 1) := by ring
    _ ≤ (|Real.log x| ^ n * Real.exp (-x)) * (x ^ (b - 1) + x ^ (c - 1)) := by gcongr
    _ = _ := by ring

/-- the integral `∫ |log x|^n x^{s-1} e^{-x}` -/
noncomputable def J (n : ℕ) (s : ℝ) : ℝ :=
  ∫ x in Ioi 0, |Real.log x| ^ n * x ^ (s - 1) * Real.exp (-x)

lemma J_nn (n : ℕ) (s : ℝ) : 0 ≤ J n s :=
  setIntegral_nonneg measurableSet_Ioi fun x hx => by
    have := Real.rpow_pos_of_pos (show (0 : ℝ) < x from hx) (s - 1); positivity

lemma hd_U (n : ℕ) : ∃ A : ℝ, 0 ≤ A ∧ ∀ ε ∈ U0, ∀ y : ℝ, 0 < y →
    |hd n ε y| ≤ A / y ^ (1 / 5 : ℝ) ∧ |hd n ε y| ≤ A / y ^ 2 := by
  obtain ⟨c, hc, hb⟩ := hd_bd_gen n
  refine ⟨c * (J n (4 / 25) + J n 1 + J n 1 + J n 3), by
    have := J_nn n (4 / 25); have := J_nn n 1; have := J_nn n 3; positivity,
    fun ε hε y hy => ⟨?_, ?_⟩⟩
  · have h1 := hb ε y hy (1 / 10) (by norm_num) (by norm_num) (by linarith [hε.1])
    have h2 := J_mono n (a := 1 / 10 * (2 - 2 * ε) + ε) (b := 4 / 25) (c := 1) (by norm_num)
      (by linarith [hε.1]) (by linarith [hε.2])
    rw [show (2 : ℝ) * (1 / 10) = 1 / 5 by norm_num] at h1
    have hJ3 := J_nn n 3
    have hy' : 0 < y ^ (1 / 5 : ℝ) := Real.rpow_pos_of_pos hy _
    calc _ ≤ _ := h1
      _ ≤ c / y ^ (1 / 5 : ℝ) * (J n (4 / 25) + J n 1) := by
          unfold J; gcongr
      _ ≤ c / y ^ (1 / 5 : ℝ) * (J n (4 / 25) + J n 1 + J n 1 + J n 3) := by
          apply mul_le_mul_of_nonneg_left _ (by positivity); linarith [J_nn n 1, J_nn n 3]
      _ = _ := by ring
  · have h1 := hb ε y hy 1 (by norm_num) (by norm_num) (by linarith [hε.2])
    have h2 := J_mono n (a := 1 * (2 - 2 * ε) + ε) (b := 1) (c := 3) (by norm_num)
      (by linarith [hε.2]) (by linarith [hε.1])
    rw [show (2 : ℝ) * 1 = 2 by norm_num, Real.rpow_two] at h1
    have hJ := J_nn n (4 / 25)
    calc _ ≤ _ := h1
      _ ≤ c / y ^ 2 * (J n 1 + J n 3) := by
          unfold J; gcongr
      _ ≤ c / y ^ 2 * (J n (4 / 25) + J n 1 + J n 1 + J n 3) := by
          apply mul_le_mul_of_nonneg_left _ (by positivity); linarith [J_nn n 1, J_nn n (4 / 25)]
      _ = _ := by ring

lemma Gamma_le {s : ℝ} (hs : 0 < s) (hs3 : s ≤ 3) : Real.Gamma s ≤ 2 / s + 2 := by
  have h13 : ∀ t ∈ Icc (1 : ℝ) 3, Real.Gamma t ≤ 2 := by
    intro t ht
    have := Real.convexOn_Gamma.le_max_of_mem_Icc (show (1 : ℝ) ∈ Ioi 0 by norm_num)
      (show (3 : ℝ) ∈ Ioi 0 by norm_num) ht
    have g3 : Real.Gamma 3 = 2 := by
      rw [show (3 : ℝ) = (2 : ℕ) + 1 by norm_num, Real.Gamma_nat_eq_factorial]; norm_num
    rw [Real.Gamma_one, g3] at this
    simpa using this
  rcases le_total 1 s with h1 | h1
  · have := h13 s ⟨h1, hs3⟩
    have : 0 ≤ 2 / s := by positivity
    linarith
  · have e := Real.Gamma_add_one hs.ne'
    have := h13 (s + 1) ⟨by linarith, by linarith⟩
    rw [e] at this
    have : Real.Gamma s ≤ 2 / s := by
      rw [le_div_iff₀ hs]; linarith
    linarith

lemma J_le (n : ℕ) {a : ℝ} (ha : 0 < a) (ha2 : a ≤ 2) :
    J n a ≤ 16 * (n.factorial : ℝ) * 2 ^ n / a ^ (n + 1) := by
  have h := (J_int n ha).2
  have g1 := Gamma_le (s := a / 2) (by positivity) (by linarith)
  have g2 := Gamma_le (s := 3 * a / 2) (by positivity) (by linarith)
  have e1 : 2 / (a / 2) + 2 ≤ 8 / a := by
    rw [div_add' _ _ _ (by positivity), div_le_div_iff₀ (by positivity) ha]; nlinarith
  have e2 : 2 / (3 * a / 2) + 2 ≤ 8 / a := by
    rw [div_add' _ _ _ (by positivity), div_le_div_iff₀ (by positivity) ha]; nlinarith
  have hf : (0 : ℝ) ≤ n.factorial / (a / 2) ^ n := by positivity
  unfold J
  calc _ ≤ _ := h
    _ ≤ (n.factorial : ℝ) / (a / 2) ^ n * (8 / a + 8 / a) := by gcongr <;> linarith
    _ = _ := by
        rw [div_pow]
        field_simp
        ring

lemma hd_log (n : ℕ) : ∃ C : ℝ, 0 ≤ C ∧ ∀ ε ∈ I0, ∀ y : ℝ, 0 < y → y ≤ 1 →
    |hd n ε y| ≤ C * (1 + |Real.log y|) ^ (n + 1) := by
  obtain ⟨c, hc, hb⟩ := hd_bd_gen n
  refine ⟨c * Real.exp 2 * (16 * (n.factorial : ℝ) * 2 ^ n), by positivity,
    fun ε hε y hy hy1 => ?_⟩
  set L := |Real.log y| with hL
  have hL0 : 0 ≤ L := abs_nonneg _
  set θ := 1 / (1 + L) with hθ
  have hθ0 : 0 < θ := by positivity
  have hθ1 : θ ≤ 1 := by rw [hθ, div_le_one (by positivity)]; linarith
  have he0 := hε.1
  have he1 := hε.2
  set a := θ * (2 - 2 * ε) + ε with ha
  have ha0 : θ ≤ a := by nlinarith
  have ha2 : a ≤ 2 := by nlinarith
  have hapos : 0 < a := by linarith
  have h1 := hb ε y hy θ hθ0.le hθ1 hapos
  have h2 := J_le n hapos ha2
  -- 1 / y^(2θ) ≤ e²
  have hy2 : 1 / y ^ (2 * θ) ≤ Real.exp 2 := by
    have hl : Real.log y ≤ 0 := Real.log_nonpos hy.le hy1
    have hLy : L = -Real.log y := by rw [hL, abs_of_nonpos hl]
    have hk : -2 ≤ Real.log y * (2 * θ) := by
      rw [hθ, hLy]
      rw [show Real.log y * (2 * (1 / (1 + -Real.log y))) =
        -2 * (-Real.log y / (1 + -Real.log y)) by field_simp]
      have : -Real.log y / (1 + -Real.log y) ≤ 1 := by
        rw [div_le_one (by linarith)]; linarith
      linarith
    rw [Real.rpow_def_of_pos hy, one_div, ← Real.exp_neg]
    exact Real.exp_le_exp.mpr (by linarith)
  have hθa : 1 / a ^ (n + 1) ≤ (1 + L) ^ (n + 1) := by
    have : 1 / a ≤ 1 + L := by
      rw [div_le_iff₀ hapos]
      have : 1 = θ * (1 + L) := by rw [hθ]; field_simp
      nlinarith
    calc 1 / a ^ (n + 1) = (1 / a) ^ (n + 1) := by rw [one_div_pow]
      _ ≤ _ := by gcongr
  have hJ : J n a ≤ 16 * (n.factorial : ℝ) * 2 ^ n * (1 + L) ^ (n + 1) := by
    calc J n a ≤ 16 * (n.factorial : ℝ) * 2 ^ n / a ^ (n + 1) := h2
      _ = 16 * (n.factorial : ℝ) * 2 ^ n * (1 / a ^ (n + 1)) := by ring
      _ ≤ _ := by gcongr
  have hJn := J_nn n a
  calc |hd n ε y| ≤ c / y ^ (2 * θ) * J n a := h1
    _ = c * (1 / y ^ (2 * θ)) * J n a := by ring
    _ ≤ c * Real.exp 2 * (16 * (n.factorial : ℝ) * 2 ^ n * (1 + L) ^ (n + 1)) := by
        gcongr
    _ = _ := by ring

end R04

theorem R04.proof_5 : R04.stmt_5 := by
  refine ⟨fun y hy => (R04.h_smooth hy).1.mono R04.I0_sub, fun j => ?_⟩
  obtain ⟨C1, hC1, h1⟩ := R04.hd_log j
  obtain ⟨A, hA, h2⟩ := R04.hd_U j
  refine ⟨C1 + A, fun y hy ε hε => ?_⟩
  have e : iteratedDerivWithin j (fun e => Z04.h e y) R04.I0 ε = R04.hd j ε y := by
    rw [iteratedDerivWithin_eq_iteratedDeriv
      (show UniqueDiffOn ℝ R04.I0 from uniqueDiffOn_Icc (by norm_num))
      (((R04.h_smooth hy).1.contDiffAt (R04.U0_open.mem_nhds (R04.I0_sub hε))).of_le
        (by exact_mod_cast le_top)) hε]
    exact (R04.h_smooth hy).2 j (R04.I0_sub hε)
  rw [e]
  refine ⟨fun hy1 => ?_, fun hy1 => ?_⟩
  · have := h1 ε hε y hy hy1
    have : 0 ≤ A * (1 + |Real.log y|) ^ (j + 1) := by positivity
    nlinarith
  · have := (h2 ε (R04.I0_sub hε) y hy).2
    have : 0 ≤ C1 / y ^ 2 := by positivity
    calc _ ≤ A / y ^ 2 := by assumption
      _ ≤ (C1 + A) / y ^ 2 := by gcongr; linarith

namespace R04

lemma iter_hasDeriv {g : ℝ → ℝ} {U : Set ℝ} (hU : IsOpen U) (hg : ContDiffOn ℝ ∞ g U) {ε : ℝ}
    (hε : ε ∈ U) (n : ℕ) : HasDerivAt (iteratedDeriv n g) (iteratedDeriv (n + 1) g ε) ε := by
  have hd : DifferentiableOn ℝ (iteratedDerivWithin n g U) U :=
    hg.differentiableOn_iteratedDerivWithin (by exact_mod_cast ENat.natCast_lt_top n) hU.uniqueDiffOn
  have heq : iteratedDeriv n g =ᶠ[𝓝 ε] iteratedDerivWithin n g U :=
    Filter.eventuallyEq_of_mem (hU.mem_nhds hε) (fun e he => (iteratedDerivWithin_of_isOpen hU he).symm)
  have hda : DifferentiableAt ℝ (iteratedDeriv n g) ε :=
    ((hd ε hε).differentiableAt (hU.mem_nhds hε)).congr_of_eventuallyEq heq
  rw [iteratedDeriv_succ]; exact hda.hasDerivAt

lemma pow_bd {f : ℝ → ℝ} {U : Set ℝ} (hU : IsOpen U) (hf : ContDiffOn ℝ ∞ f U) {ε : ℝ}
    (hε : ε ∈ U) {n : ℕ} {A : ℝ} (hA : ∀ i ≤ n, |iteratedDerivWithin i f U ε| ≤ A) :
    ∀ m : ℕ, ∀ i ≤ n, |iteratedDerivWithin i (fun e => f e ^ (m + 1)) U ε| ≤
      (2 ^ n) ^ m * A ^ (m + 1) := by
  have hA0 : 0 ≤ A := (abs_nonneg _).trans (hA 0 (Nat.zero_le _))
  intro m
  induction m with
  | zero => intro i hi; simpa using hA i hi
  | succ m ih =>
    intro i hi
    have e : (fun e => f e ^ (m + 1 + 1)) = fun e => f e ^ (m + 1) * f e := by
      funext e; ring
    rw [e]
    have hpm : ContDiffOn ℝ ∞ (fun e => f e ^ (m + 1)) U := hf.pow _
    have := norm_iteratedFDerivWithin_mul_le hpm hf hU.uniqueDiffOn hε (n := i)
      (by exact_mod_cast le_top)
    simp only [norm_iteratedFDerivWithin_eq_norm_iteratedDerivWithin, Real.norm_eq_abs] at this
    have hs : ∑ j ∈ Finset.range (i + 1), (i.choose j : ℝ) = 2 ^ i := by
      exact_mod_cast Nat.sum_range_choose i
    calc _ ≤ _ := this
      _ ≤ ∑ j ∈ Finset.range (i + 1), (i.choose j : ℝ) * ((2 ^ n) ^ m * A ^ (m + 1)) * A :=
          Finset.sum_le_sum fun j hj => by
            have hj' : j ≤ i := Nat.lt_succ_iff.mp (Finset.mem_range.mp hj)
            have a1 := ih j (hj'.trans hi)
            have a2 := hA (i - j) ((Nat.sub_le i j).trans hi)
            gcongr
      _ = 2 ^ i * ((2 ^ n) ^ m * A ^ (m + 1)) * A := by
          rw [← Finset.sum_mul, ← Finset.sum_mul, hs]
      _ ≤ 2 ^ n * ((2 ^ n) ^ m * A ^ (m + 1)) * A := by
          gcongr
          · norm_num
      _ = _ := by ring

/-- the weight `min(y^{-1/5}, y^{-2})` -/
noncomputable def wt (y : ℝ) : ℝ := min (1 / y ^ (1 / 5 : ℝ)) (1 / y ^ 2)

lemma wt_nn {y : ℝ} (hy : 0 < y) : 0 ≤ wt y := by
  unfold wt; have := Real.rpow_pos_of_pos hy (1 / 5); positivity

lemma wt_int {k : ℕ} (hk1 : 1 ≤ k) (hk : k ≤ 3) : IntegrableOn (fun y => wt y ^ k) (Ioi 0) := by
  have hmeas : Measurable (fun y => wt y ^ k) := by unfold wt; fun_prop
  rw [← Ioc_union_Ioi_eq_Ioi (zero_le_one' ℝ)]
  refine IntegrableOn.union ?_ ?_
  · have hr : IntegrableOn (fun y : ℝ => y ^ (-(k : ℝ) / 5)) (Ioc 0 1) := by
      have := intervalIntegral.intervalIntegrable_rpow' (a := 0) (b := 1)
        (r := -(k : ℝ) / 5) (by have : (k : ℝ) ≤ 3 := (by exact_mod_cast hk); linarith)
      exact (intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one).mp this
    refine Integrable.mono' hr hmeas.aestronglyMeasurable ?_
    refine (ae_restrict_iff' measurableSet_Ioc).mpr (Eventually.of_forall fun y hy => ?_)
    have hy0 : 0 < y := hy.1
    rw [Real.norm_eq_abs, abs_of_nonneg (pow_nonneg (wt_nn hy0) _)]
    calc wt y ^ k ≤ (1 / y ^ (1 / 5 : ℝ)) ^ k := by
          gcongr
          · exact wt_nn hy0
          · exact min_le_left _ _
      _ = _ := by
          rw [one_div, ← Real.rpow_neg hy0.le, ← Real.rpow_natCast, ← Real.rpow_mul hy0.le]
          ring_nf
  · have hr : IntegrableOn (fun y : ℝ => y ^ (-(2 * (k : ℝ)))) (Ioi 1) :=
      integrableOn_Ioi_rpow_of_lt (by have : (1 : ℝ) ≤ k := (by exact_mod_cast hk1); linarith)
        one_pos
    refine Integrable.mono' hr hmeas.aestronglyMeasurable ?_
    refine (ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun y hy => ?_)
    have hy0 : (0 : ℝ) < y := lt_trans one_pos hy
    rw [Real.norm_eq_abs, abs_of_nonneg (pow_nonneg (wt_nn hy0) _)]
    calc wt y ^ k ≤ (1 / y ^ 2) ^ k := by
          gcongr
          · exact wt_nn hy0
          · exact min_le_right _ _
      _ = _ := by
          rw [one_div, ← Real.rpow_natCast y 2, ← Real.rpow_neg hy0.le, ← Real.rpow_natCast,
            ← Real.rpow_mul hy0.le]
          push_cast; ring_nf

lemma Hm_facts (m : ℕ) (hm : m ≤ 2) :
    (∀ ε ∈ U0, IntegrableOn (fun y => Z04.h ε y ^ (m + 1)) (Ioi 0)) ∧
      ContDiffOn ℝ ∞ (fun ε => ∫ y in Ioi (0 : ℝ), Z04.h ε y ^ (m + 1)) U0 := by
  set E : ℕ → ℝ → ℝ → ℝ := fun n ε y => iteratedDeriv n (fun e => Z04.h e y ^ (m + 1)) ε with hE
  have hm0 : ∀ ε ∈ U0, AEStronglyMeasurable (E 0 ε) (volume.restrict (Ioi 0)) := by
    intro ε hε
    have : E 0 ε = fun y => Z04.h ε y ^ (m + 1) := by funext y; simp [hE]
    rw [this]
    exact (Z04.h_int (ε := ε) (by linarith [hε.2])).aestronglyMeasurable.pow _
  have hdv : ∀ᵐ y ∂(volume.restrict (Ioi (0 : ℝ))), ∀ n, ∀ ε ∈ U0,
      HasDerivAt (fun e => E n e y) (E (n + 1) ε y) ε := by
    refine (ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun y hy n ε hε => ?_)
    exact iter_hasDeriv U0_open ((h_smooth hy).1.pow _) hε n
  choose A hA0 hAb using hd_U
  have hb : ∀ n, ∃ B : ℝ → ℝ, Integrable B (volume.restrict (Ioi 0)) ∧
      ∀ᵐ y ∂(volume.restrict (Ioi (0 : ℝ))), ∀ ε ∈ U0, |E n ε y| ≤ B y := by
    intro n
    set A' := ∑ i ∈ Finset.range (n + 1), A i with hA'
    have hA'0 : 0 ≤ A' := Finset.sum_nonneg fun i _ => hA0 i
    refine ⟨fun y => (2 ^ n) ^ m * A' ^ (m + 1) * wt y ^ (m + 1),
      (wt_int (by omega) (by omega)).const_mul _, ?_⟩
    refine (ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun y hy ε hε => ?_)
    have hy0 : (0 : ℝ) < y := hy
    have hsm := (h_smooth hy0)
    have key : ∀ i ≤ n, |iteratedDerivWithin i (fun e => Z04.h e y) U0 ε| ≤ A' * wt y := by
      intro i hi
      rw [iteratedDerivWithin_of_isOpen U0_open hε, hsm.2 i hε]
      obtain ⟨b1, b2⟩ := hAb i ε hε y hy0
      have hAi : A i ≤ A' := Finset.single_le_sum (fun j _ => hA0 j)
        (Finset.mem_range.mpr (by omega))
      calc |hd i ε y| ≤ A i * wt y := by
            unfold wt
            rw [mul_min_of_nonneg _ _ (hA0 i)]
            refine le_min ?_ ?_
            · exact b1.trans (le_of_eq (by ring))
            · exact b2.trans (le_of_eq (by ring))
        _ ≤ A' * wt y := by gcongr; exact wt_nn hy0
    have := pow_bd U0_open hsm.1 hε key m n le_rfl
    rw [iteratedDerivWithin_of_isOpen U0_open hε] at this
    calc |E n ε y| ≤ (2 ^ n) ^ m * (A' * wt y) ^ (m + 1) := this
      _ = _ := by ring
  have hch := integral_chain U0_open E hm0 hdv hb
  have hcs := chain_smooth U0_open (fun n ε => ∫ y, E n ε y ∂(volume.restrict (Ioi 0))) hch
  have e0 : (fun ε => ∫ y, E 0 ε y ∂(volume.restrict (Ioi 0))) =
      fun ε => ∫ y in Ioi (0 : ℝ), Z04.h ε y ^ (m + 1) := by
    funext ε; simp [hE]
  refine ⟨fun ε hε => ?_, by rw [← e0]; exact hcs.1⟩
  obtain ⟨B, hBi, hB⟩ := hb 0
  refine Integrable.mono' hBi ?_ ?_
  · have := hm0 ε hε
    simpa [hE] using this
  · filter_upwards [hB] with y hy
    have := hy ε hε
    simpa [hE, Real.norm_eq_abs] using this

lemma Q_cd : ContDiffOn ℝ ∞ Z04.Q U0 := by
  have hI := (Hm_facts 2 le_rfl).2
  have hG : ContDiffOn ℝ ∞ (fun e => (Real.Gamma (1 + 2 * e))⁻¹) U0 := fun e he =>
    ((Gamma_comp_cd (show 0 < 1 + 2 * e by linarith [he.1])).inv
      (Real.Gamma_pos_of_pos (show 0 < 1 + 2 * e by linarith [he.1])).ne').contDiffWithinAt
  refine ((contDiffOn_const (c := -(1 / Real.sqrt π))).mul (hG.mul hI)).congr (fun e _ => ?_)
  unfold Z04.Q
  rw [show (2 + 1 : ℕ) = 3 from rfl]
  simp only [one_div, mul_inv]
  ring

end R04

theorem R04.proof_3 : R04.stmt_3 :=
  (R04.Q_cd.of_le (by exact_mod_cast le_top)).mono R04.I0_sub

theorem R04.proof_4 : R04.stmt_4 :=
  (R04.Q_cd.of_le (WithTop.coe_le_coe.mpr le_top)).mono R04.I0_sub

theorem R04.proof_6 : R04.stmt_6 := by
  intro m hm1 hm3
  obtain ⟨k, rfl⟩ : ∃ k, m = k + 1 := ⟨m - 1, by omega⟩
  have hf := R04.Hm_facts k (by omega)
  refine ⟨fun ε hε => hf.1 ε (R04.I0_sub hε), ?_⟩
  have : Z04.H (k + 1) = fun ε => 1 / Real.sqrt π * ∫ y in Ioi (0 : ℝ), Z04.h ε y ^ (k + 1) := by
    funext ε; rfl
  rw [this]
  exact (contDiffOn_const.mul hf.2).mono R04.I0_sub

theorem R04.proof_7 : R04.stmt_7 := R04.Q_cd.mono R04.I0_sub

namespace R04

lemma P_series {p : ℝ} (hp : 1 < p) :
    HasSum (fun k => M k / p ^ (k + 1)) (W02.Pint p / 8 - 1) := by
  have h := Complex.reCLM.hasSum (W02.proof_6 p hp)
  convert h using 1
  · funext k
    simp only [Complex.reCLM_apply]
    rw [← Complex.ofReal_pow, Complex.div_ofReal_re]
    rfl
  · simp

/-- the coefficients of `P(1/ε)/8 - 1 = Σ M_k ε^{k+1}` -/
noncomputable def mc (k : ℕ) : ℝ := if k = 0 then 0 else M (k - 1)

lemma AE_f : AE (fun ε => W02.Pint (1 / ε) / 8 - 1) mc := by
  have ht : Tendsto (fun k => |M k / 2 ^ (k + 1)|) atTop (𝓝 0) := by
    simpa using (P_series one_lt_two).summable.tendsto_atTop_zero.abs
  obtain ⟨B, hB⟩ := ht.bddAbove_range
  have hB' : ∀ k, |M k / 2 ^ (k + 1)| ≤ B := fun k => hB ⟨k, rfl⟩
  have hB0 : 0 ≤ B := (abs_nonneg _).trans (hB' 0)
  intro K
  refine ⟨B * 2 ^ (K + 1) * 2, 1 / 4, by norm_num, fun ε hε hε4 => ?_⟩
  have hp : 1 < 1 / ε := by rw [lt_div_iff₀ hε]; linarith
  have hs := P_series hp
  set t : ℕ → ℝ := fun k => M k / (1 / ε) ^ (k + 1) with htdef
  have ht' : ∀ k, t k = M k * ε ^ (k + 1) := fun k => by
    simp only [htdef]; rw [one_div_pow, div_div_eq_mul_div, div_one]
  have hT := (hasSum_nat_add_iff' K).mpr hs
  have esum : ∑ k ∈ Finset.range (K + 1), mc k * ε ^ k = ∑ i ∈ Finset.range K, t i := by
    rw [Finset.sum_range_succ']
    have hm0 : mc 0 = 0 := by simp [mc]
    rw [hm0, zero_mul, add_zero]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [ht']; simp [mc]
  show |W02.Pint (1 / ε) / 8 - 1 - ∑ k ∈ Finset.range (K + 1), mc k * ε ^ k| ≤ _
  rw [esum, ← hT.tsum_eq, ← Real.norm_eq_abs]
  refine (tsum_of_norm_bounded ((hasSum_geometric_two).mul_left (B * 2 ^ (K + 1) * ε ^ (K + 1)))
    fun n => ?_).trans (le_of_eq (by ring))
  rw [Real.norm_eq_abs, ht']
  have h1 := hB' (n + K)
  have e : M (n + K) * ε ^ (n + K + 1) =
      (M (n + K) / 2 ^ (n + K + 1)) * ((2 * ε) ^ (K + 1) * (2 * ε) ^ n) := by
    rw [← pow_add, show K + 1 + n = n + K + 1 by ring, mul_pow]; field_simp
  rw [e, abs_mul, abs_of_nonneg (by positivity : (0 : ℝ) ≤ (2 * ε) ^ (K + 1) * (2 * ε) ^ n)]
  have h2 : (2 * ε) ^ n ≤ (1 / 2) ^ n := pow_le_pow_left₀ (by positivity) (by linarith) n
  calc |M (n + K) / 2 ^ (n + K + 1)| * ((2 * ε) ^ (K + 1) * (2 * ε) ^ n)
      ≤ B * ((2 * ε) ^ (K + 1) * (1 / 2) ^ n) := by gcongr
    _ = _ := by rw [mul_pow]; ring

/-- the Taylor coefficients of `Q` at `0` (two-sided) -/
noncomputable def qd (k : ℕ) : ℝ := iteratedDeriv k Z04.Q 0 / (k.factorial : ℝ)

lemma Qc_eq (k : ℕ) : Qc k = qd k := by
  unfold Qc qd
  rw [iteratedDerivWithin_eq_iteratedDeriv
    (show UniqueDiffOn ℝ I0 from uniqueDiffOn_Icc (by norm_num))
    ((Q_cd.contDiffAt (U0_open.mem_nhds (I0_sub ⟨le_rfl, by norm_num⟩))).of_le
      (by exact_mod_cast le_top)) ⟨le_rfl, by norm_num⟩]

lemma AE_taylor {f : ℝ → ℝ} {U : Set ℝ} (hU : IsOpen U) {δ : ℝ} (hδ : 0 < δ)
    (hsub : Icc 0 δ ⊆ U) (hf : ContDiffOn ℝ ∞ f U) :
    AE f (fun k => iteratedDeriv k f 0 / (k.factorial : ℝ)) := by
  intro K
  obtain ⟨C, hC⟩ := taylor_bd hU hδ hsub hf K
  exact ⟨C, δ, hδ, fun ε hε hεδ => hC ε ⟨hε.le, hεδ⟩⟩

lemma AE_G : AE Z04.G g :=
  AE_taylor (U := Ioo (-1 / 2) 1) isOpen_Ioo (δ := 1 / 4) (by norm_num)
    (fun x hx => ⟨by linarith [hx.1], by linarith [hx.2]⟩) proof_2.1

lemma AE_Q : AE Z04.Q qd :=
  AE_taylor U0_open (δ := 1 / 4) (by norm_num) I0_sub Q_cd

/-- `Φ(ε) = 24G + 48 G (P(1/ε)/8 - 1) - 8 ε² Q`, equal to `S(1/ε)` -/
noncomputable def Φ (ε : ℝ) : ℝ :=
  24 * Z04.G ε + 48 * (Z04.G ε * (W02.Pint (1 / ε) / 8 - 1)) + -8 * (ε ^ 2 * Z04.Q ε)

lemma AE_Φ : AE Φ s := by
  have h := ((AE_G.const_mul 24).add ((AE_G.mul AE_f).const_mul 48)).add
    ((AE_sq.mul AE_Q).const_mul (-8))
  have e : ∀ k, 24 * g k + 48 * (∑ i ∈ Finset.range (k + 1), g i * mc (k - i)) +
      -8 * (∑ i ∈ Finset.range (k + 1), (if i = 2 then (1 : ℝ) else 0) * qd (k - i)) = s k := by
    intro k
    have e1 : ∑ i ∈ Finset.range (k + 1), g i * mc (k - i) =
        ∑ i ∈ Finset.range k, g i * M (k - 1 - i) := by
      rw [Finset.sum_range_succ, Nat.sub_self]
      have hm0 : mc 0 = 0 := by simp [mc]
      rw [hm0, mul_zero, add_zero]
      refine Finset.sum_congr rfl fun i hi => ?_
      have hik : i < k := Finset.mem_range.mp hi
      have : mc (k - i) = M (k - 1 - i) := by
        unfold mc; rw [if_neg (by omega), show k - i - 1 = k - 1 - i by omega]
      rw [this]
    have e2 : ∑ i ∈ Finset.range (k + 1), (if i = 2 then (1 : ℝ) else 0) * qd (k - i) =
        if 2 ≤ k then Qc (k - 2) else 0 := by
      simp only [ite_mul, one_mul, zero_mul]
      rw [Finset.sum_ite_eq' (Finset.range (k + 1)) 2 (fun i => qd (k - i))]
      simp only [Finset.mem_range, Qc_eq]
      by_cases hk : 2 ≤ k
      · rw [if_pos (by omega), if_pos hk]
      · rw [if_neg (by omega), if_neg hk]
    rw [e1, e2, s]; ring
  intro K
  obtain ⟨C, δ, hδ, hb⟩ := h K
  refine ⟨C, δ, hδ, fun ε hε hεδ => ?_⟩
  have := hb ε hε hεδ
  simp only [e] at this
  exact this

lemma S_eq_Φ {p : ℝ} (hp : 1 ≤ p) : Y05.S p = Φ (1 / p) := by
  rw [(Z04.proof_6 p hp).2.2.2.2, Φ, one_div_one_div]
  ring

end R04

theorem R04.proof_8 : R04.stmt_8 := by
  intro K
  obtain ⟨C, δ, hδ, hb⟩ := R04.AE_Φ K
  refine ⟨C, max 1 (1 / δ), fun p hp => ?_⟩
  have hp1 : 1 ≤ p := le_trans (le_max_left _ _) hp
  have hp0 : 0 < p := by linarith
  have hpδ : 1 / p ≤ δ := by
    rw [div_le_iff₀ hp0]
    have := le_trans (le_max_right _ _) hp
    rw [div_le_iff₀ hδ] at this
    linarith
  have := hb (1 / p) (by positivity) hpδ
  rw [← R04.S_eq_Φ hp1] at this
  have e : ∑ k ∈ Finset.range (K + 1), R04.s k * (1 / p) ^ k =
      ∑ k ∈ Finset.range (K + 1), R04.s k / p ^ k :=
    Finset.sum_congr rfl fun k _ => by rw [one_div_pow, mul_one_div]
  rw [e, one_div_pow, mul_one_div] at this
  exact this

theorem R04.proof_9 : R04.stmt_9 := by
  obtain ⟨_, g0, g1, g2, g3⟩ := R04.proof_2
  have q0 : R04.Qc 0 = Z04.Q 0 := by simp [R04.Qc]
  have q1 : R04.Qc 1 = derivWithin Z04.Q R04.I0 0 := by
    simp [R04.Qc, iteratedDerivWithin_one]
  refine ⟨?_, ?_, ?_, ?_⟩
  · simp [R04.s, g0]
  · simp [R04.s, g0, g1]
  · simp only [R04.s, Finset.sum_range_succ, Finset.sum_range_zero, g0, g1, g2, if_pos (le_refl 2),
      Nat.sub_self, q0]
    norm_num
    ring
  · simp only [R04.s, Finset.sum_range_succ, Finset.sum_range_zero, g0, g1, g2, g3,
      if_pos (show 2 ≤ 3 by norm_num), q1]
    norm_num
    ring

theorem R04.proof_10 : R04.stmt_10 := by
  have hd : DifferentiableAt ℝ Z04.Q 0 :=
    (R04.Q_cd.contDiffAt (R04.U0_open.mem_nhds (R04.I0_sub ⟨le_rfl, by norm_num⟩))).differentiableAt
      (by simp)
  have e : derivWithin Z04.Q R04.I0 0 = deriv Z04.Q 0 :=
    hd.derivWithin ((uniqueDiffOn_Icc (by norm_num : (0 : ℝ) < 1 / 4)) 0 ⟨le_rfl, by norm_num⟩)
  obtain ⟨h0, h1, h2, h3⟩ := R04.proof_9
  refine ⟨R04.proof_7, deriv Z04.Q 0, hd.hasDerivAt.hasDerivWithinAt, R04.s, R04.proof_8,
    h0, h1, h2, ?_⟩
  rw [h3, e]
