import Mathlib
import Zeta3.Wave3.W01_GraphLength
/-!
# Y03 -- the graph of `x u + f(x)^n v` at every angle

Target id: Y03.
Paper label: `thm:fnangle` (paper/ieee/main.tex, lines 1849-1855), with the setting and the
definitions of `u, v, θ, h_c, S_θ, k_1(θ)` in main.tex lines 1833-1847.
Detailed proof: paper/ieee/proofs/sec_g7_fnangle.tex (setting lines 5-10, theorem lines 13-19,
Lemma `g7:lem:phi` lines 58-72, Lemma `g7:lem:law` lines 74-92, Lemma `g7:lem:layer`
lines 96-107, proof of the theorem lines 109-131).
The orthogonal case (`θ = π/2`, `thm:fn`) is `Zeta3.Wave3.W01_GraphLength`; its definitions
`W01.Hyp`, `W01.D`, `W01.yd`, `W01.levelMeasure`, `W01.L`, `W01.φ` are reused here, not redefined.

Quoted verbatim (main.tex 1833):

  Take the curve $x\,u+f(x)^n\,v$ of Section~\ref{sec:corners}-A, where $u,v$ are unit vectors
  with $u\cdot v=c=\cos\theta$, $0<\theta<\pi$, and $f$ is as in Theorem~\ref{thm:fn}. With
  $y=f^n$ and $h_c$ as in Section~\ref{sec:corners}-A, the length lost is
  $2-L_n=\int_0^1h_c(y')\,dx$.

Quoted verbatim (main.tex 1834-1837):

  \[
  2-L_n=\int_0^\infty h_c'(s)\,m_y(s)\,ds=(1-\cos\theta)\,\mathbb E\bigl[m_y(S_\theta)\bigr].
  \]
  Here $h_c'(s)=1-(s+c)/\sqrt{1+2cs+s^2}$ is a density of total mass $1-\cos\theta$, and
  $S_\theta$ has this density divided by its mass:
  $\mathbb P(S_\theta>s)=\bigl(\sqrt{1+2cs+s^2}-s-c\bigr)/(1-c)$.

Quoted verbatim (main.tex 1841-1843, the value of `k_1`):

  ... the cumulants of $\log(1/S_\theta)$ are $k_j(\theta)=(-1)^jj!\,a_j$. ...
  k_1&=1+\cot^2\tfrac\theta2\,\log\cos^2\tfrac\theta2,

Quoted verbatim (main.tex 1849-1857):

  \begin{theorem}[The graph of $x\,u+f(x)^n\,v$ at every angle$^\dagger$]\label{thm:fnangle}
  With $u,v,\theta$ and $f$ as above and $y=f^n$, $2-L_n=(1-\cos\theta)\,\mathbb E[m_y(S_\theta)]$,
  and
  \begin{multline*}
  2-L_n=\frac{(1-\cos\theta)\bigl(\log n+\log f'(1)+k_1(\theta)\bigr)}{f'(1)\,(n-1)}\\
  +O\bigl(n^{-2}\log^2n\bigr).
  \end{multline*}
  \end{theorem}
  ... At $\theta=\frac\pi2$ this is Theorem~\ref{thm:fn} ...

Quoted verbatim (sec_g7_fnangle.tex 6, 8, 10, the setting):

  Put $y=y_n=f^n$, let $L_n$ be the length of the curve $x\mapsto x\,u+y(x)\,v$, $0\le x\le1$,
  and let $m_y(s)=|\{x\in[0,1]:y'(x)>s\}|$. Put
  $h_c(s)=1+s-\sqrt{1+2cs+s^2},\qquad \varphi_\theta(s)=h_c'(s)=1-\frac{s+c}{\sqrt{1+2cs+s^2}}$
  ... The number $k_1(\theta)$ of the statement is the mean of $\log(1/S_\theta)$ ...
  $O(\cdot)$ refers to $n\to\infty$ with constants that depend on $f$ and $\theta$.

Quoted verbatim (sec_g7_fnangle.tex 59-61, Lemma `g7:lem:phi` (a),(b)):

  (a) $\varphi_\theta(s)=\dfrac{1-c^2}{r(s)\,(r(s)+s+c)}>0$, ... $\varphi_\theta(0)=1-c$. So
  $0<\varphi_\theta\le1-c$, and $\int_0^\infty\varphi_\theta\dd s=1-c$.
  (b) For $s\ge1$, $\varphi_\theta(s)\le C_\theta s^{-2}$ with $C_\theta=\frac{1+|c|}{2(1-|c|)}$.

Quoted verbatim (sec_g7_fnangle.tex 75, 79, Lemma `g7:lem:law` (a),(c), the part used):

  $\E[|\log S_\theta|^q]<\infty$ for every $q\ge0$. ...
  In particular $k_1=-a_1=1+\cot^2\frac\theta2\log\cos^2\frac\theta2$

Quoted verbatim (sec_g7_fnangle.tex 99, Lemma `g7:lem:layer`):

  $2-L=\int_0^1h_c(y'(x))\dd x=\int_0^\infty\varphi_\theta(s)\,m_y(s)\dd s
  =(1-\cos\theta)\,\E[m_y(S_\theta)]$.

## Formalization choices

* `u, v` live in the Euclidean plane `EuclideanSpace ℝ (Fin 2)`, with `‖u‖ = ‖v‖ = 1` and
  `inner ℝ u v = cos θ`, `0 < θ < π`. (Nothing depends on the dimension.)
* `f`, `n`, `y_n' = W01.yd f n` (the derivative within `[0,1]`), `m_y = W01.levelMeasure` and the
  hypotheses `W01.Hyp f` (`f ∈ C^1[0,1]` nondecreasing, `f(0)=0`, `f(1)=1`, `f'(1)>0`, `f ∈ C^2`
  near `1`) are exactly those of W01 (`thm:fn`), as the paper says ("$f$ is as in
  Theorem~\ref{thm:fn}"). `n : ℕ`, `f x ^ n` is the natural power.
* Length. `L u v f n = ∫_0^1 ‖u + y_n'(x) v‖ dx`, the integral of the speed of the `C^1` curve
  `x ↦ x u + y_n(x) v` (the supplement, line 104: "the length of a $C^1$ arc is the integral of
  the speed"). The integrand is continuous on `[0,1]`, so this is never the "non-integrable ⇒ 0"
  convention. At `θ = π/2` (`inner u v = 0`) this is `W01.L f n` (stmt_4).
* `φθ c s = 1 - (s+c)/√(1+2cs+s²)` is the paper's `h_c'(s)`; for `|c| < 1` the radicand is
  `(s+c)² + 1 - c² > 0`, so `Real.sqrt` is the honest square root.
* `E[m_y(S_θ)]` is written out as `∫_0^∞ (φθ(s)/(1-c)) m_y(s) ds`, the expectation for the density
  `φθ/(1-c)` of `S_θ` (`1 - c > 0`); the integrability is stated with it.
* `k_1(θ) = 1 + cot²(θ/2) log cos²(θ/2)` (`Real.cot`; `cos(θ/2) > 0` for `0 < θ < π`). Its meaning
  as the mean of `log(1/S_θ)` is stmt_2: `∫_0^∞ φθ(s) log s ds = -(1-c) k_1(θ)`.
* The `O(n^{-2} log² n)` is explicit: `∃ C N, ∀ n ≥ N, |error| ≤ C (log n)² / n²`, with `C, N`
  allowed to depend on `θ, u, v, f` (the supplement: "constants that depend on $f$ and $\theta$").
* The layer-cake lemma (stmt_3) is stated for `y ∈ C^1[0,1]` nondecreasing with `y(1)-y(0)=1`
  (the hypotheses of W01.stmt_3), which covers `y = f^n`; the supplement's Lemma `g7:lem:layer`
  allows `y ∈ C^1[0,1)` only, a generality the theorem does not use.

## Check at `θ = π/2`

`cos(π/2) = 0`, so `1 - cos θ = 1`, `φθ 0 = W01.φ`, `L u v f n = W01.L f n` for orthonormal `u,v`,
and `k_1(π/2) = 1 + cot²(π/4) log cos²(π/4) = 1 + log(1/2) = 1 - log 2` (stmt_4). Hence stmt_6 at
`θ = π/2` is literally `W01.stmt_6`, and stmt_5 at `θ = π/2` is W01's layer-cake identity for
`f^n` (density `W01.φ`, mass 1).

## Statements

* stmt_1: Lemma `g7:lem:phi` (a),(b) and the tail of `S_θ`: `φθ = (1-c²)/(r(r+s+c))`, `0 < φθ ≤ 1-c`,
  `φθ(s) ≤ C_θ/s²` for `s ≥ 1`, `∫_s^∞ φθ = √(1+2cs+s²) - s - c` for `s ≥ 0`, so total mass `1-c`.
* stmt_2: Lemma `g7:lem:law`, the part used: `φθ log`, `φθ log²` integrable on `(0,∞)` and
  `E[log(1/S_θ)] = k_1(θ)`, i.e. `∫ φθ log = -(1-c) k_1(θ)`.
* stmt_3: Lemma `g7:lem:layer`: `2 - L = ∫_0^1 h_c(y') = ∫_0^∞ φθ m_y` (integrable).
* stmt_4: the reduction at `θ = π/2`: `φθ 0 = W01.φ`, `k_1(π/2) = 1 - log 2`, and the length at
  orthogonal unit vectors is `W01.L`.
* stmt_5: the first claim of `thm:fnangle`: `2 - L_n = (1 - cos θ) E[m_{y_n}(S_θ)]`, for `n ≥ 1`.
* stmt_6: the second claim of `thm:fnangle`, with explicit `C, N`.

## NOT COVERED

* The Mellin transform `E[S_θ^w] = (πw/sin πw) Λ_w(θ)/(1-cos θ)`, the representation
  `S_θ = cos²(θ/2)/U - cos θ - sin²(θ/2) U`, the cumulants `k_2, k_3`, the geometric mean
  `e^{a_1}`, the expansion `ε J_1 + ε² J_2 + ε³ J_3 + O(ε⁴ log⁴ n)` and the contact law in the
  text after the theorem (not part of `thm:fnangle`).
* The identification "integral of the speed = length (supremum of polygonal lengths)" is
  explained above, not stated.
* Uniformity of the constants in `θ` (the paper and the supplement do not claim it).
-/
-- BEGIN STATEMENTS
namespace Y03

open MeasureTheory Set Filter

/-- `h_c(s) = 1 + s - √(1+2cs+s²)` (sec_g7_fnangle.tex line 8). -/
noncomputable def hc (c s : ℝ) : ℝ := 1 + s - Real.sqrt (1 + 2 * c * s + s ^ 2)

/-- `φ_θ(s) = h_c'(s) = 1 - (s+c)/√(1+2cs+s²)`, with `c = cos θ` (sec_g7_fnangle.tex line 8). -/
noncomputable def φθ (c s : ℝ) : ℝ := 1 - (s + c) / Real.sqrt (1 + 2 * c * s + s ^ 2)

/-- `k_1(θ) = 1 + cot²(θ/2) log cos²(θ/2)` (main.tex line 1843). -/
noncomputable def k1 (θ : ℝ) : ℝ :=
  1 + Real.cot (θ / 2) ^ 2 * Real.log (Real.cos (θ / 2) ^ 2)

/-- `L_n`: the length `∫_0^1 ‖u + y_n'(x) v‖ dx` of the curve `x ↦ x u + f(x)^n v`. -/
noncomputable def L (u v : EuclideanSpace ℝ (Fin 2)) (f : ℝ → ℝ) (n : ℕ) : ℝ :=
  ∫ x in (0 : ℝ)..1, ‖u + W01.yd f n x • v‖

/-- Lemma `g7:lem:phi` (a),(b) and `P(S_θ > s) = (√(1+2cs+s²) - s - c)/(1-c)`:
for `0 < θ < π`, `c = cos θ`, `s ≥ 0`: `φ_θ(s) = (1-c²)/(r(r+s+c))`, `0 < φ_θ(s) ≤ 1-c`,
`∫_s^∞ φ_θ = r(s) - s - c` (an honest integral; at `s = 0` it is the mass `1-c`), and
`φ_θ(s) ≤ C_θ s^{-2}` for `s ≥ 1`, `C_θ = (1+|c|)/(2(1-|c|))`. -/
def stmt_1 : Prop :=
  ∀ θ : ℝ, 0 < θ → θ < Real.pi →
    let c := Real.cos θ
    (∀ s : ℝ, 0 ≤ s →
      φθ c s = (1 - c ^ 2) / (Real.sqrt (1 + 2 * c * s + s ^ 2) *
        (Real.sqrt (1 + 2 * c * s + s ^ 2) + s + c)) ∧
      0 < φθ c s ∧ φθ c s ≤ 1 - c ∧
      IntegrableOn (φθ c) (Ioi s) ∧
      ∫ t in Ioi s, φθ c t = Real.sqrt (1 + 2 * c * s + s ^ 2) - s - c) ∧
    (∀ s : ℝ, 1 ≤ s → φθ c s ≤ (1 + |c|) / (2 * (1 - |c|)) / s ^ 2)

/-- Lemma `g7:lem:law` (a),(c), the part the theorem uses: `E[log² S_θ] < ∞` and
`E[log(1/S_θ)] = k_1(θ)`, written for the unnormalised density `φ_θ` of mass `1 - cos θ`. -/
def stmt_2 : Prop :=
  ∀ θ : ℝ, 0 < θ → θ < Real.pi →
    let c := Real.cos θ
    IntegrableOn (fun s => φθ c s * Real.log s) (Ioi 0) ∧
    IntegrableOn (fun s => φθ c s * Real.log s ^ 2) (Ioi 0) ∧
    ∫ s in Ioi (0 : ℝ), φθ c s * Real.log s = -((1 - c) * k1 θ)

/-- Lemma `g7:lem:layer` (layer cake at angle `θ`): for unit `u, v` with `u·v = cos θ`, and
`y ∈ C^1[0,1]` nondecreasing with `y(1) - y(0) = 1`, the length `∫_0^1 ‖u + y' v‖` satisfies
`2 - L = ∫_0^1 h_c(y') = ∫_0^∞ φ_θ(s) m_y(s) ds` (the last an honest integral). -/
def stmt_3 : Prop :=
  ∀ θ : ℝ, 0 < θ → θ < Real.pi →
    ∀ u v : EuclideanSpace ℝ (Fin 2), ‖u‖ = 1 → ‖v‖ = 1 → inner ℝ u v = Real.cos θ →
    ∀ y : ℝ → ℝ, ContDiffOn ℝ 1 y (Icc 0 1) → MonotoneOn y (Icc 0 1) → y 1 - y 0 = 1 →
      IntegrableOn (fun s => φθ (Real.cos θ) s * W01.levelMeasure (W01.D y) s) (Ioi 0) ∧
      2 - ∫ x in (0 : ℝ)..1, ‖u + W01.D y x • v‖ =
        ∫ x in (0 : ℝ)..1, hc (Real.cos θ) (W01.D y x) ∧
      2 - ∫ x in (0 : ℝ)..1, ‖u + W01.D y x • v‖ =
        ∫ s in Ioi (0 : ℝ), φθ (Real.cos θ) s * W01.levelMeasure (W01.D y) s

/-- The case `θ = π/2` ("At $\theta=\frac\pi2$ this is Theorem~\ref{thm:fn}"): the density is
W01's `φ`, `k_1(π/2) = 1 - log 2`, and for orthonormal `u, v` the length is W01's `L_n`. -/
def stmt_4 : Prop :=
  φθ (Real.cos (Real.pi / 2)) = W01.φ ∧ k1 (Real.pi / 2) = 1 - Real.log 2 ∧
  ∀ u v : EuclideanSpace ℝ (Fin 2), ‖u‖ = 1 → ‖v‖ = 1 → inner ℝ u v = 0 →
    ∀ f : ℝ → ℝ, ∀ n : ℕ, L u v f n = W01.L f n

/-- `thm:fnangle`, first claim: `2 - L_n = (1 - cos θ) E[m_y(S_θ)]` for `y = f^n`, with
`E[m_y(S_θ)] = ∫_0^∞ (φ_θ(s)/(1 - cos θ)) m_y(s) ds` an honest integral (`n ≥ 1`). -/
def stmt_5 : Prop :=
  ∀ θ : ℝ, 0 < θ → θ < Real.pi →
    ∀ u v : EuclideanSpace ℝ (Fin 2), ‖u‖ = 1 → ‖v‖ = 1 → inner ℝ u v = Real.cos θ →
    ∀ f : ℝ → ℝ, W01.Hyp f → ∀ n : ℕ, 1 ≤ n →
      IntegrableOn
        (fun s => φθ (Real.cos θ) s / (1 - Real.cos θ) * W01.levelMeasure (W01.yd f n) s)
        (Ioi 0) ∧
      2 - L u v f n = (1 - Real.cos θ) *
        ∫ s in Ioi (0 : ℝ),
          φθ (Real.cos θ) s / (1 - Real.cos θ) * W01.levelMeasure (W01.yd f n) s

/-- `thm:fnangle`, second claim:
`2 - L_n = (1 - cos θ)(log n + log f'(1) + k_1(θ))/(f'(1)(n-1)) + O(n^{-2} log² n)`,
with the `O` explicit (`C, N` depend on `θ, u, v, f`). -/
def stmt_6 : Prop :=
  ∀ θ : ℝ, 0 < θ → θ < Real.pi →
    ∀ u v : EuclideanSpace ℝ (Fin 2), ‖u‖ = 1 → ‖v‖ = 1 → inner ℝ u v = Real.cos θ →
    ∀ f : ℝ → ℝ, W01.Hyp f →
      ∃ C : ℝ, ∃ N : ℕ, ∀ n : ℕ, N ≤ n →
        |(2 - L u v f n) - (1 - Real.cos θ) * (Real.log n + Real.log (W01.D f 1) + k1 θ) /
            (W01.D f 1 * ((n : ℝ) - 1))|
          ≤ C * Real.log n ^ 2 / (n : ℝ) ^ 2

end Y03
-- END STATEMENTS

namespace Y03
open MeasureTheory Set Filter Topology

lemma rad_pos (c s : ℝ) (hc : c ^ 2 < 1) : 0 < 1 + 2 * c * s + s ^ 2 := by
  nlinarith [sq_nonneg (s + c)]

lemma r_sq (c s : ℝ) (hc : c ^ 2 < 1) :
    Real.sqrt (1 + 2 * c * s + s ^ 2) ^ 2 = 1 + 2 * c * s + s ^ 2 :=
  Real.sq_sqrt (rad_pos c s hc).le

lemma r_pos (c s : ℝ) (hc : c ^ 2 < 1) : 0 < Real.sqrt (1 + 2 * c * s + s ^ 2) :=
  Real.sqrt_pos.2 (rad_pos c s hc)

lemma r_gt (c s : ℝ) (hc : c ^ 2 < 1) : |s + c| < Real.sqrt (1 + 2 * c * s + s ^ 2) := by
  apply Real.lt_sqrt_of_sq_lt
  rw [sq_abs]; nlinarith

lemma rsc_pos (c s : ℝ) (hc : c ^ 2 < 1) : 0 < Real.sqrt (1 + 2 * c * s + s ^ 2) + s + c := by
  have := r_gt c s hc
  have := neg_abs_le (s + c)
  linarith

lemma phi_eq (c s : ℝ) (hc : c ^ 2 < 1) :
    φθ c s = (1 - c ^ 2) / (Real.sqrt (1 + 2 * c * s + s ^ 2) *
        (Real.sqrt (1 + 2 * c * s + s ^ 2) + s + c)) := by
  have hr := r_pos c s hc
  have h2 := r_sq c s hc
  have h3 := rsc_pos c s hc
  have e : 1 - c ^ 2 = (Real.sqrt (1 + 2 * c * s + s ^ 2) - (s + c)) *
      (Real.sqrt (1 + 2 * c * s + s ^ 2) + s + c) := by linear_combination -h2
  unfold φθ
  rw [e, mul_div_mul_right _ _ h3.ne', sub_div, div_self hr.ne']

lemma phi_pos (c s : ℝ) (hc : c ^ 2 < 1) : 0 < φθ c s := by
  rw [phi_eq c s hc]
  have := r_pos c s hc
  have := rsc_pos c s hc
  have : 0 < 1 - c ^ 2 := by linarith
  positivity

lemma phi_le (c s : ℝ) (hc : c ^ 2 < 1) (hs : 0 ≤ s) : φθ c s ≤ 1 - c := by
  have hr := r_pos c s hc
  have h2 := r_sq c s hc
  set r := Real.sqrt (1 + 2 * c * s + s ^ 2)
  unfold φθ
  have key : c * r ≤ s + c := by
    have hd : (s + c) ^ 2 - c ^ 2 * r ^ 2 = (1 - c ^ 2) * s * (s + 2 * c) := by
      rw [h2]; ring
    rcases le_or_gt 0 c with hc0 | hc0
    · have hsc : 0 ≤ s + c := by linarith
      have : (c * r) ^ 2 ≤ (s + c) ^ 2 := by nlinarith [mul_nonneg (mul_nonneg (by linarith : (0:ℝ) ≤ 1 - c ^ 2) hs) (by linarith : (0:ℝ) ≤ s + 2 * c)]
      by_contra hcon
      push Not at hcon
      nlinarith [mul_nonneg hc0 hr.le]
    · rcases le_or_gt 0 (s + c) with hsc | hsc
      · nlinarith
      · have h1 : 0 ≤ (1 - c ^ 2) * s * (-(s + 2 * c)) :=
          mul_nonneg (mul_nonneg (by linarith) hs) (by linarith)
        have : (s + c) ^ 2 ≤ (c * r) ^ 2 := by nlinarith
        nlinarith [mul_pos (neg_pos.2 hc0) hr]
  rw [sub_le_sub_iff_left, le_div_iff₀ hr]
  linarith

lemma phi_le_inv (c s : ℝ) (hc : c ^ 2 < 1) (hs : 1 ≤ s) :
    φθ c s ≤ (1 + |c|) / (2 * (1 - |c|)) / s ^ 2 := by
  have hr := r_pos c s hc
  have h2 := r_sq c s hc
  have hca : |c| < 1 := by rw [← abs_one]; exact sq_lt_sq.1 (by simpa using hc)
  have hc0 := abs_nonneg c
  have hcl := neg_abs_le c
  rw [phi_eq c s hc]
  set r := Real.sqrt (1 + 2 * c * s + s ^ 2)
  have hrs : s - |c| ≤ r := by
    apply Real.le_sqrt_of_sq_le
    have : c ^ 2 = |c| ^ 2 := (sq_abs c).symm
    nlinarith [le_abs_self c]
  have hA : (1 - |c|) * s ≤ s - |c| := by nlinarith
  have hB : 2 * ((1 - |c|) * s) ≤ r + s + c := by linarith
  have hpos : 0 < (1 - |c|) * s := by apply mul_pos <;> linarith
  have hprod : (1 - |c|) * s * (2 * ((1 - |c|) * s)) ≤ r * (r + s + c) :=
    mul_le_mul (hA.trans hrs) hB (by linarith) hr.le
  have h1c : 1 - c ^ 2 = (1 - |c|) * (1 + |c|) := by rw [← sq_abs c]; ring
  rw [h1c, div_div, div_le_div_iff₀ (by nlinarith) (by positivity)]
  have : 0 < 1 + |c| := by linarith
  nlinarith [mul_le_mul_of_nonneg_left hprod (by positivity : (0:ℝ) ≤ (1 - |c|) * (1 + |c|))]


lemma hasDerivAt_r (c s : ℝ) (hc : c ^ 2 < 1) :
    HasDerivAt (fun t => Real.sqrt (1 + 2 * c * t + t ^ 2))
      ((s + c) / Real.sqrt (1 + 2 * c * s + s ^ 2)) s := by
  have h1 : HasDerivAt (fun t : ℝ => 1 + 2 * c * t + t ^ 2) (2 * c + 2 * s) s := by
    have := ((hasDerivAt_id' s).const_mul (2 * c)).const_add 1 |>.add (hasDerivAt_pow 2 s)
    convert this using 1; simp
  have := h1.sqrt (rad_pos c s hc).ne'
  convert this using 1
  field_simp; ring

lemma hasDerivAt_G0 (c s : ℝ) (hc : c ^ 2 < 1) :
    HasDerivAt (fun t => t - Real.sqrt (1 + 2 * c * t + t ^ 2)) (φθ c s) s := by
  have := (hasDerivAt_id' s).sub (hasDerivAt_r c s hc)
  convert this using 1
  rfl

lemma G0_eq (c t : ℝ) (hc : c ^ 2 < 1) :
    t - Real.sqrt (1 + 2 * c * t + t ^ 2) =
      -c - (1 - c ^ 2) / (Real.sqrt (1 + 2 * c * t + t ^ 2) + t + c) := by
  have h3 := rsc_pos c t hc
  have h2 := r_sq c t hc
  set r := Real.sqrt (1 + 2 * c * t + t ^ 2)
  have : (1 - c ^ 2) / (r + t + c) = r - t - c := by
    rw [div_eq_iff h3.ne']; linear_combination -h2
  rw [this]; ring

lemma tendsto_G0 (c : ℝ) (hc : c ^ 2 < 1) :
    Tendsto (fun t : ℝ => t - Real.sqrt (1 + 2 * c * t + t ^ 2)) atTop (𝓝 (-c)) := by
  have h2 : Tendsto (fun t : ℝ => Real.sqrt (1 + 2 * c * t + t ^ 2) + t + c) atTop atTop := by
    refine tendsto_atTop_mono' _ ?_ (tendsto_atTop_add_const_right _ c tendsto_id)
    exact Eventually.of_forall fun t => by
      simp only [id]; linarith [Real.sqrt_nonneg (1 + 2 * c * t + t ^ 2)]
  have := ((h2.inv_tendsto_atTop).const_mul (1 - c ^ 2)).const_sub (-c)
  simp only [mul_zero, sub_zero] at this
  refine this.congr' (Eventually.of_forall fun t => ?_)
  simp only [Pi.inv_apply]
  rw [G0_eq c t hc, div_eq_mul_inv]

lemma phi_continuous (c : ℝ) (hc : c ^ 2 < 1) : Continuous (φθ c) := by
  show Continuous (fun s : ℝ => 1 - (s + c) / Real.sqrt (1 + 2 * c * s + s ^ 2))
  apply continuous_const.sub
  apply (continuous_id.add continuous_const).div (by fun_prop)
  intro x; exact (r_pos c x hc).ne'

lemma phi_integrableOn (c s : ℝ) (hc : c ^ 2 < 1) : IntegrableOn (φθ c) (Ioi s) :=
  integrableOn_Ioi_deriv_of_nonneg' (fun x _ => hasDerivAt_G0 c x hc)
    (fun x _ => (phi_pos c x hc).le) (tendsto_G0 c hc)

lemma phi_integral (c s : ℝ) (hc : c ^ 2 < 1) :
    ∫ t in Ioi s, φθ c t = Real.sqrt (1 + 2 * c * s + s ^ 2) - s - c := by
  rw [integral_Ioi_of_hasDerivAt_of_nonneg' (fun x _ => hasDerivAt_G0 c x hc)
    (fun x _ => (phi_pos c x hc).le) (tendsto_G0 c hc)]
  ring


lemma alg_H (c s r : ℝ) (hs : 0 < s) (hr : 0 < r) (hr2 : r ^ 2 = 1 + 2 * c * s + s ^ 2)
    (h1 : 0 < 1 + c * s + r) (h2 : 0 < s + c + r) :
    (1 + s - r) / s = 1 - (s + c) / r - c * ((1 + (s + c) / r) / (s + c + r)) +
      (c * 1 + (s + c) / r) / (1 + c * s + r) := by
  have key : (r - 1) * (1 + c * s + r) = s * (c * r + s + c) := by linear_combination hr2
  have e1 : c * ((1 + (s + c) / r) / (s + c + r)) = c / r := by
    field_simp; ring
  have e2 : (c * 1 + (s + c) / r) / (1 + c * s + r) = (r - 1) / (s * r) := by
    rw [div_eq_div_iff h1.ne' (by positivity)]
    field_simp
    linear_combination r * key - (1 + r) * hr2
  rw [e1, e2]
  field_simp
  linear_combination -hr2

lemma r_ge_one_sub (c s : ℝ) (hc : c ^ 2 < 1) (hs : 0 ≤ s) :
    1 - s ≤ Real.sqrt (1 + 2 * c * s + s ^ 2) := by
  apply Real.le_sqrt_of_sq_le; nlinarith

lemma r_le_one_add (c s : ℝ) (hc : c ^ 2 < 1) (hs : 0 ≤ s) :
    Real.sqrt (1 + 2 * c * s + s ^ 2) ≤ 1 + s := by
  rw [Real.sqrt_le_iff]; constructor <;> nlinarith

lemma one_cs_r_pos (c s : ℝ) (hc : c ^ 2 < 1) (hs : 0 ≤ s) :
    0 < 1 + c * s + Real.sqrt (1 + 2 * c * s + s ^ 2) := by
  have hr := r_pos c s hc
  have h2 := r_sq c s hc
  rcases lt_or_ge 0 (1 + c * s) with h | h
  · linarith
  · have hs0 : 0 < s := by
      rcases hs.lt_or_eq with h' | h'
      · exact h'
      · subst h'; linarith
    have : (1 + c * s) ^ 2 < Real.sqrt (1 + 2 * c * s + s ^ 2) ^ 2 := by
      rw [h2]; nlinarith [mul_pos hs0 hs0]
    nlinarith

/-- `H' = h_c(s)/s`. -/
noncomputable def H (c s : ℝ) : ℝ :=
  s - Real.sqrt (1 + 2 * c * s + s ^ 2)
    - c * Real.log (s + c + Real.sqrt (1 + 2 * c * s + s ^ 2))
    + Real.log (1 + c * s + Real.sqrt (1 + 2 * c * s + s ^ 2))

lemma hasDerivAt_H (c s : ℝ) (hcc : c ^ 2 < 1) (hs : 0 < s) :
    HasDerivAt (H c) (Y03.hc c s / s) s := by
  have hR := hasDerivAt_r c s hcc
  have hp1 := rsc_pos c s hcc
  have hp2 := one_cs_r_pos c s hcc hs.le
  have hA : HasDerivAt (fun t => t + c + Real.sqrt (1 + 2 * c * t + t ^ 2))
      (1 + (s + c) / Real.sqrt (1 + 2 * c * s + s ^ 2)) s := by
    have := ((hasDerivAt_id' s).add_const c).add hR
    convert this using 1
  have hB : HasDerivAt (fun t => 1 + c * t + Real.sqrt (1 + 2 * c * t + t ^ 2))
      (c * 1 + (s + c) / Real.sqrt (1 + 2 * c * s + s ^ 2)) s := by
    have := (((hasDerivAt_id' s).const_mul c).const_add 1).add hR
    convert this using 1
  have hA' := (hA.log (by linarith)).const_mul c
  have hB' := hB.log (by linarith)
  have := (((hasDerivAt_id' s).sub hR).sub hA').add hB'
  convert this using 1
  · funext t; simp only [H]; rfl
  · unfold Y03.hc
    have h2 := r_sq c s hcc
    have hr := r_pos c s hcc
    have hp1' : 0 < s + c + Real.sqrt (1 + 2 * c * s + s ^ 2) := by linarith
    generalize Real.sqrt (1 + 2 * c * s + s ^ 2) = r at *
    rw [alg_H c s r hs hr h2 hp2 hp1']


/-- antiderivative of `φθ log`. -/
noncomputable def Gl (c s : ℝ) : ℝ := hc c s * Real.log s - H c s

lemma hasDerivAt_hc (c s : ℝ) (hcc : c ^ 2 < 1) : HasDerivAt (hc c) (φθ c s) s := by
  have := ((hasDerivAt_id' s).const_add 1).sub (hasDerivAt_r c s hcc)
  convert this using 1
  · funext t; simp [hc]
  · rfl

lemma hasDerivAt_Gl (c s : ℝ) (hcc : c ^ 2 < 1) (hs : 0 < s) :
    HasDerivAt (Gl c) (φθ c s * Real.log s) s := by
  have := ((hasDerivAt_hc c s hcc).mul (Real.hasDerivAt_log hs.ne')).sub (hasDerivAt_H c s hcc hs)
  convert this using 1
  · funext t; simp [Gl]
  · rw [div_eq_mul_inv]; ring

lemma Gl_zero (c : ℝ) : Gl c 0 = 1 + c * Real.log (1 + c) - Real.log 2 := by
  unfold Gl H; norm_num
  rw [add_comm c 1]; ring_nf

lemma hc_bounds (c s : ℝ) (hcc : c ^ 2 < 1) (hs : 0 ≤ s) : 0 ≤ hc c s ∧ hc c s ≤ 2 * s := by
  have := r_ge_one_sub c s hcc hs
  have := r_le_one_add c s hcc hs
  unfold hc; constructor <;> linarith

lemma Gl_cont0 (c : ℝ) (hcc : c ^ 2 < 1) : ContinuousWithinAt (Gl c) (Ici 0) 0 := by
  have hA : ContinuousWithinAt (fun s : ℝ => hc c s * Real.log s) (Ici 0) 0 := by
    show Tendsto _ _ _
    simp only [Real.log_zero, mul_zero]
    have hB : Tendsto (fun s : ℝ => 2 * (s * Real.log s)) (𝓝[Ici 0] 0) (𝓝 0) := by
      have := ((Real.continuous_mul_log.tendsto 0).mono_left
        (nhdsWithin_le_nhds (s := Ici (0 : ℝ)))).const_mul 2
      simpa using this
    refine squeeze_zero_norm' ?_ hB.norm
    filter_upwards [self_mem_nhdsWithin] with s hs
    have hs0 : (0 : ℝ) ≤ s := hs
    obtain ⟨h1, h2⟩ := hc_bounds c s hcc hs0
    simp only [norm_mul, Real.norm_eq_abs]
    rw [abs_of_nonneg h1, abs_of_nonneg (by norm_num : (0:ℝ) ≤ 2), abs_of_nonneg hs0]
    nlinarith [abs_nonneg (Real.log s)]
  have hr : Continuous (fun s : ℝ => Real.sqrt (1 + 2 * c * s + s ^ 2)) := by fun_prop
  have hH : ContinuousAt (H c) 0 := by
    unfold H
    have h1 : ContinuousAt (fun s : ℝ => s + c + Real.sqrt (1 + 2 * c * s + s ^ 2)) 0 :=
      (by fun_prop : Continuous fun s : ℝ => s + c + Real.sqrt (1 + 2 * c * s + s ^ 2)).continuousAt
    have h2 : ContinuousAt (fun s : ℝ => 1 + c * s + Real.sqrt (1 + 2 * c * s + s ^ 2)) 0 :=
      (by fun_prop : Continuous fun s : ℝ => 1 + c * s + Real.sqrt (1 + 2 * c * s + s ^ 2)).continuousAt
    have hp1 := rsc_pos c 0 hcc
    have hp2 := one_cs_r_pos c 0 hcc le_rfl
    exact ((continuous_id.sub hr).continuousAt.sub
      ((h1.log (by norm_num at hp1 ⊢ <;> linarith)).const_mul c)).add
      (h2.log (by norm_num at hp2 ⊢))
  exact hA.sub hH.continuousWithinAt


lemma Gl_eq (c s : ℝ) (hcc : c ^ 2 < 1) (hs : 0 < s) :
    Gl c s = (s + c - Real.sqrt (1 + 2 * c * s + s ^ 2)) * Real.log s
      + (Real.sqrt (1 + 2 * c * s + s ^ 2) - s)
      + c * Real.log ((s + c + Real.sqrt (1 + 2 * c * s + s ^ 2)) / s)
      - Real.log ((1 + c * s + Real.sqrt (1 + 2 * c * s + s ^ 2)) / s) := by
  have hp1 := rsc_pos c s hcc
  have hp2 := one_cs_r_pos c s hcc hs.le
  unfold Gl H hc
  rw [Real.log_div (by linarith) hs.ne', Real.log_div (by linarith) hs.ne']
  ring

lemma tendsto_r_sub (c : ℝ) (hcc : c ^ 2 < 1) :
    Tendsto (fun t : ℝ => Real.sqrt (1 + 2 * c * t + t ^ 2) - t) atTop (𝓝 c) := by
  have := (tendsto_G0 c hcc).neg
  rw [neg_neg] at this
  exact this.congr (fun t => by ring)

lemma tendsto_ratio (c d : ℝ) (hcc : c ^ 2 < 1) :
    Tendsto (fun t : ℝ => (d + t + Real.sqrt (1 + 2 * c * t + t ^ 2)) / t) atTop (𝓝 2) := by
  have h1 : Tendsto (fun t : ℝ => (d + (Real.sqrt (1 + 2 * c * t + t ^ 2) - t)) / t) atTop
      (𝓝 0) := ((tendsto_r_sub c hcc).const_add d).div_atTop tendsto_id
  have := h1.const_add 2
  rw [add_zero] at this
  refine this.congr' ?_
  filter_upwards [eventually_gt_atTop 0] with t ht
  field_simp; ring

lemma Gl_tendsto (c : ℝ) (hcc : c ^ 2 < 1) :
    Tendsto (Gl c) atTop (𝓝 (c + c * Real.log 2 - Real.log (1 + c))) := by
  have hc1 : 0 < 1 + c := by nlinarith
  -- piece 1
  have h1 : Tendsto (fun s : ℝ => (s + c - Real.sqrt (1 + 2 * c * s + s ^ 2)) * Real.log s)
      atTop (𝓝 0) := by
    have hL : Tendsto (fun s : ℝ => Real.log s ^ 1 / (1 * s + 0)) atTop (𝓝 0) :=
      Real.tendsto_pow_log_div_mul_add_atTop 1 0 1 one_ne_zero
    refine squeeze_zero_norm' ?_ hL
    filter_upwards [eventually_ge_atTop 2] with s hs
    have hs0 : 0 < s := by linarith
    have hca : |c| < 1 := by rw [← abs_one]; exact sq_lt_sq.1 (by simpa using hcc)
    have hcl := neg_abs_le c
    have hcu := le_abs_self c
    have hl := Real.log_nonneg (by linarith : (1:ℝ) ≤ s)
    have hrr := r_gt c s hcc
    have hsc := le_abs_self (s + c)
    have e : s + c - Real.sqrt (1 + 2 * c * s + s ^ 2) =
        -((1 - c ^ 2) / (Real.sqrt (1 + 2 * c * s + s ^ 2) + s + c)) := by
      have := G0_eq c s hcc; linarith
    have hp := rsc_pos c s hcc
    have hbd : (1 - c ^ 2) / (Real.sqrt (1 + 2 * c * s + s ^ 2) + s + c) ≤ 1 / s := by
      rw [div_le_div_iff₀ hp hs0]; nlinarith
    have hnn : 0 ≤ (1 - c ^ 2) / (Real.sqrt (1 + 2 * c * s + s ^ 2) + s + c) :=
      div_nonneg (by linarith) hp.le
    rw [e, norm_mul, norm_neg, Real.norm_eq_abs, Real.norm_eq_abs, abs_of_nonneg hnn,
      abs_of_nonneg hl]
    simp only [pow_one, one_mul, add_zero]
    calc _ ≤ 1 / s * Real.log s := mul_le_mul_of_nonneg_right hbd hl
      _ = Real.log s / s := by ring
  have h2 := tendsto_r_sub c hcc
  have h3 : Tendsto (fun s : ℝ => c * Real.log ((s + c + Real.sqrt (1 + 2 * c * s + s ^ 2)) / s))
      atTop (𝓝 (c * Real.log 2)) := by
    have := ((Real.continuousAt_log (by norm_num : (2:ℝ) ≠ 0)).tendsto.comp
      (tendsto_ratio c c hcc)).const_mul c
    refine this.congr' ?_
    filter_upwards with s
    simp only [Function.comp]; ring_nf
  have h4 : Tendsto (fun s : ℝ => Real.log ((1 + c * s + Real.sqrt (1 + 2 * c * s + s ^ 2)) / s))
      atTop (𝓝 (Real.log (1 + c))) := by
    have hq : Tendsto (fun t : ℝ => (1 + (Real.sqrt (1 + 2 * c * t + t ^ 2) - t)) / t) atTop
        (𝓝 0) := ((tendsto_r_sub c hcc).const_add 1).div_atTop tendsto_id
    have hq2 := hq.const_add (1 + c)
    rw [add_zero] at hq2
    have := (Real.continuousAt_log hc1.ne').tendsto.comp hq2
    refine this.congr' ?_
    filter_upwards [eventually_gt_atTop 0] with t ht
    simp only [Function.comp]
    congr 1
    field_simp; ring
  have := ((h1.add h2).add h3).sub h4
  rw [zero_add] at this
  refine this.congr' ?_
  filter_upwards [eventually_gt_atTop 0] with s hs
  rw [Gl_eq c s hcc hs]


/-- the constant `C_θ`-type bound `φθ ≤ K/s²` for all `s > 0`. -/
noncomputable def Kc (c : ℝ) : ℝ := max 2 ((1 + |c|) / (2 * (1 - |c|)))

lemma phi_le_K (c s : ℝ) (hcc : c ^ 2 < 1) (hs : 0 < s) : φθ c s ≤ Kc c / s ^ 2 := by
  rcases le_or_gt 1 s with h | h
  · calc φθ c s ≤ (1 + |c|) / (2 * (1 - |c|)) / s ^ 2 := phi_le_inv c s hcc h
      _ ≤ Kc c / s ^ 2 := by gcongr; exact le_max_right _ _
  · have h1 := phi_le c s hcc hs.le
    have hc1 : 1 - c ≤ 2 := by nlinarith
    have hK : 2 ≤ Kc c := le_max_left _ _
    have : 2 ≤ Kc c / s ^ 2 := by
      rw [le_div_iff₀ (by positivity)]
      have : s ^ 2 ≤ 1 := by nlinarith
      nlinarith
    unfold Kc at *; linarith

lemma Kc_pos (c : ℝ) : 0 < Kc c := lt_of_lt_of_le (by norm_num) (le_max_left _ _)

lemma phi_log_cont (c : ℝ) (hcc : c ^ 2 < 1) :
    ContinuousOn (fun s => φθ c s * Real.log s ^ 2) (Ioi 0) := by
  apply ContinuousOn.mul
  · exact (phi_continuous c hcc).continuousOn
  · exact (Real.continuousOn_log.mono (fun x hx => ne_of_gt hx)).pow 2

lemma phi_log2_integrable (c : ℝ) (hcc : c ^ 2 < 1) :
    IntegrableOn (fun s => φθ c s * Real.log s ^ 2) (Ioi 0) := by
  rw [← Ioc_union_Ioi_eq_Ioi (zero_le_one' ℝ)]
  apply IntegrableOn.union
  · have hint : IntegrableOn (fun s : ℝ => 32 * s ^ (-1 / 2 : ℝ)) (Ioc 0 1) :=
      W01.rpow_half_integrable.const_mul 32
    refine hint.mono' (((phi_log_cont c hcc).mono Ioc_subset_Ioi_self).aestronglyMeasurable
      measurableSet_Ioc) ?_
    refine (ae_restrict_iff' measurableSet_Ioc).2 (Eventually.of_forall fun s hs => ?_)
    have h1 := phi_pos c s hcc
    have h2 := phi_le c s hcc hs.1.le
    have hc1 : 1 - c ≤ 2 := by nlinarith
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
    have := W01.log_sq_le_small s hs.1 hs.2
    have h0 : 0 ≤ s ^ (-1 / 2 : ℝ) := Real.rpow_nonneg hs.1.le _
    nlinarith [sq_nonneg (Real.log s)]
  · have hint : IntegrableOn (fun s : ℝ => 16 * Kc c * s ^ (-3 / 2 : ℝ)) (Ioi 1) :=
      (integrableOn_Ioi_rpow_of_lt (by norm_num) zero_lt_one).const_mul _
    refine hint.mono' (((phi_log_cont c hcc).mono (Ioi_subset_Ioi zero_le_one)).aestronglyMeasurable
      measurableSet_Ioi) ?_
    refine (ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun s hs => ?_)
    have hs' : (0 : ℝ) < s := lt_trans zero_lt_one hs
    have h1 := phi_pos c s hcc
    have h2 := phi_le_K c s hcc hs'
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
    have := W01.log_sq_le_large s hs.le
    obtain ⟨hu, hsu, -, h3⟩ := W01.rpow_quarter_facts s hs'
    rw [h3]
    set u := s ^ (1 / 4 : ℝ)
    have hK := Kc_pos c
    calc φθ c s * Real.log s ^ 2 ≤ Kc c / s ^ 2 * (16 * u ^ 2) :=
          mul_le_mul h2 this (sq_nonneg _) (by positivity)
      _ = 16 * Kc c * (u ^ 6)⁻¹ := by rw [hsu]; field_simp

lemma phi_log_integrable (c : ℝ) (hcc : c ^ 2 < 1) :
    IntegrableOn (fun s => φθ c s * Real.log s) (Ioi 0) := by
  refine ((phi_integrableOn c 0 hcc).add (phi_log2_integrable c hcc)).mono' ?_ ?_
  · apply ContinuousOn.aestronglyMeasurable _ measurableSet_Ioi
    apply ContinuousOn.mul
    · exact (phi_continuous c hcc).continuousOn
    · exact Real.continuousOn_log.mono (fun x hx => ne_of_gt hx)
  · refine (ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun s hs => ?_)
    have h1 := phi_pos c s hcc
    rw [Real.norm_eq_abs, abs_mul, abs_of_pos h1]
    have : |Real.log s| ≤ 1 + Real.log s ^ 2 := by
      rcases abs_cases (Real.log s) with ⟨h, -⟩ | ⟨h, -⟩ <;> rw [h] <;>
        nlinarith [sq_nonneg (Real.log s - 1), sq_nonneg (Real.log s + 1)]
    show _ ≤ φθ c s + φθ c s * Real.log s ^ 2
    nlinarith

lemma phi_log_integral (c : ℝ) (hcc : c ^ 2 < 1) :
    ∫ s in Ioi (0 : ℝ), φθ c s * Real.log s =
      c - 1 + (1 + c) * Real.log 2 - (1 + c) * Real.log (1 + c) := by
  rw [integral_Ioi_of_hasDerivAt_of_tendsto (Gl_cont0 c hcc) (fun x hx => hasDerivAt_Gl c x hcc hx)
    (phi_log_integrable c hcc) (Gl_tendsto c hcc), Gl_zero]
  ring

lemma cos_sq_lt (θ : ℝ) (h0 : 0 < θ) (h1 : θ < Real.pi) : Real.cos θ ^ 2 < 1 := by
  have := Real.sin_pos_of_pos_of_lt_pi h0 h1
  have := Real.sin_sq_add_cos_sq θ
  nlinarith

lemma k1_eq (θ : ℝ) (h0 : 0 < θ) (h1 : θ < Real.pi) :
    -((1 - Real.cos θ) * k1 θ) = Real.cos θ - 1 + (1 + Real.cos θ) * Real.log 2
      - (1 + Real.cos θ) * Real.log (1 + Real.cos θ) := by
  have hs : 0 < Real.sin (θ / 2) :=
    Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith)
  have hcs : Real.cos (θ / 2) ^ 2 = 1 / 2 + Real.cos θ / 2 := by
    rw [Real.cos_sq]; ring_nf
  have hss : Real.sin (θ / 2) ^ 2 = 1 / 2 - Real.cos θ / 2 := by
    rw [Real.sin_sq, hcs]; ring
  have hcc := cos_sq_lt θ h0 h1
  have hc1 : 0 < 1 + Real.cos θ := by nlinarith
  unfold k1
  rw [Real.cot_eq_cos_div_sin, div_pow, hcs, hss]
  have : (1 / 2 + Real.cos θ / 2) = (1 + Real.cos θ) / 2 := by ring
  rw [this, Real.log_div hc1.ne' (by norm_num)]
  have h1c : (1 / 2 - Real.cos θ / 2) ≠ 0 := by nlinarith
  have h1c' : 1 - Real.cos θ ≠ 0 := by nlinarith
  field_simp
  ring


lemma norm_uv (u v : EuclideanSpace ℝ (Fin 2)) (c : ℝ) (hu : ‖u‖ = 1) (hv : ‖v‖ = 1)
    (huv : inner ℝ u v = c) (t : ℝ) : ‖u + t • v‖ = Real.sqrt (1 + 2 * c * t + t ^ 2) := by
  have h := norm_add_sq_real u (t • v)
  rw [inner_smul_right, huv, norm_smul, hu, hv, Real.norm_eq_abs, mul_one, sq_abs] at h
  have e : 1 + 2 * c * t + t ^ 2 = ‖u + t • v‖ ^ 2 := by rw [h]; ring
  rw [e, Real.sqrt_sq (norm_nonneg _)]

lemma level_integrable (c : ℝ) (hcc : c ^ 2 < 1) (g : ℝ → ℝ) :
    IntegrableOn (fun s => φθ c s * W01.levelMeasure g s) (Ioi 0) := by
  refine (phi_integrableOn c 0 hcc).mono' ?_ ?_
  · exact ((phi_continuous c hcc).measurable.mul
      (W01.levelMeasure_antitone g).measurable).aestronglyMeasurable
  · refine (ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun s hs => ?_)
    have h1 := phi_pos c s hcc
    obtain ⟨h2, h3⟩ := W01.levelMeasure_le_one g s
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
    nlinarith

lemma integral_phi_eq (c : ℝ) (hcc : c ^ 2 < 1) (t : ℝ) :
    ∫ s in (0 : ℝ)..t, φθ c s = hc c t := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt (fun x _ => hasDerivAt_G0 c x hcc)
    ((phi_continuous c hcc).intervalIntegrable _ _)]
  unfold hc; norm_num; ring

lemma layer (c : ℝ) (hcc : c ^ 2 < 1) (y : ℝ → ℝ) (hy : ContDiffOn ℝ 1 y (Icc 0 1))
    (hm : MonotoneOn y (Icc 0 1)) (h1 : y 1 - y 0 = 1) :
    (2 - ∫ x in (0 : ℝ)..1, Real.sqrt (1 + 2 * c * W01.D y x + W01.D y x ^ 2) =
      ∫ x in (0 : ℝ)..1, hc c (W01.D y x)) ∧
    ∫ x in (0 : ℝ)..1, hc c (W01.D y x) =
      ∫ s in Ioi (0 : ℝ), φθ c s * W01.levelMeasure (W01.D y) s := by
  set g := W01.D y with hg_def
  have hgc : ContinuousOn g (Icc 0 1) :=
    hy.continuousOn_derivWithin (uniqueDiffOn_Icc zero_lt_one) le_rfl
  have hg0 : ∀ x, 0 ≤ g x := fun x => hm.derivWithin_nonneg
  have hint1 : ∫ x in (0 : ℝ)..1, g x = 1 := by
    have := intervalIntegral.integral_derivWithin_Icc_of_contDiffOn_Icc (a := (0 : ℝ)) (b := 1)
      hy zero_le_one
    rw [h1] at this
    exact this
  have hsc : ContinuousOn (fun x => Real.sqrt (1 + 2 * c * g x + g x ^ 2)) (Icc 0 1) :=
    ((continuousOn_const.add (continuousOn_const.mul hgc)).add (hgc.pow 2)).sqrt
  set F : ℝ → ℝ := fun x => hc c (g x) with hF
  have hF0 : ∀ x, 0 ≤ F x := fun x => (hc_bounds c (g x) hcc (hg0 x)).1
  have hFc : ContinuousOn F (Icc 0 1) := (continuousOn_const.add hgc).sub hsc
  have hgi : IntervalIntegrable g volume 0 1 :=
    (hgc.mono (by rw [uIcc_of_le zero_le_one])).intervalIntegrable
  have hsi : IntervalIntegrable (fun x => Real.sqrt (1 + 2 * c * g x + g x ^ 2)) volume 0 1 :=
    (hsc.mono (by rw [uIcc_of_le zero_le_one])).intervalIntegrable
  have hL : 2 - ∫ x in (0 : ℝ)..1, Real.sqrt (1 + 2 * c * g x + g x ^ 2) =
      ∫ x in (0 : ℝ)..1, F x := by
    simp only [hF, hc]
    rw [intervalIntegral.integral_sub ((intervalIntegrable_const).add hgi) hsi,
      intervalIntegral.integral_add intervalIntegrable_const hgi, hint1]
    norm_num
  refine ⟨hL, ?_⟩
  rw [intervalIntegral.integral_of_le zero_le_one, ← integral_Icc_eq_integral_Ioc]
  set μ : Measure ℝ := volume.restrict (Icc (0 : ℝ) 1)
  have hcake := lintegral_comp_eq_lintegral_meas_lt_mul μ (f := g) (g := φθ c)
    (Eventually.of_forall hg0) (hgc.aemeasurable measurableSet_Icc)
    (fun t _ => (phi_continuous c hcc).intervalIntegrable _ _)
    ((ae_restrict_iff' measurableSet_Ioi).2
      (Eventually.of_forall fun t _ => (phi_pos c t hcc).le))
  have hFeq : ∀ x, F x = ∫ s in (0 : ℝ)..g x, φθ c s := by
    intro x; rw [integral_phi_eq c hcc]
  rw [integral_eq_lintegral_of_nonneg_ae (Eventually.of_forall hF0)
    (hFc.aestronglyMeasurable measurableSet_Icc)]
  rw [integral_eq_lintegral_of_nonneg_ae
    ((ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun t _ =>
      mul_nonneg (phi_pos c t hcc).le (W01.levelMeasure_le_one g t).1))
    (level_integrable c hcc g).aestronglyMeasurable]
  congr 1
  simp_rw [hFeq]
  rw [hcake]
  apply setLIntegral_congr_fun measurableSet_Ioi
  intro t ht
  have hphi := (phi_pos c t hcc).le
  show μ {a | t < g a} * ENNReal.ofReal (φθ c t) = ENNReal.ofReal (φθ c t * W01.levelMeasure g t)
  rw [ENNReal.ofReal_mul hphi, mul_comm]
  congr 1
  unfold W01.levelMeasure
  rw [ENNReal.ofReal_toReal (ne_top_of_le_ne_top (b := volume (Icc (0 : ℝ) 1)) (by simp)
      (measure_mono (s := {x : ℝ | x ∈ Icc (0 : ℝ) 1 ∧ t < g x}) (fun x hx => hx.1))),
    Measure.restrict_apply' measurableSet_Icc]
  congr 1
  ext x; simp [and_comm]

end Y03

namespace Y03
open MeasureTheory Set Filter Topology

/-- The properties of a density `ψ` used in the asymptotic argument of W01. -/
structure Dens (ψ : ℝ → ℝ) (k K : ℝ) : Prop where
  cont : Continuous ψ
  pos : ∀ s, 0 < s → 0 < ψ s
  le1 : ∀ s, 0 < s → ψ s ≤ 1
  Kpos : 0 < K
  tail : ∀ s, 0 < s → ψ s ≤ K / s ^ 2
  int : IntegrableOn ψ (Ioi 0)
  int1 : ∫ s in Ioi (0:ℝ), ψ s = 1
  lint : IntegrableOn (fun s => ψ s * Real.log s) (Ioi 0)
  l2int : IntegrableOn (fun s => ψ s * Real.log s ^ 2) (Ioi 0)
  lval : ∫ s in Ioi (0:ℝ), ψ s * Real.log s = -k

lemma hasDerivAt_Ft (K A s : ℝ) (hs : 0 < s) :
    HasDerivAt (fun t => -K * (Real.log t - Real.log A + 1) / t)
      (K * (Real.log s - Real.log A) / s ^ 2) s := by
  have h1 := (((Real.hasDerivAt_log hs.ne').sub_const (Real.log A)).add_const 1).const_mul (-K)
  have := h1.div (hasDerivAt_id' s) hs.ne'
  convert this using 1
  field_simp
  ring

lemma tendsto_Ft (K A : ℝ) :
    Tendsto (fun t => -K * (Real.log t - Real.log A + 1) / t) atTop (𝓝 0) := by
  have hL : Tendsto (fun s : ℝ => Real.log s ^ 1 / (1 * s + 0)) atTop (𝓝 0) :=
    Real.tendsto_pow_log_div_mul_add_atTop 1 0 1 one_ne_zero
  have h2 : Tendsto (fun s : ℝ => (1 - Real.log A) / s) atTop (𝓝 0) :=
    tendsto_const_nhds.div_atTop tendsto_id
  have := (hL.add h2).const_mul (-K)
  rw [add_zero, mul_zero] at this
  refine this.congr' ?_
  filter_upwards [eventually_gt_atTop 0] with s hs
  simp only [pow_one, one_mul, add_zero]
  field_simp
  ring

lemma tail_bound {ψ : ℝ → ℝ} {k K : ℝ} (hD : Dens ψ k K) (A : ℝ) (hA : 0 < A) :
    IntegrableOn (fun s => ψ s * (Real.log s - Real.log A)) (Ioi A) ∧
      ∫ s in Ioi A, ψ s * (Real.log s - Real.log A) ≤ K / A := by
  have hsub : Ioi A ⊆ Ioi (0:ℝ) := Ioi_subset_Ioi hA.le
  have i1 : IntegrableOn (fun s => ψ s * Real.log s) (Ioi A) := hD.lint.mono_set hsub
  have i2 : IntegrableOn (fun s => ψ s * Real.log A) (Ioi A) :=
    (hD.int.mono_set hsub).mul_const _
  have hint : IntegrableOn (fun s => ψ s * (Real.log s - Real.log A)) (Ioi A) := by
    have := i1.sub i2
    refine this.congr_fun (fun s _ => ?_) measurableSet_Ioi
    simp only [Pi.sub_apply]; ring
  refine ⟨hint, ?_⟩
  have hder : ∀ x ∈ Ici A, HasDerivAt (fun t => -K * (Real.log t - Real.log A + 1) / t)
      (K * (Real.log x - Real.log A) / x ^ 2) x :=
    fun x hx => hasDerivAt_Ft K A x (lt_of_lt_of_le hA hx)
  have hnn : ∀ x ∈ Ioi A, 0 ≤ K * (Real.log x - Real.log A) / x ^ 2 := by
    intro x hx
    have : Real.log A ≤ Real.log x := Real.log_le_log hA (le_of_lt hx)
    have := hD.Kpos
    apply div_nonneg _ (sq_nonneg _)
    apply mul_nonneg this.le; linarith
  have gi := integrableOn_Ioi_deriv_of_nonneg' hder hnn (tendsto_Ft K A)
  have gv := integral_Ioi_of_hasDerivAt_of_nonneg' hder hnn (tendsto_Ft K A)
  have : ∫ s in Ioi A, ψ s * (Real.log s - Real.log A) ≤
      ∫ s in Ioi A, K * (Real.log s - Real.log A) / s ^ 2 := by
    apply setIntegral_mono_on hint gi measurableSet_Ioi
    intro s hs
    have hs0 : 0 < s := lt_trans hA hs
    have h0 : 0 ≤ Real.log s - Real.log A := by
      have := Real.log_le_log hA (le_of_lt hs); linarith
    calc ψ s * (Real.log s - Real.log A) ≤ K / s ^ 2 * (Real.log s - Real.log A) :=
          mul_le_mul_of_nonneg_right (hD.tail s hs0) h0
      _ = _ := by ring
  rw [gv] at this
  have e : (0:ℝ) - -K * (Real.log A - Real.log A + 1) / A = K / A := by ring
  linarith

/-- The `T` identity. -/
lemma T_eq {ψ : ℝ → ℝ} {k K : ℝ} (hD : Dens ψ k K) (n a : ℝ) (hn : 1 < n) (ha : 0 < a) :
    IntegrableOn (fun s => ψ s * (1 / (n - 1) * (Real.log (n / s) + Real.log a) / a)) (Ioi 0) ∧
    (Real.log n + Real.log a + k) / (a * (n - 1)) =
      ∫ s in Ioi (0:ℝ), ψ s * (1 / (n - 1) * (Real.log (n / s) + Real.log a) / a) := by
  set c1 := 1 / (n - 1) * (Real.log n + Real.log a) / a
  set c2 := 1 / (n - 1) / a
  have hcongr : EqOn (fun s => ψ s * (1 / (n - 1) * (Real.log (n / s) + Real.log a) / a))
      (fun s => c1 * ψ s - c2 * (ψ s * Real.log s)) (Ioi 0) := by
    intro s hs
    simp only
    rw [Real.log_div (by linarith) (ne_of_gt hs)]
    simp only [c1, c2]; ring
  have hi : IntegrableOn (fun s => c1 * ψ s - c2 * (ψ s * Real.log s)) (Ioi 0) :=
    (hD.int.const_mul c1).sub (hD.lint.const_mul c2)
  refine ⟨hi.congr_fun hcongr.symm measurableSet_Ioi, ?_⟩
  rw [setIntegral_congr_fun measurableSet_Ioi hcongr, integral_sub
    (hD.int.const_mul c1) (hD.lint.const_mul c2), integral_const_mul,
    integral_const_mul, hD.int1, hD.lval]
  have : n - 1 ≠ 0 := by linarith
  simp only [c1, c2]
  field_simp
  ring

/-- The three dominating functions. -/
noncomputable def B1 (ψ : ℝ → ℝ) (n C1 μ : ℝ) (s : ℝ) : ℝ :=
  3 * C1 * μ ^ 2 * (ψ s * (1 + Real.log n ^ 2) + ψ s * Real.log s ^ 2)
noncomputable def B3 (ψ : ℝ → ℝ) (n a μ : ℝ) (s : ℝ) : ℝ :=
  (Ioi (n * a)).indicator (fun s => μ / a * (ψ s * (Real.log s - Real.log (n * a)))) s

lemma pt_bound {ψ : ℝ → ℝ} {k K : ℝ} (hD : Dens ψ k K) (n a μ C1 σ mv s : ℝ) (hs : 0 < s) (ha : 0 < a) (hμ0 : 0 < μ) (hμ1 : μ ≤ 1)
    (hn1 : 1 ≤ n) (hC1 : 0 ≤ C1) (hσ0 : 0 < σ) (hσ1 : σ ≤ 1) (hm0 : 0 ≤ mv) (hm1 : mv ≤ 1)
    (hmA : n * a ≤ s → mv = 0)
    (hwin : σ ≤ s → s < n * a →
      |mv - μ * (Real.log (n / s) + Real.log a) / a| ≤ C1 * μ ^ 2 * (1 + |Real.log (n / s)|) ^ 2) :
    ‖ψ s * (mv - μ * (Real.log (n / s) + Real.log a) / a)‖ ≤
      B1 ψ n C1 μ s + W01.B2 n (1 + (1 + |Real.log a|) / a) σ s + B3 ψ n a μ s := by
  have hφ0 := hD.pos s hs
  have hφ1 := hD.le1 s hs
  have hln : 0 ≤ Real.log n := Real.log_nonneg hn1
  have hlam : Real.log (n / s) = Real.log n - Real.log s :=
    Real.log_div (by linarith) hs.ne'
  have hB1 : 0 ≤ B1 ψ n C1 μ s := by unfold B1; positivity
  have hB2 : 0 ≤ W01.B2 n (1 + (1 + |Real.log a|) / a) σ s := by
    unfold W01.B2; apply indicator_nonneg; intro x hx
    have : 0 < x := hx.1
    positivity
  have hB3 : 0 ≤ B3 ψ n a μ s := by
    unfold B3; apply indicator_nonneg; intro x hx
    have hx' : n * a < x := hx
    have : Real.log (n * a) ≤ Real.log x := Real.log_le_log (by positivity) hx'.le
    have : 0 < ψ x := hD.pos x (by nlinarith)
    have : 0 ≤ Real.log x - Real.log (n * a) := by linarith
    positivity
  rw [Real.norm_eq_abs, abs_mul, abs_of_pos hφ0]
  rcases lt_or_ge s σ with h1 | h1
  · -- small s
    have hsI : s ∈ Ioc (0:ℝ) 1 := ⟨hs, by linarith⟩
    have hB2v : W01.B2 n (1 + (1 + |Real.log a|) / a) σ s = 5 * (1 + (1 + |Real.log a|) / a) *
        (1 + Real.log n) * σ ^ (1 / 4 : ℝ) * s ^ (-1 / 2 : ℝ) := by
      unfold W01.B2; rw [indicator_of_mem hsI]
    have hsm := W01.small_bound s σ hs h1.le hσ1
    set K1 := 1 + (1 + |Real.log a|) / a
    set ls := |Real.log s|
    have hh : |μ * (Real.log (n / s) + Real.log a) / a| ≤ (Real.log n + ls + |Real.log a|) / a := by
      rw [abs_div, abs_of_pos ha, abs_mul, abs_of_pos hμ0]
      apply div_le_div_of_nonneg_right _ ha.le
      have : |Real.log (n / s) + Real.log a| ≤ Real.log n + ls + |Real.log a| := by
        rw [hlam]
        calc _ ≤ |Real.log n - Real.log s| + |Real.log a| := abs_add_le _ _
          _ ≤ |Real.log n| + |Real.log s| + |Real.log a| := by gcongr; exact abs_sub _ _
          _ = _ := by rw [abs_of_nonneg hln]
      calc μ * |Real.log (n / s) + Real.log a| ≤ 1 * (Real.log n + ls + |Real.log a|) :=
            mul_le_mul hμ1 this (abs_nonneg _) zero_le_one
        _ = _ := one_mul _
    have hdiff : |mv - μ * (Real.log (n / s) + Real.log a) / a| ≤
        1 + (Real.log n + ls + |Real.log a|) / a := by
      calc _ ≤ |mv| + |μ * (Real.log (n / s) + Real.log a) / a| := abs_sub _ _
        _ ≤ _ := by rw [abs_of_nonneg hm0]; linarith
    have hK : 1 + (Real.log n + ls + |Real.log a|) / a ≤ K1 * (1 + Real.log n) * (1 + ls) := by
      have hls : 0 ≤ ls := abs_nonneg _
      have hla : 0 ≤ |Real.log a| := abs_nonneg _
      have e : K1 * (1 + Real.log n) * (1 + ls) = (1 + Real.log n) * (1 + ls) +
          ((1 + |Real.log a|) * (1 + Real.log n) * (1 + ls)) / a := by
        simp only [K1]; field_simp
      rw [e]
      have h2 : Real.log n + ls + |Real.log a| ≤ (1 + |Real.log a|) * (1 + Real.log n) * (1 + ls) := by
        nlinarith [mul_nonneg hln hls, mul_nonneg hla hln, mul_nonneg hla hls,
          mul_nonneg (mul_nonneg hla hln) hls]
      have h3 : (Real.log n + ls + |Real.log a|) / a ≤
          ((1 + |Real.log a|) * (1 + Real.log n) * (1 + ls)) / a :=
        div_le_div_of_nonneg_right h2 ha.le
      nlinarith [mul_nonneg hln hls]
    have hK1 : 0 ≤ K1 := by positivity
    calc ψ s * |mv - μ * (Real.log (n / s) + Real.log a) / a| ≤
          1 * (K1 * (1 + Real.log n) * (1 + ls)) :=
          mul_le_mul hφ1 (hdiff.trans hK) (abs_nonneg _) zero_le_one
      _ ≤ K1 * (1 + Real.log n) * (5 * σ ^ (1 / 4 : ℝ) * s ^ (-1 / 2 : ℝ)) := by
          rw [one_mul]; gcongr
      _ = W01.B2 n K1 σ s := by rw [hB2v]; ring
      _ ≤ _ := by linarith
  rcases lt_or_ge s (n * a) with h2 | h2
  · -- window
    have hw := hwin h1 h2
    have hl : |Real.log (n / s)| ≤ Real.log n + |Real.log s| := by
      rw [hlam]
      calc _ ≤ |Real.log n| + |Real.log s| := abs_sub _ _
        _ = _ := by rw [abs_of_nonneg hln]
    have hsq : (1 + |Real.log (n / s)|) ^ 2 ≤ 3 * (1 + Real.log n ^ 2 + Real.log s ^ 2) := by
      have h0 : 0 ≤ 1 + |Real.log (n / s)| := by positivity
      have : (1 + |Real.log (n / s)|) ^ 2 ≤ (1 + Real.log n + |Real.log s|) ^ 2 :=
        pow_le_pow_left₀ h0 (by linarith) 2
      have h3 : (1 + Real.log n + |Real.log s|) ^ 2 ≤ 3 * (1 + Real.log n ^ 2 + |Real.log s| ^ 2) := by
        nlinarith [sq_nonneg (1 - Real.log n), sq_nonneg (1 - |Real.log s|),
          sq_nonneg (Real.log n - |Real.log s|)]
      rw [sq_abs] at h3
      linarith
    have : ψ s * |mv - μ * (Real.log (n / s) + Real.log a) / a| ≤ B1 ψ n C1 μ s := by
      calc _ ≤ ψ s * (C1 * μ ^ 2 * (1 + |Real.log (n / s)|) ^ 2) :=
            mul_le_mul_of_nonneg_left hw hφ0.le
        _ ≤ ψ s * (C1 * μ ^ 2 * (3 * (1 + Real.log n ^ 2 + Real.log s ^ 2))) := by gcongr
        _ = B1 ψ n C1 μ s := by unfold B1; ring
    linarith
  · -- large s
    have hm := hmA h2
    rw [hm, zero_sub, abs_neg]
    have hlog : Real.log (n / s) + Real.log a = -(Real.log s - Real.log (n * a)) := by
      rw [hlam, Real.log_mul (by linarith) ha.ne']; ring
    rw [hlog]
    rcases eq_or_lt_of_le h2 with h3 | h3
    · have h0 : Real.log s - Real.log (n * a) = 0 := by rw [h3, sub_self]
      rw [h0, neg_zero, mul_zero, zero_div, abs_zero, mul_zero]
      linarith
    · have hB3v : B3 ψ n a μ s = μ / a * (ψ s * (Real.log s - Real.log (n * a))) := by
        unfold B3; rw [indicator_of_mem (show s ∈ Ioi (n * a) from h3)]
      have hpos : 0 ≤ Real.log s - Real.log (n * a) := by
        have := Real.log_le_log (by positivity) h3.le; linarith
      rw [abs_div, abs_of_pos ha, abs_mul, abs_of_pos hμ0, abs_neg, abs_of_nonneg hpos]
      have : ψ s * (μ * (Real.log s - Real.log (n * a)) / a) = B3 ψ n a μ s := by
        rw [hB3v]; field_simp
      linarith


lemma B_integrals {ψ : ℝ → ℝ} {k K : ℝ} (hD : Dens ψ k K) (n C1 μ K1 σ a : ℝ) (hn : 1 ≤ n)
    (ha : 0 < a) (hμ : 0 ≤ μ) :
    IntegrableOn (fun s => B1 ψ n C1 μ s + W01.B2 n K1 σ s + B3 ψ n a μ s) (Ioi 0) ∧
    ∫ s in Ioi (0:ℝ), (B1 ψ n C1 μ s + W01.B2 n K1 σ s + B3 ψ n a μ s) ≤
      3 * C1 * μ ^ 2 * (1 + Real.log n ^ 2 + ∫ s in Ioi (0:ℝ), ψ s * Real.log s ^ 2) +
      5 * K1 * (1 + Real.log n) * σ ^ (1 / 4 : ℝ) * (∫ s in Ioc (0:ℝ) 1, s ^ (-1 / 2 : ℝ)) +
      μ / a * (K / (n * a)) := by
  have hA : 0 < n * a := by nlinarith
  obtain ⟨hti, htb⟩ := tail_bound hD (n * a) hA
  have iB1 : IntegrableOn (B1 ψ n C1 μ) (Ioi 0) := by
    have := ((hD.int.mul_const (1 + Real.log n ^ 2)).add hD.l2int).const_mul (3 * C1 * μ ^ 2)
    exact IntegrableOn.congr_fun this (fun s _ => by simp [B1]) measurableSet_Ioi
  have ig2 : IntegrableOn (fun s : ℝ => 5 * K1 * (1 + Real.log n) * σ ^ (1 / 4 : ℝ) *
      s ^ (-1 / 2 : ℝ)) (Ioc 0 1) := W01.rpow_half_integrable.const_mul _
  have iB2 : IntegrableOn (W01.B2 n K1 σ) (Ioi 0) :=
    (ig2.integrable_indicator measurableSet_Ioc).integrableOn
  have ig3 : IntegrableOn (fun s => μ / a * (ψ s * (Real.log s - Real.log (n * a))))
      (Ioi (n * a)) := hti.const_mul _
  have iB3 : IntegrableOn (B3 ψ n a μ) (Ioi 0) :=
    (ig3.integrable_indicator measurableSet_Ioi).integrableOn
  refine ⟨(iB1.add iB2).add iB3, ?_⟩
  rw [integral_add (f := fun s => B1 ψ n C1 μ s + W01.B2 n K1 σ s) (g := B3 ψ n a μ)
    (iB1.add iB2) iB3, integral_add (f := B1 ψ n C1 μ) (g := W01.B2 n K1 σ) iB1 iB2]
  have e1 : ∫ s in Ioi (0:ℝ), B1 ψ n C1 μ s =
      3 * C1 * μ ^ 2 * (1 + Real.log n ^ 2 + ∫ s in Ioi (0:ℝ), ψ s * Real.log s ^ 2) := by
    unfold B1
    rw [integral_const_mul, integral_add (hD.int.mul_const _) hD.l2int,
      integral_mul_const, hD.int1]
    ring
  have e2 : ∫ s in Ioi (0:ℝ), W01.B2 n K1 σ s =
      5 * K1 * (1 + Real.log n) * σ ^ (1 / 4 : ℝ) * (∫ s in Ioc (0:ℝ) 1, s ^ (-1 / 2 : ℝ)) := by
    unfold W01.B2
    rw [setIntegral_indicator measurableSet_Ioc,
      inter_eq_right.2 Ioc_subset_Ioi_self, integral_const_mul]
  have e3 : ∫ s in Ioi (0:ℝ), B3 ψ n a μ s ≤ μ / a * (K / (n * a)) := by
    unfold B3
    rw [setIntegral_indicator measurableSet_Ioi, inter_eq_right.2 (Ioi_subset_Ioi hA.le),
      integral_const_mul]
    exact mul_le_mul_of_nonneg_left htb (by positivity)
  rw [e1, e2]
  linarith

lemma final_arith (n ln μ C1 M2 K1 q I0 a K : ℝ) (hn : 3 ≤ n) (hln : 1 ≤ ln) (hμ0 : 0 ≤ μ)
    (hμ : μ ≤ 2 / n) (hC1 : 0 ≤ C1) (hM2 : 0 ≤ M2) (hK1 : 0 ≤ K1) (hq0 : 0 ≤ q)
    (hq : q ≤ 1 / n ^ 2) (hI0 : 0 ≤ I0) (ha : 0 < a) (hK : 0 < K) :
    3 * C1 * μ ^ 2 * (1 + ln ^ 2 + M2) + 5 * K1 * (1 + ln) * q * I0 + μ / a * (K / (n * a)) ≤
      (12 * C1 * (2 + M2) + 10 * K1 * I0 + 2 * K / a ^ 2) * ln ^ 2 / n ^ 2 := by
  have hn0 : 0 < n := by linarith
  set t := ln ^ 2 with ht
  have ht1 : 1 ≤ t := by nlinarith
  have hlt : ln ≤ t := by nlinarith
  have hμ2 : μ ^ 2 ≤ 4 / n ^ 2 := by
    calc μ ^ 2 ≤ (2 / n) ^ 2 := pow_le_pow_left₀ hμ0 hμ 2
      _ = 4 / n ^ 2 := by ring
  have T1 : 3 * C1 * μ ^ 2 * (1 + t + M2) ≤ 12 * C1 * (2 + M2) * t / n ^ 2 := by
    have h1 : 1 + t + M2 ≤ (2 + M2) * t := by nlinarith
    calc 3 * C1 * μ ^ 2 * (1 + t + M2) ≤ 3 * C1 * (4 / n ^ 2) * ((2 + M2) * t) := by
          gcongr
      _ = 12 * C1 * (2 + M2) * t / n ^ 2 := by ring
  have T2 : 5 * K1 * (1 + ln) * q * I0 ≤ 10 * K1 * I0 * t / n ^ 2 := by
    calc 5 * K1 * (1 + ln) * q * I0 ≤ 5 * K1 * (2 * t) * (1 / n ^ 2) * I0 := by
          gcongr; linarith
      _ = 10 * K1 * I0 * t / n ^ 2 := by ring
  have T3 : μ / a * (K / (n * a)) ≤ 2 * K / a ^ 2 * t / n ^ 2 := by
    calc μ / a * (K / (n * a)) ≤ (2 / n) / a * (K / (n * a)) := by gcongr
      _ = 2 * K / a ^ 2 * 1 / n ^ 2 := by field_simp
      _ ≤ 2 * K / a ^ 2 * t / n ^ 2 := by gcongr
  calc _ ≤ 12 * C1 * (2 + M2) * t / n ^ 2 + 10 * K1 * I0 * t / n ^ 2 + 2 * K / a ^ 2 * t / n ^ 2 := by
        linarith
    _ = _ := by ring

lemma dens_level_integrable {ψ : ℝ → ℝ} {k K : ℝ} (hD : Dens ψ k K) (g : ℝ → ℝ) :
    IntegrableOn (fun s => ψ s * W01.levelMeasure g s) (Ioi 0) := by
  refine hD.int.mono' ?_ ?_
  · exact (hD.cont.measurable.mul
      (W01.levelMeasure_antitone g).measurable).aestronglyMeasurable
  · refine (ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun s hs => ?_)
    have h1 := hD.pos s hs
    obtain ⟨h2, h3⟩ := W01.levelMeasure_le_one g s
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
    nlinarith

lemma main_gen {ψ : ℝ → ℝ} {k K : ℝ} (hD : Dens ψ k K) (f : ℝ → ℝ) (hf : W01.Hyp f) :
    ∃ C : ℝ, ∃ N : ℕ, ∀ n : ℕ, N ≤ n →
      |(∫ s in Ioi (0:ℝ), ψ s * W01.levelMeasure (W01.yd f n) s) -
          (Real.log n + Real.log (W01.D f 1) + k) / (W01.D f 1 * ((n : ℝ) - 1))|
        ≤ C * Real.log n ^ 2 / (n : ℝ) ^ 2 := by
  have hf' := hf
  obtain ⟨hC1, hmono0, h0, h1, -, -⟩ := hf'
  obtain ⟨η, θ₀, F₁, β, g2, hS⟩ := W01.setup' f hf
  obtain ⟨n₀, hn₀, hshape⟩ := W01.shape f η θ₀ F₁ β g2 hS
  have hS' := hS
  obtain ⟨ha, -, -, -, -, -, hθh, hθ1, haF, hβ, hf01, -⟩ := hS'
  set a := W01.D f 1 with ha_def
  have hF : 0 < F₁ := lt_of_lt_of_le ha haF
  have hθ0 : 0 < θ₀ := by linarith
  set KK := 2 / a * (1 + |Real.log F₁|) with hKK
  set C1 := 1 / a * ((β + 2 * (a + β) ^ 2) * KK ^ 2 + 2 * β / a * KK) with hC1def
  have hK0 : 0 ≤ KK := by positivity
  have hC1nn : 0 ≤ C1 := by
    apply mul_nonneg (by positivity)
    apply add_nonneg
    · exact mul_nonneg (add_nonneg hβ (by positivity)) (by positivity)
    · exact mul_nonneg (div_nonneg (by linarith) ha.le) hK0
  set K1 := 1 + (1 + |Real.log a|) / a with hK1
  have hK1nn : 0 ≤ K1 := by positivity
  set M2 := ∫ s in Ioi (0:ℝ), ψ s * Real.log s ^ 2 with hM2def
  have hM2 : 0 ≤ M2 := setIntegral_nonneg measurableSet_Ioi
    (fun s hs => mul_nonneg (hD.pos s hs).le (sq_nonneg _))
  set I0 := ∫ s in Ioc (0:ℝ) 1, s ^ (-1 / 2 : ℝ) with hI0def
  have hI0 : 0 ≤ I0 := setIntegral_nonneg measurableSet_Ioc
    (fun s hs => Real.rpow_nonneg hs.1.le _)
  have hθabs : |θ₀| < 1 := by rw [abs_of_pos hθ0]; exact hθ1
  obtain ⟨N3, hN3⟩ := eventually_atTop.1
    ((tendsto_pow_const_mul_const_pow_of_abs_lt_one 9 hθabs).eventually
      (gt_mem_nhds (div_pos hθ0 hF)))
  refine ⟨12 * C1 * (2 + M2) + 10 * K1 * I0 + 2 * K / a ^ 2, max n₀ (max N3 3), fun n hn => ?_⟩
  have hnn₀ : n₀ ≤ n := le_trans (le_max_left _ _) hn
  have hnN3 : N3 ≤ n := le_trans (le_trans (le_max_left _ _) (le_max_right _ _)) hn
  have hn3 : 3 ≤ n := le_trans (le_trans (le_max_right _ _) (le_max_right _ _)) hn
  have hnR : (3:ℝ) ≤ n := by exact_mod_cast hn3
  have hn0 : (0:ℝ) < n := by linarith
  set μ := 1 / ((n : ℝ) - 1) with hμdef
  have hμ0 : 0 < μ := by apply div_pos one_pos; linarith
  have hμ1 : μ ≤ 1 := by rw [div_le_one (by linarith)]; linarith
  have hμ2 : μ ≤ 2 / n := by
    rw [hμdef, div_le_div_iff₀ (by linarith) hn0]; linarith
  have hln : 1 ≤ Real.log n := by
    rw [← Real.log_exp 1]
    apply Real.log_le_log (Real.exp_pos 1)
    have := Real.exp_one_lt_d9; linarith
  -- σ
  set σ := (n : ℝ) * F₁ * θ₀ ^ (n - 1) with hσdef
  have hσ0 : 0 < σ := by positivity
  have hσn : σ * (n : ℝ) ^ 8 < 1 := by
    have h9 := hN3 n hnN3
    have hp : θ₀ ^ n = θ₀ ^ (n - 1) * θ₀ := by
      rw [← pow_succ]; congr 1; omega
    rw [hp, lt_div_iff₀ hF] at h9
    have : (n : ℝ) ^ 9 * θ₀ ^ (n - 1) * F₁ < 1 := by
      have h' : (n : ℝ) ^ 9 * θ₀ ^ (n - 1) * F₁ * θ₀ < 1 * θ₀ := by nlinarith
      exact lt_of_mul_lt_mul_right h' hθ0.le
    calc σ * (n : ℝ) ^ 8 = (n : ℝ) ^ 9 * θ₀ ^ (n - 1) * F₁ := by rw [hσdef]; ring
      _ < 1 := this
  have hn8 : (1:ℝ) ≤ (n : ℝ) ^ 8 := one_le_pow₀ (by linarith)
  have hσ1 : σ ≤ 1 := by nlinarith
  have hσq : σ ^ (1 / 4 : ℝ) ≤ 1 / (n : ℝ) ^ 2 := by
    have hle : σ ≤ (1 / (n : ℝ) ^ 2) ^ 4 := by
      rw [div_pow, one_pow, ← pow_mul, le_div_iff₀ (by positivity)]; norm_num; linarith
    calc σ ^ (1 / 4 : ℝ) ≤ ((1 / (n : ℝ) ^ 2) ^ 4) ^ (1 / 4 : ℝ) :=
          Real.rpow_le_rpow hσ0.le hle (by norm_num)
      _ = 1 / (n : ℝ) ^ 2 := by
          have := Real.pow_rpow_inv_natCast (x := 1 / (n : ℝ) ^ 2) (n := 4) (by positivity)
            (by norm_num)
          have e : (1 / 4 : ℝ) = ((4 : ℕ) : ℝ)⁻¹ := by norm_num
          rw [e, this]
  have hσq0 : 0 ≤ σ ^ (1 / 4 : ℝ) := Real.rpow_nonneg hσ0.le _
  obtain ⟨hTi, hT⟩ := T_eq hD (n : ℝ) a (by linarith) ha
  have hmi := dens_level_integrable hD (W01.yd f n)
  rw [hT, ← integral_sub hmi hTi]
  obtain ⟨hBi, hBint⟩ := B_integrals hD (n : ℝ) C1 μ K1 σ a (by linarith) ha hμ0.le
  obtain ⟨-, -, -, -, -, hmA, -, -⟩ := hshape n hnn₀
  have hbound : ‖∫ s in Ioi (0:ℝ), (ψ s * W01.levelMeasure (W01.yd f n) s -
      ψ s * (1 / ((n : ℝ) - 1) * (Real.log (n / s) + Real.log a) / a))‖ ≤
      ∫ s in Ioi (0:ℝ), (B1 ψ n C1 μ s + W01.B2 n K1 σ s + B3 ψ n a μ s) := by
    apply norm_integral_le_of_norm_le hBi
    refine (ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun s hs => ?_)
    rw [← mul_sub]
    exact pt_bound hD (n : ℝ) a μ C1 σ (W01.levelMeasure (W01.yd f n) s) s hs ha hμ0 hμ1 (by linarith)
      hC1nn hσ0 hσ1 (W01.levelMeasure_le_one _ s).1 (W01.levelMeasure_le_one _ s).2 (hmA s)
      (fun h1 h2 => W01.window f η θ₀ F₁ β g2 hS n₀ hn₀ hshape n hnn₀ s h1 h2)
  rw [Real.norm_eq_abs] at hbound
  refine le_trans hbound (le_trans hBint ?_)
  exact final_arith (n : ℝ) (Real.log n) μ C1 M2 K1 (σ ^ (1 / 4 : ℝ)) I0 a K hnR hln hμ0.le hμ2
    hC1nn hM2 hK1nn hσq0 hσq hI0 ha hD.Kpos

end Y03


namespace Y03
open MeasureTheory Set Filter Topology

lemma k1_pi_half : k1 (Real.pi / 2) = 1 - Real.log 2 := by
  unfold k1
  have e : Real.pi / 2 / 2 = Real.pi / 4 := by ring
  rw [e, Real.cot_eq_cos_div_sin, Real.cos_pi_div_four, Real.sin_pi_div_four]
  have h2 : (0:ℝ) < Real.sqrt 2 := by positivity
  have hs : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  rw [div_self (by positivity), div_pow, hs]
  norm_num
  rw [one_div, Real.log_inv]; ring

lemma fpow_facts (f : ℝ → ℝ) (hf : W01.Hyp f) (n : ℕ) (hn : 1 ≤ n) :
    ContDiffOn ℝ 1 (fun t => f t ^ n) (Icc 0 1) ∧ MonotoneOn (fun t => f t ^ n) (Icc 0 1) ∧
      (fun t => f t ^ n) 1 - (fun t => f t ^ n) 0 = 1 := by
  obtain ⟨hC1, hmono0, h0, h1, -, -⟩ := hf
  refine ⟨hC1.pow n, ?_, ?_⟩
  · intro x hx y hy hxy
    have hx0 : 0 ≤ f x := by
      have := hmono0 (left_mem_Icc.2 zero_le_one) hx hx.1; rwa [h0] at this
    exact pow_le_pow_left₀ hx0 (hmono0 hx hy hxy) n
  · simp only [h0, h1, one_pow]; rw [zero_pow (by omega)]; ring

lemma p1 : stmt_1 := by
  intro θ h0 h1
  have hcc := cos_sq_lt θ h0 h1
  exact ⟨fun s hs => ⟨phi_eq _ s hcc, phi_pos _ s hcc, phi_le _ s hcc hs,
    phi_integrableOn _ s hcc, phi_integral _ s hcc⟩,
    fun s hs => phi_le_inv _ s hcc hs⟩

lemma p2 : stmt_2 := by
  intro θ h0 h1
  have hcc := cos_sq_lt θ h0 h1
  exact ⟨phi_log_integrable _ hcc, phi_log2_integrable _ hcc,
    by rw [phi_log_integral _ hcc, k1_eq θ h0 h1]⟩

lemma p3 : stmt_3 := by
  intro θ h0 h1 u v hu hv huv y hy hm hy1
  have hcc := cos_sq_lt θ h0 h1
  obtain ⟨e1, e2⟩ := layer _ hcc y hy hm hy1
  have hn : ∀ x, ‖u + W01.D y x • v‖ =
      Real.sqrt (1 + 2 * Real.cos θ * W01.D y x + W01.D y x ^ 2) :=
    fun x => norm_uv u v _ hu hv huv _
  simp_rw [hn]
  exact ⟨level_integrable _ hcc _, e1, e1.trans e2⟩

lemma p4 : stmt_4 := by
  refine ⟨?_, k1_pi_half, ?_⟩
  · funext s; simp [φθ, W01.φ, Real.cos_pi_div_two]
  · intro u v hu hv huv f n
    unfold L W01.L
    congr 1; funext x
    rw [norm_uv u v 0 hu hv huv]
    ring_nf

lemma p5 : stmt_5 := by
  intro θ h0 h1 u v hu hv huv f hf n hn
  have hcc := cos_sq_lt θ h0 h1
  have hc1 : 0 < 1 - Real.cos θ := by nlinarith
  obtain ⟨hy, hm, hy1⟩ := fpow_facts f hf n hn
  obtain ⟨hI, -, hE⟩ := p3 θ h0 h1 u v hu hv huv _ hy hm hy1
  have hcongr : (fun s => φθ (Real.cos θ) s / (1 - Real.cos θ) *
      W01.levelMeasure (W01.yd f n) s) = fun s => (1 - Real.cos θ)⁻¹ *
      (φθ (Real.cos θ) s * W01.levelMeasure (W01.D (fun t => f t ^ n)) s) := by
    funext s; unfold W01.yd; field_simp
  refine ⟨?_, ?_⟩
  · rw [hcongr]; exact hI.const_mul _
  · rw [hcongr, integral_const_mul, ← mul_assoc, mul_inv_cancel₀ hc1.ne', one_mul]
    exact hE

lemma p6 : stmt_6 := by
  intro θ h0 h1 u v hu hv huv f hf
  have hcc := cos_sq_lt θ h0 h1
  set c := Real.cos θ with hcdef
  have hc1 : 0 < 1 - c := by nlinarith
  have hD : Dens (fun s => φθ c s / (1 - c)) (k1 θ) (Kc c / (1 - c)) := by
    obtain ⟨hl, -, hlv⟩ := p2 θ h0 h1
    refine ⟨(phi_continuous c hcc).div_const _, fun s _ => div_pos (phi_pos c s hcc) hc1,
      fun s hs => ?_, div_pos (Kc_pos c) hc1, fun s hs => ?_,
      (phi_integrableOn c 0 hcc).div_const _, ?_, ?_,
      IntegrableOn.congr_fun ((phi_log2_integrable c hcc).div_const (1 - c)) (fun s _ => by ring) measurableSet_Ioi,
      ?_⟩
    · rw [div_le_one hc1]; exact phi_le c s hcc hs.le
    · rw [div_right_comm (Kc c)]
      exact div_le_div_of_nonneg_right (phi_le_K c s hcc hs) hc1.le
    · rw [integral_div, phi_integral c 0 hcc]
      simp only [mul_zero, add_zero, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true,
        zero_pow, Real.sqrt_one]
      field_simp; ring
    · exact IntegrableOn.congr_fun (hl.div_const (1 - c)) (fun s _ => by ring) measurableSet_Ioi
    · have : (fun s => φθ c s / (1 - c) * Real.log s) =
          fun s => (φθ c s * Real.log s) / (1 - c) := by funext s; ring
      rw [this, integral_div, hlv, ← hcdef]
      field_simp
  obtain ⟨C, N, hCN⟩ := main_gen hD f hf
  refine ⟨(1 - c) * C, max N 1, fun n hn => ?_⟩
  have hnN : N ≤ n := le_trans (le_max_left _ _) hn
  have hn1 : 1 ≤ n := le_trans (le_max_right _ _) hn
  obtain ⟨-, h5⟩ := p5 θ h0 h1 u v hu hv huv f hf n hn1
  have := hCN n hnN
  have key : 2 - L u v f n - (1 - c) * (Real.log n + Real.log (W01.D f 1) + k1 θ) /
      (W01.D f 1 * ((n : ℝ) - 1)) = (1 - c) * ((∫ s in Ioi (0:ℝ), φθ c s / (1 - c) *
        W01.levelMeasure (W01.yd f n) s) - (Real.log n + Real.log (W01.D f 1) + k1 θ) /
      (W01.D f 1 * ((n : ℝ) - 1))) := by rw [h5]; ring
  rw [key, abs_mul, abs_of_pos hc1]
  calc _ ≤ (1 - c) * (C * Real.log n ^ 2 / (n : ℝ) ^ 2) := mul_le_mul_of_nonneg_left this hc1.le
    _ = _ := by ring

end Y03

theorem Y03.proof_1 : Y03.stmt_1 := Y03.p1
theorem Y03.proof_2 : Y03.stmt_2 := Y03.p2
theorem Y03.proof_3 : Y03.stmt_3 := Y03.p3
theorem Y03.proof_4 : Y03.stmt_4 := Y03.p4
theorem Y03.proof_5 : Y03.stmt_5 := Y03.p5
theorem Y03.proof_6 : Y03.stmt_6 := Y03.p6
