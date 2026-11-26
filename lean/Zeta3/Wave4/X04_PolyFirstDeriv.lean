import Mathlib
import Zeta3.Wave2.U03_R1Positive
/-!
# X04 — `P'(1) = P(1) R_1` for the smoothed regular `2m`-gon (paper label `thm:poly`)

Target id: X04.
Sources: `paper/ieee/main.tex`, Subsection "Polygons at the Diamond End", lines 1645–1655;
detailed proof `paper/ieee/proofs/sec_g4_poly.tex` (label `g4:thm:poly`, Lemmas
`g4:lem:polygon`, `g4:lem:edge`, `g4:lem:arc`, `g4:lem:smooth`, `g4:lem:first`, and
"Proof of Theorem 23").

Setup quoted verbatim (main.tex line 1645):
"The diamond is the case $m=2$ of a regular $2m$-gon. Smooth it the same way: with
$n_i=(\cos\frac{i\pi}m,\sin\frac{i\pi}m)$, let $K_p=\{\sum_{i=1}^m|\langle n_i,x\rangle|^p\le1\}$.
At $p=1$ this is a regular $2m$-gon of perimeter $P(1)=4m(\sec\frac{\pi}{2m}-\cos\frac{\pi}{2m})$,
and $m=2$ is the squircle."

Theorem quoted verbatim (main.tex lines 1647–1653):
"\begin{theorem}[First derivative at the regular polygon$^\dagger$]\label{thm:poly}
$P'(1)=P(1)\,R_1$, where $V=\cot\frac{\pi}{2m}$ and
\[
R_1=\frac1V\sum_{i=1}^{m-1}\sin\frac{i\pi}{m}\,\log\frac{V}{\sin(i\pi/m)} .
\]
In particular $P'(1)=0$ only for $m=2$, and for the hexagon ($m=3$) $P'(1)=2\sqrt3\log2$.
\end{theorem}"

Supplement (sec_g4_poly.tex, after Remark g4:rem:analytic's lead-in) quoted verbatim:
"For $p<1$ we call $P(p)$ the length of the closed curve $\partial K_p$ given by the formula of
Lemma~\ref{g4:lem:arc}. For $p\ge1$ this is the perimeter of the convex body $K_p$, so $P'(1)$ is
the right derivative, and it equals the derivative of the extension."

Lemma g4:lem:arc quoted verbatim:
"The map $x\mapsto x/F_p(x)^{1/p}$ is a homeomorphism of $\partial K_1$ onto $\partial K_p$ that
fixes every ray. The curve $x_p(t)=y(t)F(t,d)^{-1/p}$ runs over the arc of $\partial K_p$ between
the rays of $V_0$ and $V_1$, and $|x_p'(t)|=G(t,d):=F^{-1-1/p}\,R(F,D;t)$. Hence
$P(p)=\sum_{k=0}^{2m-1}\int_0^1G_k(t,d)\dd t$, and $G(t,0)=L$, so $P(1)=\sum_eL_e$."
with $R(F,D;t)=\sqrt{(LF-\sigma(t)D)^2+h^2D^2}$, $d=p-1$, $D=\sum_ja_j'(t)e^{d\lambda_j(t)}$.

Lemma g4:lem:smooth quoted verbatim:
"For every $k\ge0$ there is a constant $K_k$, depending on the polygon and the weights, such that
for $|d|\le\frac14$ and $0<t<1$
$|\pa_d^kG(t,d)|\le K_k\,t^{-1/4}(1-t)^{-1/4}\bigl(1+|\log t|+|\log(1-t)|\bigr)^k$."

Lemma g4:lem:first quoted verbatim:
"$P'(1)=\sum_e\Bigl(\ip{\tau_e}{V_0^e}\Sone(V_0^e)-\ip{\tau_e}{V_1^e}\Sone(V_1^e)\Bigr).$
If $\Sone(C_k)=s_0$ for all vertices, then $P'(1)=-s_0\,P(1)$."
with $\Sone(x)=\sum_ja_j(x)\log\frac{a_j(x)}{w_j}$, $0\log0=0$; in its proof
"$\pa_dG|_{d=0}=-L\Sone-\sigma\Sone'=-(\sigma\Sone)'$".

Proof of Theorem 23 quoted (excerpts): "the vertices are $C_k=\varrho^kC_0$ with
$C_0=\rho_0(0,1)$", "$V:=\sum_{i=1}^{m-1}\sin\frac{i\pi}m=\cot\frac\pi{2m}$", "$\rho_0=1/V=
\tan\frac\pi{2m}$", "$\Sone(C_k)=\Sone(C_0)=s_0$ for all $k$", "$s_0=\sum_{i=1}^{m-1}
\frac{\sin(i\pi/m)}V\log\frac{\sin(i\pi/m)}V=-R_1$".

## Formalization

* Points are `Plane := EuclideanSpace ℝ (Fin 2)` (Euclidean distance). `dot u x` is the
  Euclidean inner product written in coordinates. `nrm m i = (cos(iπ/m), sin(iπ/m))`.
* `Fp m p x = ∑_{i=1}^m |⟨n_i,x⟩|^p` (real `rpow` of a nonnegative base), `Kp m p = {Fp ≤ 1}`.
  All statements are for integers `m ≥ 2` (the paper's setting).
* `vertex m k = tan(π/2m)·(−sin(kπ/m), cos(kπ/m))` is `C_k = ϱ^k C_0`, `C_0 = ρ_0 (0,1)`,
  `ρ_0 = tan(π/2m)` (supplement, proof of Theorem 23). That these are the vertices of `K_1` is
  stmt_1, so nothing rests on this choice silently.
* `edgePt m k t = C_k + t (C_{k+1} − C_k)` is `y(t)`; `arc m p k t = Fp(y(t))^{-1/p} • y(t)`
  is the supplement's `x_p(t)`; `loop m p s = arc m p ⌊s⌋ (s − ⌊s⌋)` on `[0, 2m]` traverses the
  `2m` arcs one after another (at `s = 2m` it is `arc (2m) 0 = C_{2m}/… = C_0/…`).
* THE DEFINITION OF `P(p)`: `perim m p` is the true Euclidean length (total variation
  `eVariationOn`) of the closed curve `loop m p` on `[0, 2m]`, as a real number via `toReal`.
  This is defined for every real `p`, in particular for `p` near `1` on both sides, as the
  supplement requires ("For $p<1$ we call $P(p)$ the length of the closed curve
  $\partial K_p$"). stmt_3 says `loop m p` is a simple closed curve whose image is the boundary
  of `K_p` (= `{Fp = 1}`), so `perim` is the length of `∂K_p` (for `p ≥ 1` the perimeter of
  the convex body). Because `toReal ⊤ = 0`, finiteness is part of the statements that use
  `perim` (stmt_5 and stmt_13 assert the length is a finite real number given by an integral).
* `side`, `tau`, `sigma`, `apo` are `L`, `τ`, `σ(t) = ⟨y(t),τ⟩`, `h = √(|V_0|²−σ_0²)` of the
  supplement. `Fe m p k t = Fp m p (y t) = F(t,d)`. The supplement's
  `D = ∑ a_j' e^{dλ_j}` equals `(∂_t F)/p` (supplement: "$F'=\dots=pD$"); we write it as
  `deriv (Fe m p k) t / p`. On `0 < t < 1` every `|⟨n_j,y(t)⟩|` is positive, so `Fe` is
  differentiable there and `deriv` is the true derivative (no `deriv`-of-nondifferentiable
  convention is used: stmt_4 asserts the true derivative of the arc has norm `G`).
  `G m p k t = F^{-1-1/p}·√((L F − σ D)² + h² D²)`.
* `S1 m x = ∑_{i=1}^m a_i log a_i`, `a_i = |⟨n_i,x⟩|`. The paper sets `0 log 0 = 0`; Lean's
  `Real.log 0 = 0` gives the same value, and this is the paper's definition, not a repair.
* `R_1` and `P(1)` are `U03.R1`, `U03.P1` (Zeta3/Wave2/U03_R1Positive.lean, identical to Wave1
  `T06.R1`, `T06.P1`): the paper's formulas with `\sec = 1/cos`, sum over `Finset.Ico 1 m`.
* Range of `p` for the integral formula: the supplement works on `|d| ≤ 1/4`, i.e.
  `p ∈ [3/4, 5/4]`; stmt_4/5 use the open interval `(3/4, 5/4)`, stmt_7 the closed one.

Statements (in the order of the supplement's proof):
* stmt_1: `K_1` is the regular `2m`-gon: `K_1 = conv{C_0,…,C_{2m−1}}`, every `C_k` is on
  `{Fp = 1}`, and all sides have length `L = 2 tan(π/2m) sin(π/2m)`.
* stmt_2: the length of the actual `p = 1` curve is `P(1) = 4m(sec(π/2m) − cos(π/2m))`.
* stmt_3: for every `p > 0`, `loop m p` is a simple closed curve with image
  `frontier K_p = {Fp = 1}` (Lemma g4:lem:arc, first sentence).
* stmt_4: speed: for `p ∈ (3/4,5/4)`, `0<t<1`, `|x_p'(t)| = G(t,p−1)` (Lemma g4:lem:arc).
* stmt_5: `P(p) = ∑_{k<2m} ∫_0^1 G_k` with each `G_k` integrable, as an equality of lengths
  (finite) (Lemma g4:lem:arc, "Hence").
* stmt_6: `G(t,0) = L` (Lemma g4:lem:arc).
* stmt_7: domination of `∂_p G` for `|p−1| ≤ 1/4` (Lemma g4:lem:smooth, `k = 1`).
* stmt_8: pointwise first variation `∂_p G|_{p=1} = −L S1(y(t)) − σ(t) (S1∘y)'(t)`
  (proof of Lemma g4:lem:first).
* stmt_9: per-edge first variation: `d/dp ∫_0^1 G_k |_{p=1} = ⟨τ,C_k⟩S1(C_k) − ⟨τ,C_{k+1}⟩S1(C_{k+1})`
  (Lemma g4:lem:first).
* stmt_10: `∑_{i=1}^{m−1} sin(iπ/m) = cot(π/2m)` (proof of Theorem 23).
* stmt_11: `S1(C_k) = −R_1` for every vertex (proof of Theorem 23).
* stmt_12: `P(1) = ∑_k L_k` telescoping form used at the end of Lemma g4:lem:first:
  `∑_k (⟨τ_k,C_k⟩ − ⟨τ_k,C_{k+1}⟩) = −P(1)`.
* stmt_13: THE THEOREM: for every integer `m ≥ 2`, the arc lengths are finite for `p` near `1`
  and `HasDerivAt perim (P(1) R_1) 1`.
* stmt_14: "P'(1)=0 only for m=2": for every `m ≥ 2` the derivative exists and is `0` iff `m = 2`.
* stmt_15: hexagon: `HasDerivAt (perim 3) (2√3 log 2) 1`.

## NOT COVERED

* "$P$ is $C^\infty$ near $p=1$" and the bounds of Lemma g4:lem:smooth for `k ≥ 2`
  (only `k = 1` is stated, stmt_7, which is what thm:poly needs).
* "for $m=2$ it is not analytic there" (Theorem `thm:zero`; a different target).
* The general first-variation formula of Lemma g4:lem:first for arbitrary centrally symmetric
  polygons and weights (only the regular `2m`-gon, power convention, is stated), and the
  "general formula" paragraph (parallelogram, affine-regular hexagon).
* The notes' alternative proof by integration by parts in polar coordinates (main.tex line 1655).
* Convexity of `K_p` for `p ≥ 1` is not stated; "perimeter" is the length of `∂K_p`.
-/
-- BEGIN STATEMENTS
namespace X04

open MeasureTheory Set Filter Real
open scoped Topology

/-- The Euclidean plane. -/
abbrev Plane := EuclideanSpace ℝ (Fin 2)

/-- Euclidean inner product `⟨u, x⟩` in coordinates. -/
def dot (u x : Plane) : ℝ := u 0 * x 0 + u 1 * x 1

/-- The unit normal `n_i = (cos(iπ/m), sin(iπ/m))`. -/
noncomputable def nrm (m i : ℕ) : Plane :=
  !₂[Real.cos ((i : ℝ) * π / (m : ℝ)), Real.sin ((i : ℝ) * π / (m : ℝ))]

/-- `F_p(x) = ∑_{i=1}^m |⟨n_i, x⟩|^p`. -/
noncomputable def Fp (m : ℕ) (p : ℝ) (x : Plane) : ℝ :=
  ∑ i ∈ Finset.Icc 1 m, |dot (nrm m i) x| ^ p

/-- `K_p = {x : ∑_{i=1}^m |⟨n_i,x⟩|^p ≤ 1}`. -/
def Kp (m : ℕ) (p : ℝ) : Set Plane := {x | Fp m p x ≤ 1}

/-- The vertex `C_k = ϱ^k C_0`, `C_0 = tan(π/2m)·(0,1)`, `ϱ` the rotation by `π/m`. -/
noncomputable def vertex (m k : ℕ) : Plane :=
  Real.tan (π / (2 * (m : ℝ))) •
    !₂[-Real.sin ((k : ℝ) * π / (m : ℝ)), Real.cos ((k : ℝ) * π / (m : ℝ))]

/-- Side length `L_k = |C_{k+1} − C_k|`. -/
noncomputable def side (m k : ℕ) : ℝ := dist (vertex m (k + 1)) (vertex m k)

/-- Unit vector `τ` along the edge `e_k`, from `C_k` to `C_{k+1}`. -/
noncomputable def tau (m k : ℕ) : Plane := (side m k)⁻¹ • (vertex m (k + 1) - vertex m k)

/-- `y(t) = C_k + t (C_{k+1} − C_k)`, the straight edge of `K_1`. -/
noncomputable def edgePt (m k : ℕ) (t : ℝ) : Plane :=
  vertex m k + t • (vertex m (k + 1) - vertex m k)

/-- `σ(t) = ⟨y(t), τ⟩`. -/
noncomputable def sigma (m k : ℕ) (t : ℝ) : ℝ := dot (edgePt m k t) (tau m k)

/-- `h`: the distance from `0` to the line of the edge, `h² = |V_0|² − σ_0²`. -/
noncomputable def apo (m k : ℕ) : ℝ :=
  Real.sqrt (‖vertex m k‖ ^ 2 - (dot (vertex m k) (tau m k)) ^ 2)

/-- `F(t, d) = F_p(y(t))`, `d = p − 1`. -/
noncomputable def Fe (m : ℕ) (p : ℝ) (k : ℕ) (t : ℝ) : ℝ := Fp m p (edgePt m k t)

/-- The supplement's length element
`G(t,d) = F^{-1-1/p} √((L F − σ D)² + h² D²)` with `D = (∂_t F)/p`. -/
noncomputable def G (m : ℕ) (p : ℝ) (k : ℕ) (t : ℝ) : ℝ :=
  Fe m p k t ^ (-1 - 1 / p) *
    Real.sqrt ((side m k * Fe m p k t - sigma m k t * (deriv (Fe m p k) t / p)) ^ 2 +
      apo m k ^ 2 * (deriv (Fe m p k) t / p) ^ 2)

/-- The arc `x_p(t) = y(t) F_p(y(t))^{-1/p}` of `∂K_p` over the edge `e_k`. -/
noncomputable def arc (m : ℕ) (p : ℝ) (k : ℕ) (t : ℝ) : Plane :=
  (Fe m p k t ^ (-(1 / p))) • edgePt m k t

/-- The closed curve `∂K_p`, traversed on `[0, 2m]`: arc `k` on `[k, k+1]`. -/
noncomputable def loop (m : ℕ) (p s : ℝ) : Plane :=
  arc m p ⌊s⌋₊ (s - ⌊s⌋₊)

/-- `P(p)`: the Euclidean length of the closed curve `∂K_p` (total variation of `loop m p`
on `[0, 2m]`), as a real number. -/
noncomputable def perim (m : ℕ) (p : ℝ) : ℝ :=
  (eVariationOn (loop m p) (Icc 0 (2 * (m : ℝ)))).toReal

/-- `Σ₁(x) = ∑_{i=1}^m a_i log a_i`, `a_i = |⟨n_i,x⟩|` (`0 log 0 = 0`). -/
noncomputable def S1 (m : ℕ) (x : Plane) : ℝ :=
  ∑ i ∈ Finset.Icc 1 m, |dot (nrm m i) x| * Real.log |dot (nrm m i) x|

/-- main.tex 1645: "At $p=1$ this is a regular $2m$-gon": `K_1` is the convex hull of the
`2m` points `C_k`, which lie on `∂K_1`, and all sides have length `2 tan(π/2m) sin(π/2m)`. -/
def stmt_1 : Prop :=
  ∀ m : ℕ, 2 ≤ m →
    Kp m 1 = convexHull ℝ (vertex m '' {k | k < 2 * m}) ∧
    (∀ k : ℕ, Fp m 1 (vertex m k) = 1) ∧
    (∀ k : ℕ, side m k = 2 * Real.tan (π / (2 * (m : ℝ))) * Real.sin (π / (2 * (m : ℝ))))

/-- main.tex 1645: "of perimeter $P(1)=4m(\sec\frac{\pi}{2m}-\cos\frac{\pi}{2m})$", for the
length of the actual `p = 1` curve. -/
def stmt_2 : Prop :=
  ∀ m : ℕ, 2 ≤ m → perim m 1 = U03.P1 m

/-- Lemma g4:lem:arc (first sentences): for `p > 0`, the `2m` arcs, traversed by `loop m p`,
form a simple closed curve whose image is `∂K_p = {F_p = 1}`. -/
def stmt_3 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∀ p : ℝ, 0 < p →
    ContinuousOn (loop m p) (Icc 0 (2 * (m : ℝ))) ∧
    loop m p (2 * (m : ℝ)) = loop m p 0 ∧
    InjOn (loop m p) (Ico 0 (2 * (m : ℝ))) ∧
    loop m p '' Icc 0 (2 * (m : ℝ)) = frontier (Kp m p) ∧
    frontier (Kp m p) = {x | Fp m p x = 1}

/-- Lemma g4:lem:arc: "$|x_p'(t)|=G(t,d)$" for `|p − 1| < 1/4`, `0 < t < 1`. -/
def stmt_4 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∀ p ∈ Ioo (3 / 4 : ℝ) (5 / 4), ∀ k : ℕ, ∀ t ∈ Ioo (0 : ℝ) 1,
    ∃ v : Plane, HasDerivAt (arc m p k) v t ∧ ‖v‖ = G m p k t

/-- Lemma g4:lem:arc: "Hence $P(p)=\sum_{k=0}^{2m-1}\int_0^1G_k(t,d)\dd t$": each `G_k` is
integrable on `[0,1]` and the (finite) length of `∂K_p` is the sum of the integrals. -/
def stmt_5 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∀ p ∈ Ioo (3 / 4 : ℝ) (5 / 4),
    (∀ k : ℕ, IntervalIntegrable (G m p k) volume 0 1) ∧
    eVariationOn (loop m p) (Icc 0 (2 * (m : ℝ))) =
      ENNReal.ofReal (∑ k ∈ Finset.range (2 * m), ∫ t in (0 : ℝ)..1, G m p k t)

/-- Lemma g4:lem:arc: "$G(t,0)=L$". -/
def stmt_6 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∀ k : ℕ, ∀ t ∈ Ioo (0 : ℝ) 1, G m 1 k t = side m k

/-- Lemma g4:lem:smooth, `k = 1`: for `|d| ≤ 1/4`, `0 < t < 1`, `G` is differentiable in `p`
and `|∂_d G| ≤ K t^{-1/4}(1-t)^{-1/4}(1+|log t|+|log(1-t)|)`. -/
def stmt_7 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∃ K : ℝ, ∀ k : ℕ, ∀ p ∈ Icc (3 / 4 : ℝ) (5 / 4), ∀ t ∈ Ioo (0 : ℝ) 1,
    ∃ g : ℝ, HasDerivAt (fun q => G m q k t) g p ∧
      |g| ≤ K * t ^ (-(1 / 4 : ℝ)) * (1 - t) ^ (-(1 / 4 : ℝ)) *
        (1 + |Real.log t| + |Real.log (1 - t)|)

/-- Proof of Lemma g4:lem:first: "$\pa_dG|_{d=0}=-L\Sone-\sigma\Sone'$" pointwise on `(0,1)`. -/
def stmt_8 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∀ k : ℕ, ∀ t ∈ Ioo (0 : ℝ) 1,
    HasDerivAt (fun p => G m p k t)
      (-(side m k * S1 m (edgePt m k t) +
          sigma m k t * deriv (fun s => S1 m (edgePt m k s)) t)) 1

/-- Lemma g4:lem:first, one edge: the length `∫_0^1 G_k` of the arc over `e_k` has
`p`-derivative `⟨τ,C_k⟩ Σ₁(C_k) − ⟨τ,C_{k+1}⟩ Σ₁(C_{k+1})` at `p = 1`. -/
def stmt_9 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∀ k : ℕ,
    HasDerivAt (fun p => ∫ t in (0 : ℝ)..1, G m p k t)
      (dot (tau m k) (vertex m k) * S1 m (vertex m k) -
        dot (tau m k) (vertex m (k + 1)) * S1 m (vertex m (k + 1))) 1

/-- Proof of Theorem 23: "$V:=\sum_{i=1}^{m-1}\sin\frac{i\pi}m=\cot\frac\pi{2m}$". -/
def stmt_10 : Prop :=
  ∀ m : ℕ, 2 ≤ m →
    ∑ i ∈ Finset.Ico 1 m, Real.sin ((i : ℝ) * π / (m : ℝ)) = Real.cot (π / (2 * (m : ℝ)))

/-- Proof of Theorem 23: "$\Sone(C_k)=\Sone(C_0)=s_0$ for all $k$" and "$s_0=-R_1$". -/
def stmt_11 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∀ k : ℕ, S1 m (vertex m k) = -U03.R1 m

/-- End of Lemma g4:lem:first: "$\sum_k(\ip{\tau_k}{C_k}-\ip{\tau_k}{C_{k+1}})=-\sum_kL_k=
-P(1)$", for the regular `2m`-gon with `P(1) = 4m(sec − cos)`. -/
def stmt_12 : Prop :=
  ∀ m : ℕ, 2 ≤ m →
    ∑ k ∈ Finset.range (2 * m),
      (dot (tau m k) (vertex m k) - dot (tau m k) (vertex m (k + 1))) = -U03.P1 m

/-- Theorem thm:poly: "$P'(1)=P(1)\,R_1$" for every integer `m ≥ 2`, for the true length
`P(p)` of `∂K_p` (finite for `p` near `1`). -/
def stmt_13 : Prop :=
  ∀ m : ℕ, 2 ≤ m →
    (∀ᶠ p in 𝓝 (1 : ℝ), eVariationOn (loop m p) (Icc 0 (2 * (m : ℝ))) ≠ ⊤) ∧
    HasDerivAt (perim m) (U03.P1 m * U03.R1 m) 1

/-- Theorem thm:poly: "In particular $P'(1)=0$ only for $m=2$": for every `m ≥ 2`, `P` is
differentiable at `1` and its derivative there vanishes iff `m = 2`. -/
def stmt_14 : Prop :=
  ∀ m : ℕ, 2 ≤ m → ∃ D : ℝ, HasDerivAt (perim m) D 1 ∧ (D = 0 ↔ m = 2)

/-- Theorem thm:poly: "for the hexagon ($m=3$) $P'(1)=2\sqrt3\log2$". -/
def stmt_15 : Prop :=
  HasDerivAt (perim 3) (2 * Real.sqrt 3 * Real.log 2) 1

end X04
-- END STATEMENTS

namespace X04
open MeasureTheory Set Filter Real

lemma sin_sum_aux (h : ℝ) (n : ℕ) :
    ∑ i ∈ Finset.range n, 2 * Real.sin h * Real.sin (2 * (i : ℝ) * h) =
      Real.cos h - Real.cos ((2 * (n : ℝ) - 1) * h) := by
  induction n with
  | zero =>
    simp only [Finset.range_zero, Finset.sum_empty, CharP.cast_eq_zero]
    rw [show (2*(0:ℝ)-1)*h = -h by ring, Real.cos_neg]; ring
  | succ n ih =>
    rw [Finset.sum_range_succ, ih, Real.two_mul_sin_mul_sin,
      show h - 2*(n:ℝ)*h = -((2*n-1)*h) by ring, Real.cos_neg,
      show h + 2*(n:ℝ)*h = (2*((n:ℝ)+1)-1)*h by ring]
    push_cast; ring

lemma sum_sin_eq_cot (m : ℕ) (hm : 2 ≤ m) :
    ∑ i ∈ Finset.Ico 1 m, Real.sin ((i : ℝ) * π / (m : ℝ)) = Real.cot (π / (2 * (m : ℝ))) := by
  have hm' : (2:ℝ) ≤ m := by exact_mod_cast hm
  set h := π / (2 * (m : ℝ)) with hh
  have hpos : 0 < h := by positivity
  have hlt : h < π := by
    rw [hh, div_lt_iff₀ (by positivity)]; nlinarith [Real.pi_pos]
  have hs : 0 < Real.sin h := Real.sin_pos_of_pos_of_lt_pi hpos hlt
  have key := sin_sum_aux h m
  have hterm : ∀ i : ℕ, Real.sin ((i : ℝ) * π / (m : ℝ)) = Real.sin (2 * (i:ℝ) * h) := by
    intro i; congr 1; rw [hh]; field_simp
  have hend : (2 * (m : ℝ) - 1) * h = π - h := by
    rw [hh]; field_simp
  rw [hend, Real.cos_pi_sub] at key
  have hIco : ∑ i ∈ Finset.Ico 1 m, Real.sin ((i : ℝ) * π / (m : ℝ)) =
      ∑ i ∈ Finset.range m, Real.sin (2 * (i:ℝ) * h) := by
    rw [Finset.sum_congr rfl (fun i _ => hterm i), Finset.range_eq_Ico,
      Finset.sum_eq_sum_Ico_succ_bot (show 0 < m by omega)]
    simp
  rw [hIco, Real.cot_eq_cos_div_sin, eq_div_iff hs.ne']
  rw [← Finset.mul_sum] at key
  linarith

end X04

namespace X04
open MeasureTheory Set Filter Real

lemma vertex0 (m k : ℕ) : vertex m k 0 =
    Real.tan (π / (2 * (m : ℝ))) * -Real.sin ((k : ℝ) * π / (m : ℝ)) := by
  simp [vertex]

lemma vertex1 (m k : ℕ) : vertex m k 1 =
    Real.tan (π / (2 * (m : ℝ))) * Real.cos ((k : ℝ) * π / (m : ℝ)) := by
  simp [vertex]

lemma nrm0 (m i : ℕ) : nrm m i 0 = Real.cos ((i : ℝ) * π / (m : ℝ)) := by simp [nrm]
lemma nrm1 (m i : ℕ) : nrm m i 1 = Real.sin ((i : ℝ) * π / (m : ℝ)) := by simp [nrm]

lemma dot_nrm_vertex (m i k : ℕ) : dot (nrm m i) (vertex m k) =
    Real.tan (π / (2 * (m : ℝ))) * Real.sin (((i : ℝ) - k) * π / (m : ℝ)) := by
  rw [dot, nrm0, nrm1, vertex0, vertex1, sub_mul, sub_div, Real.sin_sub]
  ring

/-- shifting sum for an `m`-periodic function. -/
lemma shift_sum (m : ℕ) (hm : 1 ≤ m) (f : ℝ → ℝ) (hf : ∀ x, f (x + m) = f x) (k : ℕ) :
    ∑ i ∈ Finset.Icc 1 m, f ((i : ℝ) - k) = ∑ i ∈ Finset.Icc 1 m, f (i : ℝ) := by
  have hT : ∀ k : ℕ, ∑ i ∈ Finset.Icc 1 m, f ((i : ℝ) - k) =
      ∑ i ∈ Finset.range m, f (1 + (i : ℝ) - k) := by
    intro k
    rw [← Finset.Ico_add_one_right_eq_Icc, Finset.sum_Ico_eq_sum_range]
    rw [Nat.add_sub_cancel]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    push_cast; ring_nf
  induction k with
  | zero => simp
  | succ k ih =>
    rw [← ih, hT, hT]
    obtain ⟨n, rfl⟩ : ∃ n, m = n + 1 := ⟨m - 1, by omega⟩
    rw [Finset.sum_range_succ', Finset.sum_range_succ]
    congr 1
    · refine Finset.sum_congr rfl (fun i _ => ?_)
      push_cast; ring_nf
    · have := hf (-(k:ℝ))
      push_cast at this ⊢
      rw [show (1:ℝ)+0-(k+1) = -k by ring, show (1:ℝ)+n-k = -k + (n+1) by ring, this]

lemma per_abs_sin (m : ℕ) (hm : 1 ≤ m) (φ : ℝ → ℝ) (x : ℝ) :
    φ |Real.sin ((x + m) * π / m)| = φ |Real.sin (x * π / m)| := by
  have hm' : (m:ℝ) ≠ 0 := by positivity
  rw [show (x + m) * π / m = x * π / m + π by field_simp, Real.sin_add_pi, abs_neg]

/-- the sum over a full period of `φ(|sin(jπ/m)|)`. -/
lemma period_sum (m : ℕ) (hm : 1 ≤ m) (φ : ℝ → ℝ) (hφ : φ 0 = 0) (k : ℕ) :
    ∑ i ∈ Finset.Icc 1 m, φ |Real.sin (((i : ℝ) - k) * π / m)| =
      ∑ j ∈ Finset.Ico 1 m, φ (Real.sin ((j : ℝ) * π / m)) := by
  have := shift_sum m hm (fun x => φ |Real.sin (x * π / m)|) (per_abs_sin m hm φ) k
  rw [this, ← Finset.Ico_add_one_right_eq_Icc, Finset.sum_Ico_succ_top hm]
  have hm' : (m:ℝ) ≠ 0 := by positivity
  rw [show (m:ℝ) * π / m = π by field_simp, Real.sin_pi, abs_zero, hφ, add_zero]
  refine Finset.sum_congr rfl (fun j hj => ?_)
  rw [Finset.mem_Ico] at hj
  rw [abs_of_nonneg]
  apply Real.sin_nonneg_of_nonneg_of_le_pi (by positivity)
  rw [div_le_iff₀ (by positivity)]
  have : (j:ℝ) ≤ m := by exact_mod_cast hj.2.le
  nlinarith [Real.pi_pos]

end X04

namespace X04
open MeasureTheory Set Filter Real

lemma ang_facts (m : ℕ) (hm : 2 ≤ m) :
    0 < π / (2 * (m : ℝ)) ∧ π / (2 * (m : ℝ)) ≤ π / 4 := by
  have hm' : (2:ℝ) ≤ m := by exact_mod_cast hm
  refine ⟨by positivity, ?_⟩
  rw [div_le_div_iff₀ (by positivity) (by norm_num)]; nlinarith [Real.pi_pos]

lemma rho_pos (m : ℕ) (hm : 2 ≤ m) : 0 < Real.tan (π / (2 * (m : ℝ))) := by
  obtain ⟨h1, h2⟩ := ang_facts m hm
  exact Real.tan_pos_of_pos_of_lt_pi_div_two h1 (by linarith [Real.pi_pos])

lemma tan_mul_cot (m : ℕ) (hm : 2 ≤ m) :
    Real.tan (π / (2 * (m : ℝ))) * Real.cot (π / (2 * (m : ℝ))) = 1 := by
  obtain ⟨h1, h2⟩ := ang_facts m hm
  have hs : 0 < Real.sin (π / (2 * (m : ℝ))) :=
    Real.sin_pos_of_pos_of_lt_pi h1 (by linarith [Real.pi_pos])
  have hc : 0 < Real.cos (π / (2 * (m : ℝ))) :=
    Real.cos_pos_of_mem_Ioo ⟨by linarith, by linarith [Real.pi_pos]⟩
  rw [Real.tan_eq_sin_div_cos, Real.cot_eq_cos_div_sin]
  field_simp

lemma abs_dot_vertex (m i k : ℕ) (hm : 2 ≤ m) : |dot (nrm m i) (vertex m k)| =
    Real.tan (π / (2 * (m : ℝ))) * |Real.sin (((i : ℝ) - k) * π / (m : ℝ))| := by
  rw [dot_nrm_vertex, abs_mul, abs_of_pos (rho_pos m hm)]

lemma Fp_one_vertex (m : ℕ) (hm : 2 ≤ m) (k : ℕ) : Fp m 1 (vertex m k) = 1 := by
  unfold Fp
  simp only [Real.rpow_one]
  rw [Finset.sum_congr rfl (fun i _ => abs_dot_vertex m i k hm)]
  rw [period_sum m (by omega) (fun a => Real.tan (π / (2 * (m : ℝ))) * a) (by simp) k,
    ← Finset.mul_sum, sum_sin_eq_cot m hm, tan_mul_cot m hm]

lemma S1_vertex (m : ℕ) (hm : 2 ≤ m) (k : ℕ) : S1 m (vertex m k) = -U03.R1 m := by
  unfold S1
  set ρ := Real.tan (π / (2 * (m : ℝ))) with hρ
  rw [Finset.sum_congr rfl (fun i _ => by rw [abs_dot_vertex m i k hm])]
  rw [period_sum m (by omega) (fun a => ρ * a * Real.log (ρ * a)) (by simp) k]
  unfold U03.R1 U03.V
  have hV : 1 / Real.cot (π / (2 * (m : ℝ))) = ρ := by
    rw [hρ, ← tan_mul_cot m hm]; have := tan_mul_cot m hm
    have hc : Real.cot (π / (2 * (m : ℝ))) ≠ 0 := by
      intro h; rw [h] at this; simp at this
    field_simp
  rw [hV, Finset.mul_sum, ← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  have : Real.cot (π / (2 * (m : ℝ))) / Real.sin ((j : ℝ) * π / m) =
      (ρ * Real.sin ((j : ℝ) * π / m))⁻¹ := by
    rw [← hV]; field_simp
  rw [this, Real.log_inv]; ring

end X04

namespace X04
open MeasureTheory Set Filter Real

lemma dist_plane (v w : Plane) : dist v w = Real.sqrt ((v 0 - w 0) ^ 2 + (v 1 - w 1) ^ 2) := by
  rw [EuclideanSpace.dist_eq, Fin.sum_univ_two, Real.dist_eq, Real.dist_eq, sq_abs, sq_abs]

lemma norm_plane (v : Plane) : ‖v‖ = Real.sqrt (v 0 ^ 2 + v 1 ^ 2) := by
  rw [EuclideanSpace.norm_eq, Fin.sum_univ_two, Real.norm_eq_abs, Real.norm_eq_abs, sq_abs, sq_abs]

lemma dot_self (u : Plane) : dot u u = ‖u‖ ^ 2 := by
  rw [norm_plane, Real.sq_sqrt (by positivity), dot]; ring

lemma side_eq (m : ℕ) (hm : 2 ≤ m) (k : ℕ) :
    side m k = 2 * Real.tan (π / (2 * (m : ℝ))) * Real.sin (π / (2 * (m : ℝ))) := by
  obtain ⟨h1, h2⟩ := ang_facts m hm
  have hρ := rho_pos m hm
  have hs : 0 < Real.sin (π / (2 * (m : ℝ))) :=
    Real.sin_pos_of_pos_of_lt_pi h1 (by linarith [Real.pi_pos])
  unfold side
  rw [dist_plane, vertex0, vertex0, vertex1, vertex1]
  set ρ := Real.tan (π / (2 * (m : ℝ)))
  set y := π / (2 * (m : ℝ)) with hy
  have hm2 : (2:ℝ) ≤ m := by exact_mod_cast hm
  have hm' : (m:ℝ) ≠ 0 := by positivity
  have hb : ((k + 1 : ℕ) : ℝ) * π / m = (k : ℝ) * π / m + 2 * y := by
    rw [hy]; push_cast; field_simp
  rw [hb, Real.sin_add, Real.cos_add]
  set S := Real.sin ((k : ℝ) * π / m)
  set C := Real.cos ((k : ℝ) * π / m)
  have hSC : S ^ 2 + C ^ 2 = 1 := Real.sin_sq_add_cos_sq _
  have h2 : Real.sin (2 * y) ^ 2 + Real.cos (2 * y) ^ 2 = 1 := Real.sin_sq_add_cos_sq _
  have hc2 : Real.cos (2 * y) = 1 - 2 * Real.sin y ^ 2 := by
    rw [Real.cos_two_mul, Real.cos_sq']; ring
  have : (ρ * -(S * Real.cos (2 * y) + C * Real.sin (2 * y)) - ρ * -S) ^ 2 +
      (ρ * (C * Real.cos (2 * y) - S * Real.sin (2 * y)) - ρ * C) ^ 2 =
      (2 * ρ * Real.sin y) ^ 2 := by
    have e : (ρ * -(S * Real.cos (2 * y) + C * Real.sin (2 * y)) - ρ * -S) ^ 2 +
      (ρ * (C * Real.cos (2 * y) - S * Real.sin (2 * y)) - ρ * C) ^ 2 =
      ρ ^ 2 * (S ^ 2 + C ^ 2) * ((Real.cos (2 * y) - 1) ^ 2 + Real.sin (2 * y) ^ 2) := by ring
    rw [e, hSC]
    have : (Real.cos (2 * y) - 1) ^ 2 + Real.sin (2 * y) ^ 2 = 2 - 2 * Real.cos (2 * y) := by
      nlinarith [h2]
    rw [this, hc2]; ring
  rw [this, Real.sqrt_sq (by positivity)]

lemma side_pos (m : ℕ) (hm : 2 ≤ m) (k : ℕ) : 0 < side m k := by
  obtain ⟨h1, h2⟩ := ang_facts m hm
  have hs : 0 < Real.sin (π / (2 * (m : ℝ))) :=
    Real.sin_pos_of_pos_of_lt_pi h1 (by linarith [Real.pi_pos])
  rw [side_eq m hm k]; have := rho_pos m hm; positivity

lemma dot_tau_diff (m : ℕ) (hm : 2 ≤ m) (k : ℕ) :
    dot (tau m k) (vertex m k) - dot (tau m k) (vertex m (k + 1)) = -side m k := by
  have hL := side_pos m hm k
  have hn : ‖vertex m (k + 1) - vertex m k‖ = side m k := by
    rw [side, dist_eq_norm]
  have hd := dot_self (vertex m (k + 1) - vertex m k)
  rw [hn] at hd
  unfold tau
  simp only [dot, PiLp.smul_apply, PiLp.sub_apply, smul_eq_mul] at hd ⊢
  field_simp
  nlinarith [hd]

lemma P1_eq (m : ℕ) (hm : 2 ≤ m) :
    U03.P1 m = 2 * (m : ℝ) * (2 * Real.tan (π / (2 * (m : ℝ))) * Real.sin (π / (2 * (m : ℝ)))) := by
  obtain ⟨h1, h2⟩ := ang_facts m hm
  have hc : 0 < Real.cos (π / (2 * (m : ℝ))) :=
    Real.cos_pos_of_mem_Ioo ⟨by linarith, by linarith [Real.pi_pos]⟩
  unfold U03.P1
  set x := π / (2 * (m : ℝ))
  have := Real.sin_sq_add_cos_sq x
  have e : 1 / Real.cos x - Real.cos x = Real.sin x ^ 2 / Real.cos x := by
    rw [eq_div_iff hc.ne', sub_mul, one_div_mul_cancel hc.ne']; nlinarith
  rw [e, Real.tan_eq_sin_div_cos]; ring

end X04

theorem X04.proof_10 : X04.stmt_10 := X04.sum_sin_eq_cot
theorem X04.proof_11 : X04.stmt_11 := X04.S1_vertex
theorem X04.proof_12 : X04.stmt_12 := by
  intro m hm
  rw [Finset.sum_congr rfl (fun k _ => X04.dot_tau_diff m hm k)]
  simp only [X04.side_eq m hm, Finset.sum_const, Finset.card_range, nsmul_eq_mul]
  rw [X04.P1_eq m hm]; push_cast; ring

namespace X04
open MeasureTheory Set Filter Real
open scoped Topology

lemma sign_pair (m : ℕ) (hm : 2 ≤ m) (j : ℤ) :
    0 ≤ Real.sin ((j : ℝ) * π / m) * Real.sin (((j : ℝ) - 1) * π / m) ∧
    (Real.sin ((j : ℝ) * π / m) ≠ 0 ∨ Real.sin (((j : ℝ) - 1) * π / m) ≠ 0) := by
  have hm2 : (2:ℝ) ≤ m := by exact_mod_cast hm
  have hm0 : (m:ℝ) ≠ 0 := by positivity
  have hmz : (0:ℤ) < m := by exact_mod_cast (by omega : 0 < m)
  set q := j / (m:ℤ)
  set r := j % (m:ℤ)
  have hr0 : 0 ≤ r := Int.emod_nonneg _ hmz.ne'
  have hrm : r < m := Int.emod_lt_of_pos _ hmz
  have hj : (j:ℝ) = (r:ℝ) + (q:ℝ) * m := by
    have := Int.emod_add_mul_ediv j m
    have h' : ((j:ℤ):ℝ) = ((r + m * q : ℤ):ℝ) := by rw [this]
    push_cast at h'; rw [h']; ring
  have e1 : (j:ℝ) * π / m = (r:ℝ) * π / m + q * π := by rw [hj]; field_simp
  have e2 : ((j:ℝ) - 1) * π / m = ((r:ℝ) - 1) * π / m + q * π := by rw [hj]; field_simp; ring
  rw [e1, e2, Real.sin_add_int_mul_pi, Real.sin_add_int_mul_pi]
  have hq : ((-1:ℝ) ^ q) ^ 2 = 1 := by
    rw [← zpow_natCast, ← zpow_mul, mul_comm, zpow_mul]; simp
  have hrR0 : (0:ℝ) ≤ r := by exact_mod_cast hr0
  have hrRm : (r:ℝ) + 1 ≤ m := by exact_mod_cast hrm
  have hpi := Real.pi_pos
  have hq0 : ((-1:ℝ) ^ q) ≠ 0 := zpow_ne_zero _ (by norm_num)
  rcases (show r = 0 ∨ 1 ≤ r by omega) with h | h
  · rw [h]
    simp only [Int.cast_zero, zero_mul, zero_div, Real.sin_zero, mul_zero, zero_mul, le_refl,
      true_and, ne_eq, not_true_eq_false, false_or]
    rw [zero_sub, neg_mul, neg_div, Real.sin_neg]
    have : 0 < Real.sin (1 * π / m) := by
      apply Real.sin_pos_of_pos_of_lt_pi (by positivity)
      rw [div_lt_iff₀ (by positivity)]; nlinarith
    exact mul_ne_zero hq0 (by linarith)
  · have hrR1 : (1:ℝ) ≤ r := by exact_mod_cast h
    have ha : 0 < Real.sin ((r:ℝ) * π / m) := by
      apply Real.sin_pos_of_pos_of_lt_pi (by positivity)
      rw [div_lt_iff₀ (by positivity)]; nlinarith
    have hb : 0 ≤ Real.sin (((r:ℝ) - 1) * π / m) := by
      apply Real.sin_nonneg_of_nonneg_of_le_pi
      · apply div_nonneg (mul_nonneg (by linarith) hpi.le) (by positivity)
      · rw [div_le_iff₀ (by positivity)]; nlinarith
    refine ⟨?_, Or.inl (mul_ne_zero hq0 ha.ne')⟩
    have : (-1:ℝ) ^ q * Real.sin ((r:ℝ) * π / m) * ((-1:ℝ) ^ q * Real.sin (((r:ℝ) - 1) * π / m))
      = ((-1:ℝ) ^ q) ^ 2 * (Real.sin ((r:ℝ) * π / m) * Real.sin (((r:ℝ) - 1) * π / m)) := by ring
    rw [this, hq, one_mul]; positivity

lemma abs_comb {a b t : ℝ} (hab : 0 ≤ a * b) (t0 : 0 ≤ t) (t1 : t ≤ 1) :
    |(1 - t) * a + t * b| = (1 - t) * |a| + t * |b| := by
  rcases le_total 0 a with ha | ha <;> rcases le_total 0 b with hb | hb
  · rw [abs_of_nonneg ha, abs_of_nonneg hb, abs_of_nonneg (by nlinarith)]
  · have : a * b = 0 := le_antisymm (mul_nonpos_of_nonneg_of_nonpos ha hb) hab
    rcases mul_eq_zero.1 this with h | h <;> subst h <;>
      simp [abs_mul, abs_of_nonneg t0, abs_of_nonneg (by linarith : (0:ℝ) ≤ 1 - t)]
  · have : a * b = 0 := le_antisymm (mul_nonpos_of_nonpos_of_nonneg ha hb) hab
    rcases mul_eq_zero.1 this with h | h <;> subst h <;>
      simp [abs_mul, abs_of_nonneg t0, abs_of_nonneg (by linarith : (0:ℝ) ≤ 1 - t)]
  · rw [abs_of_nonpos ha, abs_of_nonpos hb, abs_of_nonpos (by nlinarith)]; ring

lemma dot_edge (m i k : ℕ) (t : ℝ) : dot (nrm m i) (edgePt m k t) =
    (1 - t) * dot (nrm m i) (vertex m k) + t * dot (nrm m i) (vertex m (k + 1)) := by
  simp only [dot, edgePt, PiLp.add_apply, PiLp.smul_apply, PiLp.sub_apply, smul_eq_mul]; ring

lemma dot_vertex_prod (m : ℕ) (hm : 2 ≤ m) (i k : ℕ) :
    0 ≤ dot (nrm m i) (vertex m k) * dot (nrm m i) (vertex m (k + 1)) ∧
    (dot (nrm m i) (vertex m k) ≠ 0 ∨ dot (nrm m i) (vertex m (k + 1)) ≠ 0) := by
  have hρ := rho_pos m hm
  obtain ⟨h1, h2⟩ := sign_pair m hm ((i:ℤ) - k)
  rw [dot_nrm_vertex, dot_nrm_vertex]
  have e : ((i:ℝ) - ((k + 1 : ℕ) : ℝ)) = (((i:ℤ) - k : ℤ) : ℝ) - 1 := by push_cast; ring
  have e' : ((i:ℝ) - k) = (((i:ℤ) - k : ℤ) : ℝ) := by push_cast; ring
  rw [e, e']
  refine ⟨?_, ?_⟩
  · have : Real.tan (π / (2 * (m:ℝ))) * Real.sin ((((i:ℤ) - k : ℤ) : ℝ) * π / m) *
        (Real.tan (π / (2 * (m:ℝ))) * Real.sin (((((i:ℤ) - k : ℤ) : ℝ) - 1) * π / m)) =
        Real.tan (π / (2 * (m:ℝ))) ^ 2 * (Real.sin ((((i:ℤ) - k : ℤ) : ℝ) * π / m) *
          Real.sin (((((i:ℤ) - k : ℤ) : ℝ) - 1) * π / m)) := by ring
    rw [this]; positivity
  · rcases h2 with h | h
    · exact Or.inl (mul_ne_zero hρ.ne' h)
    · exact Or.inr (mul_ne_zero hρ.ne' h)

lemma abs_dot_edge (m : ℕ) (hm : 2 ≤ m) (i k : ℕ) {t : ℝ} (t0 : 0 ≤ t) (t1 : t ≤ 1) :
    |dot (nrm m i) (edgePt m k t)| =
      (1 - t) * |dot (nrm m i) (vertex m k)| + t * |dot (nrm m i) (vertex m (k + 1))| := by
  rw [dot_edge]; exact abs_comb (dot_vertex_prod m hm i k).1 t0 t1

lemma abs_dot_edge_pos (m : ℕ) (hm : 2 ≤ m) (i k : ℕ) {t : ℝ} (t0 : 0 < t) (t1 : t < 1) :
    0 < |dot (nrm m i) (edgePt m k t)| := by
  rw [abs_dot_edge m hm i k t0.le t1.le]
  rcases (dot_vertex_prod m hm i k).2 with h | h
  · have := abs_pos.2 h
    have : 0 ≤ t * |dot (nrm m i) (vertex m (k + 1))| := by positivity
    nlinarith
  · have := abs_pos.2 h
    have : 0 ≤ (1 - t) * |dot (nrm m i) (vertex m k)| :=
      mul_nonneg (by linarith) (abs_nonneg _)
    nlinarith

lemma Fe_one (m : ℕ) (hm : 2 ≤ m) (k : ℕ) {t : ℝ} (t0 : 0 ≤ t) (t1 : t ≤ 1) :
    Fe m 1 k t = 1 := by
  have h0 := Fp_one_vertex m hm k
  have h1 := Fp_one_vertex m hm (k + 1)
  unfold Fe Fp at *
  simp only [Real.rpow_one] at *
  rw [Finset.sum_congr rfl (fun i _ => abs_dot_edge m hm i k t0 t1), Finset.sum_add_distrib,
    ← Finset.mul_sum, ← Finset.mul_sum, h0, h1]; ring

lemma deriv_Fe_one (m : ℕ) (hm : 2 ≤ m) (k : ℕ) {t : ℝ} (t0 : 0 < t) (t1 : t < 1) :
    deriv (Fe m 1 k) t = 0 := by
  have : Fe m 1 k =ᶠ[𝓝 t] fun _ => (1:ℝ) := by
    filter_upwards [Ioo_mem_nhds t0 t1] with s hs using Fe_one m hm k hs.1.le hs.2.le
  rw [this.deriv_eq]; simp

end X04

theorem X04.proof_6 : X04.stmt_6 := by
  intro m hm k t ht
  unfold X04.G
  rw [X04.deriv_Fe_one m hm k ht.1 ht.2, X04.Fe_one m hm k ht.1.le ht.2.le]
  norm_num
  rw [Real.sqrt_sq (X04.side_pos m hm k).le]

namespace X04
open MeasureTheory Set Filter Real
open scoped Topology

lemma var_affine (A v : Plane) (a : ℝ) :
    eVariationOn (fun s => A + (s - a) • v) (Icc a (a + 1)) = ENNReal.ofReal ‖v‖ := by
  apply le_antisymm
  · have hL : LipschitzWith ‖v‖₊ (fun s : ℝ => A + (s - a) • v) := by
      apply LipschitzWith.of_dist_le_mul
      intro x y
      rw [dist_eq_norm, dist_eq_norm]
      have : A + (x - a) • v - (A + (y - a) • v) = (x - y) • v := by
        rw [sub_smul, sub_smul, sub_smul]; abel
      rw [this, norm_smul, mul_comm]; simp
    have := hL.lipschitzOnWith (s := univ) |>.comp_eVariationOn_le (g := id)
      (s := Icc a (a + 1)) (mapsTo_univ _ _)
    simp only [Function.comp_def, id, eVariationOn_id_Icc] at this
    convert this using 1
    simp only [add_sub_cancel_left, ENNReal.ofReal_one, mul_one]
    rw [ofReal_norm]; rfl
  · have := eVariationOn.edist_le (fun s => A + (s - a) • v) (s := Icc a (a + 1))
      (left_mem_Icc.2 (by linarith)) (right_mem_Icc.2 (by linarith))
    convert this using 1
    rw [eDist_eq_eNorm_sub]
    simp [ofReal_norm]

lemma var_sum {E : Type*} [PseudoEMetricSpace E] (f : ℝ → E) (n : ℕ) :
    eVariationOn f (Icc 0 (n : ℝ)) =
      ∑ k ∈ Finset.range n, eVariationOn f (Icc (k : ℝ) (k + 1)) := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [Finset.sum_range_succ, ← ih]
    have := eVariationOn.Icc_add_Icc f (s := univ) (a := 0) (b := n) (c := n + 1)
      (by positivity) (by linarith) (mem_univ _)
    simp only [univ_inter] at this
    rw [this]; push_cast; rfl

lemma loop_one_eq (m : ℕ) (hm : 2 ≤ m) (k : ℕ) :
    EqOn (loop m 1) (fun s => vertex m k + (s - k) • (vertex m (k + 1) - vertex m k))
      (Icc (k : ℝ) (k + 1)) := by
  intro s hs
  rcases eq_or_lt_of_le hs.2 with h | h
  · subst h
    have : ⌊(k : ℝ) + 1⌋₊ = k + 1 := by exact_mod_cast Nat.floor_natCast (k + 1)
    simp only [loop, this, arc]
    push_cast
    rw [sub_self, Fe_one m hm (k + 1) le_rfl zero_le_one]
    simp [edgePt]
  · have hfl : ⌊s⌋₊ = k := by
      rw [Nat.floor_eq_iff (by linarith [hs.1, (Nat.cast_nonneg k : (0:ℝ) ≤ k)])]
      exact ⟨hs.1, h⟩
    have hF : Fe m 1 k (s - k) = 1 := Fe_one m hm k (by linarith [hs.1]) (by linarith)
    simp only [loop, hfl, arc, hF, Real.one_rpow, one_smul, edgePt]

lemma perim_one (m : ℕ) (hm : 2 ≤ m) : perim m 1 = U03.P1 m := by
  unfold perim
  have h2m : (2 * (m : ℝ)) = ((2 * m : ℕ) : ℝ) := by push_cast; ring
  rw [h2m, var_sum]
  rw [Finset.sum_congr rfl (fun k _ => by
    rw [eVariationOn.congr (loop_one_eq m hm k), var_affine])]
  rw [ENNReal.toReal_sum (fun k _ => ENNReal.ofReal_ne_top)]
  simp only [ENNReal.toReal_ofReal (norm_nonneg _)]
  have : ∀ k : ℕ, ‖vertex m (k + 1) - vertex m k‖ = side m k := fun k => by
    rw [side, dist_eq_norm]
  simp only [this, side_eq m hm, Finset.sum_const, Finset.card_range, nsmul_eq_mul]
  rw [P1_eq m hm]; push_cast; ring

end X04

theorem X04.proof_2 : X04.stmt_2 := X04.perim_one

namespace X04
open MeasureTheory Set Filter Real
open scoped Topology

/-- `a_i(s) = ⟨n_i, y(s)⟩` on edge `k`. -/
noncomputable def aE (m k i : ℕ) (s : ℝ) : ℝ := dot (nrm m i) (edgePt m k s)

/-- derivative of `|a_i|` at an interior point. -/
noncomputable def cE (m k i : ℕ) (t : ℝ) : ℝ :=
  (SignType.sign (aE m k i t) : ℝ) * (dot (nrm m i) (vertex m (k + 1)) - dot (nrm m i) (vertex m k))

lemma aE_hasDeriv (m k i : ℕ) (t : ℝ) :
    HasDerivAt (aE m k i) (dot (nrm m i) (vertex m (k + 1)) - dot (nrm m i) (vertex m k)) t := by
  have e : aE m k i = fun s => dot (nrm m i) (vertex m k) +
      s * (dot (nrm m i) (vertex m (k + 1)) - dot (nrm m i) (vertex m k)) := by
    funext s; rw [aE, dot_edge]; ring
  rw [e]
  simpa using ((hasDerivAt_id' t).mul_const
    (dot (nrm m i) (vertex m (k + 1)) - dot (nrm m i) (vertex m k))).const_add
      (dot (nrm m i) (vertex m k))

lemma aE_ne (m : ℕ) (hm : 2 ≤ m) (k i : ℕ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) : aE m k i t ≠ 0 :=
  abs_pos.1 (abs_dot_edge_pos m hm i k ht.1 ht.2)

lemma absE_hasDeriv (m : ℕ) (hm : 2 ≤ m) (k i : ℕ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    HasDerivAt (fun s => |aE m k i s|) (cE m k i t) t := by
  have := (hasDerivAt_abs (aE_ne m hm k i ht)).comp t (aE_hasDeriv m k i t)
  simpa [cE, Function.comp_def, mul_comm] using this

lemma Fe_eq (m : ℕ) (q : ℝ) (k : ℕ) :
    Fe m q k = fun s => ∑ i ∈ Finset.Icc 1 m, |aE m k i s| ^ q := rfl

lemma Fe_hasDeriv (m : ℕ) (hm : 2 ≤ m) (q : ℝ) (k : ℕ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    HasDerivAt (Fe m q k)
      (∑ i ∈ Finset.Icc 1 m, cE m k i t * q * |aE m k i t| ^ (q - 1)) t := by
  rw [Fe_eq]
  exact HasDerivAt.fun_sum (fun i _ => (absE_hasDeriv m hm k i ht).rpow_const
    (Or.inl (abs_pos.2 (aE_ne m hm k i ht)).ne'))

lemma sum_cE (m : ℕ) (hm : 2 ≤ m) (k : ℕ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    ∑ i ∈ Finset.Icc 1 m, cE m k i t = 0 := by
  have h := (Fe_hasDeriv m hm 1 k ht).deriv
  rw [deriv_Fe_one m hm k ht.1 ht.2] at h
  simpa using h.symm

lemma S1_edge_hasDeriv (m : ℕ) (hm : 2 ≤ m) (k : ℕ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    HasDerivAt (fun s => S1 m (edgePt m k s))
      (∑ i ∈ Finset.Icc 1 m, cE m k i t * Real.log |aE m k i t|) t := by
  have e : (fun s => S1 m (edgePt m k s)) =
      fun s => ∑ i ∈ Finset.Icc 1 m, |aE m k i s| * Real.log |aE m k i s| := rfl
  rw [e]
  have h := HasDerivAt.fun_sum (u := Finset.Icc 1 m) (fun i _ =>
    (absE_hasDeriv m hm k i ht).mul ((absE_hasDeriv m hm k i ht).log
      (abs_pos.2 (aE_ne m hm k i ht)).ne'))
  convert h using 1
  have : ∀ i ∈ Finset.Icc 1 m, cE m k i t * Real.log |aE m k i t| + |aE m k i t| *
      (cE m k i t / |aE m k i t|) = cE m k i t * Real.log |aE m k i t| + cE m k i t := by
    intro i _
    rw [mul_div_cancel₀ _ (abs_pos.2 (aE_ne m hm k i ht)).ne']
  rw [Finset.sum_congr rfl this, Finset.sum_add_distrib, sum_cE m hm k ht, add_zero]

lemma stmt8_core (m : ℕ) (hm : 2 ≤ m) (k : ℕ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    HasDerivAt (fun p => G m p k t)
      (-(side m k * S1 m (edgePt m k t) +
          sigma m k t * deriv (fun s => S1 m (edgePt m k s)) t)) 1 := by
  set I := Finset.Icc 1 m
  set A : ℕ → ℝ := fun i => |aE m k i t| with hAdef
  set c : ℕ → ℝ := fun i => cE m k i t with hcdef
  set L := side m k
  set σ := sigma m k t
  set h := apo m k
  have hA : ∀ i, 0 < A i := fun i => abs_pos.2 (aE_ne m hm k i ht)
  have hL : 0 < L := side_pos m hm k
  have hG : (fun p => G m p k t) = fun q : ℝ =>
      (∑ i ∈ I, A i ^ q) ^ (-1 - 1 / q) *
        Real.sqrt ((L * ∑ i ∈ I, A i ^ q - σ * ((∑ i ∈ I, c i * q * A i ^ (q - 1)) / q)) ^ 2 +
          h ^ 2 * ((∑ i ∈ I, c i * q * A i ^ (q - 1)) / q) ^ 2) := by
    funext q
    unfold G
    rw [(Fe_hasDeriv m hm q k ht).deriv]
    rfl
  have hF1 : ∑ i ∈ I, A i ^ (1:ℝ) = 1 := by
    have := Fe_one m hm k ht.1.le ht.2.le
    rw [Fe_eq] at this; simpa using this
  have hc0 : ∑ i ∈ I, c i = 0 := sum_cE m hm k ht
  have hS : S1 m (edgePt m k t) = ∑ i ∈ I, A i * Real.log (A i) := rfl
  have hS' : deriv (fun s => S1 m (edgePt m k s)) t = ∑ i ∈ I, c i * Real.log (A i) :=
    (S1_edge_hasDeriv m hm k ht).deriv
  -- derivatives
  have dF : HasDerivAt (fun q : ℝ => ∑ i ∈ I, A i ^ q) (∑ i ∈ I, A i ^ (1:ℝ) * Real.log (A i)) (1:ℝ) :=
    HasDerivAt.fun_sum (fun i _ => (hasStrictDerivAt_const_rpow (hA i) 1).hasDerivAt)
  have dN : HasDerivAt (fun q : ℝ => ∑ i ∈ I, c i * q * A i ^ (q - 1))
      (∑ i ∈ I, (c i * 1 * A i ^ ((1:ℝ) - 1) +
        c i * 1 * (A i ^ ((1:ℝ) - 1) * Real.log (A i) * 1))) (1:ℝ) := by
    apply HasDerivAt.fun_sum
    intro i _
    have h1 : HasDerivAt (fun q : ℝ => c i * q) (c i * 1) 1 := (hasDerivAt_id 1).const_mul (c i)
    have h2 : HasDerivAt (fun q : ℝ => A i ^ (q - 1)) (A i ^ ((1:ℝ) - 1) * Real.log (A i) * 1) 1 := by
      have := ((hasDerivAt_id' (1:ℝ)).sub_const 1).const_rpow (hA i)
      convert this using 1; ring
    exact h1.mul h2
  have hF1' : ∑ i ∈ I, A i = 1 := by simpa using hF1
  set S := ∑ i ∈ I, A i * Real.log (A i)
  set Dp := ∑ i ∈ I, c i * Real.log (A i)
  have dF' : HasDerivAt (fun q : ℝ => ∑ i ∈ I, A i ^ q) S 1 := by
    convert dF using 1; simp [S]
  have dD : HasDerivAt (fun q : ℝ => (∑ i ∈ I, c i * q * A i ^ (q - 1)) / q) Dp 1 := by
    convert dN.fun_div (hasDerivAt_id' (1:ℝ)) one_ne_zero using 1
    simp only [sub_self, Real.rpow_zero, mul_one, one_mul, one_pow, div_one]
    rw [Finset.sum_add_distrib, hc0]; simp [Dp]
  have hD1 : (∑ i ∈ I, c i * (1:ℝ) * A i ^ ((1:ℝ) - 1)) / 1 = 0 := by simp [hc0]
  have dE : HasDerivAt (fun q : ℝ => -1 - 1 / q) 1 1 := by
    have := (hasDerivAt_inv (x := (1:ℝ)) one_ne_zero).const_sub (-1)
    convert this using 1
    · funext q; simp [one_div]
    · simp
  have dFE : HasDerivAt (fun q : ℝ => (∑ i ∈ I, A i ^ q) ^ (-1 - 1 / q)) (-2 * S) 1 := by
    convert dF'.rpow dE (by rw [hF1]; norm_num) using 1
    rw [hF1]; simp; ring
  have dX : HasDerivAt (fun q : ℝ =>
      (L * ∑ i ∈ I, A i ^ q - σ * ((∑ i ∈ I, c i * q * A i ^ (q - 1)) / q)) ^ 2 +
        h ^ 2 * ((∑ i ∈ I, c i * q * A i ^ (q - 1)) / q) ^ 2) (2 * L * (L * S - σ * Dp)) 1 := by
    convert (((dF'.const_mul L).sub (dD.const_mul σ)).fun_pow 2).add
      ((dD.fun_pow 2).const_mul (h ^ 2)) using 1
    rw [hD1]; simp only [Pi.sub_apply]; rw [hF1, hD1]; push_cast; ring
  have hX1 : (L * ∑ i ∈ I, A i ^ (1:ℝ) - σ * ((∑ i ∈ I, c i * (1:ℝ) * A i ^ ((1:ℝ) - 1)) / 1)) ^ 2
      + h ^ 2 * ((∑ i ∈ I, c i * (1:ℝ) * A i ^ ((1:ℝ) - 1)) / 1) ^ 2 = L ^ 2 := by
    rw [hF1, hD1]; ring
  have dSq := dX.sqrt (by rw [hX1]; positivity)
  have dP := dFE.mul dSq
  rw [hG]
  convert dP using 1
  rw [hX1, Real.sqrt_sq hL.le, hS, hS', hF1]
  field_simp
  ring

end X04

theorem X04.proof_8 : X04.stmt_8 := fun m hm k _ ht => X04.stmt8_core m hm k ht

namespace X04
open MeasureTheory Set Filter Real
open scoped Topology

lemma Fe_pos (m : ℕ) (hm : 2 ≤ m) (p : ℝ) (k : ℕ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    0 < Fe m p k t := by
  rw [Fe_eq]
  apply Finset.sum_pos (fun i _ => Real.rpow_pos_of_pos (abs_pos.2 (aE_ne m hm k i ht)) _)
  exact ⟨1, Finset.mem_Icc.2 ⟨le_rfl, by omega⟩⟩

lemma edgePt_hasDeriv (m k : ℕ) (t : ℝ) :
    HasDerivAt (edgePt m k) (vertex m (k + 1) - vertex m k) t := by
  have := ((hasDerivAt_id' t).smul_const (vertex m (k + 1) - vertex m k)).const_add (vertex m k)
  have e : edgePt m k = fun s => vertex m k + s • (vertex m (k + 1) - vertex m k) := rfl
  rw [e]; simpa using this

lemma speed_identity (C0 C1 τ0 τ1 t F D L : ℝ) (hτ : τ0 ^ 2 + τ1 ^ 2 = 1) :
    (F * (L * τ0) - D * (C0 + t * (L * τ0))) ^ 2 + (F * (L * τ1) - D * (C1 + t * (L * τ1))) ^ 2 =
      (L * F - ((C0 + t * (L * τ0)) * τ0 + (C1 + t * (L * τ1)) * τ1) * D) ^ 2 +
        ((C0 ^ 2 + C1 ^ 2) - (C0 * τ0 + C1 * τ1) ^ 2) * D ^ 2 := by
  linear_combination (L ^ 2 * F ^ 2 - D ^ 2 * (2 * t * L * (C0 * τ0 + C1 * τ1) +
    t ^ 2 * L ^ 2 * (τ0 ^ 2 + τ1 ^ 2))) * hτ

end X04

namespace X04
open MeasureTheory Set Filter Real
open scoped Topology

lemma tau_facts (m : ℕ) (hm : 2 ≤ m) (k : ℕ) :
    (∀ i, (vertex m (k + 1) - vertex m k) i = side m k * tau m k i) ∧
    tau m k 0 ^ 2 + tau m k 1 ^ 2 = 1 := by
  have hL := side_pos m hm k
  have hn : side m k = ‖vertex m (k + 1) - vertex m k‖ := by rw [side, dist_eq_norm]
  have hL2 : side m k ^ 2 = (vertex m (k + 1) - vertex m k) 0 ^ 2 +
      (vertex m (k + 1) - vertex m k) 1 ^ 2 := by
    rw [hn, norm_plane, Real.sq_sqrt (by positivity)]
  refine ⟨fun i => ?_, ?_⟩
  · simp only [tau, PiLp.smul_apply, smul_eq_mul]; field_simp
  · simp only [tau, PiLp.smul_apply, smul_eq_mul]
    rw [mul_pow, mul_pow, ← mul_add, ← hL2]; field_simp

lemma apo_sq (m : ℕ) (hm : 2 ≤ m) (k : ℕ) :
    apo m k ^ 2 = (vertex m k 0 ^ 2 + vertex m k 1 ^ 2) -
      (vertex m k 0 * tau m k 0 + vertex m k 1 * tau m k 1) ^ 2 := by
  have hτ := (tau_facts m hm k).2
  have hC : ‖vertex m k‖ ^ 2 = vertex m k 0 ^ 2 + vertex m k 1 ^ 2 := by
    rw [norm_plane, Real.sq_sqrt (by positivity)]
  unfold apo
  rw [Real.sq_sqrt, hC]
  · rfl
  · rw [hC, dot]
    nlinarith [sq_nonneg (vertex m k 0 * tau m k 1 - vertex m k 1 * tau m k 0)]

lemma speed (m : ℕ) (hm : 2 ≤ m) (p : ℝ) (k : ℕ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    ∃ v : Plane, HasDerivAt (arc m p k) v t ∧ ‖v‖ = G m p k t := by
  have hF := Fe_pos m hm p k ht
  obtain ⟨F', hF'⟩ : ∃ F', HasDerivAt (Fe m p k) F' t := ⟨_, Fe_hasDeriv m hm p k ht⟩
  have hd := hF'.deriv
  have hA := (hF'.rpow_const (p := -(1 / p)) (Or.inl hF.ne')).smul (edgePt_hasDeriv m k t)
  refine ⟨_, hA, ?_⟩
  obtain ⟨hΔ, hτ⟩ := tau_facts m hm k
  have hapo := apo_sq m hm k
  set F := Fe m p k t
  set E := F ^ (-(1 / p) - 1) with hE
  have hEpos : 0 < E := Real.rpow_pos_of_pos hF _
  have hE1 : F ^ (-(1 / p)) = E * F := by
    rw [hE, Real.rpow_sub_one hF.ne']; field_simp
  have hE2 : F ^ (-1 - 1 / p) = E := by rw [hE]; congr 1; ring
  unfold G
  rw [hd, hE2, norm_plane]
  simp only [PiLp.add_apply, PiLp.smul_apply, smul_eq_mul, hE1]
  set D := F' / p
  set L := side m k
  have key : ∀ i, (E * F * (vertex m (k + 1) - vertex m k) i + F' * -(1 / p) * E * edgePt m k t i)
      = E * (F * (L * tau m k i) - D * (vertex m k i + t * (L * tau m k i))) := by
    intro i
    rw [← hΔ i]
    simp only [edgePt, PiLp.add_apply, PiLp.smul_apply, smul_eq_mul, D]
    ring
  rw [key 0, key 1, mul_pow, mul_pow, ← mul_add,
    speed_identity _ _ _ _ t F D L hτ, ← hapo]
  have hs : sigma m k t = (vertex m k 0 + t * (L * tau m k 0)) * tau m k 0 +
      (vertex m k 1 + t * (L * tau m k 1)) * tau m k 1 := by
    rw [sigma, dot, ← hΔ 0, ← hΔ 1]
    simp only [edgePt, PiLp.add_apply, PiLp.smul_apply, smul_eq_mul]
  rw [← hs, Real.sqrt_mul (by positivity), Real.sqrt_sq hEpos.le]

end X04

theorem X04.proof_4 : X04.stmt_4 := fun m hm p _ k _ ht => X04.speed m hm p k ht

namespace X04
open MeasureTheory Set Filter Real
open scoped Topology

lemma plane_ext' {v w : Plane} (h0 : v 0 = w 0) (h1 : v 1 = w 1) : v = w := by
  ext i; fin_cases i
  · exact h0
  · exact h1

lemma polar (x : Plane) : ∃ r θ : ℝ, 0 ≤ r ∧ 0 ≤ θ ∧ θ < 2 * π ∧
    x 0 = r * -Real.sin θ ∧ x 1 = r * Real.cos θ := by
  set z : ℂ := ⟨x 0, x 1⟩
  have h0 : ‖z‖ * Real.cos (Complex.arg z) = x 0 := Complex.norm_mul_cos_arg z
  have h1 : ‖z‖ * Real.sin (Complex.arg z) = x 1 := Complex.norm_mul_sin_arg z
  have ha := Complex.neg_pi_lt_arg z
  have hb := Complex.arg_le_pi z
  have hpi := Real.pi_pos
  by_cases hc : Complex.arg z - π / 2 < 0
  · refine ⟨‖z‖, Complex.arg z - π / 2 + 2 * π, norm_nonneg _, by linarith, by linarith, ?_, ?_⟩
    · rw [← h0, Real.sin_add_two_pi, Real.sin_sub_pi_div_two]; ring
    · rw [← h1, Real.cos_add_two_pi, Real.cos_sub_pi_div_two]
  · refine ⟨‖z‖, Complex.arg z - π / 2, norm_nonneg _, by linarith, by linarith, ?_, ?_⟩
    · rw [← h0, Real.sin_sub_pi_div_two]; ring
    · rw [← h1, Real.cos_sub_pi_div_two]

lemma abs_comb2 {a b u w : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) (huw : 0 ≤ u * w) :
    |a * u + b * w| = a * |u| + b * |w| := by
  rcases le_total 0 u with hu | hu <;> rcases le_total 0 w with hw | hw
  · rw [abs_of_nonneg hu, abs_of_nonneg hw, abs_of_nonneg (by positivity)]
  · have : u * w = 0 := le_antisymm (mul_nonpos_of_nonneg_of_nonpos hu hw) huw
    rcases mul_eq_zero.1 this with h | h <;> subst h <;> simp [abs_mul, abs_of_nonneg ha, abs_of_nonneg hb]
  · have : u * w = 0 := le_antisymm (mul_nonpos_of_nonpos_of_nonneg hu hw) huw
    rcases mul_eq_zero.1 this with h | h <;> subst h <;> simp [abs_mul, abs_of_nonneg ha, abs_of_nonneg hb]
  · rw [abs_of_nonpos hu, abs_of_nonpos hw, abs_of_nonpos (by nlinarith)]; ring

lemma Fp_one_comb (m : ℕ) (hm : 2 ≤ m) (k : ℕ) {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) :
    Fp m 1 (a • vertex m k + b • vertex m (k + 1)) = a + b := by
  have h0 := Fp_one_vertex m hm k
  have h1 := Fp_one_vertex m hm (k + 1)
  unfold Fp at *
  simp only [Real.rpow_one] at *
  have : ∀ i ∈ Finset.Icc 1 m, |dot (nrm m i) (a • vertex m k + b • vertex m (k + 1))| =
      a * |dot (nrm m i) (vertex m k)| + b * |dot (nrm m i) (vertex m (k + 1))| := by
    intro i _
    have e : dot (nrm m i) (a • vertex m k + b • vertex m (k + 1)) =
        a * dot (nrm m i) (vertex m k) + b * dot (nrm m i) (vertex m (k + 1)) := by
      simp only [dot, PiLp.add_apply, PiLp.smul_apply, smul_eq_mul]; ring
    rw [e]; exact abs_comb2 ha hb (dot_vertex_prod m hm i k).1
  rw [Finset.sum_congr rfl this, Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum,
    h0, h1, mul_one, mul_one]

lemma vertex_add_two_m (m : ℕ) (hm : 2 ≤ m) (j : ℕ) : vertex m (j + 2 * m) = vertex m j := by
  have hm0 : (m:ℝ) ≠ 0 := by have : (2:ℝ) ≤ m := by exact_mod_cast hm
                             positivity
  have e : ((j + 2 * m : ℕ) : ℝ) * π / m = (j : ℝ) * π / m + 2 * π := by
    push_cast; field_simp
  unfold vertex; rw [e, Real.sin_add_two_pi, Real.cos_add_two_pi]

lemma vertex_add_m (m : ℕ) (hm : 2 ≤ m) : vertex m m = -vertex m 0 := by
  have hm0 : (m:ℝ) ≠ 0 := by have : (2:ℝ) ≤ m := by exact_mod_cast hm
                             positivity
  have e : (m:ℝ) * π / m = π := by field_simp
  apply plane_ext' <;> simp only [PiLp.neg_apply, vertex0, vertex1, e] <;> simp

/-- Every point is a nonnegative combination of two consecutive vertices. -/
lemma sector (m : ℕ) (hm : 2 ≤ m) (x : Plane) : ∃ k < 2 * m, ∃ a b : ℝ, 0 ≤ a ∧ 0 ≤ b ∧
    x = a • vertex m k + b • vertex m (k + 1) := by
  obtain ⟨r, θ, hr, hθ0, hθ1, hx0, hx1⟩ := polar x
  have hm2 : (2:ℝ) ≤ m := by exact_mod_cast hm
  have hmpos : (0:ℝ) < m := by linarith
  have hpi := Real.pi_pos
  set β := π / m with hβ
  have hβpos : 0 < β := by positivity
  have hβle : β ≤ π / 2 := by
    rw [hβ, div_le_div_iff₀ hmpos (by norm_num)]; nlinarith
  obtain ⟨k, hk⟩ : ∃ k, k = ⌊θ / β⌋₊ := ⟨_, rfl⟩
  have hk1 : (k : ℝ) ≤ θ / β := hk ▸ Nat.floor_le (by positivity)
  have hk2 : θ / β < k + 1 := hk ▸ Nat.lt_floor_add_one _
  have hk2m : k < 2 * m := by
    have : θ / β < 2 * m := by
      rw [div_lt_iff₀ hβpos, hβ]; field_simp; linarith
    have : (k:ℝ) < 2 * m := lt_of_le_of_lt hk1 this
    exact_mod_cast this
  obtain ⟨δ, hδ⟩ : ∃ δ, δ = θ - k * β := ⟨_, rfl⟩
  have hδ0 : 0 ≤ δ := by
    have := (le_div_iff₀ hβpos).1 hk1; linarith
  have hδ1 : δ < β := by
    have := (div_lt_iff₀ hβpos).1 hk2; linarith
  have hsβ : 0 < Real.sin β := Real.sin_pos_of_pos_of_lt_pi hβpos (by linarith)
  have hs1 : 0 ≤ Real.sin (β - δ) := Real.sin_nonneg_of_nonneg_of_le_pi (by linarith) (by linarith)
  have hs2 : 0 ≤ Real.sin δ := Real.sin_nonneg_of_nonneg_of_le_pi hδ0 (by linarith)
  have hρ := rho_pos m hm
  refine ⟨k, hk2m, r * Real.sin (β - δ) / (Real.tan (π / (2 * (m : ℝ))) * Real.sin β),
    r * Real.sin δ / (Real.tan (π / (2 * (m : ℝ))) * Real.sin β),
    by positivity, by positivity, ?_⟩
  have hα : (k : ℝ) * π / m = k * β := by rw [hβ]; ring
  have hα1 : ((k + 1 : ℕ) : ℝ) * π / m = k * β + β := by rw [hβ]; push_cast; ring
  have hθ : θ = k * β + δ := by rw [hδ]; ring
  rw [hθ] at hx0 hx1
  apply plane_ext'
  · rw [hx0]
    simp only [PiLp.add_apply, PiLp.smul_apply, smul_eq_mul, vertex0, hα, hα1]
    rw [Real.sin_add (k * β) δ, Real.sin_add (k * β) β, Real.sin_sub β δ]
    field_simp
    ring
  · rw [hx1]
    simp only [PiLp.add_apply, PiLp.smul_apply, smul_eq_mul, vertex1, hα, hα1]
    rw [Real.cos_add (k * β) δ, Real.cos_add (k * β) β, Real.sin_sub β δ]
    field_simp
    ring

lemma Kp_convex (m : ℕ) : Convex ℝ (Kp m 1) := by
  intro x hx y hy a b ha hb hab
  simp only [Kp, mem_ofPred_eq, Fp, Real.rpow_one] at *
  calc ∑ i ∈ Finset.Icc 1 m, |dot (nrm m i) (a • x + b • y)|
      ≤ ∑ i ∈ Finset.Icc 1 m, (a * |dot (nrm m i) x| + b * |dot (nrm m i) y|) := by
        apply Finset.sum_le_sum; intro i _
        have e : dot (nrm m i) (a • x + b • y) = a * dot (nrm m i) x + b * dot (nrm m i) y := by
          simp only [dot, PiLp.add_apply, PiLp.smul_apply, smul_eq_mul]; ring
        rw [e]
        calc |a * dot (nrm m i) x + b * dot (nrm m i) y|
            ≤ |a * dot (nrm m i) x| + |b * dot (nrm m i) y| := abs_add_le _ _
          _ = _ := by rw [abs_mul, abs_mul, abs_of_nonneg ha, abs_of_nonneg hb]
    _ = a * ∑ i ∈ Finset.Icc 1 m, |dot (nrm m i) x| + b * ∑ i ∈ Finset.Icc 1 m, |dot (nrm m i) y| := by
        rw [Finset.sum_add_distrib, Finset.mul_sum, Finset.mul_sum]
    _ ≤ a * 1 + b * 1 := by gcongr
    _ = 1 := by rw [mul_one, mul_one, hab]

lemma Kp_hull (m : ℕ) (hm : 2 ≤ m) :
    Kp m 1 = convexHull ℝ (vertex m '' {k | k < 2 * m}) := by
  set H := convexHull ℝ (vertex m '' {k | k < 2 * m})
  have hH : Convex ℝ H := convex_convexHull ℝ _
  apply le_antisymm
  · intro x hx
    have mem : ∀ k, k < 2 * m → vertex m k ∈ H := fun k hk => subset_convexHull ℝ _ ⟨k, hk, rfl⟩
    have mem1 : ∀ k, k < 2 * m → vertex m (k + 1) ∈ H := by
      intro k hk
      rcases (show k + 1 < 2 * m ∨ k + 1 = 2 * m by omega) with h | h
      · exact mem _ h
      · rw [h, show 2 * m = 0 + 2 * m by ring, vertex_add_two_m m hm 0]; exact mem 0 (by omega)
    have h0 : (0 : Plane) ∈ H := by
      have := hH (mem 0 (by omega)) (mem m (by omega)) (by norm_num : (0:ℝ) ≤ 1/2)
        (by norm_num : (0:ℝ) ≤ 1/2) (by norm_num)
      rw [vertex_add_m m hm] at this
      simpa [smul_neg] using this
    obtain ⟨k, hk, a, b, ha, hb, rfl⟩ := sector m hm x
    have hab : a + b ≤ 1 := by
      have := Fp_one_comb m hm k ha hb
      simp only [Kp, mem_ofPred_eq] at hx; linarith
    rcases eq_or_lt_of_le (add_nonneg ha hb) with h | h
    · have ha0 : a = 0 := by linarith
      have hb0 : b = 0 := by linarith
      subst ha0 hb0; simpa using h0
    · have hy := hH (mem k hk) (mem1 k hk) (div_nonneg ha h.le) (div_nonneg hb h.le)
        (by field_simp)
      have := hH hy h0 h.le (by linarith : (0:ℝ) ≤ 1 - (a + b)) (by ring)
      convert this using 1
      rw [smul_zero, add_zero, smul_add, smul_smul, smul_smul]
      congr 2 <;> field_simp
  · apply convexHull_min _ (Kp_convex m)
    rintro _ ⟨k, _, rfl⟩
    simp only [Kp, mem_ofPred_eq, Fp_one_vertex m hm k, le_refl]

end X04

theorem X04.proof_1 : X04.stmt_1 := fun m hm =>
  ⟨X04.Kp_hull m hm, X04.Fp_one_vertex m hm, X04.side_eq m hm⟩


namespace X04
open MeasureTheory Set Filter Real
open scoped Topology

lemma rho_le_one (m : ℕ) (hm : 2 ≤ m) : Real.tan (π / (2 * (m : ℝ))) ≤ 1 := by
  obtain ⟨h1, h2⟩ := ang_facts m hm
  rw [← Real.tan_pi_div_four]
  rcases eq_or_lt_of_le h2 with h | h
  · rw [h]
  · exact (Real.tan_lt_tan_of_lt_of_lt_pi_div_two (by linarith) (by linarith [Real.pi_pos]) h).le

lemma abs_dot_vertex_le (m : ℕ) (hm : 2 ≤ m) (i k : ℕ) : |dot (nrm m i) (vertex m k)| ≤ 1 := by
  rw [abs_dot_vertex m i k hm]
  have := rho_le_one m hm
  have := rho_pos m hm
  have := Real.abs_sin_le_one (((i : ℝ) - k) * π / (m : ℝ))
  nlinarith [abs_nonneg (Real.sin (((i : ℝ) - k) * π / (m : ℝ)))]

/-- the smallest nonzero `|⟨n_i, C_k⟩|`. -/
noncomputable def mu (m : ℕ) : ℝ := Real.tan (π / (2 * (m : ℝ))) * Real.sin (π / m)

lemma sin_pi_div_pos (m : ℕ) (hm : 2 ≤ m) : 0 < Real.sin (π / m) := by
  have hm2 : (2:ℝ) ≤ m := by exact_mod_cast hm
  apply Real.sin_pos_of_pos_of_lt_pi (by positivity)
  rw [div_lt_iff₀ (by positivity)]; nlinarith [Real.pi_pos]

lemma mu_pos (m : ℕ) (hm : 2 ≤ m) : 0 < mu m := mul_pos (rho_pos m hm) (sin_pi_div_pos m hm)

lemma mu_le_one (m : ℕ) (hm : 2 ≤ m) : mu m ≤ 1 := by
  unfold mu
  have := rho_le_one m hm
  have := rho_pos m hm
  have := Real.sin_le_one (π / m)
  have := sin_pi_div_pos m hm
  nlinarith

lemma abs_sin_int_ge (m : ℕ) (hm : 2 ≤ m) (j : ℤ) (hj : Real.sin ((j : ℝ) * π / m) ≠ 0) :
    Real.sin (π / m) ≤ |Real.sin ((j : ℝ) * π / m)| := by
  have hm2 : (2:ℝ) ≤ m := by exact_mod_cast hm
  have hm0 : (m:ℝ) ≠ 0 := by positivity
  have hmz : (0:ℤ) < m := by exact_mod_cast (by omega : 0 < m)
  set q := j / (m:ℤ)
  set r := j % (m:ℤ)
  have hr0 : 0 ≤ r := Int.emod_nonneg _ hmz.ne'
  have hrm : r < m := Int.emod_lt_of_pos _ hmz
  have hj' : (j:ℝ) = (r:ℝ) + (q:ℝ) * m := by
    have := Int.emod_add_mul_ediv j m
    have h' : ((j:ℤ):ℝ) = ((r + m * q : ℤ):ℝ) := by rw [this]
    push_cast at h'; rw [h']; ring
  have e1 : (j:ℝ) * π / m = (r:ℝ) * π / m + q * π := by rw [hj']; field_simp
  rw [e1, Real.sin_add_int_mul_pi] at hj ⊢
  rw [abs_mul, abs_zpow, abs_neg, abs_one, one_zpow, one_mul]
  have hpi := Real.pi_pos
  rcases (show r = 0 ∨ 1 ≤ r by omega) with h | h
  · exfalso; apply hj; rw [h]; simp
  · have hrR1 : (1:ℝ) ≤ r := by exact_mod_cast h
    have hrRm : (r:ℝ) + 1 ≤ m := by exact_mod_cast hrm
    have hx1 : π / m ≤ (r:ℝ) * π / m := by
      rw [div_le_div_iff_of_pos_right (by positivity)]; nlinarith
    have hx2 : (r:ℝ) * π / m ≤ π - π / m := by
      rw [div_le_iff₀ (by positivity), sub_mul, div_mul_cancel₀ _ hm0]; nlinarith
    have hpm : π / m ≤ π / 2 := by
      rw [div_le_div_iff₀ (by positivity) (by norm_num)]; nlinarith
    have hpm0 : 0 < π / m := by positivity
    have hs0 : 0 ≤ Real.sin ((r:ℝ) * π / m) :=
      Real.sin_nonneg_of_nonneg_of_le_pi (by positivity) (by linarith)
    rw [abs_of_nonneg hs0]
    rcases le_total ((r:ℝ) * π / m) (π / 2) with h3 | h3
    · exact Real.sin_le_sin_of_le_of_le_pi_div_two (by linarith) h3 hx1
    · rw [← Real.sin_pi_sub ((r:ℝ) * π / m)]
      exact Real.sin_le_sin_of_le_of_le_pi_div_two (by linarith) (by linarith) (by linarith)

lemma abs_dot_ge (m : ℕ) (hm : 2 ≤ m) (i k : ℕ) (h : dot (nrm m i) (vertex m k) ≠ 0) :
    mu m ≤ |dot (nrm m i) (vertex m k)| := by
  rw [abs_dot_vertex m i k hm]
  have e : ((i : ℝ) - k) = (((i:ℤ) - k : ℤ) : ℝ) := by push_cast; ring
  rw [dot_nrm_vertex, e] at h
  rw [e]
  have h' : Real.sin ((((i:ℤ) - k : ℤ) : ℝ) * π / m) ≠ 0 := fun h0 => h (by rw [h0, mul_zero])
  exact mul_le_mul_of_nonneg_left (abs_sin_int_ge m hm _ h') (rho_pos m hm).le

lemma A_lower (m : ℕ) (hm : 2 ≤ m) (k i : ℕ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    mu m * (t * (1 - t)) ≤ |aE m k i t| := by
  have hmu := mu_pos m hm
  have t0 := ht.1; have t1 := ht.2
  rw [aE, abs_dot_edge m hm i k t0.le t1.le]
  have hα := abs_nonneg (dot (nrm m i) (vertex m k))
  have hβ := abs_nonneg (dot (nrm m i) (vertex m (k + 1)))
  have hmu1 := mu_le_one m hm
  rcases (dot_vertex_prod m hm i k).2 with h | h
  · have := abs_dot_ge m hm i k h
    have a1 := mul_le_mul_of_nonneg_left this (sub_nonneg.2 t1.le)
    have a2 : mu m * (t * (1 - t)) ≤ (1 - t) * mu m := by
      nlinarith [mul_nonneg (mul_nonneg hmu.le (sub_nonneg.2 t1.le)) (sub_nonneg.2 t1.le)]
    nlinarith [mul_nonneg t0.le hβ]
  · have := abs_dot_ge m hm i (k + 1) h
    have a1 := mul_le_mul_of_nonneg_left this t0.le
    have a2 : mu m * (t * (1 - t)) ≤ t * mu m := by
      nlinarith [mul_nonneg (mul_nonneg hmu.le t0.le) t0.le]
    nlinarith [mul_nonneg (sub_nonneg.2 t1.le) hα]

lemma A_sum (m : ℕ) (hm : 2 ≤ m) (k : ℕ) {t : ℝ} (t0 : 0 ≤ t) (t1 : t ≤ 1) :
    ∑ i ∈ Finset.Icc 1 m, |aE m k i t| = 1 := by
  have := Fe_one m hm k t0 t1
  rw [Fe_eq] at this; simpa using this

lemma A_le_one (m : ℕ) (hm : 2 ≤ m) (k i : ℕ) (hi : i ∈ Finset.Icc 1 m) {t : ℝ}
    (t0 : 0 ≤ t) (t1 : t ≤ 1) : |aE m k i t| ≤ 1 := by
  rw [← A_sum m hm k t0 t1]
  exact Finset.single_le_sum (f := fun i => |aE m k i t|) (fun j _ => abs_nonneg _) hi

lemma cE_abs_le (m : ℕ) (hm : 2 ≤ m) (k i : ℕ) (t : ℝ) : |cE m k i t| ≤ 2 := by
  unfold cE
  rw [abs_mul]
  have h1 : |(SignType.sign (aE m k i t) : ℝ)| ≤ 1 := by
    rcases lt_trichotomy (aE m k i t) 0 with h | h | h
    · simp [sign_neg h]
    · simp [h]
    · simp [sign_pos h]
  have h2 : |dot (nrm m i) (vertex m (k + 1)) - dot (nrm m i) (vertex m k)| ≤ 2 := by
    have := abs_dot_vertex_le m hm i (k + 1)
    have := abs_dot_vertex_le m hm i k
    calc _ ≤ |dot (nrm m i) (vertex m (k + 1))| + |dot (nrm m i) (vertex m k)| := abs_sub _ _
      _ ≤ 2 := by linarith
  nlinarith [abs_nonneg (SignType.sign (aE m k i t) : ℝ),
    abs_nonneg (dot (nrm m i) (vertex m (k + 1)) - dot (nrm m i) (vertex m k))]

lemma vertex_sq (m : ℕ) (k : ℕ) :
    vertex m k 0 ^ 2 + vertex m k 1 ^ 2 = Real.tan (π / (2 * (m : ℝ))) ^ 2 := by
  rw [vertex0, vertex1]
  have := Real.sin_sq_add_cos_sq ((k : ℝ) * π / (m : ℝ))
  nlinarith

lemma side_le_two (m : ℕ) (hm : 2 ≤ m) (k : ℕ) : side m k ≤ 2 := by
  rw [side_eq m hm]
  have := rho_le_one m hm
  have := rho_pos m hm
  have := Real.sin_le_one (π / (2 * (m : ℝ)))
  obtain ⟨h1, h2⟩ := ang_facts m hm
  have : 0 < Real.sin (π / (2 * (m : ℝ))) := Real.sin_pos_of_pos_of_lt_pi h1 (by linarith [Real.pi_pos])
  nlinarith

lemma cross_vertex (m : ℕ) (hm : 2 ≤ m) (k : ℕ) :
    vertex m k 0 * vertex m (k + 1) 1 - vertex m k 1 * vertex m (k + 1) 0 =
      Real.tan (π / (2 * (m : ℝ))) ^ 2 * Real.sin (π / m) := by
  rw [vertex0, vertex0, vertex1, vertex1]
  have hb : ((k + 1 : ℕ) : ℝ) * π / m = (k : ℝ) * π / m + π / m := by push_cast; ring
  rw [hb, Real.sin_add, Real.cos_add]
  have := Real.sin_sq_add_cos_sq ((k : ℝ) * π / (m : ℝ))
  linear_combination (Real.tan (π / (2 * (m : ℝ))) ^ 2 * Real.sin (π / m)) * this

lemma apo_pos (m : ℕ) (hm : 2 ≤ m) (k : ℕ) : 0 < apo m k := by
  have hL := side_pos m hm k
  obtain ⟨hΔ, hτ⟩ := tau_facts m hm k
  have e : apo m k ^ 2 = (vertex m k 0 * tau m k 1 - vertex m k 1 * tau m k 0) ^ 2 := by
    rw [apo_sq m hm k]; linear_combination (-(vertex m k 0 ^ 2 + vertex m k 1 ^ 2)) * hτ
  have hc : vertex m k 0 * tau m k 1 - vertex m k 1 * tau m k 0 ≠ 0 := by
    have h0 := hΔ 0; have h1 := hΔ 1
    simp only [PiLp.sub_apply] at h0 h1
    have hcr := cross_vertex m hm k
    have hpos : 0 < Real.tan (π / (2 * (m : ℝ))) ^ 2 * Real.sin (π / m) := by
      have := rho_pos m hm; have := sin_pi_div_pos m hm; positivity
    intro hz
    have : side m k * (vertex m k 0 * tau m k 1 - vertex m k 1 * tau m k 0) =
        Real.tan (π / (2 * (m : ℝ))) ^ 2 * Real.sin (π / m) := by
      rw [← hcr]; linear_combination (-vertex m k 0) * h1 + vertex m k 1 * h0
    rw [hz, mul_zero] at this; linarith
  have : 0 < apo m k ^ 2 := by rw [e]; positivity
  exact lt_of_le_of_ne (Real.sqrt_nonneg _) (fun h => by rw [← h] at this; simp at this)

lemma apo_le_one (m : ℕ) (hm : 2 ≤ m) (k : ℕ) : apo m k ≤ 1 := by
  have h := apo_sq m hm k
  rw [vertex_sq] at h
  have := rho_le_one m hm
  have := rho_pos m hm
  have ha := (apo_pos m hm k).le
  nlinarith [sq_nonneg (vertex m k 0 * tau m k 0 + vertex m k 1 * tau m k 1)]

lemma sigma_abs_le (m : ℕ) (hm : 2 ≤ m) (k : ℕ) {t : ℝ} (t0 : 0 ≤ t) (t1 : t ≤ 1) :
    |sigma m k t| ≤ 2 := by
  obtain ⟨_, hτ⟩ := tau_facts m hm k
  have hτ0 : |tau m k 0| ≤ 1 := by
    rw [abs_le]; constructor <;> nlinarith [sq_nonneg (tau m k 1)]
  have hτ1 : |tau m k 1| ≤ 1 := by
    rw [abs_le]; constructor <;> nlinarith [sq_nonneg (tau m k 0)]
  have hρ := rho_le_one m hm
  have hv : ∀ j, ∀ c : Fin 2, |vertex m j c| ≤ 1 := by
    intro j c
    have := vertex_sq m j
    have hρ0 := rho_pos m hm
    fin_cases c <;> simp only [Fin.zero_eta, Fin.mk_one] <;> rw [abs_le] <;> constructor <;>
      nlinarith [sq_nonneg (vertex m j 0), sq_nonneg (vertex m j 1)]
  have hy : ∀ c : Fin 2, |edgePt m k t c| ≤ 1 := by
    intro c
    have e : edgePt m k t c = (1 - t) * vertex m k c + t * vertex m (k + 1) c := by
      simp only [edgePt, PiLp.add_apply, PiLp.smul_apply, PiLp.sub_apply, smul_eq_mul]; ring
    rw [e]
    calc |(1 - t) * vertex m k c + t * vertex m (k + 1) c|
        ≤ |(1 - t) * vertex m k c| + |t * vertex m (k + 1) c| := abs_add_le _ _
      _ = (1 - t) * |vertex m k c| + t * |vertex m (k + 1) c| := by
          rw [abs_mul, abs_mul, abs_of_nonneg (by linarith), abs_of_nonneg t0]
      _ ≤ (1 - t) * 1 + t * 1 := by
          gcongr <;> first | linarith | exact hv _ _
      _ = 1 := by ring
  unfold sigma dot
  calc |edgePt m k t 0 * tau m k 0 + edgePt m k t 1 * tau m k 1|
      ≤ |edgePt m k t 0 * tau m k 0| + |edgePt m k t 1 * tau m k 1| := abs_add_le _ _
    _ = |edgePt m k t 0| * |tau m k 0| + |edgePt m k t 1| * |tau m k 1| := by rw [abs_mul, abs_mul]
    _ ≤ 1 * 1 + 1 * 1 := by
        gcongr <;> first | exact hy 0 | exact hy 1
    _ = 2 := by norm_num

end X04

namespace X04
open MeasureTheory Set Filter Real
open scoped Topology

lemma cs2 (a b c d : ℝ) : |a * b + c * d| ≤ Real.sqrt (a ^ 2 + c ^ 2) * Real.sqrt (b ^ 2 + d ^ 2) := by
  rw [← Real.sqrt_mul (by positivity)]
  apply Real.abs_le_sqrt
  nlinarith [sq_nonneg (a * d - b * c)]

lemma sqrt_two_le (a b : ℝ) : Real.sqrt (a ^ 2 + b ^ 2) ≤ |a| + |b| := by
  rw [Real.sqrt_le_left (by positivity)]
  nlinarith [sq_abs a, sq_abs b, abs_nonneg a, abs_nonneg b, mul_nonneg (abs_nonneg a) (abs_nonneg b)]

/-- abstract form of `G` as a function of `q` at a fixed interior point. -/
lemma Gq_deriv (I : Finset ℕ) (A c : ℕ → ℝ) (hA : ∀ i ∈ I, 0 < A i) (L σ h p : ℝ)
    (hL : 0 < L) (hh : 0 < h) (hp : p ≠ 0) (hF : 0 < ∑ i ∈ I, A i ^ p) :
    ∃ g, HasDerivAt (fun q => (∑ i ∈ I, A i ^ q) ^ (-1 - 1 / q) *
        Real.sqrt ((L * ∑ i ∈ I, A i ^ q - σ * ∑ i ∈ I, c i * A i ^ (q - 1)) ^ 2 +
          h ^ 2 * (∑ i ∈ I, c i * A i ^ (q - 1)) ^ 2)) g p ∧
      |g| ≤ (|∑ i ∈ I, A i ^ p * Real.log (A i)| * |-1 - 1 / p| *
              (∑ i ∈ I, A i ^ p) ^ (-1 - 1 / p - 1) +
            |1 / p ^ 2| * (∑ i ∈ I, A i ^ p) ^ (-1 - 1 / p) * |Real.log (∑ i ∈ I, A i ^ p)|) *
          (L * ∑ i ∈ I, A i ^ p + |σ| * |∑ i ∈ I, c i * A i ^ (p - 1)| +
            h * |∑ i ∈ I, c i * A i ^ (p - 1)|) +
        (∑ i ∈ I, A i ^ p) ^ (-1 - 1 / p) *
          (L * |∑ i ∈ I, A i ^ p * Real.log (A i)| +
            |σ| * |∑ i ∈ I, c i * A i ^ (p - 1) * Real.log (A i)| +
            h * |∑ i ∈ I, c i * A i ^ (p - 1) * Real.log (A i)|) := by
  set F := ∑ i ∈ I, A i ^ p
  set F1 := ∑ i ∈ I, A i ^ p * Real.log (A i)
  set D := ∑ i ∈ I, c i * A i ^ (p - 1)
  set D1 := ∑ i ∈ I, c i * A i ^ (p - 1) * Real.log (A i)
  have dF : HasDerivAt (fun q : ℝ => ∑ i ∈ I, A i ^ q) F1 p :=
    HasDerivAt.fun_sum (fun i hi => (hasStrictDerivAt_const_rpow (hA i hi) p).hasDerivAt)
  have dD : HasDerivAt (fun q : ℝ => ∑ i ∈ I, c i * A i ^ (q - 1)) D1 p := by
    apply HasDerivAt.fun_sum
    intro i hi
    have := (((hasDerivAt_id' p).sub_const 1).const_rpow (hA i hi)).const_mul (c i)
    convert this using 1; ring
  have dE : HasDerivAt (fun q : ℝ => -1 - 1 / q) (1 / p ^ 2) p := by
    have := (hasDerivAt_inv hp).const_sub (-1)
    convert this using 1
    · funext q; simp [one_div]
    · field_simp
  have dFE := dF.rpow dE hF
  set u := L * F - σ * D
  set X := u ^ 2 + h ^ 2 * D ^ 2 with hX
  have hXpos : 0 < X := by
    rcases eq_or_ne D 0 with h0 | h0
    · have : u = L * F := by simp [u, h0]
      have hu : u ^ 2 = X := by rw [hX, h0]; ring
      rw [← hu, this]; positivity
    · have : 0 < h ^ 2 * D ^ 2 := by positivity
      have := sq_nonneg u
      linarith
  have dX : HasDerivAt (fun q : ℝ =>
      (L * ∑ i ∈ I, A i ^ q - σ * ∑ i ∈ I, c i * A i ^ (q - 1)) ^ 2 +
        h ^ 2 * (∑ i ∈ I, c i * A i ^ (q - 1)) ^ 2)
      (2 * u * (L * F1 - σ * D1) + h ^ 2 * (2 * D * D1)) p := by
    convert (((dF.const_mul L).sub (dD.const_mul σ)).fun_pow 2).add
      ((dD.fun_pow 2).const_mul (h ^ 2)) using 1
    simp only [Pi.sub_apply, u, F, D]; push_cast; ring
  have dS := dX.sqrt hXpos.ne'
  refine ⟨_, dFE.mul dS, ?_⟩
  have hsX : Real.sqrt X ≤ L * F + |σ| * |D| + h * |D| := by
    have e : X = u ^ 2 + (h * D) ^ 2 := by rw [hX]; ring
    rw [e]
    refine (sqrt_two_le _ _).trans ?_
    rw [abs_mul, abs_of_pos hh]
    have : |u| ≤ L * F + |σ| * |D| := by
      calc |u| ≤ |L * F| + |σ * D| := abs_sub _ _
        _ = L * F + |σ| * |D| := by rw [abs_mul, abs_mul, abs_of_pos hL, abs_of_pos hF]
    linarith
  have hsq : 0 < Real.sqrt X := Real.sqrt_pos.2 hXpos
  have hX' : |2 * u * (L * F1 - σ * D1) + h ^ 2 * (2 * D * D1)| / (2 * Real.sqrt X) ≤
      L * |F1| + |σ| * |D1| + h * |D1| := by
    rw [div_le_iff₀ (by positivity)]
    have e : 2 * u * (L * F1 - σ * D1) + h ^ 2 * (2 * D * D1) =
        2 * (u * (L * F1 - σ * D1) + (h * D) * (h * D1)) := by ring
    have eX : X = u ^ 2 + (h * D) ^ 2 := by rw [hX]; ring
    rw [e, abs_mul, abs_two]
    have c := cs2 u (L * F1 - σ * D1) (h * D) (h * D1)
    rw [← eX] at c
    have hs2 := sqrt_two_le (L * F1 - σ * D1) (h * D1)
    have hs3 : |L * F1 - σ * D1| + |h * D1| ≤ L * |F1| + |σ| * |D1| + h * |D1| := by
      rw [abs_mul h, abs_of_pos hh]
      have : |L * F1 - σ * D1| ≤ L * |F1| + |σ| * |D1| := by
        calc |L * F1 - σ * D1| ≤ |L * F1| + |σ * D1| := abs_sub _ _
          _ = _ := by rw [abs_mul, abs_mul, abs_of_pos hL]
      linarith
    nlinarith
  have hFe := Real.rpow_nonneg hF.le (-1 - 1 / p)
  have hFe1 := Real.rpow_nonneg hF.le (-1 - 1 / p - 1)
  have h1 : |F1 * (-1 - 1 / p) * (∑ i ∈ I, A i ^ p) ^ (-1 - 1 / p - 1) +
      1 / p ^ 2 * (∑ i ∈ I, A i ^ p) ^ (-1 - 1 / p) * Real.log (∑ i ∈ I, A i ^ p)| ≤
      |F1| * |-1 - 1 / p| * F ^ (-1 - 1 / p - 1) + |1 / p ^ 2| * F ^ (-1 - 1 / p) * |Real.log F| := by
    refine (abs_add_le _ _).trans ?_
    rw [abs_mul, abs_mul, abs_mul, abs_mul, abs_of_nonneg hFe, abs_of_nonneg hFe1]
  calc _ ≤ |F1 * (-1 - 1 / p) * (∑ i ∈ I, A i ^ p) ^ (-1 - 1 / p - 1) +
      1 / p ^ 2 * (∑ i ∈ I, A i ^ p) ^ (-1 - 1 / p) * Real.log (∑ i ∈ I, A i ^ p)| * Real.sqrt X +
      F ^ (-1 - 1 / p) * (|2 * u * (L * F1 - σ * D1) + h ^ 2 * (2 * D * D1)| / (2 * Real.sqrt X)) := by
        refine (abs_add_le _ _).trans (le_of_eq ?_)
        rw [abs_mul, abs_mul, abs_div, abs_of_nonneg hFe, abs_of_pos (by positivity : 0 < 2 * Real.sqrt X),
          abs_of_nonneg (Real.sqrt_nonneg _)]
    _ ≤ _ := by
        gcongr

end X04

namespace X04
open MeasureTheory Set Filter Real
open scoped Topology

lemma combine {m L σ h F F1 D D1 e e' lF E E1 P Q : ℝ} (hm : 1 ≤ m) (hP : 1 ≤ P) (hQ : 1 ≤ Q)
    (hL0 : 0 ≤ L) (hL : L ≤ 2) (hσ : |σ| ≤ 2) (hh0 : 0 ≤ h) (hh : h ≤ 1) (hF0 : 0 ≤ F)
    (hF : F ≤ m) (hE0 : 0 ≤ E) (hE : E ≤ m ^ 3) (hE10 : 0 ≤ E1) (hE1 : E1 ≤ m ^ 4)
    (he : |e| ≤ 3) (he' : |e'| ≤ 2) (hlF : |lF| ≤ m) (hF1 : |F1| ≤ m * Q) (hD : |D| ≤ 2 * m * P)
    (hD1 : |D1| ≤ 2 * m * P * Q) :
    (|F1| * |e| * E1 + |e'| * E * |lF|) * (L * F + |σ| * |D| + h * |D|) +
      E * (L * |F1| + |σ| * |D1| + h * |D1|) ≤ 100 * m ^ 6 * P * Q := by
  have hm0 : 0 ≤ m := by linarith
  calc _ ≤ (m * Q * 3 * m ^ 4 + 2 * m ^ 3 * m) * (2 * m + 2 * (2 * m * P) + 1 * (2 * m * P)) +
        m ^ 3 * (2 * (m * Q) + 2 * (2 * m * P * Q) + 1 * (2 * m * P * Q)) := by
        gcongr
    _ ≤ 100 * m ^ 6 * P * Q := by
        have h1 : m ^ 4 ≤ m ^ 6 := pow_le_pow_right₀ hm (by norm_num)
        have h2 : m ^ 5 ≤ m ^ 6 := pow_le_pow_right₀ hm (by norm_num)
        have hPQ : 1 ≤ P * Q := by nlinarith
        have h6 : 0 ≤ m ^ 6 := by positivity
        have a1 : m ^ 6 * Q ≤ m ^ 6 * (P * Q) :=
          mul_le_mul_of_nonneg_left (le_mul_of_one_le_left (by linarith) hP) h6
        have a2 : m ^ 5 ≤ m ^ 6 * (P * Q) := h2.trans (le_mul_of_one_le_right h6 hPQ)
        have a3 : m ^ 5 * P ≤ m ^ 6 * (P * Q) :=
          (mul_le_mul_of_nonneg_right h2 (by linarith)).trans
            (mul_le_mul_of_nonneg_left (le_mul_of_one_le_right (by linarith) hQ) h6)
        have a4 : m ^ 4 * Q ≤ m ^ 6 * (P * Q) :=
          (mul_le_mul_of_nonneg_right h1 (by linarith)).trans a1
        have a5 : m ^ 4 * (P * Q) ≤ m ^ 6 * (P * Q) :=
          mul_le_mul_of_nonneg_right h1 (by linarith)
        nlinarith

lemma rpow_le_pow_of_mem {x m : ℝ} (n : ℕ) {e : ℝ} (hm : 1 ≤ m) (hx0 : 1 / m ≤ x) (hx1 : x ≤ m)
    (he0 : -(n : ℝ) ≤ e) (he1 : e ≤ 0) : x ^ e ≤ m ^ n := by
  have hmpos : 0 < m := by linarith
  have hx : 0 < x := lt_of_lt_of_le (by positivity) hx0
  rcases le_total 1 x with h | h
  · exact (Real.rpow_le_one_of_one_le_of_nonpos h he1).trans (one_le_pow₀ hm)
  · calc x ^ e ≤ x ^ (-(n : ℝ)) := Real.rpow_le_rpow_of_exponent_ge hx h he0
      _ = (x ^ n)⁻¹ := by rw [Real.rpow_neg hx.le, Real.rpow_natCast]
      _ ≤ ((1 / m) ^ n)⁻¹ := by
          gcongr
      _ = m ^ n := by rw [one_div, inv_pow, inv_inv]

lemma G_eq_form (m : ℕ) (hm : 2 ≤ m) (k : ℕ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) {q : ℝ} (hq : q ≠ 0) :
    G m q k t = (∑ i ∈ Finset.Icc 1 m, |aE m k i t| ^ q) ^ (-1 - 1 / q) *
        Real.sqrt ((side m k * ∑ i ∈ Finset.Icc 1 m, |aE m k i t| ^ q -
          sigma m k t * ∑ i ∈ Finset.Icc 1 m, cE m k i t * |aE m k i t| ^ (q - 1)) ^ 2 +
          apo m k ^ 2 * (∑ i ∈ Finset.Icc 1 m, cE m k i t * |aE m k i t| ^ (q - 1)) ^ 2) := by
  unfold G
  rw [(Fe_hasDeriv m hm q k ht).deriv]
  have e : (∑ i ∈ Finset.Icc 1 m, cE m k i t * q * |aE m k i t| ^ (q - 1)) / q =
      ∑ i ∈ Finset.Icc 1 m, cE m k i t * |aE m k i t| ^ (q - 1) := by
    rw [Finset.sum_div]; refine Finset.sum_congr rfl (fun i _ => ?_); field_simp
  rw [e]; rfl

lemma card_Icc_one (m : ℕ) : ((Finset.Icc 1 m).card : ℝ) = m := by
  rw [Nat.card_Icc]; simp

set_option maxHeartbeats 2000000 in
lemma stmt7_core (m : ℕ) (hm : 2 ≤ m) (k : ℕ) {p : ℝ} (hp : p ∈ Icc (3 / 4 : ℝ) (5 / 4))
    {t : ℝ} (ht : t ∈ Ioo (0 : ℝ) 1) :
    ∃ g : ℝ, HasDerivAt (fun q => G m q k t) g p ∧
      |g| ≤ (100 * (m : ℝ) ^ 6 * mu m ^ (-(1 / 4 : ℝ)) * (1 + |Real.log (mu m)|)) *
        t ^ (-(1 / 4 : ℝ)) * (1 - t) ^ (-(1 / 4 : ℝ)) * (1 + |Real.log t| + |Real.log (1 - t)|) := by
  set I := Finset.Icc 1 m
  set A : ℕ → ℝ := fun i => |aE m k i t| with hAdef
  set c : ℕ → ℝ := fun i => cE m k i t with hcdef
  have hp0 : 0 < p := by linarith [hp.1]
  have t0 := ht.1
  have t1 : 0 < 1 - t := by linarith [ht.2]
  have hA : ∀ i ∈ I, 0 < A i := fun i _ => abs_pos.2 (aE_ne m hm k i ht)
  have hA1 : ∀ i ∈ I, A i ≤ 1 := fun i hi => A_le_one m hm k i hi ht.1.le ht.2.le
  have hAsum : ∑ i ∈ I, A i = 1 := A_sum m hm k ht.1.le ht.2.le
  have hmu := mu_pos m hm
  have hmu1 := mu_le_one m hm
  have hAlow : ∀ i ∈ I, mu m * (t * (1 - t)) ≤ A i := fun i _ => A_lower m hm k i ht
  have hm1 : (1:ℝ) ≤ m := by exact_mod_cast (by omega : 1 ≤ m)
  have hm0 : (0:ℝ) < m := by linarith
  -- F bounds
  have hFup : ∑ i ∈ I, A i ^ p ≤ m := by
    calc ∑ i ∈ I, A i ^ p ≤ ∑ i ∈ I, (1:ℝ) :=
          Finset.sum_le_sum (fun i hi => Real.rpow_le_one (hA i hi).le (hA1 i hi) hp0.le)
      _ = m := by rw [Finset.sum_const, nsmul_eq_mul, card_Icc_one, mul_one]
  have hFlow : 1 / (m : ℝ) ≤ ∑ i ∈ I, A i ^ p := by
    have h2 : ∑ i ∈ I, A i ^ 2 ≤ ∑ i ∈ I, A i ^ p := by
      apply Finset.sum_le_sum; intro i hi
      rw [← Real.rpow_two]
      exact Real.rpow_le_rpow_of_exponent_ge (hA i hi) (hA1 i hi) (by linarith [hp.2])
    have h3 := sq_sum_le_card_mul_sum_sq (s := I) (f := A)
    rw [hAsum, card_Icc_one] at h3
    rw [div_le_iff₀ hm0]; nlinarith
  have hFpos : 0 < ∑ i ∈ I, A i ^ p := lt_of_lt_of_le (by positivity) hFlow
  obtain ⟨g, hg, hgb⟩ := Gq_deriv I A c hA (side m k) (sigma m k t) (apo m k) p (side_pos m hm k)
    (apo_pos m hm k) hp0.ne' hFpos
  refine ⟨g, ?_, ?_⟩
  · refine hg.congr_of_eventuallyEq ?_
    filter_upwards [lt_mem_nhds hp0] with q hq
    exact G_eq_form m hm k ht hq.ne'
  refine hgb.trans ?_
  -- the parameters
  set W := t ^ (-(1 / 4 : ℝ)) * (1 - t) ^ (-(1 / 4 : ℝ))
  set Λ := 1 + |Real.log t| + |Real.log (1 - t)|
  set M1 := mu m ^ (-(1 / 4 : ℝ))
  set M2 := 1 + |Real.log (mu m)|
  have hW : 1 ≤ W := by
    have a := Real.one_le_rpow_of_pos_of_le_one_of_nonpos t0 ht.2.le (by norm_num : -(1/4:ℝ) ≤ 0)
    have b := Real.one_le_rpow_of_pos_of_le_one_of_nonpos t1 (by linarith) (by norm_num : -(1/4:ℝ) ≤ 0)
    nlinarith
  have hM1 : 1 ≤ M1 := Real.one_le_rpow_of_pos_of_le_one_of_nonpos hmu hmu1 (by norm_num)
  have hΛ : 1 ≤ Λ := by have := abs_nonneg (Real.log t); have := abs_nonneg (Real.log (1 - t)); linarith
  have hM2 : 1 ≤ M2 := by have := abs_nonneg (Real.log (mu m)); linarith
  have hlogA : ∀ i ∈ I, |Real.log (A i)| ≤ M2 * Λ := by
    intro i hi
    have hl0 : Real.log (A i) ≤ 0 := Real.log_nonpos (hA i hi).le (hA1 i hi)
    have hl1 : Real.log (mu m * (t * (1 - t))) ≤ Real.log (A i) :=
      Real.log_le_log (by positivity) (hAlow i hi)
    rw [Real.log_mul hmu.ne' (by positivity), Real.log_mul t0.ne' t1.ne'] at hl1
    rw [abs_of_nonpos hl0]
    have := neg_abs_le (Real.log (mu m)); have := neg_abs_le (Real.log t)
    have := neg_abs_le (Real.log (1 - t))
    have : |Real.log (mu m)| * (|Real.log t| + |Real.log (1 - t)|) ≥ 0 := by positivity
    nlinarith
  have hpowA : ∀ i ∈ I, A i ^ (p - 1) ≤ M1 * W := by
    intro i hi
    calc A i ^ (p - 1) ≤ A i ^ (-(1 / 4 : ℝ)) :=
          Real.rpow_le_rpow_of_exponent_ge (hA i hi) (hA1 i hi) (by linarith [hp.1])
      _ ≤ (mu m * (t * (1 - t))) ^ (-(1 / 4 : ℝ)) :=
          Real.rpow_le_rpow_of_nonpos (by positivity) (hAlow i hi) (by norm_num)
      _ = M1 * W := by
          rw [Real.mul_rpow hmu.le (by positivity), Real.mul_rpow t0.le t1.le]
  have hcard : ∀ x : ℝ, ∑ i ∈ I, x = m * x := fun x => by
    rw [Finset.sum_const, nsmul_eq_mul, card_Icc_one]
  have hF1 : |∑ i ∈ I, A i ^ p * Real.log (A i)| ≤ m * (M2 * Λ) := by
    refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
    rw [← hcard]
    apply Finset.sum_le_sum; intro i hi
    rw [abs_mul, abs_of_nonneg (Real.rpow_nonneg (hA i hi).le _)]
    have := Real.rpow_le_one (hA i hi).le (hA1 i hi) hp0.le
    have := hlogA i hi
    nlinarith [abs_nonneg (Real.log (A i)), Real.rpow_nonneg (hA i hi).le p]
  have hD : |∑ i ∈ I, c i * A i ^ (p - 1)| ≤ 2 * m * (M1 * W) := by
    refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
    rw [show 2 * (m:ℝ) * (M1 * W) = m * (2 * (M1 * W)) by ring, ← hcard]
    apply Finset.sum_le_sum; intro i hi
    rw [abs_mul, abs_of_nonneg (Real.rpow_nonneg (hA i hi).le _)]
    have := cE_abs_le m hm k i t
    have := hpowA i hi
    nlinarith [abs_nonneg (c i), Real.rpow_nonneg (hA i hi).le (p - 1)]
  have hD1 : |∑ i ∈ I, c i * A i ^ (p - 1) * Real.log (A i)| ≤ 2 * m * (M1 * W) * (M2 * Λ) := by
    refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
    rw [show 2 * (m:ℝ) * (M1 * W) * (M2 * Λ) = m * (2 * (M1 * W) * (M2 * Λ)) by ring, ← hcard]
    apply Finset.sum_le_sum; intro i hi
    rw [abs_mul, abs_mul, abs_of_nonneg (Real.rpow_nonneg (hA i hi).le _)]
    have h1 := cE_abs_le m hm k i t
    have h2 := hpowA i hi
    have h3 := hlogA i hi
    have h4 := Real.rpow_nonneg (hA i hi).le (p - 1)
    calc |c i| * A i ^ (p - 1) * |Real.log (A i)| ≤ 2 * (M1 * W) * (M2 * Λ) := by gcongr
      _ = _ := rfl
  have hE : (∑ i ∈ I, A i ^ p) ^ (-1 - 1 / p) ≤ (m:ℝ) ^ 3 := by
    apply rpow_le_pow_of_mem 3 hm1 hFlow hFup
    · have : 1 / p ≤ 4 / 3 := by rw [div_le_iff₀ hp0]; linarith [hp.1]
      push_cast; linarith
    · have : 0 < 1 / p := by positivity
      linarith
  have hE1 : (∑ i ∈ I, A i ^ p) ^ (-1 - 1 / p - 1) ≤ (m:ℝ) ^ 4 := by
    apply rpow_le_pow_of_mem 4 hm1 hFlow hFup
    · have : 1 / p ≤ 4 / 3 := by rw [div_le_iff₀ hp0]; linarith [hp.1]
      push_cast; linarith
    · have : 0 < 1 / p := by positivity
      linarith
  have he : |-1 - 1 / p| ≤ 3 := by
    have : 1 / p ≤ 4 / 3 := by rw [div_le_iff₀ hp0]; linarith [hp.1]
    have : 0 < 1 / p := by positivity
    rw [abs_le]; constructor <;> linarith
  have he' : |1 / p ^ 2| ≤ 2 := by
    rw [abs_of_pos (by positivity), div_le_iff₀ (by positivity)]; nlinarith [hp.1]
  have hlF : |Real.log (∑ i ∈ I, A i ^ p)| ≤ m := by
    rw [abs_le]; constructor
    · have := Real.log_le_log (by positivity) hFlow
      rw [one_div, Real.log_inv] at this
      have := Real.log_le_sub_one_of_pos hm0
      linarith
    · have := Real.log_le_log hFpos hFup
      have := Real.log_le_sub_one_of_pos hm0
      linarith
  have key := combine (m := (m:ℝ)) (P := M1 * W) (Q := M2 * Λ) hm1 (by nlinarith) (by nlinarith)
    (side_pos m hm k).le (side_le_two m hm k) (sigma_abs_le m hm k ht.1.le ht.2.le)
    (apo_pos m hm k).le (apo_le_one m hm k) hFpos.le hFup (Real.rpow_nonneg hFpos.le _) hE
    (Real.rpow_nonneg hFpos.le _) hE1 he he' hlF hF1 hD hD1
  refine key.trans (le_of_eq ?_)
  ring

end X04

theorem X04.proof_7 : X04.stmt_7 := by
  intro m hm
  exact ⟨_, fun k p hp t ht => X04.stmt7_core m hm k hp ht⟩

namespace X04
open MeasureTheory Set Filter Real
open scoped Topology

lemma aE_cont (m k i : ℕ) : Continuous (aE m k i) := by
  have e : aE m k i = fun s => dot (nrm m i) (vertex m k) +
      s * (dot (nrm m i) (vertex m (k + 1)) - dot (nrm m i) (vertex m k)) := by
    funext s; rw [aE, dot_edge]; ring
  rw [e]; fun_prop

lemma Fe_meas (m : ℕ) (q : ℝ) (k : ℕ) : Measurable (Fe m q k) := by
  rw [Fe_eq]
  exact Finset.measurable_sum _ (fun i _ =>
    ((aE_cont m k i).measurable.abs).pow_const q)

lemma sigma_cont (m k : ℕ) : Continuous (sigma m k) := by
  have e : sigma m k = fun t => (vertex m k + t • (vertex m (k + 1) - vertex m k)) 0 * tau m k 0 +
      (vertex m k + t • (vertex m (k + 1) - vertex m k)) 1 * tau m k 1 := rfl
  rw [e]
  simp only [PiLp.add_apply, PiLp.smul_apply, smul_eq_mul]
  fun_prop

lemma G_meas (m : ℕ) (q : ℝ) (k : ℕ) : Measurable (G m q k) := by
  have h1 := Fe_meas m q k
  have h2 : Measurable (deriv (Fe m q k)) := measurable_deriv _
  have h3 := (sigma_cont m k).measurable
  unfold G
  exact (h1.pow_const _).mul (((((h1.const_mul _).sub (h3.mul (h2.div_const _))).pow_const 2).add
    ((h2.div_const _).pow_const 2 |>.const_mul _)).sqrt)

lemma S1e_eq (m k : ℕ) : (fun s => S1 m (edgePt m k s)) =
    fun s => ∑ i ∈ Finset.Icc 1 m, |aE m k i s| * Real.log |aE m k i s| := rfl

lemma S1e_cont (m k : ℕ) : Continuous (fun s => S1 m (edgePt m k s)) := by
  rw [S1e_eq]
  exact continuous_finsetSum _ (fun i _ => Real.continuous_mul_log.comp (aE_cont m k i).abs)

lemma ae_Ioo {P : ℝ → Prop} (h : ∀ t ∈ Ioo (0:ℝ) 1, P t) :
    ∀ᵐ t ∂(volume.restrict (uIoc (0:ℝ) 1)), P t := by
  rw [uIoc_of_le zero_le_one, ← Measure.restrict_congr_set Ioo_ae_eq_Ioc]
  filter_upwards [ae_restrict_mem measurableSet_Ioo] with t ht using h t ht

lemma ae_Ioo' {P : ℝ → Prop} (h : ∀ t ∈ Ioo (0:ℝ) 1, P t) :
    ∀ᵐ t ∂volume, t ∈ uIoc (0:ℝ) 1 → P t := by
  filter_upwards [compl_mem_ae_iff.2 (measure_singleton (1:ℝ))] with t ht h01
  rw [uIoc_of_le zero_le_one] at h01
  exact h t ⟨h01.1, lt_of_le_of_ne h01.2 ht⟩

lemma half_bound {a b : ℝ} (ha : 0 < a) (ha1 : a ≤ 1) (hb : 1 / 2 ≤ b) (hb1 : b ≤ 1) :
    a ^ (-(1 / 4 : ℝ)) * b ^ (-(1 / 4 : ℝ)) * (1 + |Real.log a| + |Real.log b|) ≤
      8 * a ^ (-(3 / 4 : ℝ)) := by
  have hb0 : 0 < b := by linarith
  have hb2 : b ^ (-(1 / 4 : ℝ)) ≤ 2 := by
    calc b ^ (-(1 / 4 : ℝ)) ≤ b ^ (-(1:ℝ)) := Real.rpow_le_rpow_of_exponent_ge hb0 hb1 (by norm_num)
      _ = 1 / b := by rw [Real.rpow_neg_one, one_div]
      _ ≤ 2 := by rw [div_le_iff₀ hb0]; linarith
  have hlb : |Real.log b| ≤ 1 := by
    rw [abs_of_nonpos (Real.log_nonpos hb0.le hb1)]
    have := Real.log_le_sub_one_of_pos (inv_pos.2 hb0)
    rw [Real.log_inv] at this
    have : b⁻¹ ≤ 2 := by rw [inv_le_comm₀ hb0 (by norm_num)]; linarith
    linarith
  have hla : |Real.log a| ≤ 2 * a ^ (-(1 / 2 : ℝ)) := by
    rw [abs_of_nonpos (Real.log_nonpos ha.le ha1)]
    have := Real.log_le_sub_one_of_pos (Real.rpow_pos_of_pos ha (-(1 / 2 : ℝ)))
    rw [Real.log_rpow ha] at this
    linarith
  have hp1 : 0 ≤ a ^ (-(1 / 4 : ℝ)) := Real.rpow_nonneg ha.le _
  have e1 : a ^ (-(1 / 4 : ℝ)) * a ^ (-(1 / 2 : ℝ)) = a ^ (-(3 / 4 : ℝ)) := by
    rw [← Real.rpow_add ha]; norm_num
  have e2 : a ^ (-(1 / 4 : ℝ)) ≤ a ^ (-(3 / 4 : ℝ)) :=
    Real.rpow_le_rpow_of_exponent_ge ha ha1 (by norm_num)
  have hq : 0 ≤ 1 + |Real.log a| + |Real.log b| := by positivity
  calc a ^ (-(1 / 4 : ℝ)) * b ^ (-(1 / 4 : ℝ)) * (1 + |Real.log a| + |Real.log b|)
      ≤ a ^ (-(1 / 4 : ℝ)) * 2 * (1 + 2 * a ^ (-(1 / 2 : ℝ)) + 1) := by gcongr
    _ = 4 * a ^ (-(1 / 4 : ℝ)) + 4 * (a ^ (-(1 / 4 : ℝ)) * a ^ (-(1 / 2 : ℝ))) := by ring
    _ ≤ 8 * a ^ (-(3 / 4 : ℝ)) := by rw [e1]; linarith

lemma bound_int (K : ℝ) : IntervalIntegrable (fun t : ℝ => K * t ^ (-(1 / 4 : ℝ)) *
    (1 - t) ^ (-(1 / 4 : ℝ)) * (1 + |Real.log t| + |Real.log (1 - t)|)) volume 0 1 := by
  have i1 : IntervalIntegrable (fun x : ℝ => x ^ (-(3 / 4 : ℝ))) volume 0 1 :=
    intervalIntegral.intervalIntegrable_rpow' (by norm_num)
  have i2 : IntervalIntegrable (fun x : ℝ => (1 - x) ^ (-(3 / 4 : ℝ))) volume 0 1 := by
    have := i1.comp_sub_left 1
    simpa using this.symm
  have i3 := ((i1.add i2).const_mul (8 * |K|))
  refine i3.mono_fun ?_ ?_
  · refine Measurable.aestronglyMeasurable ?_
    fun_prop
  · refine ae_Ioo (fun t ht => ?_)
    have t0 := ht.1; have t1 : 0 < 1 - t := by linarith [ht.2]
    have hW : 0 ≤ t ^ (-(1 / 4 : ℝ)) * (1 - t) ^ (-(1 / 4 : ℝ)) * (1 + |Real.log t| + |Real.log (1 - t)|) := by
      positivity
    have hsum : t ^ (-(1 / 4 : ℝ)) * (1 - t) ^ (-(1 / 4 : ℝ)) * (1 + |Real.log t| + |Real.log (1 - t)|)
        ≤ 8 * (t ^ (-(3 / 4 : ℝ)) + (1 - t) ^ (-(3 / 4 : ℝ))) := by
      have p1 := Real.rpow_nonneg t0.le (-(3 / 4 : ℝ))
      have p2 := Real.rpow_nonneg t1.le (-(3 / 4 : ℝ))
      rcases le_total t (1 / 2) with h | h
      · have := half_bound t0 ht.2.le (by linarith) (by linarith : 1 - t ≤ 1)
        linarith
      · have := half_bound t1 (by linarith) h ht.2.le
        have e : (1 - t) ^ (-(1 / 4 : ℝ)) * t ^ (-(1 / 4 : ℝ)) * (1 + |Real.log (1 - t)| + |Real.log t|)
            = t ^ (-(1 / 4 : ℝ)) * (1 - t) ^ (-(1 / 4 : ℝ)) * (1 + |Real.log t| + |Real.log (1 - t)|) := by
          ring
        linarith
    simp only [Real.norm_eq_abs]
    rw [abs_of_nonneg (by positivity : (0:ℝ) ≤ 8 * |K| * (t ^ (-(3 / 4 : ℝ)) + (1 - t) ^ (-(3 / 4 : ℝ))))]
    have e : K * t ^ (-(1 / 4 : ℝ)) * (1 - t) ^ (-(1 / 4 : ℝ)) * (1 + |Real.log t| + |Real.log (1 - t)|) =
        K * (t ^ (-(1 / 4 : ℝ)) * (1 - t) ^ (-(1 / 4 : ℝ)) * (1 + |Real.log t| + |Real.log (1 - t)|)) := by ring
    rw [e, abs_mul, abs_of_nonneg hW]
    calc |K| * (t ^ (-(1 / 4 : ℝ)) * (1 - t) ^ (-(1 / 4 : ℝ)) * (1 + |Real.log t| + |Real.log (1 - t)|))
        ≤ |K| * (8 * (t ^ (-(3 / 4 : ℝ)) + (1 - t) ^ (-(3 / 4 : ℝ)))) := by gcongr
      _ = _ := by ring

end X04

namespace X04
open MeasureTheory Set Filter Real
open scoped Topology

lemma sigma_hasDeriv (m : ℕ) (hm : 2 ≤ m) (k : ℕ) (t : ℝ) :
    HasDerivAt (sigma m k) (side m k) t := by
  obtain ⟨hΔ, hτ⟩ := tau_facts m hm k
  have e : sigma m k = fun s => dot (vertex m k) (tau m k) +
      s * ((vertex m (k + 1) - vertex m k) 0 * tau m k 0 +
        (vertex m (k + 1) - vertex m k) 1 * tau m k 1) := by
    funext s
    simp only [sigma, dot, edgePt, PiLp.add_apply, PiLp.smul_apply, PiLp.sub_apply, smul_eq_mul]
    ring
  have hv : (vertex m (k + 1) - vertex m k) 0 * tau m k 0 +
      (vertex m (k + 1) - vertex m k) 1 * tau m k 1 = side m k := by
    rw [hΔ 0, hΔ 1]; linear_combination side m k * hτ
  rw [e, hv]
  simpa using ((hasDerivAt_id' t).mul_const (side m k)).const_add (dot (vertex m k) (tau m k))

lemma edgePt_zero (m k : ℕ) : edgePt m k 0 = vertex m k := by simp [edgePt]
lemma edgePt_one (m k : ℕ) : edgePt m k 1 = vertex m (k + 1) := by simp [edgePt]

lemma dot_comm' (u v : Plane) : dot u v = dot v u := by unfold dot; ring

/-- the explicit constant of stmt_7. -/
noncomputable def K7 (m : ℕ) : ℝ :=
  100 * (m : ℝ) ^ 6 * mu m ^ (-(1 / 4 : ℝ)) * (1 + |Real.log (mu m)|)

lemma stmt9_core (m : ℕ) (hm : 2 ≤ m) (k : ℕ) :
    HasDerivAt (fun p => ∫ t in (0 : ℝ)..1, G m p k t)
      (dot (tau m k) (vertex m k) * S1 m (vertex m k) -
        dot (tau m k) (vertex m (k + 1)) * S1 m (vertex m (k + 1))) 1 := by
  set E : ℝ → ℝ := fun t => -(side m k * S1 m (edgePt m k t) +
      sigma m k t * deriv (fun s => S1 m (edgePt m k s)) t) with hEdef
  have hs : Ioo (3 / 4 : ℝ) (5 / 4) ∈ 𝓝 (1 : ℝ) := Ioo_mem_nhds (by norm_num) (by norm_num)
  have hEq : ∀ t ∈ Ioo (0:ℝ) 1, E t = deriv (fun q => G m q k t) 1 := fun t ht =>
    (stmt8_core m hm k ht).deriv.symm
  have hEmeas : Measurable E := by
    have h1 := (S1e_cont m k).measurable
    have h2 : Measurable (deriv (fun s => S1 m (edgePt m k s))) := measurable_deriv _
    have h3 := (sigma_cont m k).measurable
    exact ((h1.const_mul _).add (h3.mul h2)).neg
  have hint1 : IntervalIntegrable (G m 1 k) volume 0 1 := by
    refine (intervalIntegrable_const (c := side m k)).congr_uIoo ?_
    intro t ht
    rw [uIoo_of_le zero_le_one] at ht
    exact (proof_6 m hm k t ht).symm
  have hmeas' : AEStronglyMeasurable (fun t => deriv (fun q => G m q k t) 1)
      (volume.restrict (uIoc (0:ℝ) 1)) :=
    hEmeas.aestronglyMeasurable.congr (ae_Ioo hEq)
  have hbound : ∀ᵐ t ∂volume, t ∈ uIoc (0:ℝ) 1 → ∀ x ∈ Ioo (3 / 4 : ℝ) (5 / 4),
      ‖deriv (fun q => G m q k t) x‖ ≤ K7 m * t ^ (-(1 / 4 : ℝ)) * (1 - t) ^ (-(1 / 4 : ℝ)) *
        (1 + |Real.log t| + |Real.log (1 - t)|) := by
    refine ae_Ioo' (fun t ht x hx => ?_)
    obtain ⟨g, hg, hgb⟩ := stmt7_core m hm k (Ioo_subset_Icc_self hx) ht
    rw [hg.deriv, Real.norm_eq_abs]; exact hgb
  have hdiff : ∀ᵐ t ∂volume, t ∈ uIoc (0:ℝ) 1 → ∀ x ∈ Ioo (3 / 4 : ℝ) (5 / 4),
      HasDerivAt (fun q => G m q k t) (deriv (fun q => G m q k t) x) x := by
    refine ae_Ioo' (fun t ht x hx => ?_)
    obtain ⟨g, hg, -⟩ := stmt7_core m hm k (Ioo_subset_Icc_self hx) ht
    rw [hg.deriv]; exact hg
  obtain ⟨hI, hD⟩ := intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (F := fun p t => G m p k t) (F' := fun p t => deriv (fun q => G m q k t) p) (x₀ := 1) hs
    (Eventually.of_forall (fun x => (G_meas m x k).aestronglyMeasurable)) hint1 hmeas'
    hbound (bound_int (K7 m)) hdiff
  refine hD.congr_deriv ?_
  have hIE : IntervalIntegrable E volume 0 1 := by
    refine hI.congr_uIoo ?_
    intro t ht; rw [uIoo_of_le zero_le_one] at ht; exact (hEq t ht).symm
  have hint : ∫ t in (0:ℝ)..1, deriv (fun q => G m q k t) 1 = ∫ t in (0:ℝ)..1, E t := by
    refine intervalIntegral.integral_congr_ae (ae_Ioo' (fun t ht => (hEq t ht).symm))
  rw [hint]
  have hFTC := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le (f := fun t => sigma m k t * S1 m (edgePt m k t))
    (f' := fun t => -E t) zero_le_one
    (((sigma_cont m k).mul (S1e_cont m k)).continuousOn)
    (fun t ht => by
      have := (sigma_hasDeriv m hm k t).mul (S1_edge_hasDeriv m hm k ht)
      convert this using 1
      simp only [hEdef, neg_neg]
      rw [(S1_edge_hasDeriv m hm k ht).deriv])
    hIE.neg
  rw [intervalIntegral.integral_neg] at hFTC
  have : ∫ t in (0:ℝ)..1, E t = sigma m k 0 * S1 m (edgePt m k 0) - sigma m k 1 * S1 m (edgePt m k 1) := by
    linarith
  rw [this, sigma, sigma, edgePt_zero, edgePt_one, dot_comm' (vertex m k), dot_comm' (vertex m (k + 1))]

end X04

theorem X04.proof_9 : X04.stmt_9 := X04.stmt9_core

namespace X04
open MeasureTheory Set Filter Real
open scoped Topology

/-! ### length of a curve as the integral of its speed (after W04, abstracted) -/

section VarInt
variable (f : ℝ → Plane) (g : ℝ → ℝ)

/-- hypotheses: `f` continuous on `[0,1]`, differentiable on `(0,1)` with speed `g`. -/
structure SpeedData : Prop where
  cont : ContinuousOn f (Icc 0 1)
  der : ∀ t ∈ Ioo (0:ℝ) 1, HasDerivAt f (deriv f t) t ∧ ‖deriv f t‖ = g t
  int : IntervalIntegrable g volume 0 1
  nonneg : ∀ t, 0 ≤ g t
  gcont : ∀ t ∈ Ioo (0:ℝ) 1, ContinuousAt g t
  meas : Measurable g

variable {f g}

lemma SpeedData.sub (h : SpeedData f g) {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) (hb : b ≤ 1) :
    IntervalIntegrable g volume a b := by
  refine h.int.mono_set ?_
  rw [uIcc_of_le zero_le_one, uIcc_of_le hab]
  exact Icc_subset_Icc ha hb

lemma SpeedData.Gcont (h : SpeedData f g) :
    ContinuousOn (fun x => ∫ t in (0:ℝ)..x, g t) (Icc 0 1) := by
  have := intervalIntegral.continuousOn_primitive_interval' h.int left_mem_uIcc
  rw [uIcc_of_le zero_le_one] at this
  exact this

lemma SpeedData.Gderiv (h : SpeedData f g) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    HasDerivAt (fun x => ∫ t in (0:ℝ)..x, g t) (g t) t :=
  intervalIntegral.integral_hasDerivAt_right (h.sub le_rfl ht.1.le ht.2.le)
    h.meas.stronglyMeasurable.stronglyMeasurableAtFilter (h.gcont t ht)

lemma SpeedData.Gsub (h : SpeedData f g) {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) (hb : b ≤ 1) :
    (∫ t in (0:ℝ)..b, g t) - (∫ t in (0:ℝ)..a, g t) = ∫ t in a..b, g t := by
  rw [intervalIntegral.integral_interval_sub_left (h.sub le_rfl (ha.trans hab) hb)
    (h.sub le_rfl ha (hab.trans hb))]

lemma SpeedData.Gmono (h : SpeedData f g) {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) (hb : b ≤ 1) :
    (∫ t in (0:ℝ)..a, g t) ≤ ∫ t in (0:ℝ)..b, g t := by
  have := h.Gsub ha hab hb
  have h2 : 0 ≤ ∫ t in a..b, g t := intervalIntegral.integral_nonneg hab fun u _ => h.nonneg u
  linarith

lemma SpeedData.up_open (h : SpeedData f g) {x y : ℝ} (hx : 0 < x) (hxy : x ≤ y) (hy : y < 1) :
    ‖f y - f x‖ ≤ (∫ t in (0:ℝ)..y, g t) - ∫ t in (0:ℝ)..x, g t := by
  have hI : ∀ z ∈ Ico x y, z ∈ Ioo (0:ℝ) 1 := fun z hz => ⟨by linarith [hz.1], by linarith [hz.2]⟩
  have hsub : Icc x y ⊆ Icc 0 1 := Icc_subset_Icc hx.le hy.le
  refine image_norm_le_of_norm_deriv_right_le_deriv_boundary' (a := x) (b := y)
    (f := fun z => f z - f x) (f' := deriv f)
    ((h.cont.mono hsub).sub continuousOn_const)
    (fun z hz => ((h.der z (hI z hz)).1.sub_const _).hasDerivWithinAt)
    (B := fun z => (∫ t in (0:ℝ)..z, g t) - ∫ t in (0:ℝ)..x, g t) (B' := g) (by simp)
    ((h.Gcont.mono hsub).sub continuousOn_const)
    (fun z hz => ((h.Gderiv (hI z hz)).sub_const _).hasDerivWithinAt)
    (fun z hz => (h.der z (hI z hz)).2.le) ⟨hxy, le_rfl⟩

lemma eps_tend' {x : ℝ} (hx : x ∈ Icc (0:ℝ) 1) :
    Tendsto (fun e : ℝ => e + (1 - 2 * e) * x) (𝓝[>] 0) (𝓝[Icc 0 1] x) := by
  apply tendsto_nhdsWithin_iff.2 ⟨?_, ?_⟩
  · have : Tendsto (fun e : ℝ => e + (1 - 2 * e) * x) (𝓝 0) (𝓝 (0 + (1 - 2 * 0) * x)) :=
      ((continuous_id.add ((continuous_const.sub (continuous_const.mul continuous_id)).mul
        continuous_const)).tendsto 0)
    simpa using this.mono_left nhdsWithin_le_nhds
  · filter_upwards [Ioo_mem_nhdsGT (show (0:ℝ) < 1/2 by norm_num)] with e he
    constructor <;> nlinarith [he.1, he.2, hx.1, hx.2]

lemma SpeedData.up_closed (h : SpeedData f g) {x y : ℝ} (hx : 0 ≤ x) (hxy : x ≤ y) (hy : y ≤ 1) :
    ‖f y - f x‖ ≤ (∫ t in (0:ℝ)..y, g t) - ∫ t in (0:ℝ)..x, g t := by
  have tx := eps_tend' ⟨hx, hxy.trans hy⟩
  have ty := eps_tend' ⟨hx.trans hxy, hy⟩
  have L : Tendsto (fun e : ℝ => ‖f (e + (1 - 2 * e) * y) - f (e + (1 - 2 * e) * x)‖)
      (𝓝[>] 0) (𝓝 ‖f y - f x‖) :=
    (((h.cont y ⟨hx.trans hxy, hy⟩).tendsto.comp ty).sub
      ((h.cont x ⟨hx, hxy.trans hy⟩).tendsto.comp tx)).norm
  have R : Tendsto (fun e : ℝ => (∫ t in (0:ℝ)..(e + (1 - 2 * e) * y), g t) -
      ∫ t in (0:ℝ)..(e + (1 - 2 * e) * x), g t)
      (𝓝[>] 0) (𝓝 ((∫ t in (0:ℝ)..y, g t) - ∫ t in (0:ℝ)..x, g t)) :=
    ((h.Gcont y ⟨hx.trans hxy, hy⟩).tendsto.comp ty).sub
      ((h.Gcont x ⟨hx, hxy.trans hy⟩).tendsto.comp tx)
  refine le_of_tendsto_of_tendsto L R ?_
  filter_upwards [Ioo_mem_nhdsGT (show (0:ℝ) < 1/2 by norm_num)] with e he
  exact h.up_open (by nlinarith [he.1, he.2]) (by nlinarith [he.1, he.2]) (by nlinarith [he.1, he.2])

lemma var_le' {E F : Type*} [PseudoEMetricSpace E] [PseudoEMetricSpace F]
    (f : ℝ → E) (h : ℝ → F) (s : Set ℝ)
    (hd : ∀ a ∈ s, ∀ b ∈ s, edist (f a) (f b) ≤ edist (h a) (h b)) :
    eVariationOn f s ≤ eVariationOn h s := by
  unfold eVariationOn
  exact iSup_mono fun p => Finset.sum_le_sum fun i _ => hd _ (p.2.2.2 _) _ (p.2.2.2 _)

lemma SpeedData.var_up (h : SpeedData f g) :
    eVariationOn f (Icc 0 1) ≤ ENNReal.ofReal (∫ t in (0:ℝ)..1, g t) := by
  have key : ∀ a ∈ Icc (0:ℝ) 1, ∀ b ∈ Icc (0:ℝ) 1, a ≤ b →
      edist (f a) (f b) ≤ edist (∫ t in (0:ℝ)..a, g t) (∫ t in (0:ℝ)..b, g t) := by
    intro a ha b hb hab
    rw [edist_dist, edist_dist, dist_eq_norm, norm_sub_rev, Real.dist_eq, abs_sub_comm,
      abs_of_nonneg (by linarith [h.Gmono ha.1 hab hb.2])]
    exact ENNReal.ofReal_le_ofReal (h.up_closed ha.1 hab hb.2)
  calc eVariationOn f (Icc 0 1) ≤ eVariationOn (fun x => ∫ t in (0:ℝ)..x, g t) (Icc 0 1) := by
        refine var_le' _ _ _ fun a ha b hb => ?_
        rcases le_total a b with hh | hh
        · exact key a ha b hb hh
        · rw [edist_comm, edist_comm (∫ t in (0:ℝ)..a, g t)]; exact key b hb a ha hh
    _ = ENNReal.ofReal (∫ t in (0:ℝ)..1, g t) := by
        have hm : MonotoneOn (fun x => ∫ t in (0:ℝ)..x, g t) (Icc 0 1) :=
          fun a ha b hb hab => h.Gmono ha.1 hab hb.2
        have := hm.eVariationOn_eq (a := 0) (b := 1) ⟨le_rfl, zero_le_one⟩ ⟨zero_le_one, le_rfl⟩
        simpa using this

lemma V_inc' (hbv : BoundedVariationOn f (Icc 0 1)) {x z : ℝ}
    (hx : 0 ≤ x) (hxz : x ≤ z) (hz : z ≤ 1) :
    ‖f z - f x‖ ≤
      variationOnFromTo f (Icc 0 1) 0 z - variationOnFromTo f (Icc 0 1) 0 x := by
  rw [variationOnFromTo.sub_right hbv.locallyBoundedVariationOn ⟨le_rfl, zero_le_one⟩
    ⟨hx.trans hxz, hz⟩ ⟨hx, hxz.trans hz⟩, variationOnFromTo.eq_of_le _ _ hxz, ← dist_eq_norm,
    dist_comm, dist_edist]
  refine ENNReal.toReal_mono (hbv.mono inter_subset_left) ?_
  exact eVariationOn.edist_le _ ⟨⟨hx, hxz.trans hz⟩, le_rfl, hxz⟩ ⟨⟨hx.trans hxz, hz⟩, hxz, le_rfl⟩

lemma SpeedData.V_cont (h : SpeedData f g) (hbv : BoundedVariationOn f (Icc 0 1)) :
    ContinuousOn (variationOnFromTo f (Icc 0 1) 0) (Icc 0 1) := fun x hx =>
  (hbv.continuousWithinAt_variationOnFromTo_iff ⟨le_rfl, zero_le_one⟩ hx).2 (h.cont x hx)

lemma SpeedData.low_open (h : SpeedData f g) (hbv : BoundedVariationOn f (Icc 0 1))
    {a b : ℝ} (ha : 0 < a) (hab : a ≤ b) (hb : b < 1) :
    (∫ t in (0:ℝ)..b, g t) - (∫ t in (0:ℝ)..a, g t) ≤
      variationOnFromTo f (Icc 0 1) 0 b - variationOnFromTo f (Icc 0 1) 0 a := by
  set V := variationOnFromTo f (Icc 0 1) 0 with hV
  set Gp : ℝ → ℝ := fun x => ∫ t in (0:ℝ)..x, g t with hGp
  set φ : ℝ → ℝ := fun z => (Gp z - Gp a) - (V z - V a) with hφ
  have hsub : Icc a b ⊆ Icc 0 1 := Icc_subset_Icc ha.le hb.le
  have hfc : ContinuousOn φ (Icc a b) :=
    (((h.Gcont).mono hsub).sub continuousOn_const).sub (((h.V_cont hbv).mono hsub).sub continuousOn_const)
  have key := image_le_of_liminf_slope_right_le_deriv_boundary (a := a) (b := b) hfc
    (B := fun _ => 0) (B' := fun _ => 0) (by simp [φ]) continuousOn_const
    (fun _ _ => hasDerivWithinAt_const _ _ _) ?_ ⟨hab, le_rfl⟩
  · simp only [φ] at key; linarith
  intro x hx r hr
  have hx01 : x ∈ Ioo (0:ℝ) 1 := ⟨by linarith [hx.1], by linarith [hx.2]⟩
  obtain ⟨hd, hn⟩ := h.der x hx01
  have t1 : Tendsto (slope f x) (𝓝[>] x) (𝓝 (deriv f x)) :=
    (hasDerivAt_iff_tendsto_slope.1 hd).mono_left (nhdsWithin_mono _ fun z hz => ne_of_gt hz)
  have t2 : Tendsto (slope Gp x) (𝓝[>] x) (𝓝 (g x)) :=
    (hasDerivAt_iff_tendsto_slope.1 (h.Gderiv hx01)).mono_left
      (nhdsWithin_mono _ fun z hz => ne_of_gt hz)
  have t3 : Tendsto (fun z => slope Gp x z - ‖slope f x z‖) (𝓝[>] x)
      (𝓝 (g x - ‖deriv f x‖)) := t2.sub t1.norm
  rw [hn, sub_self] at t3
  refine Filter.Eventually.frequently ?_
  filter_upwards [t3 (Iio_mem_nhds hr), self_mem_nhdsWithin,
    nhdsWithin_le_nhds (Iio_mem_nhds hx01.2)] with z hz hzx hz1
  simp only [mem_preimage, mem_Iio] at hz
  have hzx' : x < z := hzx
  have hpos : 0 < z - x := by linarith
  have hVi := V_inc' hbv hx01.1.le hzx'.le (le_of_lt hz1)
  refine lt_of_le_of_lt ?_ hz
  rw [slope_def_field, slope_def_field, slope_def_module, norm_smul, Real.norm_eq_abs,
    abs_of_pos (inv_pos.2 hpos)]
  simp only [φ]
  rw [div_le_iff₀ hpos]
  have : (z - x)⁻¹ * ‖f z - f x‖ * (z - x) = ‖f z - f x‖ := by
    field_simp
  rw [sub_mul, div_mul_cancel₀ _ hpos.ne', this]
  linarith

lemma SpeedData.var_eq (h : SpeedData f g) :
    eVariationOn f (Icc 0 1) = ENNReal.ofReal (∫ t in (0:ℝ)..1, g t) := by
  have hup := h.var_up
  have hne : eVariationOn f (Icc 0 1) ≠ ⊤ := ne_top_of_le_ne_top ENNReal.ofReal_ne_top hup
  have hbv : BoundedVariationOn f (Icc 0 1) := hne
  set W := (eVariationOn f (Icc 0 1)).toReal with hW
  have hab : ∀ a b : ℝ, 0 < a → a ≤ b → b < 1 →
      (∫ t in (0:ℝ)..b, g t) - (∫ t in (0:ℝ)..a, g t) ≤ W := by
    intro a b ha hab hb
    have h1 := h.low_open hbv ha hab hb
    rw [variationOnFromTo.sub_right hbv.locallyBoundedVariationOn ⟨le_rfl, zero_le_one⟩
      ⟨by linarith, hb.le⟩ ⟨ha.le, by linarith⟩, variationOnFromTo.eq_of_le _ _ hab] at h1
    exact h1.trans (ENNReal.toReal_mono hne (eVariationOn.mono _ inter_subset_left))
  have t0 := eps_tend' (x := 0) ⟨le_rfl, zero_le_one⟩
  have t1 := eps_tend' (x := 1) ⟨zero_le_one, le_rfl⟩
  have R : Tendsto (fun e : ℝ => (∫ t in (0:ℝ)..(e + (1 - 2 * e) * 1), g t) -
      ∫ t in (0:ℝ)..(e + (1 - 2 * e) * 0), g t)
      (𝓝[>] 0) (𝓝 ((∫ t in (0:ℝ)..1, g t) - ∫ t in (0:ℝ)..0, g t)) :=
    ((h.Gcont 1 ⟨zero_le_one, le_rfl⟩).tendsto.comp t1).sub
      ((h.Gcont 0 ⟨le_rfl, zero_le_one⟩).tendsto.comp t0)
  have hG : (∫ t in (0:ℝ)..1, g t) ≤ W := by
    have := le_of_tendsto (b := W) R ?_
    · simpa using this
    filter_upwards [Ioo_mem_nhdsGT (show (0:ℝ) < 1/2 by norm_num)] with e he
    exact hab _ _ (by nlinarith [he.1, he.2]) (by nlinarith [he.1, he.2]) (by nlinarith [he.1, he.2])
  refine le_antisymm hup ?_
  rw [← ENNReal.ofReal_toReal hne]
  exact ENNReal.ofReal_le_ofReal hG

end VarInt

end X04

namespace X04
open MeasureTheory Set Filter Real
open scoped Topology

lemma Fe_pos_Icc (m : ℕ) (hm : 2 ≤ m) (p : ℝ) (k : ℕ) {t : ℝ} (ht : t ∈ Icc (0:ℝ) 1) :
    0 < Fe m p k t := by
  have hs := A_sum m hm k ht.1 ht.2
  obtain ⟨i, hi, hpos⟩ : ∃ i ∈ Finset.Icc 1 m, 0 < |aE m k i t| := by
    by_contra hc
    simp only [not_exists, not_and, not_lt] at hc
    have : ∑ i ∈ Finset.Icc 1 m, |aE m k i t| = 0 :=
      Finset.sum_eq_zero (fun i hi => le_antisymm (hc i hi) (abs_nonneg _))
    linarith
  rw [Fe_eq]
  exact Finset.sum_pos' (fun j _ => Real.rpow_nonneg (abs_nonneg _) _)
    ⟨i, hi, Real.rpow_pos_of_pos hpos _⟩

lemma Fe_cont (m : ℕ) {p : ℝ} (hp : 0 < p) (k : ℕ) : Continuous (Fe m p k) := by
  rw [Fe_eq]
  exact continuous_finsetSum _ (fun i _ =>
    (aE_cont m k i).abs.rpow_const (fun _ => Or.inr hp.le))

lemma edgePt_cont (m k : ℕ) : Continuous (edgePt m k) := by
  have e : edgePt m k = fun s => vertex m k + s • (vertex m (k + 1) - vertex m k) := rfl
  rw [e]; fun_prop

lemma arc_contOn (m : ℕ) (hm : 2 ≤ m) {p : ℝ} (hp : 0 < p) (k : ℕ) :
    ContinuousOn (arc m p k) (Icc 0 1) := by
  have h1 : ContinuousOn (fun t => Fe m p k t ^ (-(1 / p))) (Icc 0 1) :=
    (Fe_cont m hp k).continuousOn.rpow_const (fun t ht => Or.inl (Fe_pos_Icc m hm p k ht).ne')
  exact h1.smul (edgePt_cont m k).continuousOn

lemma cE_eventually (m : ℕ) (hm : 2 ≤ m) (k i : ℕ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    cE m k i =ᶠ[𝓝 t] fun _ => cE m k i t := by
  have hne := aE_ne m hm k i ht
  have hc := (aE_cont m k i).continuousAt (x := t)
  rcases lt_or_gt_of_ne hne with h | h
  · filter_upwards [hc.eventually (gt_mem_nhds h)] with s hs
    simp only [cE, sign_neg hs, sign_neg h]
  · filter_upwards [hc.eventually (lt_mem_nhds h)] with s hs
    simp only [cE, sign_pos hs, sign_pos h]

lemma G_contAt (m : ℕ) (hm : 2 ≤ m) {p : ℝ} (hp : p ≠ 0) (k : ℕ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    ContinuousAt (G m p k) t := by
  have hev : G m p k =ᶠ[𝓝 t] fun s => (∑ i ∈ Finset.Icc 1 m, |aE m k i s| ^ p) ^ (-1 - 1 / p) *
        Real.sqrt ((side m k * ∑ i ∈ Finset.Icc 1 m, |aE m k i s| ^ p -
          sigma m k s * ∑ i ∈ Finset.Icc 1 m, cE m k i s * |aE m k i s| ^ (p - 1)) ^ 2 +
          apo m k ^ 2 * (∑ i ∈ Finset.Icc 1 m, cE m k i s * |aE m k i s| ^ (p - 1)) ^ 2) := by
    filter_upwards [Ioo_mem_nhds ht.1 ht.2] with s hs using G_eq_form m hm k hs hp
  refine ContinuousAt.congr ?_ hev.symm
  have hA : ∀ i ∈ Finset.Icc 1 m, ContinuousAt (fun s => |aE m k i s|) t :=
    fun i _ => (aE_cont m k i).abs.continuousAt
  have hA0 : ∀ i, |aE m k i t| ≠ 0 := fun i => (abs_pos.2 (aE_ne m hm k i ht)).ne'
  have hF : ContinuousAt (fun s => ∑ i ∈ Finset.Icc 1 m, |aE m k i s| ^ p) t :=
    tendsto_finsetSum _ (fun i hi => (hA i hi).rpow_const (Or.inl (hA0 i)))
  have hD : ContinuousAt (fun s => ∑ i ∈ Finset.Icc 1 m, cE m k i s * |aE m k i s| ^ (p - 1)) t :=
    tendsto_finsetSum _ (fun i hi =>
      (continuousAt_const.congr (cE_eventually m hm k i ht).symm).mul
        ((hA i hi).rpow_const (Or.inl (hA0 i))))
  have hFt : (∑ i ∈ Finset.Icc 1 m, |aE m k i t| ^ p) ≠ 0 := by
    have := Fe_pos m hm p k ht; rw [Fe_eq] at this; exact this.ne'
  exact (hF.rpow_const (Or.inl hFt)).mul
    ((((continuousAt_const.mul hF).sub ((sigma_cont m k).continuousAt.mul hD)).pow 2).add
      (continuousAt_const.mul (hD.pow 2))).sqrt

lemma G_nonneg (m : ℕ) (p : ℝ) (k : ℕ) (t : ℝ) : 0 ≤ G m p k t := by
  unfold G
  have : 0 ≤ Fe m p k t := by
    rw [Fe_eq]; exact Finset.sum_nonneg (fun i _ => Real.rpow_nonneg (abs_nonneg _) _)
  exact mul_nonneg (Real.rpow_nonneg this _) (Real.sqrt_nonneg _)

/-- size bound for `G` on `|p−1| ≤ 1/4`. -/
lemma G_bound (m : ℕ) (hm : 2 ≤ m) (k : ℕ) {p : ℝ} (hp : p ∈ Icc (3 / 4 : ℝ) (5 / 4))
    {t : ℝ} (ht : t ∈ Ioo (0 : ℝ) 1) :
    G m p k t ≤ (8 * (m : ℝ) ^ 4 * mu m ^ (-(1 / 4 : ℝ))) *
        t ^ (-(1 / 4 : ℝ)) * (1 - t) ^ (-(1 / 4 : ℝ)) * (1 + |Real.log t| + |Real.log (1 - t)|) := by
  set I := Finset.Icc 1 m
  set A : ℕ → ℝ := fun i => |aE m k i t| with hAdef
  set c : ℕ → ℝ := fun i => cE m k i t with hcdef
  have hp0 : 0 < p := by linarith [hp.1]
  have t0 := ht.1
  have t1 : 0 < 1 - t := by linarith [ht.2]
  have hA : ∀ i ∈ I, 0 < A i := fun i _ => abs_pos.2 (aE_ne m hm k i ht)
  have hA1 : ∀ i ∈ I, A i ≤ 1 := fun i hi => A_le_one m hm k i hi ht.1.le ht.2.le
  have hAsum : ∑ i ∈ I, A i = 1 := A_sum m hm k ht.1.le ht.2.le
  have hmu := mu_pos m hm
  have hmu1 := mu_le_one m hm
  have hAlow : ∀ i ∈ I, mu m * (t * (1 - t)) ≤ A i := fun i _ => A_lower m hm k i ht
  have hm1 : (1:ℝ) ≤ m := by exact_mod_cast (by omega : 1 ≤ m)
  have hm0 : (0:ℝ) < m := by linarith
  have hFup : ∑ i ∈ I, A i ^ p ≤ m := by
    calc ∑ i ∈ I, A i ^ p ≤ ∑ i ∈ I, (1:ℝ) :=
          Finset.sum_le_sum (fun i hi => Real.rpow_le_one (hA i hi).le (hA1 i hi) hp0.le)
      _ = m := by rw [Finset.sum_const, nsmul_eq_mul, card_Icc_one, mul_one]
  have hFlow : 1 / (m : ℝ) ≤ ∑ i ∈ I, A i ^ p := by
    have h2 : ∑ i ∈ I, A i ^ 2 ≤ ∑ i ∈ I, A i ^ p := by
      apply Finset.sum_le_sum; intro i hi
      rw [← Real.rpow_two]
      exact Real.rpow_le_rpow_of_exponent_ge (hA i hi) (hA1 i hi) (by linarith [hp.2])
    have h3 := sq_sum_le_card_mul_sum_sq (s := I) (f := A)
    rw [hAsum, card_Icc_one] at h3
    rw [div_le_iff₀ hm0]; nlinarith
  have hFpos : 0 < ∑ i ∈ I, A i ^ p := lt_of_lt_of_le (by positivity) hFlow
  set W := t ^ (-(1 / 4 : ℝ)) * (1 - t) ^ (-(1 / 4 : ℝ))
  set Λ := 1 + |Real.log t| + |Real.log (1 - t)|
  set M1 := mu m ^ (-(1 / 4 : ℝ))
  have hW : 1 ≤ W := by
    have a := Real.one_le_rpow_of_pos_of_le_one_of_nonpos t0 ht.2.le (by norm_num : -(1/4:ℝ) ≤ 0)
    have b := Real.one_le_rpow_of_pos_of_le_one_of_nonpos t1 (by linarith) (by norm_num : -(1/4:ℝ) ≤ 0)
    nlinarith
  have hM1 : 1 ≤ M1 := Real.one_le_rpow_of_pos_of_le_one_of_nonpos hmu hmu1 (by norm_num)
  have hΛ : 1 ≤ Λ := by have := abs_nonneg (Real.log t); have := abs_nonneg (Real.log (1 - t)); linarith
  have hpowA : ∀ i ∈ I, A i ^ (p - 1) ≤ M1 * W := by
    intro i hi
    calc A i ^ (p - 1) ≤ A i ^ (-(1 / 4 : ℝ)) :=
          Real.rpow_le_rpow_of_exponent_ge (hA i hi) (hA1 i hi) (by linarith [hp.1])
      _ ≤ (mu m * (t * (1 - t))) ^ (-(1 / 4 : ℝ)) :=
          Real.rpow_le_rpow_of_nonpos (by positivity) (hAlow i hi) (by norm_num)
      _ = M1 * W := by
          rw [Real.mul_rpow hmu.le (by positivity), Real.mul_rpow t0.le t1.le]
  have hcard : ∀ x : ℝ, ∑ i ∈ I, x = m * x := fun x => by
    rw [Finset.sum_const, nsmul_eq_mul, card_Icc_one]
  have hD : |∑ i ∈ I, c i * A i ^ (p - 1)| ≤ 2 * m * (M1 * W) := by
    refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
    rw [show 2 * (m:ℝ) * (M1 * W) = m * (2 * (M1 * W)) by ring, ← hcard]
    apply Finset.sum_le_sum; intro i hi
    rw [abs_mul, abs_of_nonneg (Real.rpow_nonneg (hA i hi).le _)]
    have := cE_abs_le m hm k i t
    have := hpowA i hi
    nlinarith [abs_nonneg (c i), Real.rpow_nonneg (hA i hi).le (p - 1)]
  have hE : (∑ i ∈ I, A i ^ p) ^ (-1 - 1 / p) ≤ (m:ℝ) ^ 3 := by
    apply rpow_le_pow_of_mem 3 hm1 hFlow hFup
    · have : 1 / p ≤ 4 / 3 := by rw [div_le_iff₀ hp0]; linarith [hp.1]
      push_cast; linarith
    · have : 0 < 1 / p := by positivity
      linarith
  rw [G_eq_form m hm k ht hp0.ne']
  set F := ∑ i ∈ I, A i ^ p
  set D := ∑ i ∈ I, c i * A i ^ (p - 1)
  set σ := sigma m k t
  have hσ := sigma_abs_le m hm k ht.1.le ht.2.le
  have hL := side_le_two m hm k
  have hL0 := side_pos m hm k
  have hh := apo_le_one m hm k
  have hh0 := apo_pos m hm k
  have hsX : Real.sqrt ((side m k * F - σ * D) ^ 2 + apo m k ^ 2 * D ^ 2) ≤
      side m k * F + |σ| * |D| + apo m k * |D| := by
    have e : (side m k * F - σ * D) ^ 2 + apo m k ^ 2 * D ^ 2 =
        (side m k * F - σ * D) ^ 2 + (apo m k * D) ^ 2 := by ring
    rw [e]
    refine (sqrt_two_le _ _).trans ?_
    rw [abs_mul (apo m k), abs_of_pos hh0]
    have : |side m k * F - σ * D| ≤ side m k * F + |σ| * |D| := by
      calc _ ≤ |side m k * F| + |σ * D| := abs_sub _ _
        _ = _ := by rw [abs_mul, abs_mul, abs_of_pos hL0, abs_of_pos hFpos]
    linarith
  have hS : side m k * F + |σ| * |D| + apo m k * |D| ≤ 8 * m * (M1 * W) := by
    have : side m k * F ≤ 2 * m := by nlinarith
    have : |σ| * |D| ≤ 2 * (2 * m * (M1 * W)) := by
      have := abs_nonneg D; have := abs_nonneg σ; nlinarith
    have : apo m k * |D| ≤ 1 * (2 * m * (M1 * W)) := by
      have := abs_nonneg D; nlinarith
    have : (m:ℝ) ≤ m * (M1 * W) := le_mul_of_one_le_right hm0.le (by nlinarith)
    linarith
  have hEpos := Real.rpow_nonneg hFpos.le (-1 - 1 / p)
  calc F ^ (-1 - 1 / p) * Real.sqrt ((side m k * F - σ * D) ^ 2 + apo m k ^ 2 * D ^ 2)
      ≤ (m:ℝ) ^ 3 * (8 * m * (M1 * W)) := by
        gcongr
        exact hsX.trans hS
    _ = 8 * (m:ℝ) ^ 4 * M1 * W * 1 := by ring
    _ ≤ 8 * (m:ℝ) ^ 4 * M1 * W * Λ := by
        gcongr
    _ = _ := by ring

lemma G_int (m : ℕ) (hm : 2 ≤ m) (k : ℕ) {p : ℝ} (hp : p ∈ Icc (3 / 4 : ℝ) (5 / 4)) :
    IntervalIntegrable (G m p k) volume 0 1 := by
  refine (bound_int (8 * (m : ℝ) ^ 4 * mu m ^ (-(1 / 4 : ℝ)))).mono_fun
    (G_meas m p k).aestronglyMeasurable ?_
  refine ae_Ioo (fun t ht => ?_)
  dsimp only
  rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_of_nonneg (G_nonneg m p k t)]
  exact (G_bound m hm k hp ht).trans (le_abs_self _)

lemma speedData (m : ℕ) (hm : 2 ≤ m) (k : ℕ) {p : ℝ} (hp : p ∈ Ioo (3 / 4 : ℝ) (5 / 4)) :
    SpeedData (arc m p k) (G m p k) where
  cont := arc_contOn m hm (by linarith [hp.1]) k
  der := fun t ht => by
    obtain ⟨v, hv, hn⟩ := speed m hm p k ht
    rw [hv.deriv]; exact ⟨hv, hn⟩
  int := G_int m hm k (Ioo_subset_Icc_self hp)
  nonneg := G_nonneg m p k
  gcont := fun t ht => G_contAt m hm (by linarith [hp.1] : p ≠ 0) k ht
  meas := G_meas m p k

lemma loop_eq (m : ℕ) (p : ℝ) (k : ℕ) :
    EqOn (loop m p) (fun s => arc m p k (s - k)) (Icc (k : ℝ) (k + 1)) := by
  intro s hs
  rcases eq_or_lt_of_le hs.2 with h | h
  · subst h
    have : ⌊(k : ℝ) + 1⌋₊ = k + 1 := by exact_mod_cast Nat.floor_natCast (k + 1)
    simp only [loop, this]
    push_cast
    rw [sub_self, show (k:ℝ) + 1 - k = 1 by ring]
    simp only [arc, Fe, edgePt_zero, edgePt_one]
  · have hfl : ⌊s⌋₊ = k := by
      rw [Nat.floor_eq_iff (by linarith [hs.1, (Nat.cast_nonneg k : (0:ℝ) ≤ k)])]
      exact ⟨hs.1, h⟩
    simp only [loop, hfl]

lemma var_shift' {E : Type*} [PseudoEMetricSpace E] (g : ℝ → E) (k : ℝ) :
    eVariationOn (fun s => g (s - k)) (Icc k (k+1)) = eVariationOn g (Icc 0 1) := by
  rw [show (fun s => g (s - k)) = g ∘ fun t => t - k from rfl,
    eVariationOn.comp_eq_of_monotoneOn g (fun t => t - k) (fun a _ b _ h => by simpa using h)]
  congr 1
  ext x
  simp only [mem_image, mem_Icc]
  constructor
  · rintro ⟨y, ⟨h1, h2⟩, rfl⟩; constructor <;> linarith
  · rintro ⟨h1, h2⟩; exact ⟨x + k, ⟨by linarith, by linarith⟩, by ring⟩

lemma stmt5_core (m : ℕ) (hm : 2 ≤ m) {p : ℝ} (hp : p ∈ Ioo (3 / 4 : ℝ) (5 / 4)) :
    (∀ k : ℕ, IntervalIntegrable (G m p k) volume 0 1) ∧
    eVariationOn (loop m p) (Icc 0 (2 * (m : ℝ))) =
      ENNReal.ofReal (∑ k ∈ Finset.range (2 * m), ∫ t in (0 : ℝ)..1, G m p k t) := by
  refine ⟨fun k => G_int m hm k (Ioo_subset_Icc_self hp), ?_⟩
  have h2m : (2 * (m : ℝ)) = ((2 * m : ℕ) : ℝ) := by push_cast; ring
  rw [h2m, var_sum, ENNReal.ofReal_sum_of_nonneg (fun k _ =>
    intervalIntegral.integral_nonneg zero_le_one (fun t _ => G_nonneg m p k t))]
  refine Finset.sum_congr rfl (fun k _ => ?_)
  rw [eVariationOn.congr (loop_eq m p k), var_shift', (speedData m hm k hp).var_eq]

end X04

theorem X04.proof_5 : X04.stmt_5 := fun m hm p hp => X04.stmt5_core m hm hp

namespace X04
open MeasureTheory Set Filter Real
open scoped Topology

lemma stmt13_core (m : ℕ) (hm : 2 ≤ m) :
    (∀ᶠ p in 𝓝 (1 : ℝ), eVariationOn (loop m p) (Icc 0 (2 * (m : ℝ))) ≠ ⊤) ∧
    HasDerivAt (perim m) (U03.P1 m * U03.R1 m) 1 := by
  have hs : Ioo (3 / 4 : ℝ) (5 / 4) ∈ 𝓝 (1 : ℝ) := Ioo_mem_nhds (by norm_num) (by norm_num)
  refine ⟨?_, ?_⟩
  · filter_upwards [hs] with p hp
    rw [(stmt5_core m hm hp).2]; exact ENNReal.ofReal_ne_top
  have hev : (fun p => ∑ k ∈ Finset.range (2 * m), ∫ t in (0 : ℝ)..1, G m p k t) =ᶠ[𝓝 1] perim m := by
    filter_upwards [hs] with p hp
    rw [perim, (stmt5_core m hm hp).2, ENNReal.toReal_ofReal]
    exact Finset.sum_nonneg (fun k _ =>
      intervalIntegral.integral_nonneg zero_le_one (fun t _ => G_nonneg m p k t))
  have hd := HasDerivAt.fun_sum (u := Finset.range (2 * m)) (fun k _ => stmt9_core m hm k)
  refine (hd.congr_of_eventuallyEq hev.symm).congr_deriv ?_
  have h12 := proof_12 m hm
  simp only [S1_vertex m hm]
  have : ∀ k ∈ Finset.range (2 * m), dot (tau m k) (vertex m k) * -U03.R1 m -
      dot (tau m k) (vertex m (k + 1)) * -U03.R1 m =
      -U03.R1 m * (dot (tau m k) (vertex m k) - dot (tau m k) (vertex m (k + 1))) :=
    fun k _ => by ring
  rw [Finset.sum_congr rfl this, ← Finset.mul_sum, h12]; ring

end X04

theorem X04.proof_13 : X04.stmt_13 := X04.stmt13_core

theorem X04.proof_14 : X04.stmt_14 := fun m hm =>
  ⟨_, (X04.stmt13_core m hm).2, U03.proof_5 m hm⟩

theorem X04.proof_15 : X04.stmt_15 := by
  have h := (X04.stmt13_core 3 (by norm_num)).2
  have e6 : Real.pi / (2 * ((3 : ℕ) : ℝ)) = Real.pi / 6 := by push_cast; ring
  have s3 : Real.sqrt 3 ≠ 0 := by positivity
  have hV : U03.V 3 = Real.sqrt 3 := by
    rw [U03.V, e6, Real.cot_eq_cos_div_sin, Real.cos_pi_div_six, Real.sin_pi_div_six]; ring
  have hs1 : Real.sin (((1 : ℕ) : ℝ) * Real.pi / ((3 : ℕ) : ℝ)) = Real.sqrt 3 / 2 := by
    rw [show ((1 : ℕ) : ℝ) * Real.pi / ((3 : ℕ) : ℝ) = Real.pi / 3 by push_cast; ring,
      Real.sin_pi_div_three]
  have hs2 : Real.sin (((2 : ℕ) : ℝ) * Real.pi / ((3 : ℕ) : ℝ)) = Real.sqrt 3 / 2 := by
    rw [show ((2 : ℕ) : ℝ) * Real.pi / ((3 : ℕ) : ℝ) = Real.pi - Real.pi / 3 by push_cast; ring,
      Real.sin_pi_sub, Real.sin_pi_div_three]
  have hR : U03.R1 3 = Real.log 2 := by
    rw [U03.R1, hV]
    rw [show Finset.Ico 1 3 = {1, 2} by decide, Finset.sum_pair (by norm_num), hs1, hs2]
    have : Real.sqrt 3 / (Real.sqrt 3 / 2) = 2 := by field_simp
    rw [this]; field_simp; ring
  have hP : U03.P1 3 = 2 * Real.sqrt 3 := by
    rw [U03.P1, e6, Real.cos_pi_div_six]
    have h3 : Real.sqrt 3 ^ 2 = 3 := Real.sq_sqrt (by norm_num)
    field_simp
    push_cast
    nlinarith [h3]
  unfold X04.stmt_15
  rw [hR, hP] at h; exact h

namespace X04
open MeasureTheory Set Filter Real
open scoped Topology

lemma dot_cont (u : Plane) : Continuous (fun x : Plane => dot u x) := by
  unfold dot; fun_prop

lemma Fp_cont (m : ℕ) {p : ℝ} (hp : 0 < p) : Continuous (Fp m p) := by
  unfold Fp
  exact continuous_finsetSum _ (fun i _ => (dot_cont _).abs.rpow_const (fun _ => Or.inr hp.le))

lemma Fp_smul (m : ℕ) {p : ℝ} (c : ℝ) (x : Plane) : Fp m p (c • x) = |c| ^ p * Fp m p x := by
  unfold Fp
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  have : dot (nrm m i) (c • x) = c * dot (nrm m i) x := by
    simp only [dot, PiLp.smul_apply, smul_eq_mul]; ring
  rw [this, abs_mul, Real.mul_rpow (abs_nonneg _) (abs_nonneg _)]

lemma Fp_pos (m : ℕ) (hm : 2 ≤ m) {p : ℝ} {x : Plane} (hx : x ≠ 0) : 0 < Fp m p x := by
  unfold Fp
  have hm2 : (2:ℝ) ≤ m := by exact_mod_cast hm
  have hm0 : (m:ℝ) ≠ 0 := by positivity
  have hne : ∃ i ∈ Finset.Icc 1 m, dot (nrm m i) x ≠ 0 := by
    by_cases h0 : x 0 = 0
    · have h1 : x 1 ≠ 0 := by
        intro h1; apply hx; exact plane_ext' (by simp [h0]) (by simp [h1])
      refine ⟨1, Finset.mem_Icc.2 ⟨le_rfl, by omega⟩, ?_⟩
      rw [dot, nrm0, nrm1, h0, mul_zero, zero_add]
      exact mul_ne_zero (by simpa using (sin_pi_div_pos m hm).ne') h1
    · refine ⟨m, Finset.mem_Icc.2 ⟨by omega, le_rfl⟩, ?_⟩
      rw [dot, nrm0, nrm1, show (m:ℝ) * π / m = π by field_simp, Real.cos_pi, Real.sin_pi]
      simpa using h0
  obtain ⟨i, hi, h⟩ := hne
  exact Finset.sum_pos' (fun j _ => Real.rpow_nonneg (abs_nonneg _) _)
    ⟨i, hi, Real.rpow_pos_of_pos (abs_pos.2 h) _⟩

lemma frontier_Kp (m : ℕ) (hm : 2 ≤ m) {p : ℝ} (hp : 0 < p) :
    frontier (Kp m p) = {x | Fp m p x = 1} := by
  apply subset_antisymm
  · exact frontier_le_subset_eq (Fp_cont m hp) continuous_const
  · intro x hx
    simp only [mem_ofPred_eq] at hx
    unfold Kp
    rw [frontier, (isClosed_le (Fp_cont m hp) continuous_const).closure_eq]
    refine ⟨le_of_eq hx, fun hint => ?_⟩
    have hK : {y | Fp m p y ≤ 1} ∈ 𝓝 x := mem_interior_iff_mem_nhds.1 hint
    have ht : Tendsto (fun c : ℝ => c • x) (𝓝 1) (𝓝 x) := by
      have := (tendsto_id (x := 𝓝 (1:ℝ))).smul_const x
      simpa using this
    have hev := (ht.eventually hK).filter_mono (nhdsWithin_le_nhds (s := Ioi (1:ℝ)))
    obtain ⟨c, hc1, hc2⟩ := (hev.and self_mem_nhdsWithin).exists
    simp only [Fp_smul, hx, mul_one] at hc1
    have hc' : 1 < c := hc2
    have : 1 < |c| ^ p := Real.one_lt_rpow (by rw [abs_of_pos (by linarith)]; exact hc') hp
    linarith

lemma Fp_arc (m : ℕ) (hm : 2 ≤ m) {p : ℝ} (hp : 0 < p) (k : ℕ) {t : ℝ} (ht : t ∈ Icc (0:ℝ) 1) :
    Fp m p (arc m p k t) = 1 := by
  have hF := Fe_pos_Icc m hm p k ht
  unfold arc
  rw [Fp_smul, abs_of_pos (Real.rpow_pos_of_pos hF _), ← Real.rpow_mul hF.le]
  show Fe m p k t ^ (-(1 / p) * p) * Fe m p k t = 1
  rw [show -(1 / p) * p = -1 by field_simp, Real.rpow_neg_one, inv_mul_cancel₀ hF.ne']

lemma loop_floor (m : ℕ) (p : ℝ) {s : ℝ} (hs : 0 ≤ s) :
    loop m p s = arc m p ⌊s⌋₊ (s - ⌊s⌋₊) := rfl

lemma frac_mem {s : ℝ} (hs : 0 ≤ s) : s - ⌊s⌋₊ ∈ Ico (0:ℝ) 1 :=
  ⟨by linarith [Nat.floor_le hs], by linarith [Nat.lt_floor_add_one s]⟩

lemma loop_contOn (m : ℕ) (hm : 2 ≤ m) {p : ℝ} (hp : 0 < p) :
    ContinuousOn (loop m p) (Icc 0 (2 * (m : ℝ))) := by
  have piece : ∀ k : ℕ, ContinuousOn (loop m p) (Icc (k : ℝ) (k + 1)) := by
    intro k
    have h1 : ContinuousOn (fun s => arc m p k (s - k)) (Icc (k : ℝ) (k + 1)) := by
      refine (arc_contOn m hm hp k).comp (continuousOn_id.sub continuousOn_const) ?_
      intro s hs; exact ⟨by simp; linarith [hs.1], by simp; linarith [hs.2]⟩
    exact h1.congr (loop_eq m p k)
  have hn : ∀ n : ℕ, ContinuousOn (loop m p) (Icc 0 (n : ℝ)) := by
    intro n
    induction n with
    | zero =>
      simp only [CharP.cast_eq_zero, Icc_self]
      exact continuousOn_singleton _ _
    | succ n ih =>
      rw [show ((n + 1 : ℕ) : ℝ) = n + 1 by push_cast; ring,
        ← Icc_union_Icc_eq_Icc (Nat.cast_nonneg n) (by linarith)]
      exact ih.union_of_isClosed (piece n) isClosed_Icc isClosed_Icc
  have := hn (2 * m)
  push_cast at this
  exact this

lemma loop_closed (m : ℕ) (hm : 2 ≤ m) (p : ℝ) :
    loop m p (2 * (m : ℝ)) = loop m p 0 := by
  have h1 : ⌊(2 * (m : ℝ))⌋₊ = 2 * m := by exact_mod_cast Nat.floor_natCast (2 * m)
  simp only [loop, h1, Nat.floor_zero, CharP.cast_eq_zero, sub_zero]
  push_cast
  rw [sub_self]
  simp only [arc, Fe, edgePt_zero]
  rw [show 2 * m = 0 + 2 * m by ring, vertex_add_two_m m hm 0]

lemma loop_image (m : ℕ) (hm : 2 ≤ m) {p : ℝ} (hp : 0 < p) :
    loop m p '' Icc 0 (2 * (m : ℝ)) = {x | Fp m p x = 1} := by
  apply subset_antisymm
  · rintro _ ⟨s, hs, rfl⟩
    rw [loop_floor m p hs.1]
    exact Fp_arc m hm hp _ (Ico_subset_Icc_self (frac_mem hs.1))
  · intro x hx
    simp only [mem_ofPred_eq] at hx
    have hx0 : x ≠ 0 := by
      rintro rfl
      have : Fp m p (0 : Plane) = 0 := by
        have := Fp_smul m (p := p) 0 (0 : Plane)
        rw [smul_zero, abs_zero, Real.zero_rpow hp.ne', zero_mul] at this; exact this
      linarith
    obtain ⟨k, hk, a, b, ha, hb, rfl⟩ := sector m hm x
    have hab : 0 < a + b := by
      rcases eq_or_lt_of_le (add_nonneg ha hb) with h | h
      · exfalso; apply hx0
        have ha0 : a = 0 := by linarith
        have hb0 : b = 0 := by linarith
        simp [ha0, hb0]
      · exact h
    set t := b / (a + b) with htdef
    have ht0 : 0 ≤ t := div_nonneg hb hab.le
    have ht1 : t ≤ 1 := by rw [div_le_one hab]; linarith
    have hy : edgePt m k t = (a + b)⁻¹ • (a • vertex m k + b • vertex m (k + 1)) := by
      apply plane_ext'
      · simp only [edgePt, PiLp.add_apply, PiLp.smul_apply, PiLp.sub_apply, smul_eq_mul, htdef]
        field_simp; ring
      · simp only [edgePt, PiLp.add_apply, PiLp.smul_apply, PiLp.sub_apply, smul_eq_mul, htdef]
        field_simp; ring
    have hFe : Fe m p k t = ((a + b)⁻¹) ^ p := by
      rw [Fe, hy, Fp_smul, hx, mul_one, abs_of_pos (inv_pos.2 hab)]
    refine ⟨k + t, ⟨by positivity, ?_⟩, ?_⟩
    · have : (k:ℝ) + 1 ≤ 2 * m := by exact_mod_cast hk
      linarith
    · rw [loop_eq m p k ⟨by linarith, by linarith⟩]
      dsimp only
      rw [add_sub_cancel_left, arc, hFe, hy,
        ← Real.rpow_mul (inv_pos.2 hab).le, show p * -(1 / p) = -1 by field_simp,
        Real.rpow_neg_one, inv_inv, smul_smul, mul_inv_cancel₀ hab.ne', one_smul]

end X04

namespace X04
open MeasureTheory Set Filter Real
open scoped Topology

/-- the unit vector `(-sin θ, cos θ)`. -/
noncomputable def eU (θ : ℝ) : Plane := !₂[-Real.sin θ, Real.cos θ]

/-- the angle of `y(t)` measured from `C_k`. -/
noncomputable def phi (m : ℕ) (t : ℝ) : ℝ :=
  Real.arctan (t * Real.sin (π / m) / (1 - t + t * Real.cos (π / m)))

lemma beta_facts (m : ℕ) (hm : 2 ≤ m) :
    0 < π / m ∧ π / m ≤ π / 2 ∧ 0 < Real.sin (π / m) ∧ 0 ≤ Real.cos (π / m) := by
  have hm2 : (2:ℝ) ≤ m := by exact_mod_cast hm
  have h1 : 0 < π / m := by positivity
  have h2 : π / m ≤ π / 2 := by
    rw [div_le_div_iff₀ (by positivity) (by norm_num)]; nlinarith [Real.pi_pos]
  exact ⟨h1, h2, sin_pi_div_pos m hm,
    Real.cos_nonneg_of_mem_Icc ⟨by linarith [Real.pi_pos], h2⟩⟩

lemma R_pos (m : ℕ) (hm : 2 ≤ m) {t : ℝ} (t1 : t < 1) (t0 : 0 ≤ t) :
    0 < 1 - t + t * Real.cos (π / m) := by
  have := (beta_facts m hm).2.2.2
  nlinarith

lemma polar_RI {R I : ℝ} (hR : 0 < R) :
    Real.sqrt (R ^ 2 + I ^ 2) * Real.cos (Real.arctan (I / R)) = R ∧
    Real.sqrt (R ^ 2 + I ^ 2) * Real.sin (Real.arctan (I / R)) = I := by
  have hs : Real.sqrt (1 + (I / R) ^ 2) = Real.sqrt (R ^ 2 + I ^ 2) / R := by
    rw [show 1 + (I / R) ^ 2 = (R ^ 2 + I ^ 2) / R ^ 2 by field_simp, Real.sqrt_div' _ (sq_nonneg R),
      Real.sqrt_sq hR.le]
  have hpos : 0 < Real.sqrt (R ^ 2 + I ^ 2) := Real.sqrt_pos.2 (by positivity)
  rw [Real.cos_arctan, Real.sin_arctan, hs]
  constructor <;> field_simp

lemma phi_nonneg (m : ℕ) (hm : 2 ≤ m) {t : ℝ} (t0 : 0 ≤ t) (t1 : t < 1) : 0 ≤ phi m t := by
  unfold phi
  rw [Real.arctan_nonneg]
  have := R_pos m hm t1 t0
  have := (beta_facts m hm).2.2.1
  positivity

lemma phi_lt_half (m : ℕ) (t : ℝ) : phi m t < π / 2 := Real.arctan_lt_pi_div_two _

lemma phi_lt (m : ℕ) (hm : 2 ≤ m) {t : ℝ} (t0 : 0 ≤ t) (t1 : t < 1) : phi m t < π / m := by
  obtain ⟨b0, b1, sb, cb⟩ := beta_facts m hm
  have hR := R_pos m hm t1 t0
  set R := 1 - t + t * Real.cos (π / m)
  set I := t * Real.sin (π / m)
  obtain ⟨hc, hsn⟩ := polar_RI (I := I) hR
  have hpos : 0 < Real.sqrt (R ^ 2 + I ^ 2) := Real.sqrt_pos.2 (by positivity)
  have hsin : 0 < Real.sin (π / m - phi m t) := by
    rw [Real.sin_sub]
    have : Real.sqrt (R ^ 2 + I ^ 2) * (Real.sin (π / m) * Real.cos (phi m t) -
        Real.cos (π / m) * Real.sin (phi m t)) = Real.sin (π / m) * (1 - t) := by
      unfold phi
      rw [mul_sub, ← mul_assoc, ← mul_assoc, mul_comm (Real.sqrt _) (Real.sin _),
        mul_comm (Real.sqrt _) (Real.cos _), mul_assoc, mul_assoc, hc, hsn]
      simp only [R, I]; ring
    have h2 : 0 < Real.sin (π / m) * (1 - t) := mul_pos sb (by linarith)
    rw [← this] at h2
    exact pos_of_mul_pos_right h2 hpos.le
  by_contra hle
  push Not at hle
  have := Real.sin_nonpos_of_nonpos_of_neg_pi_le (by linarith : π / m - phi m t ≤ 0)
    (by linarith [phi_nonneg m hm t0 t1, phi_lt_half m t, Real.pi_pos])
  linarith

lemma phi_inj (m : ℕ) (hm : 2 ≤ m) {t t' : ℝ} (t0 : 0 ≤ t) (t1 : t < 1) (t0' : 0 ≤ t') (t1' : t' < 1)
    (h : phi m t = phi m t') : t = t' := by
  obtain ⟨b0, b1, sb, cb⟩ := beta_facts m hm
  have hR := R_pos m hm t1 t0
  have hR' := R_pos m hm t1' t0'
  have := Real.arctan_injective h
  rw [div_eq_div_iff hR.ne' hR'.ne'] at this
  have e : Real.sin (π / m) * (t - t') = 0 := by linear_combination this
  rcases mul_eq_zero.1 e with h1 | h1
  · linarith
  · linarith

lemma edge_polar (m : ℕ) (hm : 2 ≤ m) (k : ℕ) {t : ℝ} (t0 : 0 ≤ t) (t1 : t < 1) :
    edgePt m k t = (Real.tan (π / (2 * (m : ℝ))) *
      Real.sqrt ((1 - t + t * Real.cos (π / m)) ^ 2 + (t * Real.sin (π / m)) ^ 2)) •
        eU ((k : ℝ) * (π / m) + phi m t) := by
  have hR := R_pos m hm t1 t0
  obtain ⟨hc, hsn⟩ := polar_RI (I := t * Real.sin (π / m)) hR
  have hk : (k : ℝ) * π / m = (k : ℝ) * (π / m) := by ring
  have hk1 : ((k + 1 : ℕ) : ℝ) * π / m = (k : ℝ) * (π / m) + π / m := by push_cast; ring
  set a := (k : ℝ) * (π / m)
  set S := Real.sqrt ((1 - t + t * Real.cos (π / m)) ^ 2 + (t * Real.sin (π / m)) ^ 2)
  apply plane_ext'
  · simp only [edgePt, PiLp.add_apply, PiLp.smul_apply, PiLp.sub_apply, smul_eq_mul, vertex0, eU]
    simp only [Matrix.cons_val_zero]
    rw [hk, hk1, Real.sin_add, Real.sin_add]
    unfold phi
    linear_combination (Real.tan (π / (2 * (m : ℝ))) * Real.sin a) * hc +
      (Real.tan (π / (2 * (m : ℝ))) * Real.cos a) * hsn
  · simp only [edgePt, PiLp.add_apply, PiLp.smul_apply, PiLp.sub_apply, smul_eq_mul, vertex1, eU]
    simp only [Matrix.cons_val_one, Matrix.cons_val_zero]
    rw [hk, hk1, Real.cos_add, Real.cos_add]
    unfold phi
    linear_combination (-Real.tan (π / (2 * (m : ℝ))) * Real.cos a) * hc +
      (Real.tan (π / (2 * (m : ℝ))) * Real.sin a) * hsn

lemma eU_inj {c c' θ θ' : ℝ} (hc : 0 < c) (hc' : 0 < c') (h0 : 0 ≤ θ) (h1 : θ < 2 * π)
    (h0' : 0 ≤ θ') (h1' : θ' < 2 * π) (h : c • eU θ = c' • eU θ') : θ = θ' := by
  have e0 := congrArg (fun v : Plane => v 0) h
  have e1 := congrArg (fun v : Plane => v 1) h
  simp only [eU, PiLp.smul_apply, smul_eq_mul, Matrix.cons_val_zero,
    Matrix.cons_val_one] at e0 e1
  have hsq : c ^ 2 = c' ^ 2 := by
    have a := Real.sin_sq_add_cos_sq θ
    have b := Real.sin_sq_add_cos_sq θ'
    have ea : c ^ 2 = (c * -Real.sin θ) ^ 2 + (c * Real.cos θ) ^ 2 := by linear_combination (-c ^ 2) * a
    have eb : c' ^ 2 = (c' * -Real.sin θ') ^ 2 + (c' * Real.cos θ') ^ 2 := by
      linear_combination (-c' ^ 2) * b
    rw [ea, eb, e0, e1]
  have hcc : c = c' := by nlinarith
  subst hcc
  have hs : Real.sin θ = Real.sin θ' := by
    have := mul_left_cancel₀ hc.ne' e0; linarith
  have hco : Real.cos θ = Real.cos θ' := mul_left_cancel₀ hc.ne' e1
  have hcos : Real.cos (θ - θ') = 1 := by
    rw [Real.cos_sub, ← hs, ← hco]; nlinarith [Real.sin_sq_add_cos_sq θ]
  have := (Real.cos_eq_one_iff_of_lt_of_lt (by linarith) (by linarith)).1 hcos
  linarith

lemma loop_polar (m : ℕ) (hm : 2 ≤ m) {p : ℝ} {s : ℝ} (hs : s ∈ Ico 0 (2 * (m : ℝ))) :
    ∃ c : ℝ, 0 < c ∧ loop m p s = c • eU ((⌊s⌋₊ : ℝ) * (π / m) + phi m (s - ⌊s⌋₊)) := by
  obtain ⟨t0, t1⟩ := frac_mem hs.1
  have hF := Fe_pos_Icc m hm p ⌊s⌋₊ ⟨t0, t1.le⟩
  have hR := R_pos m hm t1 t0
  refine ⟨Fe m p ⌊s⌋₊ (s - ⌊s⌋₊) ^ (-(1 / p)) * (Real.tan (π / (2 * (m : ℝ))) *
      Real.sqrt ((1 - (s - ⌊s⌋₊) + (s - ⌊s⌋₊) * Real.cos (π / m)) ^ 2 +
        ((s - ⌊s⌋₊) * Real.sin (π / m)) ^ 2)), ?_, ?_⟩
  · have := rho_pos m hm
    have : 0 < Real.sqrt ((1 - (s - ⌊s⌋₊) + (s - ⌊s⌋₊) * Real.cos (π / m)) ^ 2 +
        ((s - ⌊s⌋₊) * Real.sin (π / m)) ^ 2) := Real.sqrt_pos.2 (by positivity)
    have := Real.rpow_pos_of_pos hF (-(1 / p))
    positivity
  · rw [loop_floor m p hs.1, arc, edge_polar m hm _ t0 t1, smul_smul]

lemma loop_inj (m : ℕ) (hm : 2 ≤ m) (p : ℝ) : InjOn (loop m p) (Ico 0 (2 * (m : ℝ))) := by
  intro s hs s' hs' heq
  obtain ⟨c, hc, hcs⟩ := loop_polar m hm (p := p) hs
  obtain ⟨c', hc', hcs'⟩ := loop_polar m hm (p := p) hs'
  obtain ⟨t0, t1⟩ := frac_mem hs.1
  obtain ⟨t0', t1'⟩ := frac_mem hs'.1
  obtain ⟨b0, -, -, -⟩ := beta_facts m hm
  set k := ⌊s⌋₊ with hkdef
  set k' := ⌊s'⌋₊ with hk'def
  have hm0 : (0:ℝ) < m := by have : (2:ℝ) ≤ m := by exact_mod_cast hm
                             linarith
  have hk : (k : ℝ) + 1 ≤ 2 * m := by
    have := Nat.floor_le hs.1
    have h2 : (k : ℝ) < 2 * m := by linarith [hs.2]
    have : k < 2 * m := by exact_mod_cast h2
    exact_mod_cast this
  have hk' : (k' : ℝ) + 1 ≤ 2 * m := by
    have := Nat.floor_le hs'.1
    have h2 : (k' : ℝ) < 2 * m := by linarith [hs'.2]
    have : k' < 2 * m := by exact_mod_cast h2
    exact_mod_cast this
  have hpi : 2 * (m : ℝ) * (π / m) = 2 * π := by field_simp
  have φ1 := phi_lt m hm t0 t1
  have φ1' := phi_lt m hm t0' t1'
  have φ0 := phi_nonneg m hm t0 t1
  have φ0' := phi_nonneg m hm t0' t1'
  have hθ1 : (k : ℝ) * (π / m) + phi m (s - k) < 2 * π := by
    have : ((k : ℝ) + 1) * (π / m) ≤ 2 * m * (π / m) := mul_le_mul_of_nonneg_right hk b0.le
    nlinarith
  have hθ1' : (k' : ℝ) * (π / m) + phi m (s' - k') < 2 * π := by
    have : ((k' : ℝ) + 1) * (π / m) ≤ 2 * m * (π / m) := mul_le_mul_of_nonneg_right hk' b0.le
    nlinarith
  rw [hcs, hcs'] at heq
  have hθ := eU_inj hc hc' (by positivity) hθ1 (by positivity) hθ1' heq
  have hkk : k = k' := by
    by_contra hne
    rcases lt_or_gt_of_ne hne with h | h
    · have : (k : ℝ) + 1 ≤ k' := by exact_mod_cast h
      have : ((k : ℝ) + 1) * (π / m) ≤ k' * (π / m) := mul_le_mul_of_nonneg_right this b0.le
      nlinarith
    · have : (k' : ℝ) + 1 ≤ k := by exact_mod_cast h
      have : ((k' : ℝ) + 1) * (π / m) ≤ k * (π / m) := mul_le_mul_of_nonneg_right this b0.le
      nlinarith
  rw [← hkk] at hθ t0' t1'
  have hphi : phi m (s - k) = phi m (s' - k) := by linarith
  have := phi_inj m hm t0 t1 t0' t1' hphi
  linarith

end X04

theorem X04.proof_3 : X04.stmt_3 := by
  intro m hm p hp
  refine ⟨X04.loop_contOn m hm hp, X04.loop_closed m hm p, X04.loop_inj m hm p, ?_,
    X04.frontier_Kp m hm hp⟩
  rw [X04.frontier_Kp m hm hp, X04.loop_image m hm hp]
