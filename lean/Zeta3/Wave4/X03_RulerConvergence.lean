import Mathlib
/-!
# X03 -- every fixed ruler `q ≥ 1` gives a convergent expansion of `L_q(p)` in powers of `1/p`

Target id: X03.
Paper label: `prop:ruler`, last sentence (paper/ieee/main.tex, lines 1384-1390).
Detailed proof: paper/ieee/proofs/sec_ruler.tex, Lemma `lem:ruler:formF` (lines 48-54),
Lemma `lem:ruler:holo` (lines 56-76), "Proof of Proposition~4" (lines 78-84).
The two exact rulers `L_1 = 8`, `L_∞ = 8·2^{-1/p}` (first sentence of `prop:ruler`) were done in
Wave1/T14_Rulers.lean and are not repeated here. This file does not import T14; `Lq` below is
the same literal formula, restated.

Quoted verbatim (main.tex 1384):

  let $L_q(p)=\frac4p\int_0^1\bigl(t^{q/p-q}+(1-t)^{q/p-q}\bigr)^{1/q}\dd t$ be the length of the
  $\ell^p$ circle measured with the $\ell^q$ norm, so that $L_2=P$ and $L_p=2\pi_p$.

Quoted verbatim (main.tex 1386-1388, the claim formalized is the last sentence):

  \begin{proposition}[Two exact rulers$^\dagger$]\label{prop:ruler}
  For every $p\ge1$, $L_1(p)=8$ and $L_\infty(p)=8\cdot2^{-1/p}$. For every fixed $q\ge1$, the
  expansion of $L_q(p)$ in powers of $1/p$ converges for every $p>1$.
  \end{proposition}

Quoted verbatim (main.tex 1390):

  and for fixed $q$ that of Lemma~\ref{lem:holo} works with $z=(b_\varepsilon/a_\varepsilon)^q$ and
  $(1+z)^{1/q}$ in place of $\sqrt{1+z}$: the integrand is $b_\varepsilon-a_\varepsilon z\,g(z)$,
  where $g(z)=((1+z)^{1/q}-1)/z$ is bounded on $|z|\le1$ and
  $|a_\varepsilon z|=|b_\varepsilon|^q|a_\varepsilon|^{1-q}\le|b_\varepsilon|$ because $q\ge1$.

Quoted verbatim (sec_ruler.tex 44-46, definitions):

  For $\operatorname{Re}\varepsilon<1$ and $0<t<1$ put $a_\varepsilon=t^{\varepsilon-1}=\exp((\varepsilon-1)\log t)$
  and $b_\varepsilon=(1-t)^{\varepsilon-1}$, and $a_\varepsilon^q:=t^{q(\varepsilon-1)}$,
  $b_\varepsilon^q:=(1-t)^{q(\varepsilon-1)}$. For real $\varepsilon\in(0,1]$ define
  F_q(\varepsilon)=\int_0^1\Bigl(a_\varepsilon+b_\varepsilon-(a_\varepsilon^q+b_\varepsilon^q)^{1/q}\Bigr)dt .

Quoted verbatim (sec_ruler.tex 49, Lemma `lem:ruler:formF`):

  For $0<\varepsilon\le1$: $L_q(1/\varepsilon)=8-4\varepsilon F_q(\varepsilon)$, that is,
  $L_q(p)/8=1-\varepsilon F_q(\varepsilon)/2$. The integral $F_q(\varepsilon)$ converges absolutely.

Quoted verbatim (sec_ruler.tex 57, Lemma `lem:ruler:holo`):

  $F_q$ extends to a holomorphic function on $\{\operatorname{Re}\varepsilon<1\}$.

Quoted verbatim (sec_ruler.tex 61-69, from the proof of `lem:ruler:holo`):

  $F_q(\varepsilon)=2\int_0^{1/2}f_\varepsilon(t)\,dt$ ...
  z=z(t,\varepsilon)=\frac{b_\varepsilon^q}{a_\varepsilon^q}=\Bigl(\frac t{1-t}\Bigr)^{q(1-\varepsilon)}
  ...
  f_\varepsilon(t)=a_\varepsilon+b_\varepsilon-a_\varepsilon(1+z)^{1/q}=b_\varepsilon-a_\varepsilon\bigl((1+z)^{1/q}-1\bigr)=b_\varepsilon-a_\varepsilon z\,g(z),
  with the principal branch $(1+z)^{1/q}=\sum_{n\ge0}\binom{1/q}nz^n$ on $|z|<1$ ...
  so $|g(z)|\le1$ on $|z|\le1$ ... Hence $|f_\varepsilon(t)|\le2|b_\varepsilon|\le2^{2-\operatorname{Re}\varepsilon}$.

Quoted verbatim (sec_ruler.tex 80-84, proof of Proposition 4):

  put $M_q(\varepsilon)=-F_q(\varepsilon)/2$. ... So its Taylor series
  $\sum_{k\ge0}K_k(q)\varepsilon^k$ at $0$ converges to $M_q(\varepsilon)$ for $|\varepsilon|<1$.
  ... Hence
  L_q(p)=8\Bigl(1+\sum_{k\ge0}\frac{K_k(q)}{p^{k+1}}\Bigr)\qquad\text{for all }p>1,
  the series converges, and, because $M_q$ is holomorphic at $0$, it is the asymptotic expansion
  of $L_q(p)$ as $p\to\infty$

## Conventions

* `q : ℝ` with `1 ≤ q` throughout (`q = ∞` is T14). `p : ℝ`, `ε = 1/p`.
* `Lq p q` is the paper's integral, literally, with `Real.rpow` (bases `t, 1-t ∈ [0,1]` are
  nonnegative; `0 ^ (negative) = 0` only affects the endpoints, a null set).
* Real side: `gR q ε t = t^(ε-1) + (1-t)^(ε-1) - (t^(q(ε-1)) + (1-t)^(q(ε-1)))^(1/q)`, using the
  supplement's `a^q := t^{q(ε-1)}`, and `FR q ε = ∫_0^1 gR q ε t`.
* Complex side (proof of lem:ruler:holo): for `ε : ℂ` and real `t ∈ (0,1/2]`,
  `aC ε t = t^(ε-1)`, `bC ε t = (1-t)^(ε-1)` (complex `cpow` of a positive real base, i.e.
  `exp((ε-1) log t)`), `zC q ε t = (t/(1-t))^(q(1-ε))`, and
  `gC q ε t = aC + bC - aC · (1+zC)^(1/q)` with Mathlib's principal-branch `cpow`. On `|z| < 1`
  this is the binomial series of the supplement (`1+z` lies in the right half plane), and at
  `t = 1/2` one has `z = 1`, so the branch matters only where the supplement also uses it.
  The extension is `FC q ε = 2 ∫_0^{1/2} gC q ε t dt`.
* `K_k(q)` (the supplement's coefficients, Taylor coefficients of `M_q = -F_q/2` at `0`) is
  `Kk q k = -(iteratedDeriv k (FC q) 0) / (2 k!)`, complex derivatives of the extension.
  Indexing: `L_q(p)/8 = 1 + Σ_{k≥0} K_k(q) p^{-(k+1)}`.
* Every interval integral whose value is asserted comes with an explicit `IntervalIntegrable`
  claim, so nothing holds through Lean's "non-integrable ↦ 0" convention. Sums are `HasSum`.

## Statements

* stmt_1 (lem:ruler:formF): for `q ≥ 1`, `p ≥ 1`, `ε = 1/p`: the integrand of `L_q` and the
  integrand of `F_q` are integrable on `[0,1]` and `L_q(p) = 8 - 4 ε F_q(ε)`.
* stmt_2 (proof of lem:ruler:holo, "|g(z)| ≤ 1 on |z| ≤ 1"): for `q ≥ 1` and `‖z‖ ≤ 1`,
  `‖(1+z)^(1/q) - 1‖ ≤ ‖z‖` (principal branch).
* stmt_3 (proof of lem:ruler:holo, the bound): `‖gC q ε t‖ ≤ 2^(2 - Re ε)` for `q ≥ 1`,
  `Re ε < 1`, `0 < t ≤ 1/2`.
* stmt_4 (proof of lem:ruler:holo, "For real ε ∈ (0,1) this is the F_q of Lemma formF", and
  `F_q = 2∫_0^{1/2}`): for `q ≥ 1` and real `ε < 1`, the real integrand is integrable on `[0,1]`
  and `FC q ε = FR q ε`.
* stmt_5 (lem:ruler:holo): for `q ≥ 1`, `gC q ε` is integrable on `[0,1/2]` when `Re ε < 1`,
  and `FC q` is holomorphic (`DifferentiableOn ℂ`) on `{Re ε < 1}`.
* stmt_6 (proof of Prop. 4, "it is the asymptotic expansion of L_q(p) as p → ∞"): for `q ≥ 1`
  and every `K`, `L_q(p)/8 - 1 - Σ_{k<K} K_k(q)/p^(k+1) = O(p^{-(K+1)})` as `p → ∞`.
* stmt_7 (prop:ruler, last sentence; the final claim): for every `q ≥ 1` and every real `p > 1`,
  `Σ_{k≥0} K_k(q)/p^(k+1)` converges (`HasSum`) to `L_q(p)/8 - 1`.

## NOT COVERED

* The first sentence of `prop:ruler` (`L_1(p) = 8`, `L_∞(p) = 8·2^{-1/p}`): done in T14.
* Lemma `lem:ruler:formula` for general `1 < q < ∞` (that `L_q(p)` is the `ℓ^q` length of the
  curve as a supremum over partitions): this file takes `L_q(p)` as the paper's integral.
* The sign/summability facts about `binom(1/q, n)` used to prove stmt_2 are proof method, not
  stated separately; likewise "uniformly on compact subsets", Morera and Fubini.
* That the coefficients `K_k(q)` are real is not stated separately.
* The remarks after the proposition (rulers depending analytically on `1/p`, the dual ruler,
  linear growth, radius exactly `1`) are other targets or not claimed by `prop:ruler`.
-/
-- BEGIN STATEMENTS
namespace X03

open MeasureTheory Real Filter Asymptotics

/-- The paper's `L_q(p)=\frac4p\int_0^1(t^{q/p-q}+(1-t)^{q/p-q})^{1/q}\,dt` (main.tex 1384). -/
noncomputable def Lq (p q : ℝ) : ℝ :=
  (4 / p) * ∫ t in (0 : ℝ)..1, (t ^ (q / p - q) + (1 - t) ^ (q / p - q)) ^ (1 / q)

/-- The integrand of `L_q(p)`. -/
noncomputable def lqIntegrand (p q t : ℝ) : ℝ :=
  (t ^ (q / p - q) + (1 - t) ^ (q / p - q)) ^ (1 / q)

/-- The integrand of `F_q` for real `ε`:
`a + b - (a^q + b^q)^(1/q)` with `a^q := t^{q(ε-1)}`, `b^q := (1-t)^{q(ε-1)}`. -/
noncomputable def gR (q ε t : ℝ) : ℝ :=
  t ^ (ε - 1) + (1 - t) ^ (ε - 1) - (t ^ (q * (ε - 1)) + (1 - t) ^ (q * (ε - 1))) ^ (1 / q)

/-- `F_q(ε) = ∫_0^1 (a + b - (a^q + b^q)^(1/q)) dt` (sec_ruler.tex 46), real `ε`. -/
noncomputable def FR (q ε : ℝ) : ℝ := ∫ t in (0 : ℝ)..1, gR q ε t

/-- `a_ε(t) = t^(ε-1)` for complex `ε` (cpow of the positive real `t`). -/
noncomputable def aC (ε : ℂ) (t : ℝ) : ℂ := (t : ℂ) ^ (ε - 1)

/-- `b_ε(t) = (1-t)^(ε-1)` for complex `ε`. -/
noncomputable def bC (ε : ℂ) (t : ℝ) : ℂ := ((1 - t : ℝ) : ℂ) ^ (ε - 1)

/-- `z = (t/(1-t))^(q(1-ε))` (sec_ruler.tex 62). -/
noncomputable def zC (q : ℝ) (ε : ℂ) (t : ℝ) : ℂ :=
  ((t / (1 - t) : ℝ) : ℂ) ^ ((q : ℂ) * (1 - ε))

/-- The complex integrand `f_ε(t) = a + b - a (1+z)^(1/q)`, principal branch. -/
noncomputable def gC (q : ℝ) (ε : ℂ) (t : ℝ) : ℂ :=
  aC ε t + bC ε t - aC ε t * (1 + zC q ε t) ^ (1 / (q : ℂ))

/-- The extension from the proof of lem:ruler:holo: `F_q(ε) = 2 ∫_0^{1/2} f_ε(t) dt`. -/
noncomputable def FC (q : ℝ) (ε : ℂ) : ℂ := 2 * ∫ t in (0 : ℝ)..(1 / 2), gC q ε t

/-- `K_k(q)`: the `k`-th Taylor coefficient at `0` of `M_q = -F_q/2`,
i.e. `-F_q^{(k)}(0) / (2 k!)`. -/
noncomputable def Kk (q : ℝ) (k : ℕ) : ℂ :=
  -(iteratedDeriv k (FC q) 0) / (2 * (Nat.factorial k : ℂ))

/-- lem:ruler:formF: "For $0<\varepsilon\le1$: $L_q(1/\varepsilon)=8-4\varepsilon F_q(\varepsilon)$
... The integral $F_q(\varepsilon)$ converges absolutely", for `q ≥ 1`, `p = 1/ε ≥ 1`, with the
integrability of both integrands. -/
def stmt_1 : Prop :=
  ∀ q : ℝ, 1 ≤ q → ∀ p : ℝ, 1 ≤ p →
    IntervalIntegrable (lqIntegrand p q) volume 0 1 ∧
      IntervalIntegrable (gR q (1 / p)) volume 0 1 ∧
      Lq p q = 8 - 4 * (1 / p) * FR q (1 / p)

/-- Proof of lem:ruler:holo: "$|g(z)|\le1$ on $|z|\le1$", `g(z) = ((1+z)^{1/q}-1)/z`, stated
without division as `‖(1+z)^(1/q) - 1‖ ≤ ‖z‖` (principal branch), for `q ≥ 1`. -/
def stmt_2 : Prop :=
  ∀ q : ℝ, 1 ≤ q → ∀ z : ℂ, ‖z‖ ≤ 1 → ‖(1 + z) ^ (1 / (q : ℂ)) - 1‖ ≤ ‖z‖

/-- Proof of lem:ruler:holo: "$|f_\varepsilon(t)|\le2|b_\varepsilon|\le2^{2-\operatorname{Re}\varepsilon}$"
for `q ≥ 1`, `Re ε < 1`, `0 < t ≤ 1/2`. -/
def stmt_3 : Prop :=
  ∀ q : ℝ, 1 ≤ q → ∀ ε : ℂ, ε.re < 1 → ∀ t : ℝ, 0 < t → t ≤ 1 / 2 →
    ‖gC q ε t‖ ≤ (2 : ℝ) ^ (2 - ε.re)

/-- Proof of lem:ruler:holo: "For real $\varepsilon\in(0,1)$ this is the $F_q$ of
Lemma formF" (with `F_q = 2∫_0^{1/2}` by symmetry): for `q ≥ 1` and real `ε < 1` the real
integrand is integrable on `[0,1]` and `FC q ε = FR q ε`. -/
def stmt_4 : Prop :=
  ∀ q : ℝ, 1 ≤ q → ∀ ε : ℝ, ε < 1 →
    IntervalIntegrable (gR q ε) volume 0 1 ∧ FC q (ε : ℂ) = ((FR q ε : ℝ) : ℂ)

/-- lem:ruler:holo: "$F_q$ extends to a holomorphic function on
$\{\operatorname{Re}\varepsilon<1\}$", with the integrability of the complex integrand. -/
def stmt_5 : Prop :=
  ∀ q : ℝ, 1 ≤ q →
    (∀ ε : ℂ, ε.re < 1 → IntervalIntegrable (gC q ε) volume 0 (1 / 2)) ∧
      DifferentiableOn ℂ (FC q) {ε : ℂ | ε.re < 1}

/-- Proof of Prop. 4: "because $M_q$ is holomorphic at $0$, it is the asymptotic expansion of
$L_q(p)$ as $p\to\infty$": for every `K`,
`L_q(p)/8 - 1 - Σ_{k<K} K_k(q)/p^(k+1) = O(p^{-(K+1)})`. -/
def stmt_6 : Prop :=
  ∀ q : ℝ, 1 ≤ q → ∀ K : ℕ,
    (fun p : ℝ => ((Lq p q / 8 - 1 : ℝ) : ℂ) -
        ∑ k ∈ Finset.range K, Kk q k / (p : ℂ) ^ (k + 1))
      =O[atTop] (fun p : ℝ => (1 / p) ^ (K + 1))

/-- prop:ruler, last sentence: "For every fixed $q\ge1$, the expansion of $L_q(p)$ in powers of
$1/p$ converges for every $p>1$", in the supplement's form
`L_q(p) = 8(1 + Σ_{k≥0} K_k(q)/p^{k+1})`, stated as `HasSum` to `L_q(p)/8 - 1`. -/
def stmt_7 : Prop :=
  ∀ q : ℝ, 1 ≤ q → ∀ p : ℝ, 1 < p →
    HasSum (fun k : ℕ => Kk q k / (p : ℂ) ^ (k + 1)) ((Lq p q / 8 - 1 : ℝ) : ℂ)

end X03
-- END STATEMENTS

open MeasureTheory Set Filter Asymptotics

namespace X03

/-! ### Statement 2: `‖(1+z)^α - 1‖ ≤ ‖z‖` for `0 < α ≤ 1`, `‖z‖ ≤ 1`.

Instead of the binomial series we compare `s ↦ (1+sz)^α - 1` with the real function
`s ↦ 1 - (1 - s‖z‖)^α` on `[0,1]` (derivative bound), and use `(1-r)^α ≥ 1-r`. -/

lemma re_one_add_smul_nonneg {z : ℂ} {s : ℝ} (hs0 : 0 ≤ s) :
    1 - s * ‖z‖ ≤ (1 + (s : ℂ) * z).re := by
  have h := Complex.abs_re_le_norm z
  have : ((s : ℂ) * z).re = s * z.re := by simp
  simp only [Complex.add_re, Complex.one_re, this]
  nlinarith [neg_abs_le z.re]

lemma binom_bound {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1) {z : ℂ} (hz : ‖z‖ ≤ 1) :
    ‖(1 + z) ^ (α : ℂ) - 1‖ ≤ ‖z‖ := by
  set r := ‖z‖ with hr
  have hr0 : 0 ≤ r := norm_nonneg z
  let f : ℝ → ℂ := fun s => (1 + (s : ℂ) * z) ^ (α : ℂ) - 1
  let f' : ℝ → ℂ := fun s => (α : ℂ) * (1 + (s : ℂ) * z) ^ ((α : ℂ) - 1) * z
  let B : ℝ → ℝ := fun s => 1 - (1 - s * r) ^ α
  let B' : ℝ → ℝ := fun s => α * (1 - s * r) ^ (α - 1) * r
  have hcont : ContinuousOn f (Icc 0 1) := by
    intro s hs
    have hre := re_one_add_smul_nonneg (z := z) hs.1
    have h1 : 0 ≤ (1 + (s : ℂ) * z).re := by nlinarith [hs.2]
    have hc : ContinuousAt (fun w : ℂ => w ^ (α : ℂ)) (1 + (s : ℂ) * z) :=
      Complex.continuousAt_cpow_const_of_re_pos (Or.inl h1) (by simpa using hα0)
    have hl : ContinuousAt (fun s : ℝ => 1 + (s : ℂ) * z) s := by fun_prop
    exact ((hc.comp (f := fun s : ℝ => 1 + (s : ℂ) * z) hl).sub
      continuousAt_const).continuousWithinAt
  have hderiv : ∀ s ∈ Ico (0 : ℝ) 1, HasDerivWithinAt f (f' s) (Ici s) s := by
    intro s hs
    have hre := re_one_add_smul_nonneg (z := z) hs.1
    rw [← hr] at hre
    have hpos : 0 < (1 + (s : ℂ) * z).re := by nlinarith [hs.2, hs.1]
    have hslit : 1 + (s : ℂ) * z ∈ Complex.slitPlane := Or.inl hpos
    have hl : HasDerivAt (fun w : ℂ => 1 + w * z) z (s : ℂ) := by
      simpa using ((hasDerivAt_id (s : ℂ)).mul_const z).const_add 1
    have h2 := (hl.cpow_const (c := (α : ℂ)) hslit).sub_const 1
    have h3 := h2.comp_ofReal
    refine HasDerivAt.hasDerivWithinAt ?_
    simpa [f, f', mul_comm, mul_assoc, mul_left_comm] using h3
  have hB : ∀ s ∈ Ico (0 : ℝ) 1, HasDerivWithinAt B (B' s) (Ici s) s := by
    intro s hs
    have hpos : 0 < 1 - s * r := by nlinarith [hs.2, hs.1]
    have hl : HasDerivAt (fun s : ℝ => 1 - s * r) (-r) s := by
      simpa using ((hasDerivAt_id s).mul_const r).const_sub 1
    have h2 := (hl.rpow_const (p := α) (Or.inl hpos.ne')).const_sub 1
    refine HasDerivAt.hasDerivWithinAt ?_
    convert h2 using 1
    simp only [B']; ring
  have hBc : ContinuousOn B (Icc 0 1) := by
    have : Continuous B := by
      simp only [B]
      exact continuous_const.sub
        ((Real.continuous_rpow_const hα0.le).comp (by fun_prop))
    exact this.continuousOn
  have hbound : ∀ s ∈ Ico (0 : ℝ) 1, ‖f' s‖ ≤ B' s := by
    intro s hs
    have hre := re_one_add_smul_nonneg (z := z) hs.1
    rw [← hr] at hre
    have hpos : 0 < 1 - s * r := by nlinarith [hs.2, hs.1]
    have hn : 1 - s * r ≤ ‖1 + (s : ℂ) * z‖ := hre.trans (Complex.re_le_norm _)
    have he : ((α : ℂ) - 1) = ((α - 1 : ℝ) : ℂ) := by push_cast; ring
    simp only [f', B', norm_mul, he, Complex.norm_cpow_real, Complex.norm_real,
      Real.norm_eq_abs, abs_of_pos hα0]
    rw [← hr]
    gcongr ?_ * _
    gcongr α * ?_
    exact Real.rpow_le_rpow_of_nonpos hpos hn (by linarith)
  have key := image_norm_le_of_norm_deriv_right_le_deriv_boundary' hcont hderiv
    (by simp [f, B]) hBc hB hbound (right_mem_Icc.2 zero_le_one)
  simp only [f, B, Complex.ofReal_one, one_mul] at key
  refine key.trans ?_
  have : 1 - r ≤ (1 - r) ^ α :=
    Real.self_le_rpow_of_le_one (by linarith) (by linarith) hα1
  linarith

lemma inv_q_cast (q : ℝ) : (1 / (q : ℂ)) = ((1 / q : ℝ) : ℂ) := by push_cast; rfl

lemma inv_q_pos {q : ℝ} (hq : 1 ≤ q) : 0 < 1 / q := by
  have : 0 < q := by linarith
  positivity

lemma inv_q_le {q : ℝ} (hq : 1 ≤ q) : 1 / q ≤ 1 := by
  rw [div_le_one (by linarith)]; exact hq

end X03

theorem X03.proof_2 : X03.stmt_2 := by
  intro q hq z hz
  rw [X03.inv_q_cast]
  exact X03.binom_bound (X03.inv_q_pos hq) (X03.inv_q_le hq) hz

namespace X03

/-! ### Real-side integrals (statement 1 preliminaries). -/

lemma rpow_iInt {ε : ℝ} (h0 : 0 < ε) :
    IntervalIntegrable (fun t : ℝ => t ^ (ε - 1)) volume 0 1 :=
  intervalIntegral.intervalIntegrable_rpow' (by linarith)

lemma rpow_iInt' {ε : ℝ} (h0 : 0 < ε) :
    IntervalIntegrable (fun t : ℝ => (1 - t) ^ (ε - 1)) volume 0 1 := by
  have := (rpow_iInt h0).comp_sub_left 1
  simpa using this.symm

lemma rpow_int {ε : ℝ} (h0 : 0 < ε) : ∫ t in (0:ℝ)..1, t ^ (ε - 1) = 1 / ε := by
  rw [integral_rpow (Or.inl (by linarith))]
  simp [Real.zero_rpow h0.ne']

lemma rpow_int' {ε : ℝ} (h0 : 0 < ε) : ∫ t in (0:ℝ)..1, (1 - t) ^ (ε - 1) = 1 / ε := by
  rw [intervalIntegral.integral_comp_sub_left (fun t : ℝ => t ^ (ε - 1)) 1]
  simpa using rpow_int h0

lemma half_rpow (ε : ℝ) : (1 / 2 : ℝ) ^ (ε - 1) = (2 : ℝ) ^ (1 - ε) := by
  rw [one_div, Real.inv_rpow (by norm_num), ← Real.rpow_neg (by norm_num)]
  ring_nf

/-! ### Statement 3: the bound on the complex integrand. -/

lemma cpow_pos {x : ℝ} (hx : 0 < x) (w : ℂ) :
    ((x : ℝ) : ℂ) ^ w = Complex.exp ((Real.log x : ℂ) * w) := by
  rw [Complex.cpow_def_of_ne_zero (by exact_mod_cast hx.ne'), Complex.ofReal_log hx.le]

lemma gC_alt (q : ℝ) (ε : ℂ) (t : ℝ) :
    gC q ε t = bC ε t - aC ε t * ((1 + zC q ε t) ^ (1 / (q : ℂ)) - 1) := by
  unfold gC; ring

lemma norm_zC {q : ℝ} {ε : ℂ} {t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    ‖zC q ε t‖ = (t / (1 - t)) ^ (q * (1 - ε.re)) := by
  unfold zC
  rw [Complex.norm_cpow_eq_rpow_re_of_pos (div_pos ht0 (by linarith))]
  congr 1
  simp

lemma gC_bound {q : ℝ} (hq : 1 ≤ q) {ε : ℂ} (hε : ε.re < 1) {t : ℝ} (ht0 : 0 < t)
    (ht : t ≤ 1 / 2) :
    ‖gC q ε t‖ ≤ 2 * ‖bC ε t‖ ∧ ‖bC ε t‖ ≤ (2 : ℝ) ^ (1 - ε.re) := by
  have ht1 : t < 1 := by linarith
  have h1 : 0 < 1 - t := by linarith
  have hu0 : 0 < t / (1 - t) := div_pos ht0 h1
  have hu1 : t / (1 - t) ≤ 1 := by rw [div_le_one h1]; linarith
  have hna : ‖aC ε t‖ = t ^ (ε.re - 1) := by
    unfold aC; rw [Complex.norm_cpow_eq_rpow_re_of_pos ht0]; simp
  have hnb : ‖bC ε t‖ = (1 - t) ^ (ε.re - 1) := by
    unfold bC; rw [Complex.norm_cpow_eq_rpow_re_of_pos h1]; simp
  have hnz := norm_zC (q := q) (ε := ε) ht0 ht1
  have hz1 : ‖zC q ε t‖ ≤ 1 := by
    rw [hnz]
    exact Real.rpow_le_one hu0.le hu1 (by nlinarith)
  have hb2 : ‖bC ε t‖ ≤ (2 : ℝ) ^ (1 - ε.re) := by
    rw [hnb, ← half_rpow]
    exact Real.rpow_le_rpow_of_nonpos (by norm_num) (by linarith) (by linarith)
  have haz : ‖aC ε t‖ * ‖zC q ε t‖ ≤ ‖bC ε t‖ := by
    rw [hna, hnb, hnz, Real.rpow_def_of_pos ht0, Real.rpow_def_of_pos hu0,
      Real.rpow_def_of_pos h1, ← Real.exp_add, Real.exp_le_exp,
      Real.log_div ht0.ne' h1.ne']
    have hl : Real.log t ≤ Real.log (1 - t) := Real.log_le_log ht0 (by linarith)
    have key : (Real.log t - Real.log (1 - t)) * ((1 - ε.re) * (q - 1)) ≤ 0 :=
      mul_nonpos_of_nonpos_of_nonneg (by linarith)
        (mul_nonneg (by linarith) (by linarith))
    nlinarith [key]
  have hB := binom_bound (inv_q_pos hq) (inv_q_le hq) hz1
  rw [← inv_q_cast] at hB
  refine ⟨?_, hb2⟩
  rw [gC_alt]
  calc ‖bC ε t - aC ε t * ((1 + zC q ε t) ^ (1 / (q : ℂ)) - 1)‖
      ≤ ‖bC ε t‖ + ‖aC ε t * ((1 + zC q ε t) ^ (1 / (q : ℂ)) - 1)‖ := norm_sub_le _ _
    _ ≤ ‖bC ε t‖ + ‖aC ε t‖ * ‖zC q ε t‖ := by
        rw [norm_mul]; gcongr
    _ ≤ 2 * ‖bC ε t‖ := by linarith

end X03

theorem X03.proof_3 : X03.stmt_3 := by
  intro q hq ε hε t ht0 ht
  obtain ⟨h1, h2⟩ := X03.gC_bound hq hε ht0 ht
  have : (2 : ℝ) ^ (2 - ε.re) = 2 * (2 : ℝ) ^ (1 - ε.re) := by
    rw [show (2 : ℝ) - ε.re = 1 + (1 - ε.re) by ring, Real.rpow_add (by norm_num), Real.rpow_one]
  rw [this]; linarith

namespace X03

/-! ### Statement 4: real `ε`. -/

lemma real_core {q : ℝ} (hq : 1 ≤ q) (ε : ℝ) {t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    (t ^ (q * (ε - 1)) + (1 - t) ^ (q * (ε - 1))) ^ (1 / q) =
      t ^ (ε - 1) * (1 + (t / (1 - t)) ^ (q * (1 - ε))) ^ (1 / q) := by
  have h1 : 0 < 1 - t := by linarith
  have hu0 : 0 < t / (1 - t) := div_pos ht0 h1
  have hq0 : q ≠ 0 := by linarith
  have hprod : t ^ (q * (ε - 1)) * (t / (1 - t)) ^ (q * (1 - ε)) = (1 - t) ^ (q * (ε - 1)) := by
    rw [Real.rpow_def_of_pos ht0, Real.rpow_def_of_pos hu0, Real.rpow_def_of_pos h1,
      ← Real.exp_add, Real.log_div ht0.ne' h1.ne']
    congr 1; ring
  have hsplit : t ^ (q * (ε - 1)) + (1 - t) ^ (q * (ε - 1)) =
      t ^ (q * (ε - 1)) * (1 + (t / (1 - t)) ^ (q * (1 - ε))) := by
    rw [mul_add, mul_one, hprod]
  rw [hsplit, Real.mul_rpow (Real.rpow_nonneg ht0.le _)
    (by positivity), ← Real.rpow_mul ht0.le]
  congr 2
  field_simp

lemma gC_ofReal {q : ℝ} (hq : 1 ≤ q) (ε : ℝ) {t : ℝ} (ht0 : 0 < t) (ht : t < 1) :
    gC q (ε : ℂ) t = ((gR q ε t : ℝ) : ℂ) := by
  have h1 : 0 ≤ 1 - t := by linarith
  have hu : 0 ≤ t / (1 - t) := div_nonneg ht0.le h1
  have ha : aC (ε : ℂ) t = ((t ^ (ε - 1) : ℝ) : ℂ) := by
    unfold aC; rw [Complex.ofReal_cpow ht0.le]; push_cast; rfl
  have hb : bC (ε : ℂ) t = (((1 - t) ^ (ε - 1) : ℝ) : ℂ) := by
    unfold bC; rw [Complex.ofReal_cpow h1]; push_cast; rfl
  have hz : zC q (ε : ℂ) t = (((t / (1 - t)) ^ (q * (1 - ε)) : ℝ) : ℂ) := by
    unfold zC; rw [Complex.ofReal_cpow hu]; push_cast; rfl
  have hz0 : 0 ≤ (t / (1 - t)) ^ (q * (1 - ε)) := Real.rpow_nonneg hu _
  have hs : (1 + (((t / (1 - t)) ^ (q * (1 - ε)) : ℝ) : ℂ)) ^ (1 / (q : ℂ)) =
      (((1 + (t / (1 - t)) ^ (q * (1 - ε))) ^ (1 / q) : ℝ) : ℂ) := by
    rw [inv_q_cast, ← Complex.ofReal_one, ← Complex.ofReal_add,
      ← Complex.ofReal_cpow (by linarith)]
  unfold gC gR
  rw [ha, hb, hz, hs, real_core hq ε ht0 ht]
  push_cast; ring

lemma gR_meas (q ε : ℝ) : Measurable (gR q ε) := by
  unfold gR; fun_prop

lemma gR_symm (q ε t : ℝ) : gR q ε (1 - t) = gR q ε t := by
  unfold gR
  simp only [sub_sub_cancel]
  ring_nf

lemma gR_bound {q : ℝ} (hq : 1 ≤ q) {ε t : ℝ} (hε : ε < 1) (ht0 : 0 < t) (ht1 : t < 1) :
    ‖gR q ε t‖ ≤ (2 : ℝ) ^ (2 - ε) := by
  have key : ∀ s : ℝ, 0 < s → s ≤ 1 / 2 → ‖gR q ε s‖ ≤ (2 : ℝ) ^ (2 - ε) := by
    intro s hs0 hs
    have := proof_3 q hq (ε : ℂ) (by simpa using hε) s hs0 hs
    rw [gC_ofReal hq ε hs0 (by linarith), Complex.norm_real] at this
    simpa using this
  rcases le_total t (1 / 2) with h | h
  · exact key t ht0 h
  · rw [← gR_symm]; exact key (1 - t) (by linarith) (by linarith)

lemma gR_iInt {q : ℝ} (hq : 1 ≤ q) {ε : ℝ} (hε : ε < 1) :
    IntervalIntegrable (gR q ε) volume 0 1 := by
  refine IntervalIntegrable.mono_fun' (g := fun _ => (2 : ℝ) ^ (2 - ε))
    intervalIntegrable_const (gR_meas q ε).aestronglyMeasurable ?_
  rw [uIoc_of_le zero_le_one, ← restrict_Ioo_eq_restrict_Ioc]
  exact ae_restrict_of_forall_mem measurableSet_Ioo fun t ht => gR_bound hq hε ht.1 ht.2

lemma FR_half {q : ℝ} (hq : 1 ≤ q) {ε : ℝ} (hε : ε < 1) :
    FR q ε = 2 * ∫ t in (0:ℝ)..(1/2), gR q ε t := by
  have hi := gR_iInt hq hε
  have h1 : IntervalIntegrable (gR q ε) volume 0 (1/2) :=
    hi.mono_set (by rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
                    exact Icc_subset_Icc le_rfl (by norm_num))
  have h2 : IntervalIntegrable (gR q ε) volume (1/2) 1 :=
    hi.mono_set (by rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
                    exact Icc_subset_Icc (by norm_num) le_rfl)
  unfold FR
  rw [← intervalIntegral.integral_add_adjacent_intervals h1 h2]
  have : ∫ t in (1/2:ℝ)..1, gR q ε t = ∫ t in (0:ℝ)..(1/2), gR q ε t := by
    have := intervalIntegral.integral_comp_sub_left (gR q ε) (a := 0) (b := 1/2) 1
    simp only [gR_symm] at this
    rw [this]; norm_num
  rw [this]; ring

end X03

theorem X03.proof_4 : X03.stmt_4 := by
  intro q hq ε hε
  refine ⟨X03.gR_iInt hq hε, ?_⟩
  unfold X03.FC
  rw [X03.FR_half hq hε]
  have : ∫ t in (0:ℝ)..(1/2), X03.gC q (ε : ℂ) t =
      ∫ t in (0:ℝ)..(1/2), ((X03.gR q ε t : ℝ) : ℂ) := by
    refine intervalIntegral.integral_congr_ae (Filter.Eventually.of_forall fun t ht => ?_)
    rw [uIoc_of_le (by norm_num)] at ht
    exact X03.gC_ofReal hq ε ht.1 (by linarith [ht.2])
  rw [this, intervalIntegral.integral_ofReal]
  push_cast; ring

namespace X03

/-! ### Statement 1. -/

lemma gR_iInt_le {q : ℝ} (hq : 1 ≤ q) {ε : ℝ} (hε : ε ≤ 1) :
    IntervalIntegrable (gR q ε) volume 0 1 := by
  rcases lt_or_eq_of_le hε with h | h
  · exact gR_iInt hq h
  · subst h
    have : gR q 1 = fun _ => (2 : ℝ) - 2 ^ (1 / q) := by
      funext t; unfold gR; simp; norm_num
    rw [this]; exact intervalIntegrable_const

end X03

theorem X03.proof_1 : X03.stmt_1 := by
  intro q hq p hp
  set ε := 1 / p with hεdef
  have hp0 : 0 < p := by linarith
  have h0 : 0 < ε := by positivity
  have h1 : ε ≤ 1 := by rw [hεdef, div_le_one hp0]; exact hp
  have ha := X03.rpow_iInt h0
  have hb := X03.rpow_iInt' h0
  have hab : IntervalIntegrable (fun t => t ^ (ε - 1) + (1 - t) ^ (ε - 1)) volume 0 1 :=
    ha.add hb
  have hg := X03.gR_iInt_le hq h1
  have heq : X03.lqIntegrand p q = fun t => (t ^ (ε - 1) + (1 - t) ^ (ε - 1)) - X03.gR q ε t := by
    funext t
    have : q / p - q = q * (ε - 1) := by rw [hεdef]; ring
    simp only [X03.lqIntegrand, X03.gR, this]
    ring
  have hP : IntervalIntegrable (X03.lqIntegrand p q) volume 0 1 := by
    rw [heq]; exact hab.sub hg
  refine ⟨hP, hg, ?_⟩
  have hL : X03.Lq p q = (4 / p) * ∫ t in (0:ℝ)..1, X03.lqIntegrand p q t := rfl
  rw [hL, heq, intervalIntegral.integral_sub hab hg, intervalIntegral.integral_add ha hb,
    X03.rpow_int h0, X03.rpow_int' h0]
  unfold X03.FR
  have : (4 / p) = 4 * ε := by rw [hεdef]; ring
  rw [this]
  field_simp
  ring

namespace X03

/-! ### Statement 5: holomorphy. -/

lemma gC_meas (q : ℝ) (ε : ℂ) : Measurable (gC q ε) := by
  unfold gC aC bC zC
  fun_prop

/-- explicit `ε`-derivative of `gC`. -/
noncomputable def gCd (q : ℝ) (e : ℂ) (t : ℝ) : ℂ :=
  aC e t * Complex.log (t : ℂ) + bC e t * Complex.log ((1 - t : ℝ) : ℂ) -
    (aC e t * Complex.log (t : ℂ) * (1 + zC q e t) ^ (1 / (q : ℂ)) +
      aC e t * ((1 / (q : ℂ)) * (1 + zC q e t) ^ ((1 / (q : ℂ)) - 1) *
        (zC q e t * Complex.log ((t / (1 - t) : ℝ) : ℂ) * ((q : ℂ) * (0 - 1)))))

lemma gCd_meas (q : ℝ) (ε : ℂ) : Measurable (gCd q ε) := by
  unfold gCd aC bC zC
  have := Complex.measurable_log
  fun_prop

lemma zC_slit {q : ℝ} (hq : 1 ≤ q) {e : ℂ} (he : e.re < 1) {t : ℝ} (ht0 : 0 < t)
    (ht : t ≤ 1 / 2) : 1 + zC q e t ∈ Complex.slitPlane := by
  have h1 : 0 < 1 - t := by linarith
  rcases eq_or_lt_of_le ht with h | h
  · subst h
    unfold zC
    norm_num
  · have hu0 : 0 < t / (1 - t) := div_pos ht0 h1
    have hu1 : t / (1 - t) < 1 := by rw [div_lt_one h1]; linarith
    have hn : ‖zC q e t‖ < 1 := by
      rw [norm_zC ht0 (by linarith)]
      exact Real.rpow_lt_one hu0.le hu1 (by nlinarith)
    exact Complex.mem_slitPlane_of_norm_lt_one hn

lemma gC_hasDerivAt {q : ℝ} (hq : 1 ≤ q) {e : ℂ} (he : e.re < 1) {t : ℝ} (ht0 : 0 < t)
    (ht : t ≤ 1 / 2) : HasDerivAt (fun x => gC q x t) (gCd q e t) e := by
  have h1 : 0 < 1 - t := by linarith
  have ha : HasDerivAt (fun x => aC x t) (aC e t * Complex.log (t : ℂ) * 1) e :=
    ((hasDerivAt_id e).sub_const 1).const_cpow (Or.inl (by exact_mod_cast ht0.ne'))
  have hb : HasDerivAt (fun x => bC x t) (bC e t * Complex.log ((1 - t : ℝ) : ℂ) * 1) e :=
    ((hasDerivAt_id e).sub_const 1).const_cpow (Or.inl (by exact_mod_cast h1.ne'))
  have hz : HasDerivAt (fun x => zC q x t)
      (zC q e t * Complex.log ((t / (1 - t) : ℝ) : ℂ) * ((q : ℂ) * (0 - 1))) e :=
    (((hasDerivAt_const e (1 : ℂ)).sub (hasDerivAt_id e)).const_mul (q : ℂ)).const_cpow
      (Or.inl (by exact_mod_cast (div_pos ht0 h1).ne'))
  have hS := ((hasDerivAt_const e (1 : ℂ)).add hz).cpow_const (c := (1 / (q : ℂ)))
    (zC_slit hq he ht0 ht)
  have := (ha.add hb).sub (ha.mul hS)
  unfold gC gCd
  convert this using 1
  simp only [Pi.add_apply]
  ring

lemma gC_iInt {q : ℝ} (hq : 1 ≤ q) {ε : ℂ} (hε : ε.re < 1) :
    IntervalIntegrable (gC q ε) volume 0 (1 / 2) := by
  refine IntervalIntegrable.mono_fun' (g := fun _ => (2 : ℝ) ^ (2 - ε.re))
    intervalIntegrable_const (gC_meas q ε).aestronglyMeasurable ?_
  rw [uIoc_of_le (by norm_num)]
  exact ae_restrict_of_forall_mem measurableSet_Ioc fun t ht => proof_3 q hq ε hε t ht.1 ht.2

lemma gCd_bound {q : ℝ} (hq : 1 ≤ q) {x₀ : ℂ} (h : x₀.re < 1) {x : ℂ}
    (hx : x ∈ Metric.ball x₀ ((1 - x₀.re) / 3)) {t : ℝ} (ht0 : 0 < t) (ht : t ≤ 1 / 2) :
    ‖gCd q x t‖ ≤ (2 : ℝ) ^ (2 - x₀.re + 2 * ((1 - x₀.re) / 3)) / ((1 - x₀.re) / 3) := by
  set δ := (1 - x₀.re) / 3 with hδ
  have hδ0 : 0 < δ := by rw [hδ]; linarith
  have hxre : |x.re - x₀.re| < δ := by
    have := Complex.abs_re_le_norm (x - x₀)
    rw [Metric.mem_ball, Complex.dist_eq] at hx
    simp only [Complex.sub_re] at this
    linarith
  have hsub : Metric.closedBall x δ ⊆ {e : ℂ | e.re < 1} := by
    intro y hy
    rw [Metric.mem_closedBall, Complex.dist_eq] at hy
    have := Complex.abs_re_le_norm (y - x)
    simp only [Complex.sub_re] at this
    simp only [Set.mem_ofPred_eq]
    linarith [(abs_le.1 (this.trans hy)).2, (abs_lt.1 hxre).2]
  have hdiff : DifferentiableOn ℂ (fun e => gC q e t) {e : ℂ | e.re < 1} := fun e he =>
    (gC_hasDerivAt hq he ht0 ht).differentiableAt.differentiableWithinAt
  have hxU : x.re < 1 := hsub (Metric.mem_closedBall_self hδ0.le)
  rw [← (gC_hasDerivAt hq hxU ht0 ht).deriv]
  apply Complex.norm_deriv_le_of_forall_mem_sphere_norm_le hδ0 (hdiff.diffContOnCl_ball hsub)
  intro y hy
  have hyU : y.re < 1 := hsub (Metric.sphere_subset_closedBall hy)
  rw [mem_sphere_iff_norm] at hy
  have hyre : x₀.re - 2 * δ ≤ y.re := by
    have := Complex.abs_re_le_norm (y - x)
    simp only [Complex.sub_re] at this
    linarith [(abs_le.1 (this.trans hy.le)).1, (abs_lt.1 hxre).1]
  refine (proof_3 q hq y hyU t ht0 ht).trans ?_
  exact Real.rpow_le_rpow_of_exponent_le (by norm_num) (by linarith)

lemma FC_hasDerivAt {q : ℝ} (hq : 1 ≤ q) {x₀ : ℂ} (h : x₀.re < 1) :
    HasDerivAt (FC q) (2 * ∫ t in (0:ℝ)..(1 / 2), gCd q x₀ t) x₀ := by
  have hδ0 : 0 < (1 - x₀.re) / 3 := by linarith
  have key := intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (μ := volume) (a := 0) (b := 1 / 2) (F := fun e t => gC q e t) (F' := gCd q)
    (x₀ := x₀) (s := Metric.ball x₀ ((1 - x₀.re) / 3))
    (bound := fun _ => (2 : ℝ) ^ (2 - x₀.re + 2 * ((1 - x₀.re) / 3)) / ((1 - x₀.re) / 3))
    (Metric.ball_mem_nhds x₀ hδ0)
    (Filter.Eventually.of_forall fun x => (gC_meas q x).aestronglyMeasurable)
    (gC_iInt hq h) ((gCd_meas q x₀).aestronglyMeasurable)
    (Filter.Eventually.of_forall fun t ht x hx => by
      rw [uIoc_of_le (by norm_num)] at ht
      exact gCd_bound hq h hx ht.1 ht.2)
    intervalIntegrable_const
    (Filter.Eventually.of_forall fun t ht x hx => by
      rw [uIoc_of_le (by norm_num)] at ht
      have hxre : x.re < 1 := by
        have := Complex.abs_re_le_norm (x - x₀)
        rw [Metric.mem_ball, Complex.dist_eq] at hx
        simp only [Complex.sub_re] at this
        linarith [(abs_le.1 (this.trans hx.le)).2]
      exact gC_hasDerivAt hq hxre ht.1 ht.2)
  exact key.2.const_mul (2 : ℂ)

end X03

theorem X03.proof_5 : X03.stmt_5 := fun _ hq =>
  ⟨fun _ hε => X03.gC_iInt hq hε, fun _ hx₀ =>
    (X03.FC_hasDerivAt hq hx₀).differentiableAt.differentiableWithinAt⟩

namespace X03

/-! ### Statements 6 and 7: the Taylor series of `FC q` at `0`. -/

lemma FC_diff_ball {q : ℝ} (hq : 1 ≤ q) : DifferentiableOn ℂ (FC q) (Metric.ball 0 1) := by
  refine (proof_5 q hq).2.mono fun e he => ?_
  rw [mem_ball_zero_iff] at he
  have := Complex.re_le_norm e
  simp only [Set.mem_ofPred_eq]; linarith

lemma Lq_FC {q : ℝ} (hq : 1 ≤ q) {p : ℝ} (hp : 1 < p) :
    ((Lq p q / 8 - 1 : ℝ) : ℂ) = -(1 / (2 * (p : ℂ))) * FC q ((1 / p : ℝ) : ℂ) := by
  have hp0 : (0 : ℝ) < p := by linarith
  have h1 : (1 / p : ℝ) < 1 := by rw [div_lt_one hp0]; exact hp
  rw [(proof_4 q hq (1 / p) h1).2, (proof_1 q hq p hp.le).2.2]
  push_cast
  field_simp
  ring

lemma Kk_term (q : ℝ) {p : ℝ} (hp0 : 0 < p) (k : ℕ) :
    Kk q k / (p : ℂ) ^ (k + 1) =
      -(1 / (2 * (p : ℂ))) * ((k.factorial : ℂ)⁻¹ • ((1 / p : ℝ) : ℂ) ^ k •
        iteratedDeriv k (FC q) 0) := by
  have hpc : (p : ℂ) ≠ 0 := by exact_mod_cast hp0.ne'
  simp only [smul_eq_mul, Kk]
  push_cast
  rw [div_pow, one_pow]
  have : (p : ℂ) ^ k ≠ 0 := pow_ne_zero _ hpc
  have : (k.factorial : ℂ) ≠ 0 := by exact_mod_cast k.factorial_ne_zero
  field_simp
  ring

lemma taylor_tail {f : ℂ → ℂ} (hf : DifferentiableOn ℂ f (Metric.ball 0 1)) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ε : ℂ, ‖ε‖ ≤ 3 / 8 → ∀ K : ℕ,
      ‖f ε - ∑ n ∈ Finset.range K, (n.factorial : ℂ)⁻¹ • ε ^ n • iteratedDeriv n f 0‖ ≤
        2 * C * (4 / 3) ^ K * ‖ε‖ ^ K := by
  set T : ℕ → ℂ → ℂ := fun n z => (n.factorial : ℂ)⁻¹ • z ^ n • iteratedDeriv n f 0 with hT
  have hsum : ∀ z : ℂ, ‖z‖ < 1 → HasSum (fun n => T n z) (f z) := by
    intro z hz
    have := Complex.hasSum_taylorSeries_on_ball hf (z := z) (by simpa using hz)
    simpa only [sub_zero] using this
  have H34 := hsum (3 / 4) (by norm_num)
  obtain ⟨C, hC⟩ := (Metric.isBounded_iff_subset_closedBall (0 : ℂ)).1
    (Metric.isBounded_range_of_tendsto _ H34.summable.tendsto_atTop_zero)
  have hCn : ∀ n, ‖T n (3 / 4)‖ ≤ C := fun n => by
    have := hC ⟨n, rfl⟩
    simpa using this
  have hC0 : 0 ≤ C := (norm_nonneg _).trans (hCn 0)
  refine ⟨C, hC0, fun ε hε K => ?_⟩
  set ρ := ‖ε‖ * (4 / 3) with hρ
  have hρ0 : 0 ≤ ρ := by positivity
  have hρ1 : ρ ≤ 1 / 2 := by rw [hρ]; linarith
  have hTn : ∀ n, ‖T n ε‖ ≤ C * ρ ^ n := by
    intro n
    have he : ‖T n ε‖ = ‖T n (3 / 4)‖ * ρ ^ n := by
      simp only [hT, smul_eq_mul, norm_mul, norm_pow, hρ, mul_pow]
      have : (3 / 4 : ℝ) ^ n * (4 / 3) ^ n = 1 := by rw [← mul_pow]; norm_num
      have h34 : ‖(3 / 4 : ℂ)‖ = 3 / 4 := by norm_num
      rw [h34]
      linear_combination (-(‖(n.factorial : ℂ)⁻¹‖ * ‖ε‖ ^ n * ‖iteratedDeriv n f 0‖)) * this
    rw [he]
    exact mul_le_mul_of_nonneg_right (hCn n) (pow_nonneg hρ0 _)
  have hεs := hsum ε (by linarith)
  have htail := (hasSum_nat_add_iff' K).2 hεs
  have hg : HasSum (fun n => C * ρ ^ (n + K)) (C * ρ ^ K * (1 - ρ)⁻¹) := by
    have := (hasSum_geometric_of_lt_one hρ0 (by linarith)).mul_left (C * ρ ^ K)
    convert this using 2 with n
    ring
  have hb := htail.norm_le_of_bounded hg (fun n => hTn (n + K))
  refine hb.trans ?_
  have h2 : (1 - ρ)⁻¹ ≤ 2 := by
    rw [inv_le_comm₀ (by linarith) (by norm_num)]; linarith
  have hρK : ρ ^ K = ‖ε‖ ^ K * (4 / 3) ^ K := by rw [hρ, mul_pow]
  calc C * ρ ^ K * (1 - ρ)⁻¹ ≤ C * ρ ^ K * 2 := by
        gcongr
    _ = 2 * C * (4 / 3) ^ K * ‖ε‖ ^ K := by rw [hρK]; ring

end X03

theorem X03.proof_6 : X03.stmt_6 := by
  intro q hq K
  obtain ⟨C, hC0, hC⟩ := X03.taylor_tail (X03.FC_diff_ball hq)
  refine Asymptotics.IsBigO.of_bound (C * (4 / 3) ^ K) ?_
  filter_upwards [Filter.eventually_ge_atTop (3 : ℝ)] with p hp
  have hp0 : (0 : ℝ) < p := by linarith
  set ε : ℂ := ((1 / p : ℝ) : ℂ) with hε
  have hεn : ‖ε‖ = 1 / p := by
    rw [hε, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (by positivity)]
  have hε38 : ‖ε‖ ≤ 3 / 8 := by
    rw [hεn, div_le_iff₀ hp0]; linarith
  have heq : ((X03.Lq p q / 8 - 1 : ℝ) : ℂ) -
      ∑ k ∈ Finset.range K, X03.Kk q k / (p : ℂ) ^ (k + 1) =
      -(1 / (2 * (p : ℂ))) * (X03.FC q ε - ∑ n ∈ Finset.range K,
        (n.factorial : ℂ)⁻¹ • ε ^ n • iteratedDeriv n (X03.FC q) 0) := by
    rw [X03.Lq_FC hq (by linarith), mul_sub, Finset.mul_sum]
    congr 1
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [X03.Kk_term q hp0 k]
  rw [heq, norm_mul]
  have h2p : ‖-(1 / (2 * (p : ℂ)))‖ = (1 / p) / 2 := by
    rw [norm_neg, norm_div, norm_one, norm_mul, Complex.norm_real, Real.norm_eq_abs,
      abs_of_pos hp0]
    norm_num; ring
  have hg : ‖(1 / p) ^ (K + 1)‖ = (1 / p) ^ (K + 1) := by
    rw [Real.norm_eq_abs, abs_of_pos (by positivity)]
  rw [h2p, hg]
  have := hC ε hε38 K
  rw [hεn] at this
  calc 1 / p / 2 * ‖X03.FC q ε - ∑ n ∈ Finset.range K,
        (n.factorial : ℂ)⁻¹ • ε ^ n • iteratedDeriv n (X03.FC q) 0‖
      ≤ 1 / p / 2 * (2 * C * (4 / 3) ^ K * (1 / p) ^ K) := by gcongr
    _ = C * (4 / 3) ^ K * (1 / p) ^ (K + 1) := by ring

theorem X03.proof_7 : X03.stmt_7 := by
  intro q hq p hp
  have hp0 : (0 : ℝ) < p := by linarith
  have h1 : (1 / p : ℝ) < 1 := by rw [div_lt_one hp0]; exact hp
  have hz : ((1 / p : ℝ) : ℂ) ∈ Metric.ball (0 : ℂ) 1 := by
    rw [mem_ball_zero_iff, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (by positivity)]
    exact h1
  have H := Complex.hasSum_taylorSeries_on_ball (X03.FC_diff_ball hq) hz
  have H2 := H.mul_left (-(1 / (2 * (p : ℂ))))
  rw [X03.Lq_FC hq hp]
  convert H2 using 1
  funext k
  rw [X03.Kk_term q hp0 k, sub_zero]
