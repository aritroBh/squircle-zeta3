import Mathlib
import Zeta3.Wave2.D00_Dilog
/-!
# D04 -- the two closed forms of `M_1(α)` (paper labels `thm:Malpha`, `tab:corner`,
`app:corner`, `prop:M1short`)

Source: `paper/ieee/main.tex`.

**Setup, lines 1466-1474** (Subsection "Smoothed Corners of Any Angle"):

> its length is $2+2u\,M(u;\alpha)$, where
> \begin{multline*}
> M(u;\alpha)=\int_0^1\Bigl[(1+g)^{-1-u}\sqrt{1+g^{2-2u}-2\cos\alpha\,g^{1-u}}\\
> -1\Bigr]g^{u-1}\dd g=\sum_{k\ge0}M_k(\alpha)u^k,
> \end{multline*}
> ... with $c=\cos\frac\alpha2$,
> \[ M_0(\alpha)=-(1+c)\log(1+c)-(1-c)\log(1-c). \]

**`thm:Malpha`, lines 1476-1478**:

> For $0<\alpha<\pi$, with $c=\cos\frac\alpha2$ and $s=\sin\frac\alpha2$,
> $M_1(\alpha)=\zeta(2)+U+V+\mathcal L$, where $U$ (19 dilogarithms), $V$ (8 dilogarithms that
> depend on $c$ only) and $\mathcal L$ (five products of logarithms) are given in
> Table~\ref{tab:corner} (Appendix~\ref{app:corner}). Every $\Li_2$ there is the real
> dilogarithm at a point of $(-\infty,1)$, every logarithm is of a positive number, and every
> coefficient lies in $\mathbb Z[\frac12]+\mathbb Z[\frac12]\,c$.

**Table `tab:corner`, lines 2307-2318** (verbatim):

> With $X_-=\frac{1+s-c}{1-c}$ and $X_+=\frac{1+s+c}{1+c}$:
> U&=-\Li_2\bigl(-\tfrac1s\bigr)+\Li_2\bigl(-\tfrac{s}{1-c}\bigr)+\Li_2\bigl(-\tfrac{s}{1+c}\bigr)
> -\Li_2\bigl(-\tfrac{1+s}{c}\bigr)-\Li_2\bigl(-\tfrac{s(1+s)}{c^2}\bigr)+\Li_2\bigl(-\tfrac{s(1+s)}{c(1+c)}\bigr)
> +c\Li_2\bigl(-\tfrac{1+s-c}{c}\bigr)\\
> &\quad+c\Li_2\bigl(-\tfrac{1+s-c}{c(1-c)}\bigr)-2c\Li_2\bigl(-\tfrac{1+s-c}{2c}\bigr)
> -2\Li_2\bigl(\tfrac1{1+s}\bigr)+\Li_2\bigl(\tfrac{c^2}{1+s}\bigr)+\Li_2\bigl(\tfrac{c}{1+s}\bigr)+\Li_2\bigl(\tfrac{s}{1+s}\bigr)-\Li_2\bigl(\tfrac{c(1-c)}{s(1+s)}\bigr)\\
> &\quad-c\Li_2\bigl(\tfrac{1-c}{1+s-c}\bigr)+c\Li_2\bigl(\tfrac{1+c}{1+s+c}\bigr)+c\Li_2\bigl(\tfrac{c}{1+s+c}\bigr)+c\Li_2\bigl(\tfrac{c(1+c)}{1+s+c}\bigr)-2c\Li_2\bigl(\tfrac{2c}{1+s+c}\bigr),\\
> V&=(1-c)\Li_2\bigl(-\tfrac1c\bigr)+\Li_2\bigl(1-\tfrac1{c^2}\bigr)-(1+c)\Li_2\bigl(-\tfrac{1-c}{c}\bigr)+2c\Li_2\bigl(-\tfrac{1-c}{2c}\bigr)-\Li_2(c^2)-(1+c)\Li_2(c)\\
> &\quad+(1-c)\Li_2\bigl(\tfrac{c}{1+c}\bigr)+2c\Li_2\bigl(\tfrac{2c}{1+c}\bigr),\\
> \mathcal L&=\log(1+s)\log\tfrac{c}{1+s}+\log\tfrac{c}{1+c}\log\tfrac{1+s}{s}+c\log4\,\log X_-+\tfrac c2\log^2X_--\tfrac c2\log X_+\,\log\bigl((1+c)(1+s+c)\bigr).

**`prop:M1short`, lines 1623-1629**:

> For $0<\alpha<\pi$, with $\lambda=\log2$ and $M_0(\alpha)$ as above,
> \begin{multline*}
> M_1(\alpha)=\zeta(2)-\lambda^2-\lambda M_0(\alpha)-(1-c)\bigl[\Li_2(c)+\Li_2(\tfrac{1-c}2)\bigr]\\
> -(1+c)\bigl[\Li_2(-c)+\Li_2(\tfrac{1+c}2)\bigr].
> \end{multline*}

## Formalization

* `Li_2` is `D00.Li2` throughout; `ζ(2)` is written `π²/6` (Mathlib: `riemannZeta_two`).
* `U` and `V` are written as lists of (coefficient, argument) pairs, one entry per table term,
  in the table's order (19 and 8 entries), and summed as `Σ coef * Li2 arg`. This makes the side
  claim "every `Li_2` argument is `< 1`" refer to exactly the arguments used in the formula.
* `log² X` is `(log X)^2`; `log 4` is `Real.log 4`.
* For `0 < α < π` we have `0 < c, s < 1`, so no denominator in the table vanishes and no
  `x/0 = 0` convention is used; all statements carry the hypothesis `0 < α < π`.
* `stmt_1`: table form `=` four-dilogarithm form (the content shared by `thm:Malpha` and
  `prop:M1short`).
* `stmt_2`: every `Li_2` argument of the table is `< 1`.
* `stmt_3`: every logarithm argument of the table (`1+s`, `c/(1+s)`, `c/(1+c)`, `(1+s)/s`, `4`,
  `X_-`, `X_+`, `(1+c)(1+s+c)`) is `> 0`.
* `stmt_4`: `thm:Malpha` itself, with `M_1(α)` read as the `u`-derivative at `u = 0` of the
  paper's integral `M(u;α)` (the coefficient of `u` in `Σ M_k(α) u^k`). The integrand
  `[(1+g)^{-1-u}√(1+g^{2-2u}-2cos α g^{1-u}) - 1] g^{u-1}` is integrable on `[0,1]` for `|u| < 1`
  (bounded by a multiple of `g^{-|u|}`), so `HasDerivAt` is not satisfied via the junk value of a
  non-integrable integral. The radicand is `≥ 0` (it is `|1 - g^{1-u}e^{iα}|²`).
* `stmt_5`: `prop:M1short` itself, in the same reading of `M_1(α)`.

Numerical sanity check of the transcription (mpmath, 30 digits): table form, four-dilogarithm form
and the appendix integral for `M_1` agree at `α = 0.7` (`0.44055722308436422152939730707`),
`α = 1.9` (`0.170563655834963449429426571892`), `α = 2π/3` (`0.125795678937973233753868197684`)
and `α = 3.0` (`0.00250222317317432384821017382`, to ~27 digits).

NOT COVERED:
* The coefficient rule "every coefficient lies in `ℤ[½] + ℤ[½] c`" and "`V` depends on `c`
  only": both are true by inspection of the definitions below (the coefficients are written as
  `±1, ±2, ±c, ±2c, 1±c, c/2, -c/2` and `vTerms` takes only `c`), so a Lean statement would be
  vacuous bookkeeping; `log 4` and the `c/2` factors of `𝓛` are not `Li_2` coefficients.
* The full expansion `M(u;α) = Σ_k M_k(α) u^k` (analyticity in `u`) and `M_0`; only the
  derivative at `u = 0` is used.
* The special value at `α = 2π/3` stated after `thm:Malpha` (a separate claim).
-/
-- BEGIN STATEMENTS
namespace D04

open Real

/-- `c = cos(α/2)`. -/
noncomputable def c (α : ℝ) : ℝ := Real.cos (α / 2)

/-- `s = sin(α/2)`. -/
noncomputable def s (α : ℝ) : ℝ := Real.sin (α / 2)

/-- `X_- = (1+s-c)/(1-c)` (table `tab:corner`). -/
noncomputable def Xm (α : ℝ) : ℝ := (1 + s α - c α) / (1 - c α)

/-- `X_+ = (1+s+c)/(1+c)` (table `tab:corner`). -/
noncomputable def Xp (α : ℝ) : ℝ := (1 + s α + c α) / (1 + c α)

/-- The 19 terms of `U` in table `tab:corner`, as (coefficient, `Li_2` argument), in order. -/
noncomputable def uTerms (α : ℝ) : List (ℝ × ℝ) :=
  let c := c α
  let s := s α
  [ (-1,     -(1 / s)),
    (1,      -(s / (1 - c))),
    (1,      -(s / (1 + c))),
    (-1,     -((1 + s) / c)),
    (-1,     -(s * (1 + s) / c ^ 2)),
    (1,      -(s * (1 + s) / (c * (1 + c)))),
    (c,      -((1 + s - c) / c)),
    (c,      -((1 + s - c) / (c * (1 - c)))),
    (-2 * c, -((1 + s - c) / (2 * c))),
    (-2,     1 / (1 + s)),
    (1,      c ^ 2 / (1 + s)),
    (1,      c / (1 + s)),
    (1,      s / (1 + s)),
    (-1,     c * (1 - c) / (s * (1 + s))),
    (-c,     (1 - c) / (1 + s - c)),
    (c,      (1 + c) / (1 + s + c)),
    (c,      c / (1 + s + c)),
    (c,      c * (1 + c) / (1 + s + c)),
    (-2 * c, 2 * c / (1 + s + c)) ]

/-- The 8 terms of `V` in table `tab:corner` (they depend on `c` only), in order. -/
noncomputable def vTerms (c : ℝ) : List (ℝ × ℝ) :=
  [ (1 - c,    -(1 / c)),
    (1,        1 - 1 / c ^ 2),
    (-(1 + c), -((1 - c) / c)),
    (2 * c,    -((1 - c) / (2 * c))),
    (-1,       c ^ 2),
    (-(1 + c), c),
    (1 - c,    c / (1 + c)),
    (2 * c,    2 * c / (1 + c)) ]

/-- `Σ coef · Li₂(arg)` over a list of terms. -/
noncomputable def liSum (l : List (ℝ × ℝ)) : ℝ :=
  (l.map (fun p => p.1 * D00.Li2 p.2)).sum

/-- `U` of table `tab:corner` (19 dilogarithms). -/
noncomputable def U (α : ℝ) : ℝ := liSum (uTerms α)

/-- `V` of table `tab:corner` (8 dilogarithms in `c`). -/
noncomputable def V (α : ℝ) : ℝ := liSum (vTerms (c α))

/-- `𝓛` of table `tab:corner` (five products of logarithms). -/
noncomputable def Lcal (α : ℝ) : ℝ :=
  let c := c α
  let s := s α
  Real.log (1 + s) * Real.log (c / (1 + s))
    + Real.log (c / (1 + c)) * Real.log ((1 + s) / s)
    + c * Real.log 4 * Real.log (Xm α)
    + c / 2 * Real.log (Xm α) ^ 2
    - c / 2 * Real.log (Xp α) * Real.log ((1 + c) * (1 + s + c))

/-- The right-hand side of `thm:Malpha`: `ζ(2) + U + V + 𝓛`. -/
noncomputable def Mtable (α : ℝ) : ℝ := π ^ 2 / 6 + U α + V α + Lcal α

/-- `M_0(α) = -(1+c)log(1+c) - (1-c)log(1-c)` (line 1474). -/
noncomputable def M0 (α : ℝ) : ℝ :=
  -(1 + c α) * Real.log (1 + c α) - (1 - c α) * Real.log (1 - c α)

/-- The right-hand side of `prop:M1short` (`λ = log 2`). -/
noncomputable def Mshort (α : ℝ) : ℝ :=
  π ^ 2 / 6 - Real.log 2 ^ 2 - Real.log 2 * M0 α
    - (1 - c α) * (D00.Li2 (c α) + D00.Li2 ((1 - c α) / 2))
    - (1 + c α) * (D00.Li2 (-(c α)) + D00.Li2 ((1 + c α) / 2))

/-- The corner integral `M(u;α) = ∫_0^1 [(1+g)^{-1-u} √(1+g^{2-2u}-2cos α g^{1-u}) - 1] g^{u-1} dg`
(lines 1468-1471). -/
noncomputable def Mint (u α : ℝ) : ℝ :=
  ∫ g in (0 : ℝ)..1,
    ((1 + g) ^ (-1 - u) * Real.sqrt (1 + g ^ (2 - 2 * u) - 2 * Real.cos α * g ^ (1 - u)) - 1)
      * g ^ (u - 1)

/-- `thm:Malpha` = `prop:M1short` on their common part: for `0 < α < π` the 19+8-dilogarithm
table form and the four-dilogarithm form of `M_1(α)` are equal. -/
def stmt_1 : Prop :=
  ∀ α : ℝ, 0 < α → α < π → Mtable α = Mshort α

/-- `thm:Malpha`: "Every `Li_2` there is the real dilogarithm at a point of `(-∞,1)`". -/
def stmt_2 : Prop :=
  ∀ α : ℝ, 0 < α → α < π → ∀ p ∈ uTerms α ++ vTerms (c α), p.2 < 1

/-- `thm:Malpha`: "every logarithm is of a positive number" (the eight logarithm arguments of
`𝓛`). -/
def stmt_3 : Prop :=
  ∀ α : ℝ, 0 < α → α < π →
    0 < 1 + s α ∧ 0 < c α / (1 + s α) ∧ 0 < c α / (1 + c α) ∧ 0 < (1 + s α) / s α ∧
    (0 : ℝ) < 4 ∧ 0 < Xm α ∧ 0 < Xp α ∧ 0 < (1 + c α) * (1 + s α + c α)

/-- `thm:Malpha`: `M_1(α) = ζ(2) + U + V + 𝓛`, with `M_1(α) = ∂_u M(u;α)` at `u = 0`. -/
def stmt_4 : Prop :=
  ∀ α : ℝ, 0 < α → α < π → HasDerivAt (fun u => Mint u α) (Mtable α) 0

/-- `prop:M1short`: `M_1(α)` equals the four-dilogarithm form, with `M_1(α) = ∂_u M(u;α)` at
`u = 0`. -/
def stmt_5 : Prop :=
  ∀ α : ℝ, 0 < α → α < π → HasDerivAt (fun u => Mint u α) (Mshort α) 0

end D04
-- END STATEMENTS

namespace D04

lemma cs_facts {α : ℝ} (h0 : 0 < α) (h1 : α < Real.pi) :
    0 < c α ∧ c α < 1 ∧ 0 < s α ∧ s α < 1 ∧ c α ^ 2 + s α ^ 2 = 1 := by
  have ha : 0 < α / 2 := by linarith
  have hb : α / 2 < Real.pi / 2 := by linarith
  have hc : 0 < c α := Real.cos_pos_of_mem_Ioo ⟨by linarith [Real.pi_pos], hb⟩
  have hs : 0 < s α := Real.sin_pos_of_pos_of_lt_pi ha (by linarith [Real.pi_pos])
  have hcs : c α ^ 2 + s α ^ 2 = 1 := by unfold c s; rw [Real.cos_sq_add_sin_sq]
  exact ⟨hc, by nlinarith, hs, by nlinarith, hcs⟩

end D04


namespace D04

open Real


noncomputable def mono (k : ℝ) (a : ℕ) (b c d : ℤ) (t : ℝ) : ℝ :=
  k * t ^ a * (1 - t) ^ b * (1 + t) ^ c * (1 + t ^ 2) ^ d

noncomputable def ldr (a : ℕ) (b c d : ℤ) (t : ℝ) : ℝ :=
  a / t - b / (1 - t) + c / (1 + t) + 2 * d * t / (1 + t ^ 2)

lemma hasDerivAt_mono (k : ℝ) (a : ℕ) (b c d : ℤ) {t : ℝ} (h0 : 0 < t) (h1 : t < 1) :
    HasDerivAt (mono k a b c d) (mono k a b c d t * ldr a b c d t) t := by
  have ht : t ≠ 0 := h0.ne'
  have hm : (1 - t) ≠ 0 := by linarith
  have hp : (1 + t) ≠ 0 := by linarith
  have hq : (1 + t ^ 2) ≠ 0 := by positivity
  have e1 := hasDerivAt_zpow (a : ℤ) t (Or.inl ht)
  have e2 := (hasDerivAt_zpow b (1 - t) (Or.inl hm)).comp t ((hasDerivAt_id t).const_sub 1)
  have e3 := (hasDerivAt_zpow c (1 + t) (Or.inl hp)).comp t ((hasDerivAt_id t).const_add 1)
  have e4 := (hasDerivAt_zpow d (1 + t ^ 2) (Or.inl hq)).comp t
    (((hasDerivAt_id t).pow 2).const_add 1)
  have := (((e1.const_mul k).mul e2).mul e3).mul e4
  have hf : mono k a b c d = ((fun y => k * y ^ (a : ℤ)) * (fun x => x ^ b) ∘ HSub.hSub 1 *
      (fun x => x ^ c) ∘ HAdd.hAdd 1 * (fun x => x ^ d) ∘ fun x => 1 + (id ^ 2) x) := by
    ext x; simp [mono, Function.comp, zpow_natCast]
  rw [hf]
  refine this.congr_deriv ?_
  simp only [mono, ldr, zpow_sub_one₀ ht, zpow_sub_one₀ hm, zpow_sub_one₀ hp, zpow_sub_one₀ hq,
    Function.comp_apply, id, zpow_natCast, Pi.mul_apply, Pi.pow_apply, Int.cast_natCast]
  norm_num
  field_simp
  ring

lemma mono_pos {k : ℝ} (hk : 0 < k) (a : ℕ) (b c d : ℤ) {t : ℝ} (h0 : 0 < t) (h1 : t < 1) :
    0 < mono k a b c d t := by
  unfold mono
  have : 0 < 1 - t := by linarith
  have : 0 < 1 + t := by linarith
  positivity

lemma mono_ne {k : ℝ} (hk : k ≠ 0) (a : ℕ) (b c d : ℤ) {t : ℝ} (h0 : 0 < t) (h1 : t < 1) :
    mono k a b c d t ≠ 0 := by
  unfold mono
  have : (1 - t) ≠ 0 := by linarith
  have : (1 + t) ≠ 0 := by linarith
  have : t ≠ 0 := h0.ne'
  have : (1 + t^2) ≠ 0 := by positivity
  exact mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero hk (pow_ne_zero _ h0.ne'))
    (zpow_ne_zero _ (by linarith))) (zpow_ne_zero _ (by linarith))) (zpow_ne_zero _ (by positivity))

lemma log_mono {k : ℝ} (hk : 0 < k) (a : ℕ) (b c d : ℤ) {t : ℝ} (h0 : 0 < t) (h1 : t < 1) :
    Real.log (mono k a b c d t) = Real.log k + a * Real.log t + b * Real.log (1 - t)
      + c * Real.log (1 + t) + d * Real.log (1 + t ^ 2) := by
  have : 0 < 1 - t := by linarith
  have : 0 < 1 + t := by linarith
  have : 0 < 1 + t ^ 2 := by positivity
  unfold mono
  rw [Real.log_mul (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
    Real.log_pow, Real.log_zpow, Real.log_zpow, Real.log_zpow]

lemma hasDerivAt_Li2_mono {k k' : ℝ} {a a' : ℕ} {b c d b' c' d' : ℤ} {t : ℝ} (h0 : 0 < t)
    (h1 : t < 1) (hk : k ≠ 0) (hk' : 0 < k')
    (h : 1 - mono k a b c d t = mono k' a' b' c' d' t) :
    HasDerivAt (fun x => D00.Li2 (mono k a b c d x))
      (-(Real.log k' + a' * Real.log t + b' * Real.log (1 - t) + c' * Real.log (1 + t)
        + d' * Real.log (1 + t ^ 2)) * ldr a b c d t) t := by
  have hlt : mono k a b c d t < 1 := by
    have := mono_pos hk' a' b' c' d' h0 h1; linarith
  have hne := mono_ne hk a b c d h0 h1
  have := (D00.hasDerivAt_Li2 hlt hne).comp t (hasDerivAt_mono k a b c d h0 h1)
  refine this.congr_deriv ?_
  rw [h, log_mono hk' _ _ _ _ h0 h1]
  field_simp

/-- Auxiliary function for the constant part. -/
noncomputable def HA (x : ℝ) : ℝ :=
    (1 : ℝ) * D00.Li2 (mono (-2 : ℝ) 1 0 0 (-1) x) +
    (-1 : ℝ) * D00.Li2 (mono (-1 : ℝ) 0 (-1) 1 0 x) +
    (-1 : ℝ) * D00.Li2 (mono (-2 : ℝ) 1 (-2) 0 0 x) +
    (1 : ℝ) * D00.Li2 (mono (-1 : ℝ) 1 (-1) 1 0 x) +
    (-2 : ℝ) * D00.Li2 (mono (1 : ℝ) 0 0 (-2) 1 x) +
    (1 : ℝ) * D00.Li2 (mono (1 : ℝ) 0 2 0 (-1) x) +
    (1 : ℝ) * D00.Li2 (mono (1 : ℝ) 0 1 (-1) 0 x) +
    (1 : ℝ) * D00.Li2 (mono (2 : ℝ) 1 0 (-2) 0 x) +
    (-1 : ℝ) * D00.Li2 (mono (1 : ℝ) 1 1 (-1) 0 x) +
    (1 : ℝ) * D00.Li2 (mono (-1 : ℝ) 0 (-1) (-1) 1 x) +
    (1 : ℝ) * D00.Li2 (mono (-4 : ℝ) 2 (-2) (-2) 0 x) +
    (-1 : ℝ) * D00.Li2 (mono (-2 : ℝ) 2 (-1) (-1) 0 x) +
    (-1 : ℝ) * D00.Li2 (mono (1 : ℝ) 0 2 2 (-2) x) +
    (1 : ℝ) * D00.Li2 (mono ((1 : ℝ) / 2) 0 1 1 0 x) +
    (1 : ℝ) * D00.Li2 (mono (1 : ℝ) 2 0 0 (-1) x) +
    (1 : ℝ) * D00.Li2 (mono (-1 : ℝ) 0 1 1 (-1) x) +
    (1 : ℝ) * D00.Li2 (mono (1 : ℝ) 0 0 0 (-1) x) +
    ((1 : ℝ) / 4) * D00.Li2 (mono (1 : ℝ) 4 0 0 0 x) +
    ((1 : ℝ) / 4) * D00.Li2 (mono (1 : ℝ) 0 1 1 1 x) +
    (4 : ℝ) * (Real.log (1 - x) * Real.log (1 + x)) +
    (-1 : ℝ) * (Real.log (1 - x) * Real.log (1 + x ^ 2)) +
    (-1 : ℝ) * (Real.log (1 - x) * Real.log 2) +
    (1 : ℝ) * (Real.log (1 + x) * Real.log (1 + x ^ 2)) +
    (-3 : ℝ) * (Real.log (1 + x) * Real.log 2) +
    ((1 : ℝ) / 2) * (Real.log (1 + x ^ 2) * Real.log (1 + x ^ 2)) +
    (1 : ℝ) * (Real.log (1 + x ^ 2) * Real.log 2) +
    ((1 : ℝ) / 2) * (Real.log 2 * Real.log 2) +
    ((-1 : ℝ) / 4) * (π ^ 2 / 6)

lemma hasDerivAt_HA {t : ℝ} (h0 : 0 < t) (h1 : t < 1) : HasDerivAt HA 0 t := by
  have hm : (1 - t) ≠ 0 := by linarith
  have hp : (1 + t) ≠ 0 := by linarith
  have hq : (1 + t ^ 2) ≠ 0 := by positivity
  have ht : t ≠ 0 := h0.ne'
  have L4 : Real.log (4:ℝ) = 2 * Real.log 2 := by
    rw [show (4:ℝ) = 2 ^ 2 by norm_num, Real.log_pow]; push_cast; ring
  have Lh : Real.log ((1:ℝ) / 2) = -Real.log 2 := by rw [one_div, Real.log_inv]
  have d1 : HasDerivAt (fun x => Real.log (1 - x)) (-1 / (1 - t)) t := ((hasDerivAt_id t).const_sub 1).log hm
  have d2 : HasDerivAt (fun x => Real.log (1 + x)) (1 / (1 + t)) t := ((hasDerivAt_id t).const_add 1).log hp
  have d3 : HasDerivAt (fun x => Real.log (1 + x ^ 2)) ((2 * t) / (1 + t ^ 2)) t := by
    have := (((hasDerivAt_id t).pow 2).const_add 1).log hq
    simpa using this
  have d4 : HasDerivAt (fun _ : ℝ => Real.log 2) 0 t := hasDerivAt_const _ _
  have l0 := hasDerivAt_Li2_mono (k := (-2 : ℝ)) (k' := (1 : ℝ)) (a := 1) (b := 0) (c := 0) (d := (-1)) (a' := 0) (b' := 0) (c' := 2) (d' := (-1)) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l1 := hasDerivAt_Li2_mono (k := (-1 : ℝ)) (k' := (2 : ℝ)) (a := 0) (b := (-1)) (c := 1) (d := 0) (a' := 0) (b' := (-1)) (c' := 0) (d' := 0) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l2 := hasDerivAt_Li2_mono (k := (-2 : ℝ)) (k' := (1 : ℝ)) (a := 1) (b := (-2)) (c := 0) (d := 0) (a' := 0) (b' := (-2)) (c' := 0) (d' := 1) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l3 := hasDerivAt_Li2_mono (k := (-1 : ℝ)) (k' := (1 : ℝ)) (a := 1) (b := (-1)) (c := 1) (d := 0) (a' := 0) (b' := (-1)) (c' := 0) (d' := 1) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l4 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (2 : ℝ)) (a := 0) (b := 0) (c := (-2)) (d := 1) (a' := 1) (b' := 0) (c' := (-2)) (d' := 0) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l5 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (2 : ℝ)) (a := 0) (b := 2) (c := 0) (d := (-1)) (a' := 1) (b' := 0) (c' := 0) (d' := (-1)) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l6 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (2 : ℝ)) (a := 0) (b := 1) (c := (-1)) (d := 0) (a' := 1) (b' := 0) (c' := (-1)) (d' := 0) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l7 := hasDerivAt_Li2_mono (k := (2 : ℝ)) (k' := (1 : ℝ)) (a := 1) (b := 0) (c := (-2)) (d := 0) (a' := 0) (b' := 0) (c' := (-2)) (d' := 1) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l8 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (1 : ℝ)) (a := 1) (b := 1) (c := (-1)) (d := 0) (a' := 0) (b' := 0) (c' := (-1)) (d' := 1) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l9 := hasDerivAt_Li2_mono (k := (-1 : ℝ)) (k' := (2 : ℝ)) (a := 0) (b := (-1)) (c := (-1)) (d := 1) (a' := 0) (b' := (-1)) (c' := (-1)) (d' := 0) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l10 := hasDerivAt_Li2_mono (k := (-4 : ℝ)) (k' := (1 : ℝ)) (a := 2) (b := (-2)) (c := (-2)) (d := 0) (a' := 0) (b' := (-2)) (c' := (-2)) (d' := 2) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l11 := hasDerivAt_Li2_mono (k := (-2 : ℝ)) (k' := (1 : ℝ)) (a := 2) (b := (-1)) (c := (-1)) (d := 0) (a' := 0) (b' := (-1)) (c' := (-1)) (d' := 1) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l12 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (4 : ℝ)) (a := 0) (b := 2) (c := 2) (d := (-2)) (a' := 2) (b' := 0) (c' := 0) (d' := (-2)) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l13 := hasDerivAt_Li2_mono (k := ((1 : ℝ) / 2)) (k' := ((1 : ℝ) / 2)) (a := 0) (b := 1) (c := 1) (d := 0) (a' := 0) (b' := 0) (c' := 0) (d' := 1) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l14 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (1 : ℝ)) (a := 2) (b := 0) (c := 0) (d := (-1)) (a' := 0) (b' := 0) (c' := 0) (d' := (-1)) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l15 := hasDerivAt_Li2_mono (k := (-1 : ℝ)) (k' := (2 : ℝ)) (a := 0) (b := 1) (c := 1) (d := (-1)) (a' := 0) (b' := 0) (c' := 0) (d' := (-1)) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l16 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (1 : ℝ)) (a := 0) (b := 0) (c := 0) (d := (-1)) (a' := 2) (b' := 0) (c' := 0) (d' := (-1)) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l17 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (1 : ℝ)) (a := 4) (b := 0) (c := 0) (d := 0) (a' := 0) (b' := 1) (c' := 1) (d' := 1) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l18 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (1 : ℝ)) (a := 0) (b := 1) (c := 1) (d := 1) (a' := 4) (b' := 0) (c' := 0) (d' := 0) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have hall := ((((((((((((((((((((((((((((l0.const_mul (1 : ℝ)).add (l1.const_mul (-1 : ℝ))).add (l2.const_mul (-1 : ℝ))).add (l3.const_mul (1 : ℝ))).add (l4.const_mul (-2 : ℝ))).add (l5.const_mul (1 : ℝ))).add (l6.const_mul (1 : ℝ))).add (l7.const_mul (1 : ℝ))).add (l8.const_mul (-1 : ℝ))).add (l9.const_mul (1 : ℝ))).add (l10.const_mul (1 : ℝ))).add (l11.const_mul (-1 : ℝ))).add (l12.const_mul (-1 : ℝ))).add (l13.const_mul (1 : ℝ))).add (l14.const_mul (1 : ℝ))).add (l15.const_mul (1 : ℝ))).add (l16.const_mul (1 : ℝ))).add (l17.const_mul ((1 : ℝ) / 4))).add (l18.const_mul ((1 : ℝ) / 4))).add ((d1.mul d2).const_mul (4 : ℝ))).add ((d1.mul d3).const_mul (-1 : ℝ))).add ((d1.mul d4).const_mul (-1 : ℝ))).add ((d2.mul d3).const_mul (1 : ℝ))).add ((d2.mul d4).const_mul (-3 : ℝ))).add ((d3.mul d3).const_mul ((1 : ℝ) / 2))).add ((d3.mul d4).const_mul (1 : ℝ))).add ((d4.mul d4).const_mul ((1 : ℝ) / 2)))).add_const (((-1 : ℝ) / 4) * (π ^ 2 / 6))
  refine hall.congr_deriv ?_
  simp only [ldr, Real.log_one, L4, Lh]
  push_cast
  field_simp
  ring

/-- Auxiliary function for the `c` part. -/
noncomputable def HB (x : ℝ) : ℝ :=
    (1 : ℝ) * D00.Li2 (mono (-2 : ℝ) 1 (-1) 0 0 x) +
    (-1 : ℝ) * D00.Li2 (mono (-1 : ℝ) 1 1 0 (-1) x) +
    (-2 : ℝ) * D00.Li2 (mono (-1 : ℝ) 1 (-1) 0 0 x) +
    (-1 : ℝ) * D00.Li2 (mono (1 : ℝ) 1 0 (-1) 0 x) +
    (1 : ℝ) * D00.Li2 (mono (1 : ℝ) 0 0 (-1) 0 x) +
    (1 : ℝ) * D00.Li2 (mono ((1 : ℝ) / 2) 0 1 0 0 x) +
    (1 : ℝ) * D00.Li2 (mono (1 : ℝ) 0 1 0 (-1) x) +
    (-2 : ℝ) * D00.Li2 (mono (1 : ℝ) 0 1 0 0 x) +
    (-1 : ℝ) * D00.Li2 (mono (-1 : ℝ) 0 (-1) (-1) 1 x) +
    (-1 : ℝ) * D00.Li2 (mono (-2 : ℝ) 2 (-1) (-1) 0 x) +
    (2 : ℝ) * D00.Li2 (mono (-1 : ℝ) 2 (-1) (-1) 0 x) +
    ((-3 : ℝ) / 2) * D00.Li2 (mono (1 : ℝ) 0 1 1 (-1) x) +
    (-1 : ℝ) * D00.Li2 (mono ((1 : ℝ) / 2) 0 1 1 0 x) +
    (2 : ℝ) * D00.Li2 (mono (1 : ℝ) 0 1 1 0 x) +
    (-1 : ℝ) * D00.Li2 (mono (1 : ℝ) 2 0 0 (-1) x) +
    (1 : ℝ) * D00.Li2 (mono (-1 : ℝ) 0 1 1 (-1) x) +
    (1 : ℝ) * D00.Li2 (mono (1 : ℝ) 0 0 0 (-1) x) +
    ((1 : ℝ) / 2) * D00.Li2 (mono (2 : ℝ) 2 0 0 (-1) x) +
    ((-1 : ℝ) / 2) * (Real.log (1 - x) * Real.log (1 - x)) +
    ((1 : ℝ) / 2) * (Real.log (1 - x) * Real.log (1 + x ^ 2)) +
    ((1 : ℝ) / 2) * (Real.log (1 - x) * Real.log 2) +
    ((1 : ℝ) / 2) * (Real.log (1 + x) * Real.log (1 + x ^ 2)) +
    ((3 : ℝ) / 2) * (Real.log (1 + x) * Real.log 2) +
    ((-1 : ℝ) / 2) * (Real.log (1 + x ^ 2) * Real.log 2) +
    ((-3 : ℝ) / 2) * (π ^ 2 / 6)

lemma hasDerivAt_HB {t : ℝ} (h0 : 0 < t) (h1 : t < 1) : HasDerivAt HB 0 t := by
  have hm : (1 - t) ≠ 0 := by linarith
  have hp : (1 + t) ≠ 0 := by linarith
  have hq : (1 + t ^ 2) ≠ 0 := by positivity
  have ht : t ≠ 0 := h0.ne'
  have L4 : Real.log (4:ℝ) = 2 * Real.log 2 := by
    rw [show (4:ℝ) = 2 ^ 2 by norm_num, Real.log_pow]; push_cast; ring
  have Lh : Real.log ((1:ℝ) / 2) = -Real.log 2 := by rw [one_div, Real.log_inv]
  have d1 : HasDerivAt (fun x => Real.log (1 - x)) (-1 / (1 - t)) t := ((hasDerivAt_id t).const_sub 1).log hm
  have d2 : HasDerivAt (fun x => Real.log (1 + x)) (1 / (1 + t)) t := ((hasDerivAt_id t).const_add 1).log hp
  have d3 : HasDerivAt (fun x => Real.log (1 + x ^ 2)) ((2 * t) / (1 + t ^ 2)) t := by
    have := (((hasDerivAt_id t).pow 2).const_add 1).log hq
    simpa using this
  have d4 : HasDerivAt (fun _ : ℝ => Real.log 2) 0 t := hasDerivAt_const _ _
  have l0 := hasDerivAt_Li2_mono (k := (-2 : ℝ)) (k' := (1 : ℝ)) (a := 1) (b := (-1)) (c := 0) (d := 0) (a' := 0) (b' := (-1)) (c' := 1) (d' := 0) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l1 := hasDerivAt_Li2_mono (k := (-1 : ℝ)) (k' := (1 : ℝ)) (a := 1) (b := 1) (c := 0) (d := (-1)) (a' := 0) (b' := 0) (c' := 1) (d' := (-1)) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l2 := hasDerivAt_Li2_mono (k := (-1 : ℝ)) (k' := (1 : ℝ)) (a := 1) (b := (-1)) (c := 0) (d := 0) (a' := 0) (b' := (-1)) (c' := 0) (d' := 0) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l3 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (1 : ℝ)) (a := 1) (b := 0) (c := (-1)) (d := 0) (a' := 0) (b' := 0) (c' := (-1)) (d' := 0) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l4 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (1 : ℝ)) (a := 0) (b := 0) (c := (-1)) (d := 0) (a' := 1) (b' := 0) (c' := (-1)) (d' := 0) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l5 := hasDerivAt_Li2_mono (k := ((1 : ℝ) / 2)) (k' := ((1 : ℝ) / 2)) (a := 0) (b := 1) (c := 0) (d := 0) (a' := 0) (b' := 0) (c' := 1) (d' := 0) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l6 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (1 : ℝ)) (a := 0) (b := 1) (c := 0) (d := (-1)) (a' := 1) (b' := 0) (c' := 1) (d' := (-1)) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l7 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (1 : ℝ)) (a := 0) (b := 1) (c := 0) (d := 0) (a' := 1) (b' := 0) (c' := 0) (d' := 0) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l8 := hasDerivAt_Li2_mono (k := (-1 : ℝ)) (k' := (2 : ℝ)) (a := 0) (b := (-1)) (c := (-1)) (d := 1) (a' := 0) (b' := (-1)) (c' := (-1)) (d' := 0) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l9 := hasDerivAt_Li2_mono (k := (-2 : ℝ)) (k' := (1 : ℝ)) (a := 2) (b := (-1)) (c := (-1)) (d := 0) (a' := 0) (b' := (-1)) (c' := (-1)) (d' := 1) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l10 := hasDerivAt_Li2_mono (k := (-1 : ℝ)) (k' := (1 : ℝ)) (a := 2) (b := (-1)) (c := (-1)) (d := 0) (a' := 0) (b' := (-1)) (c' := (-1)) (d' := 0) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l11 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (2 : ℝ)) (a := 0) (b := 1) (c := 1) (d := (-1)) (a' := 2) (b' := 0) (c' := 0) (d' := (-1)) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l12 := hasDerivAt_Li2_mono (k := ((1 : ℝ) / 2)) (k' := ((1 : ℝ) / 2)) (a := 0) (b := 1) (c := 1) (d := 0) (a' := 0) (b' := 0) (c' := 0) (d' := 1) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l13 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (1 : ℝ)) (a := 0) (b := 1) (c := 1) (d := 0) (a' := 2) (b' := 0) (c' := 0) (d' := 0) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l14 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (1 : ℝ)) (a := 2) (b := 0) (c := 0) (d := (-1)) (a' := 0) (b' := 0) (c' := 0) (d' := (-1)) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l15 := hasDerivAt_Li2_mono (k := (-1 : ℝ)) (k' := (2 : ℝ)) (a := 0) (b := 1) (c := 1) (d := (-1)) (a' := 0) (b' := 0) (c' := 0) (d' := (-1)) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l16 := hasDerivAt_Li2_mono (k := (1 : ℝ)) (k' := (1 : ℝ)) (a := 0) (b := 0) (c := 0) (d := (-1)) (a' := 2) (b' := 0) (c' := 0) (d' := (-1)) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have l17 := hasDerivAt_Li2_mono (k := (2 : ℝ)) (k' := (1 : ℝ)) (a := 2) (b := 0) (c := 0) (d := (-1)) (a' := 0) (b' := 1) (c' := 1) (d' := (-1)) h0 h1 (by norm_num) (by norm_num) (by simp only [mono, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul]; field_simp; ring)
  have hall := (((((((((((((((((((((((((l0.const_mul (1 : ℝ)).add (l1.const_mul (-1 : ℝ))).add (l2.const_mul (-2 : ℝ))).add (l3.const_mul (-1 : ℝ))).add (l4.const_mul (1 : ℝ))).add (l5.const_mul (1 : ℝ))).add (l6.const_mul (1 : ℝ))).add (l7.const_mul (-2 : ℝ))).add (l8.const_mul (-1 : ℝ))).add (l9.const_mul (-1 : ℝ))).add (l10.const_mul (2 : ℝ))).add (l11.const_mul ((-3 : ℝ) / 2))).add (l12.const_mul (-1 : ℝ))).add (l13.const_mul (2 : ℝ))).add (l14.const_mul (-1 : ℝ))).add (l15.const_mul (1 : ℝ))).add (l16.const_mul (1 : ℝ))).add (l17.const_mul ((1 : ℝ) / 2))).add ((d1.mul d1).const_mul ((-1 : ℝ) / 2))).add ((d1.mul d3).const_mul ((1 : ℝ) / 2))).add ((d1.mul d4).const_mul ((1 : ℝ) / 2))).add ((d2.mul d3).const_mul ((1 : ℝ) / 2))).add ((d2.mul d4).const_mul ((3 : ℝ) / 2))).add ((d3.mul d4).const_mul ((-1 : ℝ) / 2)))).add_const (((-3 : ℝ) / 2) * (π ^ 2 / 6))
  refine hall.congr_deriv ?_
  simp only [ldr, Real.log_one, L4, Lh]
  push_cast
  field_simp
  ring

lemma continuousAt_mono (k : ℝ) (a : ℕ) (b c d : ℤ) {x : ℝ} (h0 : -1 < x) (h1 : x < 1) :
    ContinuousAt (mono k a b c d) x := by
  have hm : (1 - x) ≠ 0 := by linarith
  have hp : (1 + x) ≠ 0 := by linarith
  have hq : (1 + x ^ 2) ≠ 0 := by positivity
  exact (((continuousAt_const.mul (continuousAt_id.pow a)).mul
    ((continuousAt_const.sub continuousAt_id).zpow₀ b (Or.inl hm))).mul
    ((continuousAt_const.add continuousAt_id).zpow₀ c (Or.inl hp))).mul
    ((continuousAt_const.add (continuousAt_id.pow 2)).zpow₀ d (Or.inl hq))

lemma continuousAt_HA {x : ℝ} (h0 : -1 < x) (h1 : x < 1) : ContinuousAt HA x := by
  have hm : (1 - x) ≠ 0 := by linarith
  have hp : (1 + x) ≠ 0 := by linarith
  have hq : (1 + x ^ 2) ≠ 0 := by positivity
  have c1 : ContinuousAt (fun x => Real.log (1 - x)) x := (continuousAt_const.sub continuousAt_id).log hm
  have c2 : ContinuousAt (fun x => Real.log (1 + x)) x := (continuousAt_const.add continuousAt_id).log hp
  have c3 : ContinuousAt (fun x => Real.log (1 + x ^ 2)) x :=
    (continuousAt_const.add (continuousAt_id.pow 2)).log hq
  have c4 : ContinuousAt (fun _ : ℝ => Real.log 2) x := continuousAt_const
  exact (((((((((((((((((((((((((((((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (-2 : ℝ) 1 0 0 (-1) h0 h1)).const_mul (1 : ℝ)).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (-1 : ℝ) 0 (-1) 1 0 h0 h1)).const_mul (-1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (-2 : ℝ) 1 (-2) 0 0 h0 h1)).const_mul (-1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (-1 : ℝ) 1 (-1) 1 0 h0 h1)).const_mul (1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 0 0 (-2) 1 h0 h1)).const_mul (-2 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 0 2 0 (-1) h0 h1)).const_mul (1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 0 1 (-1) 0 h0 h1)).const_mul (1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (2 : ℝ) 1 0 (-2) 0 h0 h1)).const_mul (1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 1 1 (-1) 0 h0 h1)).const_mul (-1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (-1 : ℝ) 0 (-1) (-1) 1 h0 h1)).const_mul (1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (-4 : ℝ) 2 (-2) (-2) 0 h0 h1)).const_mul (1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (-2 : ℝ) 2 (-1) (-1) 0 h0 h1)).const_mul (-1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 0 2 2 (-2) h0 h1)).const_mul (-1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono ((1 : ℝ) / 2) 0 1 1 0 h0 h1)).const_mul (1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 2 0 0 (-1) h0 h1)).const_mul (1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (-1 : ℝ) 0 1 1 (-1) h0 h1)).const_mul (1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 0 0 0 (-1) h0 h1)).const_mul (1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 4 0 0 0 h0 h1)).const_mul ((1 : ℝ) / 4))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 0 1 1 1 h0 h1)).const_mul ((1 : ℝ) / 4))).add
      ((c1.mul c2).const_mul (4 : ℝ))).add
      ((c1.mul c3).const_mul (-1 : ℝ))).add
      ((c1.mul c4).const_mul (-1 : ℝ))).add
      ((c2.mul c3).const_mul (1 : ℝ))).add
      ((c2.mul c4).const_mul (-3 : ℝ))).add
      ((c3.mul c3).const_mul ((1 : ℝ) / 2))).add
      ((c3.mul c4).const_mul (1 : ℝ))).add
      ((c4.mul c4).const_mul ((1 : ℝ) / 2)))).add continuousAt_const

lemma HA_zero : HA 0 = 0 := by
  have hh : D00.Li2 (2⁻¹) = π ^ 2 / 12 - (1 / 2) * Real.log 2 ^ 2 := by
    rw [← D00.Li2_half]; norm_num
  norm_num [HA, mono, hh, D00.Li2_half]
  ring

lemma HA_eq_zero {t : ℝ} (h0 : 0 < t) (h1 : t < 1) : HA t = 0 := by
  have := D00.eq_of_deriv_zero (f := HA) (a := 0) (b := t) h0.le
    (fun x hx => (continuousAt_HA (by linarith [hx.1]) (by linarith [hx.2])).continuousWithinAt)
    (fun x hx => hasDerivAt_HA hx.1 (by linarith [hx.2]))
  rw [← this, HA_zero]

lemma continuousAt_HB {x : ℝ} (h0 : -1 < x) (h1 : x < 1) : ContinuousAt HB x := by
  have hm : (1 - x) ≠ 0 := by linarith
  have hp : (1 + x) ≠ 0 := by linarith
  have hq : (1 + x ^ 2) ≠ 0 := by positivity
  have c1 : ContinuousAt (fun x => Real.log (1 - x)) x := (continuousAt_const.sub continuousAt_id).log hm
  have c2 : ContinuousAt (fun x => Real.log (1 + x)) x := (continuousAt_const.add continuousAt_id).log hp
  have c3 : ContinuousAt (fun x => Real.log (1 + x ^ 2)) x :=
    (continuousAt_const.add (continuousAt_id.pow 2)).log hq
  have c4 : ContinuousAt (fun _ : ℝ => Real.log 2) x := continuousAt_const
  exact ((((((((((((((((((((((((((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (-2 : ℝ) 1 (-1) 0 0 h0 h1)).const_mul (1 : ℝ)).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (-1 : ℝ) 1 1 0 (-1) h0 h1)).const_mul (-1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (-1 : ℝ) 1 (-1) 0 0 h0 h1)).const_mul (-2 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 1 0 (-1) 0 h0 h1)).const_mul (-1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 0 0 (-1) 0 h0 h1)).const_mul (1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono ((1 : ℝ) / 2) 0 1 0 0 h0 h1)).const_mul (1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 0 1 0 (-1) h0 h1)).const_mul (1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 0 1 0 0 h0 h1)).const_mul (-2 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (-1 : ℝ) 0 (-1) (-1) 1 h0 h1)).const_mul (-1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (-2 : ℝ) 2 (-1) (-1) 0 h0 h1)).const_mul (-1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (-1 : ℝ) 2 (-1) (-1) 0 h0 h1)).const_mul (2 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 0 1 1 (-1) h0 h1)).const_mul ((-3 : ℝ) / 2))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono ((1 : ℝ) / 2) 0 1 1 0 h0 h1)).const_mul (-1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 0 1 1 0 h0 h1)).const_mul (2 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 2 0 0 (-1) h0 h1)).const_mul (-1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (-1 : ℝ) 0 1 1 (-1) h0 h1)).const_mul (1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (1 : ℝ) 0 0 0 (-1) h0 h1)).const_mul (1 : ℝ))).add
      ((D00.continuous_Li2.continuousAt.comp (continuousAt_mono (2 : ℝ) 2 0 0 (-1) h0 h1)).const_mul ((1 : ℝ) / 2))).add
      ((c1.mul c1).const_mul ((-1 : ℝ) / 2))).add
      ((c1.mul c3).const_mul ((1 : ℝ) / 2))).add
      ((c1.mul c4).const_mul ((1 : ℝ) / 2))).add
      ((c2.mul c3).const_mul ((1 : ℝ) / 2))).add
      ((c2.mul c4).const_mul ((3 : ℝ) / 2))).add
      ((c3.mul c4).const_mul ((-1 : ℝ) / 2)))).add continuousAt_const

lemma HB_zero : HB 0 = 0 := by
  have hh : D00.Li2 (2⁻¹) = π ^ 2 / 12 - (1 / 2) * Real.log 2 ^ 2 := by
    rw [← D00.Li2_half]; norm_num
  norm_num [HB, mono, hh, D00.Li2_half]
  ring

lemma HB_eq_zero {t : ℝ} (h0 : 0 < t) (h1 : t < 1) : HB t = 0 := by
  have := D00.eq_of_deriv_zero (f := HB) (a := 0) (b := t) h0.le
    (fun x hx => (continuousAt_HB (by linarith [hx.1]) (by linarith [hx.2])).continuousWithinAt)
    (fun x hx => hasDerivAt_HB hx.1 (by linarith [hx.2]))
  rw [← this, HB_zero]

noncomputable def monoz (k : ℝ) (a b c d : ℤ) (t : ℝ) : ℝ :=
  k * t ^ a * (1 - t) ^ b * (1 + t) ^ c * (1 + t ^ 2) ^ d

lemma monoz_pos {k : ℝ} (hk : 0 < k) (a b c d : ℤ) {t : ℝ} (h0 : 0 < t) (h1 : t < 1) :
    0 < monoz k a b c d t := by
  unfold monoz
  have : 0 < 1 - t := by linarith
  have : 0 < 1 + t := by linarith
  positivity

lemma log_monoz {k : ℝ} (hk : 0 < k) (a b c d : ℤ) {t : ℝ} (h0 : 0 < t) (h1 : t < 1) :
    Real.log (monoz k a b c d t) = Real.log k + a * Real.log t + b * Real.log (1 - t)
      + c * Real.log (1 + t) + d * Real.log (1 + t ^ 2) := by
  have : 0 < 1 - t := by linarith
  have : 0 < 1 + t := by linarith
  have : 0 < 1 + t ^ 2 := by positivity
  unfold monoz
  rw [Real.log_mul (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
    Real.log_zpow, Real.log_zpow, Real.log_zpow, Real.log_zpow]

lemma tan_param {α : ℝ} (h0 : 0 < α) (h1 : α < Real.pi) :
    ∃ t : ℝ, 0 < t ∧ t < 1 ∧ c α = (1 - t ^ 2) / (1 + t ^ 2) ∧ s α = 2 * t / (1 + t ^ 2) := by
  obtain ⟨hc0, -, -, -, -⟩ := cs_facts h0 h1
  have hcos : 0 < Real.cos (α / 4) :=
    Real.cos_pos_of_mem_Ioo ⟨by linarith [Real.pi_pos], by linarith [Real.pi_pos]⟩
  have hsin : 0 < Real.sin (α / 4) := Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith)
  have hsc := Real.sin_sq_add_cos_sq (α / 4)
  have hc : c α = (1 - (Real.sin (α / 4) / Real.cos (α / 4)) ^ 2) /
      (1 + (Real.sin (α / 4) / Real.cos (α / 4)) ^ 2) := by
    unfold c
    rw [show α / 2 = 2 * (α / 4) by ring, Real.cos_two_mul]
    field_simp
    nlinarith
  have hs : s α = 2 * (Real.sin (α / 4) / Real.cos (α / 4)) /
      (1 + (Real.sin (α / 4) / Real.cos (α / 4)) ^ 2) := by
    unfold s
    rw [show α / 2 = 2 * (α / 4) by ring, Real.sin_two_mul]
    field_simp
    nlinarith
  refine ⟨Real.sin (α / 4) / Real.cos (α / 4), by positivity, ?_, hc, hs⟩
  set t := Real.sin (α / 4) / Real.cos (α / 4)
  have ht : 0 < t := by positivity
  rw [hc] at hc0
  have : 0 < 1 - t ^ 2 := by
    have hp : 0 < 1 + t ^ 2 := by positivity
    have := mul_pos hc0 hp
    rwa [div_mul_cancel₀ _ hp.ne'] at this
  nlinarith

end D04

theorem D04.proof_1 : D04.stmt_1 := by
  intro α hα0 hα1
  obtain ⟨t, h0, h1, hc, hs⟩ := D04.tan_param hα0 hα1
  simp only [D04.Mtable, D04.Mshort, D04.U, D04.V, D04.Lcal, D04.M0, D04.liSum, D04.uTerms,
    D04.vTerms, D04.Xm, D04.Xp, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  generalize D04.c α = c at hc ⊢
  generalize D04.s α = s at hs ⊢
  have ht : t ≠ 0 := h0.ne'
  have hm : (1 - t) ≠ 0 := by linarith
  have hp : (1 + t) ≠ 0 := by linarith
  have hq : (1 + t ^ 2) ≠ 0 := by positivity
  have hm2 : (1 - t ^ 2) ≠ 0 := by nlinarith
  have L4 : Real.log (4:ℝ) = 2 * Real.log 2 := by
    rw [show (4:ℝ) = 2 ^ 2 by norm_num, Real.log_pow]; push_cast; ring
  have Lh : Real.log ((1:ℝ) / 2) = -Real.log 2 := by rw [one_div, Real.log_inv]
  have e1 : 1 + s - c = 2 * t * (1 + t) / (1 + t ^ 2) := by rw [hc, hs]; field_simp <;> ring
  have e2 : 1 + s + c = 2 * (1 + t) / (1 + t ^ 2) := by rw [hc, hs]; field_simp <;> ring
  have e3 : 1 + s = (1 + t) ^ 2 / (1 + t ^ 2) := by rw [hs]; field_simp <;> ring
  have e4 : 1 - c = 2 * t ^ 2 / (1 + t ^ 2) := by rw [hc]; field_simp <;> ring
  have e5 : 1 + c = 2 / (1 + t ^ 2) := by rw [hc]; field_simp <;> ring
  have i0 : D00.Li2 (-(1 / s)) + D00.Li2 (D04.mono (-2 : ℝ) 1 0 0 (-1) t) = -(Real.pi ^ 2 / 6) - 1 / 2 * (-1 * Real.log 2 + -1 * Real.log t + Real.log (1 + t ^ 2)) ^ 2 := by
    have q1 : -(1 / s) = -D04.monoz ((1 : ℝ) / 2) (-1) 0 0 1 t := by (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
    have q2 : 1 / (-(1 / s)) = D04.mono (-2 : ℝ) 1 0 0 (-1) t := by (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
    have q3 := D04.monoz_pos (k := ((1 : ℝ) / 2)) (by norm_num) (-1) 0 0 1 h0 h1
    have q4 := D00.inversion_neg (x := -(1 / s)) (by rw [q1]; linarith)
    rw [q2] at q4; rw [q4, q1, neg_neg, D04.log_monoz (by norm_num) _ _ _ _ h0 h1]
    (try simp only [Real.log_one, L4, Lh]); push_cast; ring
  have i1 : D00.Li2 (-(s / (1 - c))) + D00.Li2 (D04.mono (-1 : ℝ) 1 0 0 0 t) = -(Real.pi ^ 2 / 6) - 1 / 2 * (-1 * Real.log t) ^ 2 := by
    have q1 : -(s / (1 - c)) = -D04.monoz (1 : ℝ) (-1) 0 0 0 t := by (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
    have q2 : 1 / (-(s / (1 - c))) = D04.mono (-1 : ℝ) 1 0 0 0 t := by (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
    have q3 := D04.monoz_pos (k := (1 : ℝ)) (by norm_num) (-1) 0 0 0 h0 h1
    have q4 := D00.inversion_neg (x := -(s / (1 - c))) (by rw [q1]; linarith)
    rw [q2] at q4; rw [q4, q1, neg_neg, D04.log_monoz (by norm_num) _ _ _ _ h0 h1]
    (try simp only [Real.log_one, L4, Lh]); push_cast; ring
  have r2 : D00.Li2 (-(s / (1 + c))) = D00.Li2 (D04.mono (-1 : ℝ) 1 0 0 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r3 : D00.Li2 (-((1 + s) / c)) = D00.Li2 (D04.mono (-1 : ℝ) 0 (-1) 1 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r4 : D00.Li2 (-(s * (1 + s) / c ^ 2)) = D00.Li2 (D04.mono (-2 : ℝ) 1 (-2) 0 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r5 : D00.Li2 (-(s * (1 + s) / (c * (1 + c)))) = D00.Li2 (D04.mono (-1 : ℝ) 1 (-1) 1 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r6 : D00.Li2 (-((1 + s - c) / c)) = D00.Li2 (D04.mono (-2 : ℝ) 1 (-1) 0 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have i7 : D00.Li2 (-((1 + s - c) / (c * (1 - c)))) + D00.Li2 (D04.mono (-1 : ℝ) 1 1 0 (-1) t) = -(Real.pi ^ 2 / 6) - 1 / 2 * (-1 * Real.log t + -1 * Real.log (1 - t) + Real.log (1 + t ^ 2)) ^ 2 := by
    have q1 : -((1 + s - c) / (c * (1 - c))) = -D04.monoz (1 : ℝ) (-1) (-1) 0 1 t := by (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
    have q2 : 1 / (-((1 + s - c) / (c * (1 - c)))) = D04.mono (-1 : ℝ) 1 1 0 (-1) t := by (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
    have q3 := D04.monoz_pos (k := (1 : ℝ)) (by norm_num) (-1) (-1) 0 1 h0 h1
    have q4 := D00.inversion_neg (x := -((1 + s - c) / (c * (1 - c)))) (by rw [q1]; linarith)
    rw [q2] at q4; rw [q4, q1, neg_neg, D04.log_monoz (by norm_num) _ _ _ _ h0 h1]
    (try simp only [Real.log_one, L4, Lh]); push_cast; ring
  have r8 : D00.Li2 (-((1 + s - c) / (2 * c))) = D00.Li2 (D04.mono (-1 : ℝ) 1 (-1) 0 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r9 : D00.Li2 (1 / (1 + s)) = D00.Li2 (D04.mono (1 : ℝ) 0 0 (-2) 1 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r10 : D00.Li2 (c ^ 2 / (1 + s)) = D00.Li2 (D04.mono (1 : ℝ) 0 2 0 (-1) t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r11 : D00.Li2 (c / (1 + s)) = D00.Li2 (D04.mono (1 : ℝ) 0 1 (-1) 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r12 : D00.Li2 (s / (1 + s)) = D00.Li2 (D04.mono (2 : ℝ) 1 0 (-2) 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r13 : D00.Li2 (c * (1 - c) / (s * (1 + s))) = D00.Li2 (D04.mono (1 : ℝ) 1 1 (-1) 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r14 : D00.Li2 ((1 - c) / (1 + s - c)) = D00.Li2 (D04.mono (1 : ℝ) 1 0 (-1) 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r15 : D00.Li2 ((1 + c) / (1 + s + c)) = D00.Li2 (D04.mono (1 : ℝ) 0 0 (-1) 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r16 : D00.Li2 (c / (1 + s + c)) = D00.Li2 (D04.mono ((1 : ℝ) / 2) 0 1 0 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r17 : D00.Li2 (c * (1 + c) / (1 + s + c)) = D00.Li2 (D04.mono (1 : ℝ) 0 1 0 (-1) t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r18 : D00.Li2 (2 * c / (1 + s + c)) = D00.Li2 (D04.mono (1 : ℝ) 0 1 0 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r19 : D00.Li2 (-(1 / c)) = D00.Li2 (D04.mono (-1 : ℝ) 0 (-1) (-1) 1 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r20 : D00.Li2 (1 - 1 / c ^ 2) = D00.Li2 (D04.mono (-4 : ℝ) 2 (-2) (-2) 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r21 : D00.Li2 (-((1 - c) / c)) = D00.Li2 (D04.mono (-2 : ℝ) 2 (-1) (-1) 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r22 : D00.Li2 (-((1 - c) / (2 * c))) = D00.Li2 (D04.mono (-1 : ℝ) 2 (-1) (-1) 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r23 : D00.Li2 (c ^ 2) = D00.Li2 (D04.mono (1 : ℝ) 0 2 2 (-2) t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r24 : D00.Li2 (c) = D00.Li2 (D04.mono (1 : ℝ) 0 1 1 (-1) t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r25 : D00.Li2 (c / (1 + c)) = D00.Li2 (D04.mono ((1 : ℝ) / 2) 0 1 1 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r26 : D00.Li2 (2 * c / (1 + c)) = D00.Li2 (D04.mono (1 : ℝ) 0 1 1 0 t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r27 : D00.Li2 ((1 - c) / 2) = D00.Li2 (D04.mono (1 : ℝ) 2 0 0 (-1) t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r28 : D00.Li2 (-c) = D00.Li2 (D04.mono (-1 : ℝ) 0 1 1 (-1) t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have r29 : D00.Li2 ((1 + c) / 2) = D00.Li2 (D04.mono (1 : ℝ) 0 0 0 (-1) t) := by
    congr 1; (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring
  have g0 : Real.log (1 + s) = (2 * Real.log (1 + t) + -1 * Real.log (1 + t ^ 2)) := by
    rw [show 1 + s = D04.monoz (1 : ℝ) 0 0 2 (-1) t by (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring, D04.log_monoz (by norm_num) _ _ _ _ h0 h1]
    (try simp only [Real.log_one, L4, Lh]); push_cast; ring
  have g1 : Real.log (c / (1 + s)) = (Real.log (1 - t) + -1 * Real.log (1 + t)) := by
    rw [show c / (1 + s) = D04.monoz (1 : ℝ) 0 1 (-1) 0 t by (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring, D04.log_monoz (by norm_num) _ _ _ _ h0 h1]
    (try simp only [Real.log_one, L4, Lh]); push_cast; ring
  have g2 : Real.log (c / (1 + c)) = (-1 * Real.log 2 + Real.log (1 - t) + Real.log (1 + t)) := by
    rw [show c / (1 + c) = D04.monoz ((1 : ℝ) / 2) 0 1 1 0 t by (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring, D04.log_monoz (by norm_num) _ _ _ _ h0 h1]
    (try simp only [Real.log_one, L4, Lh]); push_cast; ring
  have g3 : Real.log ((1 + s) / s) = (-1 * Real.log 2 + -1 * Real.log t + 2 * Real.log (1 + t)) := by
    rw [show (1 + s) / s = D04.monoz ((1 : ℝ) / 2) (-1) 0 2 0 t by (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring, D04.log_monoz (by norm_num) _ _ _ _ h0 h1]
    (try simp only [Real.log_one, L4, Lh]); push_cast; ring
  have g4 : Real.log ((1 + s - c) / (1 - c)) = (-1 * Real.log t + Real.log (1 + t)) := by
    rw [show (1 + s - c) / (1 - c) = D04.monoz (1 : ℝ) (-1) 0 1 0 t by (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring, D04.log_monoz (by norm_num) _ _ _ _ h0 h1]
    (try simp only [Real.log_one, L4, Lh]); push_cast; ring
  have g5 : Real.log ((1 + s + c) / (1 + c)) = (Real.log (1 + t)) := by
    rw [show (1 + s + c) / (1 + c) = D04.monoz (1 : ℝ) 0 0 1 0 t by (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring, D04.log_monoz (by norm_num) _ _ _ _ h0 h1]
    (try simp only [Real.log_one, L4, Lh]); push_cast; ring
  have g6 : Real.log ((1 + c) * (1 + s + c)) = (2 * Real.log 2 + Real.log (1 + t) + -2 * Real.log (1 + t ^ 2)) := by
    rw [show (1 + c) * (1 + s + c) = D04.monoz (4 : ℝ) 0 0 1 (-2) t by (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring, D04.log_monoz (by norm_num) _ _ _ _ h0 h1]
    (try simp only [Real.log_one, L4, Lh]); push_cast; ring
  have g7 : Real.log (1 + c) = (1 * Real.log 2 + -1 * Real.log (1 + t ^ 2)) := by
    rw [show 1 + c = D04.monoz (2 : ℝ) 0 0 0 (-1) t by (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring, D04.log_monoz (by norm_num) _ _ _ _ h0 h1]
    (try simp only [Real.log_one, L4, Lh]); push_cast; ring
  have g8 : Real.log (1 - c) = (1 * Real.log 2 + 2 * Real.log t + -1 * Real.log (1 + t ^ 2)) := by
    rw [show 1 - c = D04.monoz (2 : ℝ) 2 0 0 (-1) t by (try rw [e1]); (try rw [e2]); (try rw [e3]); (try rw [e4]); (try rw [e5]); (try rw [hc]); (try rw [hs]); simp only [D04.mono, D04.monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring, D04.log_monoz (by norm_num) _ _ _ _ h0 h1]
    (try simp only [Real.log_one, L4, Lh]); push_cast; ring
  rw [r2, r3, r4, r5, r6, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17, r18, r19, r20, r21, r22, r23, r24, r25, r26, r27, r28, r29, g0, g1, g2, g3, g4, g5, g6, g7, g8, L4]
  have hA := D04.HA_eq_zero h0 h1
  have hB := D04.HB_eq_zero h0 h1
  unfold D04.HA at hA
  unfold D04.HB at hB
  have rA : D00.Li2 (D04.mono (1 : ℝ) 4 0 0 0 t) + D00.Li2 (D04.mono (1 : ℝ) 0 1 1 1 t) =
      Real.pi ^ 2 / 6 - (4 * Real.log t) * (Real.log (1 - t) + Real.log (1 + t) + Real.log (1 + t ^ 2)) := by
    have q : 1 - D04.mono (1 : ℝ) 4 0 0 0 t = D04.mono (1 : ℝ) 0 1 1 1 t := by
      simp only [D04.mono, zpow_one, pow_zero, zpow_zero, mul_one, one_mul]; ring
    have p1 := D04.mono_pos (k := (1 : ℝ)) (by norm_num) 4 0 0 0 h0 h1
    have p2 := D04.mono_pos (k := (1 : ℝ)) (by norm_num) 0 1 1 1 h0 h1
    have := D00.reflection p1 (by linarith)
    rw [q, D04.log_mono (by norm_num) _ _ _ _ h0 h1, D04.log_mono (by norm_num) _ _ _ _ h0 h1] at this
    rw [this]; simp only [Real.log_one]; push_cast; ring
  have rB : D00.Li2 (D04.mono (2 : ℝ) 2 0 0 (-1) t) + D00.Li2 (D04.mono (1 : ℝ) 0 1 1 (-1) t) =
      Real.pi ^ 2 / 6 - (Real.log 2 + 2 * Real.log t - Real.log (1 + t ^ 2)) *
        (Real.log (1 - t) + Real.log (1 + t) - Real.log (1 + t ^ 2)) := by
    have q : 1 - D04.mono (2 : ℝ) 2 0 0 (-1) t = D04.mono (1 : ℝ) 0 1 1 (-1) t := by
      simp only [D04.mono, zpow_neg, zpow_one, pow_zero, zpow_zero, mul_one, one_mul]
      field_simp <;> ring
    have p1 := D04.mono_pos (k := (2 : ℝ)) (by norm_num) 2 0 0 (-1) h0 h1
    have p2 := D04.mono_pos (k := (1 : ℝ)) (by norm_num) 0 1 1 (-1) h0 h1
    have := D00.reflection p1 (by linarith)
    rw [q, D04.log_mono (by norm_num) _ _ _ _ h0 h1, D04.log_mono (by norm_num) _ _ _ _ h0 h1] at this
    rw [this]; simp only [Real.log_one]; push_cast; ring
  linear_combination hA + c * hB - (1 / 4) * rA - (c / 2) * rB - i0 + i1 + c * i7

theorem D04.proof_2 : D04.stmt_2 := by
  intro α h0 h1 p hp
  obtain ⟨hc0, hc1, hs0, hs1, hcs⟩ := D04.cs_facts h0 h1
  simp only [D04.uTerms, D04.vTerms, List.mem_append, List.mem_cons, List.not_mem_nil,
    or_false] at hp
  generalize D04.c α = c at *
  generalize D04.s α = s at *
  have h1s : 0 < 1 + s := by linarith
  have h1c : 0 < 1 - c := by linarith
  have h2c : 0 < 1 + c := by linarith
  have hd : 0 < 1 + s - c := by linarith
  have hd2 : 0 < 1 + s + c := by linarith
  rcases hp with (h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h) | (h | h | h | h | h | h | h | h) <;> subst h <;> dsimp only
  · have : 0 < 1 / s := by positivity
    linarith
  · have : 0 < s / (1 - c) := by positivity
    linarith
  · have : 0 < s / (1 + c) := by positivity
    linarith
  · have : 0 < (1 + s) / c := by positivity
    linarith
  · have : 0 < s * (1 + s) / c ^ 2 := by positivity
    linarith
  · have : 0 < s * (1 + s) / (c * (1 + c)) := by positivity
    linarith
  · have : 0 < (1 + s - c) / c := by positivity
    linarith
  · have : 0 < (1 + s - c) / (c * (1 - c)) := by positivity
    linarith
  · have : 0 < (1 + s - c) / (2 * c) := by positivity
    linarith
  · rw [div_lt_one h1s]; linarith
  · rw [div_lt_one h1s]; nlinarith
  · rw [div_lt_one h1s]; linarith
  · rw [div_lt_one h1s]; linarith
  · rw [div_lt_one (by positivity)]; nlinarith [mul_pos hs0 hs0, mul_pos hc0 hs0, mul_pos hs0 hc0]
  · rw [div_lt_one hd]; linarith
  · rw [div_lt_one hd2]; linarith
  · rw [div_lt_one hd2]; linarith
  · rw [div_lt_one hd2]; nlinarith
  · rw [div_lt_one hd2]; linarith
  · have : 0 < 1 / c := by positivity
    linarith
  · have : 0 < 1 / c ^ 2 := by positivity
    have : 1 < 1 / c ^ 2 := by
      rw [lt_div_iff₀ (by positivity)]; nlinarith
    linarith
  · have : 0 < (1 - c) / c := by positivity
    linarith
  · have : 0 < (1 - c) / (2 * c) := by positivity
    linarith
  · nlinarith
  · linarith
  · rw [div_lt_one h2c]; linarith
  · rw [div_lt_one h2c]; linarith

theorem D04.proof_3 : D04.stmt_3 := by
  intro α h0 h1
  obtain ⟨hc0, hc1, hs0, hs1, hcs⟩ := D04.cs_facts h0 h1
  have hX : ∀ x y : ℝ, 0 < x → 0 < y → 0 < x / y := fun x y hx hy => div_pos hx hy
  refine ⟨by linarith, hX _ _ hc0 (by linarith), hX _ _ hc0 (by linarith),
    hX _ _ (by linarith) hs0, by norm_num, ?_, ?_, by positivity⟩
  · unfold D04.Xm
    apply hX
    · nlinarith
    · linarith
  · unfold D04.Xp
    apply hX <;> linarith


namespace D04
open Real

noncomputable def PL (ρ σ t : ℝ) : ℝ :=
  -D00.Li2 (-(t - ρ) / (ρ - σ)) + Real.log (ρ - σ) * Real.log (t - ρ)
noncomputable def PL1 (σ t : ℝ) : ℝ := -D00.Li2 (-(t - 1) / (1 - σ))
noncomputable def PU (ρ σ t : ℝ) : ℝ :=
  D00.Li2 ((σ - ρ) / (t - ρ)) + Real.log (t - ρ) ^ 2 / 2
noncomputable def PD (ρ t : ℝ) : ℝ := Real.log (t - ρ) ^ 2 / 2

lemma dPL {ρ σ t : ℝ} (h1 : σ < ρ) (h2 : ρ < t) :
    HasDerivAt (PL ρ σ) (Real.log (t - σ) / (t - ρ)) t := by
  have hx : -(t - ρ) / (ρ - σ) < 0 := div_neg_of_neg_of_pos (by linarith) (by linarith)
  have hd : HasDerivAt (fun x => -(x - ρ) / (ρ - σ)) (-1 / (ρ - σ)) t := by
    have := (((hasDerivAt_id t).sub_const ρ).neg).div_const (ρ - σ)
    simpa using this
  have e1 := ((D00.hasDerivAt_Li2 (by linarith) hx.ne).comp t hd).neg
  have e2 := (((hasDerivAt_id t).sub_const ρ).log (by simp; linarith)).const_mul (Real.log (ρ - σ))
  have := e1.add e2
  unfold PL
  refine this.congr_deriv ?_
  have hq : 1 - -(t - ρ) / (ρ - σ) = (t - σ) / (ρ - σ) := by field_simp; ring
  rw [hq, Real.log_div (by linarith) (by linarith)]
  simp only [id]
  field_simp
  ring

lemma dPL1 {σ t : ℝ} (h1 : σ < 1) (h2 : 1 < t) :
    HasDerivAt (PL1 σ) ((Real.log (t - σ) - Real.log (1 - σ)) / (t - 1)) t := by
  have hx : -(t - 1) / (1 - σ) < 0 := div_neg_of_neg_of_pos (by linarith) (by linarith)
  have hd : HasDerivAt (fun x => -(x - 1) / (1 - σ)) (-1 / (1 - σ)) t := by
    have := (((hasDerivAt_id t).sub_const 1).neg).div_const (1 - σ)
    simpa using this
  have e1 := ((D00.hasDerivAt_Li2 (by linarith) hx.ne).comp t hd).neg
  unfold PL1
  refine e1.congr_deriv ?_
  have hq : 1 - -(t - 1) / (1 - σ) = (t - σ) / (1 - σ) := by field_simp; ring
  rw [hq, Real.log_div (by linarith) (by linarith)]
  have : t - 1 ≠ 0 := by linarith
  have : 1 - σ ≠ 0 := by linarith
  field_simp

lemma dPU {ρ σ t : ℝ} (h1 : ρ < σ) (h2 : σ < t) :
    HasDerivAt (PU ρ σ) (Real.log (t - σ) / (t - ρ)) t := by
  have hx0 : 0 < (σ - ρ) / (t - ρ) := div_pos (by linarith) (by linarith)
  have hx1 : (σ - ρ) / (t - ρ) < 1 := by rw [div_lt_one (by linarith)]; linarith
  have htr : t - ρ ≠ 0 := by linarith
  have hd : HasDerivAt (fun x => (σ - ρ) / (x - ρ)) (-(σ - ρ) / (t - ρ) ^ 2) t := by
    have := ((hasDerivAt_id t).sub_const ρ).inv htr
    have := this.const_mul (σ - ρ)
    convert this using 1
    · ext x; simp [div_eq_mul_inv]
    · simp; ring
  have e1 := (D00.hasDerivAt_Li2 hx1 hx0.ne').comp t hd
  have e2 := ((((hasDerivAt_id t).sub_const ρ).log htr).pow 2).div_const 2
  have := e1.add e2
  unfold PU
  refine this.congr_deriv ?_
  have hq : 1 - (σ - ρ) / (t - ρ) = (t - σ) / (t - ρ) := by field_simp; ring
  rw [hq, Real.log_div (by linarith) (by linarith)]
  simp only [id]
  have : σ - ρ ≠ 0 := by linarith
  field_simp
  ring

lemma dPD {ρ t : ℝ} (h : ρ < t) :
    HasDerivAt (PD ρ) (Real.log (t - ρ) / (t - ρ)) t := by
  have htr : t - ρ ≠ 0 := by linarith
  have e2 := ((((hasDerivAt_id t).sub_const ρ).log htr).pow 2).div_const 2
  unfold PD
  refine e2.congr_deriv ?_
  simp only [id]
  field_simp
  ring

lemma dLg {ρ t : ℝ} (h : ρ < t) :
    HasDerivAt (fun x => Real.log (x - ρ)) (1 / (t - ρ)) t := by
  have := ((hasDerivAt_id t).sub_const ρ).log (by simp; linarith)
  simpa using this


noncomputable def kap (c : ℝ) : ℝ := 2 * c ^ 2 - 1
noncomputable def r1 (c : ℝ) : ℝ := 2 * c - 1
noncomputable def r2 (c : ℝ) : ℝ := -1 - 2 * c

noncomputable def H (c t : ℝ) : ℝ :=
  PL1 (kap c) t - PL1 (r1 c) t - PL1 (r2 c) t
  - 2 * PU (-1) 1 t - 2 * PD (-1) t + PU (-1) (kap c) t + PU (-1) (r1 c) t + PL (-1) (r2 c) t
  + PU (kap c) 1 t + PL (kap c) (-1) t - PU (kap c) (r1 c) t - PL (kap c) (r2 c) t
  - c * PU (r1 c) 1 t - c * PL (r1 c) (-1) t - c * PL (r1 c) (kap c) t + 2 * c * PD (r1 c) t
  + 2 * c * PL (r1 c) (r2 c) t
  + c * PU (r2 c) 1 t + c * PU (r2 c) (-1) t + c * PU (r2 c) (kap c) t
  - 2 * c * PU (r2 c) (r1 c) t - 2 * c * PD (r2 c) t
  + Real.log 2 * Real.log (t - -1) - c * Real.log 2 * Real.log (t - r1 c)
  + c * Real.log 2 * Real.log (t - r2 c)

noncomputable def cA (c t : ℝ) : ℝ :=
  -2 / (t - -1) + 1 / (t - kap c) - c / (t - r1 c) + c / (t - r2 c)
noncomputable def cB (c t : ℝ) : ℝ :=
  -1 / (t - 1) + 1 / (t - -1) - 1 / (t - kap c) + 2 * c / (t - r1 c) - 2 * c / (t - r2 c)

lemma order_facts {c : ℝ} (h0 : 0 < c) (h1 : c < 1) :
    r2 c < -1 ∧ -1 < kap c ∧ kap c < r1 c ∧ r1 c < 1 ∧ -1 < r1 c := by
  unfold r2 kap r1; refine ⟨by linarith, by nlinarith, by nlinarith, by linarith, by linarith⟩

lemma hasDerivAt_H {c t : ℝ} (h0 : 0 < c) (h1 : c < 1) (ht : 1 < t) :
    HasDerivAt (H c)
      (cA c t * (Real.log (t - 1) + Real.log (t - -1) - Real.log 2 - Real.log (t - kap c))
        + cB c t * (Real.log (t - r1 c) + Real.log (t - r2 c) - Real.log 2
          - Real.log (t - kap c))) t := by
  obtain ⟨o1, o2, o3, o4, o5⟩ := order_facts h0 h1
  have q1 := dPL1 (σ := kap c) (by linarith) ht
  have q2 := dPL1 (σ := r1 c) (by linarith) ht
  have q3 := dPL1 (σ := r2 c) (by linarith) ht
  have a1 := dPU (ρ := -1) (σ := 1) (t := t) (by norm_num) ht
  have a2 := dPD (ρ := -1) (t := t) (by linarith)
  have a3 := dPU (ρ := -1) (σ := kap c) (t := t) o2 (by linarith)
  have a4 := dPU (ρ := -1) (σ := r1 c) (t := t) o5 (by linarith)
  have a5 := dPL (ρ := -1) (σ := r2 c) (t := t) o1 (by linarith)
  have b1 := dPU (ρ := kap c) (σ := 1) (t := t) (by linarith) ht
  have b2 := dPL (ρ := kap c) (σ := -1) (t := t) o2 (by linarith)
  have b3 := dPU (ρ := kap c) (σ := r1 c) (t := t) o3 (by linarith)
  have b4 := dPL (ρ := kap c) (σ := r2 c) (t := t) (by linarith) (by linarith)
  have c1 := dPU (ρ := r1 c) (σ := 1) (t := t) o4 ht
  have c2 := dPL (ρ := r1 c) (σ := -1) (t := t) o5 (by linarith)
  have c3 := dPL (ρ := r1 c) (σ := kap c) (t := t) o3 (by linarith)
  have c4 := dPD (ρ := r1 c) (t := t) (by linarith)
  have c5 := dPL (ρ := r1 c) (σ := r2 c) (t := t) (by linarith) (by linarith)
  have d1 := dPU (ρ := r2 c) (σ := 1) (t := t) (by linarith) ht
  have d2 := dPU (ρ := r2 c) (σ := -1) (t := t) o1 (by linarith)
  have d3 := dPU (ρ := r2 c) (σ := kap c) (t := t) (by linarith) (by linarith)
  have d4 := dPU (ρ := r2 c) (σ := r1 c) (t := t) (by linarith) (by linarith)
  have d5 := dPD (ρ := r2 c) (t := t) (by linarith)
  have l1 := (dLg (ρ := -1) (t := t) (by linarith)).const_mul (Real.log 2)
  have l2 := (dLg (ρ := r1 c) (t := t) (by linarith)).const_mul (c * Real.log 2)
  have l3 := (dLg (ρ := r2 c) (t := t) (by linarith)).const_mul (c * Real.log 2)
  have hall := q1
  replace hall := hall.sub q2
  replace hall := hall.sub q3
  replace hall := hall.sub (a1.const_mul 2)
  replace hall := hall.sub (a2.const_mul 2)
  replace hall := hall.add a3
  replace hall := hall.add a4
  replace hall := hall.add a5
  replace hall := hall.add b1
  replace hall := hall.add b2
  replace hall := hall.sub b3
  replace hall := hall.sub b4
  replace hall := hall.sub (c1.const_mul c)
  replace hall := hall.sub (c2.const_mul c)
  replace hall := hall.sub (c3.const_mul c)
  replace hall := hall.add (c4.const_mul (2 * c))
  replace hall := hall.add (c5.const_mul (2 * c))
  replace hall := hall.add (d1.const_mul c)
  replace hall := hall.add (d2.const_mul c)
  replace hall := hall.add (d3.const_mul c)
  replace hall := hall.sub (d4.const_mul (2 * c))
  replace hall := hall.sub (d5.const_mul (2 * c))
  replace hall := hall.add l1
  replace hall := hall.sub l2
  replace hall := hall.add l3
  have hH : H c = fun t => H c t := rfl
  rw [hH]
  unfold H
  refine hall.congr_deriv ?_
  have e1 : Real.log (1 - kap c) = Real.log 2 + Real.log (1 - c) + Real.log (1 + c) := by
    rw [show 1 - kap c = 2 * (1 - c) * (1 + c) by unfold kap; ring,
      Real.log_mul (by positivity) (by linarith), Real.log_mul (by norm_num) (by linarith)]
  have e2 : Real.log (1 - r1 c) = Real.log 2 + Real.log (1 - c) := by
    rw [show 1 - r1 c = 2 * (1 - c) by unfold r1; ring, Real.log_mul (by norm_num) (by linarith)]
  have e3 : Real.log (1 - r2 c) = Real.log 2 + Real.log (1 + c) := by
    rw [show 1 - r2 c = 2 * (1 + c) by unfold r2; ring, Real.log_mul (by norm_num) (by linarith)]
  rw [e1, e2, e3]
  unfold cA cB
  ring


lemma alg_AB {c t : ℝ} (h0 : 0 < c) (h1 : c < 1) (ht : 1 < t) :
    let g := (t ^ 2 - 1) / (2 * (t - kap c))
    let R := (t ^ 2 - 2 * kap c * t + 1) / (2 * (t - kap c))
    cA c t * ((t - kap c) / R) = -(g - kap c) / ((1 + g) * R) + (R / (1 + g) - 1) / g ∧
    cB c t * ((t - kap c) / R) = -R / ((1 + g) * g) := by
  intro g R
  obtain ⟨o1, o2, o3, o4, o5⟩ := order_facts h0 h1
  have hk : t - kap c ≠ 0 := by linarith
  have hk' : 0 < t - kap c := by linarith
  have hR0 : 0 < t ^ 2 - 2 * kap c * t + 1 := by nlinarith
  have hg0 : 0 < t ^ 2 - 1 := by nlinarith
  have hp : t - r1 c ≠ 0 := by linarith
  have hq : t - r2 c ≠ 0 := by linarith
  have hm : t - -1 ≠ 0 := by linarith
  have h1' : t - 1 ≠ 0 := by linarith
  have hgpos : 0 < g := by positivity
  have hRpos : 0 < R := by positivity
  have h1g : 1 + g = (t - r1 c) * (t - r2 c) / (2 * (t - kap c)) := by
    simp only [g]; unfold r1 r2 kap at *; field_simp; ring
  have h1gpos : 0 < 1 + g := by linarith
  rw [h1g]
  simp only [g, R]
  unfold cA cB
  unfold r1 r2 kap at *
  have hpq : (t - (2 * c - 1)) * (t - (-1 - 2 * c)) ≠ 0 := mul_ne_zero hp hq
  have hR1 : t ^ 2 - 2 * (2 * c ^ 2 - 1) * t + 1 ≠ 0 := hR0.ne'
  have hg1 : t ^ 2 - 1 ≠ 0 := hg0.ne'
  have hR2 : t * (t - 2 * (2 * c ^ 2 - 1)) + 1 ≠ 0 := by nlinarith
  constructor
  · field_simp
    ring
  · field_simp
    ring

end D04

namespace D04
open Real

lemma H_end {α : ℝ} (hα0 : 0 < α) (hα1 : α < π) :
    H (c α) (1 + 2 * s α) - H (c α) 1 = Mtable α := by
  obtain ⟨τ, h0, h1, hc, hs⟩ := tan_param hα0 hα1
  simp only [Mtable, U, V, Lcal, liSum, uTerms, vTerms, Xm, Xp, List.map_cons, List.map_nil,
    List.sum_cons, List.sum_nil]
  generalize c α = c at hc ⊢
  generalize s α = s at hs ⊢
  simp only [H, PL, PL1, PU, PD, kap, r1, r2]
  have ht : τ ≠ 0 := h0.ne'
  have hm : (1 - τ) ≠ 0 := by linarith
  have hp : (1 + τ) ≠ 0 := by linarith
  have hq : (1 + τ ^ 2) ≠ 0 := by positivity
  have hm2 : (1 - τ ^ 2) ≠ 0 := by nlinarith
  have hc' : c = (1 - τ) * (1 + τ) / (1 + τ ^ 2) := by rw [hc]; ring
  have E0 : (1 + 2 * s - 1 : ℝ) = 4 * τ / (1 + τ ^ 2) := by simp only [hc, hs]; field_simp; ring
  have E1 : (1 - (2 * c ^ 2 - 1) : ℝ) = 8 * τ ^ 2 / (1 + τ ^ 2) ^ 2 := by simp only [hc, hs]; field_simp; ring
  have E2 : (1 - (2 * c - 1) : ℝ) = 4 * τ ^ 2 / (1 + τ ^ 2) := by simp only [hc, hs]; field_simp; ring
  have E3 : (1 - (-1 - 2 * c) : ℝ) = 4 / (1 + τ ^ 2) := by simp only [hc, hs]; field_simp; ring
  have E4 : (1 + 2 * s - -1 : ℝ) = 2 * (1 + τ) ^ 2 / (1 + τ ^ 2) := by simp only [hc, hs]; field_simp; ring
  have E5 : (2 * c ^ 2 - 1 - -1 : ℝ) = 2 * (1 - τ) ^ 2 * (1 + τ) ^ 2 / (1 + τ ^ 2) ^ 2 := by simp only [hc, hs]; field_simp; ring
  have E6 : (2 * c - 1 - -1 : ℝ) = 2 * (1 - τ) * (1 + τ) / (1 + τ ^ 2) := by simp only [hc, hs]; field_simp; ring
  have E7 : (-1 - (-1 - 2 * c) : ℝ) = 2 * (1 - τ) * (1 + τ) / (1 + τ ^ 2) := by simp only [hc, hs]; field_simp; ring
  have E8 : (1 + 2 * s - (2 * c ^ 2 - 1) : ℝ) = 4 * τ * (1 + τ) ^ 2 / (1 + τ ^ 2) ^ 2 := by simp only [hc, hs]; field_simp; ring
  have E9 : (2 * c ^ 2 - 1 - (-1 - 2 * c) : ℝ) = 4 * (1 - τ) * (1 + τ) / (1 + τ ^ 2) ^ 2 := by simp only [hc, hs]; field_simp; ring
  have E10 : (2 * c - 1 - (2 * c ^ 2 - 1) : ℝ) = 4 * τ ^ 2 * (1 - τ) * (1 + τ) / (1 + τ ^ 2) ^ 2 := by simp only [hc, hs]; field_simp; ring
  have E11 : (1 + 2 * s - (2 * c - 1) : ℝ) = 4 * τ * (1 + τ) / (1 + τ ^ 2) := by simp only [hc, hs]; field_simp; ring
  have E12 : (2 * c - 1 - (-1 - 2 * c) : ℝ) = 4 * (1 - τ) * (1 + τ) / (1 + τ ^ 2) := by simp only [hc, hs]; field_simp; ring
  have E13 : (1 + 2 * s - (-1 - 2 * c) : ℝ) = 4 * (1 + τ) / (1 + τ ^ 2) := by simp only [hc, hs]; field_simp; ring
  have E14 : (1 + s - c : ℝ) = 2 * τ * (1 + τ) / (1 + τ ^ 2) := by simp only [hc, hs]; field_simp; ring
  have E15 : (1 + s + c : ℝ) = 2 * (1 + τ) / (1 + τ ^ 2) := by simp only [hc, hs]; field_simp; ring
  have E16 : (1 - 1 / c ^ 2 : ℝ) = -(4 * τ ^ 2) / ((1 - τ) ^ 2 * (1 + τ) ^ 2) := by simp only [hc, hs]; field_simp; ring
  have E17 : (1 + s : ℝ) = (1 + τ) ^ 2 / (1 + τ ^ 2) := by simp only [hc, hs]; field_simp; ring
  have E18 : (1 - c : ℝ) = 2 * τ ^ 2 / (1 + τ ^ 2) := by simp only [hc, hs]; field_simp; ring
  have E19 : (1 + c : ℝ) = 2 / (1 + τ ^ 2) := by simp only [hc, hs]; field_simp; ring
  have Q0 : D00.Li2 (-(1 + 2 * s - 1) / (1 - (2 * c ^ 2 - 1))) = D00.Li2 (-(1 / s)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q1 : D00.Li2 (-(1 + 2 * s - 1) / (1 - (2 * c - 1))) = D00.Li2 (-(s / (1 - c))) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q2 : D00.Li2 (-(1 + 2 * s - 1) / (1 - (-1 - 2 * c))) = D00.Li2 (-(s / (1 + c))) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q3 : D00.Li2 ((1 - -1) / (1 + 2 * s - -1)) = D00.Li2 (1 / (1 + s)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q4 : D00.Li2 ((2 * c ^ 2 - 1 - -1) / (1 + 2 * s - -1)) = D00.Li2 (c ^ 2 / (1 + s)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q5 : D00.Li2 ((2 * c - 1 - -1) / (1 + 2 * s - -1)) = D00.Li2 (c / (1 + s)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q6 : D00.Li2 (-(1 + 2 * s - -1) / (-1 - (-1 - 2 * c))) = D00.Li2 (-((1 + s) / c)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q7 : D00.Li2 ((1 - (2 * c ^ 2 - 1)) / (1 + 2 * s - (2 * c ^ 2 - 1))) = D00.Li2 (s / (1 + s)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q8 : D00.Li2 (-(1 + 2 * s - (2 * c ^ 2 - 1)) / (2 * c ^ 2 - 1 - -1)) = D00.Li2 (-(s * (1 + s) / c ^ 2)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q9 : D00.Li2 ((2 * c - 1 - (2 * c ^ 2 - 1)) / (1 + 2 * s - (2 * c ^ 2 - 1))) = D00.Li2 (c * (1 - c) / (s * (1 + s))) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q10 : D00.Li2 (-(1 + 2 * s - (2 * c ^ 2 - 1)) / (2 * c ^ 2 - 1 - (-1 - 2 * c))) = D00.Li2 (-(s * (1 + s) / (c * (1 + c)))) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q11 : D00.Li2 ((1 - (2 * c - 1)) / (1 + 2 * s - (2 * c - 1))) = D00.Li2 ((1 - c) / (1 + s - c)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q12 : D00.Li2 (-(1 + 2 * s - (2 * c - 1)) / (2 * c - 1 - -1)) = D00.Li2 (-((1 + s - c) / c)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q13 : D00.Li2 (-(1 + 2 * s - (2 * c - 1)) / (2 * c - 1 - (2 * c ^ 2 - 1))) = D00.Li2 (-((1 + s - c) / (c * (1 - c)))) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q14 : D00.Li2 (-(1 + 2 * s - (2 * c - 1)) / (2 * c - 1 - (-1 - 2 * c))) = D00.Li2 (-((1 + s - c) / (2 * c))) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q15 : D00.Li2 ((1 - (-1 - 2 * c)) / (1 + 2 * s - (-1 - 2 * c))) = D00.Li2 ((1 + c) / (1 + s + c)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q16 : D00.Li2 ((-1 - (-1 - 2 * c)) / (1 + 2 * s - (-1 - 2 * c))) = D00.Li2 (c / (1 + s + c)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q17 : D00.Li2 ((2 * c ^ 2 - 1 - (-1 - 2 * c)) / (1 + 2 * s - (-1 - 2 * c))) = D00.Li2 (c * (1 + c) / (1 + s + c)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q18 : D00.Li2 ((2 * c - 1 - (-1 - 2 * c)) / (1 + 2 * s - (-1 - 2 * c))) = D00.Li2 (2 * c / (1 + s + c)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q19 : D00.Li2 (-(1 - 1) / (1 - (2 * c ^ 2 - 1))) = 0 := by norm_num
  have Q20 : D00.Li2 (-(1 - 1) / (1 - (2 * c - 1))) = 0 := by norm_num
  have Q21 : D00.Li2 (-(1 - 1) / (1 - (-1 - 2 * c))) = 0 := by norm_num
  have Q22 : D00.Li2 ((1 - -1) / (1 - -1)) = π ^ 2 / 6 := by
    rw [show ((1 - -1) / (1 - -1)) = (1 : ℝ) by (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp]; exact D00.Li2_one
  have Q23 : D00.Li2 ((2 * c ^ 2 - 1 - -1) / (1 - -1)) = D00.Li2 (c ^ 2) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q24 : D00.Li2 ((2 * c - 1 - -1) / (1 - -1)) = D00.Li2 c := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q25 : D00.Li2 (-(1 - -1) / (-1 - (-1 - 2 * c))) = D00.Li2 (-(1 / c)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q26 : D00.Li2 ((1 - (2 * c ^ 2 - 1)) / (1 - (2 * c ^ 2 - 1))) = π ^ 2 / 6 := by
    rw [show ((1 - (2 * c ^ 2 - 1)) / (1 - (2 * c ^ 2 - 1))) = (1 : ℝ) by (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp]; exact D00.Li2_one
  have Q27 : D00.Li2 (-(1 - (2 * c ^ 2 - 1)) / (2 * c ^ 2 - 1 - -1)) = D00.Li2 (1 - 1 / c ^ 2) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q28 : D00.Li2 ((2 * c - 1 - (2 * c ^ 2 - 1)) / (1 - (2 * c ^ 2 - 1))) = D00.Li2 (c / (1 + c)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q29 : D00.Li2 (-(1 - (2 * c ^ 2 - 1)) / (2 * c ^ 2 - 1 - (-1 - 2 * c))) = D00.Li2 (-((1 - c) / c)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q30 : D00.Li2 ((1 - (2 * c - 1)) / (1 - (2 * c - 1))) = π ^ 2 / 6 := by
    rw [show ((1 - (2 * c - 1)) / (1 - (2 * c - 1))) = (1 : ℝ) by (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp]; exact D00.Li2_one
  have Q31 : D00.Li2 (-(1 - (2 * c - 1)) / (2 * c - 1 - -1)) = D00.Li2 (-((1 - c) / c)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q32 : D00.Li2 (-(1 - (2 * c - 1)) / (2 * c - 1 - (2 * c ^ 2 - 1))) = D00.Li2 (-(1 / c)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q33 : D00.Li2 (-(1 - (2 * c - 1)) / (2 * c - 1 - (-1 - 2 * c))) = D00.Li2 (-((1 - c) / (2 * c))) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q34 : D00.Li2 ((1 - (-1 - 2 * c)) / (1 - (-1 - 2 * c))) = π ^ 2 / 6 := by
    rw [show ((1 - (-1 - 2 * c)) / (1 - (-1 - 2 * c))) = (1 : ℝ) by (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp]; exact D00.Li2_one
  have Q35 : D00.Li2 ((-1 - (-1 - 2 * c)) / (1 - (-1 - 2 * c))) = D00.Li2 (c / (1 + c)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q36 : D00.Li2 ((2 * c ^ 2 - 1 - (-1 - 2 * c)) / (1 - (-1 - 2 * c))) = D00.Li2 c := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have Q37 : D00.Li2 ((2 * c - 1 - (-1 - 2 * c)) / (1 - (-1 - 2 * c))) = D00.Li2 (2 * c / (1 + c)) := by
    congr 1; (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); field_simp <;> ring
  have G0 : Real.log (1 + 2 * s - -1) = 1 * Real.log 2 + 0 * Real.log τ + 0 * Real.log (1 - τ) + 2 * Real.log (1 + τ) + -1 * Real.log (1 + τ ^ 2) := by
    rw [show (1 + 2 * s - -1) = monoz ((2 : ℝ) ^ (1 : ℤ)) (0) (0) (2) (-1) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G1 : Real.log (-1 - (-1 - 2 * c)) = 1 * Real.log 2 + 0 * Real.log τ + 1 * Real.log (1 - τ) + 1 * Real.log (1 + τ) + -1 * Real.log (1 + τ ^ 2) := by
    rw [show (-1 - (-1 - 2 * c)) = monoz ((2 : ℝ) ^ (1 : ℤ)) (0) (1) (1) (-1) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G2 : Real.log (1 + 2 * s - (2 * c ^ 2 - 1)) = 2 * Real.log 2 + 1 * Real.log τ + 0 * Real.log (1 - τ) + 2 * Real.log (1 + τ) + -2 * Real.log (1 + τ ^ 2) := by
    rw [show (1 + 2 * s - (2 * c ^ 2 - 1)) = monoz ((2 : ℝ) ^ (2 : ℤ)) (1) (0) (2) (-2) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G3 : Real.log (2 * c ^ 2 - 1 - -1) = 1 * Real.log 2 + 0 * Real.log τ + 2 * Real.log (1 - τ) + 2 * Real.log (1 + τ) + -2 * Real.log (1 + τ ^ 2) := by
    rw [show (2 * c ^ 2 - 1 - -1) = monoz ((2 : ℝ) ^ (1 : ℤ)) (0) (2) (2) (-2) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G4 : Real.log (2 * c ^ 2 - 1 - (-1 - 2 * c)) = 2 * Real.log 2 + 0 * Real.log τ + 1 * Real.log (1 - τ) + 1 * Real.log (1 + τ) + -2 * Real.log (1 + τ ^ 2) := by
    rw [show (2 * c ^ 2 - 1 - (-1 - 2 * c)) = monoz ((2 : ℝ) ^ (2 : ℤ)) (0) (1) (1) (-2) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G5 : Real.log (1 + 2 * s - (2 * c - 1)) = 2 * Real.log 2 + 1 * Real.log τ + 0 * Real.log (1 - τ) + 1 * Real.log (1 + τ) + -1 * Real.log (1 + τ ^ 2) := by
    rw [show (1 + 2 * s - (2 * c - 1)) = monoz ((2 : ℝ) ^ (2 : ℤ)) (1) (0) (1) (-1) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G6 : Real.log (2 * c - 1 - -1) = 1 * Real.log 2 + 0 * Real.log τ + 1 * Real.log (1 - τ) + 1 * Real.log (1 + τ) + -1 * Real.log (1 + τ ^ 2) := by
    rw [show (2 * c - 1 - -1) = monoz ((2 : ℝ) ^ (1 : ℤ)) (0) (1) (1) (-1) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G7 : Real.log (2 * c - 1 - (2 * c ^ 2 - 1)) = 2 * Real.log 2 + 2 * Real.log τ + 1 * Real.log (1 - τ) + 1 * Real.log (1 + τ) + -2 * Real.log (1 + τ ^ 2) := by
    rw [show (2 * c - 1 - (2 * c ^ 2 - 1)) = monoz ((2 : ℝ) ^ (2 : ℤ)) (2) (1) (1) (-2) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G8 : Real.log (2 * c - 1 - (-1 - 2 * c)) = 2 * Real.log 2 + 0 * Real.log τ + 1 * Real.log (1 - τ) + 1 * Real.log (1 + τ) + -1 * Real.log (1 + τ ^ 2) := by
    rw [show (2 * c - 1 - (-1 - 2 * c)) = monoz ((2 : ℝ) ^ (2 : ℤ)) (0) (1) (1) (-1) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G9 : Real.log (1 + 2 * s - (-1 - 2 * c)) = 2 * Real.log 2 + 0 * Real.log τ + 0 * Real.log (1 - τ) + 1 * Real.log (1 + τ) + -1 * Real.log (1 + τ ^ 2) := by
    rw [show (1 + 2 * s - (-1 - 2 * c)) = monoz ((2 : ℝ) ^ (2 : ℤ)) (0) (0) (1) (-1) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G11 : Real.log (1 - -1) = 1 * Real.log 2 + 0 * Real.log τ + 0 * Real.log (1 - τ) + 0 * Real.log (1 + τ) + 0 * Real.log (1 + τ ^ 2) := by
    rw [show (1 - -1) = monoz ((2 : ℝ) ^ (1 : ℤ)) (0) (0) (0) (0) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G12 : Real.log (1 - (2 * c ^ 2 - 1)) = 3 * Real.log 2 + 2 * Real.log τ + 0 * Real.log (1 - τ) + 0 * Real.log (1 + τ) + -2 * Real.log (1 + τ ^ 2) := by
    rw [show (1 - (2 * c ^ 2 - 1)) = monoz ((2 : ℝ) ^ (3 : ℤ)) (2) (0) (0) (-2) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G13 : Real.log (1 - (2 * c - 1)) = 2 * Real.log 2 + 2 * Real.log τ + 0 * Real.log (1 - τ) + 0 * Real.log (1 + τ) + -1 * Real.log (1 + τ ^ 2) := by
    rw [show (1 - (2 * c - 1)) = monoz ((2 : ℝ) ^ (2 : ℤ)) (2) (0) (0) (-1) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G14 : Real.log (1 - (-1 - 2 * c)) = 2 * Real.log 2 + 0 * Real.log τ + 0 * Real.log (1 - τ) + 0 * Real.log (1 + τ) + -1 * Real.log (1 + τ ^ 2) := by
    rw [show (1 - (-1 - 2 * c)) = monoz ((2 : ℝ) ^ (2 : ℤ)) (0) (0) (0) (-1) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G15 : Real.log (1 + s) = 0 * Real.log 2 + 0 * Real.log τ + 0 * Real.log (1 - τ) + 2 * Real.log (1 + τ) + -1 * Real.log (1 + τ ^ 2) := by
    rw [show (1 + s) = monoz ((2 : ℝ) ^ (0 : ℤ)) (0) (0) (2) (-1) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G16 : Real.log (c / (1 + s)) = 0 * Real.log 2 + 0 * Real.log τ + 1 * Real.log (1 - τ) + -1 * Real.log (1 + τ) + 0 * Real.log (1 + τ ^ 2) := by
    rw [show (c / (1 + s)) = monoz ((2 : ℝ) ^ (0 : ℤ)) (0) (1) (-1) (0) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G17 : Real.log (c / (1 + c)) = -1 * Real.log 2 + 0 * Real.log τ + 1 * Real.log (1 - τ) + 1 * Real.log (1 + τ) + 0 * Real.log (1 + τ ^ 2) := by
    rw [show (c / (1 + c)) = monoz ((2 : ℝ) ^ (-1 : ℤ)) (0) (1) (1) (0) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G18 : Real.log ((1 + s) / s) = -1 * Real.log 2 + -1 * Real.log τ + 0 * Real.log (1 - τ) + 2 * Real.log (1 + τ) + 0 * Real.log (1 + τ ^ 2) := by
    rw [show ((1 + s) / s) = monoz ((2 : ℝ) ^ (-1 : ℤ)) (-1) (0) (2) (0) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G19 : Real.log 4 = 2 * Real.log 2 + 0 * Real.log τ + 0 * Real.log (1 - τ) + 0 * Real.log (1 + τ) + 0 * Real.log (1 + τ ^ 2) := by
    rw [show 4 = monoz ((2 : ℝ) ^ (2 : ℤ)) (0) (0) (0) (0) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G20 : Real.log ((1 + s - c) / (1 - c)) = 0 * Real.log 2 + -1 * Real.log τ + 0 * Real.log (1 - τ) + 1 * Real.log (1 + τ) + 0 * Real.log (1 + τ ^ 2) := by
    rw [show ((1 + s - c) / (1 - c)) = monoz ((2 : ℝ) ^ (0 : ℤ)) (-1) (0) (1) (0) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G21 : Real.log ((1 + s + c) / (1 + c)) = 0 * Real.log 2 + 0 * Real.log τ + 0 * Real.log (1 - τ) + 1 * Real.log (1 + τ) + 0 * Real.log (1 + τ ^ 2) := by
    rw [show ((1 + s + c) / (1 + c)) = monoz ((2 : ℝ) ^ (0 : ℤ)) (0) (0) (1) (0) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  have G22 : Real.log ((1 + c) * (1 + s + c)) = 2 * Real.log 2 + 0 * Real.log τ + 0 * Real.log (1 - τ) + 1 * Real.log (1 + τ) + -2 * Real.log (1 + τ ^ 2) := by
    rw [show ((1 + c) * (1 + s + c)) = monoz ((2 : ℝ) ^ (2 : ℤ)) (0) (0) (1) (-2) τ by
      (try rw [E0]); (try rw [E1]); (try rw [E2]); (try rw [E3]); (try rw [E4]); (try rw [E5]); (try rw [E6]); (try rw [E7]); (try rw [E8]); (try rw [E9]); (try rw [E10]); (try rw [E11]); (try rw [E12]); (try rw [E13]); (try rw [E14]); (try rw [E15]); (try rw [E16]); (try rw [E17]); (try rw [E18]); (try rw [E19]); (try rw [hc']); (try rw [hs]); simp only [monoz, zpow_neg, zpow_ofNat, zpow_zero, pow_zero, zpow_one, pow_one, mul_one, one_mul, zpow_natCast]; field_simp <;> ring,
      log_monoz (by positivity) _ _ _ _ h0 h1, Real.log_zpow]; push_cast; ring
  rw [Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12, Q13, Q14, Q15, Q16, Q17, Q18, Q19, Q20, Q21, Q22, Q23, Q24, Q25, Q26, Q27, Q28, Q29, Q30, Q31, Q32, Q33, Q34, Q35, Q36, Q37, G0, G1, G2, G3, G4, G5, G6, G7, G8, G9, G11, G12, G13, G14, G15, G16, G17, G18, G19, G20, G21, G22]
  ring

end D04
namespace D04
open Real

noncomputable def Fi (α u g : ℝ) : ℝ :=
  ((1 + g) ^ (-1 - u) * Real.sqrt (1 + g ^ (2 - 2 * u) - 2 * Real.cos α * g ^ (1 - u)) - 1)
    * g ^ (u - 1)

noncomputable def Fd (α u g : ℝ) : ℝ :=
  -Real.log (1 + g) * (1 + g) ^ (-1 - u)
      * Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u)) / g ^ (1 - u)
    - Real.log g * (1 + g) ^ (-1 - u) * (g ^ (1 - u) - Real.cos α)
      / Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u))
    + ((1 + g) ^ (-1 - u) * Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u)) - 1)
      * Real.log g / g ^ (1 - u)

lemma Fi_eq {α u g : ℝ} (hg : 0 < g) :
    Fi α u g = ((1 + g) ^ (-1 - u)
      * Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u)) - 1) / g ^ (1 - u) := by
  unfold Fi
  have h1 : g ^ (2 - 2 * u) = (g ^ (1 - u)) ^ 2 := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hg.le]; congr 1; push_cast; ring
  have h2 : g ^ (u - 1) = (g ^ (1 - u))⁻¹ := by
    rw [← Real.rpow_neg hg.le]; congr 1; ring
  rw [h1, h2, div_eq_mul_inv]

lemma rad_pos {α x : ℝ} (hα0 : 0 < α) (hα1 : α < π) :
    0 < 1 + x ^ 2 - 2 * Real.cos α * x := by
  have hs : 0 < Real.sin α := Real.sin_pos_of_pos_of_lt_pi hα0 hα1
  have := Real.sin_sq_add_cos_sq α
  nlinarith [sq_nonneg (x - Real.cos α)]

lemma hasDerivAt_Fi {α u g : ℝ} (hα0 : 0 < α) (hα1 : α < π) (hg : 0 < g) :
    HasDerivAt (fun u => Fi α u g) (Fd α u g) u := by
  have hfun : (fun u => Fi α u g) = fun u => ((1 + g) ^ (-1 - u)
      * Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u)) - 1) / g ^ (1 - u) := by
    funext u; exact Fi_eq hg
  rw [hfun]
  have hh : HasDerivAt (fun u : ℝ => g ^ (1 - u)) (Real.log g * (-1) * g ^ (1 - u)) u :=
    ((hasDerivAt_id u).const_sub 1).const_rpow hg |>.congr_deriv (by simp)
  have hp : HasDerivAt (fun u : ℝ => (1 + g) ^ (-1 - u))
      (Real.log (1 + g) * (-1) * (1 + g) ^ (-1 - u)) u :=
    ((hasDerivAt_id u).const_sub (-1)).const_rpow (by linarith) |>.congr_deriv (by simp)
  have hr := ((hh.pow 2).const_add 1).sub (hh.const_mul (2 * Real.cos α))
  have hpos := rad_pos (x := g ^ (1 - u)) hα0 hα1
  have hsq := hr.sqrt (by exact hpos.ne')
  have hgu : 0 < g ^ (1 - u) := Real.rpow_pos_of_pos hg _
  have := ((hp.mul hsq).sub_const 1).div hh hgu.ne'
  refine this.congr_deriv ?_
  unfold Fd
  have hsr : 0 < Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u)) :=
    Real.sqrt_pos.mpr hpos
  simp only [Pi.add_apply, Pi.sub_apply, Pi.mul_apply, Pi.pow_apply, Nat.cast_ofNat]
  field_simp
  ring



lemma Fd_bound {α u g : ℝ} (hα0 : 0 < α) (hα1 : α < π) (hg0 : 0 < g) (hg1 : g ≤ 1)
    (hu0 : -(1 / 4) < u) (hu1 : u < 1 / 4) :
    ‖Fd α u g‖ ≤ 18 * g ^ (-(1 / 2) : ℝ) := by
  set κ := Real.cos α with hκ
  have hκ1 : κ ≤ 1 := Real.cos_le_one α
  have hκ2 : -1 ≤ κ := Real.neg_one_le_cos α
  set h := g ^ (1 - u) with hh
  set p := (1 + g) ^ (-1 - u) with hp
  set w := g ^ (-(1 / 4) : ℝ) with hw
  have hh0 : 0 < h := Real.rpow_pos_of_pos hg0 _
  have hh1 : h ≤ 1 := Real.rpow_le_one hg0.le hg1 (by linarith)
  have hw1 : 1 ≤ w := Real.one_le_rpow_of_pos_of_le_one_of_nonpos hg0 hg1 (by norm_num)
  have hgw : g ≤ h * w := by
    have e : g = h * g ^ u := by
      rw [hh, ← Real.rpow_add hg0]; norm_num
    have : g ^ u ≤ w := Real.rpow_le_rpow_of_exponent_ge hg0 hg1 hu0.le
    rw [e]; exact mul_le_mul_of_nonneg_left this hh0.le
  have hww : w * w = g ^ (-(1 / 2) : ℝ) := by
    rw [hw, ← Real.rpow_add hg0]; norm_num
  have hlog : |Real.log g| ≤ 4 * w := by
    have := Real.abs_log_mul_self_rpow_lt g (1 / 4) hg0 hg1 (by norm_num)
    have e : g ^ (1 / 4 : ℝ) * w = 1 := by
      rw [hw, ← Real.rpow_add hg0]; norm_num
    have hq : 0 < g ^ (1 / 4 : ℝ) := Real.rpow_pos_of_pos hg0 _
    rw [abs_mul, abs_of_pos hq] at this
    nlinarith [abs_nonneg (Real.log g)]
  have hp0 : 0 < p := Real.rpow_pos_of_pos (by linarith) _
  have hp1 : p ≤ 1 := Real.rpow_le_one_of_one_le_of_nonpos (by linarith) (by linarith)
  have hp2 : 1 - 2 * g ≤ p := by
    have : (1 + g) ^ (-2 : ℝ) ≤ p := Real.rpow_le_rpow_of_exponent_le (by linarith) (by linarith)
    have e : (1 + g) ^ (-2 : ℝ) = 1 / (1 + g) ^ 2 := by
      rw [Real.rpow_neg (by linarith), one_div]; norm_num
    rw [e] at this
    have : 1 - 2 * g ≤ 1 / (1 + g) ^ 2 := by
      rw [le_div_iff₀ (by positivity)]; nlinarith [sq_nonneg g]
    linarith
  set r := 1 + h ^ 2 - 2 * κ * h with hr
  have hr0 : 0 < r := rad_pos hα0 hα1
  set S := Real.sqrt r with hS
  have hS0 : 0 < S := Real.sqrt_pos.mpr hr0
  have hS1 : S ≤ 1 + h := by
    rw [hS, show 1 + h = Real.sqrt ((1 + h) ^ 2) by rw [Real.sqrt_sq (by linarith)]]
    exact Real.sqrt_le_sqrt (by nlinarith)
  have hS2 : 1 - h ≤ S := by
    rw [hS, show 1 - h = Real.sqrt ((1 - h) ^ 2) by rw [Real.sqrt_sq (by linarith)]]
    exact Real.sqrt_le_sqrt (by nlinarith)
  have hS3 : |h - κ| ≤ S := Real.abs_le_sqrt (by
    have := Real.sin_sq_add_cos_sq α
    nlinarith [sq_nonneg (Real.sin α)])
  have hL0 : 0 ≤ Real.log (1 + g) := Real.log_nonneg (by linarith)
  have hL1 : Real.log (1 + g) ≤ g := by
    have := Real.log_le_sub_one_of_pos (show 0 < 1 + g by linarith); linarith
  have hww1 : w ≤ w * w := by nlinarith
  have T1 : |-Real.log (1 + g) * p * S / h| ≤ 2 * w := by
    have n0 : 0 ≤ Real.log (1 + g) * p := mul_nonneg hL0 hp0.le
    rw [abs_div, abs_of_pos hh0, div_le_iff₀ hh0, abs_of_nonpos (by
      have : 0 ≤ Real.log (1 + g) * p * S := mul_nonneg n0 hS0.le
      linarith)]
    have a1 : Real.log (1 + g) * p ≤ g :=
      (mul_le_of_le_one_right hL0 hp1).trans hL1
    have a2 : Real.log (1 + g) * p * S ≤ Real.log (1 + g) * p * 2 :=
      mul_le_mul_of_nonneg_left (by linarith) n0
    have : -Real.log (1 + g) * p * S = -(Real.log (1 + g) * p * S) := by ring
    rw [this, neg_neg]
    linarith
  have T2 : |Real.log g * p * (h - κ) / S| ≤ 4 * w := by
    rw [abs_div, abs_of_pos hS0, div_le_iff₀ hS0, abs_mul, abs_mul, abs_of_pos hp0]
    have b1 : |Real.log g| * p ≤ 4 * w :=
      (mul_le_of_le_one_right (abs_nonneg _) hp1).trans hlog
    have b0 : 0 ≤ |Real.log g| * p := mul_nonneg (abs_nonneg _) hp0.le
    have b2 : |Real.log g| * p * |h - κ| ≤ |Real.log g| * p * S :=
      mul_le_mul_of_nonneg_left hS3 b0
    have b3 : |Real.log g| * p * S ≤ 4 * w * S := mul_le_mul_of_nonneg_right b1 hS0.le
    linarith
  have hpS : |p * S - 1| ≤ h + 2 * g := by
    rw [abs_le]; constructor
    · have b1 : p * (1 - h) ≤ p * S := mul_le_mul_of_nonneg_left hS2 hp0.le
      have b2 : (1 - 2 * g) * (1 - h) ≤ p * (1 - h) :=
        mul_le_mul_of_nonneg_right hp2 (by linarith)
      have b3 : (1 - 2 * g) * (1 - h) = 1 - h - 2 * g + 2 * (g * h) := by ring
      have b4 : 0 ≤ g * h := mul_nonneg hg0.le hh0.le
      linarith
    · have : p * S ≤ 1 * (1 + h) := mul_le_mul hp1 hS1 hS0.le (by norm_num)
      linarith
  have T3 : |(p * S - 1) * Real.log g / h| ≤ 12 * (w * w) := by
    rw [abs_div, abs_of_pos hh0, div_le_iff₀ hh0, abs_mul]
    have c1 : |p * S - 1| ≤ h * (1 + 2 * w) := by
      have : h * (1 + 2 * w) = h + 2 * (h * w) := by ring
      linarith
    have c2 : |p * S - 1| * |Real.log g| ≤ h * (1 + 2 * w) * (4 * w) :=
      mul_le_mul c1 hlog (abs_nonneg _) (by positivity)
    have c3 : h * (1 + 2 * w) * (4 * w) = 4 * (h * w) + 8 * (h * (w * w)) := by ring
    have c4 : h * w ≤ h * (w * w) := mul_le_mul_of_nonneg_left hww1 hh0.le
    have c5 : 12 * (w * w) * h = 12 * (h * (w * w)) := by ring
    linarith
  unfold Fd
  rw [← hκ, ← hh, ← hp, ← hr, ← hS, Real.norm_eq_abs, ← hww]
  calc |-Real.log (1 + g) * p * S / h - Real.log g * p * (h - κ) / S
        + (p * S - 1) * Real.log g / h|
      ≤ |-Real.log (1 + g) * p * S / h| + |Real.log g * p * (h - κ) / S|
        + |(p * S - 1) * Real.log g / h| := by
        refine (abs_add_le _ _).trans ?_
        gcongr
        exact abs_sub _ _
    _ ≤ 2 * w + 4 * w + 12 * (w * w) := by gcongr
    _ ≤ 18 * (w * w) := by linarith

end D04

namespace D04
open Real

lemma alg_AB' {c t g R : ℝ} (h0 : 0 < c) (h1 : c < 1) (ht : 1 < t)
    (hg : g = (t ^ 2 - 1) / (2 * (t - kap c)))
    (hR : R = (t ^ 2 - 2 * kap c * t + 1) / (2 * (t - kap c))) :
    cA c t * ((t - kap c) / R) = -(g - kap c) / ((1 + g) * R) + (R / (1 + g) - 1) / g ∧
    cB c t * ((t - kap c) / R) = -R / ((1 + g) * g) := by
  subst hg hR; exact alg_AB h0 h1 ht

lemma cPL1 (σ : ℝ) : Continuous (PL1 σ) := by
  unfold PL1; exact (D00.continuous_Li2.comp (by fun_prop)).neg

lemma cPL {ρ σ t : ℝ} (h : ρ < t) : ContinuousAt (PL ρ σ) t := by
  unfold PL
  have hl : ContinuousAt (fun x => Real.log (x - ρ)) t :=
    (continuousAt_id.sub continuousAt_const).log (by simp; linarith)
  exact (D00.continuous_Li2.continuousAt.comp (f := fun x => -(x - ρ) / (ρ - σ))
    (by fun_prop)).neg.add (continuousAt_const.mul hl)

lemma cPU {ρ σ t : ℝ} (h : ρ < t) : ContinuousAt (PU ρ σ) t := by
  unfold PU
  have hl : ContinuousAt (fun x => Real.log (x - ρ)) t :=
    (continuousAt_id.sub continuousAt_const).log (by simp; linarith)
  have hd : ContinuousAt (fun x => (σ - ρ) / (x - ρ)) t :=
    continuousAt_const.div (continuousAt_id.sub continuousAt_const) (by show t - ρ ≠ 0; linarith)
  exact (D00.continuous_Li2.continuousAt.comp (f := fun x => (σ - ρ) / (x - ρ)) hd).add
    ((hl.pow 2).div_const 2)

lemma cPD {ρ t : ℝ} (h : ρ < t) : ContinuousAt (PD ρ) t := by
  unfold PD
  have hl : ContinuousAt (fun x => Real.log (x - ρ)) t :=
    (continuousAt_id.sub continuousAt_const).log (by simp; linarith)
  exact (hl.pow 2).div_const 2

lemma cLg {ρ t : ℝ} (h : ρ < t) : ContinuousAt (fun x => Real.log (x - ρ)) t :=
  (continuousAt_id.sub continuousAt_const).log (by simp; linarith)

lemma continuousAt_H {c t : ℝ} (h0 : 0 < c) (h1 : c < 1) (ht : 1 ≤ t) :
    ContinuousAt (H c) t := by
  obtain ⟨o1, o2, o3, o4, o5⟩ := order_facts h0 h1
  have q1 := (cPL1 (kap c)).continuousAt (x := t)
  have q2 := (cPL1 (r1 c)).continuousAt (x := t)
  have q3 := (cPL1 (r2 c)).continuousAt (x := t)
  have a1 := cPU (ρ := -1) (σ := 1) (t := t) (by linarith)
  have a2 := cPD (ρ := -1) (t := t) (by linarith)
  have a3 := cPU (ρ := -1) (σ := kap c) (t := t) (by linarith)
  have a4 := cPU (ρ := -1) (σ := r1 c) (t := t) (by linarith)
  have a5 := cPL (ρ := -1) (σ := r2 c) (t := t) (by linarith)
  have b1 := cPU (ρ := kap c) (σ := 1) (t := t) (by linarith)
  have b2 := cPL (ρ := kap c) (σ := -1) (t := t) (by linarith)
  have b3 := cPU (ρ := kap c) (σ := r1 c) (t := t) (by linarith)
  have b4 := cPL (ρ := kap c) (σ := r2 c) (t := t) (by linarith)
  have c1 := cPU (ρ := r1 c) (σ := 1) (t := t) (by linarith)
  have c2 := cPL (ρ := r1 c) (σ := -1) (t := t) (by linarith)
  have c3 := cPL (ρ := r1 c) (σ := kap c) (t := t) (by linarith)
  have c4 := cPD (ρ := r1 c) (t := t) (by linarith)
  have c5 := cPL (ρ := r1 c) (σ := r2 c) (t := t) (by linarith)
  have d1 := cPU (ρ := r2 c) (σ := 1) (t := t) (by linarith)
  have d2 := cPU (ρ := r2 c) (σ := -1) (t := t) (by linarith)
  have d3 := cPU (ρ := r2 c) (σ := kap c) (t := t) (by linarith)
  have d4 := cPU (ρ := r2 c) (σ := r1 c) (t := t) (by linarith)
  have d5 := cPD (ρ := r2 c) (t := t) (by linarith)
  have l1 := cLg (ρ := -1) (t := t) (by linarith)
  have l2 := cLg (ρ := r1 c) (t := t) (by linarith)
  have l3 := cLg (ρ := r2 c) (t := t) (by linarith)
  have hH : H c = fun t => H c t := rfl
  rw [hH]
  unfold H
  fun_prop


lemma Fi0_bound {α g : ℝ} (hα0 : 0 < α) (hα1 : α < π) (hg0 : 0 < g) (hg1 : g ≤ 1) :
    ‖Fi α 0 g‖ ≤ 3 := by
  rw [Fi_eq hg0]
  simp only [sub_zero, Real.rpow_one, Real.rpow_neg_one]
  have hr0 := rad_pos (x := g) hα0 hα1
  have hκ1 : Real.cos α ≤ 1 := Real.cos_le_one α
  set S := Real.sqrt (1 + g ^ 2 - 2 * Real.cos α * g) with hS
  have hS0 : 0 < S := Real.sqrt_pos.mpr hr0
  have hS1 : S ≤ 1 + g := by
    rw [hS, show 1 + g = Real.sqrt ((1 + g) ^ 2) by rw [Real.sqrt_sq (by linarith)]]
    exact Real.sqrt_le_sqrt (by nlinarith [Real.neg_one_le_cos α])
  have hS2 : 1 - g ≤ S := by
    rw [hS, show 1 - g = Real.sqrt ((1 - g) ^ 2) by rw [Real.sqrt_sq (by linarith)]]
    exact Real.sqrt_le_sqrt (by nlinarith)
  have hg1' : 0 < 1 + g := by linarith
  have e1 : (1 + g)⁻¹ * S ≤ 1 := by
    rw [inv_mul_le_iff₀ hg1']; linarith
  have e2 : 1 - 2 * g ≤ (1 + g)⁻¹ * S := by
    rw [le_inv_mul_iff₀ hg1']; nlinarith
  rw [Real.norm_eq_abs, abs_div, abs_of_pos hg0, div_le_iff₀ hg0, abs_le]
  constructor <;> nlinarith

lemma hasDerivAt_Mint {α : ℝ} (hα0 : 0 < α) (hα1 : α < π) :
    IntervalIntegrable (Fd α 0) MeasureTheory.volume 0 1 ∧
      HasDerivAt (fun u => Mint u α) (∫ g in (0 : ℝ)..1, Fd α 0 g) 0 := by
  have hmeas : ∀ u, MeasureTheory.AEStronglyMeasurable (Fi α u)
      (MeasureTheory.volume.restrict (Set.uIoc (0 : ℝ) 1)) := fun u =>
    (by unfold Fi; fun_prop : Measurable (Fi α u)).aestronglyMeasurable
  have hmeas' : MeasureTheory.AEStronglyMeasurable (Fd α 0)
      (MeasureTheory.volume.restrict (Set.uIoc (0 : ℝ) 1)) :=
    (by unfold Fd; fun_prop : Measurable (Fd α 0)).aestronglyMeasurable
  have hint0 : IntervalIntegrable (Fi α 0) MeasureTheory.volume 0 1 := by
    refine (intervalIntegrable_const (c := (3 : ℝ))).mono_fun' (hmeas 0) ?_
    refine (MeasureTheory.ae_restrict_iff' measurableSet_uIoc).mpr
      (Filter.Eventually.of_forall fun g hg => ?_)
    rw [Set.uIoc_of_le zero_le_one] at hg
    exact Fi0_bound hα0 hα1 hg.1 hg.2
  have key := intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (μ := MeasureTheory.volume) (a := 0) (b := 1)
    (F := fun u g => Fi α u g) (F' := fun u g => Fd α u g) (x₀ := 0)
    (bound := fun g => 18 * g ^ (-(1 / 2) : ℝ)) (s := Metric.ball (0 : ℝ) (1 / 4))
    (Metric.ball_mem_nhds _ (by norm_num))
    (Filter.Eventually.of_forall (fun u => hmeas u)) hint0 hmeas'
    (Filter.Eventually.of_forall fun g hg u hu => by
      rw [Set.uIoc_of_le zero_le_one] at hg
      rw [Metric.mem_ball, Real.dist_eq, sub_zero, abs_lt] at hu
      exact Fd_bound hα0 hα1 hg.1 hg.2 hu.1 hu.2)
    ((intervalIntegral.intervalIntegrable_rpow' (by norm_num)).const_mul 18)
    (Filter.Eventually.of_forall fun g hg u _ => by
      rw [Set.uIoc_of_le zero_le_one] at hg
      exact hasDerivAt_Fi hα0 hα1 hg.1)
  exact key


lemma cos_eq_kap {α : ℝ} : Real.cos α = kap (c α) := by
  unfold kap c
  rw [← Real.cos_two_mul]; ring_nf

lemma integral_Fd0 {α : ℝ} (hα0 : 0 < α) (hα1 : α < π)
    (hint : IntervalIntegrable (Fd α 0) MeasureTheory.volume 0 1) :
    ∫ g in (0 : ℝ)..1, Fd α 0 g = H (c α) (1 + 2 * s α) - H (c α) 1 := by
  obtain ⟨hc0, hc1, hs0, hs1, hcs⟩ := cs_facts hα0 hα1
  have hκ : Real.cos α = kap (c α) := cos_eq_kap
  have hκ1 : Real.cos α < 1 := by
    rw [hκ]; unfold kap; nlinarith
  have hκ2 : -1 ≤ Real.cos α := Real.neg_one_le_cos α
  have hT1 : ∀ g ∈ Set.Icc (0 : ℝ) 1, 1 ≤ g + Real.sqrt (1 + g ^ 2 - 2 * Real.cos α * g) := by
    intro g hg
    have : 1 - g ≤ Real.sqrt (1 + g ^ 2 - 2 * Real.cos α * g) := by
      rw [show 1 - g = Real.sqrt ((1 - g) ^ 2) by rw [Real.sqrt_sq (by linarith [hg.2])]]
      exact Real.sqrt_le_sqrt (by nlinarith [hg.1])
    linarith
  have hcont : ContinuousOn (fun g => H (c α) (g + Real.sqrt (1 + g ^ 2 - 2 * Real.cos α * g)))
      (Set.Icc 0 1) := by
    intro g hg
    exact ((continuousAt_H hc0 hc1 (hT1 g hg)).comp
      (f := fun g => g + Real.sqrt (1 + g ^ 2 - 2 * Real.cos α * g)) (by fun_prop)).continuousWithinAt
  have hderiv : ∀ g ∈ Set.Ioo (0 : ℝ) 1,
      HasDerivAt (fun g => H (c α) (g + Real.sqrt (1 + g ^ 2 - 2 * Real.cos α * g)))
        (Fd α 0 g) g := by
    intro g ⟨hg0, hg1⟩
    have hr0 : 0 < 1 + g ^ 2 - 2 * Real.cos α * g := rad_pos hα0 hα1
    have hTd : HasDerivAt (fun x => x + Real.sqrt (1 + x ^ 2 - 2 * Real.cos α * x))
        (1 + (2 * g - 2 * Real.cos α) / (2 * Real.sqrt (1 + g ^ 2 - 2 * Real.cos α * g))) g := by
      have h1 := ((((hasDerivAt_id g).pow 2).const_add 1).sub
        ((hasDerivAt_id g).const_mul (2 * Real.cos α))).sqrt (by exact hr0.ne')
      have h2 := (hasDerivAt_id g).add h1
      refine h2.congr_deriv ?_
      simp only [id, Nat.cast_ofNat, Pi.add_apply, Pi.sub_apply, Pi.pow_apply]; ring
    set R := Real.sqrt (1 + g ^ 2 - 2 * Real.cos α * g) with hRdef
    have hR0 : 0 < R := Real.sqrt_pos.mpr hr0
    have hRsq : R ^ 2 = 1 + g ^ 2 - 2 * Real.cos α * g := Real.sq_sqrt hr0.le
    have hRg : 1 - g < R := by
      rw [hRdef, show 1 - g = Real.sqrt ((1 - g) ^ 2) by rw [Real.sqrt_sq (by linarith)]]
      exact Real.sqrt_lt_sqrt (sq_nonneg _) (by nlinarith)
    set t := g + R with htdef
    have ht1 : 1 < t := by linarith
    have hHd := (hasDerivAt_H hc0 hc1 ht1).comp g hTd
    refine hHd.congr_deriv ?_
    rw [hκ] at hRsq
    obtain ⟨o1, o2, o3, o4, o5⟩ := order_facts hc0 hc1
    have htk : 0 < t - kap (c α) := by linarith
    have hT' : 1 + (2 * g - 2 * Real.cos α) / (2 * R) = (t - kap (c α)) / R := by
      rw [hκ, htdef]; field_simp; ring
    have hg_t : g = (t ^ 2 - 1) / (2 * (t - kap (c α))) := by
      rw [eq_div_iff (by positivity), htdef]; linear_combination (-1 : ℝ) * hRsq
    have hR_t : R = (t ^ 2 - 2 * kap (c α) * t + 1) / (2 * (t - kap (c α))) := by
      rw [eq_div_iff (by positivity), htdef]; linear_combination hRsq
    obtain ⟨eA, eB⟩ := alg_AB' hc0 hc1 ht1 hg_t hR_t
    have hlg : Real.log g = Real.log (t - 1) + Real.log (t - -1) - Real.log 2
        - Real.log (t - kap (c α)) := by
      rw [show g = (t - 1) * (t - -1) / (2 * (t - kap (c α))) by rw [hg_t]; ring,
        Real.log_div (by
          have : 0 < (t - 1) * (t - -1) := mul_pos (by linarith) (by linarith)
          exact this.ne') (by positivity),
        Real.log_mul (by linarith) (by linarith), Real.log_mul (by norm_num) htk.ne']
      ring
    have hl1g : Real.log (1 + g) = Real.log (t - r1 (c α)) + Real.log (t - r2 (c α))
        - Real.log 2 - Real.log (t - kap (c α)) := by
      rw [show 1 + g = (t - r1 (c α)) * (t - r2 (c α)) / (2 * (t - kap (c α))) by
          rw [hg_t]; unfold r1 r2 kap at *; field_simp; ring,
        Real.log_div (by
          have : 0 < (t - r1 (c α)) * (t - r2 (c α)) := mul_pos (by linarith) (by linarith)
          exact this.ne') (by positivity),
        Real.log_mul (by linarith) (by linarith), Real.log_mul (by norm_num) htk.ne']
      ring
    rw [hT']
    calc (cA (c α) t * (Real.log (t - 1) + Real.log (t - -1) - Real.log 2
            - Real.log (t - kap (c α)))
          + cB (c α) t * (Real.log (t - r1 (c α)) + Real.log (t - r2 (c α)) - Real.log 2
            - Real.log (t - kap (c α)))) * ((t - kap (c α)) / R)
        = cA (c α) t * ((t - kap (c α)) / R) * (Real.log (t - 1) + Real.log (t - -1)
            - Real.log 2 - Real.log (t - kap (c α)))
          + cB (c α) t * ((t - kap (c α)) / R) * (Real.log (t - r1 (c α))
            + Real.log (t - r2 (c α)) - Real.log 2 - Real.log (t - kap (c α))) := by ring
      _ = (-(g - kap (c α)) / ((1 + g) * R) + (R / (1 + g) - 1) / g) * Real.log g
          + -R / ((1 + g) * g) * Real.log (1 + g) := by rw [eA, eB, hlg, hl1g]
      _ = Fd α 0 g := by
          unfold Fd
          simp only [sub_zero, Real.rpow_one, Real.rpow_neg_one]
          rw [← hRdef, hκ]
          have : 0 < 1 + g := by linarith
          field_simp
          ring
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one hcont hderiv hint]
  have e2 : Real.sqrt (1 + 1 ^ 2 - 2 * Real.cos α * 1) = 2 * s α := by
    rw [show 1 + 1 ^ 2 - 2 * Real.cos α * 1 = (2 * s α) ^ 2 by
      rw [hκ]; unfold kap; nlinarith, Real.sqrt_sq (by linarith)]
  simp only [e2]
  norm_num

end D04

theorem D04.proof_4 : D04.stmt_4 := by
  intro α h0 h1
  obtain ⟨hint, hd⟩ := D04.hasDerivAt_Mint h0 h1
  rw [D04.integral_Fd0 h0 h1 hint, D04.H_end h0 h1] at hd
  exact hd

theorem D04.proof_5 : D04.stmt_5 := by
  intro α h0 h1
  rw [← D04.proof_1 α h0 h1]
  exact D04.proof_4 α h0 h1
