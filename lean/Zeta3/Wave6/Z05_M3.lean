import Mathlib
import Zeta3.Wave3.W02_Holo
import Zeta3.Wave2.D00_Dilog
import Zeta3.Wave5.Y01_M2
/-!
# Z05 -- the coefficient `M_3 = W_4 + W_3' + T` of the squircle expansion (Theorem `thm:M3`)

Target id: Z05.
Paper labels: `thm:M3` (main.tex 957-966), `eq:W3p` = (T2) of Table `tab:closed`
(main.tex 931-934), `eq:T` (main.tex 903-905), `thm:W4R` (main.tex 980-982) with
`eq:R` = (T3) of Table `tab:closed` (main.tex 935-946).
Detailed proof: paper/ieee/proofs/sec_g1_M3.tex (Lemmas `g1:lem:red3`, `g1:lem:pole`,
`g1:lem:W3p`, proof in `g1:sub:M3:proof`), with Lemma `g1:lem:Mk` of sec_g1_M2.tex and
Lemma `g1:lem:init` of sec_g1_wt2.tex.

Quoted verbatim (main.tex 957-966, Theorem thm:M3):

  \begin{theorem}[$M_3$$^\dagger$]\label{thm:M3}
  $M_3=W_4+W_3'+T$, with $T$ as in \eqref{eq:T}, $W_3'$ of weight three given by \eqref{eq:W3p}
  in Table~\ref{tab:closed}, and
  \begin{multline*}
  W_4=\int_b^a\Bigl[\ell_1^3\Bigl(\frac{8-7\sqrt2}{48v}-\frac{1}{3(v+a)}\Bigr)
  +\frac{\sqrt2}{16}\frac{\ell_1^2\ell_2}{v}\Bigr]dv,
  \end{multline*}
  where $\ell_1=\log\frac{(v-b)(v+a)}{4v}$ and $\ell_2=\log\frac{(a-v)(v+b)}{4v}$. Numerically
  $M_3=0.02277052149178657745\ldots$ and $W_4=0.79578505983224081132\ldots$
  \end{theorem}

Quoted verbatim (main.tex 903-905, eq:T):

  T=\frac{\pi^2+6\lambda^2}{24}=\frac{\zeta(2)-\Li_2(\frac12)}{2}

Quoted verbatim (main.tex 931-934, eq:W3p = (T2)):

  W_3'&=\Bigl(\frac{7\sqrt2}{32}-\frac3{16}\Bigr)\zeta(3)+(1+\sqrt2)\Li_3(b)+(1-\sqrt2)\Li_3(-b)-2\sqrt2\,\Li_3\bigl(\tfrac1{\sqrt2}\bigr)
  +(2L-\sqrt2\lambda)\Li_2(b)\notag\\
  &\quad+\Bigl(\frac{\sqrt2}{12}-\frac18\Bigr)\pi^2L-\frac{1+\sqrt2}{24}\pi^2\lambda+\frac{1+\sqrt2}{24}\lambda^3+\frac{\sqrt2}{4}\lambda^2L
  -\frac{1+\sqrt2}{2}\lambda L^2+\frac{5-\sqrt2}{6}L^3\tag{T2}\label{eq:W3p}

Quoted verbatim (main.tex 935-946, eq:R = (T3)):

  \mathcal R&=(1-573\sqrt2)\Li_4(b)-(1+573\sqrt2)\Li_4(-b)+56\sqrt2\,\Li_4(-b^2)
  +17\sqrt2\,\Li_4\bigl(\tfrac1{\sqrt2}\bigr)-\Bigl(\frac{37}{4}+\frac{17\sqrt2}{16}\Bigr)\Li_4\bigl(\tfrac12\bigr)\notag\\
  &\quad+(10+86\sqrt2)\Li_4\bigl(1-\tfrac1{\sqrt2}\bigr)-(10-86\sqrt2)\Li_4(2-\sqrt2)\notag\\
  &\quad+L\Bigl[(40+7\sqrt2)\Li_3(b)-(40-7\sqrt2)\Li_3(-b)-72\Li_3\bigl(\tfrac1{\sqrt2}\bigr)+\Bigl(\frac{63}{8}-\frac{3\sqrt2}{16}\Bigr)\zeta(3)\Bigr]\notag\\
  &\quad+\lambda\Bigl[(6+13\sqrt2)\Li_3(b)+(6-13\sqrt2)\Li_3(-b)-\frac{47\sqrt2}{2}\Li_3\bigl(\tfrac1{\sqrt2}\bigr)+\Bigl(\frac{329\sqrt2}{128}-\frac38\Bigr)\zeta(3)\Bigr]
  +30\Li_2(b)^2\notag\\
  &\quad+\Li_2(b)\Bigl[\Bigl(30+\frac{7\sqrt2}{2}\Bigr)L^2-24\lambda L-\frac{19\sqrt2}{8}\lambda^2+\Bigl(\frac{\sqrt2}{24}-\frac{19}{4}\Bigr)\pi^2\Bigr]
  +\Bigl(\frac{103}{384}-\frac{4963\sqrt2}{5760}\Bigr)\pi^4\notag\\
  &\quad+\pi^2\Bigl[\Bigl(\frac{67}{24}+\frac{185\sqrt2}{96}\Bigr)L^2+\Bigl(\frac{41\sqrt2}{12}-\frac{35}{12}\Bigr)\lambda L+\Bigl(\frac5{12}-\frac{1105\sqrt2}{384}\Bigr)\lambda^2\Bigr]
  +\Bigl(\frac{13}{24}-\frac{25\sqrt2}{24}\Bigr)L^4-\Bigl(\frac{73}{6}+\frac{29\sqrt2}{3}\Bigr)L^3\lambda\notag\\
  &\quad+\Bigl(\frac{63}{8}+\frac{67\sqrt2}{16}\Bigr)L^2\lambda^2+\Bigl(\frac{41}{24}+\frac{17\sqrt2}{6}\Bigr)L\lambda^3+\Bigl(\frac{137\sqrt2}{192}-\frac{61}{192}\Bigr)\lambda^4\tag{T3}\label{eq:R}

Quoted verbatim (main.tex 980-982, Theorem thm:W4R):

  \begin{theorem}[$M_3$ in closed form]\label{thm:W4R}
  $W_4=\mathcal R$; in particular $c_0=c_1=0$ in Theorem~\ref{thm:W4}. Hence
  $M_3=\mathcal R+W_3'+T$ is a closed form in classical polylogarithms, with $\mathcal R$,
  $W_3'$ and $T$ as in \eqref{eq:R}, \eqref{eq:W3p} and \eqref{eq:T}.
  \end{theorem}

Quoted verbatim (supplement sec_g1_M3.tex, integral form of `W_3'`, eq. g1:W3p):

  W_3'=\int_b^a\bigl[\ell_1^2\,\Omega_{20}+\ell_1\ell_2\,\Omega_{11}\bigr]dv,\\
  \Omega_{20}=\frac{1-\sqrt2}{4v}+\frac1{4(v-a)}-\frac1{4(v+a)}-\frac1{2(v+b)},\qquad
  \Omega_{11}=\frac{1+\sqrt2}{4v}-\frac1{2(v+a)} .

Quoted verbatim (supplement, Section g1:sub:M3:red):

  g_{11}(v)&=-\frac{\sqrt2\,(v^2+2v-1)}{4(v^2+1)}=-\frac{\sqrt2\,(v-b)(v+a)}{4(v^2+1)},\\
  G_{30}&=-\frac{2\sqrt2\,v^3}{3(v^2+1)^3}+\frac{\sqrt2\,(v-1)}{4(v^2+1)}+\frac{1+3\sqrt2}{24},\\
  G_{21}&=\frac{2\sqrt2\,v^3}{(v^2+1)^3}-\frac{\sqrt2\,(3v-1)}{4(v^2+1)}+\frac{1-\sqrt2}8 .
  ... and $G_{30}(b)=G_{21}(b)=g_{11}(b)=g_{20}(b)=0$

  R_{30}=G_{30}'-\frac1{3(v+a)}+\frac{8-7\sqrt2}{48v},\qquad R_{21}=G_{21}'+\frac{\sqrt2}{16v}.

  -3G_{30}D_1-G_{21}D_2=\Omega_{20}+g_{20}',\qquad -2G_{21}D_1=\Omega_{11}+g_{11}'

  So $M_3=W_4+\int_b^a[\ell_1^2(\Omega_{20}+g_{20}')+\ell_1\ell_2(\Omega_{11}+g_{11}')]dv$

  ... the degree-two parts contribute ... $\int\ell_1\Rwt=T$.

(`g_{20} = a/4 - √2/(2(v²+1))`, sec_g1_wt2.tex; `g_{20}' = √2 v/(v²+1)²`, Lemma g1:lem:k2.)

Lemma g1:lem:Mk (sec_g1_M2.tex): `M_k=-\frac12\int_0^1\varphi_k(t)\,dt`,
`\varphi_k(t)=[\varepsilon^k]\Phi_\varepsilon(t)`, with
`Φ_ε = A e^x + B e^y - S√(1+u)`, `x = εℓ₁`, `y = εℓ₂`, `A = cS = 1/t`, `B = c'S = 1/(1-t)`,
`u = c²(e^{2x}-1) + c'²(e^{2y}-1)`, `c = (1-t)/q`, `c' = t/q`, `S = q/(t(1-t))`.
The supplement does not print `φ₃`; it is the `ε³` coefficient of this expression, computed
here by hand (`√(1+u) = 1 + u/2 - u²/8 + u³/16 - ...`):
`φ₃ = (1/t)(1/6 - 2c/3 + c³ - c⁵/2) ℓ₁³ + (1/(1-t))(1/6 - 2c'/3 + c'³ - c'⁵/2) ℓ₂³
      + t(1-t)/q³ ((1 - 3c²/2) ℓ₁²ℓ₂ + (1 - 3c'²/2) ℓ₁ℓ₂²)`
(same shape as Y01's `φ₂`; checked numerically below).

## Conventions

* `M_3` is `W02.Mk 3 = -F'''(0)/(2·3!)` (complex, imported, not redefined); real numbers are
  cast to `ℂ`.
* Reused from Zeta3.Wave5.Y01_M2 (not redefined): `a = 1+√2`, `b = √2-1`, `L = log a`,
  `λ = log 2` (`Y01.lam`), `ζ(3) = Y01.z3 = (riemannZeta 3).re` (Y01.stmt_4, proved, says
  `riemannZeta 3 = ↑z3` and `HasSum (fun n => 1/(n+1)^3) z3`), the trilogarithm
  `Y01.Li3 x = ∫_0^x Li₂(t)/t dt` (Y01.stmt_1-3: integrable, the series `Σ x^k/k³` for `|x|≤1`,
  derivative `Li₂(x)/x`), `Y01.tv`, `Y01.l1`, `Y01.l2`, `Y01.c`, `Y01.c'`, `Y01.g20'`,
  `Y01.Rwt`, `Y01.Tint` (`T = ∫_b^a ℓ₁ R_wt dv`). `Li₂ = D00.Li2`.
* `Y01.l1 v = log t(v)` with `t(v) = (v²+2v-1)/(4v)`, and `(v-b)(v+a) = v²+2v-1`, so
  `Y01.l1` is literally the paper's `ℓ₁ = log((v-b)(v+a)/(4v))`. `Y01.l2 v = log(1 - t(v))` and
  `4v - (v²+2v-1) = 1+2v-v² = (a-v)(v+b)`, so `Y01.l2` is the paper's
  `ℓ₂ = log((a-v)(v+b)/(4v))`. On `(b,a)` both arguments lie in `(0,1)` (genuine logarithms).
* `T = (π²+6λ²)/24` is Y01.stmt_8 (proved there: `Y01.Tint = (π²+6λ²)/24`), not restated.
* `1/(v-a)` in `Ω₂₀` and `R_wt` is Lean's `0` at `v = a` only, a null set.
* Every integral whose value is asserted comes with an explicit `IntervalIntegrable` claim.
* The tetralogarithm (needed only for thm:W4R) is defined as `Li4 x = ∫_0^x Li₃(t)/t dt`;
  stmt_10 records integrability and the paper's series `Σ x^k/k⁴` for `|x| ≤ 1`. All seven
  arguments of (T3) lie in `(-1,1)`.

## Statements

* stmt_1: the `W_4` integrand is integrable on `[b,a]` ("All integrals converge").
* stmt_2: Lemma g1:lem:Mk at `k = 3`: `φ₃` integrable on `[0,1]`, `M_3 = -½∫_0^1 φ₃`.
* stmt_3: Step 1 of g1:lem:red3 (Euler substitution and folding):
  `M_3 = ∫_b^a [ℓ₁³ R₃₀ + ℓ₁²ℓ₂ R₂₁] dv`, with `R₃₀ = G₃₀' - 1/(3(v+a)) + (8-7√2)/(48v)`,
  `R₂₁ = G₂₁' + √2/(16v)`.
* stmt_4: the kernel facts of g1:sub:M3:red: derivatives of `G₃₀, G₂₁, g₂₀, g₁₁, ℓ₁, ℓ₂`,
  vanishing at `b`, and the two identities `-3G₃₀D₁-G₂₁D₂ = Ω₂₀+g₂₀'`, `-2G₂₁D₁ = Ω₁₁+g₁₁'`
  on `(b,a)`.
* stmt_5: Step 2 of g1:lem:red3 (degree-three integration by parts):
  `M_3 = W_4 + ∫_b^a [ℓ₁²(Ω₂₀+g₂₀') + ℓ₁ℓ₂(Ω₁₁+g₁₁')] dv`.
* stmt_6: Step 3 of g1:lem:red3 (degree-two integration by parts and folding):
  `∫_b^a [ℓ₁² g₂₀' + ℓ₁ℓ₂ g₁₁'] dv = T` (`= Y01.Tint`).
* stmt_7: Lemma g1:lem:red3: `M_3 = W_4 + W_3'(integral) + T`, all integrals convergent.
* stmt_8: Lemma g1:lem:W3p: the integral g1:W3p equals (T2).
* stmt_9: Theorem thm:M3: `M_3 = W_4 + (T2) + (π²+6λ²)/24`.
* stmt_10: properties of `Li4` relied on for (T3): integrability of `Li₃(t)/t` on `[0,x]`,
  `x ≤ 1`, and `Li4 x = Σ_{k≥1} x^k/k⁴` for `|x| ≤ 1`.
* stmt_11: Theorem thm:W4R, the "Hence" step only, as an implication: if `W_4 = 𝓡` (T3), then
  `M_3 = 𝓡 + W_3' + T`. The identity `W_4 = 𝓡` itself is NOT COVERED: the paper proves it by a
  computer certificate over 5184 level-8 words, which this file does not reproduce.

Numerical check before freezing (mpmath, 45 digits): `W_4` (quadrature) =
0.795785059832240811320359819629726762546...; (T2) = -1.304361308532061199148508733130571...;
`W_4 + (T2) + T` = 0.022770521491786577456730509742327800250...; `-½∫φ₃` (stmt_2) and the
integral of stmt_3 agree with this to `10⁻⁴⁰`; the right side of stmt_5 to `2·10⁻⁴¹`; the
integral of stmt_8 agrees with (T2) to `10⁻⁴²`; `Tint - T` is `4·10⁻⁴⁵`; independently,
`-F'''(0)/(2·3!)` by Cauchy's formula on `|ε| = 1/2` (64 nodes, 30 digits, from W02's `Φ`)
gives 0.0227705214917865774567305100647, agreeing to `3·10⁻²⁸` (truncation of the Cauchy sum).
The identities of stmt_4 were checked exactly with sympy. (T3) agrees with `W_4` to
`2·10⁻⁴¹`, and the integral definition of `Li4` agrees with mpmath's `polylog(4,·)` at `±b`.

## NOT COVERED

* The printed decimals of `M_3` and `W_4` (not required for this target).
* "of weight three": a property of the displayed expression (T2), not of the number.
* Lemma g1:lem:pole and the raw table of `W_3'` (18 trilogarithms, 18 dilogarithms, Table
  g1:tab:rawW3p) and its reduction by the tables g1:tab:dilog / g1:tab:trilog: computer
  generated tables, not reproduced; stmt_8 states their combined conclusion. The trilogarithm
  functional equations and Spence--Kummer that this reduction uses are Y01.stmt_9 and
  Y01.stmt_10 (proved there) and are not restated.
* thm:W4R's identity `W_4 = 𝓡` itself (a computer certificate in the paper; `stmt_11` only
  assumes it). The statement block was changed by the planner on 2026-11-26 after the audit: the
  audited `stmt_11` (`W4 = Rclosed`) was removed as unproved and the audited `stmt_12` became the
  implication `stmt_11`; `stmt_1` to `stmt_10` and all definitions are unchanged.
* thm:W4R's "in particular `c_0 = c_1 = 0` in Theorem thm:W4" (thm:W4 is a separate target)
  and the description of the certificate (level-8 words and relations).
-/
-- BEGIN STATEMENTS
namespace Z05

open MeasureTheory Real

/-- `φ₃ = [ε³]Φ_ε` (Lemma g1:lem:Mk at `k = 3`, expanded as in the docstring). -/
noncomputable def phi3 (t : ℝ) : ℝ :=
  1 / t * (1 / 6 - 2 * Y01.c t / 3 + Y01.c t ^ 3 - Y01.c t ^ 5 / 2) * Real.log t ^ 3
    + 1 / (1 - t) * (1 / 6 - 2 * Y01.c' t / 3 + Y01.c' t ^ 3 - Y01.c' t ^ 5 / 2)
      * Real.log (1 - t) ^ 3
    + t * (1 - t) / W02.q t ^ 3
      * ((1 - 3 * Y01.c t ^ 2 / 2) * Real.log t ^ 2 * Real.log (1 - t)
        + (1 - 3 * Y01.c' t ^ 2 / 2) * Real.log t * Real.log (1 - t) ^ 2)

/-- The integrand of `W_4` (Theorem thm:M3). -/
noncomputable def W4integrand (v : ℝ) : ℝ :=
  Y01.l1 v ^ 3 * ((8 - 7 * Real.sqrt 2) / (48 * v) - 1 / (3 * (v + Y01.a)))
    + Real.sqrt 2 / 16 * (Y01.l1 v ^ 2 * Y01.l2 v / v)

/-- `W_4 = ∫_b^a [ℓ₁³((8-7√2)/(48v) - 1/(3(v+a))) + (√2/16) ℓ₁²ℓ₂/v] dv`. -/
noncomputable def W4 : ℝ := ∫ v in Y01.b..Y01.a, W4integrand v

/-- The right-hand side of eq:W3p = (T2) of Table tab:closed. -/
noncomputable def W3pclosed : ℝ :=
  (7 * Real.sqrt 2 / 32 - 3 / 16) * Y01.z3
    + (1 + Real.sqrt 2) * Y01.Li3 Y01.b + (1 - Real.sqrt 2) * Y01.Li3 (-Y01.b)
    - 2 * Real.sqrt 2 * Y01.Li3 (1 / Real.sqrt 2)
    + (2 * Y01.L - Real.sqrt 2 * Y01.lam) * D00.Li2 Y01.b
    + (Real.sqrt 2 / 12 - 1 / 8) * π ^ 2 * Y01.L
    - (1 + Real.sqrt 2) / 24 * π ^ 2 * Y01.lam
    + (1 + Real.sqrt 2) / 24 * Y01.lam ^ 3
    + Real.sqrt 2 / 4 * Y01.lam ^ 2 * Y01.L
    - (1 + Real.sqrt 2) / 2 * Y01.lam * Y01.L ^ 2
    + (5 - Real.sqrt 2) / 6 * Y01.L ^ 3

/-- `Ω₂₀ = (1-√2)/(4v) + 1/(4(v-a)) - 1/(4(v+a)) - 1/(2(v+b))`. -/
noncomputable def Om20 (v : ℝ) : ℝ :=
  (1 - Real.sqrt 2) / (4 * v) + 1 / (4 * (v - Y01.a)) - 1 / (4 * (v + Y01.a))
    - 1 / (2 * (v + Y01.b))

/-- `Ω₁₁ = (1+√2)/(4v) - 1/(2(v+a))`. -/
noncomputable def Om11 (v : ℝ) : ℝ := (1 + Real.sqrt 2) / (4 * v) - 1 / (2 * (v + Y01.a))

/-- The integrand of the integral form g1:W3p of `W_3'`. -/
noncomputable def W3pintegrand (v : ℝ) : ℝ :=
  Y01.l1 v ^ 2 * Om20 v + Y01.l1 v * Y01.l2 v * Om11 v

/-- `W_3'` as the integral g1:W3p. -/
noncomputable def W3pint : ℝ := ∫ v in Y01.b..Y01.a, W3pintegrand v

/-- `G₃₀ = -2√2v³/(3(v²+1)³) + √2(v-1)/(4(v²+1)) + (1+3√2)/24`. -/
noncomputable def G30 (v : ℝ) : ℝ :=
  -(2 * Real.sqrt 2 * v ^ 3) / (3 * (v ^ 2 + 1) ^ 3) + Real.sqrt 2 * (v - 1) / (4 * (v ^ 2 + 1))
    + (1 + 3 * Real.sqrt 2) / 24

/-- `G₂₁ = 2√2v³/(v²+1)³ - √2(3v-1)/(4(v²+1)) + (1-√2)/8`. -/
noncomputable def G21 (v : ℝ) : ℝ :=
  2 * Real.sqrt 2 * v ^ 3 / (v ^ 2 + 1) ^ 3 - Real.sqrt 2 * (3 * v - 1) / (4 * (v ^ 2 + 1))
    + (1 - Real.sqrt 2) / 8

/-- `g₂₀ = a/4 - √2/(2(v²+1))` (sec_g1_wt2.tex; `g₂₀' = Y01.g20'`). -/
noncomputable def g20 (v : ℝ) : ℝ := Y01.a / 4 - Real.sqrt 2 / (2 * (v ^ 2 + 1))

/-- `g₁₁ = -√2(v²+2v-1)/(4(v²+1))`. -/
noncomputable def g11 (v : ℝ) : ℝ := -(Real.sqrt 2 * (v ^ 2 + 2 * v - 1)) / (4 * (v ^ 2 + 1))

/-- `G₃₀'`, written out: `-√2(v⁶-2v⁵-7v⁴-4v³+7v²-2v-1)/(4(v²+1)⁴)`. -/
noncomputable def G30d (v : ℝ) : ℝ :=
  -(Real.sqrt 2 * (v ^ 6 - 2 * v ^ 5 - 7 * v ^ 4 - 4 * v ^ 3 + 7 * v ^ 2 - 2 * v - 1))
    / (4 * (v ^ 2 + 1) ^ 4)

/-- `G₂₁'`, written out: `√2(3v⁶-2v⁵-21v⁴-4v³+21v²-2v-3)/(4(v²+1)⁴)`. -/
noncomputable def G21d (v : ℝ) : ℝ :=
  Real.sqrt 2 * (3 * v ^ 6 - 2 * v ^ 5 - 21 * v ^ 4 - 4 * v ^ 3 + 21 * v ^ 2 - 2 * v - 3)
    / (4 * (v ^ 2 + 1) ^ 4)

/-- `g₁₁'`, written out: `√2(v²-2v-1)/(2(v²+1)²)`. -/
noncomputable def g11d (v : ℝ) : ℝ := Real.sqrt 2 * (v ^ 2 - 2 * v - 1) / (2 * (v ^ 2 + 1) ^ 2)

/-- `D₁ = ℓ₁' = (2v+2)/(v²+2v-1) - 1/v`. -/
noncomputable def D1 (v : ℝ) : ℝ := (2 * v + 2) / (v ^ 2 + 2 * v - 1) - 1 / v

/-- `D₂ = ℓ₂' = (2-2v)/(1+2v-v²) - 1/v`. -/
noncomputable def D2 (v : ℝ) : ℝ := (2 - 2 * v) / (1 + 2 * v - v ^ 2) - 1 / v

/-- `R₃₀ = G₃₀' - 1/(3(v+a)) + (8-7√2)/(48v)`. -/
noncomputable def R30 (v : ℝ) : ℝ := G30d v - 1 / (3 * (v + Y01.a)) + (8 - 7 * Real.sqrt 2) / (48 * v)

/-- `R₂₁ = G₂₁' + √2/(16v)`. -/
noncomputable def R21 (v : ℝ) : ℝ := G21d v + Real.sqrt 2 / (16 * v)

/-- The real tetralogarithm `Li₄(x) = ∫_0^x Li₃(t)/t dt` (only for thm:W4R). -/
noncomputable def Li4 (x : ℝ) : ℝ := ∫ t in (0 : ℝ)..x, Y01.Li3 t / t

/-- The right-hand side of eq:R = (T3) of Table tab:closed. -/
noncomputable def Rclosed : ℝ :=
  (1 - 573 * Real.sqrt 2) * Li4 Y01.b - (1 + 573 * Real.sqrt 2) * Li4 (-Y01.b)
    + 56 * Real.sqrt 2 * Li4 (-Y01.b ^ 2)
    + 17 * Real.sqrt 2 * Li4 (1 / Real.sqrt 2)
    - (37 / 4 + 17 * Real.sqrt 2 / 16) * Li4 (1 / 2)
    + (10 + 86 * Real.sqrt 2) * Li4 (1 - 1 / Real.sqrt 2)
    - (10 - 86 * Real.sqrt 2) * Li4 (2 - Real.sqrt 2)
    + Y01.L * ((40 + 7 * Real.sqrt 2) * Y01.Li3 Y01.b - (40 - 7 * Real.sqrt 2) * Y01.Li3 (-Y01.b)
        - 72 * Y01.Li3 (1 / Real.sqrt 2) + (63 / 8 - 3 * Real.sqrt 2 / 16) * Y01.z3)
    + Y01.lam * ((6 + 13 * Real.sqrt 2) * Y01.Li3 Y01.b + (6 - 13 * Real.sqrt 2) * Y01.Li3 (-Y01.b)
        - 47 * Real.sqrt 2 / 2 * Y01.Li3 (1 / Real.sqrt 2)
        + (329 * Real.sqrt 2 / 128 - 3 / 8) * Y01.z3)
    + 30 * D00.Li2 Y01.b ^ 2
    + D00.Li2 Y01.b * ((30 + 7 * Real.sqrt 2 / 2) * Y01.L ^ 2 - 24 * Y01.lam * Y01.L
        - 19 * Real.sqrt 2 / 8 * Y01.lam ^ 2 + (Real.sqrt 2 / 24 - 19 / 4) * π ^ 2)
    + (103 / 384 - 4963 * Real.sqrt 2 / 5760) * π ^ 4
    + π ^ 2 * ((67 / 24 + 185 * Real.sqrt 2 / 96) * Y01.L ^ 2
        + (41 * Real.sqrt 2 / 12 - 35 / 12) * Y01.lam * Y01.L
        + (5 / 12 - 1105 * Real.sqrt 2 / 384) * Y01.lam ^ 2)
    + (13 / 24 - 25 * Real.sqrt 2 / 24) * Y01.L ^ 4
    - (73 / 6 + 29 * Real.sqrt 2 / 3) * Y01.L ^ 3 * Y01.lam
    + (63 / 8 + 67 * Real.sqrt 2 / 16) * Y01.L ^ 2 * Y01.lam ^ 2
    + (41 / 24 + 17 * Real.sqrt 2 / 6) * Y01.L * Y01.lam ^ 3
    + (137 * Real.sqrt 2 / 192 - 61 / 192) * Y01.lam ^ 4

/-- Supplement, Lemma g1:lem:red3, "All integrals converge": the integrand of `W_4` is
integrable on `[b,a]`. -/
def stmt_1 : Prop := IntervalIntegrable W4integrand volume Y01.b Y01.a

/-- Supplement, Lemma g1:lem:Mk at `k = 3` ("$M_k=-\frac12\int_0^1\varphi_k(t)\,dt$"), with
`φ₃` written out: `φ₃` is integrable on `[0,1]` and `M_3 = -½ ∫_0^1 φ₃(t) dt`. -/
def stmt_2 : Prop :=
  IntervalIntegrable phi3 volume 0 1 ∧
    W02.Mk 3 = (((-(1 / 2) * ∫ t in (0 : ℝ)..1, phi3 t) : ℝ) : ℂ)

/-- Supplement, Lemma g1:lem:red3, Step 1: "$M_3=\int_b^a[\ell_1^3R_{30}+\ell_1^2\ell_2R_{21}]dv$"
with the folded kernels `R₃₀ = G₃₀' - 1/(3(v+a)) + (8-7√2)/(48v)`, `R₂₁ = G₂₁' + √2/(16v)`. -/
def stmt_3 : Prop :=
  IntervalIntegrable (fun v => Y01.l1 v ^ 3 * R30 v + Y01.l1 v ^ 2 * Y01.l2 v * R21 v)
      volume Y01.b Y01.a ∧
    W02.Mk 3 =
      (((∫ v in Y01.b..Y01.a, (Y01.l1 v ^ 3 * R30 v + Y01.l1 v ^ 2 * Y01.l2 v * R21 v)) : ℝ) : ℂ)

/-- Supplement, Section g1:sub:M3:red: `G₃₀' , G₂₁', g₂₀', g₁₁'` are the displayed rational
functions; "$G_{30}(b)=G_{21}(b)=g_{11}(b)=g_{20}(b)=0$"; `ℓ₁' = D₁`, `ℓ₂' = D₂` on `(b,a)`; and
"$-3G_{30}D_1-G_{21}D_2=\Omega_{20}+g_{20}',\ -2G_{21}D_1=\Omega_{11}+g_{11}'$" on `(b,a)`. -/
def stmt_4 : Prop :=
  (∀ v : ℝ, HasDerivAt G30 (G30d v) v ∧ HasDerivAt G21 (G21d v) v ∧
      HasDerivAt g20 (Y01.g20' v) v ∧ HasDerivAt g11 (g11d v) v) ∧
    G30 Y01.b = 0 ∧ G21 Y01.b = 0 ∧ g11 Y01.b = 0 ∧ g20 Y01.b = 0 ∧
    ∀ v : ℝ, Y01.b < v → v < Y01.a →
      HasDerivAt Y01.l1 (D1 v) v ∧ HasDerivAt Y01.l2 (D2 v) v ∧
      -3 * G30 v * D1 v - G21 v * D2 v = Om20 v + Y01.g20' v ∧
      -2 * G21 v * D1 v = Om11 v + g11d v

/-- Supplement, Lemma g1:lem:red3, Step 2 (integration by parts at degree three):
"$M_3=W_4+\int_b^a[\ell_1^2(\Omega_{20}+g_{20}')+\ell_1\ell_2(\Omega_{11}+g_{11}')]dv$". -/
def stmt_5 : Prop :=
  IntervalIntegrable
      (fun v => Y01.l1 v ^ 2 * (Om20 v + Y01.g20' v) + Y01.l1 v * Y01.l2 v * (Om11 v + g11d v))
      volume Y01.b Y01.a ∧
    W02.Mk 3 = (((W4 + ∫ v in Y01.b..Y01.a,
      (Y01.l1 v ^ 2 * (Om20 v + Y01.g20' v) + Y01.l1 v * Y01.l2 v * (Om11 v + g11d v))) : ℝ) : ℂ)

/-- Supplement, Lemma g1:lem:red3, Step 3 (integration by parts at degree two and folding):
"the degree-two parts contribute ... $\int\ell_1\Rwt=T$", that is
`∫_b^a [ℓ₁² g₂₀' + ℓ₁ℓ₂ g₁₁'] dv = ∫_b^a ℓ₁ R_wt dv`. -/
def stmt_6 : Prop :=
  IntervalIntegrable (fun v => Y01.l1 v ^ 2 * Y01.g20' v + Y01.l1 v * Y01.l2 v * g11d v)
      volume Y01.b Y01.a ∧
    ∫ v in Y01.b..Y01.a, (Y01.l1 v ^ 2 * Y01.g20' v + Y01.l1 v * Y01.l2 v * g11d v) = Y01.Tint

/-- Supplement, Lemma g1:lem:red3: "$M_3=W_4+W_3'+T$ with $W_4$, $T$ as in the theorem and
Lemma~\ref{g1:lem:red2}, and $W_3'$ the integral \eqref{g1:W3p}. All integrals converge." -/
def stmt_7 : Prop :=
  IntervalIntegrable W4integrand volume Y01.b Y01.a ∧
    IntervalIntegrable W3pintegrand volume Y01.b Y01.a ∧
    IntervalIntegrable (fun v => Y01.l1 v * Y01.Rwt v) volume Y01.b Y01.a ∧
    W02.Mk 3 = (((W4 + W3pint + Y01.Tint) : ℝ) : ℂ)

/-- Supplement, Lemma g1:lem:W3p: the integral \eqref{g1:W3p} "equals the right-hand side of
\eqref{g1:T2}". -/
def stmt_8 : Prop := W3pint = W3pclosed

/-- Theorem thm:M3: "$M_3=W_4+W_3'+T$", with `W_3'` = (T2) and `T = (π²+6λ²)/24` (eq:T), for
`M_3 = W02.Mk 3`. -/
def stmt_9 : Prop :=
  W02.Mk 3 = (((W4 + W3pclosed + (π ^ 2 + 6 * Y01.lam ^ 2) / 24) : ℝ) : ℂ)

/-- Properties of `Li4` relied on (main.tex: "$\Li_n(z)=\sum_{k\ge1}z^k/k^n$" at `n = 4`):
`Li₃(t)/t` is integrable on `[0,x]` for `x ≤ 1`, and `Li4 x = Σ_{k≥1} x^k/k⁴` for `|x| ≤ 1`. -/
def stmt_10 : Prop :=
  (∀ x : ℝ, x ≤ 1 → IntervalIntegrable (fun t : ℝ => Y01.Li3 t / t) volume 0 x) ∧
    ∀ x : ℝ, |x| ≤ 1 → HasSum (fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) ^ 4) (Li4 x)

/-- Theorem thm:W4R, the step "Hence $M_3=\mathcal R+W_3'+T$", with (T3), (T2), eq:T, stated as
an implication from the certificate identity "$W_4=\mathcal R$", which is not proved here. -/
def stmt_11 : Prop :=
  W4 = Rclosed →
    W02.Mk 3 = (((Rclosed + W3pclosed + (π ^ 2 + 6 * Y01.lam ^ 2) / 24) : ℝ) : ℂ)

end Z05
-- END STATEMENTS
set_option linter.unusedSimpArgs false

namespace Z05.P4

open MeasureTheory Real Filter Topology Set

lemma s2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)

lemma hG30 (v : ℝ) : HasDerivAt G30 (G30d v) v := by
  have hp : (v ^ 2 + 1) ≠ 0 := by positivity
  have e1 : HasDerivAt (fun w : ℝ => w ^ 2 + 1) (2 * v) v := by
    simpa using (hasDerivAt_pow 2 v).add_const 1
  have e3 : HasDerivAt (fun w : ℝ => w ^ 3) (3 * v ^ 2) v := by simpa using hasDerivAt_pow 3 v
  have f1 : HasDerivAt (fun w : ℝ => -(2 * Real.sqrt 2 * w ^ 3) / (3 * (w ^ 2 + 1) ^ 3))
      ((-(2 * Real.sqrt 2 * (3 * v ^ 2)) * (3 * (v ^ 2 + 1) ^ 3)
        - (-(2 * Real.sqrt 2 * v ^ 3)) * (3 * (3 * (v ^ 2 + 1) ^ 2 * (2 * v))))
        / (3 * (v ^ 2 + 1) ^ 3) ^ 2) v := by
    apply HasDerivAt.div
    · exact (e3.const_mul (2 * Real.sqrt 2)).neg
    · exact (e1.pow 3).const_mul 3
    · exact (by positivity : (3 * (v ^ 2 + 1) ^ 3) ≠ 0)
  have f2 : HasDerivAt (fun w : ℝ => Real.sqrt 2 * (w - 1) / (4 * (w ^ 2 + 1)))
      ((Real.sqrt 2 * 1 * (4 * (v ^ 2 + 1)) - Real.sqrt 2 * (v - 1) * (4 * (2 * v)))
        / (4 * (v ^ 2 + 1)) ^ 2) v := by
    apply HasDerivAt.div
    · exact ((hasDerivAt_id v).sub_const 1).const_mul _
    · exact e1.const_mul 4
    · exact (by positivity : (4 * (v ^ 2 + 1)) ≠ 0)
  have := (f1.add f2).add_const ((1 + 3 * Real.sqrt 2) / 24)
  unfold G30
  convert this using 1
  unfold G30d
  field_simp
  ring

lemma hG21 (v : ℝ) : HasDerivAt G21 (G21d v) v := by
  have hp : (v ^ 2 + 1) ≠ 0 := by positivity
  have e1 : HasDerivAt (fun w : ℝ => w ^ 2 + 1) (2 * v) v := by
    simpa using (hasDerivAt_pow 2 v).add_const 1
  have e3 : HasDerivAt (fun w : ℝ => w ^ 3) (3 * v ^ 2) v := by simpa using hasDerivAt_pow 3 v
  have f1 : HasDerivAt (fun w : ℝ => 2 * Real.sqrt 2 * w ^ 3 / (w ^ 2 + 1) ^ 3)
      (((2 * Real.sqrt 2 * (3 * v ^ 2)) * ((v ^ 2 + 1) ^ 3)
        - (2 * Real.sqrt 2 * v ^ 3) * (3 * (v ^ 2 + 1) ^ 2 * (2 * v)))
        / ((v ^ 2 + 1) ^ 3) ^ 2) v := by
    apply HasDerivAt.div
    · exact (e3.const_mul (2 * Real.sqrt 2))
    · exact (e1.pow 3)
    · exact (by positivity : ((v ^ 2 + 1) ^ 3) ≠ 0)
  have f2 : HasDerivAt (fun w : ℝ => Real.sqrt 2 * (3 * w - 1) / (4 * (w ^ 2 + 1)))
      ((Real.sqrt 2 * (3 * 1) * (4 * (v ^ 2 + 1)) - Real.sqrt 2 * (3 * v - 1) * (4 * (2 * v)))
        / (4 * (v ^ 2 + 1)) ^ 2) v := by
    apply HasDerivAt.div
    · exact (((hasDerivAt_id v).const_mul 3).sub_const 1).const_mul _
    · exact e1.const_mul 4
    · exact (by positivity : (4 * (v ^ 2 + 1)) ≠ 0)
  have := (f1.sub f2).add_const ((1 - Real.sqrt 2) / 8)
  unfold G21
  convert this using 1
  unfold G21d
  field_simp
  ring

lemma hg20 (v : ℝ) : HasDerivAt g20 (Y01.g20' v) v := by
  have h := Y01.PC.hasDerivAt_g20 v
  have : g20 = Y01.PC.g20 := by funext w; rfl
  rw [this]; exact h

lemma hg11 (v : ℝ) : HasDerivAt g11 (g11d v) v := by
  have hp : (v ^ 2 + 1) ≠ 0 := by positivity
  have e1 : HasDerivAt (fun w : ℝ => w ^ 2 + 1) (2 * v) v := by
    simpa using (hasDerivAt_pow 2 v).add_const 1
  have e2 : HasDerivAt (fun w : ℝ => w ^ 2 + 2 * w - 1) (2 * v + 2) v := by
    have := ((hasDerivAt_pow 2 v).add ((hasDerivAt_id v).const_mul 2)).sub_const 1
    convert this using 1; norm_num; ring
  have f1 : HasDerivAt (fun w : ℝ => -(Real.sqrt 2 * (w ^ 2 + 2 * w - 1)) / (4 * (w ^ 2 + 1)))
      ((-(Real.sqrt 2 * (2 * v + 2)) * (4 * (v ^ 2 + 1))
        - (-(Real.sqrt 2 * (v ^ 2 + 2 * v - 1))) * (4 * (2 * v))) / (4 * (v ^ 2 + 1)) ^ 2) v := by
    apply HasDerivAt.div
    · exact (e2.const_mul _).neg
    · exact e1.const_mul 4
    · exact (by positivity : (4 * (v ^ 2 + 1)) ≠ 0)
  unfold g11
  convert f1 using 1
  unfold g11d
  field_simp
  ring

lemma b_sq : Y01.b ^ 2 + 1 = 4 - 2 * Real.sqrt 2 := by
  unfold Y01.b; nlinarith [s2]

lemma vals : G30 Y01.b = 0 ∧ G21 Y01.b = 0 ∧ g11 Y01.b = 0 ∧ g20 Y01.b = 0 := by
  have hs := s2
  have hp : (4 - 2 * Real.sqrt 2) ≠ 0 := by nlinarith [Y01.PC.s2lt]
  refine ⟨?_, ?_, ?_, ?_⟩
  · unfold G30; rw [b_sq]; field_simp; unfold Y01.b; ring_nf
    simp only [Y01.PE1.s2_3, Y01.PE1.s2_4, Y01.PE1.s2_5, Y01.PE1.s2_6, Y01.PE1.s2_7, hs]
    ring_nf
  · unfold G21; rw [b_sq]; field_simp; unfold Y01.b; ring_nf
    simp only [Y01.PE1.s2_3, Y01.PE1.s2_4, Y01.PE1.s2_5, Y01.PE1.s2_6, Y01.PE1.s2_7, hs]
    ring_nf
  · unfold g11; rw [b_sq]; field_simp; unfold Y01.b; ring_nf
    simp only [Y01.PE1.s2_3, Y01.PE1.s2_4, Y01.PE1.s2_5, Y01.PE1.s2_6, Y01.PE1.s2_7, hs]
    ring_nf
  · exact Y01.PC.g20_b

lemma fac1 (v : ℝ) : v ^ 2 + 2 * v - 1 = (v - Y01.b) * (v + Y01.a) := by
  unfold Y01.a Y01.b; nlinarith [s2]

lemma fac2 (v : ℝ) : 1 + 2 * v - v ^ 2 = -((v - Y01.a) * (v + Y01.b)) := by
  unfold Y01.a Y01.b; nlinarith [s2]

lemma D1_eq {v : ℝ} (hv : v ∈ Ioo Y01.b Y01.a) : D1 v = Y01.PC.D1 v := by
  have hv0 : 0 < v := Y01.PC.b_pos.trans hv.1
  have h1 : v - Y01.b ≠ 0 := by linarith [hv.1]
  have h2 : v + Y01.a ≠ 0 := by linarith [Y01.PC.a_pos]
  unfold D1 Y01.PC.D1; rw [fac1]
  have : (2 * v + 2) = (v - Y01.b) + (v + Y01.a) := by linarith [Y01.PC.hamb]
  rw [this]; field_simp; ring

lemma D2_eq {v : ℝ} (hv : v ∈ Ioo Y01.b Y01.a) : D2 v = Y01.PC.D2 v := by
  have hv0 : 0 < v := Y01.PC.b_pos.trans hv.1
  have h1 : v - Y01.a ≠ 0 := by linarith [hv.2]
  have h2 : v + Y01.b ≠ 0 := by linarith [Y01.PC.b_pos]
  unfold D2 Y01.PC.D2; rw [fac2]
  have : (2 - 2 * v) = -((v - Y01.a) + (v + Y01.b)) := by linarith [Y01.PC.hamb]
  rw [this]; field_simp; ring

lemma id1 {v : ℝ} (hv : v ∈ Ioo Y01.b Y01.a) :
    -3 * G30 v * D1 v - G21 v * D2 v = Om20 v + Y01.g20' v := by
  have hs := s2
  have hv0 : 0 < v := Y01.PC.b_pos.trans hv.1
  have h1 : v - Y01.b ≠ 0 := by linarith [hv.1]
  have h2 : v + Y01.a ≠ 0 := by linarith [Y01.PC.a_pos]
  have h3 : v - Y01.a ≠ 0 := by linarith [hv.2]
  have h4 : v + Y01.b ≠ 0 := by linarith [Y01.PC.b_pos]
  have h5 : v ^ 2 + 1 ≠ 0 := by positivity
  rw [D1_eq hv, D2_eq hv]
  unfold G30 G21 Y01.PC.D1 Y01.PC.D2 Om20 Y01.g20'
  rw [div_add_div _ _ h1 h2, div_add_div _ _ h3 h4]
  field_simp
  unfold Y01.a Y01.b
  ring_nf
  simp only [Y01.PE1.s2_3, Y01.PE1.s2_4, Y01.PE1.s2_5, Y01.PE1.s2_6, Y01.PE1.s2_7, hs]
  ring_nf

lemma id2 {v : ℝ} (hv : v ∈ Ioo Y01.b Y01.a) :
    -2 * G21 v * D1 v = Om11 v + g11d v := by
  have hs := s2
  have hv0 : 0 < v := Y01.PC.b_pos.trans hv.1
  have h1 : v - Y01.b ≠ 0 := by linarith [hv.1]
  have h2 : v + Y01.a ≠ 0 := by linarith [Y01.PC.a_pos]
  have h5 : v ^ 2 + 1 ≠ 0 := by positivity
  rw [D1_eq hv]
  unfold G21 Y01.PC.D1 Om11 g11d
  field_simp
  unfold Y01.a Y01.b
  ring_nf
  simp only [Y01.PE1.s2_3, Y01.PE1.s2_4, Y01.PE1.s2_5, Y01.PE1.s2_6, Y01.PE1.s2_7, hs]
  ring_nf


theorem main : stmt_4 := by
  refine ⟨fun v => ⟨hG30 v, hG21 v, hg20 v, hg11 v⟩, vals.1, vals.2.1, vals.2.2.1, vals.2.2.2,
    fun v hb ha => ?_⟩
  have hv : v ∈ Ioo Y01.b Y01.a := ⟨hb, ha⟩
  refine ⟨?_, ?_, id1 hv, id2 hv⟩
  · rw [D1_eq hv]; exact Y01.PC.hasDerivAt_l1 hv
  · rw [D2_eq hv]; exact Y01.PC.hasDerivAt_l2 hv

end Z05.P4

namespace Z05.P1

open MeasureTheory Real Filter Topology Set

lemma logcube_le {x y : ℝ} (hy : 0 < y) (hyx : y ≤ x) (hx1 : x ≤ 1) :
    |Real.log x| ^ 3 ≤ 216 * y ^ (-(1 / 2 : ℝ)) := by
  have hx0 : 0 < x := hy.trans_le hyx
  have h := neg_rpow_div_le_log hx0.le (ε := 1 / 6) (by norm_num)
  have hl : Real.log x ≤ 0 := Real.log_nonpos hx0.le hx1
  have hp : 0 ≤ x ^ (-(1 / 6 : ℝ)) := Real.rpow_nonneg hx0.le _
  have hcube : (x ^ (-(1 / 6 : ℝ))) ^ 3 = x ^ (-(1 / 2 : ℝ)) := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hx0.le]; norm_num
  have hmono : x ^ (-(1 / 2 : ℝ)) ≤ y ^ (-(1 / 2 : ℝ)) :=
    Real.rpow_le_rpow_of_nonpos hy hyx (by norm_num)
  have h6 : |Real.log x| ≤ 6 * x ^ (-(1 / 6 : ℝ)) := by
    rw [abs_of_nonpos hl]
    have : -x ^ (-(1 / 6 : ℝ)) / (1 / 6) = -(6 * x ^ (-(1 / 6 : ℝ))) := by ring
    linarith
  have := pow_le_pow_left₀ (abs_nonneg _) h6 3
  rw [mul_pow, hcube] at this
  nlinarith

noncomputable def Phi3 (v : ℝ) : ℝ :=
  216 * (4 : ℝ) ^ (1 / 2 : ℝ) * ((v - Y01.b) ^ (-(1 / 2 : ℝ)) + (Y01.a - v) ^ (-(1 / 2 : ℝ)))

lemma ii_Phi3 : IntervalIntegrable Phi3 volume Y01.b Y01.a := by
  have h1 := (intervalIntegral.intervalIntegrable_rpow' (r := -(1 / 2 : ℝ)) (by norm_num)
    (a := 0) (b := Y01.a - Y01.b)).comp_sub_right Y01.b
  have h2 := (intervalIntegral.intervalIntegrable_rpow' (r := -(1 / 2 : ℝ)) (by norm_num)
    (a := Y01.a - Y01.b) (b := 0)).comp_sub_left Y01.a
  simp only [zero_add, sub_add_cancel, sub_sub_cancel, sub_zero] at h1 h2
  exact (h1.add h2).const_mul _

lemma l1cube_le {v : ℝ} (hv : v ∈ Ioo Y01.b Y01.a) :
    |Y01.l1 v| ^ 3 ≤ 216 * (4 : ℝ) ^ (1 / 2 : ℝ) * (v - Y01.b) ^ (-(1 / 2 : ℝ)) := by
  have h1 : 0 < v - Y01.b := by linarith [hv.1]
  have := logcube_le (by positivity) (Y01.PC.tv_ge hv) (Y01.PC.tv_lt_one hv).le
  rw [Real.div_rpow h1.le (by norm_num), Real.rpow_neg (by norm_num : (0:ℝ) ≤ 4),
    div_eq_mul_inv, inv_inv] at this
  unfold Y01.l1; linarith

lemma l2cube_le {v : ℝ} (hv : v ∈ Ioo Y01.b Y01.a) :
    |Y01.l2 v| ^ 3 ≤ 216 * (4 : ℝ) ^ (1 / 2 : ℝ) * (Y01.a - v) ^ (-(1 / 2 : ℝ)) := by
  have h1 : 0 < Y01.a - v := by linarith [hv.2]
  have := logcube_le (by positivity) (Y01.PC.one_sub_tv_ge hv) (by linarith [Y01.PC.tv_pos hv])
  rw [Real.div_rpow h1.le (by norm_num), Real.rpow_neg (by norm_num : (0:ℝ) ≤ 4),
    div_eq_mul_inv, inv_inv] at this
  unfold Y01.l2; linarith

lemma Phi3_ge {v : ℝ} (hv : v ∈ Ioo Y01.b Y01.a) :
    |Y01.l1 v| ^ 3 + |Y01.l2 v| ^ 3 ≤ Phi3 v := by
  have := l1cube_le hv; have := l2cube_le hv
  unfold Phi3; linarith

lemma cubic_ineq {x y : ℝ} (hx : 0 ≤ x) (hy : 0 ≤ y) : x ^ 2 * y ≤ x ^ 3 + y ^ 3 := by
  rcases le_total y x with h | h
  · nlinarith [mul_le_mul_of_nonneg_left h (sq_nonneg x), pow_nonneg hy 3]
  · nlinarith [mul_le_mul_of_nonneg_right (pow_le_pow_left₀ hx h 2) hy, pow_nonneg hx 3]

lemma ii_l1cube : IntervalIntegrable (fun v => Y01.l1 v ^ 3) volume Y01.b Y01.a :=
  Y01.PC.ii_of_bound Y01.PC.b_lt_a.le ii_Phi3 (Y01.PC.meas_l1.pow_const 3).aestronglyMeasurable
    (fun v hv => by
      have := Phi3_ge hv
      rw [Real.norm_eq_abs, abs_pow]
      linarith [pow_nonneg (abs_nonneg (Y01.l2 v)) 3])

lemma ii_l2cube : IntervalIntegrable (fun v => Y01.l2 v ^ 3) volume Y01.b Y01.a :=
  Y01.PC.ii_of_bound Y01.PC.b_lt_a.le ii_Phi3 (Y01.PC.meas_l2.pow_const 3).aestronglyMeasurable
    (fun v hv => by
      have := Phi3_ge hv
      rw [Real.norm_eq_abs, abs_pow]
      linarith [pow_nonneg (abs_nonneg (Y01.l1 v)) 3])

lemma ii_l1sql2 : IntervalIntegrable (fun v => Y01.l1 v ^ 2 * Y01.l2 v) volume Y01.b Y01.a :=
  Y01.PC.ii_of_bound Y01.PC.b_lt_a.le ii_Phi3
    ((Y01.PC.meas_l1.pow_const 2).mul Y01.PC.meas_l2).aestronglyMeasurable
    (fun v hv => by
      have := Phi3_ge hv
      rw [Real.norm_eq_abs, abs_mul, abs_pow]
      linarith [cubic_ineq (abs_nonneg (Y01.l1 v)) (abs_nonneg (Y01.l2 v))])

lemma ii_l1l2sq : IntervalIntegrable (fun v => Y01.l1 v * Y01.l2 v ^ 2) volume Y01.b Y01.a :=
  Y01.PC.ii_of_bound Y01.PC.b_lt_a.le ii_Phi3
    (Y01.PC.meas_l1.mul (Y01.PC.meas_l2.pow_const 2)).aestronglyMeasurable
    (fun v hv => by
      have := Phi3_ge hv
      rw [Real.norm_eq_abs, abs_mul, abs_pow, mul_comm]
      linarith [cubic_ineq (abs_nonneg (Y01.l2 v)) (abs_nonneg (Y01.l1 v))])

theorem main : stmt_1 := by
  have e : W4integrand = fun v => Y01.l1 v ^ 3 *
      ((8 - 7 * Real.sqrt 2) / (48 * v) - 1 / (3 * (v + Y01.a)))
      + Y01.l1 v ^ 2 * Y01.l2 v * (Real.sqrt 2 / (16 * v)) := by
    funext v; unfold W4integrand; ring
  unfold stmt_1; rw [e]
  have hA : ContinuousOn (fun v : ℝ => (8 - 7 * Real.sqrt 2) / (48 * v) - 1 / (3 * (v + Y01.a)))
      (uIcc Y01.b Y01.a) := by
    intro x hx
    have := Y01.PC.pos_of_mem hx; have := Y01.PC.a_pos
    apply ContinuousAt.continuousWithinAt
    have h1 : 48 * x ≠ 0 := by positivity
    have h2 : 3 * (x + Y01.a) ≠ 0 := by positivity
    fun_prop (disch := assumption)
  have hB : ContinuousOn (fun v : ℝ => Real.sqrt 2 / (16 * v)) (uIcc Y01.b Y01.a) := by
    intro x hx
    have := Y01.PC.pos_of_mem hx
    apply ContinuousAt.continuousWithinAt
    have h1 : 16 * x ≠ 0 := by positivity
    fun_prop (disch := assumption)
  exact (ii_l1cube.mul_continuousOn hA).add (ii_l1sql2.mul_continuousOn hB)

end Z05.P1

namespace Z05.P6

open MeasureTheory Real Filter Topology Set

/-- `g₂₀' + g₁₁'`, odd under the fold. -/
noncomputable def kk (v : ℝ) : ℝ := Y01.g20' v + g11d v

lemma kk_cont : Continuous kk := by
  show Continuous fun v => Real.sqrt 2 * v / (v ^ 2 + 1) ^ 2
    + Real.sqrt 2 * (v ^ 2 - 2 * v - 1) / (2 * (v ^ 2 + 1) ^ 2)
  refine Continuous.add ?_ ?_
  · exact Continuous.div (f := fun v => Real.sqrt 2 * v) (g := fun v => (v ^ 2 + 1) ^ 2)
      (by fun_prop) (by fun_prop) (fun x => by positivity)
  · exact Continuous.div (f := fun v => Real.sqrt 2 * (v ^ 2 - 2 * v - 1))
      (g := fun v => 2 * (v ^ 2 + 1) ^ 2) (by fun_prop) (by fun_prop) (fun x => by positivity)

noncomputable def ff (v : ℝ) : ℝ := Y01.l1 v * Y01.l2 v * kk v

lemma ii_ff : IntervalIntegrable ff volume Y01.b Y01.a :=
  Y01.PC.ii_l1l2.mul_continuousOn kk_cont.continuousOn

lemma l1_inv {u : ℝ} (hu : u ≠ 0) : Y01.l1 u⁻¹ = Y01.l2 u := by
  have := Y01.PC.l2_inv (inv_ne_zero hu)
  rw [inv_inv] at this; exact this.symm

lemma fold_pt {u : ℝ} (hu : u ∈ Ioo Y01.b Y01.a) : |-(u ^ 2)⁻¹| • ff u⁻¹ = -ff u := by
  have hu0 : 0 < u := Y01.PC.b_pos.trans hu.1
  unfold ff
  rw [l1_inv hu0.ne', Y01.PC.l2_inv hu0.ne', abs_neg, abs_of_pos (by positivity), smul_eq_mul]
  unfold kk Y01.g20' g11d
  field_simp
  ring

lemma int_ff : ∫ v in Y01.b..Y01.a, ff v = 0 := by
  have himg := integral_image_eq_integral_abs_deriv_smul (s := Ioo Y01.b Y01.a)
    (f := fun u : ℝ => u⁻¹) (f' := fun u => -(u ^ 2)⁻¹) measurableSet_Ioo
    (fun x hx => (hasDerivAt_inv (Y01.PC.b_pos.trans hx.1).ne').hasDerivWithinAt)
    inv_injective.injOn ff
  rw [Y01.PC.image_inv] at himg
  have e : ∫ v in Ioo Y01.b Y01.a, ff v = -∫ v in Ioo Y01.b Y01.a, ff v := by
    conv_lhs => rw [himg]
    rw [setIntegral_congr_fun measurableSet_Ioo (fun u hu => fold_pt hu), integral_neg]
  rw [intervalIntegral.integral_of_le Y01.PC.b_lt_a.le, integral_Ioc_eq_integral_Ioo]
  linarith

theorem main : stmt_6 := by
  have hpt : ∀ v, Y01.l1 v ^ 2 * Y01.g20' v + Y01.l1 v * Y01.l2 v * g11d v
      = Y01.PC.hfun v + ff v := by
    intro v; unfold Y01.PC.hfun ff kk; ring
  have hI : IntervalIntegrable
      (fun v => Y01.l1 v ^ 2 * Y01.g20' v + Y01.l1 v * Y01.l2 v * g11d v) volume Y01.b Y01.a := by
    simp_rw [hpt]; exact Y01.PC.ii_hfun.add ii_ff
  refine ⟨hI, ?_⟩
  simp_rw [hpt]
  rw [intervalIntegral.integral_add Y01.PC.ii_hfun ii_ff, Y01.PC.int_hfun, int_ff, add_zero]

end Z05.P6

namespace Z05.P10

open MeasureTheory Real Filter Topology Set

/-- The continuous extension of `Li₃(t)/t` (value `1` at `0`). -/
noncomputable def H4 (t : ℝ) : ℝ := if t = 0 then 1 else Y01.Li3 t / t

lemma H4_of_ne {t : ℝ} (h : t ≠ 0) : H4 t = Y01.Li3 t / t := by simp [H4, h]

lemma continuousAt_H4_zero : ContinuousAt H4 0 := by
  rw [← continuousWithinAt_compl_self, ContinuousWithinAt]
  have hd := Y01.hasDerivAt_Li3_H 0
  have hH : Y01.H 0 = 1 := by simp [Y01.H]
  rw [hH, hasDerivAt_iff_tendsto_slope] at hd
  have h0 : H4 0 = 1 := by simp [H4]
  rw [h0]
  apply hd.congr'
  filter_upwards [self_mem_nhdsWithin] with t ht
  rw [H4_of_ne ht, slope_def_field, Y01.Li3_zero]
  simp

lemma continuous_H4 : Continuous H4 := by
  rw [continuous_iff_continuousAt]
  intro t
  rcases eq_or_ne t 0 with rfl | h0
  · exact continuousAt_H4_zero
  · have hc : ContinuousAt (fun s => Y01.Li3 s / s) t :=
      (Y01.continuous_Li3.continuousAt).div continuousAt_id h0
    apply hc.congr
    filter_upwards [isOpen_ne.mem_nhds h0] with s hs
    exact (H4_of_ne hs).symm

lemma Li4_eq_H (x : ℝ) : Li4 x = ∫ t in (0:ℝ)..x, H4 t := by
  unfold Li4
  apply intervalIntegral.integral_congr_ae
  filter_upwards [(Set.countable_singleton (0:ℝ)).ae_notMem volume] with t ht _
  exact (H4_of_ne ht).symm

lemma hasDerivAt_Li4_H (x : ℝ) : HasDerivAt Li4 (H4 x) x := by
  have e : Li4 = fun x => ∫ t in (0:ℝ)..x, H4 t := funext Li4_eq_H
  rw [e]
  exact intervalIntegral.integral_hasDerivAt_right (continuous_H4.intervalIntegrable _ _)
    (continuous_H4.stronglyMeasurableAtFilter _ _) continuous_H4.continuousAt

lemma continuous_Li4 : Continuous Li4 :=
  continuous_iff_continuousAt.2 fun x => (hasDerivAt_Li4_H x).continuousAt

lemma hasDerivAt_Li4 {x : ℝ} (h0 : x ≠ 0) : HasDerivAt Li4 (Y01.Li3 x / x) x := by
  have := hasDerivAt_Li4_H x
  rwa [H4_of_ne h0] at this

lemma Li4_zero : Li4 0 = 0 := by simp [Li4]

/-- The power series `Σ x^(n+1)/(n+1)⁴`. -/
noncomputable def S4 (x : ℝ) : ℝ := ∑' n : ℕ, x ^ (n + 1) / ((n : ℝ) + 1) ^ 4

lemma summable_inv_four : Summable (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 4) := by
  have := (summable_nat_add_iff 1).mpr (Real.summable_one_div_nat_pow.mpr (by norm_num : 1 < 4))
  simpa using this

lemma term_bound4 {x : ℝ} (hx : |x| ≤ 1) (n : ℕ) :
    ‖x ^ (n + 1) / ((n : ℝ) + 1) ^ 4‖ ≤ 1 / ((n : ℝ) + 1) ^ 4 := by
  rw [Real.norm_eq_abs, abs_div, abs_pow, abs_of_pos (by positivity : (0:ℝ) < ((n:ℝ) + 1) ^ 4)]
  apply div_le_div_of_nonneg_right _ (by positivity)
  exact pow_le_one₀ (abs_nonneg x) hx

lemma summable_S4 {x : ℝ} (hx : |x| ≤ 1) :
    Summable (fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) ^ 4) :=
  Summable.of_norm_bounded summable_inv_four (term_bound4 hx)

lemma continuousOn_S4 : ContinuousOn S4 (Icc (-1) 1) := by
  apply continuousOn_tsum (fun n => (continuous_pow (n + 1)).continuousOn.div_const _)
    summable_inv_four
  intro n x hx
  exact term_bound4 (abs_le.mpr hx) n

lemma S4_zero : S4 0 = 0 := by simp [S4]

lemma hasDerivAt_S4 {y : ℝ} (hy : |y| < 1) (h0 : y ≠ 0) :
    HasDerivAt S4 (Y01.Li3 y / y) y := by
  set r := (|y| + 1) / 2 with hr
  have hr1 : r < 1 := by rw [hr]; linarith
  have hyr : |y| < r := by rw [hr]; linarith
  have hr0 : 0 ≤ r := by rw [hr]; positivity
  have ht : y ∈ Ioo (-r) r := ⟨by linarith [neg_abs_le y], by linarith [le_abs_self y]⟩
  have hd := hasDerivAt_tsum_of_isPreconnected (t := Ioo (-r) r) (y₀ := 0)
    (g := fun (n : ℕ) (z : ℝ) => z ^ (n + 1) / ((n : ℝ) + 1) ^ 4)
    (g' := fun (n : ℕ) (z : ℝ) => z ^ n / ((n : ℝ) + 1) ^ 3) (u := fun n : ℕ => r ^ n)
    (summable_geometric_of_lt_one hr0 hr1) isOpen_Ioo isPreconnected_Ioo
    (by
      intro n z _
      have := (hasDerivAt_pow (n + 1) z).div_const (((n : ℝ) + 1) ^ 4)
      refine this.congr_deriv ?_
      push_cast
      field_simp)
    (by
      intro n z hz
      rw [Real.norm_eq_abs, abs_div, abs_pow,
        abs_of_pos (by positivity : (0:ℝ) < ((n:ℝ) + 1) ^ 3)]
      have hz' : |z| ≤ r := abs_le.mpr ⟨hz.1.le, hz.2.le⟩
      calc |z| ^ n / ((n:ℝ) + 1) ^ 3 ≤ |z| ^ n := by
            apply div_le_self (by positivity)
            have : (0:ℝ) ≤ n := Nat.cast_nonneg n
            have : (1:ℝ) ≤ (n + 1) := by linarith
            exact one_le_pow₀ this
        _ ≤ r ^ n := pow_le_pow_left₀ (abs_nonneg z) hz' n)
    ⟨by linarith [abs_nonneg y, hr], by linarith [abs_nonneg y, hr]⟩
    (by simp)
    ht
  have hsum : HasSum (fun n : ℕ => y ^ n / ((n : ℝ) + 1) ^ 3) (Y01.Li3 y / y) := by
    have := (Y01.Li3_hasSum hy.le).div_const y
    convert this using 1
    funext n
    field_simp
    ring
  rw [hsum.tsum_eq] at hd
  exact hd

lemma Li4_eq_S4 {x : ℝ} (hx : |x| ≤ 1) : Li4 x = S4 x := by
  have key : (fun z => Li4 z - S4 z) 0 = (fun z => Li4 z - S4 z) x := by
    apply D00.eq_of_deriv_zero' (f := fun z => Li4 z - S4 z)
    · apply continuous_Li4.continuousOn.sub
      apply continuousOn_S4.mono
      intro z hz
      rw [mem_uIcc] at hz
      have hx' := abs_le.mp hx
      constructor <;> rcases hz with h | h <;> linarith [h.1, h.2]
    · intro z hz
      obtain ⟨hz0, hzx⟩ := D00.mem_uIoo_zero hz
      have := (hasDerivAt_Li4 hz0).sub (hasDerivAt_S4 (lt_of_lt_of_le hzx hx) hz0)
      rw [sub_self] at this
      exact this
  simp [S4_zero, Li4_zero] at key
  linarith

theorem main : stmt_10 := by
  refine ⟨fun x _ => ?_, fun x hx => ?_⟩
  · apply (continuous_H4.intervalIntegrable 0 x).congr_ae
    apply MeasureTheory.ae_restrict_of_ae
    filter_upwards [(Set.countable_singleton (0:ℝ)).ae_notMem MeasureTheory.volume] with t ht
    exact H4_of_ne ht
  · rw [Li4_eq_S4 hx]
    exact (summable_S4 hx).hasSum

end Z05.P10

namespace Z05.P2

open MeasureTheory Real Filter Topology Set

noncomputable def Lt (t : ℝ) : ℂ := Complex.log (t : ℂ)
noncomputable def Ls (t : ℝ) : ℂ := Complex.log ((1 - t : ℝ) : ℂ)
noncomputable def Lq (t : ℝ) : ℂ := Complex.log ((t / (1 - t) : ℝ) : ℂ)
noncomputable def Sx (e : ℂ) (t : ℝ) : ℂ := (1 + W02.zC e t) ^ (1 / 2 : ℂ)
noncomputable def Tx (e : ℂ) (t : ℝ) : ℂ := (1 + W02.zC e t) ^ ((1 / 2 : ℂ) - 1)
noncomputable def Ux (e : ℂ) (t : ℝ) : ℂ := (1 + W02.zC e t) ^ ((1 / 2 : ℂ) - 1 - 1)
noncomputable def Vx (e : ℂ) (t : ℝ) : ℂ := (1 + W02.zC e t) ^ ((1 / 2 : ℂ) - 1 - 1 - 1)
noncomputable def Wx (e : ℂ) (t : ℝ) : ℂ := W02.zC e t * Lq t * (-2)

/-- `Y01.PA.gCdd` in structured form. -/
noncomputable def gC2 (e : ℂ) (t : ℝ) : ℂ :=
  W02.aC e t * Lt t ^ 2 + W02.bC e t * Ls t ^ 2 - W02.aC e t * Lt t ^ 2 * Sx e t
    - W02.aC e t * Lt t * Tx e t * Wx e t + 1 / 4 * W02.aC e t * Ux e t * Wx e t ^ 2
    + Lq t * W02.aC e t * Tx e t * Wx e t

lemma gCdd_eq (e : ℂ) (t : ℝ) : Y01.PA.gCdd e t = gC2 e t := by
  unfold Y01.PA.gCdd gC2 Sx Tx Ux Wx Lt Ls Lq
  ring

/-- the shape of the third `ε`-derivative, in terms of its pieces. -/
noncomputable def F3 (A B x y l S T U V W : ℂ) : ℂ :=
  A * x ^ 3 + B * y ^ 3
    - (A * x ^ 3 * S + A * x ^ 2 * (1 / 2) * T * W)
    - (A * x ^ 2 * T * W + A * x * (-(1 / 2) * U * W) * W + A * x * T * (-2 * l * W))
    + 1 / 4 * (A * x * U * W ^ 2 + A * (-(3 / 2) * V * W) * W ^ 2 + A * U * (2 * W * (-2 * l * W)))
    + l * (A * x * T * W + A * (-(1 / 2) * U * W) * W + A * T * (-2 * l * W))

/-- the explicit `ε`-derivative of `gC2`. -/
noncomputable def gC3 (e : ℂ) (t : ℝ) : ℂ :=
  F3 (W02.aC e t) (W02.bC e t) (Lt t) (Ls t) (Lq t) (Sx e t) (Tx e t) (Ux e t) (Vx e t) (Wx e t)

lemma gC2_hasDerivAt {e : ℂ} (he : e.re < 1) {t : ℝ} (ht0 : 0 < t) (ht : t ≤ 1 / 2) :
    HasDerivAt (fun x => gC2 x t) (gC3 e t) e := by
  have h1 : 0 < 1 - t := by linarith
  have ha : HasDerivAt (fun x => W02.aC x t) (W02.aC e t * Lt t * 1) e :=
    ((hasDerivAt_id e).sub_const 1).const_cpow (Or.inl (by exact_mod_cast ht0.ne'))
  have hb : HasDerivAt (fun x => W02.bC x t) (W02.bC e t * Ls t * 1) e :=
    ((hasDerivAt_id e).sub_const 1).const_cpow (Or.inl (by exact_mod_cast h1.ne'))
  have hz : HasDerivAt (fun x => W02.zC x t) (W02.zC e t * Lq t * (0 - 2 * 1)) e :=
    ((hasDerivAt_const e (2 : ℂ)).sub ((hasDerivAt_id e).const_mul 2)).const_cpow
      (Or.inl (by exact_mod_cast (div_pos ht0 h1).ne'))
  have hS := ((hasDerivAt_const e (1 : ℂ)).add hz).cpow_const (c := (1 / 2 : ℂ))
    (W02.zC_slit he ht0 ht)
  have hT := ((hasDerivAt_const e (1 : ℂ)).add hz).cpow_const (c := (1 / 2 : ℂ) - 1)
    (W02.zC_slit he ht0 ht)
  have hU := ((hasDerivAt_const e (1 : ℂ)).add hz).cpow_const (c := (1 / 2 : ℂ) - 1 - 1)
    (W02.zC_slit he ht0 ht)
  have hW : HasDerivAt (fun x => W02.zC x t * Lq t * (-2))
      (W02.zC e t * Lq t * (0 - 2 * 1) * Lq t * (-2)) e :=
    (hz.mul_const (Lq t)).mul_const (-2)
  have h := (((((ha.mul_const (Lt t ^ 2)).add (hb.mul_const (Ls t ^ 2))).sub
      ((ha.mul_const (Lt t ^ 2)).mul hS)).sub (((ha.mul_const (Lt t)).mul hT).mul hW)).add
      (((ha.const_mul (1 / 4 : ℂ)).mul hU).mul (hW.pow 2))).add
      (((ha.const_mul (Lq t)).mul hT).mul hW)
  unfold gC2 Sx Tx Ux Wx
  convert h using 1
  unfold gC3 F3 Sx Tx Ux Vx Wx
  simp only [Pi.add_apply, Pi.mul_apply, Pi.sub_apply, Pi.pow_apply, Nat.cast_ofNat]
  ring

end Z05.P2

namespace Z05.P2

open MeasureTheory Real Filter Topology Set

lemma alg3 (t s Q x y : ℂ) (ht : t ≠ 0) (hs : s ≠ 0) (hQ : Q ≠ 0) (hQs : Q + s ≠ 0)
    (h1 : Q ^ 2 = t ^ 2 + s ^ 2) :
    F3 t⁻¹ s⁻¹ x y (x - y) (Q / s) (Q / s)⁻¹ ((Q / s) ^ 3)⁻¹ ((Q / s) ^ 5)⁻¹
        ((t / s) ^ 2 * (x - y) * (-2)) =
      6 * (1 / t * (1 / 6 - 2 * (s / Q) / 3 + (s / Q) ^ 3 - (s / Q) ^ 5 / 2) * x ^ 3
        + 1 / s * (1 / 6 - 2 * (t / Q) / 3 + (t / Q) ^ 3 - (t / Q) ^ 5 / 2) * y ^ 3
        + t * s / Q ^ 3 * ((1 - 3 * (s / Q) ^ 2 / 2) * x ^ 2 * y
          + (1 - 3 * (t / Q) ^ 2 / 2) * x * y ^ 2)) := by
  obtain ⟨u, hu⟩ : ∃ u, u = t / (Q + s) := ⟨_, rfl⟩
  have e1 : (1 + u ^ 2) * (Q + s) = 2 * Q := by
    rw [hu]; field_simp; linear_combination -h1
  have e2 : (1 - u ^ 2) * (Q + s) = 2 * s := by
    rw [hu]; field_simp; linear_combination h1
  have hu0 : u ≠ 0 := by rw [hu]; exact div_ne_zero ht hQs
  have hp : 1 + u ^ 2 ≠ 0 := by
    intro h; rw [h, zero_mul] at e1; exact hQ (by linear_combination -e1 / 2)
  have hm : 1 - u ^ 2 ≠ 0 := by
    intro h; rw [h, zero_mul] at e2; exact hs (by linear_combination -e2 / 2)
  have ht' : t = 2 * u * Q / (1 + u ^ 2) := by
    rw [eq_div_iff hp]
    have : t = u * (Q + s) := by rw [hu]; field_simp
    linear_combination (1 + u ^ 2) * this + u * e1
  have hs' : s = (1 - u ^ 2) * Q / (1 + u ^ 2) := by
    rw [eq_div_iff hp]
    linear_combination (1 - u ^ 2) / 2 * e1 - (1 + u ^ 2) / 2 * e2
  rw [ht', hs']
  unfold F3
  field_simp
  ring

lemma cpow_m52 (z : ℝ) (hz : 0 ≤ z) :
    (1 + (z : ℂ)) ^ ((1 / 2 : ℂ) - 1 - 1 - 1) = ((((Real.sqrt (1 + z)) ^ 5)⁻¹ : ℝ) : ℂ) := by
  have h1 : 0 ≤ 1 + z := by linarith
  rw [show (1 / 2 : ℂ) - 1 - 1 - 1 = ((-(5 / 2)) : ℝ) by push_cast; ring,
    ← Complex.ofReal_one, ← Complex.ofReal_add, ← Complex.ofReal_cpow h1]
  congr 1
  rw [Real.rpow_neg h1, Real.sqrt_eq_rpow, ← Real.rpow_natCast, ← Real.rpow_mul h1]
  norm_num

lemma gC3_zero {t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    gC3 0 t = ((6 * phi3 t : ℝ) : ℂ) := by
  have h1 : 0 < 1 - t := by linarith
  have hq := W02.q_pos t
  have hq2 := W02.q_sq t
  have hz0 : 0 ≤ (t / (1 - t)) ^ (2 - 2 * (0:ℝ)) :=
    Real.rpow_nonneg (div_nonneg ht0.le h1.le) _
  have ha : W02.aR 0 t = t⁻¹ := by
    unfold W02.aR; rw [show (0:ℝ) - 1 = -1 by ring, Real.rpow_neg_one]
  have hb : W02.bR 0 t = (1 - t)⁻¹ := by
    unfold W02.bR; rw [show (0:ℝ) - 1 = -1 by ring, Real.rpow_neg_one]
  have hz : (t / (1 - t)) ^ (2 - 2 * (0:ℝ)) = (t / (1 - t)) ^ 2 := by
    rw [show (2:ℝ) - 2 * 0 = ((2:ℕ):ℝ) by norm_num, Real.rpow_natCast]
  have hs : Real.sqrt (1 + (t / (1 - t)) ^ 2) = W02.q t / (1 - t) := by
    rw [show 1 + (t / (1 - t)) ^ 2 = (W02.q t / (1 - t)) ^ 2 by
      rw [div_pow, div_pow, hq2]; field_simp; ring]
    exact Real.sqrt_sq (by positivity)
  have hl : Real.log (t / (1 - t)) = Real.log t - Real.log (1 - t) :=
    Real.log_div ht0.ne' h1.ne'
  rw [← Complex.ofReal_zero]
  unfold gC3 Sx Tx Ux Vx Wx Lt Ls Lq
  rw [W02.aC_ofReal 0 ht0.le, W02.bC_ofReal 0 ht1.le, W02.zC_ofReal 0 ht0.le ht1,
    W02.sqrt_cpow _ hz0, W02.sqrt_inv_cpow _ hz0, Y01.PA.cpow_m32 _ hz0, cpow_m52 _ hz0,
    ← Complex.ofReal_log ht0.le, ← Complex.ofReal_log h1.le,
    ← Complex.ofReal_log (div_nonneg ht0.le h1.le), ha, hb, hz, hs, hl]
  have key := alg3 (t : ℂ) (1 - (t : ℂ)) (W02.q t : ℂ) (Real.log t : ℂ) (Real.log (1 - t) : ℂ)
    (by exact_mod_cast ht0.ne') (by exact_mod_cast h1.ne') (by exact_mod_cast hq.ne')
    (by exact_mod_cast (by linarith : W02.q t + (1 - t) ≠ 0))
    (by exact_mod_cast hq2)
  unfold phi3 Y01.c Y01.c'
  push_cast
  unfold F3 at key ⊢
  linear_combination key

end Z05.P2

namespace Z05.P2

open MeasureTheory Real Filter Topology Set

lemma gCdd_hasDerivAt {e : ℂ} (he : e.re < 1) {t : ℝ} (ht0 : 0 < t) (ht : t ≤ 1 / 2) :
    HasDerivAt (fun x => Y01.PA.gCdd x t) (gC3 e t) e := by
  have : (fun x => Y01.PA.gCdd x t) = fun x => gC2 x t := funext fun x => gCdd_eq x t
  rw [this]; exact gC2_hasDerivAt he ht0 ht

lemma gC3_bound {x : ℂ} (hx : x ∈ Metric.ball (0 : ℂ) (1 / 12)) {t : ℝ} (ht0 : 0 < t)
    (ht : t ≤ 1 / 2) : ‖gC3 x t‖ ≤ Y01.PA.Bd / (1 / 6) / (1 / 12) := by
  have hsub : Metric.closedBall x (1 / 12) ⊆ Metric.ball (0 : ℂ) (1 / 6) := by
    intro y hy
    rw [Metric.mem_closedBall, dist_eq_norm] at hy
    rw [mem_ball_zero_iff] at hx ⊢
    calc ‖y‖ = ‖(y - x) + x‖ := by ring_nf
      _ ≤ ‖y - x‖ + ‖x‖ := norm_add_le _ _
      _ < 1 / 6 := by linarith
  have h63 : Metric.ball (0 : ℂ) (1 / 6) ⊆ Metric.ball (0 : ℂ) (1 / 3) :=
    Metric.ball_subset_ball (by norm_num)
  have hdiff : DifferentiableOn ℂ (fun e => Y01.PA.gCdd e t) (Metric.ball (0 : ℂ) (1 / 6)) :=
    fun e he => (gCdd_hasDerivAt (Y01.PA.re_lt_of_ball (h63 he)) ht0
      ht).differentiableAt.differentiableWithinAt
  have hx6 : x ∈ Metric.ball (0 : ℂ) (1 / 6) := hsub (Metric.mem_closedBall_self (by norm_num))
  rw [← (gCdd_hasDerivAt (Y01.PA.re_lt_of_ball (h63 hx6)) ht0 ht).deriv]
  apply Complex.norm_deriv_le_of_forall_mem_sphere_norm_le (by norm_num)
    (hdiff.diffContOnCl_ball hsub)
  intro y hy
  exact Y01.PA.gCdd_bound (hsub (Metric.sphere_subset_closedBall hy)) ht0 ht

lemma gCd_iInt {x : ℂ} (hx : x ∈ Metric.ball (0 : ℂ) (1 / 3)) :
    IntervalIntegrable (W02.gCd x) volume 0 (1 / 2) := by
  refine IntervalIntegrable.mono_fun' (g := fun _ => Y01.PA.Bd)
    intervalIntegrable_const (W02.gCd_meas x).aestronglyMeasurable ?_
  rw [uIoc_of_le (by norm_num)]
  exact ae_restrict_of_forall_mem measurableSet_Ioc fun t ht =>
    Y01.PA.gCd_bound0 hx ht.1 ht.2

lemma key2x {x₀ : ℂ} (hx₀ : x₀ ∈ Metric.ball (0 : ℂ) (1 / 6)) :
    HasDerivAt (fun x => ∫ t in (0 : ℝ)..(1 / 2), W02.gCd x t)
      (∫ t in (0 : ℝ)..(1 / 2), Y01.PA.gCdd x₀ t) x₀ := by
  have h6 : Metric.ball (0 : ℂ) (1 / 6) ⊆ Metric.ball (0 : ℂ) (1 / 3) :=
    Metric.ball_subset_ball (by norm_num)
  exact (intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (μ := volume) (a := 0) (b := 1 / 2) (F := fun e t => W02.gCd e t) (F' := Y01.PA.gCdd)
    (x₀ := x₀) (s := Metric.ball (0 : ℂ) (1 / 6)) (bound := fun _ => Y01.PA.Bd / (1 / 6))
    (Metric.isOpen_ball.mem_nhds hx₀)
    (Filter.Eventually.of_forall fun x => (W02.gCd_meas x).aestronglyMeasurable)
    (gCd_iInt (h6 hx₀)) ((Y01.PA.gCdd_meas x₀).aestronglyMeasurable)
    (Filter.Eventually.of_forall fun t ht x hx => by
      rw [uIoc_of_le (by norm_num)] at ht
      exact Y01.PA.gCdd_bound hx ht.1 ht.2)
    intervalIntegrable_const
    (Filter.Eventually.of_forall fun t ht x hx => by
      rw [uIoc_of_le (by norm_num)] at ht
      exact Y01.PA.gCd_hasDerivAt (Y01.PA.re_lt_of_ball (h6 hx)) ht.1 ht.2)).2

lemma iter2_ev : iteratedDeriv 2 W02.FC =ᶠ[𝓝 0]
    fun x => 2 * ∫ t in (0 : ℝ)..(1 / 2), Y01.PA.gCdd x t := by
  have h6 : Metric.ball (0 : ℂ) (1 / 6) ⊆ Metric.ball (0 : ℂ) (1 / 3) :=
    Metric.ball_subset_ball (by norm_num)
  filter_upwards [Metric.ball_mem_nhds (0 : ℂ) (by norm_num : (0:ℝ) < 1 / 6)] with x hx
  have ho : IsOpen {e : ℂ | e.re < 1} := isOpen_lt Complex.continuous_re continuous_const
  have hev : deriv W02.FC =ᶠ[𝓝 x] fun y => 2 * ∫ t in (0 : ℝ)..(1 / 2), W02.gCd y t := by
    filter_upwards [ho.mem_nhds (Y01.PA.re_lt_of_ball (h6 hx))] with y hy
    exact (W02.FC_hasDerivAt hy).deriv
  rw [iteratedDeriv_succ, iteratedDeriv_one, hev.deriv_eq]
  exact ((key2x hx).const_mul (2 : ℂ)).deriv

lemma phi3_meas : Measurable phi3 := by
  unfold phi3 Y01.c Y01.c' W02.q
  fun_prop

lemma gC3_meas : AEStronglyMeasurable (gC3 0) (volume.restrict (Set.uIoc (0 : ℝ) (1 / 2))) := by
  have hm : AEStronglyMeasurable (fun t : ℝ => (((6 * phi3 t : ℝ)) : ℂ))
      (volume.restrict (Set.uIoc (0 : ℝ) (1 / 2))) :=
    (Complex.measurable_ofReal.comp (phi3_meas.const_mul 6)).aestronglyMeasurable
  refine hm.congr ?_
  rw [uIoc_of_le (by norm_num)]
  refine (ae_restrict_iff' measurableSet_Ioc).mpr (Filter.Eventually.of_forall fun t ht => ?_)
  exact (gC3_zero ht.1 (by linarith [ht.2])).symm

lemma key3 : HasDerivAt (fun x => ∫ t in (0 : ℝ)..(1 / 2), Y01.PA.gCdd x t)
    (∫ t in (0 : ℝ)..(1 / 2), gC3 0 t) 0 := by
  have h12 : Metric.ball (0 : ℂ) (1 / 12) ⊆ Metric.ball (0 : ℂ) (1 / 3) :=
    Metric.ball_subset_ball (by norm_num)
  exact (intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (μ := volume) (a := 0) (b := 1 / 2) (F := fun e t => Y01.PA.gCdd e t) (F' := gC3)
    (x₀ := 0) (s := Metric.ball (0 : ℂ) (1 / 12))
    (bound := fun _ => Y01.PA.Bd / (1 / 6) / (1 / 12))
    (Metric.ball_mem_nhds 0 (by norm_num))
    (Filter.Eventually.of_forall fun x => (Y01.PA.gCdd_meas x).aestronglyMeasurable)
    Y01.PA.key2.1 gC3_meas
    (Filter.Eventually.of_forall fun t ht x hx => by
      rw [uIoc_of_le (by norm_num)] at ht
      exact gC3_bound hx ht.1 ht.2)
    intervalIntegrable_const
    (Filter.Eventually.of_forall fun t ht x hx => by
      rw [uIoc_of_le (by norm_num)] at ht
      exact gCdd_hasDerivAt (Y01.PA.re_lt_of_ball (h12 hx)) ht.1 ht.2)).2

lemma iter3 : iteratedDeriv 3 W02.FC 0 = 2 * ∫ t in (0 : ℝ)..(1 / 2), gC3 0 t := by
  rw [iteratedDeriv_succ, iter2_ev.deriv_eq]
  exact (key3.const_mul (2 : ℂ)).deriv

lemma phi3_symm (t : ℝ) : phi3 (1 - t) = phi3 t := by
  unfold phi3
  rw [Y01.PB.c_one_sub, show Y01.c' (1 - t) = Y01.c t by
    unfold Y01.c Y01.c'; rw [W02.q_symm], W02.q_symm, sub_sub_cancel]
  ring

lemma phi3_bound_half {t : ℝ} (ht0 : 0 < t) (ht : t ≤ 1 / 2) :
    ‖phi3 t‖ ≤ Y01.PA.Bd / (1 / 6) / (1 / 12) / 6 := by
  have h := gC3_bound (x := 0) (by simp) ht0 ht
  rw [gC3_zero ht0 (by linarith), Complex.norm_real, Real.norm_eq_abs, abs_mul] at h
  rw [Real.norm_eq_abs]
  norm_num at h ⊢
  linarith

lemma phi3_bound {t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    ‖phi3 t‖ ≤ Y01.PA.Bd / (1 / 6) / (1 / 12) / 6 := by
  rcases le_total t (1 / 2) with h | h
  · exact phi3_bound_half ht0 h
  · rw [← phi3_symm t]
    exact phi3_bound_half (by linarith) (by linarith)

lemma phi3_iInt : IntervalIntegrable phi3 volume 0 1 := by
  refine IntervalIntegrable.mono_fun' (g := fun _ => Y01.PA.Bd / (1 / 6) / (1 / 12) / 6)
    intervalIntegrable_const phi3_meas.aestronglyMeasurable ?_
  rw [uIoc_of_le zero_le_one, ← restrict_Ioo_eq_restrict_Ioc]
  exact ae_restrict_of_forall_mem measurableSet_Ioo fun t ht => phi3_bound ht.1 ht.2

lemma phi3_half : ∫ t in (0:ℝ)..1, phi3 t = 2 * ∫ t in (0:ℝ)..(1/2), phi3 t := by
  have hi := phi3_iInt
  have h1 : IntervalIntegrable phi3 volume 0 (1/2) :=
    hi.mono_set (by rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
                    exact Icc_subset_Icc le_rfl (by norm_num))
  have h2 : IntervalIntegrable phi3 volume (1/2) 1 :=
    hi.mono_set (by rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
                    exact Icc_subset_Icc (by norm_num) le_rfl)
  rw [← intervalIntegral.integral_add_adjacent_intervals h1 h2]
  have : ∫ t in (1/2:ℝ)..1, phi3 t = ∫ t in (0:ℝ)..(1/2), phi3 t := by
    have := intervalIntegral.integral_comp_sub_left phi3 (a := 0) (b := 1/2) 1
    simp only [phi3_symm] at this
    rw [this]; norm_num
  rw [this]; ring

theorem main : stmt_2 := by
  refine ⟨phi3_iInt, ?_⟩
  unfold W02.Mk
  rw [iter3]
  have : ∫ t in (0:ℝ)..(1 / 2), gC3 0 t = ∫ t in (0:ℝ)..(1 / 2), (((6 * phi3 t : ℝ)) : ℂ) := by
    refine intervalIntegral.integral_congr_ae (Filter.Eventually.of_forall fun t ht => ?_)
    rw [uIoc_of_le (by norm_num)] at ht
    exact gC3_zero ht.1 (by linarith [ht.2])
  rw [this, intervalIntegral.integral_ofReal, intervalIntegral.integral_const_mul, phi3_half]
  simp only [Nat.factorial, Nat.cast_ofNat]
  push_cast
  ring

end Z05.P2

namespace Z05.P3

open MeasureTheory Real Filter Topology Set

/-- the `ℓ₁³` part of `φ₃`. -/
noncomputable def f3 (t : ℝ) : ℝ :=
  1 / t * (1 / 6 - 2 * Y01.c t / 3 + Y01.c t ^ 3 - Y01.c t ^ 5 / 2) * Real.log t ^ 3

/-- the `ℓ₁²ℓ₂` part of `φ₃`. -/
noncomputable def g3 (t : ℝ) : ℝ :=
  t * (1 - t) / W02.q t ^ 3 * ((1 - 3 * Y01.c t ^ 2 / 2) * Real.log t ^ 2 * Real.log (1 - t))

lemma phi3_eq (t : ℝ) : phi3 t = f3 t + f3 (1 - t) + g3 t + g3 (1 - t) := by
  unfold phi3 f3 g3
  rw [Y01.PB.c_one_sub, W02.q_symm, sub_sub_cancel]
  ring

/-- `t log³ t` written through `x = t^{1/3}`. -/
lemma tlog3_eq {t : ℝ} (ht : 0 ≤ t) :
    t * Real.log t ^ 3 = 27 * (t ^ (1 / 3 : ℝ) * Real.log (t ^ (1 / 3 : ℝ))) ^ 3 := by
  rcases eq_or_lt_of_le ht with h | h
  · subst h; simp
  have h3 : (t ^ (1 / 3 : ℝ)) ^ 3 = t := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul h.le]; norm_num
  rw [Real.log_rpow h, mul_pow, h3]
  ring

lemma cont_tlog3 : ContinuousOn (fun t : ℝ => t * Real.log t ^ 3) (Icc 0 1) := by
  have hc : Continuous fun t : ℝ => 27 * (t ^ (1 / 3 : ℝ) * Real.log (t ^ (1 / 3 : ℝ))) ^ 3 :=
    ((Real.continuous_mul_log.comp (Real.continuous_rpow_const (by norm_num))).pow 3).const_mul 27
  exact hc.continuousOn.congr (fun t ht => tlog3_eq ht.1)

lemma tlog2_eq {t : ℝ} (ht : 0 ≤ t) :
    t * Real.log t ^ 2 = 4 * (Real.sqrt t * Real.log (Real.sqrt t)) ^ 2 := by
  rw [Real.log_sqrt ht, mul_pow, Real.sq_sqrt ht]; ring

lemma cont_tlog2 : ContinuousOn (fun t : ℝ => t * Real.log t ^ 2) (Icc 0 1) := by
  have hc : Continuous fun t : ℝ => 4 * (Real.sqrt t * Real.log (Real.sqrt t)) ^ 2 :=
    ((Real.continuous_mul_log.comp Real.continuous_sqrt).pow 2).const_mul 4
  exact hc.continuousOn.congr (fun t ht => tlog2_eq ht.1)

lemma f3_eq {t : ℝ} (ht : 0 ≤ t) :
    f3 t = -(t * Real.log t ^ 3) *
      (-(Y01.c t ^ 4) / 2 - Y01.c t ^ 3 / 2 + Y01.c t ^ 2 / 2 + Y01.c t / 2 - 1 / 6)
      / (W02.q t * (1 - t + W02.q t)) := by
  rcases eq_or_lt_of_le ht with h | h
  · subst h; simp [f3]
  have hq := W02.q_pos t
  have hd := Y01.PB.den_pos t
  have h1 : 1 / 6 - 2 * Y01.c t / 3 + Y01.c t ^ 3 - Y01.c t ^ 5 / 2
      = (Y01.c t - 1) * (-(Y01.c t ^ 4) / 2 - Y01.c t ^ 3 / 2 + Y01.c t ^ 2 / 2 + Y01.c t / 2
        - 1 / 6) := by ring
  unfold f3
  rw [h1, Y01.PB.c_sub_one]
  field_simp

lemma f3_cont : ContinuousOn f3 (Icc 0 1) := by
  have hc := Y01.PB.c_cont
  have hq := Y01.PB.q_cont
  have h : ContinuousOn (fun t => -(t * Real.log t ^ 3) *
      (-(Y01.c t ^ 4) / 2 - Y01.c t ^ 3 / 2 + Y01.c t ^ 2 / 2 + Y01.c t / 2 - 1 / 6)
      / (W02.q t * (1 - t + W02.q t))) (Icc 0 1) := by
    refine ContinuousOn.div ?_ ?_ ?_
    · exact cont_tlog3.neg.mul (by fun_prop : Continuous fun t => -(Y01.c t ^ 4) / 2
        - Y01.c t ^ 3 / 2 + Y01.c t ^ 2 / 2 + Y01.c t / 2 - 1 / 6).continuousOn
    · exact (by fun_prop : Continuous fun t => W02.q t * (1 - t + W02.q t)).continuousOn
    · intro t _
      exact (mul_pos (W02.q_pos t) (Y01.PB.den_pos t)).ne'
  exact h.congr (fun t ht => f3_eq ht.1)

lemma g3_eq (t : ℝ) : g3 t = (t * Real.log t ^ 2) * ((1 - t) * Real.log (1 - t)) *
    (1 - 3 * Y01.c t ^ 2 / 2) / W02.q t ^ 3 := by
  unfold g3; ring

lemma g3_cont : ContinuousOn g3 (Icc 0 1) := by
  have hc := Y01.PB.c_cont
  have hq := Y01.PB.q_cont
  have h : ContinuousOn (fun t => (t * Real.log t ^ 2) * ((1 - t) * Real.log (1 - t)) *
      (1 - 3 * Y01.c t ^ 2 / 2) / W02.q t ^ 3) (Icc 0 1) := by
    refine ContinuousOn.div ?_ ?_ ?_
    · refine (cont_tlog2.mul ?_).mul ?_
      · exact (Real.continuous_mul_log.comp (continuous_const.sub continuous_id)).continuousOn
      · exact (by fun_prop : Continuous fun t => 1 - 3 * Y01.c t ^ 2 / 2).continuousOn
    · exact (hq.pow 3).continuousOn
    · intro t _
      exact (pow_pos (W02.q_pos t) 3).ne'
  exact h.congr (fun t ht => g3_eq t)

lemma comp_one_sub {F : ℝ → ℝ} (hF : ContinuousOn F (Icc 0 1)) :
    ContinuousOn (fun t => F (1 - t)) (Icc 0 1) := by
  refine hF.comp (continuousOn_const.sub continuousOn_id) ?_
  intro x hx
  exact ⟨by linarith [hx.2], by linarith [hx.1]⟩

lemma iint {F : ℝ → ℝ} (hF : ContinuousOn F (Icc 0 1)) : IntervalIntegrable F volume 0 1 :=
  (hF.mono (by rw [uIcc_of_le zero_le_one])).intervalIntegrable

lemma int_one_sub (F : ℝ → ℝ) : ∫ t in (0 : ℝ)..1, F (1 - t) = ∫ t in (0 : ℝ)..1, F t := by
  rw [intervalIntegral.integral_comp_sub_left F 1]; norm_num

/-- the folded integrand in `t`. -/
noncomputable def G3 (t : ℝ) : ℝ := 2 * f3 t + 2 * g3 t

lemma G3_cont : ContinuousOn G3 (Icc 0 1) :=
  (continuousOn_const.mul f3_cont).add (continuousOn_const.mul g3_cont)

lemma phi3_int_eq : ∫ t in (0 : ℝ)..1, phi3 t = ∫ t in (0 : ℝ)..1, G3 t := by
  have e1 : (fun t => phi3 t) = fun t => f3 t + f3 (1 - t) + g3 t + g3 (1 - t) := funext phi3_eq
  have i1 := iint f3_cont
  have i2 := iint (comp_one_sub f3_cont)
  have i3 := iint g3_cont
  have i4 := iint (comp_one_sub g3_cont)
  rw [e1, intervalIntegral.integral_add ((i1.add i2).add i3) i4,
    intervalIntegral.integral_add (i1.add i2) i3, intervalIntegral.integral_add i1 i2,
    int_one_sub f3, int_one_sub g3]
  unfold G3
  rw [intervalIntegral.integral_add (i1.const_mul 2) (i3.const_mul 2),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul]
  ring

end Z05.P3

namespace Z05.P3

open MeasureTheory Real Filter Topology Set

lemma R30_eq {v : ℝ} (hv : 0 < v) (hb : Y01.b < v) :
    R30 v = -(1 / Y01.tv v * (1 / 6 - 2 * Y01.c (Y01.tv v) / 3 + Y01.c (Y01.tv v) ^ 3
      - Y01.c (Y01.tv v) ^ 5 / 2)) * Y01.PB.tv' v := by
  have hs := Y01.PB.s_sq
  have h1 := Y01.PB.tpoly_pos hs Y01.PB.hr1 hb
  have h2 : 0 < v + (1 + Real.sqrt 2) := by linarith [Y01.PB.s_pos]
  have h3 : 0 < v ^ 2 + 1 := by positivity
  have h1' : v * (v + 2) - 1 ≠ 0 := by nlinarith
  rw [Y01.PB.c_tv hv]
  unfold R30 G30d Y01.tv Y01.PB.tv' Y01.a
  field_simp
  ring_nf
  simp only [Y01.PE1.s2_3, Y01.PE1.s2_4, Y01.PE1.s2_5, Y01.PE1.s2_6, Y01.PE1.s2_7, Y01.PE1.s2_8,
    Y01.PE1.s2_9, Y01.PE1.s2_10, Y01.PE1.s2_11, Y01.PE1.s2_12, hs]
  ring_nf

lemma R21_eq {v : ℝ} (hv : 0 < v) :
    R21 v = -(Y01.tv v * (1 - Y01.tv v) / W02.q (Y01.tv v) ^ 3
      * (1 - 3 * Y01.c (Y01.tv v) ^ 2 / 2)) * Y01.PB.tv' v := by
  have hs := Y01.PB.s_sq
  have h3 : 0 < v ^ 2 + 1 := by positivity
  have h4 := Y01.PB.s_pos
  rw [Y01.PB.c_tv hv, Y01.PB.q_tv hv]
  unfold R21 G21d Y01.tv Y01.PB.tv'
  field_simp
  ring_nf
  simp only [Y01.PE1.s2_3, Y01.PE1.s2_4, Y01.PE1.s2_5, Y01.PE1.s2_6, Y01.PE1.s2_7, Y01.PE1.s2_8,
    Y01.PE1.s2_9, Y01.PE1.s2_10, Y01.PE1.s2_11, Y01.PE1.s2_12, hs]
  ring_nf

lemma F_eq_inner {v : ℝ} (hv : 0 < v) (hb : Y01.b < v) :
    Y01.l1 v ^ 3 * R30 v + Y01.l1 v ^ 2 * Y01.l2 v * R21 v
      = -(1 / 2) * (G3 (Y01.tv v) * Y01.PB.tv' v) := by
  rw [R30_eq hv hb, R21_eq hv]
  unfold Y01.l1 Y01.l2 G3 f3 g3
  ring

lemma F_eq {v : ℝ} (hv : v ∈ Icc Y01.b Y01.a) :
    Y01.l1 v ^ 3 * R30 v + Y01.l1 v ^ 2 * Y01.l2 v * R21 v
      = -(1 / 2) * (G3 (Y01.tv v) * Y01.PB.tv' v) := by
  rcases eq_or_lt_of_le hv.1 with h | h
  · subst h
    simp [Y01.l1, Y01.l2, Y01.PB.tv_b, G3, f3, g3]
  rcases eq_or_lt_of_le hv.2 with h' | h'
  · subst h'
    simp [Y01.l1, Y01.l2, Y01.PB.tv_a, G3, f3, g3]
  exact F_eq_inner (Y01.PB.b_pos.trans h) h

lemma inner_cont : ContinuousOn (fun v => -(1 / 2) * (G3 (Y01.tv v) * Y01.PB.tv' v))
    (Icc Y01.b Y01.a) :=
  continuousOn_const.mul ((G3_cont.comp Y01.PB.tv_cont (fun _ hx => Y01.PB.tv_mem hx)).mul
    Y01.PB.tv'_cont)

lemma subst : ∫ v in Y01.b..Y01.a, G3 (Y01.tv v) * Y01.PB.tv' v = ∫ t in (0 : ℝ)..1, G3 t := by
  have h : ∀ x ∈ uIcc Y01.b Y01.a, HasDerivAt Y01.tv (Y01.PB.tv' x) x := by
    intro x hx
    rw [uIcc_of_le Y01.PB.b_le_a] at hx
    exact Y01.PB.tv_hasDerivAt (Y01.PB.b_pos.trans_le hx.1).ne'
  have h' : ContinuousOn Y01.PB.tv' (uIcc Y01.b Y01.a) := by
    rw [uIcc_of_le Y01.PB.b_le_a]; exact Y01.PB.tv'_cont
  have hg : ContinuousOn G3 (Y01.tv '' uIcc Y01.b Y01.a) := by
    apply G3_cont.mono
    rintro _ ⟨v, hv, rfl⟩
    rw [uIcc_of_le Y01.PB.b_le_a] at hv
    exact Y01.PB.tv_mem hv
  have := intervalIntegral.integral_comp_mul_deriv' h h' hg
  rw [Y01.PB.tv_b, Y01.PB.tv_a] at this
  exact this

lemma F_int_eq : ∫ v in Y01.b..Y01.a, (Y01.l1 v ^ 3 * R30 v + Y01.l1 v ^ 2 * Y01.l2 v * R21 v)
    = -(1 / 2) * ∫ t in (0 : ℝ)..1, phi3 t := by
  rw [phi3_int_eq, ← subst, ← intervalIntegral.integral_const_mul]
  apply intervalIntegral.integral_congr
  intro v hv
  rw [uIcc_of_le Y01.PB.b_le_a] at hv
  exact F_eq hv

lemma F_intble : IntervalIntegrable
    (fun v => Y01.l1 v ^ 3 * R30 v + Y01.l1 v ^ 2 * Y01.l2 v * R21 v) volume Y01.b Y01.a := by
  apply ContinuousOn.intervalIntegrable
  rw [uIcc_of_le Y01.PB.b_le_a]
  exact inner_cont.congr (fun v hv => F_eq hv)

theorem main : stmt_2 → stmt_3 := by
  rintro ⟨_, hM⟩
  refine ⟨F_intble, ?_⟩
  rw [F_int_eq]
  exact hM

end Z05.P3

namespace Z05.P5

open MeasureTheory Real Filter Topology Set

/-- the integrand of the degree-two remainder in stmt_5. -/
noncomputable def Tg (v : ℝ) : ℝ :=
  Y01.l1 v ^ 2 * (Om20 v + Y01.g20' v) + Y01.l1 v * Y01.l2 v * (Om11 v + g11d v)

/-- `Ω₂₀ + g₂₀'` without the pole at `a`. -/
noncomputable def Om20r (v : ℝ) : ℝ :=
  (1 - Real.sqrt 2) / (4 * v) - 1 / (4 * (v + Y01.a)) - 1 / (2 * (v + Y01.b)) + Y01.g20' v

lemma Tg_split (v : ℝ) : Tg v = Y01.l1 v ^ 2 * Om20r v
    + Y01.l1 v * (Y01.l1 v * (1 / (2 * (v - Y01.a)))) * (1 / 2)
    + Y01.l1 v * Y01.l2 v * (Om11 v + g11d v) := by
  unfold Tg Om20r Om20
  rw [show 1 / (4 * (v - Y01.a)) = (1 / (2 * (v - Y01.a))) * (1 / 2) by
    rw [div_mul_div_comm, one_mul]; ring_nf]
  ring

lemma contOn_Om20r : ContinuousOn Om20r (uIcc Y01.b Y01.a) := by
  intro x hx
  have := Y01.PC.pos_of_mem hx; have := Y01.PC.a_pos; have := Y01.PC.b_pos
  apply ContinuousAt.continuousWithinAt
  have h1 : 4 * x ≠ 0 := by positivity
  have h2 : 4 * (x + Y01.a) ≠ 0 := by positivity
  have h3 : 2 * (x + Y01.b) ≠ 0 := by positivity
  have h4 : (x ^ 2 + 1) ^ 2 ≠ 0 := by positivity
  unfold Om20r Y01.g20'
  fun_prop (disch := assumption)

lemma contOn_Om11 : ContinuousOn (fun v => Om11 v + g11d v) (uIcc Y01.b Y01.a) := by
  intro x hx
  have := Y01.PC.pos_of_mem hx; have := Y01.PC.a_pos
  apply ContinuousAt.continuousWithinAt
  have h1 : 4 * x ≠ 0 := by positivity
  have h2 : 2 * (x + Y01.a) ≠ 0 := by positivity
  have h4 : 2 * (x ^ 2 + 1) ^ 2 ≠ 0 := by positivity
  unfold Om11 g11d
  fun_prop (disch := assumption)

lemma ii_pole : IntervalIntegrable
    (fun v => Y01.l1 v * (Y01.l1 v * (1 / (2 * (v - Y01.a))))) volume Y01.b Y01.a := by
  have hsub1 : uIcc Y01.b 1 ⊆ uIcc Y01.b Y01.a := by
    rw [uIcc_of_le Y01.PC.b_lt_one.le, uIcc_of_le Y01.PC.b_lt_a.le]
    exact Icc_subset_Icc le_rfl Y01.PC.one_lt_a.le
  have hsub2 : uIcc 1 Y01.a ⊆ uIcc Y01.b Y01.a := by
    rw [uIcc_of_le Y01.PC.one_lt_a.le, uIcc_of_le Y01.PC.b_lt_a.le]
    exact Icc_subset_Icc Y01.PC.b_lt_one.le le_rfl
  have i1 : IntervalIntegrable (fun v => Y01.l1 v * (Y01.l1 v * (1 / (2 * (v - Y01.a)))))
      volume Y01.b 1 := by
    have e : (fun v => Y01.l1 v * (Y01.l1 v * (1 / (2 * (v - Y01.a)))))
        = fun v => Y01.l1 v ^ 2 * (1 / (2 * (v - Y01.a))) := by funext v; ring
    rw [e]
    refine (Y01.PC.ii_l1sq.mono_set hsub1).mul_continuousOn ?_
    rw [uIcc_of_le Y01.PC.b_lt_one.le]
    intro x hx
    have hxa : 2 * (x - Y01.a) ≠ 0 := by
      have := Y01.PC.one_lt_a; have := hx.2
      intro h; linarith
    apply ContinuousAt.continuousWithinAt
    fun_prop (disch := assumption)
  have i2 : IntervalIntegrable (fun v => Y01.l1 v * (Y01.l1 v * (1 / (2 * (v - Y01.a)))))
      volume 1 Y01.a :=
    Y01.PC.ii_of_bound Y01.PC.one_lt_a.le
      ((Y01.PC.ii_l1.mono_set hsub2).norm)
      (Y01.PC.meas_l1.mul (Y01.PC.meas_l1.mul (by fun_prop))).aestronglyMeasurable
      (fun v hv => by
        have hb := Y01.PC.l1_pole_bound hv
        rw [norm_mul]
        calc ‖Y01.l1 v‖ * ‖Y01.l1 v * (1 / (2 * (v - Y01.a)))‖ ≤ ‖Y01.l1 v‖ * 1 :=
              mul_le_mul_of_nonneg_left hb (norm_nonneg _)
          _ = ‖Y01.l1 v‖ := mul_one _)
  exact i1.trans i2

lemma ii_Tg : IntervalIntegrable Tg volume Y01.b Y01.a := by
  have e : Tg = fun v => Y01.l1 v ^ 2 * Om20r v
      + Y01.l1 v * (Y01.l1 v * (1 / (2 * (v - Y01.a)))) * (1 / 2)
      + Y01.l1 v * Y01.l2 v * (Om11 v + g11d v) := funext Tg_split
  rw [e]
  exact ((Y01.PC.ii_l1sq.mul_continuousOn contOn_Om20r).add (ii_pole.mul_const _)).add
    (Y01.PC.ii_l1l2.mul_continuousOn contOn_Om11)

lemma cont_G30d : Continuous G30d := by
  show Continuous fun v => -(Real.sqrt 2 * (v ^ 6 - 2 * v ^ 5 - 7 * v ^ 4 - 4 * v ^ 3 + 7 * v ^ 2
    - 2 * v - 1)) / (4 * (v ^ 2 + 1) ^ 4)
  exact Continuous.div (f := fun v => -(Real.sqrt 2 * (v ^ 6 - 2 * v ^ 5 - 7 * v ^ 4 - 4 * v ^ 3
    + 7 * v ^ 2 - 2 * v - 1))) (g := fun v => 4 * (v ^ 2 + 1) ^ 4)
    (by fun_prop) (by fun_prop) (fun x => by positivity)

lemma cont_G21d : Continuous G21d := by
  show Continuous fun v => Real.sqrt 2 * (3 * v ^ 6 - 2 * v ^ 5 - 21 * v ^ 4 - 4 * v ^ 3
    + 21 * v ^ 2 - 2 * v - 3) / (4 * (v ^ 2 + 1) ^ 4)
  exact Continuous.div (f := fun v => Real.sqrt 2 * (3 * v ^ 6 - 2 * v ^ 5 - 21 * v ^ 4
    - 4 * v ^ 3 + 21 * v ^ 2 - 2 * v - 3)) (g := fun v => 4 * (v ^ 2 + 1) ^ 4)
    (by fun_prop) (by fun_prop) (fun x => by positivity)

/-- the polynomial-kernel part of `M_3`. -/
noncomputable def Gp (v : ℝ) : ℝ := Y01.l1 v ^ 3 * G30d v + Y01.l1 v ^ 2 * Y01.l2 v * G21d v

lemma ii_Gp : IntervalIntegrable Gp volume Y01.b Y01.a :=
  (Z05.P1.ii_l1cube.mul_continuousOn cont_G30d.continuousOn).add
    (Z05.P1.ii_l1sql2.mul_continuousOn cont_G21d.continuousOn)

/-- the primitive `H = G₃₀ℓ₁³ + G₂₁ℓ₁²ℓ₂`. -/
noncomputable def H (v : ℝ) : ℝ := G30 v * Y01.l1 v ^ 3 + G21 v * (Y01.l1 v * (Y01.l1 v * Y01.l2 v))

lemma hasDerivAt_H {v : ℝ} (hv : v ∈ Ioo Y01.b Y01.a) : HasDerivAt H (Gp v - Tg v) v := by
  obtain ⟨hd, -, -, -, -, hrest⟩ := Z05.P4.main
  obtain ⟨d30, d21, -, -⟩ := hd v
  obtain ⟨dl1, dl2, i1, i2⟩ := hrest v hv.1 hv.2
  have h := (d30.mul (dl1.pow 3)).add (d21.mul (dl1.mul (dl1.mul dl2)))
  refine h.congr_deriv ?_
  unfold Gp Tg
  simp only [Pi.add_apply, Pi.mul_apply, Pi.pow_apply, Nat.cast_ofNat]
  linear_combination (-(Y01.l1 v ^ 2)) * i1 - (Y01.l1 v * Y01.l2 v) * i2

lemma H_b : H Y01.b = 0 := by
  unfold H; rw [Z05.P4.vals.1, Z05.P4.vals.2.1]; ring

lemma H_a : H Y01.a = 0 := by
  unfold H; rw [Y01.PC.l1_a]; ring

/-- a function with derivative bounded on `[b,a]` and vanishing at `b` is `O(t(v))`. -/
lemma lin_bound {F F' : ℝ → ℝ} (hF : ∀ v, HasDerivAt F (F' v) v) (hF' : Continuous F')
    (hb : F Y01.b = 0) : ∃ C, 0 ≤ C ∧ ∀ v ∈ Ioo Y01.b Y01.a, |F v| ≤ C * Y01.tv v := by
  obtain ⟨C, hC⟩ := isCompact_Icc.exists_bound_of_continuousOn
    (hF'.continuousOn (s := Icc Y01.b Y01.a))
  refine ⟨4 * max C 0, by positivity, fun v hv => ?_⟩
  have hvI : v ∈ Icc Y01.b Y01.a := ⟨hv.1.le, hv.2.le⟩
  have hbI : Y01.b ∈ Icc Y01.b Y01.a := ⟨le_rfl, Y01.PC.b_lt_a.le⟩
  have := (convex_Icc Y01.b Y01.a).norm_image_sub_le_of_norm_hasDerivWithin_le
    (f := F) (f' := F') (C := max C 0) (fun x _ => (hF x).hasDerivWithinAt)
    (fun x hx => (hC x hx).trans (le_max_left _ _)) hbI hvI
  rw [hb, sub_zero, Real.norm_eq_abs, Real.norm_eq_abs,
    abs_of_pos (show (0:ℝ) < v - Y01.b by linarith [hv.1])] at this
  have ht := Y01.PC.tv_ge hv
  have hm : 0 ≤ max C 0 := le_max_right _ _
  nlinarith

lemma contAt_tv_b : ContinuousAt Y01.tv Y01.b := Y01.PC.contAt_tv Y01.PC.b_pos.ne'

lemma tendsto_K {K : ℝ → ℝ} (hK : ContinuousOn K (Icc 0 1)) (hK0 : K 0 = 0) :
    Tendsto (fun v => K (Y01.tv v)) (𝓝[>] Y01.b) (𝓝 0) := by
  have hT : Tendsto Y01.tv (𝓝[>] Y01.b) (𝓝[Icc 0 1] 0) := by
    rw [tendsto_nhdsWithin_iff]
    refine ⟨?_, ?_⟩
    · have := contAt_tv_b.tendsto
      rw [Y01.PC.tv_b] at this
      exact this.mono_left nhdsWithin_le_nhds
    · filter_upwards [Ioo_mem_nhdsGT Y01.PC.b_lt_a] with v hv using
        ⟨(Y01.PC.tv_pos hv).le, (Y01.PC.tv_lt_one hv).le⟩
  have := (hK 0 ⟨le_rfl, zero_le_one⟩).tendsto.comp hT
  rwa [hK0] at this

lemma tendsto_H_b : Tendsto H (𝓝[>] Y01.b) (𝓝 (H Y01.b)) := by
  obtain ⟨hd, -⟩ := Z05.P4.main
  obtain ⟨C1, hC1, h1⟩ := lin_bound (fun v => (hd v).1) cont_G30d Z05.P4.vals.1
  obtain ⟨C2, hC2, h2⟩ := lin_bound (fun v => (hd v).2.1) cont_G21d Z05.P4.vals.2.1
  have k3 := tendsto_K Z05.P3.cont_tlog3 (by simp)
  have k2 := tendsto_K Z05.P3.cont_tlog2 (by simp)
  have hl2 := Y01.PC.contAt_l2_b.tendsto.mono_left (nhdsWithin_le_nhds (s := Ioi Y01.b))
  rw [Y01.PC.l2_b] at hl2
  rw [H_b]
  have hbd : Tendsto (fun v => C1 * |Y01.tv v * Real.log (Y01.tv v) ^ 3|
      + C2 * (|Y01.tv v * Real.log (Y01.tv v) ^ 2| * |Y01.l2 v|)) (𝓝[>] Y01.b) (𝓝 0) := by
    have := ((k3.abs.const_mul C1).add ((k2.abs.mul hl2.abs).const_mul C2))
    simpa using this
  refine squeeze_zero_norm' ?_ hbd
  filter_upwards [Ioo_mem_nhdsGT Y01.PC.b_lt_a] with v hv
  have ht := Y01.PC.tv_pos hv
  have e1 := h1 v hv
  have e2 := h2 v hv
  unfold H Y01.l1
  rw [Real.norm_eq_abs]
  calc |G30 v * Real.log (Y01.tv v) ^ 3 + G21 v * (Real.log (Y01.tv v) *
        (Real.log (Y01.tv v) * Y01.l2 v))|
      ≤ |G30 v| * |Real.log (Y01.tv v)| ^ 3
        + |G21 v| * (|Real.log (Y01.tv v)| ^ 2 * |Y01.l2 v|) := by
        refine (abs_add_le _ _).trans (le_of_eq ?_)
        rw [abs_mul, abs_mul, abs_pow, abs_mul, abs_mul]; ring
    _ ≤ (C1 * Y01.tv v) * |Real.log (Y01.tv v)| ^ 3
        + (C2 * Y01.tv v) * (|Real.log (Y01.tv v)| ^ 2 * |Y01.l2 v|) := by
        gcongr
    _ = C1 * |Y01.tv v * Real.log (Y01.tv v) ^ 3|
        + C2 * (|Y01.tv v * Real.log (Y01.tv v) ^ 2| * |Y01.l2 v|) := by
        rw [abs_mul, abs_mul, abs_of_pos ht, abs_pow, abs_pow]; ring

lemma tendsto_l1l2_a : Tendsto (fun v => Y01.l1 v * Y01.l2 v) (𝓝[<] Y01.a) (𝓝 0) := by
  have hS : Tendsto (fun v => (1 - Y01.tv v) * Real.log (1 - Y01.tv v)) (𝓝[<] Y01.a) (𝓝 0) := by
    have hc : Tendsto (fun v => 1 - Y01.tv v) (𝓝 Y01.a) (𝓝 0) := by
      have h : Tendsto (fun v => 1 - Y01.tv v) (𝓝 Y01.a) (𝓝 (1 - Y01.tv Y01.a)) :=
        tendsto_const_nhds.sub (Y01.PC.contAt_tv Y01.PC.a_pos.ne').tendsto
      rwa [Y01.PC.tv_a, sub_self] at h
    have := (Real.continuous_mul_log.tendsto 0).comp hc
    simp only [Function.comp_def, zero_mul] at this
    exact this.mono_left nhdsWithin_le_nhds
  refine squeeze_zero_norm' ?_ (by simpa using (hS.norm.const_mul 2))
  filter_upwards [Ioo_mem_nhdsLT Y01.PC.one_lt_a] with v hv
  have hv' : v ∈ Ioo Y01.b Y01.a := ⟨Y01.PC.b_lt_one.trans hv.1, hv.2⟩
  have h1 := Y01.PC.l1_abs_le hv
  have h0 : 0 ≤ 1 - Y01.tv v := (Y01.PC.one_sub_tv_pos hv').le
  simp only [Real.norm_eq_abs, abs_mul]
  rw [abs_of_nonneg h0]
  unfold Y01.l2
  nlinarith [abs_nonneg (Real.log (1 - Y01.tv v))]

lemma tendsto_H_a : Tendsto H (𝓝[<] Y01.a) (𝓝 (H Y01.a)) := by
  obtain ⟨hd, -⟩ := Z05.P4.main
  have g1 := (hd Y01.a).1.continuousAt.tendsto.mono_left (nhdsWithin_le_nhds (s := Iio Y01.a))
  have g2 := (hd Y01.a).2.1.continuousAt.tendsto.mono_left (nhdsWithin_le_nhds (s := Iio Y01.a))
  have hl1 := Y01.PC.contAt_l1_a.tendsto.mono_left (nhdsWithin_le_nhds (s := Iio Y01.a))
  rw [Y01.PC.l1_a] at hl1
  have hlim := (g1.mul (hl1.pow 3)).add (g2.mul (hl1.mul tendsto_l1l2_a))
  rw [H_a]
  have e0 : G30 Y01.a * (0:ℝ) ^ 3 + G21 Y01.a * (0 * 0) = 0 := by ring
  rw [e0] at hlim
  exact hlim

lemma contOn_H : ContinuousOn H (Icc Y01.b Y01.a) := by
  intro x hx
  rcases eq_or_lt_of_le hx.1 with h | h
  · subst h
    exact (continuousWithinAt_Icc_iff_Ici Y01.PC.b_lt_a).2
      (continuousWithinAt_Ioi_iff_Ici.1 tendsto_H_b)
  rcases eq_or_lt_of_le hx.2 with h' | h'
  · subst h'
    exact (continuousWithinAt_Icc_iff_Iic Y01.PC.b_lt_a).2
      (continuousWithinAt_Iio_iff_Iic.1 tendsto_H_a)
  exact (hasDerivAt_H ⟨h, h'⟩).continuousAt.continuousWithinAt

lemma int_Hd : ∫ v in Y01.b..Y01.a, (Gp v - Tg v) = 0 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le Y01.PC.b_lt_a.le contOn_H
    (fun x hx => hasDerivAt_H hx) (ii_Gp.sub ii_Tg), H_a, H_b, sub_zero]

theorem main : stmt_3 → stmt_5 := by
  rintro ⟨hS, hM⟩
  refine ⟨ii_Tg, ?_⟩
  have hpt : ∀ v, Y01.l1 v ^ 3 * R30 v + Y01.l1 v ^ 2 * Y01.l2 v * R21 v
      = W4integrand v + Gp v := by
    intro v; unfold R30 R21 W4integrand Gp; ring
  rw [hM]
  congr 1
  have hW := Z05.P1.main
  simp_rw [hpt]
  rw [intervalIntegral.integral_add hW ii_Gp]
  have h0 := int_Hd
  rw [intervalIntegral.integral_sub ii_Gp ii_Tg] at h0
  unfold W4
  have : ∫ v in Y01.b..Y01.a, Gp v = ∫ v in Y01.b..Y01.a, Tg v := by linarith
  rw [this]
  rfl

end Z05.P5

namespace Z05.P7

open MeasureTheory Real Filter Topology Set

theorem main : stmt_5 → stmt_7 := by
  rintro ⟨hT, hM⟩
  obtain ⟨h6i, h6⟩ := Z05.P6.main
  have hpt : ∀ v, Y01.l1 v ^ 2 * (Om20 v + Y01.g20' v) + Y01.l1 v * Y01.l2 v * (Om11 v + g11d v)
      = W3pintegrand v + (Y01.l1 v ^ 2 * Y01.g20' v + Y01.l1 v * Y01.l2 v * g11d v) := by
    intro v; unfold W3pintegrand; ring
  have hW3 : IntervalIntegrable W3pintegrand volume Y01.b Y01.a := by
    have e : W3pintegrand = fun v =>
        (Y01.l1 v ^ 2 * (Om20 v + Y01.g20' v) + Y01.l1 v * Y01.l2 v * (Om11 v + g11d v))
        - (Y01.l1 v ^ 2 * Y01.g20' v + Y01.l1 v * Y01.l2 v * g11d v) := by
      funext v; rw [hpt]; ring
    rw [e]; exact hT.sub h6i
  refine ⟨Z05.P1.main, hW3, Y01.PC.ii_l1Rwt, ?_⟩
  rw [hM]
  congr 1
  simp_rw [hpt]
  rw [intervalIntegral.integral_add hW3 h6i, h6]
  unfold W3pint
  ring

end Z05.P7

namespace Z05.P9

/-- Theorem thm:M3 from Lemma g1:lem:red3 (stmt_7), Lemma g1:lem:W3p (stmt_8) and
Lemma g1:lem:T (`Y01.stmt_8`). -/
theorem main : stmt_7 → stmt_8 → stmt_9 := by
  rintro ⟨-, -, -, hM⟩ h8
  unfold stmt_9
  rw [hM, ← h8, Y01.proof_8.1]

end Z05.P9

namespace Z05.P12

/-- Theorem thm:W4R, "Hence": from stmt_9 and stmt_11. -/
theorem main : stmt_9 → stmt_11 := by
  intro h9 h11
  rw [h9, h11]

end Z05.P12


/-!
# Part P8: `W3pint = W3pclosed` (Lemma g1:lem:W3p)

Exact evaluation as in the supplement's script `g1_M2_raw.py`, with the primitives of Y01.PE1:
`ℓ₁ = ℓ_b + ℓ_{-a} - ℓ_0 - 2λ`, `ℓ₂ = ℓ_a + ℓ_{-b} - ℓ_0 - 2λ` (`ℓ_x = log|v-x|`), every product
`ℓ_xℓ_y` polarised with the endpoint letter in the numerator; only pairs `{b, a}` are split at `1`.
The endpoint pole `ℓ₁²/(v-a)` (Lemma g1:lem:pole) uses `m_x = ℓ_x - log(a-x)` and the new blocks
`fpa`, `fpg` with primitives `FPA = P₂(s)`, `FPG = P₂(s₁) - Φ(mb)`, which need `P2a` continuous
at `s = 1` (`contAt_P2a_one`).  The resulting raw value (atoms `Li₂(u)`, `Li₃(u)`) is reduced
with the tables of Y01.PE2 plus (I3) at `-a²` and `-√2`.  The block/value lemmas `zq*`, `zl*`,
`zb*`, `zv*` and the sums `zS*`, `zU*` were generated by an exact `ℚ(√2)` script.
-/
set_option linter.unusedTactic false
set_option linter.unusedVariables false

namespace Z05.P8

open Real MeasureTheory Set Filter Topology Y01 Y01.PE1

lemma abs_log_le {s : ℝ} (h : 1 / 2 ≤ s) : |Real.log s| ≤ 2 * |s - 1| := by
  have hs : 0 < s := by linarith
  have h1 := Real.log_le_sub_one_of_pos hs
  have h2 : 1 - 1 / s ≤ Real.log s := by
    have := Real.log_le_sub_one_of_pos (inv_pos.mpr hs)
    rw [Real.log_inv] at this
    rw [one_div]; linarith
  rw [abs_le]
  constructor
  · have : 1 - 1 / s ≥ -(2 * |s - 1|) := by
      rw [ge_iff_le, show 1 - 1 / s = (s - 1) / s by field_simp]
      rw [neg_le, ← neg_div, div_le_iff₀ hs]
      have := neg_abs_le (s - 1)
      have h3 : -(s - 1) ≤ |s - 1| := neg_le_abs _
      nlinarith [abs_nonneg (s - 1)]
    linarith
  · have := le_abs_self (s - 1)
    nlinarith [abs_nonneg (s - 1)]

lemma tendsto_logsq_log1 :
    Tendsto (fun s : ℝ => Real.log s ^ 2 * Real.log (1 - s)) (𝓝 1) (𝓝 0) := by
  have hc : Continuous (fun s : ℝ => 4 * |1 - s| * |(1 - s) * Real.log (1 - s)|) :=
    ((continuous_const.mul (continuous_const.sub continuous_id).abs).mul
      (Real.continuous_mul_log.comp (continuous_const.sub continuous_id)).abs)
  have h0 : Tendsto (fun s : ℝ => 4 * |1 - s| * |(1 - s) * Real.log (1 - s)|) (𝓝 1) (𝓝 0) := by
    have := hc.tendsto 1
    simpa using this
  apply squeeze_zero_norm' _ h0
  have hev : ∀ᶠ s in 𝓝 (1 : ℝ), 1 / 2 ≤ s := by
    exact Ici_mem_nhds (by norm_num)
  filter_upwards [hev] with s hs
  rw [Real.norm_eq_abs, abs_mul, abs_pow]
  have hl := abs_log_le hs
  have hl2 : |Real.log s| ^ 2 ≤ 4 * |1 - s| ^ 2 := by
    rw [abs_sub_comm] at hl
    nlinarith [abs_nonneg (Real.log s), abs_nonneg (1 - s)]
  have : 4 * |1 - s| * |(1 - s) * Real.log (1 - s)|
      = 4 * |1 - s| ^ 2 * |Real.log (1 - s)| := by rw [abs_mul]; ring
  rw [this]
  exact mul_le_mul_of_nonneg_right hl2 (abs_nonneg _)

lemma contAt_P2a_one : ContinuousAt P2a 1 := by
  have e : P2a = fun s => Real.log s ^ 2 * Real.log (1 - s)
      + (2 * Real.log s * D00.Li2 s - 2 * Li3 s) := by
    funext s; unfold P2a; ring
  rw [e]
  have h1 : ContinuousAt (fun s : ℝ => Real.log s ^ 2 * Real.log (1 - s)) 1 := by
    rw [ContinuousAt]; simpa using tendsto_logsq_log1
  have h2 : ContinuousAt (fun s : ℝ => 2 * Real.log s * D00.Li2 s - 2 * Li3 s) 1 :=
    ((continuousAt_const.mul (Real.continuousAt_log one_ne_zero)).mul
      D00.continuous_Li2.continuousAt).sub (continuousAt_const.mul continuous_Li3.continuousAt)
  exact h1.add h2

lemma continuous_P2a : Continuous P2a := by
  rw [continuous_iff_continuousAt]
  intro s
  rcases eq_or_ne s 1 with rfl | h
  · exact contAt_P2a_one
  · exact contAt_P2a h

lemma ftc_neg {F f : ℝ → ℝ} {p q : ℝ} (hpq : p ≤ q) (hc : ContinuousOn F (Icc p q))
    (hd : ∀ v ∈ Ioo p q, HasDerivAt F (f v) v) (hf : ∀ v ∈ Ioo p q, f v ≤ 0) :
    IntervalIntegrable f volume p q ∧ ∫ v in p..q, f v = F q - F p := by
  have h := ftc_pos (F := fun v => -F v) (f := fun v => -f v) hpq hc.neg
    (fun v hv => (hd v hv).neg) (fun v hv => by linarith [hf v hv])
  refine ⟨by have h1 := h.1.neg; convert h1 using 1; funext v; simp, ?_⟩
  have h2 := h.2
  rw [intervalIntegral.integral_neg] at h2
  linarith

noncomputable def fpa (x v : ℝ) : ℝ := (Real.log (v - x) - Real.log (a - x)) ^ 2 / (v - a)
noncomputable def fpg (x y v : ℝ) : ℝ :=
  ((Real.log (v - x) - Real.log (a - x)) - (Real.log (v - y) - Real.log (a - y))) ^ 2 / (v - a)
noncomputable def FPA (x v : ℝ) : ℝ := P2a (sA x a v)
noncomputable def FPG (x y v : ℝ) : ℝ :=
  P2a (s1 x y a v) - Phi P1p P2a (Real.log ((a - y) / (a - x))) (mb x y v)

lemma dFPA {x v : ℝ} (hvx : 0 < v - x) (hva : v - a ≠ 0) (hxa : 0 < a - x) :
    HasDerivAt (FPA x) (fpa x v) v := by
  have hxz : a - x ≠ 0 := hxa.ne'
  have n0 := sA_ne_zero hxz hvx.ne'
  have n1 := sA_ne_one hxz hva
  have := (hasDerivAt_P2a n0 n1).comp v (hasDerivAt_sA x a v)
  refine this.congr_deriv ?_
  rw [sA, Real.log_div hvx.ne' hxz]
  unfold fpa
  have e : (v - x) / (a - x) - 1 = (v - a) / (a - x) := by field_simp; ring
  rw [e]
  field_simp

lemma dFPG {x y v : ℝ} (hxy : x - y ≠ 0) (hvx : 0 < v - x) (hvy : 0 < v - y) (hva : v - a ≠ 0)
    (hxa : 0 < a - x) (hya : 0 < a - y) : HasDerivAt (FPG x y) (fpg x y v) v := by
  have hzx : a - x ≠ 0 := hxa.ne'
  have hzy : a - y ≠ 0 := hya.ne'
  have n0 := s1_ne_zero hzx hzy hvx.ne' hvy.ne'
  have n1 := s1_ne_one hxy hzx hvy.ne' hva
  have m0 := mb_ne_zero hvx.ne' hvy.ne'
  have m1 := mb_ne_one hxy hvy.ne'
  have dA := (hasDerivAt_P2a n0 n1).comp v (hasDerivAt_s1 (x := x) (y := y) hzx hvy.ne')
  have dB := (hasDerivAt_Phi (P1 := P1p) (P2 := P2a) (k := Real.log ((a - y) / (a - x))) m1
    (hasDerivAt_P1p m0 m1) (hasDerivAt_P2a m0 m1)).comp v (hasDerivAt_mb (x := x) hvy.ne')
  have := dA.sub dB
  refine this.congr_deriv ?_
  rw [log_s1 hzx hzy hvx.ne' hvy.ne', log_mb hvx.ne' hvy.ne', Real.log_div hzy hzx]
  unfold fpg
  have e1 : s1 x y a v - 1 = (v - a) * (x - y) / ((v - y) * (a - x)) := by
    rw [s1]; field_simp; ring
  have e2 : mb x y v - 1 = (y - x) / (v - y) := by
    rw [mb]; field_simp; ring
  rw [e1, e2]
  have hyx : y - x ≠ 0 := by intro h; apply hxy; linarith
  field_simp
  ring

lemma intPA {x p q : ℝ} (hpq : p ≤ q) (hqa : q ≤ a) (hxa : 0 < a - x)
    (hx : ∀ v ∈ Ioo p q, 0 < v - x) :
    IntervalIntegrable (fpa x) volume p q ∧ ∫ v in p..q, fpa x v = FPA x q - FPA x p :=
  ftc_neg hpq (continuous_P2a.comp (continuous_sA x a)).continuousOn
    (fun v hv => dFPA (hx v hv) (by linarith [hv.2]) hxa)
    (fun v hv => div_nonpos_of_nonneg_of_nonpos (sq_nonneg _) (by linarith [hv.2]))

lemma intPG {x y p q : ℝ} (hpq : p ≤ q) (hqa : q ≤ a) (hxy : x - y ≠ 0) (hxa : 0 < a - x)
    (hya : 0 < a - y) (hy : ∀ v ∈ Icc p q, 0 < v - y) (hx : ∀ v ∈ Ioo p q, 0 < v - x) :
    IntervalIntegrable (fpg x y) volume p q ∧ ∫ v in p..q, fpg x y v = FPG x y q - FPG x y p := by
  have hc : ContinuousOn (FPG x y) (Icc p q) := by
    have c1 : ContinuousOn (fun v => P2a (s1 x y a v)) (Icc p q) :=
      continuous_P2a.comp_continuousOn (contOn_s1 hxa.ne' (fun v hv => (hy v hv).ne'))
    have c2 : ContinuousOn (fun v => Phi P1p P2a (Real.log ((a - y) / (a - x))) (mb x y v))
        (Icc p q) :=
      contOn_comp (contOn_mb (fun v hv => (hy v hv).ne'))
        (fun v hv => cPhiM (mb_ne_one hxy (hy v hv).ne'))
    exact c1.sub c2
  exact ftc_neg hpq hc
    (fun v hv => dFPG hxy (hx v hv) (hy v (Ioo_subset_Icc_self hv)) (by linarith [hv.2]) hxa hya)
    (fun v hv => div_nonpos_of_nonneg_of_nonpos (sq_nonneg _) (by linarith [hv.2]))

lemma zq0 : sA a 0 b = 2 * b := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq1 : (1 : ℝ) - (2 * b) = b ^ 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq2 : sA a 0 a = 0 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq3 : s1 0 (-b) (-a) b = b := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq4 : mb 0 (-b) b = 1 / 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq5 : (1 : ℝ) - (b) = 2 - Real.sqrt 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq6 : (1 : ℝ) - (1 / 2) = 1 / 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq7 : s1 0 (-b) (-a) a = 1 / Real.sqrt 2 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq8 : mb 0 (-b) a = (2 + Real.sqrt 2) / 4 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq9 : (1 : ℝ) - ((2 + Real.sqrt 2) / 4) = (2 - Real.sqrt 2) / 4 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq10 : mb a (-a) b = -(1 / Real.sqrt 2) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq11 : (1 : ℝ) - (-(1 / Real.sqrt 2)) = (Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 1 * a ^ 0)) := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq12 : mb a (-a) a = 0 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq13 : sA (-b) 0 b = (Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq14 : (1 : ℝ) - (Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) = -1 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq15 : (Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0))⁻¹ = 1 / 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq16 : sA (-b) 0 a = (Real.sqrt 2 ^ 3 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq17 : (1 : ℝ) - (Real.sqrt 2 ^ 3 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) = -(a ^ 2) := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq18 : (Real.sqrt 2 ^ 3 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0))⁻¹ = (2 - Real.sqrt 2) / 4 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq19 : (1 : ℝ) - ((2 - Real.sqrt 2) / 4) = (2 + Real.sqrt 2) / 4 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq20 : mb b (-a) b = 0 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq21 : mb b (-a) a = b := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq22 : s1 a (-a) 0 b = 1 / Real.sqrt 2 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq23 : s1 a (-a) 0 a = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq24 : sA b 0 b = 0 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq25 : sA b 0 a = -(2 * a) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq26 : (1 : ℝ) - (-(2 * a)) = (Real.sqrt 2 ^ 0 * a ^ 2 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq27 : sA (-a) 0 b = (Real.sqrt 2 ^ 3 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 1)) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq28 : (1 : ℝ) - (Real.sqrt 2 ^ 3 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 1)) = -(b ^ 2) := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq29 : sA (-a) 0 a = (Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq30 : s1 (-b) (-a) 0 b = (Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 1 * a ^ 0)) := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq31 : mb (-b) (-a) b = 1 - 1 / Real.sqrt 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq32 : (1 : ℝ) - (Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 1 * a ^ 0)) = -(1 / Real.sqrt 2) := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq33 : (Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 1 * a ^ 0))⁻¹ = 2 - Real.sqrt 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq34 : (1 : ℝ) - (2 - Real.sqrt 2) = b := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq35 : (1 : ℝ) - (1 - 1 / Real.sqrt 2) = 1 / Real.sqrt 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq36 : s1 (-b) (-a) 0 a = (Real.sqrt 2 ^ 1 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq37 : mb (-b) (-a) a = 2 - Real.sqrt 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq38 : (1 : ℝ) - (Real.sqrt 2 ^ 1 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) = -a := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq39 : (Real.sqrt 2 ^ 1 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0))⁻¹ = 1 - 1 / Real.sqrt 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq40 : s1 b 0 (-a) b = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq41 : mb b 0 b = 0 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq42 : s1 b 0 (-a) a = 1 / Real.sqrt 2 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq43 : mb b 0 a = 2 * b := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq44 : s1 b (-b) 0 b = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq45 : mb b (-b) b = 0 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq46 : s1 b (-b) 0 a = -(1 / Real.sqrt 2) := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq47 : mb b (-b) a = 1 / Real.sqrt 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq48 : mb 0 (-a) b = (2 - Real.sqrt 2) / 4 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq49 : mb 0 (-a) a = 1 / 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq50 : mb (-a) 0 b = (Real.sqrt 2 ^ 3 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq51 : mb (-a) 0 a = (Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq52 : s1 b (-a) (-b) b = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq53 : s1 b (-a) (-b) a = -1 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq54 : (1 : ℝ) - (-1) = (Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq55 : (Real.sqrt 2 ^ 3 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 1))⁻¹ = (2 + Real.sqrt 2) / 4 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq56 : s1 b 0 (-b) b = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq57 : s1 b 0 (-b) a = b := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq58 : s1 0 (-a) (-b) b = -(1 / Real.sqrt 2) := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq59 : s1 0 (-a) (-b) a = -a := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq60 : (1 : ℝ) - (-a) = (Real.sqrt 2 ^ 1 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq61 : sA a (-a) b = b := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq62 : sA a (-a) a = 0 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq63 : sA (-b) (-a) b = -b := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq64 : (1 : ℝ) - (-b) = (Real.sqrt 2 ^ 1 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq65 : sA (-b) (-a) a = -Real.sqrt 2 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq66 : (1 : ℝ) - (-Real.sqrt 2) = (Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq67 : s1 a 0 (-a) b = -a := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq68 : mb a 0 b = -(2 * a) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq69 : s1 a 0 (-a) a = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq70 : mb a 0 a = 0 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq71 : sA b (-a) b = 0 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq72 : sA b (-a) a = -(1 / Real.sqrt 2) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq73 : sA 0 (-b) b = -1 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq74 : sA 0 (-b) a = -(a ^ 2) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq75 : (1 : ℝ) - (-(a ^ 2)) = (Real.sqrt 2 ^ 3 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq76 : sA b (-b) b = 0 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq77 : sA b (-b) a = -a := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq78 : sA (-a) (-b) b = (Real.sqrt 2 ^ 1 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq79 : (1 : ℝ) - (Real.sqrt 2 ^ 1 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) = -b := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq80 : sA (-a) (-b) a = (Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq81 : (1 : ℝ) - (Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) = -Real.sqrt 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq82 : s1 b (-a) 0 b = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq83 : s1 b (-a) 0 a = -a := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq84 : sA 0 (-a) b = -(b ^ 2) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq85 : (1 : ℝ) - (-(b ^ 2)) = (Real.sqrt 2 ^ 3 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 1)) := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq86 : sA 0 (-a) a = -1 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq87 : s1 b (-b) (-a) b = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq88 : s1 b (-b) (-a) a = 1 / 2 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq89 : mb (-b) 0 b = (Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq90 : mb (-b) 0 a = (Real.sqrt 2 ^ 3 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 1)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq91 : (Real.sqrt 2 ^ 1 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0))⁻¹ = 1 / Real.sqrt 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq92 : (Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0))⁻¹ = b := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq93 : sA b a b = 0 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq94 : sA b a a = 1 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq95 : sA (-a) a b = 2 - Real.sqrt 2 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq96 : sA (-a) a a = 1 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq97 : sA 0 a b = b ^ 2 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq98 : (1 : ℝ) - (b ^ 2) = 2 * b := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq99 : sA 0 a a = 1 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq100 : s1 b (-a) a b = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq101 : s1 b (-a) a a = 1 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq102 : s1 b 0 a b = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq103 : s1 b 0 a a = 1 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq104 : s1 0 (-a) a b = 1 - 1 / Real.sqrt 2 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq105 : s1 0 (-a) a a = 1 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq106 : s1 b a (-a) b = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq107 : mb b a b = 0 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq108 : s1 b a (-a) 1 = -(1 / Real.sqrt 2) := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq109 : mb b a 1 = -b := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq110 : s1 b a 0 b = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq111 : s1 b a 0 1 = -a := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq112 : s1 a b (-a) 1 = -Real.sqrt 2 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq113 : mb a b 1 = -a := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq114 : s1 a b (-a) a = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq115 : mb a b a = 0 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq116 : s1 a b 0 1 = -b := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zq117 : s1 a b 0 a = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)

lemma zl0 : Real.log (b - 0) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
  logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl1 : Real.log (a - 0) = (0 : ℝ) * lam + (1 : ℝ) * L :=
  logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl2 : Real.log (0 - a) = (0 : ℝ) * lam + (1 : ℝ) * L :=
  logN 0 1 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl3 : Real.log (b ^ 2) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
  logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl4 : Real.log (((-a) - 0) / ((-a) - (-b))) = (-1 : ℝ) * lam + (1 : ℝ) * L :=
  logP 0 1 2 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl5 : Real.log (2 - Real.sqrt 2) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
  logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl6 : Real.log (b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
  logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl7 : Real.log (1 / 2) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
  logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl8 : Real.log (1 - 1 / Real.sqrt 2) = (-1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
  logP 0 0 1 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl9 : Real.log (1 / Real.sqrt 2) = (-1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
  logP 0 0 1 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl10 : Real.log ((2 + Real.sqrt 2) / 4) = (-3 / 2 : ℝ) * lam + (1 : ℝ) * L :=
  logP 0 1 3 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl11 : Real.log ((2 - Real.sqrt 2) / 4) = (-3 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
  logP 0 0 3 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl12 : Real.log (-(1 / Real.sqrt 2)) = (-1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
  logN 0 0 1 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl13 : Real.log ((Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 1 * a ^ 0))) = (-1 / 2 : ℝ) * lam + (1 : ℝ) * L :=
  logP 0 1 1 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl14 : Real.log (0 - (-b)) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
  logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl15 : Real.log (-1) = (0 : ℝ) * lam + (0 : ℝ) * L :=
  logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl16 : Real.log (-(a ^ 2)) = (0 : ℝ) * lam + (2 : ℝ) * L :=
  logN 0 2 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl17 : Real.log ((0 - a) / (0 - (-a))) = (0 : ℝ) * lam + (0 : ℝ) * L :=
  logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl18 : Real.log (0 - b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
  logN 0 0 0 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl19 : Real.log ((Real.sqrt 2 ^ 0 * a ^ 2 / (Real.sqrt 2 ^ 0 * a ^ 0))) = (0 : ℝ) * lam + (2 : ℝ) * L :=
  logP 0 2 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl20 : Real.log (-(2 * a)) = (1 : ℝ) * lam + (1 : ℝ) * L :=
  logN 2 1 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl21 : Real.log (b - (-a)) = (3 / 2 : ℝ) * lam + (0 : ℝ) * L :=
  logP 3 0 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl22 : Real.log (a - (-a)) = (1 : ℝ) * lam + (1 : ℝ) * L :=
  logP 2 1 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl23 : Real.log (0 - (-a)) = (0 : ℝ) * lam + (1 : ℝ) * L :=
  logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl24 : Real.log (-(b ^ 2)) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
  logN 0 0 0 2 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl25 : Real.log ((0 - (-b)) / (0 - (-a))) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
  logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl26 : Real.log (-a) = (0 : ℝ) * lam + (1 : ℝ) * L :=
  logN 0 1 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl27 : Real.log (b - (-b)) = (1 : ℝ) * lam + (-1 : ℝ) * L :=
  logP 2 0 0 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl28 : Real.log (a - (-b)) = (3 / 2 : ℝ) * lam + (0 : ℝ) * L :=
  logP 3 0 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl29 : Real.log (((-a) - b) / ((-a) - 0)) = (3 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
  logP 3 0 0 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl30 : Real.log (2 * b) = (1 : ℝ) * lam + (-1 : ℝ) * L :=
  logP 2 0 0 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl31 : Real.log ((0 - b) / (0 - (-b))) = (0 : ℝ) * lam + (0 : ℝ) * L :=
  logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl32 : Real.log (((-b) - b) / ((-b) - (-a))) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
  logN 0 0 0 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl33 : Real.log ((Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0))) = (1 : ℝ) * lam + (0 : ℝ) * L :=
  logP 2 0 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl34 : Real.log (((-b) - b) / ((-b) - 0)) = (1 : ℝ) * lam + (0 : ℝ) * L :=
  logP 2 0 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl35 : Real.log (((-b) - 0) / ((-b) - (-a))) = (-1 : ℝ) * lam + (-1 : ℝ) * L :=
  logN 0 0 2 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl36 : Real.log ((Real.sqrt 2 ^ 1 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0))) = (1 / 2 : ℝ) * lam + (1 : ℝ) * L :=
  logP 1 1 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl37 : Real.log ((-a) - a) = (1 : ℝ) * lam + (1 : ℝ) * L :=
  logN 2 1 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl38 : Real.log ((-a) - (-b)) = (1 : ℝ) * lam + (0 : ℝ) * L :=
  logN 2 0 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl39 : Real.log ((Real.sqrt 2 ^ 1 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0))) = (1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
  logP 1 0 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl40 : Real.log (-b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
  logN 0 0 0 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl41 : Real.log ((Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0))) = (0 : ℝ) * lam + (1 : ℝ) * L :=
  logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl42 : Real.log (-Real.sqrt 2) = (1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
  logN 1 0 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl43 : Real.log (((-a) - a) / ((-a) - 0)) = (1 : ℝ) * lam + (0 : ℝ) * L :=
  logP 2 0 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl44 : Real.log ((-a) - b) = (3 / 2 : ℝ) * lam + (0 : ℝ) * L :=
  logN 3 0 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl45 : Real.log ((-b) - 0) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
  logN 0 0 0 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl46 : Real.log ((Real.sqrt 2 ^ 3 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0))) = (3 / 2 : ℝ) * lam + (1 : ℝ) * L :=
  logP 3 1 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl47 : Real.log ((-b) - b) = (1 : ℝ) * lam + (-1 : ℝ) * L :=
  logN 2 0 0 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl48 : Real.log ((-b) - (-a)) = (1 : ℝ) * lam + (0 : ℝ) * L :=
  logP 2 0 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl49 : Real.log ((0 - b) / (0 - (-a))) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
  logN 0 0 0 2 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl50 : Real.log ((-a) - 0) = (0 : ℝ) * lam + (1 : ℝ) * L :=
  logN 0 1 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl51 : Real.log ((Real.sqrt 2 ^ 3 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 1))) = (3 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
  logP 3 0 0 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl52 : Real.log (((-a) - b) / ((-a) - (-b))) = (1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
  logP 1 0 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl53 : Real.log ((a - (-a)) / (a - b)) = (0 : ℝ) * lam + (1 : ℝ) * L :=
  logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl54 : Real.log ((a - 0) / (a - b)) = (-1 : ℝ) * lam + (1 : ℝ) * L :=
  logP 0 1 2 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl55 : Real.log ((a - (-a)) / (a - 0)) = (1 : ℝ) * lam + (0 : ℝ) * L :=
  logP 2 0 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl56 : Real.log (((-a) - b) / ((-a) - a)) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
  logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl57 : Real.log ((0 - b) / (0 - a)) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
  logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl58 : Real.log (((-a) - a) / ((-a) - b)) = (-1 / 2 : ℝ) * lam + (1 : ℝ) * L :=
  logP 0 1 1 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zl59 : Real.log ((0 - a) / (0 - b)) = (0 : ℝ) * lam + (2 : ℝ) * L :=
  logP 0 2 0 0 _ _ (by norm_num) (by norm_num)
    (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))

lemma zb0 : IntervalIntegrable (flin 0 0) volume b a ∧
    ∫ v in b..a, flin 0 0 v = FH 0 a - FH 0 b := by
  have hx : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intH hba hx

lemma zb1 : IntervalIntegrable (flin a 0) volume b a ∧
    ∫ v in b..a, flin a 0 v = FL P1p a 0 a - FL P1p a 0 b := by
  have hxz : ((0 : ℝ) - a : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1p hba hxz hz hx (fun v hv => cPsiM (sA_ne_one hxz (hz v hv).ne'))

lemma zb2 : IntervalIntegrable (fpr 0 (-b) (-a)) volume b a ∧
    ∫ v in b..a, fpr 0 (-b) (-a) v = FG P1p P2a P2a 0 (-b) (-a) a - FG P1p P2a P2a 0 (-b) (-a) b := by
  have hxy : ((0 : ℝ) - (-b) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (((-a) : ℝ) - 0 : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (((-a) : ℝ) - (-b) : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1p isP2a isP2a hba hxy hzx hzy hz hy hx
    (fun v hv => cPhiM (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma zb3 : IntervalIntegrable (fpr a (-a) (-a)) volume b a ∧
    ∫ v in b..a, fpr a (-a) (-a) v = FIV P2a a (-a) a - FIV P2a a (-a) b := by
  have hxy : ((a : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hy : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIV isP2a hba hxy hy hx (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv).ne'))

lemma zb4 : IntervalIntegrable (fsq (-b) 0) volume b a ∧
    ∫ v in b..a, fsq (-b) 0 v = FA P1p P2g (-b) 0 a - FA P1p P2g (-b) 0 b := by
  have hxz : ((0 : ℝ) - (-b) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx2 : ∀ v ∈ Icc b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1p isP2g hba hxz hz hx (fun v hv => cPhiG (sA_ne_zero hxz (hx2 v hv)) (sA_ne_one hxz (hz v hv).ne'))

lemma zb5 : IntervalIntegrable (fpr b (-a) (-a)) volume b a ∧
    ∫ v in b..a, fpr b (-a) (-a) v = FIV P2a b (-a) a - FIV P2a b (-a) b := by
  have hxy : ((b : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hy : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIV isP2a hba hxy hy hx (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv).ne'))

lemma zb6 : IntervalIntegrable (fpr a (-a) 0) volume b a ∧
    ∫ v in b..a, fpr a (-a) 0 v = FG P1p P2a P2a a (-a) 0 a - FG P1p P2a P2a a (-a) 0 b := by
  have hxy : ((a : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : ((0 : ℝ) - a : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : ((0 : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b a, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1p isP2a isP2a hba hxy hzx hzy hz hy hx
    (fun v hv => cPhiM (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma zb7 : IntervalIntegrable (fsq b 0) volume b a ∧
    ∫ v in b..a, fsq b 0 v = FA P1n P2a b 0 a - FA P1n P2a b 0 b := by
  have hxz : ((0 : ℝ) - b : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1n isP2a hba hxz hz hx (fun v hv => cPhiN (sA_ne_one hxz (hz v hv).ne'))

lemma zb8 : IntervalIntegrable (fcst (-a)) volume b a ∧
    ∫ v in b..a, fcst (-a) v = FK (-a) a - FK (-a) b := by
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intK hba hz

lemma zb9 : IntervalIntegrable (flin (-a) 0) volume b a ∧
    ∫ v in b..a, flin (-a) 0 v = FL P1p (-a) 0 a - FL P1p (-a) 0 b := by
  have hxz : ((0 : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1p hba hxz hz hx (fun v hv => cPsiM (sA_ne_one hxz (hz v hv).ne'))

lemma zb10 : IntervalIntegrable (flin b 0) volume b a ∧
    ∫ v in b..a, flin b 0 v = FL P1n b 0 a - FL P1n b 0 b := by
  have hxz : ((0 : ℝ) - b : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1n hba hxz hz hx (fun v hv => cPsiN (sA_ne_one hxz (hz v hv).ne'))

lemma zb11 : IntervalIntegrable (fpr (-b) (-a) 0) volume b a ∧
    ∫ v in b..a, fpr (-b) (-a) 0 v = FG P1p P2g P2a (-b) (-a) 0 a - FG P1p P2g P2a (-b) (-a) 0 b := by
  have hxy : (((-b) : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : ((0 : ℝ) - (-b) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : ((0 : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx2 : ∀ v ∈ Icc b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1p isP2g isP2a hba hxy hzx hzy hz hy hx
    (fun v hv => cPhiG (s1_ne_zero hzx hzy (hx2 v hv) (hy v hv)) (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma zb12 : IntervalIntegrable (fcst (-b)) volume b a ∧
    ∫ v in b..a, fcst (-b) v = FK (-b) a - FK (-b) b := by
  have hz : ∀ v ∈ Icc b a, 0 < v - (-b) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intK hba hz

lemma zb13 : IntervalIntegrable (fsq 0 0) volume b a ∧
    ∫ v in b..a, fsq 0 0 v = FC 0 a - FC 0 b := by
  have hx : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intC hba hx

lemma zb14 : IntervalIntegrable (fpr b 0 (-a)) volume b a ∧
    ∫ v in b..a, fpr b 0 (-a) v = FG P1p P2a P2a b 0 (-a) a - FG P1p P2a P2a b 0 (-a) b := by
  have hxy : ((b : ℝ) - 0 : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (((-a) : ℝ) - b : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (((-a) : ℝ) - 0 : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1p isP2a isP2a hba hxy hzx hzy hz hy hx
    (fun v hv => cPhiM (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma zb15 : IntervalIntegrable (fpr b (-b) 0) volume b a ∧
    ∫ v in b..a, fpr b (-b) 0 v = FG P1n P2a P2a b (-b) 0 a - FG P1n P2a P2a b (-b) 0 b := by
  have hxy : ((b : ℝ) - (-b) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : ((0 : ℝ) - b : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : ((0 : ℝ) - (-b) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1n isP2a isP2a hba hxy hzx hzy hz hy hx
    (fun v hv => cPhiN (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma zb16 : IntervalIntegrable (fpr 0 (-a) (-a)) volume b a ∧
    ∫ v in b..a, fpr 0 (-a) (-a) v = FIV P2a 0 (-a) a - FIV P2a 0 (-a) b := by
  have hxy : ((0 : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hy : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIV isP2a hba hxy hy hx (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv).ne'))

lemma zb17 : IntervalIntegrable (fpr (-a) 0 0) volume b a ∧
    ∫ v in b..a, fpr (-a) 0 0 v = FIV P2g (-a) 0 a - FIV P2g (-a) 0 b := by
  have hxy : (((-a) : ℝ) - 0 : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hy : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx2 : ∀ v ∈ Icc b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIV isP2g hba hxy hy hx (fun v hv => contAt_P2g (mb_ne_zero (hx2 v hv) (hy v hv).ne') (mb_ne_one hxy (hy v hv).ne'))

lemma zb18 : IntervalIntegrable (fpr b (-a) (-b)) volume b a ∧
    ∫ v in b..a, fpr b (-a) (-b) v = FG P1n P2a P2a b (-a) (-b) a - FG P1n P2a P2a b (-a) (-b) b := by
  have hxy : ((b : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (((-b) : ℝ) - b : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (((-b) : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-b) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1n isP2a isP2a hba hxy hzx hzy hz hy hx
    (fun v hv => cPhiN (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma zb19 : IntervalIntegrable (flin (-a) (-a)) volume b a ∧
    ∫ v in b..a, flin (-a) (-a) v = FH (-a) a - FH (-a) b := by
  have hx : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intH hba hx

lemma zb20 : IntervalIntegrable (fsq (-a) 0) volume b a ∧
    ∫ v in b..a, fsq (-a) 0 v = FA P1p P2g (-a) 0 a - FA P1p P2g (-a) 0 b := by
  have hxz : ((0 : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx2 : ∀ v ∈ Icc b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1p isP2g hba hxz hz hx (fun v hv => cPhiG (sA_ne_zero hxz (hx2 v hv)) (sA_ne_one hxz (hz v hv).ne'))

lemma zb21 : IntervalIntegrable (fpr b 0 (-b)) volume b a ∧
    ∫ v in b..a, fpr b 0 (-b) v = FG P1p P2a P2a b 0 (-b) a - FG P1p P2a P2a b 0 (-b) b := by
  have hxy : ((b : ℝ) - 0 : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (((-b) : ℝ) - b : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (((-b) : ℝ) - 0 : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-b) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1p isP2a isP2a hba hxy hzx hzy hz hy hx
    (fun v hv => cPhiM (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma zb22 : IntervalIntegrable (fpr 0 (-a) (-b)) volume b a ∧
    ∫ v in b..a, fpr 0 (-a) (-b) v = FG P1n P2a P2a 0 (-a) (-b) a - FG P1n P2a P2a 0 (-a) (-b) b := by
  have hxy : ((0 : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (((-b) : ℝ) - 0 : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (((-b) : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-b) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1n isP2a isP2a hba hxy hzx hzy hz hy hx
    (fun v hv => cPhiN (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma zb23 : IntervalIntegrable (flin a (-a)) volume b a ∧
    ∫ v in b..a, flin a (-a) v = FL P1p a (-a) a - FL P1p a (-a) b := by
  have hxz : (((-a) : ℝ) - a : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1p hba hxz hz hx (fun v hv => cPsiM (sA_ne_one hxz (hz v hv).ne'))

lemma zb24 : IntervalIntegrable (fsq (-b) (-a)) volume b a ∧
    ∫ v in b..a, fsq (-b) (-a) v = FA P1n P2a (-b) (-a) a - FA P1n P2a (-b) (-a) b := by
  have hxz : (((-a) : ℝ) - (-b) : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1n isP2a hba hxz hz hx (fun v hv => cPhiN (sA_ne_one hxz (hz v hv).ne'))

lemma zb25 : IntervalIntegrable (flin (-b) (-a)) volume b a ∧
    ∫ v in b..a, flin (-b) (-a) v = FL P1n (-b) (-a) a - FL P1n (-b) (-a) b := by
  have hxz : (((-a) : ℝ) - (-b) : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1n hba hxz hz hx (fun v hv => cPsiN (sA_ne_one hxz (hz v hv).ne'))

lemma zb26 : IntervalIntegrable (fsq a 0) volume b a ∧
    ∫ v in b..a, fsq a 0 v = FA P1p P2a a 0 a - FA P1p P2a a 0 b := by
  have hxz : ((0 : ℝ) - a : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1p isP2a hba hxz hz hx (fun v hv => cPhiM (sA_ne_one hxz (hz v hv).ne'))

lemma zb27 : IntervalIntegrable (fpr a 0 (-a)) volume b a ∧
    ∫ v in b..a, fpr a 0 (-a) v = FG P1n P2a P2a a 0 (-a) a - FG P1n P2a P2a a 0 (-a) b := by
  have hxy : ((a : ℝ) - 0 : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (((-a) : ℝ) - a : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (((-a) : ℝ) - 0 : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b a, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1n isP2a isP2a hba hxy hzx hzy hz hy hx
    (fun v hv => cPhiN (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma zb28 : IntervalIntegrable (fpr a 0 0) volume b a ∧
    ∫ v in b..a, fpr a 0 0 v = FIV P2a a 0 a - FIV P2a a 0 b := by
  have hxy : ((a : ℝ) - 0 : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hy : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIV isP2a hba hxy hy hx (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv).ne'))

lemma zb29 : IntervalIntegrable (flin b (-a)) volume b a ∧
    ∫ v in b..a, flin b (-a) v = FL P1n b (-a) a - FL P1n b (-a) b := by
  have hxz : (((-a) : ℝ) - b : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1n hba hxz hz hx (fun v hv => cPsiN (sA_ne_one hxz (hz v hv).ne'))

lemma zb30 : IntervalIntegrable (fsq b (-a)) volume b a ∧
    ∫ v in b..a, fsq b (-a) v = FA P1n P2a b (-a) a - FA P1n P2a b (-a) b := by
  have hxz : (((-a) : ℝ) - b : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1n isP2a hba hxz hz hx (fun v hv => cPhiN (sA_ne_one hxz (hz v hv).ne'))

lemma zb31 : IntervalIntegrable (flin 0 (-b)) volume b a ∧
    ∫ v in b..a, flin 0 (-b) v = FL P1n 0 (-b) a - FL P1n 0 (-b) b := by
  have hxz : (((-b) : ℝ) - 0 : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-b) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1n hba hxz hz hx (fun v hv => cPsiN (sA_ne_one hxz (hz v hv).ne'))

lemma zb32 : IntervalIntegrable (fsq (-a) (-a)) volume b a ∧
    ∫ v in b..a, fsq (-a) (-a) v = FC (-a) a - FC (-a) b := by
  have hx : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intC hba hx

lemma zb33 : IntervalIntegrable (fpr b 0 0) volume b a ∧
    ∫ v in b..a, fpr b 0 0 v = FIV P2a b 0 a - FIV P2a b 0 b := by
  have hxy : ((b : ℝ) - 0 : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hy : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIV isP2a hba hxy hy hx (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv).ne'))

lemma zb34 : IntervalIntegrable (fsq b (-b)) volume b a ∧
    ∫ v in b..a, fsq b (-b) v = FA P1n P2a b (-b) a - FA P1n P2a b (-b) b := by
  have hxz : (((-b) : ℝ) - b : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-b) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1n isP2a hba hxz hz hx (fun v hv => cPhiN (sA_ne_one hxz (hz v hv).ne'))

lemma zb35 : IntervalIntegrable (fcst 0) volume b a ∧
    ∫ v in b..a, fcst 0 v = FK 0 a - FK 0 b := by
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intK hba hz

lemma zb36 : IntervalIntegrable (flin (-b) 0) volume b a ∧
    ∫ v in b..a, flin (-b) 0 v = FL P1p (-b) 0 a - FL P1p (-b) 0 b := by
  have hxz : ((0 : ℝ) - (-b) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1p hba hxz hz hx (fun v hv => cPsiM (sA_ne_one hxz (hz v hv).ne'))

lemma zb37 : IntervalIntegrable (flin (-a) (-b)) volume b a ∧
    ∫ v in b..a, flin (-a) (-b) v = FL P1p (-a) (-b) a - FL P1p (-a) (-b) b := by
  have hxz : (((-b) : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-b) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1p hba hxz hz hx (fun v hv => cPsiM (sA_ne_one hxz (hz v hv).ne'))

lemma zb38 : IntervalIntegrable (flin b (-b)) volume b a ∧
    ∫ v in b..a, flin b (-b) v = FL P1n b (-b) a - FL P1n b (-b) b := by
  have hxz : (((-b) : ℝ) - b : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-b) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1n hba hxz hz hx (fun v hv => cPsiN (sA_ne_one hxz (hz v hv).ne'))

lemma zb39 : IntervalIntegrable (fpr b (-a) 0) volume b a ∧
    ∫ v in b..a, fpr b (-a) 0 v = FG P1n P2a P2a b (-a) 0 a - FG P1n P2a P2a b (-a) 0 b := by
  have hxy : ((b : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : ((0 : ℝ) - b : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : ((0 : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1n isP2a isP2a hba hxy hzx hzy hz hy hx
    (fun v hv => cPhiN (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma zb40 : IntervalIntegrable (fsq 0 (-a)) volume b a ∧
    ∫ v in b..a, fsq 0 (-a) v = FA P1n P2a 0 (-a) a - FA P1n P2a 0 (-a) b := by
  have hxz : (((-a) : ℝ) - 0 : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1n isP2a hba hxz hz hx (fun v hv => cPhiN (sA_ne_one hxz (hz v hv).ne'))

lemma zb41 : IntervalIntegrable (flin 0 (-a)) volume b a ∧
    ∫ v in b..a, flin 0 (-a) v = FL P1n 0 (-a) a - FL P1n 0 (-a) b := by
  have hxz : (((-a) : ℝ) - 0 : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1n hba hxz hz hx (fun v hv => cPsiN (sA_ne_one hxz (hz v hv).ne'))

lemma zb42 : IntervalIntegrable (fpr (-b) (-a) (-a)) volume b a ∧
    ∫ v in b..a, fpr (-b) (-a) (-a) v = FIV P2a (-b) (-a) a - FIV P2a (-b) (-a) b := by
  have hxy : (((-b) : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hy : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIV isP2a hba hxy hy hx (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv).ne'))

lemma zb43 : IntervalIntegrable (fpr b (-b) (-a)) volume b a ∧
    ∫ v in b..a, fpr b (-b) (-a) v = FG P1p P2a P2a b (-b) (-a) a - FG P1p P2a P2a b (-b) (-a) b := by
  have hxy : ((b : ℝ) - (-b) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (((-a) : ℝ) - b : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (((-a) : ℝ) - (-b) : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1p isP2a isP2a hba hxy hzx hzy hz hy hx
    (fun v hv => cPhiM (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma zb44 : IntervalIntegrable (fsq a (-a)) volume b a ∧
    ∫ v in b..a, fsq a (-a) v = FA P1p P2a a (-a) a - FA P1p P2a a (-a) b := by
  have hxz : (((-a) : ℝ) - a : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1p isP2a hba hxz hz hx (fun v hv => cPhiM (sA_ne_one hxz (hz v hv).ne'))

lemma zb45 : IntervalIntegrable (fpr (-b) 0 0) volume b a ∧
    ∫ v in b..a, fpr (-b) 0 0 v = FIV P2g (-b) 0 a - FIV P2g (-b) 0 b := by
  have hxy : (((-b) : ℝ) - 0 : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hy : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx2 : ∀ v ∈ Icc b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIV isP2g hba hxy hy hx (fun v hv => contAt_P2g (mb_ne_zero (hx2 v hv) (hy v hv).ne') (mb_ne_one hxy (hy v hv).ne'))

lemma zb46 : IntervalIntegrable (fsq (-a) (-b)) volume b a ∧
    ∫ v in b..a, fsq (-a) (-b) v = FA P1p P2g (-a) (-b) a - FA P1p P2g (-a) (-b) b := by
  have hxz : (((-b) : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-b) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx2 : ∀ v ∈ Icc b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1p isP2g hba hxz hz hx (fun v hv => cPhiG (sA_ne_zero hxz (hx2 v hv)) (sA_ne_one hxz (hz v hv).ne'))

lemma zb47 : IntervalIntegrable (fsq 0 (-b)) volume b a ∧
    ∫ v in b..a, fsq 0 (-b) v = FA P1n P2a 0 (-b) a - FA P1n P2a 0 (-b) b := by
  have hxz : (((-b) : ℝ) - 0 : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-b) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1n isP2a hba hxz hz hx (fun v hv => cPhiN (sA_ne_one hxz (hz v hv).ne'))

lemma zb48 : IntervalIntegrable (fpa b) volume b a ∧
    ∫ v in b..a, fpa b v = FPA b a - FPA b b := by
  have hxa : (0 : ℝ) < a - b := by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, 0 < v - b := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intPA hba le_rfl hxa hx

lemma zb49 : IntervalIntegrable (fpa (-a)) volume b a ∧
    ∫ v in b..a, fpa (-a) v = FPA (-a) a - FPA (-a) b := by
  have hxa : (0 : ℝ) < a - (-a) := by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intPA hba le_rfl hxa hx

lemma zb50 : IntervalIntegrable (fpa 0) volume b a ∧
    ∫ v in b..a, fpa 0 v = FPA 0 a - FPA 0 b := by
  have hxa : (0 : ℝ) < a - 0 := by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intPA hba le_rfl hxa hx

lemma zb51 : IntervalIntegrable (fpg b (-a)) volume b a ∧
    ∫ v in b..a, fpg b (-a) v = FPG b (-a) a - FPG b (-a) b := by
  have hxy : ((b : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hxa : (0 : ℝ) < a - b := by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi]
  have hya : (0 : ℝ) < a - (-a) := by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, 0 < v - b := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intPG hba le_rfl hxy hxa hya hy hx

lemma zb52 : IntervalIntegrable (fpg b 0) volume b a ∧
    ∫ v in b..a, fpg b 0 v = FPG b 0 a - FPG b 0 b := by
  have hxy : ((b : ℝ) - 0 : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hxa : (0 : ℝ) < a - b := by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi]
  have hya : (0 : ℝ) < a - 0 := by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, 0 < v - b := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intPG hba le_rfl hxy hxa hya hy hx

lemma zb53 : IntervalIntegrable (fpg 0 (-a)) volume b a ∧
    ∫ v in b..a, fpg 0 (-a) v = FPG 0 (-a) a - FPG 0 (-a) b := by
  have hxy : ((0 : ℝ) - (-a) : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hxa : (0 : ℝ) < a - 0 := by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi]
  have hya : (0 : ℝ) < a - (-a) := by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intPG hba le_rfl hxy hxa hya hy hx

lemma zb54 : IntervalIntegrable (fpr b a (-a)) volume b 1 ∧
    ∫ v in b..1, fpr b a (-a) v = FG P1n P2a P2a b a (-a) 1 - FG P1n P2a P2a b a (-a) b := by
  have hxy : ((b : ℝ) - a : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (((-a) : ℝ) - b : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (((-a) : ℝ) - a : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b 1, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b 1, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b 1, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1n isP2a isP2a hb1 hxy hzx hzy hz hy hx
    (fun v hv => cPhiN (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma zb55 : IntervalIntegrable (fpr b a 0) volume b 1 ∧
    ∫ v in b..1, fpr b a 0 v = FG P1n P2a P2a b a 0 1 - FG P1n P2a P2a b a 0 b := by
  have hxy : ((b : ℝ) - a : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : ((0 : ℝ) - b : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : ((0 : ℝ) - a : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b 1, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b 1, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b 1, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1n isP2a isP2a hb1 hxy hzx hzy hz hy hx
    (fun v hv => cPhiN (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma zb56 : IntervalIntegrable (fpr a b (-a)) volume 1 a ∧
    ∫ v in 1..a, fpr a b (-a) v = FG P1n P2a P2a a b (-a) a - FG P1n P2a P2a a b (-a) 1 := by
  have hxy : ((a : ℝ) - b : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (((-a) : ℝ) - a : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (((-a) : ℝ) - b : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc 1 a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc 1 a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo 1 a, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1n isP2a isP2a h1a hxy hzx hzy hz hy hx
    (fun v hv => cPhiN (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma zb57 : IntervalIntegrable (fpr a b 0) volume 1 a ∧
    ∫ v in 1..a, fpr a b 0 v = FG P1n P2a P2a a b 0 a - FG P1n P2a P2a a b 0 1 := by
  have hxy : ((a : ℝ) - b : ℝ) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : ((0 : ℝ) - a : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : ((0 : ℝ) - b : ℝ) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc 1 a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc 1 a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo 1 a, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1n isP2a isP2a h1a hxy hzx hzy hz hy hx
    (fun v hv => cPhiN (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma zv0 : FH 0 b = (1 / 2 : ℝ) * L ^ 2 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zl0]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv1 : FH 0 a = (1 / 2 : ℝ) * L ^ 2 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zl1]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv2 : FL P1p a 0 b = (-2 : ℝ) * L ^ 2 + (-1 : ℝ) * D00.Li2 (b ^ 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq1, zq0, zl2, zl3]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv3 : FL P1p a 0 a = (-1 / 6 : ℝ) * Real.pi ^ 2 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq2, zl2]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv4 : FG P1p P2a P2a 0 (-b) (-a) b = (3 / 2 : ℝ) * lam ^ 3 + (-2 : ℝ) * lam ^ 2 * L + (3 : ℝ) * lam * L ^ 2 + (2 : ℝ) * lam * D00.Li2 (2 - Real.sqrt 2) + (2 : ℝ) * lam * D00.Li2 (1 / 2) + (-2 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (2 - Real.sqrt 2) + (-2 : ℝ) * L * D00.Li2 (b) + (-2 : ℝ) * Li3 (b) + (2 : ℝ) * Li3 (1 / 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq6, zq5, zq4, zq3, zl4, zl7, zl5, zl6]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv5 : FG P1p P2a P2a 0 (-b) (-a) a = (11 / 4 : ℝ) * lam ^ 3 + (-5 / 2 : ℝ) * lam ^ 2 * L + (2 : ℝ) * lam * D00.Li2 (1 - 1 / Real.sqrt 2) + (-1 : ℝ) * lam * D00.Li2 (1 / Real.sqrt 2) + (3 : ℝ) * lam * D00.Li2 ((2 + Real.sqrt 2) / 4) + (-2 : ℝ) * L * D00.Li2 (1 - 1 / Real.sqrt 2) + (-2 : ℝ) * L * D00.Li2 ((2 + Real.sqrt 2) / 4) + (-2 : ℝ) * Li3 (1 / Real.sqrt 2) + (2 : ℝ) * Li3 ((2 + Real.sqrt 2) / 4) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq9, zq8, zq7, zl4, zl10, zl11, zl8, zl9]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv6 : FIV P2a a (-a) b = (1 / 8 : ℝ) * lam ^ 3 + (-1 / 4 : ℝ) * lam ^ 2 * L + (1 : ℝ) * lam * D00.Li2 (-(1 / Real.sqrt 2)) + (2 : ℝ) * Li3 (-(1 / Real.sqrt 2)) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq11, zq10, zl13, zl12]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv7 : FIV P2a a (-a) a = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq12]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv8 : FA P1p P2g (-b) 0 b = (-2 / 3 : ℝ) * lam ^ 3 + (-2 : ℝ) * lam * D00.Li2 (1 / 2) + (2 : ℝ) * L * D00.Li2 (-1) + (-2 : ℝ) * Li3 (1 / 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq6, zq14, zq15, zq13, zl15, zl14, zl7]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv9 : FA P1p P2g (-b) 0 a = (-9 / 4 : ℝ) * lam ^ 3 + (3 : ℝ) * lam * L ^ 2 + (-3 : ℝ) * lam * D00.Li2 ((2 - Real.sqrt 2) / 4) + (10 / 3 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 (-(a ^ 2)) + (-2 : ℝ) * L * D00.Li2 ((2 - Real.sqrt 2) / 4) + (-2 : ℝ) * Li3 ((2 - Real.sqrt 2) / 4) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq19, zq17, zq18, zq16, zl10, zl11, zl16, zl14]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv10 : FIV P2a b (-a) b = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq20]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv11 : FIV P2a b (-a) a = (-1 / 2 : ℝ) * lam * L ^ 2 + (1 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 (b) + (2 : ℝ) * Li3 (b) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq5, zq21, zl5, zl6]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv12 : FG P1p P2a P2a a (-a) 0 b = (-1 / 2 : ℝ) * lam ^ 2 * L + (-1 : ℝ) * lam * D00.Li2 (1 / Real.sqrt 2) + (1 : ℝ) * lam * D00.Li2 (-(1 / Real.sqrt 2)) + (-2 : ℝ) * Li3 (1 / Real.sqrt 2) + (2 : ℝ) * Li3 (-(1 / Real.sqrt 2)) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq11, zq10, zq22, zl17, zl13, zl12, zl8, zl9]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv13 : FG P1p P2a P2a a (-a) 0 a = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq12, zq23, zl17]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv14 : FA P1n P2a b 0 b = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq24, zl18]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv15 : FA P1n P2a b 0 a = (2 : ℝ) * lam ^ 2 * L + (2 : ℝ) * lam * D00.Li2 (-(2 * a)) + (-2 : ℝ) * Li3 (-(2 * a)) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq26, zq25, zl19, zl20, zl18]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv16 : FK (-a) b = (3 / 2 : ℝ) * lam := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zl21]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv17 : FK (-a) a = (1 : ℝ) * lam + (1 : ℝ) * L := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zl22]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv18 : FL P1p (-a) 0 b = (-2 : ℝ) * L ^ 2 + (-1 : ℝ) * D00.Li2 (-(b ^ 2)) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq28, zq27, zl24, zl23]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv19 : FL P1p (-a) 0 a = (-1 : ℝ) * D00.Li2 (-1) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq14, zq29, zl15, zl23]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv20 : FL P1n b 0 b = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq24, zl18]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv21 : FL P1n b 0 a = (2 : ℝ) * lam * L + (1 : ℝ) * D00.Li2 (-(2 * a)) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq26, zq25, zl19, zl20, zl18]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv22 : FG P1p P2g P2a (-b) (-a) 0 b = (1 / 12 : ℝ) * lam ^ 3 + (1 / 2 : ℝ) * lam ^ 2 * L + (-1 : ℝ) * lam * L ^ 2 + (1 : ℝ) * lam * D00.Li2 (2 - Real.sqrt 2) + (1 : ℝ) * lam * D00.Li2 (1 - 1 / Real.sqrt 2) + (-2 / 3 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (2 - Real.sqrt 2) + (2 : ℝ) * L * D00.Li2 (1 - 1 / Real.sqrt 2) + (4 : ℝ) * L * D00.Li2 (-(1 / Real.sqrt 2)) + (-2 : ℝ) * Li3 (2 - Real.sqrt 2) + (2 : ℝ) * Li3 (1 - 1 / Real.sqrt 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq35, zq34, zq32, zq33, zq31, zq30, zl25, zl12, zl8, zl9, zl5, zl6]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv23 : FG P1p P2g P2a (-b) (-a) 0 a = (-1 / 12 : ℝ) * lam ^ 3 + (-1 : ℝ) * lam * L ^ 2 + (-1 : ℝ) * lam * D00.Li2 (2 - Real.sqrt 2) + (-1 : ℝ) * lam * D00.Li2 (1 - 1 / Real.sqrt 2) + (16 / 3 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 (2 - Real.sqrt 2) + (-2 : ℝ) * L * D00.Li2 (1 - 1 / Real.sqrt 2) + (4 : ℝ) * L * D00.Li2 (-a) + (2 : ℝ) * Li3 (2 - Real.sqrt 2) + (-2 : ℝ) * Li3 (1 - 1 / Real.sqrt 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq35, zq34, zq38, zq39, zq37, zq36, zl25, zl26, zl8, zl9, zl5, zl6]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv24 : FK (-b) b = (1 : ℝ) * lam + (-1 : ℝ) * L := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zl27]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv25 : FK (-b) a = (3 / 2 : ℝ) * lam := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zl28]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv26 : FC 0 b = (-1 / 3 : ℝ) * L ^ 3 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zl0]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv27 : FC 0 a = (1 / 3 : ℝ) * L ^ 3 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zl1]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv28 : FG P1p P2a P2a b 0 (-a) b = (-1 / 2 : ℝ) * lam * Real.pi ^ 2 + (1 / 3 : ℝ) * L * Real.pi ^ 2 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq41, zq40, zl29]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv29 : FG P1p P2a P2a b 0 (-a) a = (-5 / 4 : ℝ) * lam ^ 3 + (1 : ℝ) * lam ^ 2 * L + (-3 / 2 : ℝ) * lam * L ^ 2 + (-3 : ℝ) * lam * D00.Li2 (1 - 1 / Real.sqrt 2) + (-1 : ℝ) * lam * D00.Li2 (1 / Real.sqrt 2) + (-2 : ℝ) * lam * D00.Li2 (2 * b) + (1 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 (1 - 1 / Real.sqrt 2) + (2 : ℝ) * L * D00.Li2 (2 * b) + (-2 : ℝ) * Li3 (1 / Real.sqrt 2) + (2 : ℝ) * Li3 (2 * b) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq1, zq43, zq42, zl29, zl8, zl9, zl30, zl3]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv30 : FG P1n P2a P2a b (-b) 0 b = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq45, zq44, zl31]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv31 : FG P1n P2a P2a b (-b) 0 a = (1 / 2 : ℝ) * lam ^ 2 * L + (1 : ℝ) * lam * D00.Li2 (1 / Real.sqrt 2) + (-1 : ℝ) * lam * D00.Li2 (-(1 / Real.sqrt 2)) + (2 : ℝ) * Li3 (1 / Real.sqrt 2) + (-2 : ℝ) * Li3 (-(1 / Real.sqrt 2)) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq11, zq47, zq46, zl31, zl13, zl12, zl8, zl9]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv32 : FIV P2a 0 (-a) b = (27 / 8 : ℝ) * lam ^ 3 + (9 / 4 : ℝ) * lam ^ 2 * L + (-3 / 2 : ℝ) * lam * L ^ 2 + (3 : ℝ) * lam * D00.Li2 ((2 - Real.sqrt 2) / 4) + (-1 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 ((2 - Real.sqrt 2) / 4) + (2 : ℝ) * Li3 ((2 - Real.sqrt 2) / 4) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq19, zq48, zl10, zl11]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv33 : FIV P2a 0 (-a) a = (1 : ℝ) * lam ^ 3 + (2 : ℝ) * lam * D00.Li2 (1 / 2) + (2 : ℝ) * Li3 (1 / 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq6, zq49, zl7]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv34 : FIV P2g (-a) 0 b = (9 / 4 : ℝ) * lam ^ 3 + (-3 : ℝ) * lam * L ^ 2 + (3 : ℝ) * lam * D00.Li2 ((2 - Real.sqrt 2) / 4) + (-4 / 3 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 ((2 - Real.sqrt 2) / 4) + (2 : ℝ) * Li3 ((2 - Real.sqrt 2) / 4) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq19, zq18, zq50, zl10, zl11]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv35 : FIV P2g (-a) 0 a = (2 / 3 : ℝ) * lam ^ 3 + (2 : ℝ) * lam * D00.Li2 (1 / 2) + (2 : ℝ) * Li3 (1 / 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq6, zq15, zq51, zl7]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv36 : FG P1n P2a P2a b (-a) (-b) b = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq20, zq52, zl32]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv37 : FG P1n P2a P2a b (-a) (-b) a = (1 / 2 : ℝ) * lam * L ^ 2 + (1 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 (b) + (-2 : ℝ) * L * D00.Li2 (-1) + (2 : ℝ) * Li3 (b) + (-2 : ℝ) * Li3 (-1) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq54, zq5, zq21, zq53, zl32, zl33, zl15, zl5, zl6]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv38 : FH (-a) b = (9 / 8 : ℝ) * lam ^ 2 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zl21]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv39 : FH (-a) a = (1 / 2 : ℝ) * lam ^ 2 + (1 : ℝ) * lam * L + (1 / 2 : ℝ) * L ^ 2 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zl22]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv40 : FA P1p P2g (-a) 0 b = (-9 / 4 : ℝ) * lam ^ 3 + (3 : ℝ) * lam * L ^ 2 + (-3 : ℝ) * lam * D00.Li2 ((2 + Real.sqrt 2) / 4) + (-10 / 3 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 ((2 + Real.sqrt 2) / 4) + (-2 : ℝ) * L * D00.Li2 (-(b ^ 2)) + (-2 : ℝ) * Li3 ((2 + Real.sqrt 2) / 4) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq9, zq28, zq55, zq27, zl10, zl11, zl24, zl23]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv41 : FA P1p P2g (-a) 0 a = (-2 / 3 : ℝ) * lam ^ 3 + (-2 : ℝ) * lam * D00.Li2 (1 / 2) + (-2 : ℝ) * L * D00.Li2 (-1) + (-2 : ℝ) * Li3 (1 / 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq6, zq14, zq15, zq29, zl15, zl23, zl7]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv42 : FG P1p P2a P2a b 0 (-b) b = (-1 / 3 : ℝ) * lam * Real.pi ^ 2 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq41, zq56, zl34]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv43 : FG P1p P2a P2a b 0 (-b) a = (1 / 2 : ℝ) * lam ^ 3 + (1 : ℝ) * lam ^ 2 * L + (-7 / 2 : ℝ) * lam * L ^ 2 + (-2 : ℝ) * lam * D00.Li2 (2 - Real.sqrt 2) + (-2 : ℝ) * lam * D00.Li2 (2 * b) + (1 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (b) + (2 : ℝ) * L * D00.Li2 (2 * b) + (-2 : ℝ) * Li3 (b) + (2 : ℝ) * Li3 (2 * b) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq1, zq5, zq43, zq57, zl34, zl30, zl5, zl6, zl3]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv44 : FG P1n P2a P2a 0 (-a) (-b) b = (9 / 4 : ℝ) * lam ^ 3 + (3 : ℝ) * lam ^ 2 * L + (1 : ℝ) * lam * L ^ 2 + (-3 : ℝ) * lam * D00.Li2 (-(1 / Real.sqrt 2)) + (3 : ℝ) * lam * D00.Li2 ((2 - Real.sqrt 2) / 4) + (-2 : ℝ) * L * D00.Li2 (-(1 / Real.sqrt 2)) + (2 : ℝ) * L * D00.Li2 ((2 - Real.sqrt 2) / 4) + (-2 : ℝ) * Li3 (-(1 / Real.sqrt 2)) + (2 : ℝ) * Li3 ((2 - Real.sqrt 2) / 4) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq19, zq11, zq48, zq58, zl35, zl10, zl11, zl13, zl12]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv45 : FG P1n P2a P2a 0 (-a) (-b) a = (3 / 2 : ℝ) * lam ^ 3 + (1 : ℝ) * lam ^ 2 * L + (2 : ℝ) * lam * D00.Li2 (1 / 2) + (-2 : ℝ) * lam * D00.Li2 (-a) + (2 : ℝ) * Li3 (1 / 2) + (-2 : ℝ) * Li3 (-a) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq60, zq6, zq49, zq59, zl35, zl36, zl26, zl7]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv46 : FL P1p a (-a) b = (1 / 2 : ℝ) * lam ^ 2 + (-1 / 2 : ℝ) * lam * L + (-1 : ℝ) * L ^ 2 + (-1 : ℝ) * D00.Li2 (2 - Real.sqrt 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq5, zq61, zl37, zl5]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv47 : FL P1p a (-a) a = (-1 / 6 : ℝ) * Real.pi ^ 2 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq62, zl37]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv48 : FA P1n P2a (-b) (-a) b = (1 / 2 : ℝ) * lam ^ 3 + (-1 : ℝ) * lam ^ 2 * L + (1 / 2 : ℝ) * lam * L ^ 2 + (2 : ℝ) * lam * D00.Li2 (-b) + (-2 : ℝ) * L * D00.Li2 (-b) + (-2 : ℝ) * Li3 (-b) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq64, zq63, zl38, zl39, zl40]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv49 : FA P1n P2a (-b) (-a) a = (9 / 4 : ℝ) * lam ^ 2 * L + (3 : ℝ) * lam * D00.Li2 (-Real.sqrt 2) + (-2 : ℝ) * Li3 (-Real.sqrt 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq66, zq65, zl38, zl41, zl42]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv50 : FL P1n (-b) (-a) b = (1 / 2 : ℝ) * lam ^ 2 + (-1 / 2 : ℝ) * lam * L + (1 : ℝ) * D00.Li2 (-b) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq64, zq63, zl38, zl39, zl40]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv51 : FL P1n (-b) (-a) a = (3 / 2 : ℝ) * lam * L + (1 : ℝ) * D00.Li2 (-Real.sqrt 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq66, zq65, zl38, zl41, zl42]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv52 : FA P1p P2a a 0 b = (-2 : ℝ) * lam ^ 2 * L + (4 : ℝ) * lam * L ^ 2 + (2 : ℝ) * lam * D00.Li2 (2 * b) + (-4 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (b ^ 2) + (-2 : ℝ) * L * D00.Li2 (2 * b) + (-2 : ℝ) * Li3 (2 * b) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq1, zq0, zl2, zl30, zl3]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv53 : FA P1p P2a a 0 a = (-1 / 3 : ℝ) * L * Real.pi ^ 2 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq2, zl2]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv54 : FG P1n P2a P2a a 0 (-a) b = (1 / 2 : ℝ) * lam ^ 3 + (-3 / 2 : ℝ) * lam * L ^ 2 + (-2 : ℝ) * lam * D00.Li2 (-(2 * a)) + (2 : ℝ) * lam * D00.Li2 (-a) + (-1 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (-(2 * a)) + (2 : ℝ) * L * D00.Li2 (-a) + (2 : ℝ) * Li3 (-(2 * a)) + (-2 : ℝ) * Li3 (-a) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq26, zq60, zq68, zq67, zl43, zl19, zl36, zl20, zl26]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv55 : FG P1n P2a P2a a 0 (-a) a = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq70, zq69, zl43]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv56 : FIV P2a a 0 b = (-2 : ℝ) * lam ^ 2 * L + (-4 : ℝ) * lam * L ^ 2 + (-2 : ℝ) * lam * D00.Li2 (-(2 * a)) + (-2 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (-(2 * a)) + (2 : ℝ) * Li3 (-(2 * a)) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq26, zq68, zl19, zl20]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv57 : FIV P2a a 0 a = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq70]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv58 : FL P1n b (-a) b = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq71, zl44]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv59 : FL P1n b (-a) a = (-1 / 2 : ℝ) * lam ^ 2 + (1 : ℝ) * lam * L + (1 : ℝ) * D00.Li2 (-(1 / Real.sqrt 2)) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq11, zq72, zl44, zl13, zl12]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv60 : FA P1n P2a b (-a) b = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq71, zl44]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv61 : FA P1n P2a b (-a) a = (-1 / 2 : ℝ) * lam ^ 3 + (1 : ℝ) * lam ^ 2 * L + (2 : ℝ) * lam * D00.Li2 (-(1 / Real.sqrt 2)) + (-2 : ℝ) * Li3 (-(1 / Real.sqrt 2)) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq11, zq72, zl44, zl13, zl12]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv62 : FL P1n 0 (-b) b = (-1 : ℝ) * lam * L + (1 : ℝ) * D00.Li2 (-1) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq54, zq73, zl45, zl33, zl15]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv63 : FL P1n 0 (-b) a = (3 / 2 : ℝ) * lam * L + (1 : ℝ) * L ^ 2 + (1 : ℝ) * D00.Li2 (-(a ^ 2)) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq75, zq74, zl45, zl46, zl16]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv64 : FC (-a) b = (9 / 8 : ℝ) * lam ^ 3 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zl21]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv65 : FC (-a) a = (1 / 3 : ℝ) * lam ^ 3 + (1 : ℝ) * lam ^ 2 * L + (1 : ℝ) * lam * L ^ 2 + (1 / 3 : ℝ) * L ^ 3 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zl22]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv66 : FIV P2a b 0 b = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq41]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv67 : FIV P2a b 0 a = (2 : ℝ) * lam ^ 2 * L + (-4 : ℝ) * lam * L ^ 2 + (-2 : ℝ) * lam * D00.Li2 (2 * b) + (2 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 (2 * b) + (2 : ℝ) * Li3 (2 * b) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq1, zq43, zl30, zl3]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv68 : FA P1n P2a b (-b) b = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq76, zl47]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv69 : FA P1n P2a b (-b) a = (1 / 2 : ℝ) * lam ^ 3 + (1 : ℝ) * lam ^ 2 * L + (2 : ℝ) * lam * D00.Li2 (-a) + (-2 : ℝ) * Li3 (-a) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq60, zq77, zl47, zl36, zl26]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv70 : FK 0 b = (-1 : ℝ) * L := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zl0]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv71 : FK 0 a = (1 : ℝ) * L := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zl1]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv72 : FL P1p (-b) 0 b = (-1 : ℝ) * D00.Li2 (-1) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq14, zq13, zl15, zl14]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv73 : FL P1p (-b) 0 a = (-2 : ℝ) * L ^ 2 + (-1 : ℝ) * D00.Li2 (-(a ^ 2)) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq17, zq16, zl16, zl14]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv74 : FL P1p (-a) (-b) b = (-1 : ℝ) * lam * L + (-1 : ℝ) * D00.Li2 (-b) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq79, zq78, zl48, zl40]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv75 : FL P1p (-a) (-b) a = (1 / 2 : ℝ) * lam ^ 2 + (-1 : ℝ) * D00.Li2 (-Real.sqrt 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq81, zq80, zl48, zl42]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv76 : FL P1n b (-b) b = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq76, zl47]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv77 : FL P1n b (-b) a = (1 / 2 : ℝ) * lam ^ 2 + (1 : ℝ) * lam * L + (1 : ℝ) * D00.Li2 (-a) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq60, zq77, zl47, zl36, zl26]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv78 : FG P1n P2a P2a b (-a) 0 b = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq20, zq82, zl49]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv79 : FG P1n P2a P2a b (-a) 0 a = (2 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 (b) + (-2 : ℝ) * L * D00.Li2 (-a) + (2 : ℝ) * Li3 (b) + (-2 : ℝ) * Li3 (-a) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq60, zq5, zq21, zq83, zl49, zl36, zl26, zl5, zl6]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv80 : FA P1n P2a 0 (-a) b = (3 / 2 : ℝ) * lam * L ^ 2 + (-1 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (-(b ^ 2)) + (-2 : ℝ) * Li3 (-(b ^ 2)) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq85, zq84, zl50, zl51, zl24]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv81 : FA P1n P2a 0 (-a) a = (1 : ℝ) * lam * L ^ 2 + (2 : ℝ) * L * D00.Li2 (-1) + (-2 : ℝ) * Li3 (-1) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq54, zq86, zl50, zl33, zl15]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv82 : FL P1n 0 (-a) b = (-3 / 2 : ℝ) * lam * L + (1 : ℝ) * L ^ 2 + (1 : ℝ) * D00.Li2 (-(b ^ 2)) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq85, zq84, zl50, zl51, zl24]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv83 : FL P1n 0 (-a) a = (1 : ℝ) * lam * L + (1 : ℝ) * D00.Li2 (-1) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq54, zq86, zl50, zl33, zl15]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv84 : FIV P2a (-b) (-a) b = (1 / 8 : ℝ) * lam ^ 3 + (1 / 2 : ℝ) * lam ^ 2 * L + (1 / 2 : ℝ) * lam * L ^ 2 + (1 : ℝ) * lam * D00.Li2 (1 - 1 / Real.sqrt 2) + (2 : ℝ) * L * D00.Li2 (1 - 1 / Real.sqrt 2) + (2 : ℝ) * Li3 (1 - 1 / Real.sqrt 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq35, zq31, zl8, zl9]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv85 : FIV P2a (-b) (-a) a = (1 / 4 : ℝ) * lam ^ 2 * L + (-1 : ℝ) * lam * L ^ 2 + (-1 : ℝ) * lam * D00.Li2 (2 - Real.sqrt 2) + (1 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 (2 - Real.sqrt 2) + (2 : ℝ) * Li3 (2 - Real.sqrt 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq34, zq37, zl5, zl6]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv86 : FG P1p P2a P2a b (-b) (-a) b = (-1 / 6 : ℝ) * lam * Real.pi ^ 2 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq45, zq87, zl52]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv87 : FG P1p P2a P2a b (-b) (-a) a = (-9 / 8 : ℝ) * lam ^ 3 + (1 / 4 : ℝ) * lam ^ 2 * L + (-3 : ℝ) * lam * D00.Li2 (1 / 2) + (1 : ℝ) * lam * D00.Li2 (1 / Real.sqrt 2) + (-2 : ℝ) * Li3 (1 / 2) + (2 : ℝ) * Li3 (1 / Real.sqrt 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq6, zq47, zq88, zl52, zl8, zl7, zl9]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv88 : FA P1p P2a a (-a) b = (1 / 2 : ℝ) * lam ^ 3 + (-1 : ℝ) * lam * L ^ 2 + (-2 : ℝ) * lam * D00.Li2 (2 - Real.sqrt 2) + (-2 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (2 - Real.sqrt 2) + (-2 : ℝ) * L * D00.Li2 (b) + (-2 : ℝ) * Li3 (b) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq5, zq61, zl37, zl5, zl6]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv89 : FA P1p P2a a (-a) a = (-1 / 3 : ℝ) * lam * Real.pi ^ 2 + (-1 / 3 : ℝ) * L * Real.pi ^ 2 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq62, zl37]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv90 : FIV P2g (-b) 0 b = (2 / 3 : ℝ) * lam ^ 3 + (2 : ℝ) * lam * D00.Li2 (1 / 2) + (2 : ℝ) * Li3 (1 / 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq6, zq15, zq89, zl7]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv91 : FIV P2g (-b) 0 a = (9 / 4 : ℝ) * lam ^ 3 + (-3 : ℝ) * lam * L ^ 2 + (3 : ℝ) * lam * D00.Li2 ((2 + Real.sqrt 2) / 4) + (4 / 3 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 ((2 + Real.sqrt 2) / 4) + (2 : ℝ) * Li3 ((2 + Real.sqrt 2) / 4) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq9, zq55, zq90, zl10, zl11]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv92 : FA P1p P2g (-a) (-b) b = (-1 / 12 : ℝ) * lam ^ 3 + (-5 / 4 : ℝ) * lam ^ 2 * L + (-1 : ℝ) * lam * D00.Li2 (1 / Real.sqrt 2) + (-2 : ℝ) * lam * D00.Li2 (-b) + (-2 : ℝ) * Li3 (1 / Real.sqrt 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq79, zq91, zq78, zl48, zl40, zl8, zl9]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv93 : FA P1p P2g (-a) (-b) a = (1 / 2 : ℝ) * lam ^ 3 + (1 / 2 : ℝ) * lam * L ^ 2 + (-2 : ℝ) * lam * D00.Li2 (-Real.sqrt 2) + (-2 / 3 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (b) + (-2 : ℝ) * Li3 (b) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq81, zq5, zq92, zq80, zl48, zl42, zl5, zl6]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv94 : FA P1n P2a 0 (-b) b = (1 : ℝ) * lam * L ^ 2 + (-2 : ℝ) * L * D00.Li2 (-1) + (-2 : ℝ) * Li3 (-1) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq54, zq73, zl45, zl33, zl15]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv95 : FA P1n P2a 0 (-b) a = (3 / 2 : ℝ) * lam * L ^ 2 + (1 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 (-(a ^ 2)) + (-2 : ℝ) * Li3 (-(a ^ 2)) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq75, zq74, zl45, zl46, zl16]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv96 : FPA b b = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq93]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv97 : FPA b a = (-2 : ℝ) * z3 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq94]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv98 : FPA (-a) b = (-1 / 4 : ℝ) * lam ^ 2 * L + (1 : ℝ) * lam * L ^ 2 + (1 : ℝ) * lam * D00.Li2 (2 - Real.sqrt 2) + (-1 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (2 - Real.sqrt 2) + (-2 : ℝ) * Li3 (2 - Real.sqrt 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq34, zq95, zl5, zl6]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv99 : FPA (-a) a = (-2 : ℝ) * z3 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq96]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv100 : FPA 0 b = (4 : ℝ) * lam * L ^ 2 + (-4 : ℝ) * L ^ 3 + (-4 : ℝ) * L * D00.Li2 (b ^ 2) + (-2 : ℝ) * Li3 (b ^ 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq98, zq97, zl30, zl3]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv101 : FPA 0 a = (-2 : ℝ) * z3 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq99]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv102 : FPG b (-a) b = (1 / 3 : ℝ) * L * Real.pi ^ 2 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq20, zq100, zl53]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv103 : FPG b (-a) a = (-1 : ℝ) * lam * L ^ 2 + (2 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 (2 - Real.sqrt 2) + (2 : ℝ) * L * D00.Li2 (b) + (-2 : ℝ) * z3 + (2 : ℝ) * Li3 (b) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq5, zq21, zq101, zl53, zl5, zl6]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv104 : FPG b 0 b = (-1 / 3 : ℝ) * lam * Real.pi ^ 2 + (1 / 3 : ℝ) * L * Real.pi ^ 2 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq41, zq102, zl54]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv105 : FPG b 0 a = (4 : ℝ) * lam ^ 2 * L + (-8 : ℝ) * lam * L ^ 2 + (-2 : ℝ) * lam * D00.Li2 (b ^ 2) + (-2 : ℝ) * lam * D00.Li2 (2 * b) + (4 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 (b ^ 2) + (2 : ℝ) * L * D00.Li2 (2 * b) + (-2 : ℝ) * z3 + (2 : ℝ) * Li3 (2 * b) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq1, zq43, zq103, zl54, zl30, zl3]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv106 : FPG 0 (-a) b = (19 / 4 : ℝ) * lam ^ 3 + (3 / 4 : ℝ) * lam ^ 2 * L + (-2 : ℝ) * lam * L ^ 2 + (-1 : ℝ) * lam * D00.Li2 (1 - 1 / Real.sqrt 2) + (2 : ℝ) * lam * D00.Li2 ((2 + Real.sqrt 2) / 4) + (3 : ℝ) * lam * D00.Li2 ((2 - Real.sqrt 2) / 4) + (-1 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (1 - 1 / Real.sqrt 2) + (2 : ℝ) * L * D00.Li2 ((2 - Real.sqrt 2) / 4) + (2 : ℝ) * Li3 ((2 - Real.sqrt 2) / 4) + (-2 : ℝ) * Li3 (1 - 1 / Real.sqrt 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq19, zq35, zq48, zq104, zl10, zl11, zl55, zl8, zl9]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv107 : FPG 0 (-a) a = (2 : ℝ) * lam ^ 3 + (4 : ℝ) * lam * D00.Li2 (1 / 2) + (-2 : ℝ) * z3 + (2 : ℝ) * Li3 (1 / 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq6, zq49, zq105, zl55, zl7]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv108 : FG P1n P2a P2a b a (-a) b = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq107, zq106, zl56]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv109 : FG P1n P2a P2a b a (-a) 1 = (-1 : ℝ) * lam * L ^ 2 + (1 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (-(1 / Real.sqrt 2)) + (2 : ℝ) * L * D00.Li2 (-b) + (-2 : ℝ) * Li3 (-(1 / Real.sqrt 2)) + (2 : ℝ) * Li3 (-b) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq11, zq64, zq109, zq108, zl56, zl13, zl39, zl12, zl40]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv110 : FG P1n P2a P2a b a 0 b = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq107, zq110, zl57]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv111 : FG P1n P2a P2a b a 0 1 = (1 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (-a) + (2 : ℝ) * L * D00.Li2 (-b) + (-2 : ℝ) * Li3 (-a) + (2 : ℝ) * Li3 (-b) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq60, zq64, zq109, zq111, zl57, zl39, zl36, zl26, zl40]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv112 : FG P1n P2a P2a a b (-a) 1 = (-1 / 2 : ℝ) * lam * L ^ 2 + (-2 : ℝ) * L * D00.Li2 (-a) + (2 : ℝ) * L * D00.Li2 (-Real.sqrt 2) + (2 : ℝ) * Li3 (-a) + (-2 : ℝ) * Li3 (-Real.sqrt 2) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq66, zq60, zq113, zq112, zl58, zl41, zl36, zl42, zl26]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv113 : FG P1n P2a P2a a b (-a) a = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq115, zq114, zl58]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv114 : FG P1n P2a P2a a b 0 1 = (-1 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (-a) + (2 : ℝ) * L * D00.Li2 (-b) + (2 : ℝ) * Li3 (-a) + (-2 : ℝ) * Li3 (-b) := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq60, zq64, zq113, zq116, zl59, zl39, zl36, zl26, zl40]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

lemma zv115 : FG P1n P2a P2a a b 0 a = 0 := by
  simp only [FA, FG, FL, FK, FH, FC, FIV, FIII, FPA, FPG, Phi, Psi, P1n, P1p, P2a, P2g, E3, zq115, zq117, zl59]
  try simp only [Real.log_zero, sub_zero, sub_self, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one', Li3_one, inv_one]
  all_goals ring

noncomputable def zS0 (v : ℝ) : ℝ :=
    ((2 : ℝ) * lam) * flin 0 0 v +
    (((-1 / 2 : ℝ) * Real.sqrt 2 * lam + (-1 / 2 : ℝ) * lam) * flin a 0 v +
    (((-1 / 4 : ℝ)) * fpr 0 (-b) (-a) v +
    (((1 / 4 : ℝ)) * fpr a (-a) (-a) v +
    (((1 / 8 : ℝ) * Real.sqrt 2 + (1 / 8 : ℝ)) * fsq (-b) 0 v +
    (((1 / 4 : ℝ)) * fpr b (-a) (-a) v +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (-1 / 8 : ℝ)) * fpr a (-a) 0 v +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (3 / 8 : ℝ)) * fsq b 0 v +
    (((-3 : ℝ) * lam ^ 2) * fcst (-a) v +
    (((1 / 2 : ℝ) * Real.sqrt 2 * lam + (-3 / 2 : ℝ) * lam) * flin (-a) 0 v +
    (((1 / 2 : ℝ) * Real.sqrt 2 * lam + (-3 / 2 : ℝ) * lam) * flin b 0 v +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (-1 / 8 : ℝ)) * fpr (-b) (-a) 0 v +
    (((-2 : ℝ) * lam ^ 2) * fcst (-b) v +
    (((-1 / 2 : ℝ)) * fsq 0 0 v +
    (((-1 / 2 : ℝ)) * fpr b 0 (-a) v +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (-1 / 8 : ℝ)) * fpr b (-b) 0 v +
    (((-1 / 2 : ℝ)) * fpr 0 (-a) (-a) v +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (3 / 8 : ℝ)) * fpr (-a) 0 0 v +
    (((1 / 2 : ℝ)) * fpr b (-a) (-b) v +
    (((2 : ℝ) * lam) * flin (-a) (-a) v +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (3 / 8 : ℝ)) * fsq (-a) 0 v +
    (((-1 / 2 : ℝ)) * fpr b 0 (-b) v +
    (((-1 / 2 : ℝ)) * fpr 0 (-a) (-b) v +
    (((1 : ℝ) * lam) * flin a (-a) v +
    (((-1 / 4 : ℝ)) * fsq (-b) (-a) v +
    (((1 : ℝ) * lam) * flin (-b) (-a) v +
    (((1 / 8 : ℝ) * Real.sqrt 2 + (1 / 8 : ℝ)) * fsq a 0 v +
    (((-1 / 4 : ℝ)) * fpr a 0 (-a) v +
    (((1 / 8 : ℝ) * Real.sqrt 2 + (1 / 8 : ℝ)) * fpr a 0 0 v +
    (((2 : ℝ) * lam) * flin b (-a) v +
    (((-1 / 2 : ℝ)) * fsq b (-a) v +
    (((-2 : ℝ) * lam) * flin 0 (-b) v +
    (((-1 / 2 : ℝ)) * fsq (-a) (-a) v +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (3 / 8 : ℝ)) * fpr b 0 0 v +
    (((-1 / 2 : ℝ)) * fsq b (-b) v +
    (((2 : ℝ) * lam ^ 2) * fcst 0 v +
    (((-1 / 2 : ℝ) * Real.sqrt 2 * lam + (-1 / 2 : ℝ) * lam) * flin (-b) 0 v +
    (((2 : ℝ) * lam) * flin (-a) (-b) v +
    (((2 : ℝ) * lam) * flin b (-b) v +
    (((1 / 4 : ℝ) * Real.sqrt 2 + (-1 / 4 : ℝ)) * fpr b (-a) 0 v +
    (((3 / 4 : ℝ)) * fsq 0 (-a) v +
    (((-3 : ℝ) * lam) * flin 0 (-a) v +
    (((1 / 4 : ℝ)) * fpr (-b) (-a) (-a) v +
    (((1 / 4 : ℝ)) * fpr b (-b) (-a) v +
    (((-1 / 4 : ℝ)) * fsq a (-a) v +
    (((1 / 8 : ℝ) * Real.sqrt 2 + (1 / 8 : ℝ)) * fpr (-b) 0 0 v +
    (((-1 / 2 : ℝ)) * fsq (-a) (-b) v +
    (((1 / 2 : ℝ)) * fsq 0 (-b) v +
    (((1 / 4 : ℝ)) * fpa b v +
    (((1 / 4 : ℝ)) * fpa (-a) v +
    (((-1 / 4 : ℝ)) * fpa 0 v +
    (((-1 / 4 : ℝ)) * fpg b (-a) v +
    (((1 / 4 : ℝ)) * fpg b 0 v +
    (((1 / 4 : ℝ)) * fpg 0 (-a) v)))))))))))))))))))))))))))))))))))))))))))))))))))))

noncomputable def zU0 : ℝ :=
    ((2 : ℝ) * lam) * (FH 0 a - FH 0 b) +
    (((-1 / 2 : ℝ) * Real.sqrt 2 * lam + (-1 / 2 : ℝ) * lam) * (FL P1p a 0 a - FL P1p a 0 b) +
    (((-1 / 4 : ℝ)) * (FG P1p P2a P2a 0 (-b) (-a) a - FG P1p P2a P2a 0 (-b) (-a) b) +
    (((1 / 4 : ℝ)) * (FIV P2a a (-a) a - FIV P2a a (-a) b) +
    (((1 / 8 : ℝ) * Real.sqrt 2 + (1 / 8 : ℝ)) * (FA P1p P2g (-b) 0 a - FA P1p P2g (-b) 0 b) +
    (((1 / 4 : ℝ)) * (FIV P2a b (-a) a - FIV P2a b (-a) b) +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (-1 / 8 : ℝ)) * (FG P1p P2a P2a a (-a) 0 a - FG P1p P2a P2a a (-a) 0 b) +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (3 / 8 : ℝ)) * (FA P1n P2a b 0 a - FA P1n P2a b 0 b) +
    (((-3 : ℝ) * lam ^ 2) * (FK (-a) a - FK (-a) b) +
    (((1 / 2 : ℝ) * Real.sqrt 2 * lam + (-3 / 2 : ℝ) * lam) * (FL P1p (-a) 0 a - FL P1p (-a) 0 b) +
    (((1 / 2 : ℝ) * Real.sqrt 2 * lam + (-3 / 2 : ℝ) * lam) * (FL P1n b 0 a - FL P1n b 0 b) +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (-1 / 8 : ℝ)) * (FG P1p P2g P2a (-b) (-a) 0 a - FG P1p P2g P2a (-b) (-a) 0 b) +
    (((-2 : ℝ) * lam ^ 2) * (FK (-b) a - FK (-b) b) +
    (((-1 / 2 : ℝ)) * (FC 0 a - FC 0 b) +
    (((-1 / 2 : ℝ)) * (FG P1p P2a P2a b 0 (-a) a - FG P1p P2a P2a b 0 (-a) b) +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (-1 / 8 : ℝ)) * (FG P1n P2a P2a b (-b) 0 a - FG P1n P2a P2a b (-b) 0 b) +
    (((-1 / 2 : ℝ)) * (FIV P2a 0 (-a) a - FIV P2a 0 (-a) b) +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (3 / 8 : ℝ)) * (FIV P2g (-a) 0 a - FIV P2g (-a) 0 b) +
    (((1 / 2 : ℝ)) * (FG P1n P2a P2a b (-a) (-b) a - FG P1n P2a P2a b (-a) (-b) b) +
    (((2 : ℝ) * lam) * (FH (-a) a - FH (-a) b) +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (3 / 8 : ℝ)) * (FA P1p P2g (-a) 0 a - FA P1p P2g (-a) 0 b) +
    (((-1 / 2 : ℝ)) * (FG P1p P2a P2a b 0 (-b) a - FG P1p P2a P2a b 0 (-b) b) +
    (((-1 / 2 : ℝ)) * (FG P1n P2a P2a 0 (-a) (-b) a - FG P1n P2a P2a 0 (-a) (-b) b) +
    (((1 : ℝ) * lam) * (FL P1p a (-a) a - FL P1p a (-a) b) +
    (((-1 / 4 : ℝ)) * (FA P1n P2a (-b) (-a) a - FA P1n P2a (-b) (-a) b) +
    (((1 : ℝ) * lam) * (FL P1n (-b) (-a) a - FL P1n (-b) (-a) b) +
    (((1 / 8 : ℝ) * Real.sqrt 2 + (1 / 8 : ℝ)) * (FA P1p P2a a 0 a - FA P1p P2a a 0 b) +
    (((-1 / 4 : ℝ)) * (FG P1n P2a P2a a 0 (-a) a - FG P1n P2a P2a a 0 (-a) b) +
    (((1 / 8 : ℝ) * Real.sqrt 2 + (1 / 8 : ℝ)) * (FIV P2a a 0 a - FIV P2a a 0 b) +
    (((2 : ℝ) * lam) * (FL P1n b (-a) a - FL P1n b (-a) b) +
    (((-1 / 2 : ℝ)) * (FA P1n P2a b (-a) a - FA P1n P2a b (-a) b) +
    (((-2 : ℝ) * lam) * (FL P1n 0 (-b) a - FL P1n 0 (-b) b) +
    (((-1 / 2 : ℝ)) * (FC (-a) a - FC (-a) b) +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (3 / 8 : ℝ)) * (FIV P2a b 0 a - FIV P2a b 0 b) +
    (((-1 / 2 : ℝ)) * (FA P1n P2a b (-b) a - FA P1n P2a b (-b) b) +
    (((2 : ℝ) * lam ^ 2) * (FK 0 a - FK 0 b) +
    (((-1 / 2 : ℝ) * Real.sqrt 2 * lam + (-1 / 2 : ℝ) * lam) * (FL P1p (-b) 0 a - FL P1p (-b) 0 b) +
    (((2 : ℝ) * lam) * (FL P1p (-a) (-b) a - FL P1p (-a) (-b) b) +
    (((2 : ℝ) * lam) * (FL P1n b (-b) a - FL P1n b (-b) b) +
    (((1 / 4 : ℝ) * Real.sqrt 2 + (-1 / 4 : ℝ)) * (FG P1n P2a P2a b (-a) 0 a - FG P1n P2a P2a b (-a) 0 b) +
    (((3 / 4 : ℝ)) * (FA P1n P2a 0 (-a) a - FA P1n P2a 0 (-a) b) +
    (((-3 : ℝ) * lam) * (FL P1n 0 (-a) a - FL P1n 0 (-a) b) +
    (((1 / 4 : ℝ)) * (FIV P2a (-b) (-a) a - FIV P2a (-b) (-a) b) +
    (((1 / 4 : ℝ)) * (FG P1p P2a P2a b (-b) (-a) a - FG P1p P2a P2a b (-b) (-a) b) +
    (((-1 / 4 : ℝ)) * (FA P1p P2a a (-a) a - FA P1p P2a a (-a) b) +
    (((1 / 8 : ℝ) * Real.sqrt 2 + (1 / 8 : ℝ)) * (FIV P2g (-b) 0 a - FIV P2g (-b) 0 b) +
    (((-1 / 2 : ℝ)) * (FA P1p P2g (-a) (-b) a - FA P1p P2g (-a) (-b) b) +
    (((1 / 2 : ℝ)) * (FA P1n P2a 0 (-b) a - FA P1n P2a 0 (-b) b) +
    (((1 / 4 : ℝ)) * (FPA b a - FPA b b) +
    (((1 / 4 : ℝ)) * (FPA (-a) a - FPA (-a) b) +
    (((-1 / 4 : ℝ)) * (FPA 0 a - FPA 0 b) +
    (((-1 / 4 : ℝ)) * (FPG b (-a) a - FPG b (-a) b) +
    (((1 / 4 : ℝ)) * (FPG b 0 a - FPG b 0 b) +
    (((1 / 4 : ℝ)) * (FPG 0 (-a) a - FPG 0 (-a) b))))))))))))))))))))))))))))))))))))))))))))))))))))))

lemma zS0_int : IntervalIntegrable zS0 volume b a ∧ ∫ v in b..a, zS0 v = zU0 := by
  have h := (lin_cons ((2 : ℝ) * lam) zb0
    (lin_cons ((-1 / 2 : ℝ) * Real.sqrt 2 * lam + (-1 / 2 : ℝ) * lam) zb1
    (lin_cons ((-1 / 4 : ℝ)) zb2
    (lin_cons ((1 / 4 : ℝ)) zb3
    (lin_cons ((1 / 8 : ℝ) * Real.sqrt 2 + (1 / 8 : ℝ)) zb4
    (lin_cons ((1 / 4 : ℝ)) zb5
    (lin_cons ((-1 / 8 : ℝ) * Real.sqrt 2 + (-1 / 8 : ℝ)) zb6
    (lin_cons ((-1 / 8 : ℝ) * Real.sqrt 2 + (3 / 8 : ℝ)) zb7
    (lin_cons ((-3 : ℝ) * lam ^ 2) zb8
    (lin_cons ((1 / 2 : ℝ) * Real.sqrt 2 * lam + (-3 / 2 : ℝ) * lam) zb9
    (lin_cons ((1 / 2 : ℝ) * Real.sqrt 2 * lam + (-3 / 2 : ℝ) * lam) zb10
    (lin_cons ((-1 / 8 : ℝ) * Real.sqrt 2 + (-1 / 8 : ℝ)) zb11
    (lin_cons ((-2 : ℝ) * lam ^ 2) zb12
    (lin_cons ((-1 / 2 : ℝ)) zb13
    (lin_cons ((-1 / 2 : ℝ)) zb14
    (lin_cons ((-1 / 8 : ℝ) * Real.sqrt 2 + (-1 / 8 : ℝ)) zb15
    (lin_cons ((-1 / 2 : ℝ)) zb16
    (lin_cons ((-1 / 8 : ℝ) * Real.sqrt 2 + (3 / 8 : ℝ)) zb17
    (lin_cons ((1 / 2 : ℝ)) zb18
    (lin_cons ((2 : ℝ) * lam) zb19
    (lin_cons ((-1 / 8 : ℝ) * Real.sqrt 2 + (3 / 8 : ℝ)) zb20
    (lin_cons ((-1 / 2 : ℝ)) zb21
    (lin_cons ((-1 / 2 : ℝ)) zb22
    (lin_cons ((1 : ℝ) * lam) zb23
    (lin_cons ((-1 / 4 : ℝ)) zb24
    (lin_cons ((1 : ℝ) * lam) zb25
    (lin_cons ((1 / 8 : ℝ) * Real.sqrt 2 + (1 / 8 : ℝ)) zb26
    (lin_cons ((-1 / 4 : ℝ)) zb27
    (lin_cons ((1 / 8 : ℝ) * Real.sqrt 2 + (1 / 8 : ℝ)) zb28
    (lin_cons ((2 : ℝ) * lam) zb29
    (lin_cons ((-1 / 2 : ℝ)) zb30
    (lin_cons ((-2 : ℝ) * lam) zb31
    (lin_cons ((-1 / 2 : ℝ)) zb32
    (lin_cons ((-1 / 8 : ℝ) * Real.sqrt 2 + (3 / 8 : ℝ)) zb33
    (lin_cons ((-1 / 2 : ℝ)) zb34
    (lin_cons ((2 : ℝ) * lam ^ 2) zb35
    (lin_cons ((-1 / 2 : ℝ) * Real.sqrt 2 * lam + (-1 / 2 : ℝ) * lam) zb36
    (lin_cons ((2 : ℝ) * lam) zb37
    (lin_cons ((2 : ℝ) * lam) zb38
    (lin_cons ((1 / 4 : ℝ) * Real.sqrt 2 + (-1 / 4 : ℝ)) zb39
    (lin_cons ((3 / 4 : ℝ)) zb40
    (lin_cons ((-3 : ℝ) * lam) zb41
    (lin_cons ((1 / 4 : ℝ)) zb42
    (lin_cons ((1 / 4 : ℝ)) zb43
    (lin_cons ((-1 / 4 : ℝ)) zb44
    (lin_cons ((1 / 8 : ℝ) * Real.sqrt 2 + (1 / 8 : ℝ)) zb45
    (lin_cons ((-1 / 2 : ℝ)) zb46
    (lin_cons ((1 / 2 : ℝ)) zb47
    (lin_cons ((1 / 4 : ℝ)) zb48
    (lin_cons ((1 / 4 : ℝ)) zb49
    (lin_cons ((-1 / 4 : ℝ)) zb50
    (lin_cons ((-1 / 4 : ℝ)) zb51
    (lin_cons ((1 / 4 : ℝ)) zb52
    (lin_one ((1 / 4 : ℝ)) zb53))))))))))))))))))))))))))))))))))))))))))))))))))))))
  exact h

noncomputable def zS1 (v : ℝ) : ℝ :=
    ((1 / 4 : ℝ)) * fpr b a (-a) v +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (-1 / 8 : ℝ)) * fpr b a 0 v)

noncomputable def zU1 : ℝ :=
    ((1 / 4 : ℝ)) * (FG P1n P2a P2a b a (-a) 1 - FG P1n P2a P2a b a (-a) b) +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (-1 / 8 : ℝ)) * (FG P1n P2a P2a b a 0 1 - FG P1n P2a P2a b a 0 b))

lemma zS1_int : IntervalIntegrable zS1 volume b 1 ∧ ∫ v in b..1, zS1 v = zU1 := by
  have h := (lin_cons ((1 / 4 : ℝ)) zb54
    (lin_one ((-1 / 8 : ℝ) * Real.sqrt 2 + (-1 / 8 : ℝ)) zb55))
  exact h

noncomputable def zS2 (v : ℝ) : ℝ :=
    ((1 / 4 : ℝ)) * fpr a b (-a) v +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (-1 / 8 : ℝ)) * fpr a b 0 v)

noncomputable def zU2 : ℝ :=
    ((1 / 4 : ℝ)) * (FG P1n P2a P2a a b (-a) a - FG P1n P2a P2a a b (-a) 1) +
    (((-1 / 8 : ℝ) * Real.sqrt 2 + (-1 / 8 : ℝ)) * (FG P1n P2a P2a a b 0 a - FG P1n P2a P2a a b 0 1))

lemma zS2_int : IntervalIntegrable zS2 volume 1 a ∧ ∫ v in 1..a, zS2 v = zU2 := by
  have h := (lin_cons ((1 / 4 : ℝ)) zb56
    (lin_one ((-1 / 8 : ℝ) * Real.sqrt 2 + (-1 / 8 : ℝ)) zb57))
  exact h


lemma la_b : Real.log (a - b) = lam := by
  rw [show a - b = 2 by simp only [Y01.a, Y01.b]; ring]; rfl

lemma la_ma : Real.log (a - -a) = lam + L := by
  rw [show a - -a = 2 * a by ring, Real.log_mul (by norm_num) nza]; rfl

lemma la_0 : Real.log (a - 0) = L := by rw [sub_zero]; rfl

lemma zpw1 : ∀ v ∈ Ioo b 1, Z05.W3pintegrand v = zS0 v + zS1 v := by
  intro v hv
  have hv' : v ∈ Ioo b a := ⟨hv.1, lt_of_lt_of_le hv.2 h1a⟩
  unfold Z05.W3pintegrand Z05.Om20 Z05.Om11 Y01.l1 Y01.l2
  rw [log_tv hv', log_one_sub_tv hv']
  unfold zS0 zS1 fsq fpr flin fcst fpa fpg
  simp only [la_b, la_ma, la_0]
  simp only [div_eq_mul_inv, mul_inv, sub_zero, sub_neg_eq_add]
  ring

lemma zpw2 : ∀ v ∈ Ioo 1 a, Z05.W3pintegrand v = zS0 v + zS2 v := by
  intro v hv
  have hv' : v ∈ Ioo b a := ⟨lt_of_le_of_lt hb1 hv.1, hv.2⟩
  unfold Z05.W3pintegrand Z05.Om20 Z05.Om11 Y01.l1 Y01.l2
  rw [log_tv hv', log_one_sub_tv hv']
  unfold zS0 zS2 fsq fpr flin fcst fpa fpg
  simp only [la_b, la_ma, la_0]
  simp only [div_eq_mul_inv, mul_inv, sub_zero, sub_neg_eq_add]
  ring

theorem main8 : Z05.W3pint = zU0 + zU1 + zU2 := by
  obtain ⟨i0, e0⟩ := zS0_int
  obtain ⟨i1, e1⟩ := zS1_int
  obtain ⟨i2, e2⟩ := zS2_int
  have i0a : IntervalIntegrable zS0 volume b 1 := i0.mono_set (by
    rw [uIcc_of_le hba, uIcc_of_le hb1]; exact Icc_subset_Icc le_rfl h1a)
  have i0b : IntervalIntegrable zS0 volume 1 a := i0.mono_set (by
    rw [uIcc_of_le hba, uIcc_of_le h1a]; exact Icc_subset_Icc hb1 le_rfl)
  have k1 := congr_Ioo (A := (∫ v in b..1, zS0 v) + zU1) hb1 zpw1
    ⟨i0a.add i1, by rw [intervalIntegral.integral_add i0a i1, e1]⟩
  have k2 := congr_Ioo (A := (∫ v in (1:ℝ)..a, zS0 v) + zU2) h1a zpw2
    ⟨i0b.add i2, by rw [intervalIntegral.integral_add i0b i2, e2]⟩
  have hs := intervalIntegral.integral_add_adjacent_intervals i0a i0b
  have hw := intervalIntegral.integral_add_adjacent_intervals k1.1 k2.1
  unfold Z05.W3pint
  rw [← hw, k1.2, k2.2]
  rw [← e0, ← hs]; ring

lemma t_ma2 : Li3 (-(a ^ 2)) = Li3 (-(b ^ 2)) - π ^ 2 / 6 * (2 * L) - 1 / 6 * (2 * L) ^ 3 := by
  have h := Y01.PE2.W_I3 (z := -(a ^ 2)) (w := -(b ^ 2)) (c := 2 * L) (by nlinarith [Y01.PE2.a_pos])
    (by rw [one_div, inv_neg, ← inv_pow, ← Y01.PE2.b_eq_inv]) (by rw [neg_neg, Y01.PE2.log_a2])
  linarith

lemma t_ms : Li3 (-Real.sqrt 2)
    = Li3 (-(1 / Real.sqrt 2)) - π ^ 2 / 6 * (lam / 2) - 1 / 6 * (lam / 2) ^ 3 := by
  have h := Y01.PE2.W_I3 (z := -Real.sqrt 2) (w := -(1 / Real.sqrt 2)) (c := lam / 2)
    (by linarith [Y01.PE2.hs0]) (by rw [one_div, one_div, inv_neg]) (by rw [neg_neg, Y01.PE2.log_s])
  linarith

set_option maxHeartbeats 4000000 in
lemma final8 : zU0 + zU1 + zU2 = Z05.W3pclosed := by
  unfold zU0 zU1 zU2
  simp only [zv0, zv1, zv2, zv3, zv4, zv5, zv6, zv7, zv8, zv9, zv10, zv11, zv12, zv13, zv14, zv15, zv16, zv17, zv18, zv19, zv20, zv21, zv22, zv23, zv24, zv25, zv26, zv27, zv28, zv29, zv30, zv31, zv32, zv33, zv34, zv35, zv36, zv37, zv38, zv39, zv40, zv41, zv42, zv43, zv44, zv45, zv46, zv47, zv48, zv49, zv50, zv51, zv52, zv53, zv54, zv55, zv56, zv57, zv58, zv59, zv60, zv61, zv62, zv63, zv64, zv65, zv66, zv67, zv68, zv69, zv70, zv71, zv72, zv73, zv74, zv75, zv76, zv77, zv78, zv79, zv80, zv81, zv82, zv83, zv84, zv85, zv86, zv87, zv88, zv89, zv90, zv91, zv92, zv93, zv94, zv95, zv96, zv97, zv98, zv99, zv100, zv101, zv102, zv103, zv104, zv105, zv106, zv107, zv108, zv109, zv110, zv111, zv112, zv113, zv114, zv115]
  simp only [t_ma2, t_ms, Y01.PE2.t_m2a, Y01.PE2.t_ma, Y01.PE2.t_m1, Y01.PE2.t_mis, Y01.PE2.t_mbh, Y01.PE2.t_mb2, Y01.PE2.t_m4, Y01.PE2.t_1mis, Y01.PE2.t_h, Y01.PE2.t_2ms, Y01.PE2.t_2b, Y01.PE2.t_p4, Y01.PE2.t_b2, Y01.PE2.d_ma2, Y01.PE2.d_m2a, Y01.PE2.d_ma, Y01.PE2.d_m1, Y01.PE2.d_mis, Y01.PE2.d_mb, Y01.PE2.d_mbh, Y01.PE2.d_mb2, Y01.PE2.d_m4, Y01.PE2.d_b2, Y01.PE2.d_1mis, Y01.PE2.d_h, Y01.PE2.d_2ms, Y01.PE2.d_is, Y01.PE2.d_2b, Y01.PE2.d_p4, Y01.PE2.d_ms]
  unfold Z05.W3pclosed
  ring

theorem main : Z05.stmt_8 := by
  unfold Z05.stmt_8
  rw [main8, final8]

end Z05.P8

theorem Z05.proof_1 : Z05.stmt_1 := Z05.P1.main
theorem Z05.proof_2 : Z05.stmt_2 := Z05.P2.main
theorem Z05.proof_3 : Z05.stmt_3 := Z05.P3.main Z05.proof_2
theorem Z05.proof_4 : Z05.stmt_4 := Z05.P4.main
theorem Z05.proof_5 : Z05.stmt_5 := Z05.P5.main Z05.proof_3
theorem Z05.proof_6 : Z05.stmt_6 := Z05.P6.main
theorem Z05.proof_7 : Z05.stmt_7 := Z05.P7.main Z05.proof_5
theorem Z05.proof_8 : Z05.stmt_8 := Z05.P8.main
theorem Z05.proof_9 : Z05.stmt_9 := Z05.P9.main Z05.proof_7 Z05.proof_8
theorem Z05.proof_10 : Z05.stmt_10 := Z05.P10.main
theorem Z05.proof_11 : Z05.stmt_11 := Z05.P12.main Z05.proof_9
