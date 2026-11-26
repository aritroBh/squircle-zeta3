import Mathlib
import Zeta3.Wave5.Y02_P2
/-!
# Z06 -- the next order of the corner asymptotics at every angle

Target id: Z06.
Paper labels: the paragraph after `thm:P2asym` (paper/ieee/main.tex line 1462) and the paragraph
after it (line 1464, the values at `θ = π/2` and the first-order constant `C(θ)`); detailed
proofs in paper/ieee/proofs/sec_g3_P2asym.tex: Lemma `g3:prop:unif` (lines 32-34, proof lines
140-164), Lemma `g3:prop:n3` (lines 170-186), Corollary `g3:cor:C3` (lines 285-305),
Lemma `g3:prop:n4` (lines 190-283, Step 0(iv),(v) at lines 222-226) and Lemma `g3:lem:a4`
(lines 309-353).

Everything about the corner is reused from `Zeta3.Wave5.Y02_P2` (not redefined):
`Y02.L θ n` (the length `L_θ(n)`), `Y02.zz θ = z = sin²(θ/2)`, `Y02.ww = w`, `Y02.ell = ℓ`,
`Y02.U = U_θ`, `Y02.a1`, `Y02.a2`, `Y02.a3`, `Y02.S12` (the closed form of `S_{1,2}(z)`),
`Y02.mu = μ_θ`, `Y02.c2 = c_2(θ)`, `Y02.R θ ε = (πε/sin πε) Λ_ε(θ)`, and
`Y02.E θ n = L_θ(n) - [1 + cos θ + 2z(μ_θ/n)^{1/(n-1)} + c_2(θ)/n²]`.
`ζ(3)` is `Y02.z3 = (riemannZeta 3).re`, the real part of Mathlib's `riemannZeta 3`
(which is real and equals `∑_{n≥1} n^{-3}` by `zeta_nat_eq_tsum_of_gt_one`; Y02 proves
`HasSum (fun n => 1/(n+1)^3) z3`). `ζ(2) = π²/6`, `ζ(4) = π⁴/90`, `Li_2 = D00.Li2`,
`Li_3 = Y02.Li3`.

Quoted verbatim (main.tex line 1462):

  The remainder is uniform in the angle: for $n\ge9$ and every $0<\theta<\pi$ it is at most
  $C\sin^2\frac\theta2\,\log n/n^3$ with an absolute constant $C$, and its leading part is
  $-2z\,a_2\log n/n^3$. Since $a_2>0$ at every angle, the logarithm is really there (proved;
  refereed once).

Quoted verbatim (main.tex line 1464):

  At $\theta=\pi/2$ we get $\mu_{\pi/2}=2/e$ and $c_2(\frac\pi2)=\pi^2/12$, which recovers
  Ahmed's Corollary~1(b)~\cite{ahmed}, $L(n)=1+(2/(en))^{1/(n-1)}+\pi^2/(12n^2)+O(\log n/n^3)$.
  The constant in the first-order defect is
  $C(\theta)=-2z\,a_1=2\sin^2\frac\theta2+4\cos^2\frac\theta2\log\cos\frac\theta2$, which
  recovers Ahmed's constant~\cite[Thm.\ on corners of arbitrary angle]{ahmed}. ...
  Notice where $\zeta(3)$ enters: in $a_3$, through $S_{1,2}$, at every angle. At
  $\theta=\pi/2$, $a_1$, $a_2$ and $a_3=\frac{\zeta(3)}4-\frac13$ are Ahmed's coefficients

Quoted verbatim (sec_g3_P2asym.tex, Lemma g3:prop:unif):

  For $n\ge9$ and every $0<\theta<\pi$ the remainder
  $E_n=L_\theta(n)-\bigl[1+\cos\theta+2z(\mu_\theta/n)^{\varepsilon}+c_2(\theta)/n^2\bigr]$
  satisfies $|E_n|\le C\,\sin^2\frac\theta2\,\frac{\log n}{n^3}$ with the absolute constant
  $C=4100$, and $E_n=-2z\,a_2\,\frac{\log n}{n^3}+\Delta_n$ with $|\Delta_n|\le9000\,z/n^3$.
  Moreover $a_2(\theta)>0$ for every $0<\theta<\pi$, so the logarithm is really there.

  (proof, line 163:) ... $T_2,T_3,T_4$ satisfy
  $|T_2|+|T_3|+|T_4|\le(8870+5.2+1.4)z\,n^{-3}\le8880\,z\,n^{-3}$, without a logarithm. ...
  $T_1=-2za_2\frac{\log n}{n^3}+2za_2\bigl[(\varepsilon^2-n^{-2})-e'\bigr]$ with
  $\bigl|2za_2[\ldots]\bigr|\le2z\cdot2.15\cdot6.4\,n^{-3}\le28\,z\,n^{-3}$.
  (Here $T_1=2za_2(q\varepsilon^2-n^{-2})$, $q=(\mu_\theta/n)^\varepsilon$, line 146.)

Quoted verbatim (sec_g3_P2asym.tex, Lemma g3:prop:n3 and Corollary g3:cor:C3):

  Put $d(\theta)=2z\bigl[a_2(2+a_1)+a_3\bigr]-z(3-z-2z^2)$. There is an absolute constant $C_3$
  such that for $n\ge9$ and $0<\theta<\pi$,
  \[ \Bigl|E_n-\frac{d(\theta)-2z\,a_2\log n}{n^3}\Bigr|\le C_3\,z\,\frac{\log^2n}{n^4}. \]
  Lemma~\ref{g3:prop:n3} holds with $C_3=3.1\cdot10^4$.

Quoted verbatim (sec_g3_P2asym.tex, Lemma g3:prop:n4 and Step 0(iv),(v)):

  Put $\lambda_n=\log n-a_1$, $p_4(z)=-\frac16z(1-z)(20z^2-2z+43)$ and
  \[ Q_4(\theta,n)=z\Bigl[a_2\bigl((\lambda_n-3)^2-3\bigr)-2a_3(\lambda_n-3)+a_2^2+2a_4\Bigr]+p_4(z). \]
  For $n\ge9$ and $0<\theta<\pi$ let $L_n=\log n$ and
  $E^{(4)}_n=E_n-\frac{d(\theta)-2za_2\log n}{n^3}-\frac{Q_4(\theta,n)}{n^4}$, with $d$ as in
  Lemma~\ref{g3:prop:n3}. Then
  \[ \text{(a)}\quad E^{(4)}_n=-\frac{z\,a_2}{3}\,\frac{\log^3n}{n^5}+\Delta'_n,\quad
     |\Delta'_n|\le2.5\cdot10^5\,z\,\frac{\log^2n}{n^5};\qquad
     \text{(b)}\quad |E^{(4)}_n|\le1.2\cdot10^5\,z\,\frac{\log^3n}{n^5}. \]
  (iv) ... $H=1+a_2\varepsilon^2+a_3\varepsilon^3+b_4\varepsilon^4+\cdots$ with
  $b_4=a_4+\frac12a_2^2$. Cauchy's estimate ... gives $|a_3|\le3.1\cdot5^3=387.5$ and
  $|b_4|\le3.1\cdot5^4=1937.5$, for every $z\in(0,1)$.
  (v) ... for $0\le\varepsilon\le\frac18$,
  $H(\varepsilon)=1+a_2\varepsilon^2+a_3\varepsilon^3+b_4\varepsilon^4+\rho_5(\varepsilon)$,
  $|\rho_5(\varepsilon)|\le3.1\cdot10^5\,\varepsilon^5$.
  (Here $H(\varepsilon)=g(\varepsilon,z)e^{-a_1\varepsilon}$ and $R_\theta(\varepsilon)=2z\,g$.)

Quoted verbatim (sec_g3_P2asym.tex, Lemma g3:lem:a4):

  For $0<z<1$, with $S_{2,2}(z)=\int_0^z\frac{S_{1,2}(t)}t\dd t$:
  $a_4=\frac12\zeta(4)-2S_{2,2}(z)+2w^2\ell\,S_{1,2}(z)+\frac14(1-U_\theta^4)$.

## Conventions

* `E θ n` is Y02's `E_n` with `ε = 1/(n-1)` as a real `rpow` of the positive base `μ_θ/n`.
  All statements are for real `n ≥ 9` (the supplement's `n` is real, Step 0 of g3:prop:n4).
* `S_{1,2}(t)` inside `S_{2,2}` is needed at every `t ∈ (0,z)`, so here it is the integral
  `S12int t = ½ ∫_0^t log²(1-s)/s ds` (the supplement's own characterization, g3:lem:sums(iv);
  Y02 stmt_10 proves it equals the closed form `Y02.S12`). The integrand is continuous on
  `(0,t]` and tends to `0` at `0` (Lean's value `0/0 = 0` at the single point `s = 0` is
  irrelevant), so both interval integrals are genuine for `0 ≤ t < 1`, and `S12int s / s` is
  bounded near `0` (`S12int s ~ s²/4`).
* `a_4` is the supplement's closed form; stmt_6 asserts it is the coefficient that makes the
  order-5 Taylor estimate of `R_θ(ε) e^{-a_1ε}` true (with `b_4 = a_4 + a_2²/2`), which pins it
  down as the `ε⁴` Taylor coefficient of `log(R_θ/2z)`.
* The Taylor statements use `0 < ε` (not `0 ≤ ε`), because Lean's `Y02.R θ 0 = 0` (the factor
  `πε/sin πε` is `0/0` there), while the paper's value is the limit `2z`.
* O-claims are stated as explicit bounds with the supplement's constants.

## Statements

* stmt_1: at `θ = π/2`: `a_1 = log 2 - 1`, `μ = 2/e`, `a_2 = π²/12 + 1/2`, `c_2 = π²/12`.
* stmt_2: at `θ = π/2`: `S_{1,2}(1/2) = ζ(3)/8 - (log 2)³/6`.
* stmt_3: at `θ = π/2`: `a_3 = ζ(3)/4 - 1/3`.
* stmt_4: `C(θ) = -2z a_1 = 2 sin²(θ/2) + 4 cos²(θ/2) log cos(θ/2)` for `0 < θ < π`.
* stmt_5: Cauchy bounds (Step 0(iv)): `|a_3| ≤ 387.5`, `|a_4 + a_2²/2| ≤ 1937.5`.
* stmt_6: order-5 Taylor (Step 0(v)): `|R e^{-a_1ε} - 2z(1+a_2ε²+a_3ε³+b_4ε⁴)| ≤ 2z·3.1·10⁵ ε⁵`.
* stmt_7: the `T_1` expansion: for `-1 < a < 0`, `n ≥ 9`,
  `|(e^a/n)^ε ε² - n^{-2} + log n/n³| ≤ 6.4/n³` (proof of g3:prop:unif, second claim).
* stmt_8: the logarithm-free part: `|E_n - 2z a_2((μ/n)^ε ε² - n^{-2})| ≤ 8880 z/n³`.
* stmt_9: main.tex 1462 / g3:prop:unif second claim:
  `|E_n + 2z a_2 log n/n³| ≤ 9000 z/n³`, `n ≥ 9`, `0 < θ < π`.
* stmt_10: g3:prop:n3 with g3:cor:C3: `|E_n - (d - 2z a_2 log n)/n³| ≤ 3.1·10⁴ z log²n/n⁴`.
* stmt_11: g3:prop:n4(b): `|E⁽⁴⁾_n| ≤ 1.2·10⁵ z log³n/n⁵`.
* stmt_12: g3:prop:n4(a): `|E⁽⁴⁾_n + (z a_2/3) log³n/n⁵| ≤ 2.5·10⁵ z log²n/n⁵`.

Numerical check (mpmath, 45 digits, scratch script): stmt_1-4 hold to `1e-45`; for
`θ ∈ {0.05, 0.5, π/3, π/2, 2π/3, 2.7, π-0.01}` and `n ∈ {9, 10, 30, 100, 10³, 10⁵}` the
scaled quantities of stmt_8, 9, 10, 12, 11 are at most `3.6, 3.6, 1.9, 4.7, 1.4` (bounds
`8880, 9000, 31000, 2.5e5, 1.2e5`); stmt_7's scaled quantity is at most `2.01` (bound `6.4`);
in stmt_5/6 the observed `|a_3| ≤ 0.34`, `|b_4| ≤ 1.9`, `|ρ_5|/ε⁵ ≤ 0.37`. `E_n` was computed
from the exact series of thm:P2 (checked against direct quadrature of `L_θ(n)` to `1e-16`).

## NOT COVERED

* The first claim of line 1462 (`|E_n| ≤ C sin²(θ/2) log n/n³`) and `a_2 > 0`: proved in Y02
  (stmt_13, stmt_9), not repeated.
* "which recovers Ahmed's Corollary 1(b)", "agrees with every constant in Table 4 of [ahmed]",
  "are Ahmed's coefficients": comparisons with an external paper.
* "the logarithm is really there" (that no `O(1/n³)` bound holds) and the analogous last
  sentence of g3:prop:n4 (`log³ n` cannot be lowered): consequences of stmt_9/stmt_12 with
  `a_2 > 0`, not stated separately.
* The last sentences of the paragraph at line 1464 (the density `h_c` and the squircle corners
  `-2M_0(π-θ) = ∫_0^∞ h_c(s) ds/(s(1+s))`): a different subsection's object.
* The decomposition `E_n = T_1 + T_2 + T_3 + T_4` and `P_1 + ... + P_6` as identities (only their
  consequences stmt_8 and the final bounds are stated).
-/
-- BEGIN STATEMENTS
namespace Z06

open Real

/-- `S_{1,2}(t) = ½ ∫_0^t log²(1-s)/s ds` (g3:lem:sums(iv)), as a function of `t`. -/
noncomputable def S12int (t : ℝ) : ℝ :=
  (1 / 2) * ∫ s in (0 : ℝ)..t, Real.log (1 - s) ^ 2 / s

/-- `S_{2,2}(z) = ∫_0^z S_{1,2}(t)/t dt` (g3:lem:a4). -/
noncomputable def S22 (z : ℝ) : ℝ := ∫ t in (0 : ℝ)..z, S12int t / t

/-- `a_4 = ½ζ(4) - 2S_{2,2}(z) + 2w²ℓ S_{1,2}(z) + ¼(1 - U_θ⁴)`, with `ζ(4) = π⁴/90`. -/
noncomputable def a4 (θ : ℝ) : ℝ :=
  π ^ 4 / 180 - 2 * S22 (Y02.zz θ) + 2 * Y02.ww θ ^ 2 * Y02.ell θ * Y02.S12 θ
    + (1 - Y02.U θ ^ 4) / 4

/-- `d(θ) = 2z[a_2(2 + a_1) + a_3] - z(3 - z - 2z²)` (g3:prop:n3). -/
noncomputable def d (θ : ℝ) : ℝ :=
  2 * Y02.zz θ * (Y02.a2 θ * (2 + Y02.a1 θ) + Y02.a3 θ)
    - Y02.zz θ * (3 - Y02.zz θ - 2 * Y02.zz θ ^ 2)

/-- `λ_n = log n - a_1` (g3:prop:n4). -/
noncomputable def lam (θ n : ℝ) : ℝ := Real.log n - Y02.a1 θ

/-- `p_4(z) = -⅙ z(1-z)(20z² - 2z + 43)`. -/
noncomputable def p4 (z : ℝ) : ℝ := -(1 / 6) * z * (1 - z) * (20 * z ^ 2 - 2 * z + 43)

/-- `Q_4(θ,n) = z[a_2((λ_n-3)² - 3) - 2a_3(λ_n-3) + a_2² + 2a_4] + p_4(z)` (g3:prop:n4). -/
noncomputable def Q4 (θ n : ℝ) : ℝ :=
  Y02.zz θ * (Y02.a2 θ * ((lam θ n - 3) ^ 2 - 3) - 2 * Y02.a3 θ * (lam θ n - 3)
      + Y02.a2 θ ^ 2 + 2 * a4 θ) + p4 (Y02.zz θ)

/-- `E⁽⁴⁾_n = E_n - (d(θ) - 2z a_2 log n)/n³ - Q_4(θ,n)/n⁴` (g3:prop:n4). -/
noncomputable def E4 (θ n : ℝ) : ℝ :=
  Y02.E θ n - (d θ - 2 * Y02.zz θ * Y02.a2 θ * Real.log n) / n ^ 3 - Q4 θ n / n ^ 4

/-- main.tex 1464: at `θ = π/2`, `a_1 = log 2 - 1`, `μ_{π/2} = 2/e`, `a_2 = π²/12 + 1/2`
(Ahmed's coefficients) and `c_2(π/2) = π²/12`. -/
def stmt_1 : Prop :=
  Y02.a1 (π / 2) = Real.log 2 - 1 ∧ Y02.mu (π / 2) = 2 / Real.exp 1 ∧
    Y02.a2 (π / 2) = π ^ 2 / 12 + 1 / 2 ∧ Y02.c2 (π / 2) = π ^ 2 / 12

/-- The value of the Nielsen polylogarithm behind `a_3(π/2)`:
`S_{1,2}(1/2) = ζ(3)/8 - (log 2)³/6`. -/
def stmt_2 : Prop :=
  Y02.S12 (π / 2) = Y02.z3 / 8 - Real.log 2 ^ 3 / 6

/-- main.tex 1464: "At $\theta=\pi/2$, ... $a_3=\frac{\zeta(3)}4-\frac13$". -/
def stmt_3 : Prop :=
  Y02.a3 (π / 2) = Y02.z3 / 4 - 1 / 3

/-- main.tex 1464: the first-order constant
`C(θ) = -2z a_1 = 2 sin²(θ/2) + 4 cos²(θ/2) log cos(θ/2)`, `0 < θ < π`. -/
def stmt_4 : Prop :=
  ∀ θ : ℝ, 0 < θ → θ < π →
    -2 * Y02.zz θ * Y02.a1 θ
      = 2 * Real.sin (θ / 2) ^ 2 + 4 * Real.cos (θ / 2) ^ 2 * Real.log (Real.cos (θ / 2))

/-- g3:prop:n4 Step 0(iv): Cauchy bounds `|a_3| ≤ 387.5` and `|b_4| = |a_4 + a_2²/2| ≤ 1937.5`,
uniformly in `0 < θ < π`. -/
def stmt_5 : Prop :=
  ∀ θ : ℝ, 0 < θ → θ < π →
    |Y02.a3 θ| ≤ 387.5 ∧ |a4 θ + Y02.a2 θ ^ 2 / 2| ≤ 1937.5

/-- g3:prop:n4 Step 0(v) with g3:lem:a4 (the `a_4` above is the next Taylor coefficient):
`R_θ(ε) e^{-a_1ε} = 2z(1 + a_2ε² + a_3ε³ + b_4ε⁴ + ρ_5)`, `b_4 = a_4 + a_2²/2`,
`|ρ_5| ≤ 3.1·10⁵ ε⁵`, for `0 < ε ≤ 1/8`. -/
def stmt_6 : Prop :=
  ∀ θ ε : ℝ, 0 < θ → θ < π → 0 < ε → ε ≤ 1 / 8 →
    |Y02.R θ ε * Real.exp (-(Y02.a1 θ * ε))
        - 2 * Y02.zz θ * (1 + Y02.a2 θ * ε ^ 2 + Y02.a3 θ * ε ^ 3
            + (a4 θ + Y02.a2 θ ^ 2 / 2) * ε ^ 4)|
      ≤ 2 * Y02.zz θ * (3.1e5 * ε ^ 5)

/-- Proof of g3:prop:unif (second claim), the `T_1` step: for `-1 < a < 0` (the range of `a_1`)
and `n ≥ 9`, with `q = (e^a/n)^ε`, `ε = 1/(n-1)`: `q ε² - n^{-2} = -log n/n³ + O(1/n³)`,
precisely `|q ε² - n^{-2} + log n/n³| ≤ 6.4/n³` (`(ε² - n^{-2}) - e'` with `2.4 + 4`). -/
def stmt_7 : Prop :=
  ∀ a n : ℝ, -1 < a → a < 0 → 9 ≤ n →
    |(Real.exp a / n) ^ (1 / (n - 1)) * (1 / (n - 1)) ^ 2 - 1 / n ^ 2 + Real.log n / n ^ 3|
      ≤ 6.4 / n ^ 3

/-- Proof of g3:prop:unif (second claim): `|T_2| + |T_3| + |T_4| ≤ 8880 z/n³`, i.e.
`|E_n - 2z a_2 (q ε² - n^{-2})| ≤ 8880 z/n³`, `q = (μ_θ/n)^ε`, for `n ≥ 9`, `0 < θ < π`. -/
def stmt_8 : Prop :=
  ∀ θ n : ℝ, 0 < θ → θ < π → 9 ≤ n →
    |Y02.E θ n - 2 * Y02.zz θ * Y02.a2 θ *
        ((Y02.mu θ / n) ^ (1 / (n - 1)) * (1 / (n - 1)) ^ 2 - 1 / n ^ 2)|
      ≤ 8880 * Y02.zz θ / n ^ 3

/-- main.tex line 1462 ("its leading part is $-2z\,a_2\log n/n^3$") / g3:prop:unif second claim:
`E_n = -2z a_2 log n/n³ + Δ_n`, `|Δ_n| ≤ 9000 z/n³`, for `n ≥ 9`, `0 < θ < π`. -/
def stmt_9 : Prop :=
  ∀ θ n : ℝ, 0 < θ → θ < π → 9 ≤ n →
    |Y02.E θ n - (-(2 * Y02.zz θ * Y02.a2 θ * Real.log n / n ^ 3))|
      ≤ 9000 * Y02.zz θ / n ^ 3

/-- g3:prop:n3 with the constant of g3:cor:C3 (the `n^{-3}` term, containing `a_3`, so `ζ(3)`):
`|E_n - (d(θ) - 2z a_2 log n)/n³| ≤ 3.1·10⁴ z log²n/n⁴`, `n ≥ 9`, `0 < θ < π`. -/
def stmt_10 : Prop :=
  ∀ θ n : ℝ, 0 < θ → θ < π → 9 ≤ n →
    |Y02.E θ n - (d θ - 2 * Y02.zz θ * Y02.a2 θ * Real.log n) / n ^ 3|
      ≤ 3.1e4 * Y02.zz θ * Real.log n ^ 2 / n ^ 4

/-- g3:prop:n4(b) (the `n^{-4}` term, containing `a_4`):
`|E⁽⁴⁾_n| ≤ 1.2·10⁵ z log³n/n⁵`, `n ≥ 9`, `0 < θ < π`. -/
def stmt_11 : Prop :=
  ∀ θ n : ℝ, 0 < θ → θ < π → 9 ≤ n →
    |E4 θ n| ≤ 1.2e5 * Y02.zz θ * Real.log n ^ 3 / n ^ 5

/-- g3:prop:n4(a): `E⁽⁴⁾_n = -(z a_2/3) log³n/n⁵ + Δ'_n`, `|Δ'_n| ≤ 2.5·10⁵ z log²n/n⁵`,
`n ≥ 9`, `0 < θ < π`. -/
def stmt_12 : Prop :=
  ∀ θ n : ℝ, 0 < θ → θ < π → 9 ≤ n →
    |E4 θ n - (-(Y02.zz θ * Y02.a2 θ / 3 * Real.log n ^ 3 / n ^ 5))|
      ≤ 2.5e5 * Y02.zz θ * Real.log n ^ 2 / n ^ 5

end Z06
-- END STATEMENTS

namespace Z06h
open Real MeasureTheory Filter Topology Set

lemma zz_half : Y02.zz (π / 2) = 1 / 2 := by
  unfold Y02.zz
  rw [show π / 2 / 2 = π / 4 by ring, Real.sin_pi_div_four, div_pow, Real.sq_sqrt (by norm_num)]
  norm_num

lemma ww_half : Y02.ww (π / 2) = 1 := by unfold Y02.ww; rw [zz_half]; norm_num

lemma ell_half : Y02.ell (π / 2) = -Real.log 2 := by
  unfold Y02.ell; rw [zz_half, show (1:ℝ) - 1 / 2 = 2⁻¹ by norm_num, Real.log_inv]

lemma U_half : Y02.U (π / 2) = -Real.log 2 := by
  unfold Y02.U; rw [ww_half, ell_half]; ring

/-! Trilogarithm tools (copied in spirit from Y01, for `Y02.Li3`). -/

noncomputable def H (t : ℝ) : ℝ := if t = 0 then 1 else D00.Li2 t / t

lemma H_of_ne {t : ℝ} (h : t ≠ 0) : H t = D00.Li2 t / t := by simp [H, h]

lemma continuousAt_H_zero : ContinuousAt H 0 := by
  rw [← continuousWithinAt_compl_self, ContinuousWithinAt]
  have hd := D00.hasDerivAt_Li2_zero
  rw [hasDerivAt_iff_tendsto_slope] at hd
  have h0 : H 0 = 1 := by simp [H]
  rw [h0]
  apply hd.congr'
  filter_upwards [self_mem_nhdsWithin] with t ht
  rw [H_of_ne ht, slope_def_field]
  simp

lemma continuous_H : Continuous H := by
  rw [continuous_iff_continuousAt]
  intro t
  rcases eq_or_ne t 0 with rfl | h0
  · exact continuousAt_H_zero
  · have hc : ContinuousAt (fun s => D00.Li2 s / s) t :=
      (D00.continuous_Li2.continuousAt).div continuousAt_id h0
    apply hc.congr
    filter_upwards [isOpen_ne.mem_nhds h0] with s hs
    exact (H_of_ne hs).symm

lemma Li3_eq_H (x : ℝ) : Y02.Li3 x = ∫ t in (0:ℝ)..x, H t := by
  unfold Y02.Li3
  apply intervalIntegral.integral_congr_ae
  filter_upwards [(Set.countable_singleton (0:ℝ)).ae_notMem volume] with t ht _
  exact (H_of_ne ht).symm

lemma hasDerivAt_Li3_H (x : ℝ) : HasDerivAt Y02.Li3 (H x) x := by
  have e : Y02.Li3 = fun x => ∫ t in (0:ℝ)..x, H t := funext Li3_eq_H
  rw [e]
  exact intervalIntegral.integral_hasDerivAt_right (continuous_H.intervalIntegrable _ _)
    (continuous_H.stronglyMeasurableAtFilter _ _) continuous_H.continuousAt

lemma continuous_Li3 : Continuous Y02.Li3 :=
  continuous_iff_continuousAt.2 fun x => (hasDerivAt_Li3_H x).continuousAt

lemma hasDerivAt_Li3 {x : ℝ} (h0 : x ≠ 0) : HasDerivAt Y02.Li3 (D00.Li2 x / x) x := by
  have := hasDerivAt_Li3_H x
  rwa [H_of_ne h0] at this

lemma Li3_zero : Y02.Li3 0 = 0 := by simp [Y02.Li3]

lemma Li3_comp {f : ℝ → ℝ} {f' x : ℝ} (hf : HasDerivAt f f' x) (h0 : f x ≠ 0) :
    HasDerivAt (fun y => Y02.Li3 (f y)) (D00.Li2 (f x) / f x * f') x :=
  (hasDerivAt_Li3 h0).comp x hf

lemma D3 {w : ℝ} (hw : |w| ≤ 1) : Y02.Li3 (w ^ 2) = 4 * Y02.Li3 w + 4 * Y02.Li3 (-w) := by
  have key := D00.eq_of_deriv_zero' (f := fun u => Y02.Li3 (u ^ 2) - 4 * Y02.Li3 u - 4 * Y02.Li3 (-u))
    (a := 0) (b := w)
    ((continuous_Li3.comp (continuous_pow 2)).sub (continuous_const.mul continuous_Li3)
      |>.sub (continuous_const.mul (continuous_Li3.comp continuous_neg))).continuousOn
    (by
      intro u hu
      obtain ⟨hu0, huw⟩ := D00.mem_uIoo_zero hu
      have hu1 : |u| ≤ 1 := (lt_of_lt_of_le huw hw).le
      have d1 := Li3_comp (hasDerivAt_pow 2 u) (pow_ne_zero 2 hu0)
      have d2 := hasDerivAt_Li3 hu0
      have d3 := Li3_comp (hasDerivAt_neg u) (neg_ne_zero.mpr hu0)
      have := (d1.sub (d2.const_mul 4)).sub (d3.const_mul 4)
      refine this.congr_deriv ?_
      have dup := D00.duplication hu1
      have e : D00.Li2 (u ^ 2) = 2 * D00.Li2 u + 2 * D00.Li2 (-u) := by linarith
      simp only [e]
      field_simp
      ring)
  simp [Li3_zero] at key
  linarith

lemma g_tendsto (z : ℝ) :
    Tendsto (fun u => Real.log u ^ 2 * Real.log (1 - u)) (𝓝[Icc z 1] 1) (𝓝 0) := by
  have hc : Continuous (fun u : ℝ => 4 * (1 - u) * |(1 - u) * Real.log (1 - u)|) :=
    (continuous_const.mul (continuous_const.sub continuous_id)).mul
      ((Real.continuous_mul_log.comp (continuous_const.sub continuous_id)).abs)
  have ht : Tendsto (fun u : ℝ => 4 * (1 - u) * |(1 - u) * Real.log (1 - u)|) (𝓝[Icc z 1] 1)
      (𝓝 0) := by
    have := (hc.tendsto 1).mono_left (nhdsWithin_le_nhds (s := Icc z 1))
    simpa using this
  apply squeeze_zero_norm' _ ht
  have h12 : ∀ᶠ u in 𝓝[Icc z 1] (1:ℝ), (1:ℝ) / 2 < u :=
    nhdsWithin_le_nhds (lt_mem_nhds (by norm_num))
  filter_upwards [self_mem_nhdsWithin, h12] with u hu h2
  have hu0 : 0 < u := by linarith
  have hu1 : u ≤ 1 := hu.2
  have hl1 : Real.log u ≤ 0 := Real.log_nonpos hu0.le hu1
  have hl2 : 1 - u⁻¹ ≤ Real.log u := Real.one_sub_inv_le_log_of_pos hu0
  have hinv : u⁻¹ ≤ 2 := by
    rw [inv_le_comm₀ hu0 (by norm_num)]; linarith
  have hl3 : -(2 * (1 - u)) ≤ Real.log u := by
    have : 1 - u⁻¹ = -((1 - u) * u⁻¹) := by
      rw [mul_comm, mul_sub, inv_mul_cancel₀ hu0.ne']; ring
    nlinarith
  have hsq : Real.log u ^ 2 ≤ 4 * (1 - u) ^ 2 := by nlinarith
  rw [Real.norm_eq_abs, abs_mul, abs_of_nonneg (sq_nonneg _), abs_mul,
    abs_of_nonneg (by linarith : (0:ℝ) ≤ 1 - u)]
  have := abs_nonneg (Real.log (1 - u))
  nlinarith

lemma Tr {z : ℝ} (hz0 : 0 < z) (hz1 : z < 1) :
    Y02.Li3 z + Y02.Li3 (1 - z) + Y02.Li3 (1 - 1 / z)
      = Y02.z3 + (1 / 6) * Real.log z ^ 3 + π ^ 2 / 6 * Real.log z
        - (1 / 2) * Real.log z ^ 2 * Real.log (1 - z) := by
  set F : ℝ → ℝ := fun u => Y02.Li3 u + Y02.Li3 (1 - u) + Y02.Li3 (1 - u⁻¹)
    - ((1 / 6) * Real.log u ^ 3 + π ^ 2 / 6 * Real.log u
        - (1 / 2) * (Real.log u ^ 2 * Real.log (1 - u))) with hF
  have hd : ∀ u : ℝ, 0 < u → u < 1 → HasDerivAt F 0 u := by
    intro u hu0 hu1
    have hu0' : u ≠ 0 := hu0.ne'
    have h1u : 1 - u ≠ 0 := by linarith
    have hA : 1 - u⁻¹ ≠ 0 := by
      intro h
      have : u⁻¹ = 1 := by linarith
      rw [inv_eq_one] at this; linarith
    have d1 := hasDerivAt_Li3 hu0'
    have d2 := Li3_comp ((hasDerivAt_id u).const_sub 1) h1u
    have d3 := Li3_comp ((hasDerivAt_inv hu0').const_sub 1) hA
    have dl := Real.hasDerivAt_log hu0'
    have dl1 := ((hasDerivAt_id u).const_sub 1).log h1u
    have := ((d1.add d2).add d3).sub
      (((dl.pow 3).const_mul (1 / 6) |>.add (dl.const_mul (π ^ 2 / 6))).sub
        (((dl.pow 2).mul dl1).const_mul (1 / 2)))
    refine this.congr_deriv ?_
    have hrefl := D00.reflection hu0 hu1
    have land := D00.landen (x := 1 - u) (by linarith)
    have e1 : (1 - u) / (1 - u - 1) = 1 - u⁻¹ := by
      rw [show (1:ℝ) - u - 1 = -u by ring, div_neg, sub_div, div_self hu0']; ring
    rw [e1, show 1 - (1 - u) = u by ring] at land
    have eA : D00.Li2 (1 - u⁻¹) = -(1 / 2) * Real.log u ^ 2 - D00.Li2 (1 - u) := by linarith
    have eB : D00.Li2 (1 - u) = π ^ 2 / 6 - Real.log u * Real.log (1 - u) - D00.Li2 u := by
      linarith
    simp only [id, Pi.pow_apply]
    rw [eA, eB]
    have hu1' : u - 1 ≠ 0 := by linarith
    field_simp
    ring
  have hc1 : ContinuousWithinAt F (Icc z 1) 1 := by
    have hA : ContinuousAt (fun u : ℝ => Y02.Li3 u + Y02.Li3 (1 - u) + Y02.Li3 (1 - u⁻¹)
        - ((1 / 6) * Real.log u ^ 3 + π ^ 2 / 6 * Real.log u)) 1 := by
      have hl : ContinuousAt Real.log 1 := Real.continuousAt_log one_ne_zero
      have hi : ContinuousAt (fun u : ℝ => u⁻¹) 1 := continuousAt_inv₀ one_ne_zero
      exact ((continuous_Li3.continuousAt.add
        (continuous_Li3.continuousAt.comp (continuousAt_const.sub continuousAt_id))).add
        (continuous_Li3.continuousAt.comp (continuousAt_const.sub hi))).sub
        (((hl.pow 3).const_mul _).add (hl.const_mul _))
    have hB : ContinuousWithinAt (fun u : ℝ => (1 / 2) * (Real.log u ^ 2 * Real.log (1 - u)))
        (Icc z 1) 1 := by
      rw [ContinuousWithinAt]
      have := (g_tendsto z).const_mul (1 / 2)
      simpa using this
    have := hA.continuousWithinAt.add hB
    convert this using 1
    funext u
    simp only [hF, Pi.add_apply]
    ring
  have key := D00.eq_of_deriv_zero (f := F) (a := z) (b := 1) hz1.le
    (fun u hu => by
      rcases eq_or_lt_of_le hu.2 with h | h
      · rw [h]; exact hc1
      · exact (hd u (by linarith [hu.1]) h).continuousAt.continuousWithinAt)
    (fun u hu => hd u (by linarith [hu.1]) hu.2)
  simp only [hF] at key
  norm_num [Li3_zero, Y02.Li3_one] at key
  rw [one_div]
  linarith

lemma Li3_m1 : Y02.Li3 (-1) = -(3 / 4) * Y02.z3 := by
  have h := D3 (w := 1) (by norm_num)
  rw [one_pow, Y02.Li3_one] at h; linarith

lemma Li3_half : Y02.Li3 (1 / 2) = 7 / 8 * Y02.z3 - π ^ 2 * Real.log 2 / 12
    + Real.log 2 ^ 3 / 6 := by
  have h := Tr (z := 1 / 2) (by norm_num) (by norm_num)
  have hl : Real.log (1 / 2) = -Real.log 2 := by rw [one_div, Real.log_inv]
  rw [show (1:ℝ) - 1 / 2 = 1 / 2 by norm_num, show (1:ℝ) - 1 / (1 / 2) = -1 by norm_num,
    Li3_m1, hl] at h
  linarith

end Z06h

namespace Z06h
open Real

/-- `log n ≤ 2√n - 2`. -/
lemma log_le_sqrt {n : ℝ} (hn : 0 < n) : Real.log n ≤ 2 * Real.sqrt n - 2 := by
  have hs : 0 < Real.sqrt n := Real.sqrt_pos.mpr hn
  have h := Real.log_le_sub_one_of_pos hs
  have e : Real.log n = 2 * Real.log (Real.sqrt n) := by
    rw [← Real.log_rpow hs, Real.rpow_two, Real.sq_sqrt hn.le]
  linarith

lemma log_le_half {n : ℝ} (hn : 0 < n) : Real.log n ≤ n / 2 := by
  have h := log_le_sqrt hn
  have hs := Real.sq_sqrt hn.le
  nlinarith [sq_nonneg (Real.sqrt n - 2)]

lemma logp1_sq {n : ℝ} (hn : 1 ≤ n) : (Real.log n + 1) ^ 2 ≤ 4 * n := by
  have hn0 : 0 < n := by linarith
  have h := log_le_sqrt hn0
  have hl : 0 ≤ Real.log n := Real.log_nonneg hn
  have hs := Real.sq_sqrt hn0.le
  have : Real.log n + 1 ≤ 2 * Real.sqrt n := by linarith
  nlinarith

lemma exp_neg_bounds {t : ℝ} (ht : 0 ≤ t) :
    1 - t ≤ Real.exp (-t) ∧ Real.exp (-t) ≤ 1 - t + t ^ 2 / 2 := by
  constructor
  · linarith [Real.add_one_le_exp (-t)]
  · -- `g(t) = 1 - t + t²/2 - e^{-t}` is nondecreasing from `0`
    have key : ∀ s : ℝ, 0 ≤ s → Real.exp (-s) ≤ 1 - s + s ^ 2 / 2 := by
      intro s hs
      have hmono : MonotoneOn (fun u : ℝ => 1 - u + u ^ 2 / 2 - Real.exp (-u)) (Set.Ici 0) := by
        apply monotoneOn_of_deriv_nonneg (convex_Ici 0)
        · fun_prop
        · fun_prop
        · intro u hu
          rw [interior_Ici] at hu
          have hd : HasDerivAt (fun u : ℝ => 1 - u + u ^ 2 / 2 - Real.exp (-u))
              (-1 + u + Real.exp (-u)) u := by
            have h := (((hasDerivAt_id' u).const_sub (1:ℝ)).add
              ((hasDerivAt_pow 2 u).div_const 2)).sub ((hasDerivAt_neg u).exp)
            convert h using 1
            push_cast; ring
          rw [hd.deriv]
          linarith [Real.add_one_le_exp (-u)]
      have := hmono (Set.mem_Ici.mpr le_rfl) (Set.mem_Ici.mpr hs) hs
      simp at this
      linarith
    exact key t ht

lemma poly1 {n : ℝ} (hn : 9 ≤ n) :
    n / 2 * (1 / (n - 1) ^ 3 - 1 / n ^ 3) + 1 / (n - 1) ^ 3 ≤ 6.4 / n ^ 3 := by
  have hn0 : 0 < n := by linarith
  have hn1 : 0 < n - 1 := by linarith
  have e : n / 2 * (1 / (n - 1) ^ 3 - 1 / n ^ 3) + 1 / (n - 1) ^ 3
      = (5 * n ^ 3 - 3 * n ^ 2 + n) / (2 * ((n - 1) ^ 3 * n ^ 3)) := by
    field_simp; ring
  rw [e, div_le_div_iff₀ (by positivity) (by positivity)]
  have h4 : (8 * n) ^ 3 ≤ (9 * (n - 1)) ^ 3 := pow_le_pow_left₀ (by linarith) (by linarith) 3
  have e1 : (8 * n) ^ 3 = 512 * n ^ 3 := by ring
  have e2 : (9 * (n - 1)) ^ 3 = 729 * (n - 1) ^ 3 := by ring
  have h6 : n ≤ n ^ 2 := by nlinarith
  have h7 : 0 ≤ n ^ 3 := by positivity
  rw [e1, e2] at h4
  have key : 5 * n ^ 3 - 3 * n ^ 2 + n ≤ 6.4 * (2 * (n - 1) ^ 3) := by linarith
  have h5 := mul_le_mul_of_nonneg_right key h7
  have e4 : 6.4 * (2 * ((n - 1) ^ 3 * n ^ 3)) = 6.4 * (2 * (n - 1) ^ 3) * n ^ 3 := by ring
  rw [e4]; exact h5

lemma poly2 {n : ℝ} (hn : 9 ≤ n) : 4 * n * n ^ 3 ≤ 3.4 * (2 * (n - 1) ^ 4) := by
  have h4 : (8 * n) ^ 4 ≤ (9 * (n - 1)) ^ 4 := pow_le_pow_left₀ (by linarith) (by linarith) 4
  have e1 : (8 * n) ^ 4 = 4096 * n ^ 4 := by ring
  have e2 : (9 * (n - 1)) ^ 4 = 6561 * (n - 1) ^ 4 := by ring
  have e3 : 4 * n * n ^ 3 = 4 * n ^ 4 := by ring
  have : 0 ≤ (n - 1) ^ 4 := by positivity
  rw [e1, e2] at h4; rw [e3]; linarith

end Z06h

namespace Z06h
open Real MeasureTheory Set

lemma intg_int {t : ℝ} (ht : t < 1) :
    IntervalIntegrable (fun s : ℝ => Real.log (1 - s) ^ 2 / s) volume 0 t := by
  have e : (fun s : ℝ => Real.log (1 - s) ^ 2 / s) = fun s => D00.G s * Real.log (1 - s) :=
    funext Y02.integrand_eq
  rw [e]
  apply ContinuousOn.intervalIntegrable
  intro s hs
  apply (Y02.continuousAt_g _).continuousWithinAt
  rcases le_total 0 t with h | h
  · rw [uIcc_of_le h] at hs; linarith [hs.2]
  · rw [uIcc_of_ge h] at hs; linarith [hs.2]

lemma S12int_nonneg {t : ℝ} (ht : 0 ≤ t) : 0 ≤ Z06.S12int t := by
  unfold Z06.S12int
  apply mul_nonneg (by norm_num)
  apply intervalIntegral.integral_nonneg ht
  intro s hs
  exact div_nonneg (sq_nonneg _) hs.1

lemma S12int_small {t : ℝ} (h0 : 0 ≤ t) (h1 : t ≤ 1 / 2) : Z06.S12int t ≤ t ^ 2 := by
  unfold Z06.S12int
  have hle : ∫ s in (0:ℝ)..t, Real.log (1 - s) ^ 2 / s ≤ ∫ s in (0:ℝ)..t, 4 * s := by
    apply intervalIntegral.integral_mono_on h0 (intg_int (by linarith))
      (by apply Continuous.intervalIntegrable; fun_prop)
    intro s hs
    rcases eq_or_lt_of_le hs.1 with h | h
    · rw [← h]; simp
    · have hs1 : 0 < 1 - s := by linarith [hs.2]
      have hlo := Real.one_sub_inv_le_log_of_pos hs1
      have hhi : Real.log (1 - s) ≤ 0 := Real.log_nonpos hs1.le (by linarith)
      have hinv : (1 - s)⁻¹ ≤ 2 := by
        rw [inv_le_comm₀ hs1 (by norm_num)]; linarith [hs.2]
      have e : 1 - (1 - s)⁻¹ = -(s * (1 - s)⁻¹) := by field_simp; ring
      have hlo' : -(2 * s) ≤ Real.log (1 - s) := by nlinarith
      rw [div_le_iff₀ h]
      nlinarith
  have hcomp : ∫ s in (0:ℝ)..t, 4 * s = 2 * t ^ 2 := by
    rw [intervalIntegral.integral_const_mul, integral_id]; ring
  rw [hcomp] at hle
  linarith

lemma Li3_nonneg {x : ℝ} (h0 : 0 ≤ x) (h1 : x < 1) : 0 ≤ Y02.Li3 x := by
  unfold Y02.Li3
  apply intervalIntegral.integral_nonneg h0
  intro t ht
  rcases eq_or_lt_of_le ht.1 with h | h
  · rw [← h]; simp
  · exact div_nonneg (Y02.Li2_pos_lt h (by linarith [ht.2])).1.le h.le

lemma z3_le : Y02.z3 ≤ 2 := by
  have h2 := D00.hasSum_Li2 (x := 1) (by norm_num)
  have hle := hasSum_le (fun n : ℕ => by
    show 1 / ((n : ℝ) + 1) ^ 3 ≤ (1:ℝ) ^ (n + 1) / ((n : ℝ) + 1) ^ 2
    rw [one_pow]
    apply one_div_le_one_div_of_le (by positivity)
    have : (1:ℝ) ≤ (n:ℝ) + 1 := by linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
    nlinarith) Y02.hasSum_z3 h2
  rw [D00.Li2_one] at hle
  have := Real.pi_lt_d2
  nlinarith [Real.pi_pos]

lemma S12int_large {z : ℝ} (h0 : 1 / 2 ≤ z) (h1 : z < 1) : Z06.S12int z ≤ Y02.z3 := by
  have hS := Y02.S12_identity (z := z) (by linarith) h1
  have e : Z06.S12int z = (1 / 2) * ∫ t in (0 : ℝ)..z, Real.log (1 - t) ^ 2 / t := rfl
  rw [e, ← hS]
  have hL3 := Li3_nonneg (x := 1 - z) (by linarith) (by linarith)
  have hl : Real.log (1 - z) ≤ 0 := Real.log_nonpos (by linarith) (by linarith)
  have hL2 := (Y02.Li2_pos_lt (z := 1 - z) (by linarith) (by linarith)).1
  have hlz : Real.log z ≤ 0 := Real.log_nonpos (by linarith) h1.le
  have p1 : Real.log (1 - z) * D00.Li2 (1 - z) ≤ 0 := mul_nonpos_of_nonpos_of_nonneg hl hL2.le
  have p2 : Real.log (1 - z) ^ 2 * Real.log z ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos (sq_nonneg _) hlz
  linarith

lemma S12int_le {t : ℝ} (h0 : 0 ≤ t) (h1 : t < 1) : Z06.S12int t ≤ 4 * t := by
  rcases le_total t (1 / 2) with h | h
  · have := S12int_small h0 h; nlinarith
  · have := S12int_large h h1; linarith [z3_le]

lemma S22_abs {z : ℝ} (h0 : 0 ≤ z) (h1 : z < 1) : |Z06.S22 z| ≤ 4 * z := by
  unfold Z06.S22
  have := intervalIntegral.norm_integral_le_of_norm_le_const (a := 0) (b := z) (C := 4)
    (f := fun t => Z06.S12int t / t) (by
      intro t ht
      rw [uIoc_of_le h0] at ht
      have ht0 : 0 < t := ht.1
      rw [Real.norm_eq_abs, abs_div, abs_of_nonneg (S12int_nonneg ht0.le), abs_of_pos ht0,
        div_le_iff₀ ht0]
      exact S12int_le ht0.le (by linarith [ht.2]))
  rw [Real.norm_eq_abs, sub_zero, abs_of_nonneg h0] at this
  linarith

lemma S12_eq (θ : ℝ) (h0 : 0 < θ) (h1 : θ < π) : Y02.S12 θ = Z06.S12int (Y02.zz θ) :=
  Y02.proof_10 θ h0 h1

lemma wS12_bounds (θ : ℝ) (h0 : 0 < θ) (h1 : θ < π) :
    0 ≤ Y02.ww θ * Y02.S12 θ ∧ Y02.ww θ * Y02.S12 θ ≤ 4 := by
  have z0 : 0 < Y02.zz θ := (Y02.zz_bounds θ h0 h1).1
  have z1 : Y02.zz θ < 1 := (Y02.zz_bounds θ h0 h1).2
  rw [S12_eq θ h0 h1]
  have hS0 := S12int_nonneg z0.le
  have hS1 := S12int_le z0.le z1
  unfold Y02.ww
  set z := Y02.zz θ
  have hS1' : (1 - z) * Z06.S12int z ≤ 4 * z := by nlinarith
  constructor
  · exact mul_nonneg (div_nonneg (by linarith) z0.le) hS0
  · rw [div_mul_eq_mul_div, div_le_iff₀ z0]
    linarith

lemma a3_abs (θ : ℝ) (h0 : 0 < θ) (h1 : θ < π) : |Y02.a3 θ| ≤ 9 := by
  obtain ⟨hw0, hw1⟩ := wS12_bounds θ h0 h1
  obtain ⟨hU0, hU1, -⟩ := Y02.proof_9 θ h0 h1
  unfold Y02.a3
  have hU3 : -1 ≤ Y02.U θ ^ 3 ∧ Y02.U θ ^ 3 ≤ 0 := by
    constructor <;> nlinarith [sq_nonneg (Y02.U θ)]
  rw [abs_le]; constructor <;> nlinarith

lemma b4_abs (θ : ℝ) (h0 : 0 < θ) (h1 : θ < π) : |Z06.a4 θ + Y02.a2 θ ^ 2 / 2| ≤ 20 := by
  have z0 : 0 < Y02.zz θ := (Y02.zz_bounds θ h0 h1).1
  have z1 : Y02.zz θ < 1 := (Y02.zz_bounds θ h0 h1).2
  obtain ⟨hw0, hw1⟩ := wS12_bounds θ h0 h1
  obtain ⟨hU0, hU1, -, -, ha0, ha1⟩ := Y02.proof_9 θ h0 h1
  have hS22 := S22_abs z0.le z1
  have e : 2 * Y02.ww θ ^ 2 * Y02.ell θ * Y02.S12 θ = 2 * Y02.U θ * (Y02.ww θ * Y02.S12 θ) := by
    unfold Y02.U; ring
  unfold Z06.a4
  rw [e]
  have hpi := Real.pi_lt_d2
  have hpi0 := Real.pi_pos
  have hp4 : π ^ 4 / 180 ≤ 1 := by
    have : π ^ 2 ≤ 10 := by nlinarith
    nlinarith
  have hp4' : 0 ≤ π ^ 4 / 180 := by positivity
  have ha2 : Y02.a2 θ ^ 2 / 2 ≤ 5 := by nlinarith
  have ha2' : 0 ≤ Y02.a2 θ ^ 2 / 2 := by positivity
  have hU4 : 0 ≤ Y02.U θ ^ 4 ∧ Y02.U θ ^ 4 ≤ 1 := by
    constructor
    · positivity
    · have : Y02.U θ ^ 2 ≤ 1 := by nlinarith
      nlinarith [sq_nonneg (Y02.U θ)]
  have hUw : |2 * Y02.U θ * (Y02.ww θ * Y02.S12 θ)| ≤ 8 := by
    rw [abs_mul, abs_mul, abs_of_nonneg hw0, abs_of_pos (by norm_num : (0:ℝ) < 2)]
    have : |Y02.U θ| ≤ 1 := by rw [abs_le]; constructor <;> linarith
    have := mul_le_mul this hw1 hw0 (by norm_num)
    linarith
  rw [abs_le] at hS22 hUw ⊢
  constructor <;> nlinarith

end Z06h


open Real in
theorem Z06h.p1 : Y02.a1 (π / 2) = Real.log 2 - 1 ∧ Y02.mu (π / 2) = 2 / Real.exp 1 ∧
    Y02.a2 (π / 2) = π ^ 2 / 12 + 1 / 2 ∧ Y02.c2 (π / 2) = π ^ 2 / 12 := by
  have ha1 : Y02.a1 (π / 2) = Real.log 2 - 1 := by
    unfold Y02.a1; rw [Z06h.U_half]; ring
  have ha2 : Y02.a2 (π / 2) = π ^ 2 / 12 + 1 / 2 := by
    unfold Y02.a2; rw [Z06h.U_half, Z06h.zz_half, D00.Li2_half]; ring
  refine ⟨ha1, ?_, ha2, ?_⟩
  · unfold Y02.mu; rw [ha1, Real.exp_sub, Real.exp_log (by norm_num)]
  · unfold Y02.c2; rw [ha2, Z06h.zz_half]; ring

open Real in
theorem Z06h.p2 : Y02.S12 (π / 2) = Y02.z3 / 8 - Real.log 2 ^ 3 / 6 := by
  unfold Y02.S12
  rw [Z06h.zz_half, Z06h.ell_half, show (1:ℝ) - 1 / 2 = 1 / 2 by norm_num, Z06h.Li3_half,
    D00.Li2_half, one_div, Real.log_inv]
  ring

open Real in
theorem Z06h.p3 : Y02.a3 (π / 2) = Y02.z3 / 4 - 1 / 3 := by
  unfold Y02.a3; rw [Z06h.p2, Z06h.ww_half, Z06h.U_half]; ring

open Real in
theorem Z06h.p4 : ∀ θ : ℝ, 0 < θ → θ < π →
    -2 * Y02.zz θ * Y02.a1 θ
      = 2 * Real.sin (θ / 2) ^ 2 + 4 * Real.cos (θ / 2) ^ 2 * Real.log (Real.cos (θ / 2)) := by
  intro θ h0 h1
  obtain ⟨z0, z1⟩ := Y02.zz_bounds θ h0 h1
  have hc : 0 < Real.cos (θ / 2) := Real.cos_pos_of_mem_Ioo ⟨by linarith, by linarith⟩
  have h1z : 1 - Y02.zz θ = Real.cos (θ / 2) ^ 2 := by
    unfold Y02.zz; linarith [Real.sin_sq_add_cos_sq (θ / 2)]
  unfold Y02.a1 Y02.U Y02.ww Y02.ell
  rw [h1z, Real.log_pow]
  have : Y02.zz θ ≠ 0 := z0.ne'
  field_simp
  unfold Y02.zz
  push_cast
  ring


theorem Z06h.p5 : Z06.stmt_5 := by
  intro θ h0 h1
  exact ⟨(Z06h.a3_abs θ h0 h1).trans (by norm_num), (Z06h.b4_abs θ h0 h1).trans (by norm_num)⟩


open Real in
theorem Z06h.p7 : ∀ a n : ℝ, -1 < a → a < 0 → 9 ≤ n →
    |(Real.exp a / n) ^ (1 / (n - 1)) * (1 / (n - 1)) ^ 2 - 1 / n ^ 2 + Real.log n / n ^ 3|
      ≤ 6.4 / n ^ 3 := by
  intro a n ha0 ha1 hn
  have hn0 : 0 < n := by linarith
  have hn1 : 0 < n - 1 := by linarith
  set L := Real.log n with hL
  set ε := 1 / (n - 1) with hε
  have hLh : L ≤ n / 2 := Z06h.log_le_half hn0
  have hL0 : 0 ≤ L := Real.log_nonneg (by linarith)
  have hL1 := Z06h.logp1_sq (n := n) (by linarith)
  set t := ε * (L - a) with ht
  have hε0 : 0 < ε := by positivity
  have ht0 : 0 ≤ t := mul_nonneg hε0.le (by linarith)
  have hq : (Real.exp a / n) ^ (1 / (n - 1)) = Real.exp (-t) := by
    rw [Real.rpow_def_of_pos (div_pos (Real.exp_pos _) hn0),
      Real.log_div (Real.exp_pos _).ne' hn0.ne', Real.log_exp, ← hε, ht, ← hL]
    congr 1; ring
  rw [hq]
  obtain ⟨hlo, hhi⟩ := Z06h.exp_neg_bounds ht0
  have hε2 : 0 ≤ ε ^ 2 := sq_nonneg _
  have lo := mul_le_mul_of_nonneg_right hlo hε2
  have hi := mul_le_mul_of_nonneg_right hhi hε2
  -- algebraic bounds, in terms of `n`
  have hεn : ε * (n - 1) = 1 := by rw [hε]; field_simp
  have B1 : ε ^ 2 - 1 / n ^ 2 ≤ 3 / n ^ 3 := by
    rw [hε, div_pow, one_pow, div_sub_div _ _ (by positivity) (by positivity),
      div_le_div_iff₀ (by positivity) (by positivity)]
    have : 0 ≤ n ^ 2 * (n ^ 2 - 5 * n + 3) := mul_nonneg (by positivity) (by nlinarith)
    nlinarith
  have B0 : 0 ≤ ε ^ 2 - 1 / n ^ 2 := by
    rw [hε, div_pow, one_pow, sub_nonneg]
    apply one_div_le_one_div_of_le (by positivity)
    nlinarith
  have hLa : L - a ≤ L + 1 := by linarith
  have hLa0 : 0 ≤ L - a := by linarith
  -- `t² ε²/2 ≤ 3.4/n³`
  have B2 : t ^ 2 * ε ^ 2 / 2 ≤ 3.4 / n ^ 3 := by
    have h1 : (L - a) ^ 2 ≤ 4 * n := le_trans (pow_le_pow_left₀ hLa0 hLa 2) hL1
    have e : t ^ 2 * ε ^ 2 / 2 = (L - a) ^ 2 / (2 * (n - 1) ^ 4) := by
      rw [ht, hε]; field_simp
    rw [e, div_le_div_iff₀ (by positivity) (by positivity)]
    have h2 : (L - a) ^ 2 * n ^ 3 ≤ 4 * n * n ^ 3 := mul_le_mul_of_nonneg_right h1 (by positivity)
    have h3 : 4 * n * n ^ 3 ≤ 3.4 * (2 * (n - 1) ^ 4) := Z06h.poly2 hn
    linarith
  -- `B = t ε² - L/n³ ∈ [0, 6.4/n³]`
  have Be : t * ε ^ 2 - L / n ^ 3 = L * (ε ^ 3 - 1 / n ^ 3) - a * ε ^ 3 := by rw [ht]; ring
  have hε3 : ε ^ 3 = 1 / (n - 1) ^ 3 := by rw [hε, div_pow, one_pow]
  have Bnn : 0 ≤ t * ε ^ 2 - L / n ^ 3 := by
    rw [Be]
    have : 1 / n ^ 3 ≤ ε ^ 3 := by
      rw [hε3]; apply one_div_le_one_div_of_le (by positivity)
      exact pow_le_pow_left₀ hn1.le (by linarith) 3
    have : 0 ≤ -a * ε ^ 3 := mul_nonneg (by linarith) (by positivity)
    nlinarith
  have Bup : t * ε ^ 2 - L / n ^ 3 ≤ 6.4 / n ^ 3 := by
    rw [Be]
    have h1 : -a * ε ^ 3 ≤ ε ^ 3 := by
      have : 0 ≤ ε ^ 3 := by positivity
      nlinarith
    have h2 : L * (ε ^ 3 - 1 / n ^ 3) ≤ n / 2 * (ε ^ 3 - 1 / n ^ 3) := by
      apply mul_le_mul_of_nonneg_right hLh
      rw [hε3, sub_nonneg]; apply one_div_le_one_div_of_le (by positivity)
      exact pow_le_pow_left₀ hn1.le (by linarith) 3
    have h3 : n / 2 * (ε ^ 3 - 1 / n ^ 3) + ε ^ 3 ≤ 6.4 / n ^ 3 := by
      rw [hε3]; exact Z06h.poly1 hn
    linarith
  have lo' : ε ^ 2 - t * ε ^ 2 ≤ Real.exp (-t) * ε ^ 2 := by
    linarith [lo, show (1 - t) * ε ^ 2 = ε ^ 2 - t * ε ^ 2 by ring]
  have hi' : Real.exp (-t) * ε ^ 2 ≤ ε ^ 2 - t * ε ^ 2 + t ^ 2 * ε ^ 2 / 2 := by
    linarith [hi, show (1 - t + t ^ 2 / 2) * ε ^ 2 = ε ^ 2 - t * ε ^ 2 + t ^ 2 * ε ^ 2 / 2 by ring]
  have h64 : 3 / n ^ 3 + 3.4 / n ^ 3 = 6.4 / n ^ 3 := by ring
  rw [abs_le]
  constructor
  · linarith
  · linarith

namespace Z06h
open Real

/-- The splitting `E_n = T_1 + T_2 + T_3 + T_4` of the proof of g3:prop:unif. -/
lemma decomp (θ n : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hn : 9 ≤ n) :
    ∃ T4 : ℝ, HasSum (fun i => Y02.serTerm θ n (i + 1)) T4 ∧
      Y02.E θ n = 2 * Y02.zz θ * Y02.a2 θ *
          ((Y02.mu θ / n) ^ (1 / (n - 1)) * (1 / (n - 1)) ^ 2 - 1 / n ^ 2)
        + 2 * Y02.zz θ * (Y02.mu θ / n) ^ (1 / (n - 1)) *
          (Y02.R θ (1 / (n - 1)) * Real.exp (-(Y02.a1 θ * (1 / (n - 1)))) / (2 * Y02.zz θ)
            - (1 + Y02.a2 θ * (1 / (n - 1)) ^ 2))
        + (Y02.serTerm θ n 0 + 2 * Y02.zz θ * (1 - Y02.zz θ) * (1 / n) ^ 2) + T4 := by
  obtain ⟨z0, z1⟩ := Y02.zz_bounds θ h0 h1
  have hzz : Y02.zz θ = Real.sin (θ / 2) ^ 2 := rfl
  rw [← hzz] at z0 z1
  obtain ⟨hε0, hε1, -, -, -⟩ := Y02.eps_facts hn
  have hlen := Y02.length_at θ n h0.le h1 (by linarith) (Y02.nonint_small hε0 (by linarith))
    (Y02.Areg_small θ (Y02.eps n) h0.le h1 hε0 (by linarith))
  have hT := (hasSum_nat_add_iff' 1).mpr hlen
  simp only [Finset.sum_range_one] at hT
  refine ⟨_, hT, ?_⟩
  have hn0 : 0 < n := by linarith
  have hε : Y02.eps n = 1 / (n - 1) := rfl
  rw [← hε]
  set z := Y02.zz θ with hzdef
  set ε := Y02.eps n with hεdef
  set X := Real.exp (Y02.a1 θ * ε) with hX
  set Y := Real.exp (-(ε * Real.log n)) with hY
  have hYn : n ^ (-ε) = Y := by
    rw [Real.rpow_def_of_pos hn0, hY]; congr 1; ring
  have hmu : (Y02.mu θ / n) ^ ε = X * Y := by
    rw [Y02.mu, Real.rpow_def_of_pos (div_pos (Real.exp_pos _) hn0),
      Real.log_div (Real.exp_pos _).ne' hn0.ne', Real.log_exp, hX, hY, ← Real.exp_add]
    congr 1; ring
  have hz0 : z ≠ 0 := z0.ne'
  have hXX : Real.exp (-(Y02.a1 θ * ε)) * X = 1 := by rw [hX, ← Real.exp_add]; simp
  have hR : π * ε / Real.sin (π * ε) * Y02.Lam ε θ = Y02.R θ ε * Real.exp (-(Y02.a1 θ * ε)) * X := by
    rw [mul_assoc, hXX, mul_one]; rfl
  have hL : Y02.L θ n = 1 + Real.cos θ + π * ε / Real.sin (π * ε) * Y02.Lam ε θ * n ^ (-ε)
      + Y02.serTerm θ n 0 + (∑' i, Y02.serTerm θ n (i + 1)) := by
    rw [hT.tsum_eq]; ring
  rw [Y02.E, show (1:ℝ) / (n - 1) = ε from hε.symm, hmu, hL, hR, hYn, Y02.c2, ← hzdef]
  field_simp
  ring

lemma q_facts (θ n : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hn : 9 ≤ n) :
    0 < (Y02.mu θ / n) ^ (1 / (n - 1)) ∧ (Y02.mu θ / n) ^ (1 / (n - 1)) ≤ 1 ∧
    (Y02.mu θ / n) ^ (1 / (n - 1)) = Real.exp (-((1 / (n - 1)) * (Real.log n - Y02.a1 θ))) := by
  have hn0 : 0 < n := by linarith
  have ha := (Y02.proof_9 θ h0 h1)
  have he : (Y02.mu θ / n) ^ (1 / (n - 1)) = Real.exp (-((1 / (n - 1)) * (Real.log n - Y02.a1 θ))) := by
    rw [Y02.mu, Real.rpow_def_of_pos (div_pos (Real.exp_pos _) hn0),
      Real.log_div (Real.exp_pos _).ne' hn0.ne', Real.log_exp]
    congr 1; ring
  refine ⟨Real.rpow_pos_of_pos (div_pos (Real.exp_pos _) hn0) _, ?_, he⟩
  rw [he, Real.exp_le_one_iff, neg_nonpos]
  have hn1 : 0 < n - 1 := by linarith
  apply mul_nonneg (div_pos one_pos hn1).le
  have := Real.log_nonneg (by linarith : (1:ℝ) ≤ n)
  linarith [ha.2.2.2.1]

lemma eps_cube {n : ℝ} (hn : 9 ≤ n) : (1 / (n - 1)) ^ 3 ≤ 729 / 512 * (1 / n) ^ 3 := by
  have hn0 : 0 < n := by linarith
  have hn1 : 0 < n - 1 := by linarith
  have h : 1 / (n - 1) ≤ 9 / 8 * (1 / n) := by
    rw [mul_one_div, div_le_div_iff₀ hn1 (by positivity)]; linarith
  have := pow_le_pow_left₀ (by positivity) h 3
  calc (1 / (n - 1)) ^ 3 ≤ (9 / 8 * (1 / n)) ^ 3 := this
    _ = 729 / 512 * (1 / n) ^ 3 := by ring

lemma p8 : ∀ θ n : ℝ, 0 < θ → θ < π → 9 ≤ n →
    |Y02.E θ n - 2 * Y02.zz θ * Y02.a2 θ *
        ((Y02.mu θ / n) ^ (1 / (n - 1)) * (1 / (n - 1)) ^ 2 - 1 / n ^ 2)|
      ≤ 8880 * Y02.zz θ / n ^ 3 := by
  intro θ n h0 h1 hn
  have z0 : 0 < Y02.zz θ := (Y02.zz_bounds θ h0 h1).1
  obtain ⟨T4, hT, hE⟩ := decomp θ n h0 h1 hn
  obtain ⟨hq0, hq1, -⟩ := q_facts θ n h0 h1 hn
  obtain ⟨hε0, hε1, -, -, -⟩ := Y02.eps_facts hn
  have hρ := Y02.proof_11 θ (Y02.eps n) h0 h1 hε0 hε1
  have t3 := Y02.T3_bound θ n h0 h1 hn
  have t4 := Y02.T4_bound θ n h0 h1 hn hT
  have hε : Y02.eps n = 1 / (n - 1) := rfl
  rw [hε] at hρ
  set z := Y02.zz θ
  set q := (Y02.mu θ / n) ^ (1 / (n - 1))
  set ε := 1 / (n - 1)
  set ρ := Y02.R θ ε * Real.exp (-(Y02.a1 θ * ε)) / (2 * z) - (1 + Y02.a2 θ * ε ^ 2) with hρdef
  have hρ' : |ρ| ≤ 3100 * ε ^ 3 := by
    have e : Y02.R θ ε * Real.exp (-(Y02.a1 θ * ε)) - 2 * z * (1 + Y02.a2 θ * ε ^ 2) = 2 * z * ρ := by
      rw [hρdef]; field_simp
    rw [e, abs_mul, abs_of_pos (by positivity)] at hρ
    exact le_of_mul_le_mul_left hρ (by positivity)
  have hc := eps_cube hn
  have t2 : |2 * z * q * ρ| ≤ 2 * z * (3100 * (729 / 512 * (1 / n) ^ 3)) := by
    rw [abs_mul, abs_of_pos (by positivity)]
    have h1 : 2 * z * q ≤ 2 * z := by nlinarith
    have h2 : |ρ| ≤ 3100 * (729 / 512 * (1 / n) ^ 3) := hρ'.trans (by linarith)
    exact mul_le_mul h1 h2 (abs_nonneg _) (by positivity)
  have e : Y02.E θ n - 2 * z * Y02.a2 θ * (q * ε ^ 2 - 1 / n ^ 2)
      = 2 * z * q * ρ + (Y02.serTerm θ n 0 + 2 * z * (1 - z) * (1 / n) ^ 2) + T4 := by
    rw [hE]; ring
  rw [e]
  have hfin : 8880 * z / n ^ 3 = 8880 * (z * (1 / n) ^ 3) := by field_simp
  rw [hfin]
  have u1 := abs_add_le (2 * z * q * ρ + (Y02.serTerm θ n 0 + 2 * z * (1 - z) * (1 / n) ^ 2)) T4
  have u2 := abs_add_le (2 * z * q * ρ) (Y02.serTerm θ n 0 + 2 * z * (1 - z) * (1 / n) ^ 2)
  have hp : 0 ≤ z * (1 / n) ^ 3 := by positivity
  nlinarith

lemma p9 (h7 : ∀ a n : ℝ, -1 < a → a < 0 → 9 ≤ n →
    |(Real.exp a / n) ^ (1 / (n - 1)) * (1 / (n - 1)) ^ 2 - 1 / n ^ 2 + Real.log n / n ^ 3|
      ≤ 6.4 / n ^ 3) : ∀ θ n : ℝ, 0 < θ → θ < π → 9 ≤ n →
    |Y02.E θ n - (-(2 * Y02.zz θ * Y02.a2 θ * Real.log n / n ^ 3))|
      ≤ 9000 * Y02.zz θ / n ^ 3 := by
  intro θ n h0 h1 hn
  have z0 : 0 < Y02.zz θ := (Y02.zz_bounds θ h0 h1).1
  obtain ⟨-, -, ha1, ha2, hb0, hb1⟩ := Y02.proof_9 θ h0 h1
  have A := p8 θ n h0 h1 hn
  have B := h7 (Y02.a1 θ) n ha1 ha2 hn
  have hmu : Y02.mu θ = Real.exp (Y02.a1 θ) := rfl
  rw [← hmu] at B
  set z := Y02.zz θ
  set a := Y02.a2 θ
  set P := (Y02.mu θ / n) ^ (1 / (n - 1)) * (1 / (n - 1)) ^ 2 - 1 / n ^ 2
  have hpi := Real.pi_lt_d2
  have ha215 : a ≤ 2.16 := by nlinarith [Real.pi_pos]
  have e : Y02.E θ n - (-(2 * z * a * Real.log n / n ^ 3))
      = (Y02.E θ n - 2 * z * a * P) + 2 * z * a * (P + Real.log n / n ^ 3) := by ring
  rw [e]
  have u := abs_add_le (Y02.E θ n - 2 * z * a * P) (2 * z * a * (P + Real.log n / n ^ 3))
  have hC : |2 * z * a * (P + Real.log n / n ^ 3)| ≤ 2 * z * 2.16 * (6.4 / n ^ 3) := by
    rw [abs_mul, abs_of_pos (by positivity)]
    apply mul_le_mul (by nlinarith) B (abs_nonneg _) (by positivity)
  have hn3 : 0 < n ^ 3 := by have : 0 < n := by linarith
                             positivity
  have e2 : 2 * z * 2.16 * (6.4 / n ^ 3) = 27.648 * z / n ^ 3 := by ring
  have e3 : 9000 * z / n ^ 3 = 8880 * z / n ^ 3 + 120 * z / n ^ 3 := by ring
  have e4 : 27.648 * z / n ^ 3 ≤ 120 * z / n ^ 3 := by
    apply div_le_div_of_nonneg_right _ hn3.le; nlinarith
  linarith

end Z06h

namespace Z06h
open Real Polynomial

lemma Lam2 (θ : ℝ) : T13.LamInt 2 θ = 2 * Y02.zz θ * (1 - Y02.zz θ) * (1 - 2 * Y02.zz θ) := by
  have h := T13.proof_2 2 (by norm_num) θ
  have hb0 := congrArg (Polynomial.eval (Real.cos θ)) (T13.P_bonnet 0)
  have hb1 := congrArg (Polynomial.eval (Real.cos θ)) (T13.P_bonnet 1)
  simp only [Polynomial.eval_smul, Polynomial.eval_sub, Polynomial.eval_mul, Polynomial.eval_X,
    smul_eq_mul, T13.P_zero, T13.P_one, Polynomial.eval_one] at hb0 hb1
  norm_num at hb0 hb1 h
  have hz : Y02.zz θ = (1 - Real.cos θ) / 2 := T13.sin_half_sq θ
  rw [h, hz]
  rw [T13.P_one, Polynomial.eval_X]
  have e3 : (T13.P 3).eval (Real.cos θ) = (5 * Real.cos θ * (T13.P 2).eval (Real.cos θ)
      - 2 * Real.cos θ) / 3 := by linarith
  have e2 : (T13.P 2).eval (Real.cos θ) = (3 * Real.cos θ * Real.cos θ - 1) / 2 := by
    rw [T13.P_one, Polynomial.eval_X] at hb0; linarith
  rw [e3, e2]
  ring

/-- polynomial facts for `n ≥ 9`, `x = 1/n`, `ε = 1/(n-1)`. -/
lemma F1 {n : ℝ} (hn : 9 ≤ n) :
    0 ≤ (1 / (n - 1)) ^ 2 - (1 / n) ^ 2 - 2 * (1 / n) ^ 3 ∧
    (1 / (n - 1)) ^ 2 - (1 / n) ^ 2 - 2 * (1 / n) ^ 3 ≤ 3.6 * (1 / n) ^ 4 := by
  have hn0 : 0 < n := by linarith
  have hn1 : 0 < n - 1 := by linarith
  have e : (1 / (n - 1)) ^ 2 - (1 / n) ^ 2 - 2 * (1 / n) ^ 3 = (3 * n - 2) / (n ^ 3 * (n - 1) ^ 2) := by
    field_simp; ring
  rw [e]
  constructor
  · apply div_nonneg (by linarith) (by positivity)
  · rw [div_le_iff₀ (by positivity)]
    have e2 : 3.6 * (1 / n) ^ 4 * (n ^ 3 * (n - 1) ^ 2) = 3.6 * (n - 1) ^ 2 / n := by
      field_simp
    rw [e2, le_div_iff₀ hn0]
    nlinarith

lemma F2 {n : ℝ} (hn : 9 ≤ n) :
    0 ≤ (1 / (n - 1)) ^ 3 - (1 / n) ^ 3 ∧ (1 / (n - 1)) ^ 3 - (1 / n) ^ 3 ≤ 4 * (1 / n) ^ 4 := by
  have hn0 : 0 < n := by linarith
  have hn1 : 0 < n - 1 := by linarith
  have e : (1 / (n - 1)) ^ 3 - (1 / n) ^ 3 = (3 * n ^ 2 - 3 * n + 1) / (n ^ 3 * (n - 1) ^ 3) := by
    field_simp; ring
  rw [e]
  constructor
  · apply div_nonneg (by nlinarith) (by positivity)
  · rw [div_le_iff₀ (by positivity)]
    have e2 : 4 * (1 / n) ^ 4 * (n ^ 3 * (n - 1) ^ 3) = 4 * (n - 1) ^ 3 / n := by
      field_simp
    rw [e2, le_div_iff₀ hn0]
    have : 0 ≤ n ^ 2 * (n - 9) := mul_nonneg (by positivity) (by linarith)
    nlinarith

lemma F3 {n : ℝ} (hn : 9 ≤ n) : (1 / (n - 1)) ^ 4 ≤ 1.61 * (1 / n) ^ 4 := by
  have hn0 : 0 < n := by linarith
  have hn1 : 0 < n - 1 := by linarith
  have h : 1 / (n - 1) ≤ 9 / 8 * (1 / n) := by
    rw [mul_one_div, div_le_div_iff₀ hn1 (by positivity)]; linarith
  have := pow_le_pow_left₀ (by positivity) h 4
  have hp : 0 ≤ (1 / n) ^ 4 := by positivity
  calc (1 / (n - 1)) ^ 4 ≤ (9 / 8 * (1 / n)) ^ 4 := this
    _ = 6561 / 4096 * (1 / n) ^ 4 := by ring
    _ ≤ 1.61 * (1 / n) ^ 4 := by nlinarith

/-- `T_3 + 4z(1-z)x³`. -/
lemma P5_bound (θ n : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hn : 9 ≤ n) :
    |Y02.serTerm θ n 0 + 2 * Y02.zz θ * (1 - Y02.zz θ) * (1 / n) ^ 2
        + 4 * Y02.zz θ * (1 - Y02.zz θ) * (1 / n) ^ 3| ≤ 10.3 * Y02.zz θ * (1 / n) ^ 4 := by
  have z0 : 0 < Y02.zz θ := (Y02.zz_bounds θ h0 h1).1
  have z1 : Y02.zz θ < 1 := (Y02.zz_bounds θ h0 h1).2
  rw [Y02.serTerm_zero θ n hn]
  have hn0 : 0 < n := by linarith
  have hn2 : 0 < n - 2 := by linarith
  set z := Y02.zz θ
  have e : -(2 * z * (1 - z)) / (n * (n - 2)) + 2 * z * (1 - z) * (1 / n) ^ 2
      + 4 * z * (1 - z) * (1 / n) ^ 3 = -(8 * z * (1 - z)) * ((1 / n) ^ 4 * (n / (n - 2))) := by
    field_simp; ring
  rw [e, abs_mul, abs_neg, abs_of_pos (by nlinarith), abs_of_pos (by positivity)]
  have h1 : n / (n - 2) ≤ 9 / 7 := by rw [div_le_div_iff₀ hn2 (by norm_num)]; linarith
  have h2 : 8 * z * (1 - z) ≤ 8 * z := by nlinarith
  have hp : 0 ≤ (1 / n) ^ 4 := by positivity
  have := mul_le_mul h2 (mul_le_mul_of_nonneg_left h1 hp) (by positivity) (by positivity)
  nlinarith

/-- the `m = 2` term. -/
lemma serTerm_one (θ n : ℝ) (hn : 9 ≤ n) :
    Y02.serTerm θ n 1 = T13.LamInt 2 θ / ((2 * n - 3) * n ^ 2) := by
  have hn1 : n - 1 ≠ 0 := by intro h; linarith
  have hn3 : 2 * n - 3 ≠ 0 := by intro h; linarith
  have hn0 : n ≠ 0 := by intro h; linarith
  simp only [Y02.serTerm, Y02.eps]
  norm_num
  have e : (2:ℝ) - (n - 1)⁻¹ = (2 * n - 3) / (n - 1) := by field_simp; ring
  rw [e]
  field_simp

lemma P6a_bound (θ n : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hn : 9 ≤ n) :
    |Y02.serTerm θ n 1 - Y02.zz θ * (1 - Y02.zz θ) * (1 - 2 * Y02.zz θ) * (1 / n) ^ 3|
      ≤ 1.8 * Y02.zz θ * (1 / n) ^ 4 := by
  have z0 : 0 < Y02.zz θ := (Y02.zz_bounds θ h0 h1).1
  have z1 : Y02.zz θ < 1 := (Y02.zz_bounds θ h0 h1).2
  rw [serTerm_one θ n hn, Lam2]
  have hn0 : 0 < n := by linarith
  have hn3 : 0 < 2 * n - 3 := by linarith
  set z := Y02.zz θ
  have e : 2 * z * (1 - z) * (1 - 2 * z) / ((2 * n - 3) * n ^ 2)
      - z * (1 - z) * (1 - 2 * z) * (1 / n) ^ 3
      = 3 * (z * (1 - z) * (1 - 2 * z)) * ((1 / n) ^ 4 * (n / (2 * n - 3))) := by
    field_simp; ring
  rw [e, abs_mul, abs_mul, abs_of_pos (by norm_num : (0:ℝ) < 3), abs_of_pos (by positivity :
    (0:ℝ) < (1 / n) ^ 4 * (n / (2 * n - 3)))]
  have h1 : n / (2 * n - 3) ≤ 0.6 := by rw [div_le_iff₀ hn3]; linarith
  have h2 : |z * (1 - z) * (1 - 2 * z)| ≤ z := by
    rw [abs_le]; constructor <;> nlinarith [mul_pos z0 (by linarith : (0:ℝ) < 1 - z)]
  have hp : 0 ≤ (1 / n) ^ 4 := by positivity
  have := mul_le_mul h2 (mul_le_mul_of_nonneg_left h1 hp) (by positivity) z0.le
  nlinarith

/-- the tail `m ≥ 3`. -/
lemma P6b_bound (θ n : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hn : 9 ≤ n) {T : ℝ}
    (hT : HasSum (fun i => Y02.serTerm θ n (i + 2)) T) : |T| ≤ 2.6 * Y02.zz θ * (1 / n) ^ 4 := by
  obtain ⟨z0, z1⟩ := Y02.zz_bounds θ h0 h1
  have hz0 : 0 < Y02.zz θ := z0
  obtain ⟨hε0, hε1, hεw, hε2, hε3⟩ := Y02.eps_facts hn
  have hn0 : 0 < n := by linarith
  set z := Y02.zz θ with hzdef
  set ε := Y02.eps n with hεdef
  set w := 1 / n with hw
  have hw0 : 0 < w := by positivity
  have hw19 : w ≤ 1 / 9 := by rw [hw, div_le_div_iff₀ hn0 (by norm_num)]; linarith
  have hw1 : w < 1 := by linarith
  have hg : HasSum (fun i : ℕ => 2 * z * ε * w ^ 3 * w ^ i) (2 * z * ε * w ^ 3 * (1 - w)⁻¹) :=
    (hasSum_geometric_of_lt_one hw0.le hw1).mul_left _
  have hc : 1 - Real.cos θ = 2 * z := by
    rw [hzdef]; have := T13.sin_half_sq θ; unfold Y02.zz; linarith
  have hb : ∀ i : ℕ, |Y02.serTerm θ n (i + 2)| ≤ 2 * z * ε * w ^ 3 * w ^ i := by
    intro i
    have hL2 := (le_min_iff.mp (Y02.stmt2_aux (i + 2 + 1) (by omega) θ)).1
    rw [hc] at hL2
    have hd : 1 ≤ ((i + 2 : ℕ) : ℝ) + 1 - ε := by
      push_cast; linarith [(Nat.cast_nonneg i : (0:ℝ) ≤ i)]
    have hd0 : (0:ℝ) < ((i + 2 : ℕ) : ℝ) + 1 - ε := by linarith
    simp only [Y02.serTerm, ← hεdef]
    rw [abs_div, abs_mul, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
      abs_div, abs_of_pos hε0, abs_of_pos hd0,
      abs_of_pos (by positivity : (0:ℝ) < n ^ (i + 2 + 1))]
    have hf : ε / (((i + 2 : ℕ) : ℝ) + 1 - ε) ≤ ε := div_le_self hε0.le hd
    have hp : 1 / n ^ (i + 2 + 1) = w ^ 3 * w ^ i := by
      rw [hw, ← pow_add, one_div_pow]; ring_nf
    rw [div_eq_mul_one_div, hp]
    have := mul_le_mul hf hL2 (abs_nonneg _) hε0.le
    calc ε / (((i + 2 : ℕ) : ℝ) + 1 - ε) * |T13.LamInt (i + 2 + 1) θ| * (w ^ 3 * w ^ i)
        ≤ ε * (2 * z) * (w ^ 3 * w ^ i) := mul_le_mul_of_nonneg_right this (by positivity)
      _ = 2 * z * ε * w ^ 3 * w ^ i := by ring
  have hup := hasSum_le (fun i => (abs_le.mp (hb i)).2) hT hg
  have hlo := hasSum_le (fun i => by have := (abs_le.mp (hb i)).1; linarith : ∀ i : ℕ,
    -Y02.serTerm θ n (i + 2) ≤ 2 * z * ε * w ^ 3 * w ^ i) hT.neg hg
  have hgeo : 2 * z * ε * w ^ 3 * (1 - w)⁻¹ ≤ 2.6 * z * w ^ 4 := by
    have h1 : (1 - w)⁻¹ ≤ 9 / 8 := by
      rw [inv_le_comm₀ (by linarith) (by norm_num)]
      linarith
    have hε98 : ε ≤ 9 / 8 * w := by
      rw [hεdef, Y02.eps, hw, mul_one_div, div_le_div_iff₀ (by linarith) hn0]; linarith
    have h2 : 2 * z * ε * w ^ 3 ≤ 2 * z * (9 / 8 * w) * w ^ 3 := by
      have : 0 ≤ 2 * z * w ^ 3 := by positivity
      have := mul_le_mul_of_nonneg_left hε98 this
      linarith
    have := mul_le_mul h2 h1 (inv_nonneg.mpr (by linarith)) (by positivity)
    have hp : 0 ≤ z * w ^ 4 := by positivity
    calc _ ≤ 2 * z * (9 / 8 * w) * w ^ 3 * (9 / 8) := this
      _ = 2.53125 * (z * w ^ 4) := by ring
      _ ≤ 2.6 * z * w ^ 4 := by nlinarith
  rw [abs_le]; constructor <;> linarith

end Z06h


namespace Z06h
open Real

lemma P1_core {ε x lam q : ℝ} (hlam : 0 ≤ lam) (hq1 : 1 - ε * lam ≤ q)
    (hq2 : q ≤ 1 - ε * lam + (ε * lam) ^ 2 / 2)
    (A0 : 0 ≤ ε ^ 2 - x ^ 2 - 2 * x ^ 3) (A1 : ε ^ 2 - x ^ 2 - 2 * x ^ 3 ≤ 3.6 * x ^ 4)
    (B0 : 0 ≤ ε ^ 3 - x ^ 3) (B1 : ε ^ 3 - x ^ 3 ≤ 4 * x ^ 4) (C : ε ^ 4 ≤ 1.61 * x ^ 4) :
    |q * ε ^ 2 - x ^ 2 - 2 * x ^ 3 + lam * x ^ 3|
      ≤ 3.6 * x ^ 4 + 4 * lam * x ^ 4 + 0.805 * lam ^ 2 * x ^ 4 := by
  set D := q - 1 + ε * lam with hD
  have D0 : 0 ≤ D := by linarith
  have D1 : D ≤ (ε * lam) ^ 2 / 2 := by linarith
  have e : q * ε ^ 2 - x ^ 2 - 2 * x ^ 3 + lam * x ^ 3
      = (ε ^ 2 - x ^ 2 - 2 * x ^ 3) - lam * (ε ^ 3 - x ^ 3) + D * ε ^ 2 := by rw [hD]; ring
  rw [e]
  have hB : lam * (ε ^ 3 - x ^ 3) ≤ lam * (4 * x ^ 4) := mul_le_mul_of_nonneg_left B1 hlam
  have hB0 : 0 ≤ lam * (ε ^ 3 - x ^ 3) := mul_nonneg hlam B0
  have hC0 : 0 ≤ D * ε ^ 2 := mul_nonneg D0 (sq_nonneg _)
  have hC1 : D * ε ^ 2 ≤ (ε * lam) ^ 2 / 2 * ε ^ 2 := mul_le_mul_of_nonneg_right D1 (sq_nonneg _)
  have hC2 : (ε * lam) ^ 2 / 2 * ε ^ 2 = lam ^ 2 / 2 * ε ^ 4 := by ring
  have hC3 : lam ^ 2 / 2 * ε ^ 4 ≤ lam ^ 2 / 2 * (1.61 * x ^ 4) :=
    mul_le_mul_of_nonneg_left C (by positivity)
  rw [abs_le]; constructor <;> nlinarith

lemma P2_core {ε x lam q : ℝ} (hlam : 0 ≤ lam) (hq1 : 1 - ε * lam ≤ q) (hq : q ≤ 1)
    (hε : 0 ≤ ε) (B0 : 0 ≤ ε ^ 3 - x ^ 3) (B1 : ε ^ 3 - x ^ 3 ≤ 4 * x ^ 4)
    (C : ε ^ 4 ≤ 1.61 * x ^ 4) :
    |q * ε ^ 3 - x ^ 3| ≤ 4 * x ^ 4 + 1.61 * lam * x ^ 4 := by
  have e : q * ε ^ 3 - x ^ 3 = (ε ^ 3 - x ^ 3) - (1 - q) * ε ^ 3 := by ring
  rw [e]
  have h0 : 0 ≤ (1 - q) * ε ^ 3 := mul_nonneg (by linarith) (by positivity)
  have h1 : (1 - q) * ε ^ 3 ≤ ε * lam * ε ^ 3 := mul_le_mul_of_nonneg_right (by linarith) (by positivity)
  have h2 : ε * lam * ε ^ 3 = lam * ε ^ 4 := by ring
  have h3 : lam * ε ^ 4 ≤ lam * (1.61 * x ^ 4) := mul_le_mul_of_nonneg_left C hlam
  rw [abs_le]; constructor <;> nlinarith

lemma scal {lam L : ℝ} (h0 : 0 ≤ lam) (h1 : lam ≤ L + 1) (hL : 2 ≤ L) :
    3.6 + 4 * lam + 0.805 * lam ^ 2 ≤ 5.72 * L ^ 2 ∧ 4 + 1.61 * lam ≤ 2.3 * L ^ 2 ∧
    1 ≤ L ^ 2 / 4 := by
  have h2 : lam ≤ 1.5 * L := by linarith
  have h3 : lam ^ 2 ≤ 2.25 * L ^ 2 := by nlinarith
  refine ⟨by nlinarith, by nlinarith, by nlinarith⟩

end Z06h


namespace Z06h
open Real

lemma a2_le (θ : ℝ) (h0 : 0 < θ) (h1 : θ < π) : Y02.a2 θ ≤ 2.16 := by
  have hb1 := (Y02.proof_9 θ h0 h1).2.2.2.2.2
  have hpi := Real.pi_lt_d2
  nlinarith [Real.pi_pos]

lemma p10 : ∀ θ n : ℝ, 0 < θ → θ < π → 9 ≤ n →
    |Y02.E θ n - (Z06.d θ - 2 * Y02.zz θ * Y02.a2 θ * Real.log n) / n ^ 3|
      ≤ 3.1e4 * Y02.zz θ * Real.log n ^ 2 / n ^ 4 := by
  intro θ n h0 h1 hn
  have z0 : 0 < Y02.zz θ := (Y02.zz_bounds θ h0 h1).1
  have z1 : Y02.zz θ < 1 := (Y02.zz_bounds θ h0 h1).2
  obtain ⟨-, -, ha1, ha2, hb0, hb1⟩ := Y02.proof_9 θ h0 h1
  obtain ⟨T4, hT, hE⟩ := decomp θ n h0 h1 hn
  obtain ⟨hq0, hq1, hqe⟩ := q_facts θ n h0 h1 hn
  obtain ⟨hε0, hε1, -, -, -⟩ := Y02.eps_facts hn
  have hρ := Y02.proof_12 θ (Y02.eps n) h0 h1 hε0 hε1
  have hε : Y02.eps n = 1 / (n - 1) := rfl
  rw [hε] at hρ hε0
  have hT2 : HasSum (fun i => Y02.serTerm θ n (i + 2)) (T4 - Y02.serTerm θ n 1) := by
    have := (hasSum_nat_add_iff' 1).mpr hT
    simpa using this
  have t6b := P6b_bound θ n h0 h1 hn hT2
  have t6a := P6a_bound θ n h0 h1 hn
  have t5 := P5_bound θ n h0 h1 hn
  obtain ⟨f1a, f1b⟩ := F1 hn
  obtain ⟨f2a, f2b⟩ := F2 hn
  have f3 := F3 hn
  have ha3 := a3_abs θ h0 h1
  have hL2 := Y02.two_le_log hn
  have hn0 : 0 < n := by linarith
  set L := Real.log n with hLdef
  set lam := L - Y02.a1 θ with hlam
  have hlam0 : 0 ≤ lam := by linarith
  have hlam1 : lam ≤ L + 1 := by linarith
  set ε := 1 / (n - 1) with hεdef
  set x := 1 / n with hx
  set z := Y02.zz θ with hz
  set q := (Y02.mu θ / n) ^ (1 / (n - 1)) with hq
  have ht0 : 0 ≤ ε * lam := mul_nonneg hε0.le hlam0
  obtain ⟨qlo, qhi⟩ := exp_neg_bounds ht0
  rw [← hqe] at qlo qhi
  have b1 := P1_core hlam0 qlo qhi f1a f1b f2a f2b f3
  have b2 := P2_core hlam0 qlo hq1 hε0.le f2a f2b f3
  obtain ⟨s1, s2, s3⟩ := scal hlam0 hlam1 hL2
  set ρ4 := Y02.R θ ε * Real.exp (-(Y02.a1 θ * ε)) / (2 * z)
    - (1 + Y02.a2 θ * ε ^ 2 + Y02.a3 θ * ε ^ 3) with hρ4
  have hρ4' : |ρ4| ≤ 31000 * ε ^ 4 := by
    have e : Y02.R θ ε * Real.exp (-(Y02.a1 θ * ε))
        - 2 * z * (1 + Y02.a2 θ * ε ^ 2 + Y02.a3 θ * ε ^ 3) = 2 * z * ρ4 := by
      rw [hρ4]; field_simp
    rw [e, abs_mul, abs_of_pos (by positivity)] at hρ
    exact le_of_mul_le_mul_left hρ (by positivity)
  set W := z * x ^ 4 with hW
  have hW0 : 0 ≤ W := by positivity
  set B1 := q * ε ^ 2 - x ^ 2 - 2 * x ^ 3 + lam * x ^ 3
  set B2 := q * ε ^ 3 - x ^ 3
  set P5 := Y02.serTerm θ n 0 + 2 * z * (1 - z) * x ^ 2 + 4 * z * (1 - z) * x ^ 3
  set P6 := Y02.serTerm θ n 1 - z * (1 - z) * (1 - 2 * z) * x ^ 3
  have key : Y02.E θ n - (Z06.d θ - 2 * z * Y02.a2 θ * L) / n ^ 3
      = 2 * z * Y02.a2 θ * B1 + 2 * z * Y02.a3 θ * B2 + 2 * z * q * ρ4 + P5 + P6
        + (T4 - Y02.serTerm θ n 1) := by
    rw [hE, Z06.d]
    simp only [B1, B2, P5, P6, hρ4, hlam, hx, hz]
    field_simp
    ring
  rw [key]
  -- the pieces
  have ha216 : Y02.a2 θ ≤ 2.16 := a2_le θ h0 h1
  have h3 : W ≤ W * L ^ 2 / 4 := by
    have := mul_le_mul_of_nonneg_left s3 hW0
    linarith [show W * (L ^ 2 / 4) = W * L ^ 2 / 4 by ring]
  have c1 : |2 * z * Y02.a2 θ * B1| ≤ 4.32 * 5.72 * (W * L ^ 2) := by
    rw [abs_mul, abs_of_pos (by positivity)]
    have hb : |B1| ≤ x ^ 4 * (5.72 * L ^ 2) := by
      have : 3.6 * x ^ 4 + 4 * lam * x ^ 4 + 0.805 * lam ^ 2 * x ^ 4
          = x ^ 4 * (3.6 + 4 * lam + 0.805 * lam ^ 2) := by ring
      rw [this] at b1
      exact b1.trans (mul_le_mul_of_nonneg_left s1 (by positivity))
    have hza : 2 * z * Y02.a2 θ ≤ 4.32 * z := by
      have := mul_le_mul_of_nonneg_left ha216 (by positivity : (0:ℝ) ≤ 2 * z)
      linarith
    have := mul_le_mul hza hb (abs_nonneg _) (by positivity)
    calc _ ≤ 4.32 * z * (x ^ 4 * (5.72 * L ^ 2)) := this
      _ = _ := by rw [hW]; ring
  have c2 : |2 * z * Y02.a3 θ * B2| ≤ 18 * 2.3 * (W * L ^ 2) := by
    rw [abs_mul, abs_mul, abs_of_pos (by positivity : (0:ℝ) < 2 * z)]
    have hb : |B2| ≤ x ^ 4 * (2.3 * L ^ 2) := by
      have : 4 * x ^ 4 + 1.61 * lam * x ^ 4 = x ^ 4 * (4 + 1.61 * lam) := by ring
      rw [this] at b2
      exact b2.trans (mul_le_mul_of_nonneg_left s2 (by positivity))
    have := mul_le_mul (mul_le_mul_of_nonneg_left ha3 (by positivity : (0:ℝ) ≤ 2 * z)) hb
      (abs_nonneg _) (by positivity)
    calc _ ≤ 2 * z * 9 * (x ^ 4 * (2.3 * L ^ 2)) := this
      _ = _ := by rw [hW]; ring
  have c3 : |2 * z * q * ρ4| ≤ 24955 * (W * L ^ 2) := by
    rw [abs_mul, abs_of_pos (by positivity : (0:ℝ) < 2 * z * q)]
    have h1 : 2 * z * q ≤ 2 * z := mul_le_of_le_one_right (by positivity) hq1
    have h2 : |ρ4| ≤ 31000 * (1.61 * x ^ 4) := hρ4'.trans (by linarith)
    have := mul_le_mul h1 h2 (abs_nonneg _) (by positivity)
    calc _ ≤ 2 * z * (31000 * (1.61 * x ^ 4)) := this
      _ = 99820 * W := by rw [hW]; ring
      _ ≤ 24955 * (W * L ^ 2) := by linarith
  have c4 : |P5| + |P6| + |T4 - Y02.serTerm θ n 1| ≤ 3.7 * (W * L ^ 2) := by
    have : 10.3 * z * x ^ 4 + 1.8 * z * x ^ 4 + 2.6 * z * x ^ 4 = 14.7 * W := by rw [hW]; ring
    linarith
  have u1 := abs_add_le (2 * z * Y02.a2 θ * B1 + 2 * z * Y02.a3 θ * B2 + 2 * z * q * ρ4 + P5 + P6)
    (T4 - Y02.serTerm θ n 1)
  have u2 := abs_add_le (2 * z * Y02.a2 θ * B1 + 2 * z * Y02.a3 θ * B2 + 2 * z * q * ρ4 + P5) P6
  have u3 := abs_add_le (2 * z * Y02.a2 θ * B1 + 2 * z * Y02.a3 θ * B2 + 2 * z * q * ρ4) P5
  have u4 := abs_add_le (2 * z * Y02.a2 θ * B1 + 2 * z * Y02.a3 θ * B2) (2 * z * q * ρ4)
  have u5 := abs_add_le (2 * z * Y02.a2 θ * B1) (2 * z * Y02.a3 θ * B2)
  have hfin : 3.1e4 * z * L ^ 2 / n ^ 4 = 31000 * (W * L ^ 2) := by
    rw [hW, hx]; field_simp; norm_num
  rw [hfin]
  have hWL : 0 ≤ W * L ^ 2 := by positivity
  linarith

end Z06h

namespace Z06h
open Real MeasureTheory Set Filter Topology

lemma S12int_eq_Iint (t : ℝ) : Z06.S12int t = (1 / 2) * Y02.Iint t := by
  unfold Z06.S12int Y02.Iint
  simp_rw [Y02.integrand_eq]

lemma hasDerivAt_S12int {w : ℝ} (hw : w < 1) :
    HasDerivAt Z06.S12int ((1 / 2) * (D00.G w * Real.log (1 - w))) w := by
  have e : Z06.S12int = fun t => (1 / 2) * Y02.Iint t := funext S12int_eq_Iint
  rw [e]
  exact (Y02.hasDerivAt_Iint hw).const_mul _

lemma S12int_zero : Z06.S12int 0 = 0 := by simp [Z06.S12int]

/-- `t ↦ S_{1,2}(t)/t` is continuous on `[0, 1)` (value `0` at `0`). -/
lemma q_contOn {y : ℝ} (hy : y < 1) : ContinuousOn (fun t => Z06.S12int t / t) (Icc 0 y) := by
  intro t ht
  rcases eq_or_lt_of_le ht.1 with h | h
  · rw [← h]
    rw [ContinuousWithinAt, S12int_zero, zero_div]
    have h12 : ∀ᶠ u in 𝓝[Icc 0 y] (0:ℝ), u < 1 / 2 :=
      nhdsWithin_le_nhds (gt_mem_nhds (by norm_num))
    have ht0 : Tendsto (fun u : ℝ => u) (𝓝[Icc 0 y] (0:ℝ)) (𝓝 0) :=
      tendsto_nhdsWithin_of_tendsto_nhds (continuous_id.tendsto' 0 0 rfl)
    apply squeeze_zero_norm' _ ht0
    filter_upwards [self_mem_nhdsWithin, h12] with u hu h2
    rcases eq_or_lt_of_le hu.1 with h' | h'
    · rw [← h']; simp
    · rw [Real.norm_eq_abs, abs_div, abs_of_nonneg (S12int_nonneg hu.1), abs_of_pos h',
        div_le_iff₀ h']
      have := S12int_small hu.1 h2.le
      nlinarith
  · have hc : ContinuousAt (fun t => Z06.S12int t / t) t :=
      (hasDerivAt_S12int (by linarith [ht.2])).continuousAt.div continuousAt_id h.ne'
    exact hc.continuousWithinAt

lemma hasDerivAt_S22 {y : ℝ} (h0 : 0 < y) (h1 : y < 1) :
    HasDerivAt Z06.S22 (Z06.S12int y / y) y := by
  have hint : IntervalIntegrable (fun t => Z06.S12int t / t) volume 0 y := by
    apply ContinuousOn.intervalIntegrable
    rw [uIcc_of_le h0.le]; exact q_contOn h1
  have hco : ContinuousOn (fun t => Z06.S12int t / t) (Ioo 0 1) := by
    intro t ht
    exact ((hasDerivAt_S12int ht.2).continuousAt.div continuousAt_id ht.1.ne').continuousWithinAt
  exact intervalIntegral.integral_hasDerivAt_right hint
    (hco.stronglyMeasurableAtFilter isOpen_Ioo y ⟨h0, h1⟩)
    (hco.continuousAt (isOpen_Ioo.mem_nhds ⟨h0, h1⟩))

lemma S22_zero : Z06.S22 0 = 0 := by simp [Z06.S22]

/-- `N(z) = Σ h_{k+1} z^{k+1}/(k+1)² = Li_{2,2}(z)`. -/
noncomputable def Nser (z : ℝ) : ℝ := ∑' k : ℕ, Y02.hk k * z ^ (k + 1) / ((k:ℝ) + 1) ^ 2

lemma hasDerivAt_Nser {y : ℝ} (hy : |y| < 1) :
    HasDerivAt Nser (∑' k : ℕ, Y02.hk k * y ^ k / ((k:ℝ) + 1)) y := by
  set r := (|y| + 1) / 2 with hr
  have hr1 : r < 1 := by rw [hr]; linarith
  have hr0 : 0 ≤ r := by rw [hr]; positivity
  have ht : y ∈ Ioo (-r) r := ⟨by linarith [neg_abs_le y], by linarith [le_abs_self y]⟩
  exact hasDerivAt_tsum_of_isPreconnected (t := Ioo (-r) r) (y₀ := 0)
    (g := fun (n : ℕ) (z : ℝ) => Y02.hk n * z ^ (n + 1) / ((n : ℝ) + 1) ^ 2)
    (g' := fun (n : ℕ) (z : ℝ) => Y02.hk n * z ^ n / ((n : ℝ) + 1)) (u := fun n : ℕ => 2 * r ^ n)
    ((summable_geometric_of_lt_one hr0 hr1).mul_left 2) isOpen_Ioo isPreconnected_Ioo
    (by
      intro n z _
      have := ((hasDerivAt_pow (n + 1) z).const_mul (Y02.hk n)).div_const (((n : ℝ) + 1) ^ 2)
      refine this.congr_deriv ?_
      push_cast
      field_simp)
    (by
      intro n z hz
      rw [Real.norm_eq_abs, abs_div, abs_mul, abs_of_nonneg (Y02.hk_nonneg n), abs_pow,
        abs_of_pos (by positivity : (0:ℝ) < (n:ℝ) + 1)]
      have hz' : |z| ≤ r := abs_le.mpr ⟨hz.1.le, hz.2.le⟩
      rw [div_le_iff₀ (by positivity)]
      have h1 := mul_le_mul (Y02.hk_le n) (pow_le_pow_left₀ (abs_nonneg z) hz' n) (by positivity)
        (by norm_num)
      have h2 : 2 * r ^ n ≤ 2 * r ^ n * ((n:ℝ) + 1) := by
        have : 0 ≤ 2 * r ^ n := by positivity
        nlinarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
      linarith)
    ⟨by linarith [abs_nonneg y, hr], by linarith [abs_nonneg y, hr]⟩
    (by simp)
    ht

lemma Nser_deriv_eq {y : ℝ} (h0 : y ≠ 0) :
    ∑' k : ℕ, Y02.hk k * y ^ k / ((k:ℝ) + 1) = Y02.Mser y / y := by
  unfold Y02.Mser
  rw [eq_div_iff h0, ← tsum_mul_right]
  congr 1; funext k; ring

lemma Nser_eq {z : ℝ} (h0 : 0 ≤ z) (h1 : z < 1) :
    Nser z = (1 / 2) * D00.Li2 z ^ 2 - 2 * Z06.S22 z := by
  set F : ℝ → ℝ := fun w => Nser w - ((1 / 2) * D00.Li2 w ^ 2 - 2 * Z06.S22 w) with hF
  have hF0 : F 0 = 0 := by simp [hF, Nser, S22_zero]
  suffices F 0 = F z by simp only [hF] at this hF0; linarith
  apply D00.eq_of_deriv_zero h0
  · intro w hw
    have hw1 : w < 1 := by linarith [hw.2]
    have hwa : |w| < 1 := by rw [abs_lt]; constructor <;> linarith [hw.1]
    have hS : ContinuousWithinAt Z06.S22 (Icc 0 z) w := by
      rcases eq_or_lt_of_le hw.1 with h | h
      · rw [← h, ContinuousWithinAt, S22_zero]
        have ht0 : Tendsto (fun u : ℝ => 4 * u) (𝓝[Icc 0 z] (0:ℝ)) (𝓝 0) := by
          have hc : Continuous (fun u : ℝ => 4 * u) := by fun_prop
          have := (hc.tendsto (0:ℝ)).mono_left (nhdsWithin_le_nhds (s := Icc 0 z))
          simpa using this
        apply squeeze_zero_norm' _ ht0
        filter_upwards [self_mem_nhdsWithin] with u hu
        rw [Real.norm_eq_abs]
        exact S22_abs hu.1 (by linarith [hu.2])
      · exact (hasDerivAt_S22 h hw1).continuousAt.continuousWithinAt
    exact ((hasDerivAt_Nser hwa).continuousAt.continuousWithinAt).sub
      (((D00.continuous_Li2.continuousAt.pow 2).const_mul _).continuousWithinAt.sub
        (hS.const_mul _))
  · intro w hw
    have hw1 : w < 1 := by linarith [hw.2]
    have hwa : |w| < 1 := by rw [abs_lt]; constructor <;> linarith [hw.1]
    have hw0 : w ≠ 0 := hw.1.ne'
    have dN := hasDerivAt_Nser hwa
    rw [Nser_deriv_eq hw0, Y02.Mser_eq hw.1.le hw1] at dN
    have dL := ((D00.hasDerivAt_Li2 hw1 hw0).pow 2).const_mul (1 / 2)
    have dS := (hasDerivAt_S22 hw.1 hw1).const_mul 2
    have hall := dN.sub (dL.sub dS)
    convert hall using 1
    rw [S12int_eq_Iint]
    field_simp
    ring

end Z06h

namespace Z06h
open Real

lemma Nser_hasSum {z : ℝ} (h0 : 0 ≤ z) (h1 : z < 1) :
    HasSum (fun k : ℕ => Y02.hk k * z ^ (k + 1) / ((k:ℝ) + 1) ^ 2) (Nser z) := by
  apply Summable.hasSum
  refine Summable.of_nonneg_of_le (fun k => by have := Y02.hk_nonneg k; positivity) (fun k => ?_)
    ((summable_geometric_of_lt_one h0 h1).mul_left 2)
  rw [div_le_iff₀ (by positivity)]
  have h2 := Y02.hk_le k
  have hz : z ^ (k + 1) ≤ z ^ k := pow_le_pow_of_le_one h0 h1.le (by omega)
  have hk1 : (1:ℝ) ≤ ((k:ℝ) + 1) ^ 2 := by nlinarith [(Nat.cast_nonneg k : (0:ℝ) ≤ k)]
  have := mul_le_mul h2 hz (by positivity) (by norm_num)
  have h3 : 0 ≤ 2 * z ^ k := by positivity
  nlinarith

/-- `f_4 = Σ h_{k+1} z^{k+1}/((k+1)²(k+2)) = Li_{2,2}(z) - f_3`. -/
lemma f4_hasSum (θ : ℝ) (h0 : 0 < θ) (h1 : θ < π) :
    HasSum (fun k => Y02.hk k * (Y02.zz θ ^ (k + 1) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2))))
      ((1 / 2) * D00.Li2 (Y02.zz θ) ^ 2 - 2 * Z06.S22 (Y02.zz θ)
        - (2 * Y02.ww θ * Y02.S12 θ + (1 + Y02.U θ) * (D00.Li2 (Y02.zz θ) - 1))) := by
  have z0 : 0 < Y02.zz θ := (Y02.zz_bounds θ h0 h1).1
  have z1 : Y02.zz θ < 1 := (Y02.zz_bounds θ h0 h1).2
  have hN := Nser_hasSum z0.le z1
  rw [Nser_eq z0.le z1] at hN
  have hf3 := Y02.f3_hasSum θ h0 h1
  convert hN.sub hf3 using 1
  funext k
  field_simp
  ring

/-- `Λ_ε = 2z(1 - aε + (a - Li₂ z)ε² + f₃ε³ + f₄ε⁴ + r)`, `|r| ≤ 3ε⁵`. -/
lemma Lam_expand4 (θ ε : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hε0 : 0 < ε) (hε1 : ε ≤ 1 / 8) :
    ∃ r : ℝ, |r| ≤ 3 * ε ^ 5 ∧
      Y02.Lam ε θ = 2 * Y02.zz θ * (1 - (1 + Y02.U θ) * ε
        + ((1 + Y02.U θ) - D00.Li2 (Y02.zz θ)) * ε ^ 2
        + (2 * Y02.ww θ * Y02.S12 θ + (1 + Y02.U θ) * (D00.Li2 (Y02.zz θ) - 1)) * ε ^ 3
        + ((1 / 2) * D00.Li2 (Y02.zz θ) ^ 2 - 2 * Z06.S22 (Y02.zz θ)
          - (2 * Y02.ww θ * Y02.S12 θ + (1 + Y02.U θ) * (D00.Li2 (Y02.zz θ) - 1))) * ε ^ 4
        + r) := by
  obtain ⟨z0, z1⟩ := Y02.zz_bounds θ h0 h1
  have hzz : Y02.zz θ = Real.sin (θ / 2) ^ 2 := rfl
  have hf3 := Y02.f3_hasSum θ h0 h1
  have hf4 := f4_hasSum θ h0 h1
  have hUb : -1 < Y02.U θ ∧ Y02.U θ < 0 := Y02.U_bounds z0 z1
  set V := 2 * Y02.ww θ * Y02.S12 θ + (1 + Y02.U θ) * (D00.Li2 (Y02.zz θ) - 1) with hV
  set W := (1 / 2) * D00.Li2 (Y02.zz θ) ^ 2 - 2 * Z06.S22 (Y02.zz θ) - V with hW
  set z := Y02.zz θ with hzdef
  have zpos : 0 < z := z0
  have hA : HasSum (fun k : ℕ => z ^ (k + 1) / (((k:ℝ) + 1) * ((k:ℝ) + 2))) (1 + Y02.U θ) :=
    Y02.sumA z0 z1
  have hB : HasSum (fun k : ℕ => z ^ (k + 1) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2)))
      (D00.Li2 z - 1 - Y02.U θ) := Y02.sumB z0 z1
  set fA := fun k : ℕ => z ^ (k + 1) / (((k:ℝ) + 1) * ((k:ℝ) + 2)) with hfA
  set fB := fun k : ℕ => z ^ (k + 1) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2)) with hfB
  set cc := fun k : ℕ => z ^ (k + 1) * (ε * ((k:ℝ) + 1 + ε) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2)))
    with hcc
  set rk := fun k : ℕ => -(cc k * (Y02.pp ε k - 1 + ε ^ 2 * Y02.hk k)) with hrk
  have hz0 : z ≠ 0 := zpos.ne'
  have hsplit : ∀ k : ℕ, T13.lamTerm ε z (k + 1) =
      2 * z * (-ε * fA k - ε ^ 2 * fB k + ε ^ 3 * (Y02.hk k * fA k)
        + ε ^ 4 * (Y02.hk k * fB k) + rk k) := by
    intro k
    rw [Y02.lamTerm_closed]
    simp only [hfA, hfB, hrk, hcc]
    field_simp
    ring
  have hsum : Summable (fun k : ℕ => T13.lamTerm ε z (k + 1)) :=
    (summable_nat_add_iff 1).mpr (Y02.lamTerm_summable ε z z0.le z1)
  have hrs : Summable rk := by
    have : Summable (fun k => T13.lamTerm ε z (k + 1) / (2 * z) + ε * fA k + ε ^ 2 * fB k
        - ε ^ 3 * (Y02.hk k * fA k) - ε ^ 4 * (Y02.hk k * fB k)) :=
      ((((hsum.div_const _).add (hA.summable.mul_left ε)).add
        (hB.summable.mul_left (ε ^ 2))).sub (hf3.summable.mul_left (ε ^ 3))).sub
        (hf4.summable.mul_left (ε ^ 4))
    refine this.congr (fun k => ?_)
    rw [hsplit]; field_simp; ring
  have hrk1 : ∀ k, |rk k| ≤ 3 * ε ^ 5 * fA k := fun k => by
    obtain ⟨hp1, hp2⟩ := Y02.pp_bounds2 (ε := ε) (by nlinarith) k
    have hk0 := Y02.hk_nonneg k
    have hk2 := Y02.hk_le k
    have hzk : 0 ≤ z ^ (k + 1) := pow_nonneg zpos.le _
    have hfA0 : 0 ≤ fA k := by simp only [hfA]; positivity
    have hfrac : ε * ((k:ℝ) + 1 + ε) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2)) ≤
        9 / 8 * ε / (((k:ℝ) + 1) * ((k:ℝ) + 2)) := by
      rw [div_le_div_iff₀ (by positivity) (by positivity)]
      have : 0 ≤ ε * ((k:ℝ) + 1) * ((k:ℝ) + 2) := by positivity
      nlinarith
    have hc0 : 0 ≤ cc k := by
      simp only [hcc]; exact mul_nonneg hzk (div_nonneg (by positivity) (by positivity))
    have hc1 : cc k ≤ 9 / 8 * ε * fA k := by
      simp only [hcc, hfA]
      calc z ^ (k + 1) * (ε * ((k:ℝ) + 1 + ε) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2)))
          ≤ z ^ (k + 1) * (9 / 8 * ε / (((k:ℝ) + 1) * ((k:ℝ) + 2))) :=
            mul_le_mul_of_nonneg_left hfrac hzk
        _ = 9 / 8 * ε * (z ^ (k + 1) / (((k:ℝ) + 1) * ((k:ℝ) + 2))) := by ring
    have hd0 : 0 ≤ Y02.pp ε k - 1 + ε ^ 2 * Y02.hk k := by linarith
    have hd1 : Y02.pp ε k - 1 + ε ^ 2 * Y02.hk k ≤ 2 * ε ^ 4 := by
      have : ε ^ 4 * Y02.hk k ^ 2 / 2 ≤ 2 * ε ^ 4 := by
        have : Y02.hk k ^ 2 ≤ 4 := by nlinarith
        have h4 : 0 ≤ ε ^ 4 := by positivity
        nlinarith
      linarith
    have t3 : 0 ≤ cc k * (Y02.pp ε k - 1 + ε ^ 2 * Y02.hk k) := mul_nonneg hc0 hd0
    have t4 : cc k * (Y02.pp ε k - 1 + ε ^ 2 * Y02.hk k) ≤ 9 / 8 * ε * fA k * (2 * ε ^ 4) :=
      mul_le_mul hc1 hd1 hd0 (by positivity)
    have t5 : 9 / 8 * ε * fA k * (2 * ε ^ 4) ≤ 3 * ε ^ 5 * fA k := by
      have : 0 ≤ ε ^ 5 * fA k := by positivity
      have e : 9 / 8 * ε * fA k * (2 * ε ^ 4) = 9 / 4 * (ε ^ 5 * fA k) := by ring
      rw [e]; nlinarith
    simp only [hrk]
    rw [abs_neg, abs_of_nonneg t3]
    linarith
  refine ⟨∑' k, rk k, ?_, ?_⟩
  · have hup := hasSum_le (fun k => (abs_le.mp (hrk1 k)).2) hrs.hasSum (hA.mul_left (3 * ε ^ 5))
    have hlo := hasSum_le (fun k => by have := (abs_le.mp (hrk1 k)).1; linarith :
      ∀ k, -rk k ≤ 3 * ε ^ 5 * fA k) hrs.hasSum.neg (hA.mul_left (3 * ε ^ 5))
    have h5 : 3 * ε ^ 5 * (1 + Y02.U θ) ≤ 3 * ε ^ 5 * 1 :=
      mul_le_mul_of_nonneg_left (by linarith [hUb.2]) (by positivity)
    have h6 : 0 ≤ 3 * ε ^ 5 * (1 + Y02.U θ) := mul_nonneg (by positivity) (by linarith [hUb.1])
    rw [abs_le]; constructor <;> linarith
  · have H : HasSum (fun k : ℕ => T13.lamTerm ε z (k + 1))
        (2 * z * (-ε * (1 + Y02.U θ) - ε ^ 2 * (D00.Li2 z - 1 - Y02.U θ) + ε ^ 3 * V
          + ε ^ 4 * W + ∑' k, rk k)) := by
      have := (((((hA.mul_left (-ε)).sub (hB.mul_left (ε ^ 2))).add
        (hf3.mul_left (ε ^ 3))).add (hf4.mul_left (ε ^ 4))).add hrs.hasSum).mul_left (2 * z)
      convert this using 1
      funext k; rw [hsplit]
    have H2 := (hasSum_nat_add_iff (f := fun k => T13.lamTerm ε z k) 1).mp H
    unfold Y02.Lam
    rw [← hzz, H2.tsum_eq]
    simp only [Finset.sum_range_one, Y02.lamTerm_zero]
    ring

end Z06h

namespace Z06h
open Real

lemma csin_alg (x : ℂ) : Complex.sin x - (x - x ^ 3 / 6 + x ^ 5 / 120) =
    (Complex.exp (-x * Complex.I) - ∑ m ∈ Finset.range 7, (-x * Complex.I) ^ m / m.factorial)
      * Complex.I / 2 -
    (Complex.exp (x * Complex.I) - ∑ m ∈ Finset.range 7, (x * Complex.I) ^ m / m.factorial)
      * Complex.I / 2 := by
  simp only [Complex.sin, Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial]
  push_cast
  linear_combination (-x - x ^ 3 * (Complex.I ^ 2 - 1) / 6
    - x ^ 5 * (Complex.I ^ 4 - Complex.I ^ 2 + 1) / 120) * Complex.I_sq

lemma csin_bound7 {x : ℂ} (hx : ‖x‖ ≤ 1) :
    ‖Complex.sin x - (x - x ^ 3 / 6 + x ^ 5 / 120)‖ ≤ ‖x‖ ^ 7 / 4410 := by
  rw [csin_alg]
  have h1 := Complex.exp_bound (x := -x * Complex.I) (n := 7) (by simpa using hx) (by norm_num)
  have h2 := Complex.exp_bound (x := x * Complex.I) (n := 7) (by simpa using hx) (by norm_num)
  have n1 : ‖-x * Complex.I‖ = ‖x‖ := by simp
  have n2 : ‖x * Complex.I‖ = ‖x‖ := by simp
  rw [n1] at h1; rw [n2] at h2
  set A := Complex.exp (-x * Complex.I) - ∑ m ∈ Finset.range 7, (-x * Complex.I) ^ m / m.factorial
  set B := Complex.exp (x * Complex.I) - ∑ m ∈ Finset.range 7, (x * Complex.I) ^ m / m.factorial
  have e : A * Complex.I / 2 - B * Complex.I / 2 = (A - B) * (Complex.I / 2) := by ring
  rw [e, norm_mul]
  have hI : ‖Complex.I / 2‖ = 1 / 2 := by simp
  rw [hI]
  have := norm_sub_le A B
  have hc : ((Nat.succ 7 : ℕ) : ℝ) * ((Nat.factorial 7 : ℕ) * (7 : ℕ) : ℝ)⁻¹ = 1 / 4410 := by
    norm_num [Nat.factorial]
  rw [hc] at h1 h2
  nlinarith [norm_nonneg x, pow_nonneg (norm_nonneg x) 7]

lemma sin_bound7 {x : ℝ} (hx : |x| ≤ 1) :
    |Real.sin x - (x - x ^ 3 / 6 + x ^ 5 / 120)| ≤ |x| ^ 7 / 4410 := by
  have h := csin_bound7 (x := (x : ℂ)) (by rw [Complex.norm_real]; exact hx)
  rw [← Complex.ofReal_sin] at h
  have e : (Real.sin x : ℂ) - ((x : ℂ) - (x : ℂ) ^ 3 / 6 + (x : ℂ) ^ 5 / 120)
      = ((Real.sin x - (x - x ^ 3 / 6 + x ^ 5 / 120) : ℝ) : ℂ) := by push_cast; ring
  rw [e, Complex.norm_real, Complex.norm_real, Real.norm_eq_abs, Real.norm_eq_abs] at h
  exact h

end Z06h

namespace Z06h
open Real

/-- `|πε/sin(πε) - (1 + π²ε²/6 + 7π⁴ε⁴/360)| ≤ ε⁵` for `0 < ε ≤ 1/8`. -/
lemma sinc_bound6 {ε : ℝ} (hε0 : 0 < ε) (hε1 : ε ≤ 1 / 8) :
    |π * ε / Real.sin (π * ε) - (1 + π ^ 2 / 6 * ε ^ 2 + 7 * π ^ 4 / 360 * ε ^ 4)| ≤ ε ^ 5 := by
  have hpi := Real.pi_gt_three
  have hpi' := Real.pi_lt_d2
  have hx0 : 0 < π * ε := by positivity
  have hx1 : π * ε ≤ 2 / 5 := by nlinarith
  have hb := sin_bound7 (x := π * ε) (by rw [abs_of_pos hx0]; linarith)
  rw [abs_of_pos hx0] at hb
  have e : 1 + π ^ 2 / 6 * ε ^ 2 + 7 * π ^ 4 / 360 * ε ^ 4
      = 1 + (π * ε) ^ 2 / 6 + 7 * (π * ε) ^ 4 / 360 := by ring
  have hx5 : (π * ε) ^ 6 ≤ 1000 * ε ^ 6 := by
    rw [mul_pow]
    have h6 := pow_lt_pow_left₀ hpi' Real.pi_pos.le (by norm_num : (6:ℕ) ≠ 0)
    have : π ^ 6 ≤ 1000 := by linarith [h6, show (3.15:ℝ) ^ 6 ≤ 1000 by norm_num]
    exact mul_le_mul_of_nonneg_right this (by positivity)
  rw [e]
  generalize π * ε = x at hx0 hx1 hb hx5 ⊢
  set sx := Real.sin x with hsx
  set δ := sx - (x - x ^ 3 / 6 + x ^ 5 / 120) with hδ
  have hδ' : |δ| ≤ x ^ 7 / 4410 := hb
  have hδ2 := abs_le.mp hδ'
  have hx2 : x ^ 2 ≤ 4 / 25 := by nlinarith
  have hx7 : x ^ 7 ≤ x ^ 3 * (4 / 25) ^ 2 := by
    have : x ^ 7 = x ^ 3 * (x ^ 2) ^ 2 := by ring
    rw [this]
    exact mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (by positivity) hx2 2) (by positivity)
  have hsx0 : 9 / 10 * x ≤ sx := by
    have h3 : x ^ 3 ≤ x * (4 / 25) := by nlinarith
    have h5 : 0 ≤ x ^ 5 := by positivity
    nlinarith
  have hsxp : 0 < sx := by linarith
  have e2 : x / sx - (1 + x ^ 2 / 6 + 7 * x ^ 4 / 360)
      = (x ^ 7 / 540 - 7 * x ^ 9 / 43200 - δ * (1 + x ^ 2 / 6 + 7 * x ^ 4 / 360)) / sx := by
    field_simp
    rw [hδ]; ring
  rw [e2, abs_div, abs_of_pos hsxp, div_le_iff₀ hsxp]
  have hq : 0 < 1 + x ^ 2 / 6 + 7 * x ^ 4 / 360 := by positivity
  have hq2 : 1 + x ^ 2 / 6 + 7 * x ^ 4 / 360 ≤ 2 := by
    have : x ^ 4 ≤ 1 := by nlinarith
    nlinarith
  have hnum : |x ^ 7 / 540 - 7 * x ^ 9 / 43200 - δ * (1 + x ^ 2 / 6 + 7 * x ^ 4 / 360)|
      ≤ x ^ 7 / 540 + 7 * x ^ 9 / 43200 + x ^ 7 / 4410 * 2 := by
    have h1 : |δ * (1 + x ^ 2 / 6 + 7 * x ^ 4 / 360)| ≤ x ^ 7 / 4410 * 2 := by
      rw [abs_mul, abs_of_pos hq]
      exact mul_le_mul hδ' hq2 hq.le (by positivity)
    have h2 := abs_sub (x ^ 7 / 540 - 7 * x ^ 9 / 43200) (δ * (1 + x ^ 2 / 6 + 7 * x ^ 4 / 360))
    have h3 := abs_sub (x ^ 7 / 540) (7 * x ^ 9 / 43200)
    rw [abs_of_pos (by positivity : (0:ℝ) < x ^ 7 / 540),
      abs_of_pos (by positivity : (0:ℝ) < 7 * x ^ 9 / 43200)] at h3
    linarith
  refine hnum.trans ?_
  have hx9 : x ^ 9 ≤ x ^ 7 := by
    have : x ^ 9 = x ^ 7 * x ^ 2 := by ring
    rw [this]; exact mul_le_of_le_one_right (by positivity) (by linarith)
  have h1 : x ^ 7 / 540 + 7 * x ^ 9 / 43200 + x ^ 7 / 4410 * 2 ≤ x ^ 7 * (3 / 1000) := by
    nlinarith [pow_pos hx0 7]
  have h2 : x ^ 7 * (3 / 1000) ≤ ε ^ 5 * (9 / 10 * x) := by
    have : x ^ 7 = x * x ^ 6 := by ring
    rw [this]
    have h6 : x * x ^ 6 ≤ x * (1000 * ε ^ 6) := mul_le_mul_of_nonneg_left hx5 hx0.le
    have h7 : ε ^ 6 ≤ ε ^ 5 / 8 := by
      have : ε ^ 6 = ε ^ 5 * ε := by ring
      rw [this]; nlinarith [pow_pos hε0 5]
    have h8 : x * (1000 * ε ^ 6) ≤ x * (1000 * (ε ^ 5 / 8)) :=
      mul_le_mul_of_nonneg_left (by linarith) hx0.le
    nlinarith [pow_pos hε0 5]
  have h3 : ε ^ 5 * (9 / 10 * x) ≤ ε ^ 5 * sx := mul_le_mul_of_nonneg_left hsx0 (by positivity)
  linarith

lemma exp_bound5 {y : ℝ} (hy : |y| ≤ 1) :
    |Real.exp y - (1 + y + y ^ 2 / 2 + y ^ 3 / 6 + y ^ 4 / 24)| ≤ |y| ^ 5 := by
  have := Real.exp_bound hy (n := 5) (by norm_num)
  have e : ∑ m ∈ Finset.range 5, y ^ m / (m.factorial : ℝ)
      = 1 + y + y ^ 2 / 2 + y ^ 3 / 6 + y ^ 4 / 24 := by
    simp [Finset.sum_range_succ, Nat.factorial]
  rw [e] at this
  refine this.trans ?_
  have : ((Nat.succ 5 : ℕ) : ℝ) / ((Nat.factorial 5 : ℕ) * (5 : ℕ)) ≤ 1 := by
    norm_num [Nat.factorial]
  exact mul_le_of_le_one_right (by positivity) (by exact_mod_cast this)

lemma abs_mul_le2 {x y X Y : ℝ} (hx : |x| ≤ X) (hy : |y| ≤ Y) : |x * y| ≤ X * Y := by
  rw [abs_mul]; exact mul_le_mul hx hy (abs_nonneg _) ((abs_nonneg _).trans hx)

/-- Truncating a product of two quartics in `ε` (`|ε| ≤ 1/8`). -/
lemma mul_trunc (u0 u1 u2 u3 u4 v0 v1 v2 v3 v4 U1 U2 U3 U4 V1 V2 V3 V4 ε : ℝ)
    (h1 : |u1| ≤ U1) (h2 : |u2| ≤ U2) (h3 : |u3| ≤ U3) (h4 : |u4| ≤ U4)
    (k1 : |v1| ≤ V1) (k2 : |v2| ≤ V2) (k3 : |v3| ≤ V3) (k4 : |v4| ≤ V4)
    (hε0 : 0 ≤ ε) (hε : ε ≤ 1 / 8) :
    |(u0 + u1 * ε + u2 * ε ^ 2 + u3 * ε ^ 3 + u4 * ε ^ 4)
        * (v0 + v1 * ε + v2 * ε ^ 2 + v3 * ε ^ 3 + v4 * ε ^ 4)
      - (u0 * v0 + (u0 * v1 + u1 * v0) * ε + (u0 * v2 + u1 * v1 + u2 * v0) * ε ^ 2
        + (u0 * v3 + u1 * v2 + u2 * v1 + u3 * v0) * ε ^ 3
        + (u0 * v4 + u1 * v3 + u2 * v2 + u3 * v1 + u4 * v0) * ε ^ 4)|
      ≤ ((U1 * V4 + U2 * V3 + U3 * V2 + U4 * V1) + (U2 * V4 + U3 * V3 + U4 * V2) / 8
        + (U3 * V4 + U4 * V3) / 64 + U4 * V4 / 512) * ε ^ 5 := by
  set A := u1 * v4 + u2 * v3 + u3 * v2 + u4 * v1
  set B := u2 * v4 + u3 * v3 + u4 * v2
  set C := u3 * v4 + u4 * v3
  set D := u4 * v4
  have e : (u0 + u1 * ε + u2 * ε ^ 2 + u3 * ε ^ 3 + u4 * ε ^ 4)
        * (v0 + v1 * ε + v2 * ε ^ 2 + v3 * ε ^ 3 + v4 * ε ^ 4)
      - (u0 * v0 + (u0 * v1 + u1 * v0) * ε + (u0 * v2 + u1 * v1 + u2 * v0) * ε ^ 2
        + (u0 * v3 + u1 * v2 + u2 * v1 + u3 * v0) * ε ^ 3
        + (u0 * v4 + u1 * v3 + u2 * v2 + u3 * v1 + u4 * v0) * ε ^ 4)
      = ε ^ 5 * (A + B * ε + C * ε ^ 2 + D * ε ^ 3) := by ring
  rw [e, abs_mul, abs_of_nonneg (by positivity), mul_comm]
  apply mul_le_mul_of_nonneg_right _ (by positivity)
  have hA : |A| ≤ U1 * V4 + U2 * V3 + U3 * V2 + U4 * V1 := by
    have := abs_add_le (u1 * v4 + u2 * v3 + u3 * v2) (u4 * v1)
    have := abs_add_le (u1 * v4 + u2 * v3) (u3 * v2)
    have := abs_add_le (u1 * v4) (u2 * v3)
    have := abs_mul_le2 h1 k4; have := abs_mul_le2 h2 k3
    have := abs_mul_le2 h3 k2; have := abs_mul_le2 h4 k1
    linarith
  have hB : |B| ≤ U2 * V4 + U3 * V3 + U4 * V2 := by
    have := abs_add_le (u2 * v4 + u3 * v3) (u4 * v2)
    have := abs_add_le (u2 * v4) (u3 * v3)
    have := abs_mul_le2 h2 k4; have := abs_mul_le2 h3 k3; have := abs_mul_le2 h4 k2
    linarith
  have hC : |C| ≤ U3 * V4 + U4 * V3 := by
    have := abs_add_le (u3 * v4) (u4 * v3)
    have := abs_mul_le2 h3 k4; have := abs_mul_le2 h4 k3
    linarith
  have hD : |D| ≤ U4 * V4 := abs_mul_le2 h4 k4
  have hε' : |ε| ≤ 1 / 8 := by rw [abs_of_nonneg hε0]; exact hε
  have hε2 : |ε ^ 2| ≤ 1 / 64 := by
    rw [abs_pow]; calc |ε| ^ 2 ≤ (1 / 8) ^ 2 := pow_le_pow_left₀ (abs_nonneg _) hε' 2
      _ = 1 / 64 := by norm_num
  have hε3 : |ε ^ 3| ≤ 1 / 512 := by
    rw [abs_pow]; calc |ε| ^ 3 ≤ (1 / 8) ^ 3 := pow_le_pow_left₀ (abs_nonneg _) hε' 3
      _ = 1 / 512 := by norm_num
  have r1 := abs_mul_le2 hB hε'
  have r2 := abs_mul_le2 hC hε2
  have r3 := abs_mul_le2 hD hε3
  have := abs_add_le (A + B * ε + C * ε ^ 2) (D * ε ^ 3)
  have := abs_add_le (A + B * ε) (C * ε ^ 2)
  have := abs_add_le A (B * ε)
  linarith

end Z06h

namespace Z06h
open Real

lemma pow_small {ε : ℝ} (h0 : 0 ≤ ε) (h1 : ε ≤ 1 / 8) :
    ε ^ 2 ≤ 1 / 64 ∧ ε ^ 3 ≤ 1 / 512 ∧ ε ^ 4 ≤ 1 / 4096 ∧ ε ^ 5 ≤ 1 / 32768 := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · calc ε ^ 2 ≤ (1 / 8) ^ 2 := pow_le_pow_left₀ h0 h1 2
      _ = 1 / 64 := by norm_num
  · calc ε ^ 3 ≤ (1 / 8) ^ 3 := pow_le_pow_left₀ h0 h1 3
      _ = 1 / 512 := by norm_num
  · calc ε ^ 4 ≤ (1 / 8) ^ 4 := pow_le_pow_left₀ h0 h1 4
      _ = 1 / 4096 := by norm_num
  · calc ε ^ 5 ≤ (1 / 8) ^ 5 := pow_le_pow_left₀ h0 h1 5
      _ = 1 / 32768 := by norm_num

/-- the algebra of Step 0(v) (order 5). -/
lemma taylor_alg4 (a L Z Y f g ε s F e : ℝ) (ha0 : 0 < a) (ha1 : a < 1) (hL0 : 0 < L) (hL1 : L < 2)
    (hZ0 : 0 < Z) (hZ1 : Z < 2) (hY0 : 0 < Y) (hY1 : Y < 2) (hf : |f| ≤ 2) (hg : |g| ≤ 4)
    (hε0 : 0 < ε) (hε1 : ε ≤ 1 / 8)
    (hs : |s - (1 + Z * ε ^ 2 + Y * ε ^ 4)| ≤ ε ^ 5)
    (hF : |F - (1 - a * ε + (a - L) * ε ^ 2 + f * ε ^ 3 + g * ε ^ 4)| ≤ 3 * ε ^ 5)
    (he : |e - (1 + a * ε + a ^ 2 * ε ^ 2 / 2 + a ^ 3 * ε ^ 3 / 6 + a ^ 4 * ε ^ 4 / 24)| ≤ ε ^ 5) :
    |s * F * e - (1 + (Z + (a - L) - a ^ 2 / 2) * ε ^ 2 + (f + a * (a - L) - a ^ 3 / 3) * ε ^ 3
      + (g + f * a + (a - L) * a ^ 2 / 2 - a ^ 4 / 8 + Z * ((a - L) - a ^ 2 / 2) + Y) * ε ^ 4)|
      ≤ 100 * ε ^ 5 := by
  obtain ⟨b, hbd⟩ : ∃ b : ℝ, b = a - L := ⟨_, rfl⟩
  rw [← hbd] at hF ⊢
  have hb : |b| ≤ 2 := by rw [hbd, abs_le]; constructor <;> linarith
  have ha : |a| ≤ 1 := by rw [abs_le]; constructor <;> linarith
  have hna : |-a| ≤ 1 := by rw [abs_neg]; exact ha
  have hε5 : 0 < ε ^ 5 := by positivity
  have ha2 : |a ^ 2 / 2| ≤ 1 / 2 := by
    rw [abs_div, abs_pow, abs_two]
    have := pow_le_one₀ (abs_nonneg a) ha (n := 2); linarith
  have ha3 : |a ^ 3 / 6| ≤ 1 / 6 := by
    rw [abs_div, abs_pow, abs_of_pos (by norm_num : (0:ℝ) < 6)]
    have := pow_le_one₀ (abs_nonneg a) ha (n := 3); linarith
  have ha4 : |a ^ 4 / 24| ≤ 1 / 24 := by
    rw [abs_div, abs_pow, abs_of_pos (by norm_num : (0:ℝ) < 24)]
    have := pow_le_one₀ (abs_nonneg a) ha (n := 4); linarith
  have T1 := mul_trunc 1 (-a) b f g 1 a (a ^ 2 / 2) (a ^ 3 / 6) (a ^ 4 / 24) 1 2 2 4 1 (1 / 2)
    (1 / 6) (1 / 24) ε hna hb hf hg ha ha2 ha3 ha4 hε0.le hε1
  obtain ⟨w2, hw2⟩ : ∃ w : ℝ, w = b - a ^ 2 / 2 := ⟨_, rfl⟩
  obtain ⟨w3, hw3⟩ : ∃ w : ℝ, w = f + a * b - a ^ 3 / 3 := ⟨_, rfl⟩
  obtain ⟨w4, hw4⟩ : ∃ w : ℝ, w = g + f * a + b * a ^ 2 / 2 - a ^ 4 / 8 := ⟨_, rfl⟩
  have E1 : 1 * 1 + (1 * a + -a * 1) * ε + (1 * (a ^ 2 / 2) + -a * a + b * 1) * ε ^ 2
      + (1 * (a ^ 3 / 6) + -a * (a ^ 2 / 2) + b * a + f * 1) * ε ^ 3
      + (1 * (a ^ 4 / 24) + -a * (a ^ 3 / 6) + b * (a ^ 2 / 2) + f * a + g * 1) * ε ^ 4
      = 1 + 0 * ε + w2 * ε ^ 2 + w3 * ε ^ 3 + w4 * ε ^ 4 := by rw [hw2, hw3, hw4]; ring
  rw [E1] at T1
  have hM1 : (1 * (1 / 24) + 2 * (1 / 6) + 2 * (1 / 2) + 4 * 1) + (2 * (1 / 24) + 2 * (1 / 6)
      + 4 * (1 / 2)) / 8 + (2 * (1 / 24) + 4 * (1 / 6)) / 64 + 4 * (1 / 24) / 512 ≤ (6:ℝ) := by
    norm_num
  have T1' : |(1 - a * ε + b * ε ^ 2 + f * ε ^ 3 + g * ε ^ 4)
      * (1 + a * ε + a ^ 2 * ε ^ 2 / 2 + a ^ 3 * ε ^ 3 / 6 + a ^ 4 * ε ^ 4 / 24)
      - (1 + 0 * ε + w2 * ε ^ 2 + w3 * ε ^ 3 + w4 * ε ^ 4)| ≤ 6 * ε ^ 5 := by
    have e : (1 - a * ε + b * ε ^ 2 + f * ε ^ 3 + g * ε ^ 4)
        * (1 + a * ε + a ^ 2 * ε ^ 2 / 2 + a ^ 3 * ε ^ 3 / 6 + a ^ 4 * ε ^ 4 / 24)
        = (1 + -a * ε + b * ε ^ 2 + f * ε ^ 3 + g * ε ^ 4)
          * (1 + a * ε + a ^ 2 / 2 * ε ^ 2 + a ^ 3 / 6 * ε ^ 3 + a ^ 4 / 24 * ε ^ 4) := by ring
    rw [e]
    exact T1.trans (mul_le_mul_of_nonneg_right hM1 hε5.le)
  -- bounds on the `w`'s
  have qa := abs_le.mp ha
  have qb := abs_le.mp hb
  have qf := abs_le.mp hf
  have qg := abs_le.mp hg
  have pa2 : 0 ≤ a ^ 2 ∧ a ^ 2 ≤ 1 := ⟨by positivity, pow_le_one₀ ha0.le ha1.le⟩
  have pa3 : 0 ≤ a ^ 3 ∧ a ^ 3 ≤ 1 := ⟨by positivity, pow_le_one₀ ha0.le ha1.le⟩
  have pa4 : 0 ≤ a ^ 4 ∧ a ^ 4 ≤ 1 := ⟨by positivity, pow_le_one₀ ha0.le ha1.le⟩
  have pab := abs_le.mp (abs_mul_le2 ha hb)
  have pfa := abs_le.mp (abs_mul_le2 hf ha)
  have pba2 : |b * a ^ 2 / 2| ≤ 2 * (1 / 2) := by
    have := abs_mul_le2 hb ha2
    rwa [show b * (a ^ 2 / 2) = b * a ^ 2 / 2 by ring] at this
  have qba2 := abs_le.mp pba2
  have hW2 : |w2| ≤ 2.5 := by rw [hw2, abs_le]; constructor <;> linarith
  have hW3 : |w3| ≤ 5 := by rw [hw3, abs_le]; constructor <;> linarith
  have hW4 : |w4| ≤ 8 := by rw [hw4, abs_le]; constructor <;> linarith
  have h00 : |(0:ℝ)| ≤ 0 := by simp
  have hZ : |Z| ≤ 2 := by rw [abs_le]; constructor <;> linarith
  have hY : |Y| ≤ 2 := by rw [abs_le]; constructor <;> linarith
  have T2 := mul_trunc 1 0 Z 0 Y 1 0 w2 w3 w4 0 2 0 2 0 2.5 5 8 ε h00 hZ h00 hY h00 hW2 hW3 hW4
    hε0.le hε1
  have hM2 : ((0:ℝ) * 8 + 2 * 5 + 0 * 2.5 + 2 * 0) + (2 * 8 + 0 * 5 + 2 * 2.5) / 8
      + (0 * 8 + 2 * 5) / 64 + 2 * 8 / 512 ≤ 13 := by norm_num
  set s0 := 1 + Z * ε ^ 2 + Y * ε ^ 4 with hs0d
  set F0 := 1 - a * ε + b * ε ^ 2 + f * ε ^ 3 + g * ε ^ 4 with hF0d
  set e0 := 1 + a * ε + a ^ 2 * ε ^ 2 / 2 + a ^ 3 * ε ^ 3 / 6 + a ^ 4 * ε ^ 4 / 24 with he0d
  set G4 := 1 + 0 * ε + w2 * ε ^ 2 + w3 * ε ^ 3 + w4 * ε ^ 4 with hG4
  set T := 1 + (Z + b - a ^ 2 / 2) * ε ^ 2 + (f + a * b - a ^ 3 / 3) * ε ^ 3
      + (g + f * a + b * a ^ 2 / 2 - a ^ 4 / 8 + Z * (b - a ^ 2 / 2) + Y) * ε ^ 4 with hT
  have T2' : |s0 * G4 - T| ≤ 13 * ε ^ 5 := by
    have e1 : s0 * G4 = (1 + 0 * ε + Z * ε ^ 2 + 0 * ε ^ 3 + Y * ε ^ 4)
        * (1 + 0 * ε + w2 * ε ^ 2 + w3 * ε ^ 3 + w4 * ε ^ 4) := by rw [hs0d, hG4]; ring
    have e2 : T = 1 * 1 + (1 * 0 + 0 * 1) * ε + (1 * w2 + 0 * 0 + Z * 1) * ε ^ 2
        + (1 * w3 + 0 * w2 + Z * 0 + 0 * 1) * ε ^ 3
        + (1 * w4 + 0 * w3 + Z * w2 + 0 * 0 + Y * 1) * ε ^ 4 := by
      rw [hT, hw2, hw3, hw4]; ring
    rw [e1, e2]
    exact T2.trans (mul_le_mul_of_nonneg_right hM2 hε5.le)
  have hε' : ε ≤ 1 / 8 := hε1
  obtain ⟨hε2, hε3, hε4, hε5'⟩ := pow_small hε0.le hε1
  have hε2' : 0 ≤ ε ^ 2 := by positivity
  have hε3' : 0 ≤ ε ^ 3 := by positivity
  have hε4' : 0 ≤ ε ^ 4 := by positivity
  have m1 := abs_le.mp (abs_mul_le2 hZ (abs_of_nonneg hε2' ▸ hε2))
  have m1' := abs_le.mp (abs_mul_le2 hY (abs_of_nonneg hε4' ▸ hε4))
  have m2 := abs_le.mp (abs_mul_le2 ha (abs_of_nonneg hε0.le ▸ hε'))
  have m3 := abs_le.mp (abs_mul_le2 hb (abs_of_nonneg hε2' ▸ hε2))
  have m4 := abs_le.mp (abs_mul_le2 hf (abs_of_nonneg hε3' ▸ hε3))
  have m5 := abs_le.mp (abs_mul_le2 hg (abs_of_nonneg hε4' ▸ hε4))
  have hs0 : |s0| ≤ 2 := by rw [hs0d, abs_le]; constructor <;> linarith
  have hF0 : |F0| ≤ 2 := by rw [hF0d, abs_le]; constructor <;> linarith
  have he0 : |e0| ≤ 2 := by
    have n2 : 0 ≤ a ^ 2 * ε ^ 2 / 2 ∧ a ^ 2 * ε ^ 2 / 2 ≤ 1 / 128 := by
      constructor
      · positivity
      · have := mul_le_mul pa2.2 hε2 hε2' (by norm_num); linarith
    have n3 : 0 ≤ a ^ 3 * ε ^ 3 / 6 ∧ a ^ 3 * ε ^ 3 / 6 ≤ 1 / 512 := by
      constructor
      · positivity
      · have := mul_le_mul pa3.2 hε3 hε3' (by norm_num); linarith
    have n4 : 0 ≤ a ^ 4 * ε ^ 4 / 24 ∧ a ^ 4 * ε ^ 4 / 24 ≤ 1 / 4096 := by
      constructor
      · positivity
      · have := mul_le_mul pa4.2 hε4 hε4' (by norm_num); linarith
    rw [he0d, abs_le]; constructor <;> linarith
  set rs := s - s0 with hrs
  set rF := F - F0 with hrF
  set re := e - e0 with hre
  have hFb : |F| ≤ 3 := by
    have : F = F0 + rF := by rw [hrF]; ring
    rw [this]; calc |F0 + rF| ≤ |F0| + |rF| := abs_add_le _ _
      _ ≤ 3 := by linarith
  have heb : |e| ≤ 3 := by
    have : e = e0 + re := by rw [hre]; ring
    rw [this]; calc |e0 + re| ≤ |e0| + |re| := abs_add_le _ _
      _ ≤ 3 := by linarith
  have hdecomp : s * F * e - T = s0 * (F0 * e0 - G4) + (s0 * G4 - T)
      + rs * F * e + s0 * rF * e + s0 * F0 * re := by
    rw [hrs, hrF, hre]; ring
  rw [hdecomp]
  have t0 : |s0 * (F0 * e0 - G4)| ≤ 2 * (6 * ε ^ 5) := abs_mul_le2 hs0 T1'
  have t2 : |rs * F * e| ≤ ε ^ 5 * 3 * 3 := abs_mul_le2 (abs_mul_le2 hs hFb) heb
  have t3 : |s0 * rF * e| ≤ 2 * (3 * ε ^ 5) * 3 := abs_mul_le2 (abs_mul_le2 hs0 hF) heb
  have t4 : |s0 * F0 * re| ≤ 2 * 2 * ε ^ 5 := abs_mul_le2 (abs_mul_le2 hs0 hF0) he
  have u1 := abs_add_le (s0 * (F0 * e0 - G4) + (s0 * G4 - T) + rs * F * e + s0 * rF * e)
    (s0 * F0 * re)
  have u2 := abs_add_le (s0 * (F0 * e0 - G4) + (s0 * G4 - T) + rs * F * e) (s0 * rF * e)
  have u3 := abs_add_le (s0 * (F0 * e0 - G4) + (s0 * G4 - T)) (rs * F * e)
  have u4 := abs_add_le (s0 * (F0 * e0 - G4)) (s0 * G4 - T)
  linarith

end Z06h

namespace Z06h
open Real

lemma f4_bound (θ : ℝ) (h0 : 0 < θ) (h1 : θ < π) :
    |(1 / 2) * D00.Li2 (Y02.zz θ) ^ 2 - 2 * Z06.S22 (Y02.zz θ)
      - (2 * Y02.ww θ * Y02.S12 θ + (1 + Y02.U θ) * (D00.Li2 (Y02.zz θ) - 1))| ≤ 4 := by
  obtain ⟨z0, z1⟩ := Y02.zz_bounds θ h0 h1
  have hf4 := f4_hasSum θ h0 h1
  have hB := Y02.sumB z0 z1
  have hUb := Y02.U_bounds z0 z1
  have hLb := Y02.Li2_pos_lt z0 z1
  have hz : 0 < Y02.zz θ := z0
  have h0' := hasSum_le (fun k => mul_nonneg (Y02.hk_nonneg k) (by positivity)) hasSum_zero hf4
  have h1' := hasSum_le (fun k => mul_le_mul_of_nonneg_right (Y02.hk_le k)
    (by positivity)) hf4 (hB.mul_left 2)
  have hU : Y02.U θ = (1 - Real.sin (θ / 2) ^ 2) / Real.sin (θ / 2) ^ 2
      * Real.log (1 - Real.sin (θ / 2) ^ 2) := rfl
  rw [← hU] at hUb
  have hpi := Real.pi_lt_d2
  have : π ^ 2 / 6 < 1.66 := by nlinarith [Real.pi_pos]
  rw [abs_le]; constructor <;> nlinarith

lemma V_bound (θ : ℝ) (h0 : 0 < θ) (h1 : θ < π) :
    |2 * Y02.ww θ * Y02.S12 θ + (1 + Y02.U θ) * (D00.Li2 (Y02.zz θ) - 1)| ≤ 2 := by
  obtain ⟨z0, z1⟩ := Y02.zz_bounds θ h0 h1
  have hf3 := Y02.f3_hasSum θ h0 h1
  have hA := Y02.sumA z0 z1
  have hUb := Y02.U_bounds z0 z1
  have hz : 0 < Y02.zz θ := z0
  have hV0 := hasSum_le (fun k => mul_nonneg (Y02.hk_nonneg k) (by positivity)) hasSum_zero hf3
  have hV1 := hasSum_le (fun k => mul_le_mul_of_nonneg_right (Y02.hk_le k)
    (by positivity)) hf3 (hA.mul_left 2)
  have hU : Y02.U θ = (1 - Real.sin (θ / 2) ^ 2) / Real.sin (θ / 2) ^ 2
      * Real.log (1 - Real.sin (θ / 2) ^ 2) := rfl
  rw [← hU] at hUb hV1
  rw [abs_le]; constructor <;> linarith [hUb.2]

lemma pi_sq_lt : π ^ 2 / 6 < 1.66 := by
  have hpi := Real.pi_lt_d2
  nlinarith [Real.pi_pos]

lemma taylor_core4 (θ ε : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hε0 : 0 < ε) (hε1 : ε ≤ 1 / 8) :
    |Y02.R θ ε * Real.exp (-(Y02.a1 θ * ε)) / (2 * Y02.zz θ)
      - (1 + Y02.a2 θ * ε ^ 2 + Y02.a3 θ * ε ^ 3 + (Z06.a4 θ + Y02.a2 θ ^ 2 / 2) * ε ^ 4)|
      ≤ 100 * ε ^ 5 := by
  obtain ⟨z0, z1⟩ := Y02.zz_bounds θ h0 h1
  have hUb : -1 < Y02.U θ ∧ Y02.U θ < 0 := Y02.U_bounds z0 z1
  have hLb : 0 < D00.Li2 (Y02.zz θ) ∧ D00.Li2 (Y02.zz θ) < π ^ 2 / 6 := Y02.Li2_pos_lt z0 z1
  have hpi := Real.pi_gt_three
  have hpi' := Real.pi_lt_d2
  obtain ⟨r, hr, hLam⟩ := Lam_expand4 θ ε h0 h1 hε0 hε1
  have hVb := V_bound θ h0 h1
  have hWb := f4_bound θ h0 h1
  set V := 2 * Y02.ww θ * Y02.S12 θ + (1 + Y02.U θ) * (D00.Li2 (Y02.zz θ) - 1) with hV
  set W := (1 / 2) * D00.Li2 (Y02.zz θ) ^ 2 - 2 * Z06.S22 (Y02.zz θ) - V with hW
  have hz0' : 0 < Y02.zz θ := z0
  have hz0 : Y02.zz θ ≠ 0 := hz0'.ne'
  set F0 := 1 - (1 + Y02.U θ) * ε + ((1 + Y02.U θ) - D00.Li2 (Y02.zz θ)) * ε ^ 2 + V * ε ^ 3
    + W * ε ^ 4 with hF0
  have e1 : Y02.R θ ε * Real.exp (-(Y02.a1 θ * ε)) / (2 * Y02.zz θ) =
      (π * ε / Real.sin (π * ε)) * (F0 + r) * Real.exp ((1 + Y02.U θ) * ε) := by
    rw [Y02.R, hLam, Y02.a1]
    field_simp
    ring_nf
  have hpoly : 1 + Y02.a2 θ * ε ^ 2 + Y02.a3 θ * ε ^ 3 + (Z06.a4 θ + Y02.a2 θ ^ 2 / 2) * ε ^ 4
      = 1 + (π ^ 2 / 6 + ((1 + Y02.U θ) - D00.Li2 (Y02.zz θ)) - (1 + Y02.U θ) ^ 2 / 2) * ε ^ 2
        + (V + (1 + Y02.U θ) * ((1 + Y02.U θ) - D00.Li2 (Y02.zz θ)) - (1 + Y02.U θ) ^ 3 / 3) * ε ^ 3
        + (W + V * (1 + Y02.U θ) + ((1 + Y02.U θ) - D00.Li2 (Y02.zz θ)) * (1 + Y02.U θ) ^ 2 / 2
          - (1 + Y02.U θ) ^ 4 / 8
          + π ^ 2 / 6 * (((1 + Y02.U θ) - D00.Li2 (Y02.zz θ)) - (1 + Y02.U θ) ^ 2 / 2)
          + 7 * π ^ 4 / 360) * ε ^ 4 := by
    rw [hW, hV, Z06.a4, Y02.a2, Y02.a3]
    unfold Y02.U
    ring
  rw [e1, hpoly]
  have hY1 : 7 * π ^ 4 / 360 < 2 := by
    have h4 := pow_lt_pow_left₀ hpi' Real.pi_pos.le (four_ne_zero)
    have : π ^ 4 < 99 := by linarith [h4, show (3.15:ℝ) ^ 4 ≤ 99 by norm_num]
    linarith
  refine taylor_alg4 (1 + Y02.U θ) (D00.Li2 (Y02.zz θ)) (π ^ 2 / 6) (7 * π ^ 4 / 360) V W ε _ _ _
    (by linarith) (by linarith) hLb.1 (by linarith [pi_sq_lt]) (by positivity)
    (by linarith [pi_sq_lt]) (by positivity) hY1 hVb hWb hε0 hε1 (sinc_bound6 hε0 hε1) ?_ ?_
  · rw [show F0 + r - F0 = r by ring]
    simpa [hF0] using hr
  · have hp0 : 0 < 1 + Y02.U θ := by linarith
    have hp1 : 1 + Y02.U θ < 1 := by linarith
    have hpe : 0 < (1 + Y02.U θ) * ε := mul_pos hp0 hε0
    have hy : |(1 + Y02.U θ) * ε| ≤ 1 := by
      rw [abs_of_pos hpe]
      have := mul_le_mul hp1.le hε1 hε0.le (by norm_num)
      linarith
    have := exp_bound5 hy
    rw [show 1 + (1 + Y02.U θ) * ε + ((1 + Y02.U θ) * ε) ^ 2 / 2 + ((1 + Y02.U θ) * ε) ^ 3 / 6
        + ((1 + Y02.U θ) * ε) ^ 4 / 24 =
      1 + (1 + Y02.U θ) * ε + (1 + Y02.U θ) ^ 2 * ε ^ 2 / 2 + (1 + Y02.U θ) ^ 3 * ε ^ 3 / 6
        + (1 + Y02.U θ) ^ 4 * ε ^ 4 / 24 by ring] at this
    refine this.trans ?_
    rw [abs_of_pos hpe, mul_pow]
    have : (1 + Y02.U θ) ^ 5 ≤ 1 := pow_le_one₀ hp0.le hp1.le
    exact mul_le_of_le_one_left (by positivity) this

lemma p6 : ∀ θ ε : ℝ, 0 < θ → θ < π → 0 < ε → ε ≤ 1 / 8 →
    |Y02.R θ ε * Real.exp (-(Y02.a1 θ * ε))
        - 2 * Y02.zz θ * (1 + Y02.a2 θ * ε ^ 2 + Y02.a3 θ * ε ^ 3
            + (Z06.a4 θ + Y02.a2 θ ^ 2 / 2) * ε ^ 4)|
      ≤ 2 * Y02.zz θ * (3.1e5 * ε ^ 5) := by
  intro θ ε h0 h1 hε0 hε1
  have z0 : 0 < Y02.zz θ := (Y02.zz_bounds θ h0 h1).1
  have H := taylor_core4 θ ε h0 h1 hε0 hε1
  have e : Y02.R θ ε * Real.exp (-(Y02.a1 θ * ε))
      - 2 * Y02.zz θ * (1 + Y02.a2 θ * ε ^ 2 + Y02.a3 θ * ε ^ 3
            + (Z06.a4 θ + Y02.a2 θ ^ 2 / 2) * ε ^ 4) =
      2 * Y02.zz θ * (Y02.R θ ε * Real.exp (-(Y02.a1 θ * ε)) / (2 * Y02.zz θ) -
        (1 + Y02.a2 θ * ε ^ 2 + Y02.a3 θ * ε ^ 3 + (Z06.a4 θ + Y02.a2 θ ^ 2 / 2) * ε ^ 4)) := by
    field_simp
  rw [e, abs_mul, abs_of_pos (by positivity)]
  apply mul_le_mul_of_nonneg_left _ (by positivity)
  have : 0 ≤ ε ^ 5 := by positivity
  linarith

end Z06h

namespace Z06h
open Real

lemma G1 {n : ℝ} (hn : 9 ≤ n) :
    0 ≤ (1 / (n - 1)) ^ 2 - (1 / n) ^ 2 - 2 * (1 / n) ^ 3 - 3 * (1 / n) ^ 4 ∧
    (1 / (n - 1)) ^ 2 - (1 / n) ^ 2 - 2 * (1 / n) ^ 3 - 3 * (1 / n) ^ 4 ≤ 5 * (1 / n) ^ 5 := by
  have hn0 : 0 < n := by linarith
  have hn1 : 0 < n - 1 := by linarith
  have e : (1 / (n - 1)) ^ 2 - (1 / n) ^ 2 - 2 * (1 / n) ^ 3 - 3 * (1 / n) ^ 4
      = (4 * n - 3) / (n ^ 4 * (n - 1) ^ 2) := by field_simp; ring
  rw [e]
  constructor
  · apply div_nonneg (by linarith) (by positivity)
  · rw [div_le_iff₀ (by positivity)]
    have e2 : 5 * (1 / n) ^ 5 * (n ^ 4 * (n - 1) ^ 2) = 5 * (n - 1) ^ 2 / n := by field_simp
    rw [e2, le_div_iff₀ hn0]
    nlinarith

lemma G2 {n : ℝ} (hn : 9 ≤ n) :
    0 ≤ (1 / (n - 1)) ^ 3 - (1 / n) ^ 3 - 3 * (1 / n) ^ 4 ∧
    (1 / (n - 1)) ^ 3 - (1 / n) ^ 3 - 3 * (1 / n) ^ 4 ≤ 8 * (1 / n) ^ 5 := by
  have hn0 : 0 < n := by linarith
  have hn1 : 0 < n - 1 := by linarith
  have e : (1 / (n - 1)) ^ 3 - (1 / n) ^ 3 - 3 * (1 / n) ^ 4
      = (6 * n ^ 2 - 8 * n + 3) / (n ^ 4 * (n - 1) ^ 3) := by field_simp; ring
  rw [e]
  constructor
  · apply div_nonneg (by nlinarith) (by positivity)
  · rw [div_le_iff₀ (by positivity)]
    have e2 : 8 * (1 / n) ^ 5 * (n ^ 4 * (n - 1) ^ 3) = 8 * (n - 1) ^ 3 / n := by field_simp
    rw [e2, le_div_iff₀ hn0]
    have : 0 ≤ n ^ 2 * (n - 9) := mul_nonneg (by positivity) (by linarith)
    nlinarith

lemma G3 {n : ℝ} (hn : 9 ≤ n) :
    0 ≤ (1 / (n - 1)) ^ 4 - (1 / n) ^ 4 ∧ (1 / (n - 1)) ^ 4 - (1 / n) ^ 4 ≤ 7 * (1 / n) ^ 5 := by
  have hn0 : 0 < n := by linarith
  have hn1 : 0 < n - 1 := by linarith
  have e : (1 / (n - 1)) ^ 4 - (1 / n) ^ 4
      = (4 * n ^ 3 - 6 * n ^ 2 + 4 * n - 1) / (n ^ 4 * (n - 1) ^ 4) := by field_simp; ring
  rw [e]
  constructor
  · apply div_nonneg (by nlinarith) (by positivity)
  · rw [div_le_iff₀ (by positivity)]
    have e2 : 7 * (1 / n) ^ 5 * (n ^ 4 * (n - 1) ^ 4) = 7 * (n - 1) ^ 4 / n := by field_simp
    rw [e2, le_div_iff₀ hn0]
    have : 0 ≤ n ^ 3 * (n - 9) := mul_nonneg (by positivity) (by linarith)
    nlinarith

lemma G4 {n : ℝ} (hn : 9 ≤ n) :
    1 / (n - 1) ≤ 1.125 * (1 / n) ∧ (1 / (n - 1)) ^ 5 ≤ 1.81 * (1 / n) ^ 5 ∧
    (1 / (n - 1)) ^ 6 ≤ 2.03 * (1 / n) ^ 6 ∧
    (1 / (n - 1)) - 1 / n = (1 / n) * (1 / (n - 1)) := by
  have hn0 : 0 < n := by linarith
  have hn1 : 0 < n - 1 := by linarith
  have h : 1 / (n - 1) ≤ 9 / 8 * (1 / n) := by
    rw [mul_one_div, div_le_div_iff₀ hn1 (by positivity)]; linarith
  have hp : 0 ≤ 1 / (n - 1) := by positivity
  refine ⟨by linarith, ?_, ?_, ?_⟩
  · have := pow_le_pow_left₀ hp h 5
    have hx : 0 ≤ (1 / n) ^ 5 := by positivity
    calc (1 / (n - 1)) ^ 5 ≤ (9 / 8 * (1 / n)) ^ 5 := this
      _ = 59049 / 32768 * (1 / n) ^ 5 := by ring
      _ ≤ 1.81 * (1 / n) ^ 5 := by nlinarith
  · have := pow_le_pow_left₀ hp h 6
    have hx : 0 ≤ (1 / n) ^ 6 := by positivity
    calc (1 / (n - 1)) ^ 6 ≤ (9 / 8 * (1 / n)) ^ 6 := this
      _ = 531441 / 262144 * (1 / n) ^ 6 := by ring
      _ ≤ 2.03 * (1 / n) ^ 6 := by nlinarith
  · field_simp; ring

lemma G5 {ε x : ℝ} (hx : 0 ≤ x) (hxε : x ≤ ε) (hd : ε - x = x * ε) :
    0 ≤ ε ^ 5 - x ^ 5 ∧ ε ^ 5 - x ^ 5 ≤ 5 * x * ε ^ 5 := by
  have hε : 0 ≤ ε := hx.trans hxε
  have e : ε ^ 5 - x ^ 5 = (ε - x) * (ε ^ 4 + ε ^ 3 * x + ε ^ 2 * x ^ 2 + ε * x ^ 3 + x ^ 4) := by
    ring
  have h1 : 0 ≤ ε - x := by linarith
  have hs : ε ^ 4 + ε ^ 3 * x + ε ^ 2 * x ^ 2 + ε * x ^ 3 + x ^ 4 ≤ 5 * ε ^ 4 := by
    have a1 : ε ^ 3 * x ≤ ε ^ 3 * ε := mul_le_mul_of_nonneg_left hxε (by positivity)
    have a2 : x ^ 2 ≤ ε ^ 2 := pow_le_pow_left₀ hx hxε 2
    have a3 : x ^ 3 ≤ ε ^ 3 := pow_le_pow_left₀ hx hxε 3
    have a4 : x ^ 4 ≤ ε ^ 4 := pow_le_pow_left₀ hx hxε 4
    have b2 : ε ^ 2 * x ^ 2 ≤ ε ^ 2 * ε ^ 2 := mul_le_mul_of_nonneg_left a2 (by positivity)
    have b3 : ε * x ^ 3 ≤ ε * ε ^ 3 := mul_le_mul_of_nonneg_left a3 hε
    nlinarith
  have hs0 : 0 ≤ ε ^ 4 + ε ^ 3 * x + ε ^ 2 * x ^ 2 + ε * x ^ 3 + x ^ 4 := by positivity
  rw [e]
  constructor
  · exact mul_nonneg h1 hs0
  · rw [hd]
    have := mul_le_mul_of_nonneg_left hs (mul_nonneg hx hε)
    nlinarith

lemma exp_neg4 {t : ℝ} (h0 : 0 ≤ t) (h1 : t ≤ 1) :
    |Real.exp (-t) - (1 - t + t ^ 2 / 2 - t ^ 3 / 6)| ≤ t ^ 4 / 10 := by
  have hy : |-t| ≤ 1 := by rw [abs_neg, abs_of_nonneg h0]; exact h1
  have := Real.exp_bound hy (n := 4) (by norm_num)
  have e : ∑ m ∈ Finset.range 4, (-t) ^ m / (m.factorial : ℝ) = 1 - t + t ^ 2 / 2 - t ^ 3 / 6 := by
    simp [Finset.sum_range_succ, Nat.factorial]; ring
  rw [e, abs_neg, abs_of_nonneg h0] at this
  refine this.trans ?_
  norm_num [Nat.factorial]
  nlinarith [pow_nonneg h0 4]

lemma p3b_aux {ε x lam : ℝ} (hx : 0 ≤ x) (hε : 0 ≤ ε) (hl0 : 0 ≤ lam) (hlx : lam * x ≤ 1)
    (D1 : ε ^ 5 - x ^ 5 ≤ 5 * x * ε ^ 5) (E5 : ε ^ 5 ≤ 1.81 * x ^ 5) :
    lam ^ 3 * (ε ^ 5 - x ^ 5) ≤ 9.05 * (lam ^ 2 * x ^ 5) := by
  have h1 : lam ^ 3 * (ε ^ 5 - x ^ 5) ≤ lam ^ 3 * (5 * x * ε ^ 5) :=
    mul_le_mul_of_nonneg_left D1 (by positivity)
  have h2 : lam ^ 3 * (5 * x * ε ^ 5) = 5 * (lam * x) * (lam ^ 2 * ε ^ 5) := by ring
  have h0 : 0 ≤ lam ^ 2 * ε ^ 5 := by positivity
  have h3 : 5 * (lam * x) * (lam ^ 2 * ε ^ 5) ≤ 5 * 1 * (lam ^ 2 * ε ^ 5) :=
    mul_le_mul_of_nonneg_right (by linarith) h0
  have h4 : lam ^ 2 * ε ^ 5 ≤ lam ^ 2 * (1.81 * x ^ 5) :=
    mul_le_mul_of_nonneg_left E5 (by positivity)
  have e : lam ^ 2 * (1.81 * x ^ 5) = 1.81 * (lam ^ 2 * x ^ 5) := by ring
  linarith

lemma p4b_aux {x lam L : ℝ} (hx : 0 ≤ x) (hL : 0 ≤ L) (hlL : L ≤ lam) (hl1 : lam ≤ L + 1) :
    0 ≤ x ^ 5 * (lam ^ 3 - L ^ 3) ∧ x ^ 5 * (lam ^ 3 - L ^ 3) ≤ 3 * (lam ^ 2 * x ^ 5) := by
  have x5 : 0 ≤ x ^ 5 := by positivity
  have e3 : lam ^ 3 - L ^ 3 = (lam - L) * (lam ^ 2 + lam * L + L ^ 2) := by ring
  have h1 : lam * L ≤ lam ^ 2 := by nlinarith
  have h1' : L ^ 2 ≤ lam ^ 2 := pow_le_pow_left₀ hL hlL 2
  have hl0 : 0 ≤ lam := hL.trans hlL
  have h2 : 0 ≤ lam ^ 2 + lam * L + L ^ 2 := by
    have := mul_nonneg hl0 hL; positivity
  have h4 : 0 ≤ lam - L := by linarith
  have h5 : (lam - L) * (lam ^ 2 + lam * L + L ^ 2) ≤ 1 * (3 * lam ^ 2) :=
    mul_le_mul (by linarith) (by linarith) h2 (by norm_num)
  have h6 : 0 ≤ (lam - L) * (lam ^ 2 + lam * L + L ^ 2) := mul_nonneg h4 h2
  rw [e3]
  constructor
  · exact mul_nonneg x5 h6
  · have := mul_le_mul_of_nonneg_left h5 x5
    linarith [show x ^ 5 * (1 * (3 * lam ^ 2)) = 3 * (lam ^ 2 * x ^ 5) by ring]

lemma p5_aux {ε x lam η : ℝ} (hx : 0 ≤ x) (hl0 : 0 ≤ lam) (hlx2 : x * lam ^ 2 ≤ 4)
    (hq : |η| ≤ (ε * lam) ^ 4 / 10) (E6 : ε ^ 6 ≤ 2.03 * x ^ 6) :
    |ε ^ 2 * η| ≤ 0.82 * (lam ^ 2 * x ^ 5) := by
  rw [abs_mul, abs_of_nonneg (sq_nonneg _)]
  have h1 : ε ^ 2 * |η| ≤ ε ^ 2 * ((ε * lam) ^ 4 / 10) := mul_le_mul_of_nonneg_left hq (sq_nonneg _)
  have h2 : ε ^ 2 * ((ε * lam) ^ 4 / 10) = ε ^ 6 * (lam ^ 4 / 10) := by ring
  have h3 : ε ^ 6 * (lam ^ 4 / 10) ≤ 2.03 * x ^ 6 * (lam ^ 4 / 10) :=
    mul_le_mul_of_nonneg_right E6 (by positivity)
  have h4 : 2.03 * x ^ 6 * (lam ^ 4 / 10) = 0.203 * (x * lam ^ 2) * (lam ^ 2 * x ^ 5) := by ring
  have h0 : 0 ≤ lam ^ 2 * x ^ 5 := by positivity
  have h5 : 0.203 * (x * lam ^ 2) * (lam ^ 2 * x ^ 5) ≤ 0.203 * 4 * (lam ^ 2 * x ^ 5) := by
    apply mul_le_mul_of_nonneg_right _ h0
    linarith
  linarith

/-- the `P_1` bracket. -/
lemma D1_core {ε x lam L q : ℝ} (hx : 0 ≤ x) (hxε : x ≤ ε) (hL : 0 ≤ L) (hlL : L ≤ lam)
    (hl1 : lam ≤ L + 1) (hlx : lam * x ≤ 1) (hlx2 : x * lam ^ 2 ≤ 4)
    (hq : |q - (1 - ε * lam + (ε * lam) ^ 2 / 2 - (ε * lam) ^ 3 / 6)| ≤ (ε * lam) ^ 4 / 10)
    (A0 : 0 ≤ ε ^ 2 - x ^ 2 - 2 * x ^ 3 - 3 * x ^ 4)
    (A1 : ε ^ 2 - x ^ 2 - 2 * x ^ 3 - 3 * x ^ 4 ≤ 5 * x ^ 5)
    (B0 : 0 ≤ ε ^ 3 - x ^ 3 - 3 * x ^ 4) (B1 : ε ^ 3 - x ^ 3 - 3 * x ^ 4 ≤ 8 * x ^ 5)
    (C0 : 0 ≤ ε ^ 4 - x ^ 4) (C1 : ε ^ 4 - x ^ 4 ≤ 7 * x ^ 5)
    (D0 : 0 ≤ ε ^ 5 - x ^ 5) (D1 : ε ^ 5 - x ^ 5 ≤ 5 * x * ε ^ 5)
    (E5 : ε ^ 5 ≤ 1.81 * x ^ 5) (E6 : ε ^ 6 ≤ 2.03 * x ^ 6) :
    |q * ε ^ 2 - x ^ 2 - (2 - lam) * x ^ 3 - (3 - 3 * lam + lam ^ 2 / 2) * x ^ 4
        + L ^ 3 * x ^ 5 / 6| ≤ x ^ 5 * (5 + 8 * lam + 7 * lam ^ 2) := by
  have hl0 : 0 ≤ lam := hL.trans hlL
  have hε : 0 ≤ ε := hx.trans hxε
  set η := q - (1 - ε * lam + (ε * lam) ^ 2 / 2 - (ε * lam) ^ 3 / 6) with hη
  have e : q * ε ^ 2 - x ^ 2 - (2 - lam) * x ^ 3 - (3 - 3 * lam + lam ^ 2 / 2) * x ^ 4
        + L ^ 3 * x ^ 5 / 6
      = (ε ^ 2 - x ^ 2 - 2 * x ^ 3 - 3 * x ^ 4) - lam * (ε ^ 3 - x ^ 3 - 3 * x ^ 4)
        + lam ^ 2 / 2 * (ε ^ 4 - x ^ 4)
        - (lam ^ 3 * (ε ^ 5 - x ^ 5) + x ^ 5 * (lam ^ 3 - L ^ 3)) / 6 + ε ^ 2 * η := by
    rw [hη]; ring
  rw [e]
  have x5 : 0 ≤ x ^ 5 := by positivity
  have p1a : 0 ≤ lam * (ε ^ 3 - x ^ 3 - 3 * x ^ 4) := mul_nonneg hl0 B0
  have p1b : lam * (ε ^ 3 - x ^ 3 - 3 * x ^ 4) ≤ 8 * (lam * x ^ 5) := by
    have := mul_le_mul_of_nonneg_left B1 hl0
    linarith [show lam * (8 * x ^ 5) = 8 * (lam * x ^ 5) by ring]
  have p2a : 0 ≤ lam ^ 2 / 2 * (ε ^ 4 - x ^ 4) := mul_nonneg (by positivity) C0
  have p2b : lam ^ 2 / 2 * (ε ^ 4 - x ^ 4) ≤ 3.5 * (lam ^ 2 * x ^ 5) := by
    have := mul_le_mul_of_nonneg_left C1 (by positivity : (0:ℝ) ≤ lam ^ 2 / 2)
    linarith [show lam ^ 2 / 2 * (7 * x ^ 5) = 3.5 * (lam ^ 2 * x ^ 5) by ring]
  have p3a : 0 ≤ lam ^ 3 * (ε ^ 5 - x ^ 5) := mul_nonneg (by positivity) D0
  have p3b := p3b_aux hx hε hl0 hlx D1 E5
  obtain ⟨p4a, p4b⟩ := p4b_aux hx hL hlL hl1
  have p5 := abs_le.mp (p5_aux (η := η) hx hl0 hlx2 hq E6)
  have hR : x ^ 5 * (5 + 8 * lam + 7 * lam ^ 2)
      = 5 * x ^ 5 + 8 * (lam * x ^ 5) + 7 * (lam ^ 2 * x ^ 5) := by ring
  rw [hR, abs_le]
  constructor <;> linarith

/-- the `P_2` bracket. -/
lemma D2_core {ε x lam q : ℝ} (hx : 0 ≤ x) (hxε : x ≤ ε) (hl0 : 0 ≤ lam)
    (hq1 : 1 - ε * lam ≤ q) (hq2 : q ≤ 1 - ε * lam + (ε * lam) ^ 2 / 2)
    (B0 : 0 ≤ ε ^ 3 - x ^ 3 - 3 * x ^ 4) (B1 : ε ^ 3 - x ^ 3 - 3 * x ^ 4 ≤ 8 * x ^ 5)
    (C0 : 0 ≤ ε ^ 4 - x ^ 4) (C1 : ε ^ 4 - x ^ 4 ≤ 7 * x ^ 5) (E5 : ε ^ 5 ≤ 1.81 * x ^ 5) :
    |q * ε ^ 3 - x ^ 3 - (3 - lam) * x ^ 4| ≤ x ^ 5 * (8 + 7 * lam + lam ^ 2) := by
  have hε : 0 ≤ ε := hx.trans hxε
  set D := q - 1 + ε * lam with hD
  have e : q * ε ^ 3 - x ^ 3 - (3 - lam) * x ^ 4
      = (ε ^ 3 - x ^ 3 - 3 * x ^ 4) - lam * (ε ^ 4 - x ^ 4) + D * ε ^ 3 := by rw [hD]; ring
  rw [e]
  have D0 : 0 ≤ D := by linarith
  have D1 : D ≤ (ε * lam) ^ 2 / 2 := by linarith
  have p1 : 0 ≤ lam * (ε ^ 4 - x ^ 4) ∧ lam * (ε ^ 4 - x ^ 4) ≤ lam * (7 * x ^ 5) :=
    ⟨mul_nonneg hl0 C0, mul_le_mul_of_nonneg_left C1 hl0⟩
  have p2 : 0 ≤ D * ε ^ 3 := mul_nonneg D0 (by positivity)
  have p3 : D * ε ^ 3 ≤ lam ^ 2 * x ^ 5 := by
    have h1 : D * ε ^ 3 ≤ (ε * lam) ^ 2 / 2 * ε ^ 3 := mul_le_mul_of_nonneg_right D1 (by positivity)
    have h2 : (ε * lam) ^ 2 / 2 * ε ^ 3 = lam ^ 2 / 2 * ε ^ 5 := by ring
    have h3 : lam ^ 2 / 2 * ε ^ 5 ≤ lam ^ 2 / 2 * (1.81 * x ^ 5) :=
      mul_le_mul_of_nonneg_left E5 (by positivity)
    have : 0 ≤ lam ^ 2 * x ^ 5 := by positivity
    nlinarith
  have x5 : 0 ≤ x ^ 5 := by positivity
  rw [abs_le]; constructor <;> nlinarith

/-- the `P_3` bracket. -/
lemma D3_core {ε x lam q : ℝ} (hx : 0 ≤ x) (hxε : x ≤ ε) (hl0 : 0 ≤ lam)
    (hq1 : 1 - ε * lam ≤ q) (hq : q ≤ 1)
    (C0 : 0 ≤ ε ^ 4 - x ^ 4) (C1 : ε ^ 4 - x ^ 4 ≤ 7 * x ^ 5) (E5 : ε ^ 5 ≤ 1.81 * x ^ 5) :
    |q * ε ^ 4 - x ^ 4| ≤ x ^ 5 * (7 + 2 * lam) := by
  have hε : 0 ≤ ε := hx.trans hxε
  have e : q * ε ^ 4 - x ^ 4 = (ε ^ 4 - x ^ 4) - (1 - q) * ε ^ 4 := by ring
  rw [e]
  have h0 : 0 ≤ (1 - q) * ε ^ 4 := mul_nonneg (by linarith) (by positivity)
  have h1 : (1 - q) * ε ^ 4 ≤ ε * lam * ε ^ 4 := mul_le_mul_of_nonneg_right (by linarith) (by positivity)
  have h2 : ε * lam * ε ^ 4 = lam * ε ^ 5 := by ring
  have h3 : lam * ε ^ 5 ≤ lam * (1.81 * x ^ 5) := mul_le_mul_of_nonneg_left E5 hl0
  have x5 : 0 ≤ x ^ 5 := by positivity
  have : 0 ≤ lam * x ^ 5 := by positivity
  rw [abs_le]; constructor <;> nlinarith

end Z06h

namespace Z06h
open Real Polynomial

lemma Lam3 (θ : ℝ) : T13.LamInt 3 θ =
    2 * Y02.zz θ * (1 - Y02.zz θ) * (1 - 5 * Y02.zz θ + 5 * Y02.zz θ ^ 2) := by
  have h := T13.proof_2 3 (by norm_num) θ
  have hb0 := congrArg (Polynomial.eval (Real.cos θ)) (T13.P_bonnet 0)
  have hb1 := congrArg (Polynomial.eval (Real.cos θ)) (T13.P_bonnet 1)
  have hb2 := congrArg (Polynomial.eval (Real.cos θ)) (T13.P_bonnet 2)
  simp only [Polynomial.eval_smul, Polynomial.eval_sub, Polynomial.eval_mul, Polynomial.eval_X,
    smul_eq_mul, T13.P_zero, T13.P_one, Polynomial.eval_one] at hb0 hb1 hb2
  norm_num at hb0 hb1 hb2 h
  have hz : Y02.zz θ = (1 - Real.cos θ) / 2 := T13.sin_half_sq θ
  rw [h, hz]
  set c := Real.cos θ
  have e2 : (T13.P 2).eval c = (3 * c * c - 1) / 2 := by
    rw [T13.P_one, Polynomial.eval_X] at hb0; linarith
  have e3 : (T13.P 3).eval c = (5 * c * (T13.P 2).eval c - 2 * c) / 3 := by linarith
  have e4 : (T13.P 4).eval c = (7 * c * (T13.P 3).eval c - 3 * (T13.P 2).eval c) / 4 := by
    linarith
  rw [e4, e3, e2]
  ring

lemma serTerm_two (θ n : ℝ) (hn : 9 ≤ n) :
    Y02.serTerm θ n 2 = -(T13.LamInt 3 θ / ((3 * n - 4) * n ^ 3)) := by
  have hn1 : n - 1 ≠ 0 := by intro h; linarith
  have hn3 : 3 * n - 4 ≠ 0 := by intro h; linarith
  have hn0 : n ≠ 0 := by intro h; linarith
  simp only [Y02.serTerm, Y02.eps]
  norm_num
  have e : (3:ℝ) - (n - 1)⁻¹ = (3 * n - 4) / (n - 1) := by field_simp; ring
  rw [e]
  field_simp

lemma D5_bound (θ n : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hn : 9 ≤ n) :
    |Y02.serTerm θ n 0 + 2 * Y02.zz θ * (1 - Y02.zz θ) * (1 / n) ^ 2
        + 4 * Y02.zz θ * (1 - Y02.zz θ) * ((1 / n) ^ 3 + 2 * (1 / n) ^ 4)|
      ≤ 21 * Y02.zz θ * (1 / n) ^ 5 := by
  have z0 : 0 < Y02.zz θ := (Y02.zz_bounds θ h0 h1).1
  have z1 : Y02.zz θ < 1 := (Y02.zz_bounds θ h0 h1).2
  rw [Y02.serTerm_zero θ n hn]
  have hn0 : 0 < n := by linarith
  have hn2 : 0 < n - 2 := by linarith
  set z := Y02.zz θ
  have e : -(2 * z * (1 - z)) / (n * (n - 2)) + 2 * z * (1 - z) * (1 / n) ^ 2
      + 4 * z * (1 - z) * ((1 / n) ^ 3 + 2 * (1 / n) ^ 4)
      = -(16 * z * (1 - z)) * ((1 / n) ^ 5 * (n / (n - 2))) := by
    field_simp; ring
  rw [e, abs_mul, abs_neg, abs_of_pos (by nlinarith), abs_of_pos (by positivity)]
  have h1 : n / (n - 2) ≤ 9 / 7 := by rw [div_le_div_iff₀ hn2 (by norm_num)]; linarith
  have h2 : 16 * z * (1 - z) ≤ 16 * z := by nlinarith
  have hp : 0 ≤ (1 / n) ^ 5 := by positivity
  have := mul_le_mul h2 (mul_le_mul_of_nonneg_left h1 hp) (by positivity) (by positivity)
  nlinarith

lemma D6a_bound (θ n : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hn : 9 ≤ n) :
    |Y02.serTerm θ n 1 - 2 * Y02.zz θ * (1 - Y02.zz θ) * (1 - 2 * Y02.zz θ) / 2 * (1 / n) ^ 3
        - 3 / 4 * (2 * Y02.zz θ * (1 - Y02.zz θ) * (1 - 2 * Y02.zz θ)) * (1 / n) ^ 4|
      ≤ 2.7 * Y02.zz θ * (1 / n) ^ 5 := by
  have z0 : 0 < Y02.zz θ := (Y02.zz_bounds θ h0 h1).1
  have z1 : Y02.zz θ < 1 := (Y02.zz_bounds θ h0 h1).2
  rw [serTerm_one θ n hn, Lam2]
  have hn0 : 0 < n := by linarith
  have hn3 : 0 < 2 * n - 3 := by linarith
  set z := Y02.zz θ
  have e : 2 * z * (1 - z) * (1 - 2 * z) / ((2 * n - 3) * n ^ 2)
      - 2 * z * (1 - z) * (1 - 2 * z) / 2 * (1 / n) ^ 3
      - 3 / 4 * (2 * z * (1 - z) * (1 - 2 * z)) * (1 / n) ^ 4
      = 9 / 2 * (z * (1 - z) * (1 - 2 * z)) * ((1 / n) ^ 5 * (n / (2 * n - 3))) := by
    field_simp; ring
  rw [e, abs_mul, abs_mul, abs_of_pos (by norm_num : (0:ℝ) < 9 / 2), abs_of_pos (by positivity :
    (0:ℝ) < (1 / n) ^ 5 * (n / (2 * n - 3)))]
  have h1 : n / (2 * n - 3) ≤ 0.6 := by rw [div_le_iff₀ hn3]; linarith
  have h2 : |z * (1 - z) * (1 - 2 * z)| ≤ z := by
    rw [abs_le]; constructor <;> nlinarith [mul_pos z0 (by linarith : (0:ℝ) < 1 - z)]
  have hp : 0 ≤ (1 / n) ^ 5 := by positivity
  have := mul_le_mul h2 (mul_le_mul_of_nonneg_left h1 hp) (by positivity) z0.le
  nlinarith

lemma D6b_bound (θ n : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hn : 9 ≤ n) :
    |Y02.serTerm θ n 2
        + 1 / 3 * (2 * Y02.zz θ * (1 - Y02.zz θ) * (1 - 5 * Y02.zz θ + 5 * Y02.zz θ ^ 2))
          * (1 / n) ^ 4| ≤ 1.1 * Y02.zz θ * (1 / n) ^ 5 := by
  have z0 : 0 < Y02.zz θ := (Y02.zz_bounds θ h0 h1).1
  have z1 : Y02.zz θ < 1 := (Y02.zz_bounds θ h0 h1).2
  rw [serTerm_two θ n hn, Lam3]
  have hn0 : 0 < n := by linarith
  have hn3 : 0 < 3 * n - 4 := by linarith
  set z := Y02.zz θ
  have e : -(2 * z * (1 - z) * (1 - 5 * z + 5 * z ^ 2) / ((3 * n - 4) * n ^ 3))
      + 1 / 3 * (2 * z * (1 - z) * (1 - 5 * z + 5 * z ^ 2)) * (1 / n) ^ 4
      = -(8 / 3 * (z * (1 - z) * (1 - 5 * z + 5 * z ^ 2))) * ((1 / n) ^ 5 * (n / (3 * n - 4))) := by
    field_simp; ring
  rw [e, abs_mul, abs_neg, abs_mul, abs_of_pos (by norm_num : (0:ℝ) < 8 / 3), abs_of_pos
    (by positivity : (0:ℝ) < (1 / n) ^ 5 * (n / (3 * n - 4)))]
  have h1 : n / (3 * n - 4) ≤ 0.4 := by rw [div_le_iff₀ hn3]; linarith
  have h2 : |z * (1 - z) * (1 - 5 * z + 5 * z ^ 2)| ≤ z := by
    have hw : 0 < z * (1 - z) := mul_pos z0 (by linarith)
    have hw1 : z * (1 - z) ≤ z := by nlinarith
    have hq : |1 - 5 * z + 5 * z ^ 2| ≤ 1 := by
      rw [abs_le]; constructor <;> nlinarith
    rw [abs_mul, abs_of_pos hw]
    calc z * (1 - z) * |1 - 5 * z + 5 * z ^ 2| ≤ z * (1 - z) * 1 :=
          mul_le_mul_of_nonneg_left hq hw.le
      _ ≤ z := by linarith
  have hp : 0 ≤ (1 / n) ^ 5 := by positivity
  have := mul_le_mul h2 (mul_le_mul_of_nonneg_left h1 hp) (by positivity) z0.le
  nlinarith

/-- the tail `m ≥ 4`. -/
lemma tail4_bound (θ n : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hn : 9 ≤ n) {T : ℝ}
    (hT : HasSum (fun i => Y02.serTerm θ n (i + 3)) T) : |T| ≤ 2.6 * Y02.zz θ * (1 / n) ^ 5 := by
  obtain ⟨z0, z1⟩ := Y02.zz_bounds θ h0 h1
  have hz0 : 0 < Y02.zz θ := z0
  obtain ⟨hε0, hε1, hεw, hε2, hε3⟩ := Y02.eps_facts hn
  have hn0 : 0 < n := by linarith
  set z := Y02.zz θ with hzdef
  set ε := Y02.eps n with hεdef
  set w := 1 / n with hw
  have hw0 : 0 < w := by positivity
  have hw19 : w ≤ 1 / 9 := by rw [hw, div_le_div_iff₀ hn0 (by norm_num)]; linarith
  have hw1 : w < 1 := by linarith
  have hg : HasSum (fun i : ℕ => 2 * z * ε * w ^ 4 * w ^ i) (2 * z * ε * w ^ 4 * (1 - w)⁻¹) :=
    (hasSum_geometric_of_lt_one hw0.le hw1).mul_left _
  have hc : 1 - Real.cos θ = 2 * z := by
    rw [hzdef]; have := T13.sin_half_sq θ; unfold Y02.zz; linarith
  have hb : ∀ i : ℕ, |Y02.serTerm θ n (i + 3)| ≤ 2 * z * ε * w ^ 4 * w ^ i := by
    intro i
    have hL2 := (le_min_iff.mp (Y02.stmt2_aux (i + 3 + 1) (by omega) θ)).1
    rw [hc] at hL2
    have hd : 1 ≤ ((i + 3 : ℕ) : ℝ) + 1 - ε := by
      push_cast; linarith [(Nat.cast_nonneg i : (0:ℝ) ≤ i)]
    have hd0 : (0:ℝ) < ((i + 3 : ℕ) : ℝ) + 1 - ε := by linarith
    simp only [Y02.serTerm, ← hεdef]
    rw [abs_div, abs_mul, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
      abs_div, abs_of_pos hε0, abs_of_pos hd0,
      abs_of_pos (by positivity : (0:ℝ) < n ^ (i + 3 + 1))]
    have hf : ε / (((i + 3 : ℕ) : ℝ) + 1 - ε) ≤ ε := div_le_self hε0.le hd
    have hp : 1 / n ^ (i + 3 + 1) = w ^ 4 * w ^ i := by
      rw [hw, ← pow_add, one_div_pow]; ring_nf
    rw [div_eq_mul_one_div, hp]
    have := mul_le_mul hf hL2 (abs_nonneg _) hε0.le
    calc ε / (((i + 3 : ℕ) : ℝ) + 1 - ε) * |T13.LamInt (i + 3 + 1) θ| * (w ^ 4 * w ^ i)
        ≤ ε * (2 * z) * (w ^ 4 * w ^ i) := mul_le_mul_of_nonneg_right this (by positivity)
      _ = 2 * z * ε * w ^ 4 * w ^ i := by ring
  have hup := hasSum_le (fun i => (abs_le.mp (hb i)).2) hT hg
  have hlo := hasSum_le (fun i => by have := (abs_le.mp (hb i)).1; linarith : ∀ i : ℕ,
    -Y02.serTerm θ n (i + 3) ≤ 2 * z * ε * w ^ 4 * w ^ i) hT.neg hg
  have hgeo : 2 * z * ε * w ^ 4 * (1 - w)⁻¹ ≤ 2.6 * z * w ^ 5 := by
    have h1 : (1 - w)⁻¹ ≤ 9 / 8 := by
      rw [inv_le_comm₀ (by linarith) (by norm_num)]
      linarith
    have hε98 : ε ≤ 9 / 8 * w := by
      rw [hεdef, Y02.eps, hw, mul_one_div, div_le_div_iff₀ (by linarith) hn0]; linarith
    have h2 : 2 * z * ε * w ^ 4 ≤ 2 * z * (9 / 8 * w) * w ^ 4 := by
      have : 0 ≤ 2 * z * w ^ 4 := by positivity
      have := mul_le_mul_of_nonneg_left hε98 this
      linarith
    have := mul_le_mul h2 h1 (inv_nonneg.mpr (by linarith)) (by positivity)
    have hp : 0 ≤ z * w ^ 5 := by positivity
    calc _ ≤ 2 * z * (9 / 8 * w) * w ^ 4 * (9 / 8) := this
      _ = 2.53125 * (z * w ^ 5) := by ring
      _ ≤ 2.6 * z * w ^ 5 := by nlinarith
  rw [abs_le]; constructor <;> linarith

end Z06h


namespace Z06h
open Real

lemma scal5 {lam L : ℝ} (h0 : 0 ≤ lam) (h1 : lam ≤ L + 1) (hL : 2 ≤ L) :
    5 + 8 * lam + 7 * lam ^ 2 ≤ 23 * L ^ 2 ∧ 8 + 7 * lam + lam ^ 2 ≤ 9.5 * L ^ 2 ∧
    7 + 2 * lam ≤ 3.25 * L ^ 2 := by
  have h2 : lam ≤ 1.5 * L := by linarith
  have h3 : lam ^ 2 ≤ 2.25 * L ^ 2 := by nlinarith
  refine ⟨by nlinarith, by nlinarith, by nlinarith⟩

set_option maxHeartbeats 1000000 in
lemma E4_core (θ n : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hn : 9 ≤ n) :
    |Z06.E4 θ n + Y02.zz θ * Y02.a2 θ / 3 * Real.log n ^ 3 / n ^ 5|
      ≤ 1000 * (Y02.zz θ * (1 / n) ^ 5 * Real.log n ^ 2) := by
  have z0 : 0 < Y02.zz θ := (Y02.zz_bounds θ h0 h1).1
  have z1 : Y02.zz θ < 1 := (Y02.zz_bounds θ h0 h1).2
  obtain ⟨-, -, ha1, ha2, hb0, hb1⟩ := Y02.proof_9 θ h0 h1
  obtain ⟨T4, hT, hE⟩ := decomp θ n h0 h1 hn
  obtain ⟨hq0, hq1, hqe⟩ := q_facts θ n h0 h1 hn
  obtain ⟨hε0, hε1, -, -, -⟩ := Y02.eps_facts hn
  have hρ := taylor_core4 θ (Y02.eps n) h0 h1 hε0 hε1
  have hε : Y02.eps n = 1 / (n - 1) := rfl
  rw [hε] at hρ hε0 hε1
  have hT6 : HasSum (fun i => Y02.serTerm θ n (i + 3))
      (T4 - Y02.serTerm θ n 1 - Y02.serTerm θ n 2) := by
    have := (hasSum_nat_add_iff' 2).mpr hT
    convert this using 1
    simp [Finset.sum_range_succ]; ring
  have t6 := tail4_bound θ n h0 h1 hn hT6
  have t5 := D5_bound θ n h0 h1 hn
  have t6a := D6a_bound θ n h0 h1 hn
  have t6b := D6b_bound θ n h0 h1 hn
  obtain ⟨g1a, g1b⟩ := G1 hn
  obtain ⟨g2a, g2b⟩ := G2 hn
  obtain ⟨g3a, g3b⟩ := G3 hn
  obtain ⟨g4a, g4b, g4c, g4d⟩ := G4 hn
  have ha3 := a3_abs θ h0 h1
  have hb4 := b4_abs θ h0 h1
  have ha216 := a2_le θ h0 h1
  have hL2 := Y02.two_le_log hn
  have hn0 : 0 < n := by linarith
  have hLh := log_le_half hn0
  have hL1sq := logp1_sq (n := n) (by linarith)
  set L := Real.log n with hLdef
  set lam := L - Y02.a1 θ with hlam
  have hlam0 : 0 ≤ lam := by linarith
  have hlamL : L ≤ lam := by linarith
  have hlam1 : lam ≤ L + 1 := by linarith
  set ε := 1 / (n - 1) with hεdef
  set x := 1 / n with hx
  have hx0 : 0 ≤ x := by positivity
  have hxε : x ≤ ε := by
    rw [hx, hεdef]; apply one_div_le_one_div_of_le (by linarith); linarith
  obtain ⟨g5a, g5b⟩ := G5 hx0 hxε g4d
  set z := Y02.zz θ with hz
  set q := (Y02.mu θ / n) ^ (1 / (n - 1)) with hq
  have ht0 : 0 ≤ ε * lam := mul_nonneg hε0.le hlam0
  have ht1 : ε * lam ≤ 1 := by
    rw [hεdef, div_mul_eq_mul_div, one_mul, div_le_one (by linarith)]; linarith
  have hlx : lam * x ≤ 1 := by
    rw [hx, mul_one_div, div_le_one hn0]; linarith
  have hlx2 : x * lam ^ 2 ≤ 4 := by
    have h1 : lam ^ 2 ≤ (L + 1) ^ 2 := pow_le_pow_left₀ hlam0 hlam1 2
    have h2 : lam ^ 2 ≤ 4 * n := h1.trans hL1sq
    rw [hx, one_div_mul_eq_div, div_le_iff₀ hn0]; linarith
  obtain ⟨qlo, qhi⟩ := exp_neg_bounds ht0
  have q4 := exp_neg4 ht0 ht1
  rw [← hqe] at qlo qhi q4
  have b1 := D1_core hx0 hxε (by linarith) hlamL hlam1 hlx hlx2 q4 g1a g1b g2a g2b g3a g3b g5a g5b
    g4b g4c
  have b2 := D2_core hx0 hxε hlam0 qlo qhi g2a g2b g3a g3b g4b
  have b3 := D3_core hx0 hxε hlam0 qlo hq1 g3a g3b g4b
  obtain ⟨s1, s2, s3⟩ := scal5 hlam0 hlam1 hL2
  set b4 := Z06.a4 θ + Y02.a2 θ ^ 2 / 2 with hb4d
  set ρ5 := Y02.R θ ε * Real.exp (-(Y02.a1 θ * ε)) / (2 * z)
    - (1 + Y02.a2 θ * ε ^ 2 + Y02.a3 θ * ε ^ 3 + b4 * ε ^ 4) with hρ5
  set W := z * x ^ 5 with hW
  have hW0 : 0 ≤ W := by positivity
  have hWL : W ≤ W * L ^ 2 / 4 := by
    have h4 : (4:ℝ) ≤ L ^ 2 := by
      have := mul_le_mul hL2 hL2 (by norm_num) (by linarith)
      linarith [show L ^ 2 = L * L by ring]
    have := mul_le_mul_of_nonneg_left h4 hW0
    linarith
  set B1 := q * ε ^ 2 - x ^ 2 - (2 - lam) * x ^ 3 - (3 - 3 * lam + lam ^ 2 / 2) * x ^ 4
    + L ^ 3 * x ^ 5 / 6
  set B2 := q * ε ^ 3 - x ^ 3 - (3 - lam) * x ^ 4
  set B3 := q * ε ^ 4 - x ^ 4
  set D5 := Y02.serTerm θ n 0 + 2 * z * (1 - z) * x ^ 2 + 4 * z * (1 - z) * (x ^ 3 + 2 * x ^ 4)
  set D6a := Y02.serTerm θ n 1 - 2 * z * (1 - z) * (1 - 2 * z) / 2 * x ^ 3
    - 3 / 4 * (2 * z * (1 - z) * (1 - 2 * z)) * x ^ 4
  set D6b := Y02.serTerm θ n 2 + 1 / 3 * (2 * z * (1 - z) * (1 - 5 * z + 5 * z ^ 2)) * x ^ 4
  set T6 := T4 - Y02.serTerm θ n 1 - Y02.serTerm θ n 2
  have key : Z06.E4 θ n + z * Y02.a2 θ / 3 * L ^ 3 / n ^ 5
      = 2 * z * Y02.a2 θ * B1 + 2 * z * Y02.a3 θ * B2 + 2 * z * b4 * B3 + 2 * z * q * ρ5
        + D5 + D6a + D6b + T6 := by
    rw [Z06.E4, hE, Z06.Q4, Z06.d, Z06.p4, Z06.lam]
    simp only [B1, B2, B3, D5, D6a, D6b, T6, hρ5, hlam, hx, hz, hb4d]
    field_simp
    ring
  rw [key]
  have c1 : |2 * z * Y02.a2 θ * B1| ≤ 4.32 * 23 * (W * L ^ 2) := by
    rw [abs_mul, abs_of_pos (by positivity)]
    have hb : |B1| ≤ x ^ 5 * (23 * L ^ 2) := b1.trans (mul_le_mul_of_nonneg_left s1 (by positivity))
    have hza : 2 * z * Y02.a2 θ ≤ 4.32 * z := by
      have := mul_le_mul_of_nonneg_left ha216 (by positivity : (0:ℝ) ≤ 2 * z)
      linarith
    have := mul_le_mul hza hb (abs_nonneg _) (by positivity)
    calc _ ≤ 4.32 * z * (x ^ 5 * (23 * L ^ 2)) := this
      _ = _ := by rw [hW]; ring
  have c2 : |2 * z * Y02.a3 θ * B2| ≤ 18 * 9.5 * (W * L ^ 2) := by
    rw [abs_mul, abs_mul, abs_of_pos (by positivity : (0:ℝ) < 2 * z)]
    have hb : |B2| ≤ x ^ 5 * (9.5 * L ^ 2) := b2.trans (mul_le_mul_of_nonneg_left s2 (by positivity))
    have := mul_le_mul (mul_le_mul_of_nonneg_left ha3 (by positivity : (0:ℝ) ≤ 2 * z)) hb
      (abs_nonneg _) (by positivity)
    calc _ ≤ 2 * z * 9 * (x ^ 5 * (9.5 * L ^ 2)) := this
      _ = _ := by rw [hW]; ring
  have c3 : |2 * z * b4 * B3| ≤ 40 * 3.25 * (W * L ^ 2) := by
    rw [abs_mul, abs_mul, abs_of_pos (by positivity : (0:ℝ) < 2 * z)]
    have hb : |B3| ≤ x ^ 5 * (3.25 * L ^ 2) := b3.trans (mul_le_mul_of_nonneg_left s3 (by positivity))
    have := mul_le_mul (mul_le_mul_of_nonneg_left hb4 (by positivity : (0:ℝ) ≤ 2 * z)) hb
      (abs_nonneg _) (by positivity)
    calc _ ≤ 2 * z * 20 * (x ^ 5 * (3.25 * L ^ 2)) := this
      _ = _ := by rw [hW]; ring
  have c4 : |2 * z * q * ρ5| ≤ 90.5 * (W * L ^ 2) := by
    rw [abs_mul, abs_of_pos (by positivity : (0:ℝ) < 2 * z * q)]
    have h1 : 2 * z * q ≤ 2 * z := mul_le_of_le_one_right (by positivity) hq1
    have h2 : |ρ5| ≤ 100 * (1.81 * x ^ 5) := hρ.trans (by linarith)
    have := mul_le_mul h1 h2 (abs_nonneg _) (by positivity)
    calc _ ≤ 2 * z * (100 * (1.81 * x ^ 5)) := this
      _ = 362 * W := by rw [hW]; ring
      _ ≤ 90.5 * (W * L ^ 2) := by linarith
  have c5 : |D5| + |D6a| + |D6b| + |T6| ≤ 6.85 * (W * L ^ 2) := by
    have : 21 * z * x ^ 5 + 2.7 * z * x ^ 5 + 1.1 * z * x ^ 5 + 2.6 * z * x ^ 5 = 27.4 * W := by
      rw [hW]; ring
    linarith
  have u1 := abs_add_le (2 * z * Y02.a2 θ * B1 + 2 * z * Y02.a3 θ * B2 + 2 * z * b4 * B3
    + 2 * z * q * ρ5 + D5 + D6a + D6b) T6
  have u2 := abs_add_le (2 * z * Y02.a2 θ * B1 + 2 * z * Y02.a3 θ * B2 + 2 * z * b4 * B3
    + 2 * z * q * ρ5 + D5 + D6a) D6b
  have u3 := abs_add_le (2 * z * Y02.a2 θ * B1 + 2 * z * Y02.a3 θ * B2 + 2 * z * b4 * B3
    + 2 * z * q * ρ5 + D5) D6a
  have u4 := abs_add_le (2 * z * Y02.a2 θ * B1 + 2 * z * Y02.a3 θ * B2 + 2 * z * b4 * B3
    + 2 * z * q * ρ5) D5
  have u5 := abs_add_le (2 * z * Y02.a2 θ * B1 + 2 * z * Y02.a3 θ * B2 + 2 * z * b4 * B3)
    (2 * z * q * ρ5)
  have u6 := abs_add_le (2 * z * Y02.a2 θ * B1 + 2 * z * Y02.a3 θ * B2) (2 * z * b4 * B3)
  have u7 := abs_add_le (2 * z * Y02.a2 θ * B1) (2 * z * Y02.a3 θ * B2)
  have hWL0 : 0 ≤ W * L ^ 2 := by positivity
  linarith

end Z06h


namespace Z06h
open Real

lemma p12 : ∀ θ n : ℝ, 0 < θ → θ < π → 9 ≤ n →
    |Z06.E4 θ n - (-(Y02.zz θ * Y02.a2 θ / 3 * Real.log n ^ 3 / n ^ 5))|
      ≤ 2.5e5 * Y02.zz θ * Real.log n ^ 2 / n ^ 5 := by
  intro θ n h0 h1 hn
  have z0 : 0 < Y02.zz θ := (Y02.zz_bounds θ h0 h1).1
  have H := E4_core θ n h0 h1 hn
  have hn0 : 0 < n := by linarith
  rw [sub_neg_eq_add]
  refine H.trans ?_
  have e : 2.5e5 * Y02.zz θ * Real.log n ^ 2 / n ^ 5
      = 2.5e5 * (Y02.zz θ * (1 / n) ^ 5 * Real.log n ^ 2) := by field_simp
  rw [e]
  have : 0 ≤ Y02.zz θ * (1 / n) ^ 5 * Real.log n ^ 2 := by positivity
  linarith

lemma p11 : ∀ θ n : ℝ, 0 < θ → θ < π → 9 ≤ n →
    |Z06.E4 θ n| ≤ 1.2e5 * Y02.zz θ * Real.log n ^ 3 / n ^ 5 := by
  intro θ n h0 h1 hn
  have z0 : 0 < Y02.zz θ := (Y02.zz_bounds θ h0 h1).1
  have H := E4_core θ n h0 h1 hn
  have ha0 := (Y02.proof_9 θ h0 h1).2.2.2.2.1
  have ha216 := a2_le θ h0 h1
  have hL2 := Y02.two_le_log hn
  have hn0 : 0 < n := by linarith
  set L := Real.log n
  set z := Y02.zz θ
  set c := z * Y02.a2 θ / 3 * L ^ 3 / n ^ 5
  set X := z * (1 / n) ^ 5 with hX
  have hX0 : 0 ≤ X := by positivity
  have hc : c = Y02.a2 θ / 3 * (X * L ^ 3) := by
    show z * Y02.a2 θ / 3 * L ^ 3 / n ^ 5 = Y02.a2 θ / 3 * (z * (1 / n) ^ 5 * L ^ 3)
    field_simp
  have hc0 : 0 ≤ c := by rw [hc]; positivity
  have hc1 : c ≤ 0.72 * (X * L ^ 3) := by
    rw [hc]; apply mul_le_mul_of_nonneg_right (by linarith) (by positivity)
  have hL3 : X * L ^ 2 * 2 ≤ X * L ^ 3 := by
    have : X * L ^ 2 * 2 ≤ X * L ^ 2 * L := mul_le_mul_of_nonneg_left hL2 (by positivity)
    linarith [show X * L ^ 2 * L = X * L ^ 3 by ring]
  have e : Z06.E4 θ n = (Z06.E4 θ n + c) - c := by ring
  rw [e]
  have u := abs_sub (Z06.E4 θ n + c) c
  rw [abs_of_nonneg hc0] at u
  have hfin : 1.2e5 * z * L ^ 3 / n ^ 5 = 1.2e5 * (X * L ^ 3) := by rw [hX]; field_simp
  rw [hfin]
  have : 1000 * (X * L ^ 2) ≤ 500 * (X * L ^ 3) := by linarith
  have h0' : 0 ≤ X * L ^ 3 := by positivity
  linarith

end Z06h

theorem Z06.proof_1 : Z06.stmt_1 := Z06h.p1
theorem Z06.proof_2 : Z06.stmt_2 := Z06h.p2
theorem Z06.proof_3 : Z06.stmt_3 := Z06h.p3
theorem Z06.proof_4 : Z06.stmt_4 := Z06h.p4
theorem Z06.proof_5 : Z06.stmt_5 := Z06h.p5
theorem Z06.proof_7 : Z06.stmt_7 := Z06h.p7
theorem Z06.proof_6 : Z06.stmt_6 := Z06h.p6
theorem Z06.proof_8 : Z06.stmt_8 := Z06h.p8
theorem Z06.proof_9 : Z06.stmt_9 := Z06h.p9 Z06h.p7
theorem Z06.proof_10 : Z06.stmt_10 := Z06h.p10
theorem Z06.proof_11 : Z06.stmt_11 := Z06h.p11
theorem Z06.proof_12 : Z06.stmt_12 := Z06h.p12
