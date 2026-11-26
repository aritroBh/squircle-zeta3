import Mathlib
import Zeta3.Wave3.W02_Holo
/-!
# Q05 -- the decay law of `M_k` (Theorem `thm:decay`, equation `eq:Mk`)

Target id: Q05.
Paper labels: `thm:decay` (paper/ieee/main.tex, lines 1287-1297), `eq:Mk` (main.tex, lines
1276-1278), the bound on `M_k^main` (main.tex, lines 1281-1284), the definition of `M_k^main`
(main.tex, line 1313). Detailed proof: paper/ieee/proofs/sec_g5_decay.tex (lemmas
`g5:lem:kernels`, `g5:lem:logit`, `g5:lem:kernel`, `g5:lem:Mk`, `g5:lem:main`, `g5:lem:W`,
`g5:lem:F`, `g5:prop:rest`, the proof of `g5:thm:decay`) and, for the main term,
`g5:prop:main` of paper/ieee/proofs/sec_g5_varpi1.tex (line 249).

Quoted verbatim (main.tex, thm:decay):

  For every $k$ with $N=k+1\ge100$,
  \Bigl|M_k+\frac{\sqrt\pi}{N}\,e^{-\sqrt{\pi N}}\cos\sqrt{\pi N}\Bigr|\le1.21\cdot10^{9}\,N^{-3/2}e^{-\sqrt{\pi N}}.
  In particular \eqref{eq:Mk} holds.

Quoted verbatim (main.tex, eq:Mk):

  M_k=-\frac{\sqrt\pi}{k}\,e^{-\sqrt{\pi k}}\Bigl[\cos\sqrt{\pi k}+O\bigl(k^{-1/2}\bigr)\Bigr].

Quoted verbatim (main.tex, lines 1281-1284, the bound on the main term):

  A computer-assisted proof gives, for all $N\ge100$,
  \Bigl|M_k^{\mathrm{main}}+\frac{\sqrt\pi}{N}\,e^{-\sqrt{\pi N}}\cos\sqrt{\pi N}\Bigr|\le254.23\,N^{-3/2}e^{-\sqrt{\pi N}}.

Quoted verbatim (main.tex, proof idea of thm:decay):

  Add the bound on $M_k^{\mathrm{main}}$ and an earlier, weaker rest bound,
  $|M_k-M_k^{\mathrm{main}}|\le10{,}600{,}000\,e^{-2\sqrt N}$ for $N\ge100$, [...] Since
  $2>\sqrt\pi$, the factor $e^{-(2-\sqrt\pi)\sqrt N}$ and
  $\sup_{u>0}u^3e^{-(2-\sqrt\pi)u}=\bigl(3/(e(2-\sqrt\pi))\bigr)^3<114.1$ give
  $1.06\cdot10^7\cdot114.1+254.23<1.21\cdot10^9$.

Quoted verbatim (main.tex, line 1313):

  Here $M_k^{\mathrm{main}}=\frac2\pi\operatorname{Im}\int_0^\infty H(s)e^{i\pi s/2}s^{-1}\bigl(\frac s{1+s}\bigr)^N ds$
  with $H(s)=\frac12B(\frac s2,\frac32)$, exactly.

Quoted verbatim (supplement, the lemmas; `\Lline` is the half line `{1/2 + r ω : r ≥ 0}`,
`ω = e^{iπ/4}`, `E(s)=e^{iπs/2}`, `K(s)=1/(2cos(πs/2))`, `Q_k(s)=s^k/(1+s)^{k+1}`):

  g5:lem:kernels. (a) $H(s)=\int_0^\infty\sqrt{1-\eu^{-2u}}\,\eu^{-su}\dd u$ for
  $\operatorname{Re}s>0$, $H$ is holomorphic there, and with $x=\operatorname{Re}s$
  |H(s)|\le\frac1x,\qquad|H(s)|\le\sqrt{\tfrac\pi2}\,x^{-3/2}.
  Moreover $H(2n+1)=\frac{\sqrt\pi}4\,\Gamma(n+\frac12)/(n+1)!$ for $n\ge0$.
  (b) $K$ is meromorphic with simple poles at the odd integers $2n+1$, with residue
  $(-1)^{n+1}/\pi$ there.
  (c) On $\Lline$: $|H|\le2$, $|E|=\eu^{-ar}$, $|K|\le2\,\eu^{-ar}$ and $|K-E|\le2\,\eu^{-3ar}$,
  with $a=\pi/(2\sqrt2)$.

  g5:lem:logit. For real $\varepsilon<1$ put $C(z)=1-\eu^{-z}/\bigl(1+\sqrt{1+\eu^{-2z}}\bigr)$ for
  $\operatorname{Re}z\ge0$. Then
  M(\varepsilon)=-\int_0^\infty\eu^{-x}(1+\eu^{-x})^{-1-\varepsilon}\,C\bigl((1-\varepsilon)x\bigr)\dd x ,
  and the integral converges absolutely.

  g5:lem:kernel. For every real $z\ge0$,
  $\;C(z)=1-\frac2\pi\operatorname{Im}\int_{\Lline}H(s)K(s)\,\eu^{-sz}\,ds$, and the integral
  converges absolutely.

  (definitions) P_k(s)=\frac1{k!}\partial_\varepsilon^k\Phi(0,s)=\frac{(-1)^k}{k!}\int_0^1\frac{t^s}{1+t}\bigl(s\log t+\log(1+t)\bigr)^kdt,
  $D(\varepsilon)=\int_0^1(1+t)^{-1-\varepsilon}dt=(1-2^{-\varepsilon})/\varepsilon$ and
  $D_k=[\varepsilon^k]D=(-1)^k(\log2)^{k+1}/(k+1)!$.

  g5:lem:Mk. For every $k\ge0$,
  M_k=-D_k+\frac2\pi\operatorname{Im}\int_{\Lline}H(s)K(s)P_k(s)\,ds ,
  and the integral converges absolutely.

  g5:lem:main. For $k\ge1$,
  $\displaystyle\int_0^\infty H(s)E(s)Q_k(s)\,ds=\int_0^{1/2}HEQ_k\,ds+\int_{\Lline}HEQ_k\,ds$.

  g5:lem:W. Let $k\ge1$ and let $s$ satisfy $\operatorname{Re}s>0$ and $0\le\arg s\le\pi/4$.
  Put $\rho=|s|$. Then
  |W_k(s)|\le\frac{2\rho^k}{|s+2|^{k+1}}+2\sum_{i=1}^k\frac{2^i}{i!}\,\frac{\rho^{k-i}}{|s+i+1|^{k-i+1}} .
  (with $W_k=P_k-Q_k$)

  g5:lem:F. For an integer $j\ge0$, a real $q\ge1$ and $\lambda>0$ put
  $\mathcal F_{j,q}(\lambda)=\int_0^\infty\eu^{-\lambda r}\rho^j|s+q|^{-j-1}dr$ with
  $s=\frac12+r\omega\in\Lline$, $\rho=|s|$. Then
  \mathcal F_{j,q}(\lambda)\le\frac4\lambda\exp\Bigl(\frac\lambda2+\frac{\lambda q}{\sqrt2}\Bigr)\exp\Bigl(-\sqrt{3\lambda jq/\sqrt2}\Bigr).

  Constants. Put $\xi=2\eu^{\pi/4}$ and
  B_1=\frac{64}{\pi a}\,\eu^{a/2+\pi/4},\qquad B_2=\frac{32}{3\pi a}\,\eu^{3a/2+3\pi/4},\qquad C_\star=2+B_2+B_1\bigl(\eu^{\pi/4}+\eu^{\xi}+\eu^{\mathrm e\xi}\bigr).

  g5:prop:rest. For every $k$ with $N=k+1\ge100$,
  |M_k-M_k^{\mathrm{main}}|\le{}&B_1\Bigl[\eu^{\pi/4}\eu^{-\sqrt{3\pi k/2}}+\xi\,\eu^{-\sqrt{3\pi(k-1)/2}}
  &+(\eu^{\xi}-1-\xi)\,\eu^{-\sqrt{9\pi(k-2)/4}}+\eu^{\mathrm e\xi}\eu^{-k/2}\Bigr]+B_2\,\eu^{-\sqrt{9\pi k/4}}+2\eu^{-(N-1)/2}.
  Consequently: (a) $|M_k-M_k^{\mathrm{main}}|\le C_\star\,\eu^{-2\sqrt N}$, and $C_\star<10{,}600{,}000$.

  g5:prop:main (sec_g5_varpi1.tex). For every integer $N\ge100$,
  $|M_k^{\mathrm{main}}+X|\le71{,}000\,N^{-2}\eu^{-\theta}$,
  with $\theta=\sqrt{\pi N}$, $X=\frac{\sqrt\pi}N\eu^{-\theta}\operatorname{Re}\bigl[\eu^{\ii\theta}(1+\varpi_1/\sqrt N)\bigr]$,
  $\varpi_1=\bigl[(24+12\pi+\pi^2)+\ii(24-12\pi+\pi^2)\bigr]/(48\sqrt\pi)$.

  Proof of g5:thm:decay (steps used): "The function $u^3\eu^{-\delta u}$ on $u>0$ has its maximum
  $(3/(\mathrm e\delta))^3$ at $u=3/\delta$", "$(3/(\mathrm e\delta))^3=114.0963\ldots<114.1$";
  "$\frac{\sqrt\pi}{N^{3/2}}\eu^{-\theta}|\operatorname{Re}[\eu^{\ii\theta}\varpi_1]|\le\sqrt\pi\,|\varpi_1|\,N^{-3/2}\eu^{-\theta}\le1.494\,N^{-3/2}\eu^{-\theta}$";
  "$\bigl(1{,}209{,}421{,}579.87+7100+1.5\bigr)N^{-3/2}\eu^{-\theta}[...]\le1.21\cdot10^9\,N^{-3/2}\eu^{-\theta}$";
  and for eq:Mk: "\Bigl|M_k+\frac{\sqrt\pi}k\eu^{-\sqrt{\pi k}}\cos\sqrt{\pi k}\Bigr|\le\bigl(1.21\cdot10^9+\sqrt\pi\cdot1.9\bigr)k^{-3/2}\eu^{-\sqrt{\pi k}}"
  for $k\ge99$.

## Conventions

* `M_k` is the real number `M k := (W02.Mk k).re`, with `W02.Mk k = -F^{(k)}(0)/(2 k!)` the
  Taylor coefficient of `M = -F/2` (W02). Realness is not assumed: Q07.stmt_2 (old stmt_10, g5:lem:Mk) states
  `W02.Mk k` itself, as a complex number, equal to a real number, so it contains `Im M_k = 0`.
* For real `ε < 1`, the paper's `M(ε) = -F(ε)/2` is `-(W02.FR ε)/2`, the real integral of the
  paper's eq:F (W02.stmt_4 shows it agrees with the holomorphic extension `W02.FC`).
* `H s = (1/2) B(s/2, 3/2)` with Mathlib's `Complex.betaIntegral` (the Euler integral, which is
  the Beta function for positive real parts). `E`, `K`, `Q k` as in the supplement; `K` is
  `1/(2 cos(πs/2))` and takes Lean's junk value `0` at the poles, which no statement uses.
* The half line `\Lline` is parametrised as `Lpt r = 1/2 + r ω`, `r > 0`, and
  `∫_\Lline f(s) ds = ∫_{r>0} f(Lpt r) ω dr` (`lineInt`). Every line, half-line or interval
  integral whose value is asserted is accompanied by an `Integrable`/`IntegrableOn` claim
  ("converges absolutely"), so nothing holds through the junk value `0`.
* `P k s` is the explicit integral of the supplement (the second expression of its definition);
  `(t:ℂ)^s` is the principal power of a positive real `t ∈ (0,1)`.
* `M_k^main` is `Mmain k = (2/π) Im ∫_{s>0} H(s) e^{iπs/2} s^{-1} (s/(1+s))^N ds` exactly as
  main.tex line 1313 (real `s` cast to `ℂ`). For `k ≥ 1` its integrability follows from stmt_11
  (g5:lem:main: integrability of `H E Q_k` on `(0,∞)` and `s^{-1}(s/(1+s))^N = Q_k(s)`).
* `N^{-3/2}` is the real power `(N:ℝ) ^ (-(3/2) : ℝ)` of a positive base; `N^{-2}` is
  `1 / (N:ℝ)^2`.
* `D k = (-1)^k (log 2)^{k+1}/(k+1)!` is the closed form the supplement gives for `[ε^k]D`.

## Statements (in the order: elementary first, then the supplement's order)

* stmt_1: the supremum: `u^3 e^{-(2-√π)u} ≤ (3/(e(2-√π)))^3` for `u > 0`, attained at
  `u = 3/(2-√π)`, and `(3/(e(2-√π)))^3 < 114.1`.
* stmt_2: the arithmetic: `1.06·10^7·114.1 + 254.23 < 1.21·10^9` (paper) and
  `10,600,000·114.1 + 7100 + 1.5 < 1.21·10^9` (supplement).
* stmt_3: `√π |ϖ₁| ≤ 1.494`.
* stmt_4: `C_⋆ < 10,600,000`.
* stmt_5: g5:lem:kernels (a).
* stmt_6: g5:lem:kernels (b).
* stmt_7: g5:lem:kernels (c).
* stmt_8: g5:lem:logit.
* stmt_11: g5:lem:main (the contour shift, with the integrability of the three integrals).
* stmt_13: g5:lem:F.

## NOT COVERED

* The statement "$K$ is meromorphic" is stated as: `K` is complex differentiable off the odd
  integers, and `(s-(2n+1))K(s) → (-1)^{n+1}/π` at each odd integer (stmt_6); this gives the
  simple poles with these residues.
* The intermediate formula `M(ε)=-D(ε)+\frac2\pi Im∫_\Lline HKΦ(ε,s)ds` for real ε, the
  identity `P_k = \frac1{k!}∂_ε^kΦ(0,s)`, `D(ε)=(1-2^{-ε})/ε`, and the second display of
  g5:lem:main (the decomposition of `M_k - M_k^main` into `D_k, U_1, U_2, U_3`, which is algebra
  from old stmt_10 (now Q07.stmt_2) and stmt_11) are not stated separately.
* g5:prop:rest (b) (`c_rest < 120,000`, used only for thm:varpi) is not stated.
* Lemmas g5:lem:kernel, g5:lem:Mk, g5:lem:W, Propositions g5:prop:rest and g5:prop:main and Theorem thm:decay itself were audited as stmt_9, 10, 12, 14, 15, 16, 18 and moved, unproved, to Q07_DecayRest.lean; the interval-arithmetic bound with constant 254.23 (audited as stmt_17) has no written proof in the supplement and is not stated any more (planner split, 2026-11-26).
-/
-- BEGIN STATEMENTS
namespace Q05

open Real MeasureTheory Set Filter

/-- `M_k`, the coefficient of `1/p^{k+1}` (W02), as a real number. -/
noncomputable def M (k : ℕ) : ℝ := (W02.Mk k).re

/-- `H(s) = (1/2) B(s/2, 3/2)`. -/
noncomputable def H (s : ℂ) : ℂ := (1 / 2) * Complex.betaIntegral (s / 2) (3 / 2)

/-- `E(s) = e^{iπs/2}`. -/
noncomputable def E (s : ℂ) : ℂ := Complex.exp (Complex.I * π * s / 2)

/-- `K(s) = 1/(2 cos(πs/2))`. -/
noncomputable def K (s : ℂ) : ℂ := 1 / (2 * Complex.cos (π * s / 2))

/-- `Q_k(s) = s^k/(1+s)^{k+1}`. -/
noncomputable def Q (k : ℕ) (s : ℂ) : ℂ := s ^ k / (1 + s) ^ (k + 1)

/-- `ω = e^{iπ/4}`. -/
noncomputable def ω : ℂ := Complex.exp (π / 4 * Complex.I)

/-- The point `s = 1/2 + r ω` of the half line `\Lline`. -/
noncomputable def Lpt (r : ℝ) : ℂ := 1 / 2 + r * ω

/-- `∫_\Lline f(s) ds = ∫_0^∞ f(1/2 + rω) ω dr`. -/
noncomputable def lineInt (f : ℂ → ℂ) : ℂ := ∫ r in Ioi (0 : ℝ), f (Lpt r) * ω

/-- `a = π/(2√2)`. -/
noncomputable def a : ℝ := π / (2 * √2)

/-- `C(z) = 1 - e^{-z}/(1+√(1+e^{-2z}))`, for real `z`. -/
noncomputable def Cf (z : ℝ) : ℝ := 1 - Real.exp (-z) / (1 + √(1 + Real.exp (-2 * z)))

/-- `P_k(s) = ((-1)^k/k!) ∫_0^1 t^s/(1+t) (s log t + log(1+t))^k dt`. -/
noncomputable def P (k : ℕ) (s : ℂ) : ℂ :=
  (-1) ^ k / (Nat.factorial k : ℂ) *
    ∫ t in (0 : ℝ)..1, (t : ℂ) ^ s / (1 + (t : ℂ)) *
      (s * (Real.log t : ℂ) + (Real.log (1 + t) : ℂ)) ^ k

/-- `D_k = (-1)^k (log 2)^{k+1}/(k+1)!`. -/
noncomputable def D (k : ℕ) : ℝ := (-1) ^ k * Real.log 2 ^ (k + 1) / (Nat.factorial (k + 1) : ℝ)

/-- `W_k = P_k - Q_k`. -/
noncomputable def W (k : ℕ) (s : ℂ) : ℂ := P k s - Q k s

/-- `M_k^main = (2/π) Im ∫_0^∞ H(s) e^{iπs/2} s^{-1} (s/(1+s))^N ds`, `N = k+1`. -/
noncomputable def Mmain (k : ℕ) : ℝ :=
  2 / π * (∫ s in Ioi (0 : ℝ), H s * E s * (s : ℂ)⁻¹ * ((s : ℂ) / (1 + s)) ^ (k + 1)).im

/-- The integrand `e^{-λr} ρ^j |s+q|^{-j-1}` of `𝓕_{j,q}(λ)`, `s = 1/2 + rω`, `ρ = |s|`. -/
noncomputable def FjqIntegrand (j : ℕ) (q l : ℝ) (r : ℝ) : ℝ :=
  Real.exp (-l * r) * ‖Lpt r‖ ^ j / ‖Lpt r + q‖ ^ (j + 1)

/-- `𝓕_{j,q}(λ) = ∫_0^∞ e^{-λr} ρ^j |s+q|^{-j-1} dr`. -/
noncomputable def Fjq (j : ℕ) (q l : ℝ) : ℝ :=
  ∫ r in Ioi (0 : ℝ), FjqIntegrand j q l r

/-- `ξ = 2 e^{π/4}`. -/
noncomputable def ξ : ℝ := 2 * Real.exp (π / 4)

/-- `B_1 = (64/(πa)) e^{a/2+π/4}`. -/
noncomputable def B1 : ℝ := 64 / (π * a) * Real.exp (a / 2 + π / 4)

/-- `B_2 = (32/(3πa)) e^{3a/2+3π/4}`. -/
noncomputable def B2 : ℝ := 32 / (3 * π * a) * Real.exp (3 * a / 2 + 3 * π / 4)

/-- `C_⋆ = 2 + B_2 + B_1 (e^{π/4} + e^ξ + e^{eξ})`. -/
noncomputable def Cstar : ℝ :=
  2 + B2 + B1 * (Real.exp (π / 4) + Real.exp ξ + Real.exp (Real.exp 1 * ξ))

/-- `ϖ₁ = [(24+12π+π²) + i(24-12π+π²)]/(48√π)`. -/
noncomputable def varpi1 : ℂ :=
  ⟨(24 + 12 * π + π ^ 2) / (48 * √π), (24 - 12 * π + π ^ 2) / (48 * √π)⟩

/-- `X = (√π/N) e^{-θ} Re[e^{iθ}(1+ϖ₁/√N)]`, `θ = √(πN)`, with `N` real. -/
noncomputable def X (N : ℝ) : ℝ :=
  √π / N * Real.exp (-√(π * N)) *
    (Complex.exp (Complex.I * (√(π * N) : ℝ)) * (1 + varpi1 / (√N : ℝ))).re

/-- Proof of g5:thm:decay (and main.tex, proof idea): "$\sup_{u>0}u^3e^{-(2-\sqrt\pi)u}=
\bigl(3/(e(2-\sqrt\pi))\bigr)^3<114.1$": the bound for every `u > 0`, the value at
`u = 3/(2-√π)`, and the numerical inequality. -/
def stmt_1 : Prop :=
  (∀ u : ℝ, 0 < u →
      u ^ 3 * Real.exp (-(2 - √π) * u) ≤ (3 / (Real.exp 1 * (2 - √π))) ^ 3) ∧
    (3 / (2 - √π)) ^ 3 * Real.exp (-(2 - √π) * (3 / (2 - √π))) =
      (3 / (Real.exp 1 * (2 - √π))) ^ 3 ∧
    (3 / (Real.exp 1 * (2 - √π))) ^ 3 < 114.1

/-- main.tex proof idea: "$1.06\cdot10^7\cdot114.1+254.23<1.21\cdot10^9$", and the sum of the
supplement's proof with `7100` (= 71,000/10) and `1.5` in place of `254.23`. -/
def stmt_2 : Prop :=
  (1.06e7 : ℝ) * 114.1 + 254.23 < 1.21e9 ∧ (10600000 : ℝ) * 114.1 + 7100 + 1.5 < 1.21e9

/-- Proof of g5:thm:decay, the last term: "$\sqrt\pi\,|\varpi_1|\le1.494$". -/
def stmt_3 : Prop := √π * ‖varpi1‖ ≤ 1.494

/-- g5:prop:rest (a), the constant: "$C_\star<10{,}600{,}000$". -/
def stmt_4 : Prop := Cstar < 10600000

/-- g5:lem:kernels (a): for `Re s > 0`, `H(s) = ∫_0^∞ √(1-e^{-2u}) e^{-su} du` (absolutely
convergent), `H` is holomorphic on `Re s > 0`, `|H(s)| ≤ 1/x`, `|H(s)| ≤ √(π/2) x^{-3/2}`
(`x = Re s`), and `H(2n+1) = (√π/4) Γ(n+1/2)/(n+1)!` for `n ≥ 0`. -/
def stmt_5 : Prop :=
  (∀ s : ℂ, 0 < s.re →
      IntegrableOn (fun u : ℝ => ((√(1 - Real.exp (-2 * u)) : ℝ) : ℂ) * Complex.exp (-s * u))
        (Ioi 0) ∧
      H s = ∫ u in Ioi (0 : ℝ), ((√(1 - Real.exp (-2 * u)) : ℝ) : ℂ) * Complex.exp (-s * u)) ∧
    DifferentiableOn ℂ H {s : ℂ | 0 < s.re} ∧
    (∀ s : ℂ, 0 < s.re → ‖H s‖ ≤ 1 / s.re ∧ ‖H s‖ ≤ √(π / 2) * s.re ^ (-(3 / 2) : ℝ)) ∧
    (∀ n : ℕ, H (2 * n + 1) =
      ((√π / 4 * Real.Gamma (n + 1 / 2) / (Nat.factorial (n + 1) : ℝ) : ℝ) : ℂ))

/-- g5:lem:kernels (b): "$K$ is meromorphic with simple poles at the odd integers $2n+1$, with
residue $(-1)^{n+1}/\pi$ there": `K` is complex differentiable at every point that is not an
odd integer, and `(s - (2n+1)) K(s) → (-1)^{n+1}/π` as `s → 2n+1`, `s ≠ 2n+1`, for every
integer `n`. -/
def stmt_6 : Prop :=
  (∀ s : ℂ, (∀ n : ℤ, s ≠ 2 * n + 1) → DifferentiableAt ℂ K s) ∧
    ∀ n : ℤ, Tendsto (fun s : ℂ => (s - (2 * n + 1)) * K s)
      (nhdsWithin ((2 * n + 1 : ℂ)) {(2 * n + 1 : ℂ)}ᶜ) (nhds ((-1 : ℂ) ^ (n + 1) / π))

/-- g5:lem:kernels (c): "On $\Lline$: $|H|\le2$, $|E|=\eu^{-ar}$, $|K|\le2\,\eu^{-ar}$ and
$|K-E|\le2\,\eu^{-3ar}$, with $a=\pi/(2\sqrt2)$." -/
def stmt_7 : Prop :=
  ∀ r : ℝ, 0 ≤ r →
    ‖H (Lpt r)‖ ≤ 2 ∧ ‖E (Lpt r)‖ = Real.exp (-a * r) ∧ ‖K (Lpt r)‖ ≤ 2 * Real.exp (-a * r) ∧
      ‖K (Lpt r) - E (Lpt r)‖ ≤ 2 * Real.exp (-3 * a * r)

/-- g5:lem:logit: for real `ε < 1`,
`M(ε) = -∫_0^∞ e^{-x}(1+e^{-x})^{-1-ε} C((1-ε)x) dx`, absolutely convergent; here
`M(ε) = -F(ε)/2` with the paper's real `F` (`W02.FR`). -/
def stmt_8 : Prop :=
  ∀ ε : ℝ, ε < 1 →
    IntegrableOn
        (fun x : ℝ => Real.exp (-x) * (1 + Real.exp (-x)) ^ (-1 - ε) * Cf ((1 - ε) * x))
        (Ioi 0) ∧
      -(W02.FR ε) / 2 =
        -∫ x in Ioi (0 : ℝ), Real.exp (-x) * (1 + Real.exp (-x)) ^ (-1 - ε) * Cf ((1 - ε) * x)

/-- g5:lem:main: "For $k\ge1$, $\int_0^\infty H(s)E(s)Q_k(s)\,ds=\int_0^{1/2}HEQ_k\,ds+
\int_{\Lline}HEQ_k\,ds$", with the absolute convergence of the integrals over `(0,∞)` and over
`\Lline`, and the identity `s^{-1}(s/(1+s))^N = Q_k(s)` linking it to `M_k^main`. -/
def stmt_11 : Prop :=
  ∀ k : ℕ, 1 ≤ k →
    IntegrableOn (fun s : ℝ => H s * E s * Q k s) (Ioi 0) ∧
      IntegrableOn (fun r : ℝ => H (Lpt r) * E (Lpt r) * Q k (Lpt r)) (Ioi 0) ∧
      (∀ s : ℝ, 0 < s → (s : ℂ)⁻¹ * ((s : ℂ) / (1 + s)) ^ (k + 1) = Q k s) ∧
      ∫ s in Ioi (0 : ℝ), H s * E s * Q k s =
        (∫ s in (0 : ℝ)..(1 / 2), H s * E s * Q k s) + lineInt fun s => H s * E s * Q k s

/-- g5:lem:F: for `j ≥ 0`, real `q ≥ 1`, `λ > 0`, the integral `𝓕_{j,q}(λ)` converges and
`𝓕_{j,q}(λ) ≤ (4/λ) exp(λ/2 + λq/√2) exp(-√(3λjq/√2))`. -/
def stmt_13 : Prop :=
  ∀ (j : ℕ) (q l : ℝ), 1 ≤ q → 0 < l →
    IntegrableOn (FjqIntegrand j q l) (Ioi 0) ∧
      Fjq j q l ≤ 4 / l * Real.exp (l / 2 + l * q / √2) * Real.exp (-√(3 * l * j * q / √2))

end Q05
-- END STATEMENTS

open Real MeasureTheory Set Filter

namespace Q05
open Real

lemma sqrt_pi_gt : (1.7724536 : ℝ) < √π := by
  rw [Real.lt_sqrt (by norm_num)]; nlinarith [Real.pi_gt_d6]

lemma sqrt_pi_lt : √π < 1.7724539 := by
  rw [Real.sqrt_lt' (by norm_num)]; nlinarith [Real.pi_lt_d20]

lemma xexp_le (x : ℝ) : x * Real.exp (-x) ≤ Real.exp (-1) := by
  have h := Real.add_one_le_exp (x - 1)
  have h2 : Real.exp (x - 1) * Real.exp (-x) = Real.exp (-1) := by
    rw [← Real.exp_add]; ring_nf
  have : 0 < Real.exp (-x) := Real.exp_pos _
  nlinarith

end Q05

theorem Q05.proof_1 : Q05.stmt_1 := by
  have hd : 0 < 2 - √π := by linarith [Q05.sqrt_pi_lt]
  have hdne : 2 - √π ≠ 0 := hd.ne'
  refine ⟨fun u hu => ?_, ?_, ?_⟩
  · set d := 2 - √π
    have h := Q05.xexp_le (d * u / 3)
    have h0 : 0 ≤ d * u / 3 * Real.exp (-(d * u / 3)) := by positivity
    have h3 := pow_le_pow_left₀ h0 h 3
    have e1 : (d * u / 3 * Real.exp (-(d * u / 3))) ^ 3 =
        (d / 3) ^ 3 * (u ^ 3 * Real.exp (-d * u)) := by
      rw [mul_pow, ← Real.exp_nat_mul]; ring_nf
    have e2 : Real.exp (-1) ^ 3 = (d / 3) ^ 3 * (3 / (Real.exp 1 * d)) ^ 3 := by
      rw [Real.exp_neg, div_pow, div_pow, mul_pow]; field_simp
    rw [e1, e2] at h3
    exact le_of_mul_le_mul_left h3 (by positivity)
  · rw [show -(2 - √π) * (3 / (2 - √π)) = -3 by field_simp]
    rw [show (-3 : ℝ) = (3 : ℕ) * (-1) by norm_num, Real.exp_nat_mul, Real.exp_neg]
    field_simp
  · have he := Real.exp_one_gt_d9
    have hs := Q05.sqrt_pi_lt
    have hlt : 3 / (Real.exp 1 * (2 - √π)) < 4.8502 := by
      have hd' : (0.2275461 : ℝ) < 2 - √π := by linarith
      rw [div_lt_iff₀ (by positivity)]
      have := mul_lt_mul'' he hd' (by norm_num) (by norm_num)
      nlinarith
    have h0 : 0 ≤ 3 / (Real.exp 1 * (2 - √π)) := by positivity
    calc (3 / (Real.exp 1 * (2 - √π))) ^ 3 ≤ (4.8502 : ℝ) ^ 3 :=
          pow_le_pow_left₀ h0 hlt.le 3
      _ < 114.1 := by norm_num

theorem Q05.proof_2 : Q05.stmt_2 := by
  unfold Q05.stmt_2; norm_num

theorem Q05.proof_3 : Q05.stmt_3 := by
  unfold Q05.stmt_3
  have hp := Real.pi_pos
  have hs : 0 < √π := Real.sqrt_pos.2 hp
  have hsq : √π ^ 2 = π := Real.sq_sqrt hp.le
  rw [Complex.norm_eq_sqrt_sq_add_sq, ← Real.sqrt_mul_self (by norm_num : (0:ℝ) ≤ 1.494),
    ← Real.sqrt_mul hp.le]
  apply Real.sqrt_le_sqrt
  simp only [Q05.varpi1]
  rw [div_pow, div_pow, mul_pow, hsq]
  have : π * (((24 + 12 * π + π ^ 2) ^ 2 / (48 ^ 2 * π)) +
      ((24 - 12 * π + π ^ 2) ^ 2 / (48 ^ 2 * π))) =
      ((24 + 12 * π + π ^ 2) ^ 2 + (24 - 12 * π + π ^ 2) ^ 2) / 2304 := by
    field_simp; ring
  rw [this, div_le_iff₀ (by norm_num)]
  nlinarith [Real.pi_lt_d6, Real.pi_gt_d6]

namespace Q05

lemma exp_le_pow (x c : ℝ) (n : ℕ) (hn : 0 < n) (h0 : 0 ≤ c) (hx : x ≤ c) (hc : c ≤ n) :
    Real.exp x ≤ ((∑ m ∈ Finset.range 10, (c / n) ^ m / (m.factorial : ℝ)) +
      (c / n) ^ 10 * (10 + 1) / ((Nat.factorial 10 : ℝ) * 10)) ^ n := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have h1 : Real.exp x ≤ Real.exp (c / n) ^ n := by
    rw [← Real.exp_nat_mul, mul_div_cancel₀ _ hn'.ne']; exact Real.exp_le_exp.2 hx
  refine h1.trans (pow_le_pow_left₀ (Real.exp_pos _).le ?_ n)
  have := Real.exp_bound' (x := c / n) (by positivity) (by rwa [div_le_one hn']) (n := 10)
    (by norm_num)
  exact_mod_cast this

end Q05

theorem Q05.proof_4 : Q05.stmt_4 := by
  unfold Q05.stmt_4 Q05.Cstar Q05.B1 Q05.B2 Q05.ξ
  have hpl : (3.1415926 : ℝ) < π := by linarith [Real.pi_gt_d20]
  have hph : π < (3.1415927 : ℝ) := by linarith [Real.pi_lt_d20]
  have hsl : (1.4142135 : ℝ) < √2 := by
    rw [Real.lt_sqrt (by norm_num)]; norm_num
  have hsh : √2 < (1.4142136 : ℝ) := by
    rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have hal : (1.110720 : ℝ) < Q05.a := by
    unfold Q05.a; rw [lt_div_iff₀ (by positivity)]; nlinarith
  have hah : Q05.a < (1.110721 : ℝ) := by
    unfold Q05.a; rw [div_lt_iff₀ (by positivity)]; nlinarith
  have hpa : (3.1415926 * 1.110720 : ℝ) < π * Q05.a := by
    have := mul_lt_mul'' hpl hal (by norm_num) (by norm_num); linarith
  have hpa0 : 0 < π * Q05.a := by nlinarith
  have hapos : 0 < Q05.a := by linarith
  have hc1 : 64 / (π * Q05.a) < 18.3412 := by
    rw [div_lt_iff₀ (by nlinarith)]; nlinarith
  have hc2 : 32 / (3 * π * Q05.a) < 3.0569 := by
    rw [div_lt_iff₀ (by nlinarith)]; nlinarith
  have he1 : Real.exp (Q05.a / 2 + π / 4) < 3.82195 :=
    lt_of_le_of_lt (Q05.exp_le_pow _ 1.34076 2 (by norm_num) (by norm_num) (by linarith)
      (by norm_num)) (by norm_num)
  have he2 : Real.exp (3 * Q05.a / 2 + 3 * π / 4) < 55.83 :=
    lt_of_le_of_lt (Q05.exp_le_pow _ 4.02228 5 (by norm_num) (by norm_num) (by linarith)
      (by norm_num)) (by norm_num)
  have he3 : Real.exp (π / 4) < 2.1933 :=
    lt_of_le_of_lt (Q05.exp_le_pow _ 0.7854 1 (by norm_num) (by norm_num) (by linarith)
      (by norm_num)) (by norm_num)
  have hxi : 2 * Real.exp (π / 4) < 4.3866 := by linarith
  have he4 : Real.exp (2 * Real.exp (π / 4)) < 80.37 :=
    lt_of_le_of_lt (Q05.exp_le_pow _ 4.3866 5 (by norm_num) (by norm_num) hxi.le
      (by norm_num)) (by norm_num)
  have hexi : Real.exp 1 * (2 * Real.exp (π / 4)) < 11.92402 := by
    have := mul_lt_mul'' Real.exp_one_lt_d9 hxi (by positivity) (by positivity); linarith
  have he5 : Real.exp (Real.exp 1 * (2 * Real.exp (π / 4))) < 150860 :=
    lt_of_le_of_lt (Q05.exp_le_pow _ 11.92402 12 (by norm_num) (by norm_num) hexi.le
      (by norm_num)) (by norm_num)
  have hB1 : 64 / (π * Q05.a) * Real.exp (Q05.a / 2 + π / 4) < 70.1 := by
    have := mul_lt_mul'' hc1 he1 (by positivity) (by positivity)
    linarith
  have hB2 : 32 / (3 * π * Q05.a) * Real.exp (3 * Q05.a / 2 + 3 * π / 4) < 171 := by
    have := mul_lt_mul'' hc2 he2 (by positivity) (by positivity)
    linarith
  have hS : Real.exp (π / 4) + Real.exp (2 * Real.exp (π / 4)) +
      Real.exp (Real.exp 1 * (2 * Real.exp (π / 4))) < 2.1933 + 80.37 + 150860 := by linarith
  have := mul_lt_mul'' hB1 hS (by positivity) (by positivity)
  linarith

namespace Q05
open Real MeasureTheory Set Filter

/-- The integrand of `g5:lem:kernels (a)`. -/
lemma image_exp_neg_two : (fun u : ℝ => Real.exp (-2 * u)) '' Ioi 0 = Ioo 0 1 := by
  ext x
  simp only [mem_image, mem_Ioi, mem_Ioo]
  constructor
  · rintro ⟨u, hu, rfl⟩
    exact ⟨Real.exp_pos _, Real.exp_lt_one_iff.2 (by linarith)⟩
  · rintro ⟨h0, h1⟩
    refine ⟨-Real.log x / 2, ?_, ?_⟩
    · have := Real.log_neg h0 h1; linarith
    · rw [show -2 * (-Real.log x / 2) = Real.log x by ring, Real.exp_log h0]

lemma beta_integrand_comp (s : ℂ) (u : ℝ) (hu : 0 ≤ u) :
    |(-2 * Real.exp (-2 * u))| • ((Real.exp (-2 * u) : ℂ) ^ (s / 2 - 1) *
      (1 - (Real.exp (-2 * u) : ℂ)) ^ ((3 / 2 : ℂ) - 1)) =
    2 * (((√(1 - Real.exp (-2 * u)) : ℝ) : ℂ) * Complex.exp (-s * u)) := by
  have he : 0 < Real.exp (-2 * u) := Real.exp_pos _
  rw [abs_of_neg (by linarith), Complex.real_smul]
  have h1 : (Real.exp (-2 * u) : ℂ) ^ (s / 2 - 1) = Complex.exp ((-2 * u) * (s / 2 - 1)) := by
    rw [Complex.cpow_def_of_ne_zero (by exact_mod_cast he.ne'), ← Complex.ofReal_log he.le,
      Real.log_exp]
    push_cast; ring_nf
  have hu1 : 0 ≤ 1 - Real.exp (-2 * u) := by
    have := Real.exp_le_one_iff.2 (by linarith : -2 * u ≤ 0); linarith
  have h2 : (1 - (Real.exp (-2 * u) : ℂ)) ^ ((3 / 2 : ℂ) - 1) =
      (((√(1 - Real.exp (-2 * u)) : ℝ) : ℂ)) := by
    rw [Real.sqrt_eq_rpow, Complex.ofReal_cpow hu1]
    push_cast; norm_num
  rw [h1, h2]
  have h3 : (((-2 * Real.exp (-2 * u)) : ℝ) : ℂ) = -2 * Complex.exp (-2 * u) := by push_cast; ring
  push_cast
  have h4 : Complex.exp (-2 * u) * Complex.exp ((-2 * u) * (s / 2 - 1)) = Complex.exp (-s * u) := by
    rw [← Complex.exp_add]; ring_nf
  linear_combination (2 * (((√(1 - Real.exp (-2 * u)) : ℝ) : ℂ))) * h4

/-- The integrand `√(1-e^{-2u}) e^{-su}`. -/
noncomputable def Hint (s : ℂ) (u : ℝ) : ℂ :=
  ((√(1 - Real.exp (-2 * u)) : ℝ) : ℂ) * Complex.exp (-s * u)

lemma H_rep (s : ℂ) (hs : 0 < s.re) :
    IntegrableOn (Hint s) (Ioi 0) ∧ H s = ∫ u in Ioi (0 : ℝ), Hint s u := by
  set g : ℝ → ℂ := fun x => (x : ℂ) ^ (s / 2 - 1) * (1 - (x : ℂ)) ^ ((3 / 2 : ℂ) - 1)
  have hder : ∀ x ∈ Ioi (0 : ℝ), HasDerivWithinAt (fun u : ℝ => Real.exp (-2 * u))
      (-2 * Real.exp (-2 * x)) (Ioi 0) x := by
    intro x _
    have := ((hasDerivAt_id x).const_mul (-2 : ℝ)).exp
    simp only [id, mul_one] at this
    exact (by simpa [mul_comm] using this : HasDerivAt (fun u : ℝ => Real.exp (-2 * u))
      (-2 * Real.exp (-2 * x)) x).hasDerivWithinAt
  have hinj : InjOn (fun u : ℝ => Real.exp (-2 * u)) (Ioi 0) := by
    intro x _ y _ h
    have := Real.exp_injective h
    linarith
  have hcv : IntervalIntegrable g volume 0 1 :=
    Complex.betaIntegral_convergent (by simpa using hs) (by norm_num)
  have hg : IntegrableOn g (Ioo 0 1) :=
    ((intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one).1 hcv).mono_set
      Ioo_subset_Ioc_self
  rw [← image_exp_neg_two] at hg
  have hI := (integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ioi hder hinj g).1 hg
  have hEq := integral_image_eq_integral_abs_deriv_smul measurableSet_Ioi hder hinj g
  rw [image_exp_neg_two] at hEq
  have hc : EqOn (fun x => |(-2 * Real.exp (-2 * x))| • g (Real.exp (-2 * x)))
      (fun u => 2 * Hint s u) (Ioi 0) := by
    intro u hu
    exact beta_integrand_comp s u (le_of_lt hu)
  have hI2 : IntegrableOn (fun u => 2 * Hint s u) (Ioi 0) := hI.congr_fun hc measurableSet_Ioi
  have hI3 : IntegrableOn (Hint s) (Ioi 0) := by
    have : IntegrableOn (fun u => 1 / 2 * (2 * Hint s u)) (Ioi 0) := hI2.const_mul (1 / 2 : ℂ)
    refine this.congr_fun (fun u _ => by ring) measurableSet_Ioi
  refine ⟨hI3, ?_⟩
  have hB : Complex.betaIntegral (s / 2) (3 / 2) = ∫ x in Ioo (0 : ℝ) 1, g x := by
    rw [Complex.betaIntegral, intervalIntegral.integral_of_le zero_le_one,
      integral_Ioc_eq_integral_Ioo]
  rw [H, hB, hEq, setIntegral_congr_fun measurableSet_Ioi hc, integral_const_mul]
  ring

lemma norm_Hint (s : ℂ) (u : ℝ) : ‖Hint s u‖ = √(1 - Real.exp (-2 * u)) * Real.exp (-(s.re * u)) := by
  rw [Hint, norm_mul, Complex.norm_exp, Complex.norm_real, Real.norm_eq_abs,
    abs_of_nonneg (Real.sqrt_nonneg _)]
  congr 2
  simp

lemma H_le_one_div (s : ℂ) (hs : 0 < s.re) : ‖H s‖ ≤ 1 / s.re := by
  obtain ⟨hI, hH⟩ := H_rep s hs
  rw [hH]
  have hg : IntegrableOn (fun u : ℝ => Real.exp (-s.re * u)) (Ioi 0) :=
    exp_neg_integrableOn_Ioi 0 hs
  calc ‖∫ u in Ioi (0 : ℝ), Hint s u‖ ≤ ∫ u in Ioi (0 : ℝ), ‖Hint s u‖ :=
        norm_integral_le_integral_norm _
    _ ≤ ∫ u in Ioi (0 : ℝ), Real.exp (-s.re * u) := by
        refine setIntegral_mono_on hI.norm hg measurableSet_Ioi (fun u hu => ?_)
        rw [norm_Hint, show -s.re * u = -(s.re * u) by ring]
        have h1 : √(1 - Real.exp (-2 * u)) ≤ 1 := by
          rw [Real.sqrt_le_one]; linarith [Real.exp_pos (-2 * u)]
        have := Real.exp_pos (-(s.re * u))
        nlinarith
    _ = 1 / s.re := by
        rw [integral_exp_mul_Ioi (by linarith) 0]; simp

lemma H_le_rpow (s : ℂ) (hs : 0 < s.re) : ‖H s‖ ≤ √(π / 2) * s.re ^ (-(3 / 2) : ℝ) := by
  obtain ⟨hI, hH⟩ := H_rep s hs
  rw [hH]
  set x := s.re
  have hg0 : IntegrableOn (fun u : ℝ => u ^ (1 / 2 : ℝ) * Real.exp (-x * u ^ (1 : ℝ))) (Ioi 0) :=
    integrableOn_rpow_mul_exp_neg_mul_rpow (by norm_num) (by norm_num) hs
  have hg : IntegrableOn (fun u : ℝ => √2 * (u ^ ((3 / 2 : ℝ) - 1) * Real.exp (-(x * u))))
      (Ioi 0) := by
    have hg1 : IntegrableOn (fun u : ℝ => √2 * (u ^ (1 / 2 : ℝ) * Real.exp (-x * u ^ (1 : ℝ))))
        (Ioi 0) := hg0.const_mul √2
    refine hg1.congr_fun (fun u _ => ?_) measurableSet_Ioi
    simp only [Real.rpow_one, neg_mul]; norm_num
  have hG : Real.Gamma (3 / 2) = √π / 2 := by
    rw [show (3 / 2 : ℝ) = 1 / 2 + 1 by norm_num, Real.Gamma_add_one (by norm_num),
      Real.Gamma_one_half_eq]; ring
  calc ‖∫ u in Ioi (0 : ℝ), Hint s u‖ ≤ ∫ u in Ioi (0 : ℝ), ‖Hint s u‖ :=
        norm_integral_le_integral_norm _
    _ ≤ ∫ u in Ioi (0 : ℝ), √2 * (u ^ ((3 / 2 : ℝ) - 1) * Real.exp (-(x * u))) := by
        refine setIntegral_mono_on hI.norm hg measurableSet_Ioi (fun u hu => ?_)
        rw [norm_Hint]
        have hu' : 0 ≤ u := le_of_lt hu
        have h1 : √(1 - Real.exp (-2 * u)) ≤ √2 * u ^ ((3 / 2 : ℝ) - 1) := by
          rw [show (3 / 2 : ℝ) - 1 = 1 / 2 by norm_num, ← Real.sqrt_eq_rpow,
            ← Real.sqrt_mul (by norm_num)]
          apply Real.sqrt_le_sqrt
          linarith [Real.add_one_le_exp (-2 * u)]
        have := Real.exp_pos (-(x * u))
        nlinarith
    _ = √(π / 2) * x ^ (-(3 / 2) : ℝ) := by
        rw [integral_const_mul, Real.integral_rpow_mul_exp_neg_mul_Ioi (by norm_num) hs, hG,
          one_div, Real.inv_rpow hs.le, ← Real.rpow_neg hs.le, Real.sqrt_div' _ (by norm_num)]
        have h2 : (0 : ℝ) < √2 := by positivity
        field_simp
        rw [Real.sq_sqrt (by norm_num)]

lemma H_eq_Gamma (s : ℂ) (hs : 0 < s.re) :
    H s = 1 / 2 * (Complex.Gamma (s / 2) * Complex.Gamma (3 / 2) / Complex.Gamma (s / 2 + 3 / 2)) := by
  rw [H, Complex.betaIntegral_eq_Gamma_mul_div _ _ (by simpa using hs) (by norm_num)]

lemma H_diff : DifferentiableOn ℂ H {s : ℂ | 0 < s.re} := by
  have hop : IsOpen {s : ℂ | 0 < s.re} := isOpen_lt continuous_const Complex.continuous_re
  intro s hs
  have hs' : 0 < s.re := hs
  refine DifferentiableAt.differentiableWithinAt ?_
  have hev : (fun z => 1 / 2 * (Complex.Gamma (z / 2) * Complex.Gamma (3 / 2) /
      Complex.Gamma (z / 2 + 3 / 2))) =ᶠ[nhds s] H := by
    filter_upwards [hop.mem_nhds hs] with z hz
    exact (H_eq_Gamma z hz).symm
  refine (Filter.EventuallyEq.differentiableAt_iff hev).1 ?_
  have h1 : DifferentiableAt ℂ Complex.Gamma (s / 2) := by
    apply Complex.differentiableAt_Gamma
    intro m hm
    have := congrArg Complex.re hm
    simp at this
    have : (0 : ℝ) ≤ m := Nat.cast_nonneg m
    linarith
  have h2 : DifferentiableAt ℂ Complex.Gamma (s / 2 + 3 / 2) := by
    apply Complex.differentiableAt_Gamma
    intro m hm
    have := congrArg Complex.re hm
    simp at this
    have : (0 : ℝ) ≤ m := Nat.cast_nonneg m
    norm_num at this
    linarith
  have h3 : Complex.Gamma (s / 2 + 3 / 2) ≠ 0 :=
    Complex.Gamma_ne_zero_of_re_pos (by simp; linarith)
  have hc1 : DifferentiableAt ℂ (fun z : ℂ => Complex.Gamma (z / 2)) s :=
    h1.comp s (by fun_prop)
  have hc2 : DifferentiableAt ℂ (fun z : ℂ => Complex.Gamma (z / 2 + 3 / 2)) s :=
    h2.comp s (by fun_prop)
  exact ((hc1.mul_const _).div hc2 h3).const_mul _

lemma H_odd (n : ℕ) : H (2 * n + 1) =
    ((√π / 4 * Real.Gamma (n + 1 / 2) / (Nat.factorial (n + 1) : ℝ) : ℝ) : ℂ) := by
  rw [H_eq_Gamma _ (by simp; positivity)]
  have e1 : ((2 * n + 1 : ℂ)) / 2 = ((n + 1 / 2 : ℝ) : ℂ) := by push_cast; ring
  have e2 : ((2 * n + 1 : ℂ)) / 2 + 3 / 2 = (((n + 1 : ℕ) : ℝ) + 1 : ℝ) := by push_cast; ring
  have e3 : (3 / 2 : ℂ) = ((1 / 2 + 1 : ℝ) : ℂ) := by push_cast; norm_num
  rw [e2, e1, e3, Complex.Gamma_ofReal, Complex.Gamma_ofReal, Complex.Gamma_ofReal,
    Real.Gamma_nat_eq_factorial, Real.Gamma_add_one (by norm_num), Real.Gamma_one_half_eq]
  push_cast
  ring

end Q05

theorem Q05.proof_5 : Q05.stmt_5 :=
  ⟨fun s hs => Q05.H_rep s hs, Q05.H_diff, fun s hs => ⟨Q05.H_le_one_div s hs, Q05.H_le_rpow s hs⟩,
    Q05.H_odd⟩

namespace Q05
open Filter

lemma cos_ne_zero_of (s : ℂ) (hs : ∀ n : ℤ, s ≠ 2 * n + 1) :
    Complex.cos (π * s / 2) ≠ 0 := by
  intro h
  obtain ⟨k, hk⟩ := Complex.cos_eq_zero_iff.1 h
  apply hs k
  have hp : (π : ℂ) ≠ 0 := Complex.ofReal_ne_zero.2 Real.pi_ne_zero
  have : (π : ℂ) * s = π * (2 * k + 1) := by linear_combination 2 * hk
  exact mul_left_cancel₀ hp this

lemma sin_odd (n : ℤ) : Complex.sin (π * (2 * n + 1 : ℂ) / 2) = (-1 : ℂ) ^ n := by
  have : (π * (2 * n + 1 : ℂ) / 2) = ((π / 2 + n * π : ℝ) : ℂ) := by push_cast; ring
  rw [this, ← Complex.ofReal_sin, Real.sin_add_int_mul_pi, Real.sin_pi_div_two, mul_one]
  push_cast; rfl

lemma K_residue (n : ℤ) : Tendsto (fun s : ℂ => (s - (2 * n + 1)) * K s)
      (nhdsWithin ((2 * n + 1 : ℂ)) {(2 * n + 1 : ℂ)}ᶜ) (nhds ((-1 : ℂ) ^ (n + 1) / π)) := by
  have hf : HasDerivAt (fun s : ℂ => Complex.cos (π * s / 2))
      (-Complex.sin (π * (2 * n + 1 : ℂ) / 2) * (π / 2)) (2 * n + 1 : ℂ) := by
    have h1 : HasDerivAt (fun s : ℂ => π * s / 2) (π / 2) (2 * n + 1 : ℂ) := by
      simpa using ((hasDerivAt_id (2 * n + 1 : ℂ)).const_mul (π : ℂ)).div_const 2
    exact (Complex.hasDerivAt_cos _).comp (2 * n + 1 : ℂ) h1
  have hf0 : Complex.cos (π * (2 * n + 1 : ℂ) / 2) = 0 := by
    rw [Complex.cos_eq_zero_iff]; exact ⟨n, by ring⟩
  have hsl := hasDerivAt_iff_tendsto_slope.1 hf
  have hsin := sin_odd n
  have hp : (π : ℂ) ≠ 0 := Complex.ofReal_ne_zero.2 Real.pi_ne_zero
  have hm : (-1 : ℂ) ^ n ≠ 0 := zpow_ne_zero _ (by norm_num)
  rw [hsin] at hsl
  have hne : (2 * (-(-1 : ℂ) ^ n * (π / 2))) ≠ 0 := by simp [hp, hm]
  have ht := (hsl.const_mul (2 : ℂ)).inv₀ hne
  have hlim : (2 * (-(-1 : ℂ) ^ n * (π / 2)))⁻¹ = (-1 : ℂ) ^ (n + 1) / π := by
    rw [zpow_add_one₀ (by norm_num : (-1 : ℂ) ≠ 0)]
    field_simp
    have h2 : ((-1 : ℂ) ^ n) ^ 2 = 1 := by
      rw [← zpow_natCast, ← zpow_mul, mul_comm, zpow_mul]; norm_num
    linear_combination h2
  rw [hlim] at ht
  refine ht.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with s hs
  have hs' : s - (2 * n + 1) ≠ 0 := sub_ne_zero.2 hs
  simp only [slope_def_field, hf0, sub_zero, K]
  field_simp

end Q05

theorem Q05.proof_6 : Q05.stmt_6 := by
  unfold Q05.stmt_6
  refine ⟨fun s hs => ?_, Q05.K_residue⟩
  · have h := Q05.cos_ne_zero_of s hs
    have h2 : (2 : ℂ) * Complex.cos (π * s / 2) ≠ 0 := mul_ne_zero two_ne_zero h
    show DifferentiableAt ℂ (fun s : ℂ => 1 / (2 * Complex.cos (π * s / 2))) s
    fun_prop (disch := exact h2)

namespace Q05
open Real MeasureTheory Set Filter

lemma sqrt2_sq : √2 * √2 = 2 := Real.mul_self_sqrt (by norm_num)

lemma ω_eq : ω = ((√2 / 2 : ℝ) : ℂ) + ((√2 / 2 : ℝ) : ℂ) * Complex.I := by
  rw [ω, show (π : ℂ) / 4 = ((π / 4 : ℝ) : ℂ) by push_cast; ring, Complex.exp_mul_I, ← Complex.ofReal_cos, ← Complex.ofReal_sin,
    Real.cos_pi_div_four, Real.sin_pi_div_four]

lemma Lpt_eq (r : ℝ) : Lpt r = ((1 / 2 + r * (√2 / 2) : ℝ) : ℂ) + ((r * (√2 / 2) : ℝ) : ℂ) * Complex.I := by
  rw [Lpt, ω_eq]; push_cast; ring

lemma Lpt_re (r : ℝ) : (Lpt r).re = 1 / 2 + r * (√2 / 2) := by
  rw [Lpt_eq]; simp

lemma a_eq : a = π / 2 * (√2 / 2) := by
  rw [a]; have := sqrt2_sq; have h2 : (0:ℝ) < √2 := by positivity
  field_simp; nlinarith

lemma w_eq (r : ℝ) : (π : ℂ) * Lpt r / 2 = ((π / 4 + a * r : ℝ) : ℂ) + ((a * r : ℝ) : ℂ) * Complex.I := by
  rw [Lpt_eq, a_eq]; push_cast; ring

lemma key_ineq (b : ℝ) (hb : 0 ≤ b) :
    3 * Real.exp (2 * b) ^ 2 + 4 ≥ 8 * Real.exp (2 * b) * Real.sin (2 * b) := by
  set Y := Real.exp (2 * b)
  have hY : 1 + 2 * b ≤ Y := by have := Real.add_one_le_exp (2 * b); linarith
  have hY0 : 0 < Y := Real.exp_pos _
  have hs1 : Real.sin (2 * b) ≤ 1 := Real.sin_le_one _
  have hs2 : Real.sin (2 * b) ≤ 2 * b := Real.sin_le (by linarith)
  rcases le_total Y 2 with h | h
  · nlinarith [mul_le_mul_of_nonneg_left hs2 hY0.le]
  · nlinarith [mul_le_mul_of_nonneg_left hs1 hY0.le]

lemma cos_lower (b : ℝ) (hb : 0 ≤ b) :
    Real.exp b / 4 ≤ ‖Complex.cos (((π / 4 + b : ℝ) : ℂ) + ((b : ℝ) : ℂ) * Complex.I)‖ := by
  rw [Complex.cos_add_mul_I, ← Complex.ofReal_cos, ← Complex.ofReal_sin, ← Complex.ofReal_cosh,
    ← Complex.ofReal_sinh]
  rw [Complex.norm_eq_sqrt_sq_add_sq]
  simp only [Complex.sub_re, Complex.sub_im, Complex.mul_re, Complex.mul_im, Complex.ofReal_re,
    Complex.ofReal_im, Complex.I_re, Complex.I_im]
  apply Real.le_sqrt_of_sq_le
  have hc : Real.cos (π / 4 + b) = √2 / 2 * (Real.cos b - Real.sin b) := by
    rw [Real.cos_add, Real.cos_pi_div_four, Real.sin_pi_div_four]; ring
  have hsn : Real.sin (π / 4 + b) = √2 / 2 * (Real.cos b + Real.sin b) := by
    rw [Real.sin_add, Real.cos_pi_div_four, Real.sin_pi_div_four]; ring
  rw [hc, hsn]
  have h2 := sqrt2_sq
  have hch := Real.cosh_sq b
  have hpy := Real.sin_sq_add_cos_sq b
  have hsh : Real.sinh b = (Real.exp b - Real.exp (-b)) / 2 := Real.sinh_eq b
  have hk := key_ineq b hb
  have hs2 : Real.sin (2 * b) = 2 * Real.sin b * Real.cos b := Real.sin_two_mul b
  have he2 : Real.exp (2 * b) = Real.exp b ^ 2 := by rw [← Real.exp_nat_mul]; ring_nf
  have hen : Real.exp (-b) * Real.exp b = 1 := by rw [← Real.exp_add]; simp
  rw [hs2, he2] at hk
  -- target: (e^b/4)^2 ≤ (√2/2 (c-s) cosh)^2 + (-(√2/2 (c+s)) sinh)^2 (in some arrangement)
  have hcs : (√2 / 2) ^ 2 = 1 / 2 := by nlinarith
  have hE : 0 < Real.exp b := Real.exp_pos b
  -- normalise
  have goal : (Real.exp b / 4) ^ 2 ≤
      (√2 / 2 * (Real.cos b - Real.sin b) * Real.cosh b) ^ 2 +
        (√2 / 2 * (Real.cos b + Real.sin b) * Real.sinh b) ^ 2 := by
    have e1 : (√2 / 2 * (Real.cos b - Real.sin b) * Real.cosh b) ^ 2 +
        (√2 / 2 * (Real.cos b + Real.sin b) * Real.sinh b) ^ 2 =
        1 / 2 * ((Real.cos b - Real.sin b) ^ 2 * Real.cosh b ^ 2 +
          (Real.cos b + Real.sin b) ^ 2 * Real.sinh b ^ 2) := by
      rw [mul_pow, mul_pow, mul_pow, mul_pow, hcs]; ring
    rw [e1, hch]
    have e2 : (Real.cos b - Real.sin b) ^ 2 * (Real.sinh b ^ 2 + 1) +
        (Real.cos b + Real.sin b) ^ 2 * Real.sinh b ^ 2 =
        (1 - 2 * Real.sin b * Real.cos b) + 2 * Real.sinh b ^ 2 := by
      nlinarith [hpy]
    rw [e2, hsh]
    set X := Real.exp b
    set Z := Real.exp (-b)
    set P := 2 * Real.sin b * Real.cos b
    -- (X/4)^2 ≤ 1/2 (1 - P + (X - Z)^2/2), with XZ = 1, 3X^4+4 ≥ 8X^2 P
    have hZ : Z = 1 / X := by field_simp; linarith
    rw [hZ]
    rw [div_pow, show (1 : ℝ) / 2 * (1 - P + 2 * ((X - 1 / X) / 2) ^ 2) =
      (2 * X ^ 2 - 2 * X ^ 2 * P + X ^ 4 - 2 * X ^ 2 + 1) / (4 * X ^ 2) by field_simp; ring]
    rw [div_le_div_iff₀ (by norm_num) (by positivity)]
    nlinarith
  linarith [goal, sq_nonneg (Real.cos b), neg_sq (√2 / 2 * (Real.cos b + Real.sin b) * Real.sinh b)]


lemma a_pos : 0 < a := by rw [a]; positivity

lemma stmt7_at (r : ℝ) (hr : 0 ≤ r) :
    ‖H (Lpt r)‖ ≤ 2 ∧ ‖E (Lpt r)‖ = Real.exp (-a * r) ∧ ‖K (Lpt r)‖ ≤ 2 * Real.exp (-a * r) ∧
      ‖K (Lpt r) - E (Lpt r)‖ ≤ 2 * Real.exp (-3 * a * r) := by
  have hre := Lpt_re r
  have hs2 : 0 ≤ r * (√2 / 2) := by positivity
  have hb : 0 ≤ a * r := mul_nonneg a_pos.le hr
  set w : ℂ := ((π / 4 + a * r : ℝ) : ℂ) + ((a * r : ℝ) : ℂ) * Complex.I with hw
  have hwL : (π : ℂ) * Lpt r / 2 = w := w_eq r
  have hEw : E (Lpt r) = Complex.exp (w * Complex.I) := by
    rw [E, ← hwL]; ring_nf
  have hKw : K (Lpt r) = 1 / (2 * Complex.cos w) := by rw [K, hwL]
  have hcl := cos_lower (a * r) hb
  have hcpos : 0 < ‖Complex.cos w‖ := lt_of_lt_of_le (by positivity) hcl
  have hnE : ‖Complex.exp (w * Complex.I)‖ = Real.exp (-a * r) := by
    rw [Complex.norm_exp, hw]; congr 1; simp
  have hnK : ‖K (Lpt r)‖ ≤ 2 * Real.exp (-a * r) := by
    rw [hKw, norm_div, norm_one, norm_mul, Complex.norm_two]
    rw [div_le_iff₀ (by positivity)]
    have : Real.exp (-a * r) * Real.exp (a * r) = 1 := by rw [← Real.exp_add]; simp
    nlinarith [Real.exp_pos (-a * r)]
  refine ⟨?_, by rw [hEw, hnE], hnK, ?_⟩
  · have h1 := H_le_one_div (Lpt r) (by rw [hre]; linarith)
    refine h1.trans ?_
    rw [hre, div_le_iff₀ (by linarith)]; linarith
  · have hc0 : Complex.cos w ≠ 0 := norm_pos_iff.1 hcpos
    have h2c : 2 * Complex.cos w = Complex.exp (w * Complex.I) + Complex.exp (-w * Complex.I) :=
      Complex.two_cos w
    have hpq : Complex.exp (w * Complex.I) * Complex.exp (-w * Complex.I) = 1 := by
      rw [← Complex.exp_add]; simp
    have hdiff : K (Lpt r) - E (Lpt r) = -(Complex.exp (w * Complex.I)) ^ 2 * K (Lpt r) := by
      rw [hKw, hEw]
      have h2 : (2 : ℂ) * Complex.cos w ≠ 0 := mul_ne_zero two_ne_zero hc0
      field_simp
      linear_combination (-(Complex.exp (w * Complex.I))) * h2c - hpq
    rw [hdiff, norm_mul, norm_neg, norm_pow, hnE]
    have he : Real.exp (-a * r) ^ 2 * (2 * Real.exp (-a * r)) = 2 * Real.exp (-3 * a * r) := by
      rw [← Real.exp_nat_mul, mul_left_comm, ← Real.exp_add]; ring_nf
    rw [← he]
    exact mul_le_mul_of_nonneg_left hnK (by positivity)

end Q05

theorem Q05.proof_7 : Q05.stmt_7 := Q05.stmt7_at

namespace Q05
open Real MeasureTheory Set Filter

lemma Lpt_im (r : ℝ) : (Lpt r).im = r * (√2 / 2) := by
  rw [Lpt_eq]; simp

lemma norm_Lpt_le (r : ℝ) (hr : 0 ≤ r) : ‖Lpt r‖ ≤ 1 / 2 + r := by
  have hω : ‖ω‖ = 1 := by
    rw [ω, Complex.norm_exp]; simp
  calc ‖Lpt r‖ ≤ ‖(1 / 2 : ℂ)‖ + ‖(r : ℂ) * ω‖ := norm_add_le _ _
    _ = 1 / 2 + r := by
      rw [norm_mul, hω, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hr]; norm_num

lemma two_sqrt (X : ℝ) : 2 * √X = √(4 * X) := by
  rw [Real.sqrt_mul (by norm_num), show (4 : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]

lemma amgm (A B t : ℝ) (hA : 0 ≤ A) (hB : 0 ≤ B) (ht : 0 < t) :
    2 * √(A * B) ≤ A * t + B / t := by
  have e : √(A * B) = √(A * t) * √(B / t) := by
    rw [← Real.sqrt_mul (by positivity)]; congr 1; field_simp
  rw [e]
  have h1 := Real.sq_sqrt (by positivity : 0 ≤ A * t)
  have h2 := Real.sq_sqrt (by positivity : 0 ≤ B / t)
  nlinarith [sq_nonneg (√(A * t) - √(B / t))]

lemma Fjq_pointwise (j : ℕ) (q l r : ℝ) (hq : 1 ≤ q) (hl : 0 < l) (hr : 0 ≤ r) :
    FjqIntegrand j q l r ≤ Real.exp (-(l / 4) * r) *
      (Real.exp (3 * l / 4 * (1 / 2 + q / √2)) * Real.exp (-√(3 * l * j * q / √2))) := by
  set ρ := ‖Lpt r‖
  set d := ‖Lpt r + q‖
  set c := q / √2
  have hs2 : (0 : ℝ) < √2 := by positivity
  have h22 := sqrt2_sq
  have hc : 0 < c := by positivity
  have hρ0 : 0 ≤ ρ := norm_nonneg _
  have hre := Lpt_re r
  have him := Lpt_im r
  -- ρ^2 = re^2 + im^2, d^2 = (re+q)^2 + im^2
  have hρ2 : ρ ^ 2 = (Lpt r).re ^ 2 + (Lpt r).im ^ 2 := by
    simp only [ρ]; rw [Complex.norm_eq_sqrt_sq_add_sq, Real.sq_sqrt (by positivity)]
  have hd2 : d ^ 2 = ((Lpt r).re + q) ^ 2 + (Lpt r).im ^ 2 := by
    simp only [d]
    rw [Complex.norm_eq_sqrt_sq_add_sq, Real.sq_sqrt (by positivity)]; simp
  have hd0 : 0 ≤ d := norm_nonneg _
  -- re ≥ ρ/√2
  have hreρ : ρ ≤ √2 * (Lpt r).re := by
    have h1 : 0 ≤ √2 * (Lpt r).re := by rw [hre]; positivity
    have h2 : ρ ^ 2 ≤ (√2 * (Lpt r).re) ^ 2 := by
      rw [hρ2, mul_pow, Real.sq_sqrt (by norm_num), hre, him]; nlinarith
    exact (pow_le_pow_iff_left₀ hρ0 h1 (by norm_num)).1 h2
  have hdρ : ρ + c ≤ d := by
    have h2 : (ρ + c) ^ 2 ≤ d ^ 2 := by
      rw [hd2]
      have hcc : c ^ 2 = q ^ 2 / 2 := by simp only [c]; rw [div_pow, Real.sq_sqrt (by norm_num)]
      have : 2 * c * ρ ≤ 2 * q * (Lpt r).re := by
        simp only [c]
        rw [show 2 * (q / √2) * ρ = 2 * q * (ρ / √2) by ring]
        have : ρ / √2 ≤ (Lpt r).re := by rw [div_le_iff₀ hs2]; linarith
        nlinarith
      nlinarith [hρ2]
    nlinarith [sq_nonneg (ρ + c - d)]
  have hd1 : 1 ≤ d := by
    have : q ^ 2 ≤ d ^ 2 := by
      have hp : 0 ≤ r * (√2 / 2) := by positivity
      rw [hd2, hre]; nlinarith [sq_nonneg (Lpt r).im]
    nlinarith
  have hdpos : 0 < d := by linarith
  -- step 1: ρ^j / d^(j+1) ≤ (ρ/(ρ+c))^j
  have st1 : ρ ^ j / d ^ (j + 1) ≤ (ρ / (ρ + c)) ^ j := by
    rw [div_pow, pow_succ]
    rw [div_le_div_iff₀ (by positivity) (by positivity)]
    have : (ρ + c) ^ j ≤ d ^ j := pow_le_pow_left₀ (by positivity) hdρ j
    have : d ^ j ≤ d ^ j * d := le_mul_of_one_le_right (by positivity) hd1
    exact mul_le_mul_of_nonneg_left (by linarith) (by positivity)
  -- step 2: ρ/(ρ+c) ≤ exp(-c/(ρ+c))
  have st2 : (ρ / (ρ + c)) ^ j ≤ Real.exp (-(j * c / (1 / 2 + r + c))) := by
    have h1 : ρ / (ρ + c) ≤ Real.exp (-(c / (ρ + c))) := by
      have := Real.add_one_le_exp (-(c / (ρ + c)))
      have e : ρ / (ρ + c) = -(c / (ρ + c)) + 1 := by field_simp; ring
      linarith
    calc (ρ / (ρ + c)) ^ j ≤ Real.exp (-(c / (ρ + c))) ^ j :=
          pow_le_pow_left₀ (by positivity) h1 j
      _ = Real.exp (-(j * (c / (ρ + c)))) := by rw [← Real.exp_nat_mul]; ring_nf
      _ ≤ Real.exp (-(j * c / (1 / 2 + r + c))) := by
          apply Real.exp_le_exp.2
          have hρle := norm_Lpt_le r hr
          have : c / (1 / 2 + r + c) ≤ c / (ρ + c) :=
            div_le_div_of_nonneg_left hc.le (by positivity) (by linarith)
          have hj : (0 : ℝ) ≤ j := Nat.cast_nonneg j
          rw [mul_div_assoc]
          nlinarith
  -- step 3: AM-GM
  set t := 1 / 2 + r + c
  have ht : 0 < t := by positivity
  have st3 : 2 * √(3 * l / 4 * (j * c)) ≤ 3 * l / 4 * t + j * c / t :=
    amgm _ _ t (by positivity) (by positivity) ht
  have hsq : 2 * √(3 * l / 4 * (j * c)) = √(3 * l * j * q / √2) := by
    rw [two_sqrt]
    congr 1; simp only [c]; ring
  rw [hsq] at st3
  unfold FjqIntegrand
  have hmain : ρ ^ j / d ^ (j + 1) ≤ Real.exp (-(j * c / t)) := st1.trans st2
  calc Real.exp (-l * r) * ρ ^ j / d ^ (j + 1) = Real.exp (-l * r) * (ρ ^ j / d ^ (j + 1)) := by
        ring
    _ ≤ Real.exp (-l * r) * Real.exp (-(j * c / t)) :=
        mul_le_mul_of_nonneg_left hmain (Real.exp_pos _).le
    _ = Real.exp (-(l / 4) * r) * Real.exp (3 * l / 4 * (1 / 2 + q / √2) -
          (3 * l / 4 * t + j * c / t)) := by
        rw [← Real.exp_add, ← Real.exp_add]; congr 1; simp only [t, c]; ring
    _ ≤ _ := by
        simp only [← Real.exp_add]
        apply Real.exp_le_exp.2
        linarith

lemma Fjq_meas (j : ℕ) (q l : ℝ) : Measurable (FjqIntegrand j q l) := by
  unfold FjqIntegrand Lpt
  fun_prop

end Q05

theorem Q05.proof_13 : Q05.stmt_13 := by
  intro j q l hq hl
  set C := Real.exp (3 * l / 4 * (1 / 2 + q / √2)) * Real.exp (-√(3 * l * j * q / √2))
  have hg : IntegrableOn (fun r : ℝ => Real.exp (-(l / 4) * r) * C) (Ioi 0) :=
    (exp_neg_integrableOn_Ioi 0 (by positivity : 0 < l / 4)).mul_const C
  have hpt : ∀ r ∈ Ioi (0 : ℝ), Q05.FjqIntegrand j q l r ≤ Real.exp (-(l / 4) * r) * C :=
    fun r hr => Q05.Fjq_pointwise j q l r hq hl (le_of_lt hr)
  have hnn : ∀ r, 0 ≤ Q05.FjqIntegrand j q l r := fun r => by
    unfold Q05.FjqIntegrand; positivity
  have hI : IntegrableOn (Q05.FjqIntegrand j q l) (Ioi 0) := by
    refine Integrable.mono' hg (Q05.Fjq_meas j q l).aestronglyMeasurable ?_
    refine (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall fun r hr => ?_)
    rw [Real.norm_eq_abs, abs_of_nonneg (hnn r)]
    exact hpt r hr
  refine ⟨hI, ?_⟩
  unfold Q05.Fjq
  calc ∫ r in Ioi (0 : ℝ), Q05.FjqIntegrand j q l r
        ≤ ∫ r in Ioi (0 : ℝ), Real.exp (-(l / 4) * r) * C :=
          setIntegral_mono_on hI hg measurableSet_Ioi hpt
    _ = 4 / l * C := by
          rw [integral_mul_const, integral_exp_mul_Ioi (by linarith) 0]
          simp
    _ ≤ 4 / l * Real.exp (l / 2 + l * q / √2) * Real.exp (-√(3 * l * j * q / √2)) := by
          rw [mul_assoc]
          apply mul_le_mul_of_nonneg_left _ (by positivity)
          apply mul_le_mul_of_nonneg_right _ (Real.exp_pos _).le
          apply Real.exp_le_exp.2
          have : 0 ≤ q / √2 := by positivity
          have e : l / 2 + l * q / √2 = l * (1 / 2 + q / √2) := by ring
          rw [e]; nlinarith

namespace Q05
open Real MeasureTheory Set Filter

lemma logit_pt (ε x : ℝ) :
    |(-(Real.exp x / (1 + Real.exp x) ^ 2))| * W02.gR ε (1 / (1 + Real.exp x)) =
      Real.exp (-x) * (1 + Real.exp (-x)) ^ (-1 - ε) * Cf ((1 - ε) * x) := by
  set E := Real.exp x with hE
  set L := 1 + E with hLdef
  have hE0 : 0 < E := Real.exp_pos x
  have hL : 0 < L := by positivity
  set A := L ^ (1 - ε) with hA
  set w := Real.exp (-(1 - ε) * x) with hw
  have hA0 : 0 < A := Real.rpow_pos_of_pos hL _
  have hw0 : 0 < w := Real.exp_pos _
  have h1 : W02.aR ε (1 / L) = A := by
    rw [W02.aR, one_div, Real.inv_rpow hL.le, ← Real.rpow_neg hL.le, hA]; ring_nf
  have hEw : E ^ (ε - 1) = w := by
    rw [hE, ← Real.exp_mul, hw]; ring_nf
  have h2 : W02.bR ε (1 / L) = w * A := by
    rw [W02.bR, show 1 - 1 / L = E / L by field_simp; ring, Real.div_rpow hE0.le hL.le, hEw, hA,
      show (1 : ℝ) - ε = -(ε - 1) by ring, Real.rpow_neg hL.le]
    field_simp
  have h3 : (1 + Real.exp (-x)) ^ (-1 - ε) * Real.exp (-x) = A * w * E / L ^ 2 := by
    have e1 : 1 + Real.exp (-x) = L / E := by rw [Real.exp_neg, ← hE]; field_simp; ring
    rw [e1, Real.div_rpow hL.le hE0.le, show (-1 - ε) = (1 - ε) - 2 by ring,
      Real.rpow_sub hL, ← hA, Real.rpow_two]
    rw [hE, ← Real.exp_mul, hw]
    have : Real.exp (-x) = Real.exp (-(1 - ε) * x) * Real.exp x * Real.exp (x * ((1 - ε) - 2)) := by
      rw [← Real.exp_add, ← Real.exp_add]; ring_nf
    rw [this]
    have := Real.exp_pos (x * ((1 - ε) - 2))
    field_simp
  have hS : √(A ^ 2 + (w * A) ^ 2) = A * √(1 + w ^ 2) := by
    rw [show A ^ 2 + (w * A) ^ 2 = A ^ 2 * (1 + w ^ 2) by ring, Real.sqrt_mul (by positivity),
      Real.sqrt_sq hA0.le]
  have hC : Cf ((1 - ε) * x) = 1 - w / (1 + √(1 + w ^ 2)) := by
    have e1 : Real.exp (-((1 - ε) * x)) = w := by rw [hw]; ring_nf
    have e2 : Real.exp (-2 * ((1 - ε) * x)) = w ^ 2 := by rw [hw, ← Real.exp_nat_mul]; ring_nf
    rw [Cf, e1, e2]
  have hgR : W02.gR ε (1 / L) = A * (1 + w - √(1 + w ^ 2)) := by
    rw [W02.gR, h1, h2, hS]; ring
  rw [hgR, hC, abs_neg, abs_of_pos (by positivity)]
  rw [show Real.exp (-x) * (1 + Real.exp (-x)) ^ (-1 - ε) =
    (1 + Real.exp (-x)) ^ (-1 - ε) * Real.exp (-x) by ring, h3]
  set S := √(1 + w ^ 2)
  have hS2 : S ^ 2 = 1 + w ^ 2 := Real.sq_sqrt (by positivity)
  have hS0 : 0 ≤ S := Real.sqrt_nonneg _
  have h1S : 0 < 1 + S := by linarith
  field_simp
  linear_combination (-1 : ℝ) * hS2

lemma image_logit : (fun x : ℝ => 1 / (1 + Real.exp x)) '' Ioi 0 = Ioo 0 (1 / 2) := by
  ext t
  simp only [mem_image, mem_Ioi, mem_Ioo]
  constructor
  · rintro ⟨x, hx, rfl⟩
    have h1 : 1 < Real.exp x := Real.one_lt_exp_iff.2 hx
    refine ⟨by positivity, ?_⟩
    rw [div_lt_div_iff₀ (by positivity) (by norm_num)]; linarith
  · rintro ⟨h0, h1⟩
    have hv : 1 < 1 / t - 1 := by
      rw [lt_sub_iff_add_lt, lt_div_iff₀ h0]; linarith
    refine ⟨Real.log (1 / t - 1), Real.log_pos hv, ?_⟩
    rw [Real.exp_log (by linarith)]
    field_simp; ring

lemma logit_deriv (x : ℝ) : HasDerivAt (fun x : ℝ => 1 / (1 + Real.exp x))
    (-(Real.exp x / (1 + Real.exp x) ^ 2)) x := by
  have h := ((Real.hasDerivAt_exp x).const_add 1).inv (by positivity)
  simp only [one_div]
  convert h using 1
  field_simp

lemma logit_inj : InjOn (fun x : ℝ => 1 / (1 + Real.exp x)) (Ioi 0) := by
  intro x _ y _ h
  simp only at h
  have h' : 1 + Real.exp x = 1 + Real.exp y := by
    rw [div_eq_div_iff (by positivity) (by positivity)] at h
    linarith
  exact Real.exp_injective (by linarith)

end Q05

theorem Q05.proof_8 : Q05.stmt_8 := by
  intro ε hε
  have hder : ∀ x ∈ Ioi (0 : ℝ), HasDerivWithinAt (fun x : ℝ => 1 / (1 + Real.exp x))
      (-(Real.exp x / (1 + Real.exp x) ^ 2)) (Ioi 0) x :=
    fun x _ => (Q05.logit_deriv x).hasDerivWithinAt
  have hi := W02.gR_iInt hε
  have hg : IntegrableOn (W02.gR ε) (Ioo 0 (1 / 2)) := by
    have := ((intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one).1 hi)
    exact this.mono_set (Ioo_subset_Ioc_self.trans (Ioc_subset_Ioc_right (by norm_num)))
  rw [← Q05.image_logit] at hg
  have hI := (integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ioi hder
    Q05.logit_inj (W02.gR ε)).1 hg
  have hEq := integral_image_eq_integral_abs_deriv_smul measurableSet_Ioi hder Q05.logit_inj
    (W02.gR ε)
  rw [Q05.image_logit] at hEq
  have hc : EqOn (fun x => |(-(Real.exp x / (1 + Real.exp x) ^ 2))| • W02.gR ε (1 / (1 + Real.exp x)))
      (fun x => Real.exp (-x) * (1 + Real.exp (-x)) ^ (-1 - ε) * Q05.Cf ((1 - ε) * x)) (Ioi 0) :=
    fun x _ => by simpa [smul_eq_mul] using Q05.logit_pt ε x
  refine ⟨hI.congr_fun hc measurableSet_Ioi, ?_⟩
  rw [W02.FR_half hε, intervalIntegral.integral_of_le (by norm_num),
    integral_Ioc_eq_integral_Ioo, hEq, setIntegral_congr_fun measurableSet_Ioi hc]
  ring

namespace Q05
open Real MeasureTheory Set Filter

/-- The integrand `H E Q_k` of g5:lem:main. -/
noncomputable def fM (k : ℕ) (s : ℂ) : ℂ := H s * E s * Q k s

lemma one_add_ne (s : ℂ) (hs : 0 ≤ s.re) : 1 + s ≠ 0 := by
  intro h; have := congrArg Complex.re h; simp at this; linarith

lemma norm_le_norm_one_add (s : ℂ) (hs : 0 ≤ s.re) : ‖s‖ ≤ ‖1 + s‖ := by
  rw [Complex.norm_eq_sqrt_sq_add_sq, Complex.norm_eq_sqrt_sq_add_sq]
  apply Real.sqrt_le_sqrt
  simp; nlinarith

lemma one_le_norm_one_add (s : ℂ) : 1 + s.re ≤ ‖1 + s‖ := by
  have := Complex.re_le_norm (1 + s); simpa using this

lemma Q_le (k : ℕ) (s : ℂ) (hs : 0 ≤ s.re) : ‖Q k s‖ ≤ 1 / ‖1 + s‖ := by
  have h0 : 0 < ‖1 + s‖ := norm_pos_iff.2 (one_add_ne s hs)
  rw [Q, norm_div, norm_pow, norm_pow, pow_succ, div_le_div_iff₀ (by positivity) h0]
  have := pow_le_pow_left₀ (norm_nonneg s) (norm_le_norm_one_add s hs) k
  nlinarith [pow_nonneg h0.le k]

lemma E_norm (s : ℂ) : ‖E s‖ = Real.exp (-(π * s.im / 2)) := by
  rw [E, Complex.norm_exp]; congr 1; simp; ring

lemma fM_diff (k : ℕ) : DifferentiableOn ℂ (fM k) {s : ℂ | 0 < s.re} := by
  intro s hs
  have hs' : 0 < s.re := hs
  have h1 : DifferentiableWithinAt ℂ H {s : ℂ | 0 < s.re} s := H_diff s hs
  have h2 : DifferentiableAt ℂ E s := by unfold E; fun_prop
  have h3 : DifferentiableAt ℂ (Q k) s := by
    unfold Q
    have : (1 + s) ^ (k + 1) ≠ 0 := pow_ne_zero _ (one_add_ne s hs'.le)
    fun_prop (disch := exact this)
  exact (h1.mul h2.differentiableWithinAt).mul h3.differentiableWithinAt

lemma fM_cont (k : ℕ) : ContinuousOn (fM k) {s : ℂ | 0 < s.re} := (fM_diff k).continuousOn

lemma Q_le' (k : ℕ) (s : ℂ) : ‖Q k s‖ = ‖s‖ ^ k / ‖1 + s‖ ^ (k + 1) := by
  rw [Q, norm_div, norm_pow, norm_pow]

lemma fM_real_int (k : ℕ) (hk : 1 ≤ k) : IntegrableOn (fun x : ℝ => fM k x) (Ioi 0) := by
  have hc : ContinuousOn (fun x : ℝ => fM k x) (Ioi 0) :=
    (fM_cont k).comp Complex.continuous_ofReal.continuousOn (fun x hx => by simpa using hx)
  rw [← Ioc_union_Ioi_eq_Ioi zero_le_one]
  refine IntegrableOn.union ?_ ?_
  · have hc1 : ContinuousOn (fun x : ℝ => fM k x) (Ioc 0 1) := hc.mono Ioc_subset_Ioi_self
    refine Integrable.mono' (integrableOn_const (C := (1 : ℝ)) (by simp))
      (hc1.aestronglyMeasurable measurableSet_Ioc) ?_
    refine (ae_restrict_iff' measurableSet_Ioc).2 (Filter.Eventually.of_forall fun x hx => ?_)
    have hx0 : 0 < x := hx.1
    have hH := H_le_one_div (x : ℂ) (by simpa using hx0)
    simp only [Complex.ofReal_re] at hH
    have hE : ‖E (x : ℂ)‖ = 1 := by rw [E_norm]; simp
    have hQ : ‖Q k (x : ℂ)‖ ≤ x := by
      rw [Q_le', Complex.norm_real, Real.norm_eq_abs, abs_of_pos hx0]
      have h1 : 1 ≤ ‖1 + (x : ℂ)‖ ^ (k + 1) := by
        apply one_le_pow₀; have := one_le_norm_one_add (x : ℂ); simp at this; linarith
      have h2 : x ^ k ≤ x := by
        calc x ^ k ≤ x ^ 1 := pow_le_pow_of_le_one hx0.le hx.2 hk
          _ = x := pow_one x
      calc x ^ k / ‖1 + (x : ℂ)‖ ^ (k + 1) ≤ x ^ k / 1 :=
            div_le_div_of_nonneg_left (by positivity) one_pos h1
        _ ≤ x := by rw [div_one]; exact h2
    rw [fM, norm_mul, norm_mul, hE, mul_one]
    calc ‖H x‖ * ‖Q k x‖ ≤ 1 / x * x :=
          mul_le_mul hH hQ (norm_nonneg _) (by positivity)
      _ = 1 := by field_simp
  · have hc1 : ContinuousOn (fun x : ℝ => fM k x) (Ioi 1) := hc.mono (Ioi_subset_Ioi zero_le_one)
    have hg : IntegrableOn (fun x : ℝ => √(π / 2) * x ^ (-(3 / 2) : ℝ)) (Ioi 1) :=
      (integrableOn_Ioi_rpow_of_lt (by norm_num) one_pos).const_mul _
    refine Integrable.mono' hg (hc1.aestronglyMeasurable measurableSet_Ioi) ?_
    refine (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall fun x hx => ?_)
    have hx0 : 0 < x := lt_trans one_pos hx
    have hH := H_le_rpow (x : ℂ) (by simpa using hx0)
    simp only [Complex.ofReal_re] at hH
    have hE : ‖E (x : ℂ)‖ = 1 := by rw [E_norm]; simp
    have hQ : ‖Q k (x : ℂ)‖ ≤ 1 := by
      refine (Q_le k x (by simp; linarith)).trans ?_
      have := one_le_norm_one_add (x : ℂ); simp at this
      rw [div_le_one (by linarith)]; linarith
    rw [fM, norm_mul, norm_mul, hE, mul_one]
    calc ‖H x‖ * ‖Q k x‖ ≤ (√(π / 2) * x ^ (-(3 / 2) : ℝ)) * 1 :=
          mul_le_mul hH hQ (norm_nonneg _) (by positivity)
      _ = _ := mul_one _

lemma Lpt_re_pos (r : ℝ) (hr : 0 ≤ r) : 0 < (Lpt r).re := by
  rw [Lpt_re]; positivity

lemma fM_line_int (k : ℕ) : IntegrableOn (fun r : ℝ => fM k (Lpt r)) (Ioi 0) := by
  have hLc : Continuous Lpt := by unfold Lpt; fun_prop
  have hc : ContinuousOn (fun r : ℝ => fM k (Lpt r)) (Ioi 0) :=
    (fM_cont k).comp hLc.continuousOn (fun r hr => Lpt_re_pos r (le_of_lt hr))
  have hg : IntegrableOn (fun r : ℝ => 2 * Real.exp (-a * r)) (Ioi 0) :=
    (exp_neg_integrableOn_Ioi 0 a_pos).const_mul 2
  refine Integrable.mono' hg (hc.aestronglyMeasurable measurableSet_Ioi) ?_
  refine (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall fun r hr => ?_)
  have hr0 : 0 ≤ r := le_of_lt hr
  obtain ⟨hH, hE, -, -⟩ := stmt7_at r hr0
  have hre := Lpt_re_pos r hr0
  have hQ : ‖Q k (Lpt r)‖ ≤ 1 := by
    refine (Q_le k _ hre.le).trans ?_
    have := one_le_norm_one_add (Lpt r)
    rw [div_le_one (by linarith)]; linarith
  rw [fM, norm_mul, norm_mul, hE]
  calc ‖H (Lpt r)‖ * Real.exp (-a * r) * ‖Q k (Lpt r)‖ ≤ 2 * Real.exp (-a * r) * 1 := by
        gcongr
    _ = _ := mul_one _

end Q05

namespace Q05
open Real MeasureTheory Set Filter

lemma mem_ball_aux (u v X : ℝ) (hX : 1 / 4 < X) (h : (u - X) ^ 2 + v ^ 2 < (X - 1 / 4) ^ 2) :
    ((u : ℂ) + (v : ℂ) * Complex.I) ∈ Metric.ball (X : ℂ) (X - 1 / 4) := by
  rw [Metric.mem_ball, dist_eq_norm, Complex.norm_eq_sqrt_sq_add_sq]
  rw [Real.sqrt_lt' (by linarith)]
  simpa using h

lemma ball_sub (X : ℝ) : Metric.ball (X : ℂ) (X - 1 / 4) ⊆ {s : ℂ | 0 < s.re} := by
  intro z hz
  rw [Metric.mem_ball, dist_eq_norm] at hz
  have := Complex.abs_re_le_norm (z - X)
  simp only [Complex.sub_re, Complex.ofReal_re] at this
  show 0 < z.re
  have := neg_abs_le (z.re - X)
  linarith

lemma sq2half (t : ℝ) : t * √2 * (√2 / 2) = t := by
  rw [show t * √2 * (√2 / 2) = t * (√2 * √2) / 2 by ring, sqrt2_sq]; ring

lemma Lpt_T (R : ℝ) : Lpt ((R - 1 / 2) * √2) = (R : ℂ) + ((R - 1 / 2 : ℝ) : ℂ) * Complex.I := by
  rw [Lpt_eq, sq2half]
  apply Complex.ext <;> simp

lemma contour_R (k : ℕ) (R : ℝ) (hR : 1 ≤ R) :
    ∫ r in (0 : ℝ)..((R - 1 / 2) * √2), fM k (Lpt r) * ω =
      (∫ x in (1 / 2 : ℝ)..R, fM k x) +
        ∫ y in (0 : ℝ)..(R - 1 / 2), fM k (R + y * Complex.I) * Complex.I := by
  set X : ℝ := 2 * R ^ 2 with hX
  have hX1 : 2 ≤ X := by nlinarith
  have hdiff : DifferentiableOn ℂ (fM k) (Metric.ball (X : ℂ) (X - 1 / 4)) :=
    (fM_diff k).mono (ball_sub X)
  obtain ⟨F, hF⟩ := hdiff.isExactOn_ball
  have hcont : ContinuousOn (fM k) (Metric.ball (X : ℂ) (X - 1 / 4)) := hdiff.continuousOn
  -- membership
  have hmR : ∀ x ∈ uIcc (1 / 2 : ℝ) R, (x : ℂ) ∈ Metric.ball (X : ℂ) (X - 1 / 4) := by
    intro x hx
    rw [uIcc_of_le (by linarith)] at hx
    have := mem_ball_aux x 0 X (by linarith) (by nlinarith [hx.1, hx.2])
    simpa using this
  have hmV : ∀ y ∈ uIcc (0 : ℝ) (R - 1 / 2),
      ((R : ℂ) + (y : ℂ) * Complex.I) ∈ Metric.ball (X : ℂ) (X - 1 / 4) := by
    intro y hy
    rw [uIcc_of_le (by linarith)] at hy
    exact mem_ball_aux R y X (by linarith) (by nlinarith [hy.1, hy.2])
  have h22 := sqrt2_sq
  have hs2 : (0 : ℝ) < √2 := by positivity
  have hmL : ∀ r ∈ uIcc (0 : ℝ) ((R - 1 / 2) * √2), Lpt r ∈ Metric.ball (X : ℂ) (X - 1 / 4) := by
    intro r hr
    have hRh0 : 0 ≤ R - 1 / 2 := by linarith
    have hRh : 0 ≤ (R - 1 / 2) * √2 := by positivity
    rw [uIcc_of_le hRh] at hr
    rw [Lpt_eq]
    set u := r * (√2 / 2)
    have hu0 : 0 ≤ u := by have := hr.1; positivity
    have hu1 : u ≤ R - 1 / 2 := by
      have := hr.2
      have : r * (√2 / 2) ≤ (R - 1 / 2) * √2 * (√2 / 2) :=
        mul_le_mul_of_nonneg_right this (by positivity)
      rw [sq2half] at this
      exact this
    apply mem_ball_aux _ _ X (by linarith)
    have h1 : u * u ≤ u * R := mul_le_mul_of_nonneg_left (by linarith) hu0
    have h2 : u * R ≤ u * R ^ 2 := mul_le_mul_of_nonneg_left (by nlinarith) hu0
    nlinarith
  -- FTC on the three segments
  have hA : ∫ x in (1 / 2 : ℝ)..R, fM k x = F R - F (1 / 2 : ℝ) := by
    apply intervalIntegral.integral_eq_sub_of_hasDerivAt
    · intro x hx; exact (hF _ (hmR x hx)).comp_ofReal
    · apply ContinuousOn.intervalIntegrable
      exact hcont.comp Complex.continuous_ofReal.continuousOn (fun x hx => hmR x hx)
  have hV : ∫ y in (0 : ℝ)..(R - 1 / 2), fM k (R + y * Complex.I) * Complex.I =
      F (R + ((R - 1 / 2 : ℝ) : ℂ) * Complex.I) - F (R + ((0 : ℝ) : ℂ) * Complex.I) := by
    apply intervalIntegral.integral_eq_sub_of_hasDerivAt (f := fun y : ℝ => F (R + y * Complex.I))
    · intro y hy
      have h1 : HasDerivAt (fun z : ℂ => (R : ℂ) + z * Complex.I) Complex.I (y : ℂ) := by
        simpa using ((hasDerivAt_id (y : ℂ)).mul_const Complex.I).const_add (R : ℂ)
      exact ((hF _ (hmV y hy)).comp (y : ℂ) h1).comp_ofReal
    · apply ContinuousOn.intervalIntegrable
      refine ContinuousOn.mul ?_ continuousOn_const
      refine hcont.comp (by fun_prop) (fun y hy => hmV y hy)
  have hB : ∫ r in (0 : ℝ)..((R - 1 / 2) * √2), fM k (Lpt r) * ω =
      F (Lpt ((R - 1 / 2) * √2)) - F (Lpt 0) := by
    apply intervalIntegral.integral_eq_sub_of_hasDerivAt (f := fun r : ℝ => F (Lpt r))
    · intro r hr
      have h1 : HasDerivAt (fun z : ℂ => (1 / 2 : ℂ) + z * ω) ω (r : ℂ) := by
        simpa using ((hasDerivAt_id (r : ℂ)).mul_const ω).const_add (1 / 2 : ℂ)
      have hFr : HasDerivAt F (fM k (Lpt r)) ((1 / 2 : ℂ) + (r : ℂ) * ω) := hF (Lpt r) (hmL r hr)
      exact (hFr.comp (r : ℂ) h1).comp_ofReal
    · apply ContinuousOn.intervalIntegrable
      refine ContinuousOn.mul ?_ continuousOn_const
      have hLc : Continuous Lpt := by unfold Lpt; fun_prop
      exact hcont.comp hLc.continuousOn (fun r hr => hmL r hr)
  rw [hA, hV, hB, Lpt_T]
  have h0 : Lpt 0 = ((1 / 2 : ℝ) : ℂ) := by rw [Lpt]; push_cast; ring
  rw [h0]
  push_cast
  rw [zero_mul, add_zero]
  ring

lemma vert_bound (k : ℕ) (R : ℝ) (hR : 1 ≤ R) :
    ‖∫ y in (0 : ℝ)..(R - 1 / 2), fM k (R + y * Complex.I) * Complex.I‖ ≤ R⁻¹ := by
  have hb : ∀ y ∈ Set.uIoc (0 : ℝ) (R - 1 / 2), ‖fM k (R + y * Complex.I) * Complex.I‖ ≤ 1 / R ^ 2 := by
    intro y hy
    rw [uIoc_of_le (by linarith)] at hy
    have hy0 : 0 ≤ y := le_of_lt hy.1
    set s : ℂ := (R : ℂ) + (y : ℂ) * Complex.I
    have hre : s.re = R := by simp [s]
    have him : s.im = y := by simp [s]
    have hH := H_le_one_div s (by rw [hre]; linarith)
    rw [hre] at hH
    have hE : ‖E s‖ ≤ 1 := by
      rw [E_norm, him]; apply Real.exp_le_one_iff.2
      have := Real.pi_pos; nlinarith
    have hQ : ‖Q k s‖ ≤ 1 / R := by
      refine (Q_le k s (by rw [hre]; linarith)).trans ?_
      have := one_le_norm_one_add s
      rw [hre] at this
      exact div_le_div_of_nonneg_left zero_le_one (by linarith) (by linarith)
    rw [norm_mul, Complex.norm_I, mul_one, fM, norm_mul, norm_mul]
    calc ‖H s‖ * ‖E s‖ * ‖Q k s‖ ≤ (1 / R) * 1 * (1 / R) := by gcongr
      _ = 1 / R ^ 2 := by ring
  have := intervalIntegral.norm_integral_le_of_norm_le_const hb
  refine this.trans ?_
  rw [sub_zero, abs_of_nonneg (by linarith)]
  have hR0 : 0 < R := by linarith
  rw [div_mul_eq_mul_div, one_mul, div_le_iff₀ (by positivity)]
  have : R⁻¹ * R ^ 2 = R := by field_simp
  rw [this]; linarith

lemma shift (k : ℕ) (hk : 1 ≤ k) :
    ∫ x in Ioi (1 / 2 : ℝ), fM k x = lineInt (fM k) := by
  have hI1 : IntegrableOn (fun x : ℝ => fM k x) (Ioi (1 / 2)) :=
    (fM_real_int k hk).mono_set (Ioi_subset_Ioi (by norm_num))
  have hI2 : IntegrableOn (fun r : ℝ => fM k (Lpt r) * ω) (Ioi 0) := by
    have := fM_line_int k
    exact this.mul_const ω
  have tA := intervalIntegral_tendsto_integral_Ioi (1 / 2 : ℝ) hI1 tendsto_id
  have hT : Tendsto (fun R : ℝ => (R - 1 / 2) * √2) atTop atTop :=
    (tendsto_atTop_add_const_right _ _ tendsto_id).atTop_mul_const (by positivity)
  have tB := intervalIntegral_tendsto_integral_Ioi (0 : ℝ) hI2 hT
  have tV : Tendsto (fun R : ℝ => ∫ y in (0 : ℝ)..(R - 1 / 2),
      fM k (R + y * Complex.I) * Complex.I) atTop (nhds 0) := by
    refine squeeze_zero_norm' ?_ tendsto_inv_atTop_zero
    filter_upwards [eventually_ge_atTop (1 : ℝ)] with R hR
    exact vert_bound k R hR
  have tAV := tA.add tV
  rw [add_zero] at tAV
  have heq : (fun R : ℝ => ∫ r in (0 : ℝ)..((R - 1 / 2) * √2), fM k (Lpt r) * ω) =ᶠ[atTop]
      (fun R : ℝ => (∫ x in (1 / 2 : ℝ)..R, fM k x) +
        ∫ y in (0 : ℝ)..(R - 1 / 2), fM k (R + y * Complex.I) * Complex.I) := by
    filter_upwards [eventually_ge_atTop (1 : ℝ)] with R hR
    exact contour_R k R hR
  exact tendsto_nhds_unique tAV (tB.congr' heq)

end Q05

theorem Q05.proof_11 : Q05.stmt_11 := by
  intro k hk
  have hI1 : IntegrableOn (fun s : ℝ => Q05.fM k s) (Ioi 0) := Q05.fM_real_int k hk
  refine ⟨hI1, Q05.fM_line_int k, fun s hs => ?_, ?_⟩
  · have hs0 : (s : ℂ) ≠ 0 := by exact_mod_cast hs.ne'
    have hs1 : (1 + (s : ℂ)) ≠ 0 := Q05.one_add_ne s (by simp; linarith)
    rw [Q05.Q, div_pow, pow_succ (s : ℂ)]
    field_simp
  · have hsplit : ∫ s in Ioi (0 : ℝ), Q05.fM k s =
        (∫ s in Ioc (0 : ℝ) (1 / 2), Q05.fM k s) + ∫ s in Ioi (1 / 2 : ℝ), Q05.fM k s := by
      rw [← Ioc_union_Ioi_eq_Ioi (by norm_num : (0 : ℝ) ≤ 1 / 2)]
      exact setIntegral_union (Ioc_disjoint_Ioi le_rfl) measurableSet_Ioi
        (hI1.mono_set (Ioc_subset_Ioi_self)) (hI1.mono_set (Ioi_subset_Ioi (by norm_num)))
    have := Q05.shift k hk
    show ∫ s in Ioi (0 : ℝ), Q05.fM k s =
      (∫ s in (0 : ℝ)..(1 / 2), Q05.fM k s) + Q05.lineInt (Q05.fM k)
    rw [hsplit, intervalIntegral.integral_of_le (by norm_num), this]
