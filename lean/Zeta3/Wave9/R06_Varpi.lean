import Mathlib
import Zeta3.Wave8.Q07_DecayRest
/-!
# R06 -- the first correction `ϖ₁` (Theorem `thm:varpi`)

Target id: R06.
Paper labels: `thm:varpi` (paper/ieee/main.tex, lines 1301-1308; `ϖ₁` defined at line 1275),
in the supplement `g5:thm:varpi` (paper/ieee/proofs/sec_g5_varpi1.tex, statement lines 6-19,
plan line 21, proof of Theorem 13 line 257-259), and `g5:prop:rest` (b)
(paper/ieee/proofs/sec_g5_decay.tex, lines 193-203, gap note line 205, proof of (b) lines
229-235).

Quoted verbatim (main.tex, thm:varpi):

  For every $k$ with $N=k+1\ge100$,
  \Bigl|M_k+\frac{\sqrt\pi}{N}\,e^{-\theta}\operatorname{Re}\Bigl[e^{i\theta}\Bigl(1+\frac{\varpi_1}{\sqrt N}\Bigr)\Bigr]\Bigr|\le C\,N^{-2}e^{-\theta},
  with $C=1{,}000{,}668$.

Quoted verbatim (main.tex line 1275, and sec_g5_varpi1.tex lines 8-10):

  \varpi_1=\bigl[(24+12\pi+\pi^2)+i\,(24-12\pi+\pi^2)\bigr]/(48\sqrt\pi)=0.84121508777-0.04501183769\,i

Quoted verbatim (sec_g5_varpi1.tex, plan, line 21):

  Put $X=\frac{\sqrt\pi}N\eu^{-\theta}\operatorname{Re}[\eu^{\ii\theta}(1+\varpi_1/\sqrt N)]$. Then
  $|M_k+X|\le|M_k-M_k^{\mathrm{main}}|+|M_k^{\mathrm{main}}+X|$. The first term is bounded by
  Proposition~\ref{g5:prop:rest}(b): at most $c_{\mathrm{rest}}N^{-2}\eu^{-\theta}$ with
  $c_{\mathrm{rest}}<120{,}000$. The second term is bounded in Proposition~\ref{g5:prop:main}
  below by $C_{\mathrm{main}}N^{-2}\eu^{-\theta}$ with $C_{\mathrm{main}}=71{,}000$. So the theorem
  holds with $191{,}000$, and a fortiori with $1{,}000{,}668$.

Quoted verbatim (sec_g5_decay.tex, g5:prop:rest (b)):

  (b) $N^2\eu^{\sqrt{\pi N}}\,|M_k-M_k^{\mathrm{main}}|\le c_{\mathrm{rest}}$ with $c_{\mathrm{rest}}<120{,}000$.

Quoted verbatim (sec_g5_decay.tex, proof of (b), the definition of `c_rest`):

  Let $p_0,\dots,p_5$ be the six terms of the displayed inequality multiplied by
  $N^2\eu^{\sqrt{\pi N}}$ [...] So these pieces decrease for $N\ge100$. [...] Therefore
  $N^2\eu^{\sqrt{\pi N}}|M_k-\Mmain|\le\sum_ip_i(100)=c_{\mathrm{rest}}$, which the script
  evaluates and bounds by $120{,}000$.

Quoted verbatim (sec_g5_varpi1.tex, proof of Theorem 13):

  Therefore $|M_k+X|\le191{,}000\,N^{-2}\eu^{-\theta}\le1{,}000{,}668\,N^{-2}\eu^{-\theta}$.

## Conventions

* `M_k` is `Q05.M k := (W02.Mk k).re`, the real part of the library's Taylor coefficient
  `W02.Mk k` of `M = -F/2` (Wave 3, W02). Q07.stmt_2 (g5:lem:Mk, proved) shows that `W02.Mk k`
  is itself a real number, so `Q05.M k` is exactly the paper's `M_k`, and nothing is lost by
  taking the real part. `M_k^main` is `Q05.Mmain k`, the integral of main.tex line 1313.
* `ϖ₁` is `Q05.varpi1`, defined (inside Q05's statement block) as the complex number with real
  part `(24+12π+π²)/(48√π)` and imaginary part `(24-12π+π²)/(48√π)`; this is the paper's
  formula, and stmt_3 records it together with the printed decimals.
* `X` is `Q05.X N = (√π/N) e^{-√(πN)} Re[e^{i√(πN)}(1+ϖ₁/√N)]` (Q05's block), the `X` of the
  supplement's plan, with `N = k+1` cast to `ℝ`. The final statement stmt_5 does not use `X`:
  it writes the expression out.
* `θ = √(πN)`; `N^{-2} e^{-θ}` is written `C / N^2 * exp(-√(πN))` (as in Q07.stmt_6).
* `restBound k` is the right-hand side of the displayed inequality of g5:prop:rest, written
  exactly as in Q07.stmt_4 (proved there: `|M_k - M_k^main| ≤ restBound k` for `N ≥ 100`).
  For `k ≥ 99` the quantities `k-1`, `k-2` are real subtractions of positive numbers.
* `crest` is the supplement's `c_rest = Σ_i p_i(100)`, i.e. `100² e^{√(100π)} restBound 99`
  (the six pieces at `N = 100`, `k = 99`). Numerically (mpmath, 40 digits) `crest =
  114,536.103195...`, the six pieces being `31918.94, 71215.51, 11125.48, 0.0017, 276.17,
  3.2·10^{-10}` (in the order of the display: the supplement lists the `B_2` piece before the
  `e^{eξ}` piece), matching the gap note's `114,536.10...`.

## Statements

* stmt_1: the elementary estimates of the proof of g5:prop:rest (b): for `N ≥ 100`,
  `N² e^{√(πN)} restBound k ≤ crest` (each piece decreases in `N`), and `crest < 120,000`.
* stmt_2: g5:prop:rest (b): for `N ≥ 100`, `N² e^{√(πN)} |M_k - M_k^main| ≤ crest`, equivalently
  `|M_k - M_k^main| ≤ crest N^{-2} e^{-θ}`, and `crest < 120,000`.
* stmt_3: the definition of `ϖ₁` and of `X` unfolded: `Re ϖ₁`, `Im ϖ₁` as printed, the printed
  decimals correct to `5·10^{-12}` (rounded to 11 places), and
  `X(N) = (√π/N)e^{-θ}[cos θ (1 + Re ϖ₁/√N) - sin θ Im ϖ₁/√N]` for `N > 0`.
* stmt_4: the supplement's intermediate bound `|M_k + X| ≤ 191,000 N^{-2} e^{-θ}`, `N ≥ 100`.
* stmt_5: thm:varpi with `C = 1,000,668`, the expression written out.

## NOT COVERED

* The proof of g5:prop:main (Regions I and II of sec_g5_varpi1.tex) is not restated: the bound
  `|M_k^main + X| ≤ 71,000 N^{-2} e^{-θ}` is Q07.stmt_6 (proved) and is imported, not restated.
* The paper's own route to `C = 1,000,668` (main.tex line 1310: interval arithmetic giving
  `10^6` for the main term, the rest bound `1030 N^{-1/2} e^{-√(2πN)}` and `667.357`) is not
  stated; the supplement replaces it by `71,000 + 120,000`, and the supplement says the rest bound
  `1030 N^{-1/2} e^{-√(2πN)}` is not proved there. The constant `1.76` for `100 ≤ N ≤ 400` is not
  stated either.
* No high-precision values of `M_k` were computed for this file; only the arithmetic of `c_rest`
  and the monotonicity of the six pieces (`k = 99..20000`) and the decimals of `ϖ₁` were checked
  numerically.
-/
-- BEGIN STATEMENTS
namespace R06

open Real Q05

/-- The right-hand side of the displayed inequality of g5:prop:rest (as in Q07.stmt_4):
`B_1[e^{π/4}e^{-√(3πk/2)} + ξ e^{-√(3π(k-1)/2)} + (e^ξ-1-ξ) e^{-√(9π(k-2)/4)} + e^{eξ}e^{-k/2}]
 + B_2 e^{-√(9πk/4)} + 2e^{-(N-1)/2}`, `N = k+1`. -/
noncomputable def restBound (k : ℕ) : ℝ :=
  B1 * (Real.exp (π / 4) * Real.exp (-√(3 * π * k / 2)) +
      ξ * Real.exp (-√(3 * π * ((k : ℝ) - 1) / 2)) +
      (Real.exp ξ - 1 - ξ) * Real.exp (-√(9 * π * ((k : ℝ) - 2) / 4)) +
      Real.exp (Real.exp 1 * ξ) * Real.exp (-(k : ℝ) / 2)) +
    B2 * Real.exp (-√(9 * π * k / 4)) + 2 * Real.exp (-(((k + 1 : ℕ) : ℝ) - 1) / 2)

/-- `c_rest = Σ_i p_i(100)`: the six pieces of the display times `N² e^{√(πN)}`, at `N = 100`
(`k = 99`). -/
noncomputable def crest : ℝ :=
  (100 : ℝ) ^ 2 * Real.exp (√(π * 100)) * restBound 99

/-- Proof of g5:prop:rest (b): "So these pieces decrease for $N\ge100$ [...]
$\sum_ip_i(100)=c_{\mathrm{rest}}$, which the script evaluates and bounds by $120{,}000$":
`N² e^{√(πN)} restBound k ≤ c_rest` for `N = k+1 ≥ 100`, and `c_rest < 120,000`. -/
def stmt_1 : Prop :=
  (∀ k : ℕ, 100 ≤ k + 1 →
      (((k + 1 : ℕ) : ℝ)) ^ 2 * Real.exp (√(π * ((k + 1 : ℕ) : ℝ))) * restBound k ≤ crest) ∧
    crest < 120000

/-- g5:prop:rest (b): "$N^2\eu^{\sqrt{\pi N}}\,|M_k-M_k^{\mathrm{main}}|\le c_{\mathrm{rest}}$
with $c_{\mathrm{rest}}<120{,}000$", for `N = k+1 ≥ 100`, also in the form used in the plan of
Theorem 13: `|M_k - M_k^main| ≤ c_rest N^{-2} e^{-θ}`. -/
def stmt_2 : Prop :=
  (∀ k : ℕ, 100 ≤ k + 1 →
      (((k + 1 : ℕ) : ℝ)) ^ 2 * Real.exp (√(π * ((k + 1 : ℕ) : ℝ))) * |M k - Mmain k| ≤ crest ∧
        |M k - Mmain k| ≤
          crest / (((k + 1 : ℕ) : ℝ)) ^ 2 * Real.exp (-√(π * ((k + 1 : ℕ) : ℝ)))) ∧
    crest < 120000

/-- The definitions of `ϖ₁` and `X` (main.tex line 1275; sec_g5_varpi1.tex lines 8-10, 21):
`ϖ₁ = [(24+12π+π²) + i(24-12π+π²)]/(48√π) = 0.84121508777 - 0.04501183769 i` (decimals rounded
to 11 places), and `X = (√π/N) e^{-θ} Re[e^{iθ}(1+ϖ₁/√N)]` in real form. -/
def stmt_3 : Prop :=
  varpi1.re = (24 + 12 * π + π ^ 2) / (48 * √π) ∧
    varpi1.im = (24 - 12 * π + π ^ 2) / (48 * √π) ∧
    |varpi1.re - 0.84121508777| < 5e-12 ∧
    |varpi1.im - (-0.04501183769)| < 5e-12 ∧
    ∀ N : ℝ, 0 < N →
      X N = √π / N * Real.exp (-√(π * N)) *
        (Real.cos (√(π * N)) * (1 + (24 + 12 * π + π ^ 2) / (48 * √π) / √N) -
          Real.sin (√(π * N)) * ((24 - 12 * π + π ^ 2) / (48 * √π) / √N))

/-- Proof of Theorem 13 (sec_g5_varpi1.tex): "Therefore
$|M_k+X|\le191{,}000\,N^{-2}\eu^{-\theta}$", for `N = k+1 ≥ 100`. -/
def stmt_4 : Prop :=
  ∀ k : ℕ, 100 ≤ k + 1 →
    |M k + X ((k + 1 : ℕ) : ℝ)| ≤
      191000 / (((k + 1 : ℕ) : ℝ)) ^ 2 * Real.exp (-√(π * ((k + 1 : ℕ) : ℝ)))

/-- thm:varpi: "For every $k$ with $N=k+1\ge100$,
$|M_k+\frac{\sqrt\pi}{N}e^{-\theta}\operatorname{Re}[e^{i\theta}(1+\frac{\varpi_1}{\sqrt N})]|\le
C N^{-2}e^{-\theta}$, with $C=1{,}000{,}668$", `θ = √(πN)`. -/
def stmt_5 : Prop :=
  ∀ k : ℕ, 100 ≤ k + 1 →
    |M k + √π / ((k + 1 : ℕ) : ℝ) * Real.exp (-√(π * ((k + 1 : ℕ) : ℝ))) *
        (Complex.exp (Complex.I * ((√(π * ((k + 1 : ℕ) : ℝ)) : ℝ) : ℂ)) *
          (1 + varpi1 / ((√((k + 1 : ℕ) : ℝ) : ℝ) : ℂ))).re| ≤
      1000668 / (((k + 1 : ℕ) : ℝ)) ^ 2 * Real.exp (-√(π * ((k + 1 : ℕ) : ℝ)))

end R06
-- END STATEMENTS

namespace R06.P1

open Real Q05

noncomputable def g (c d N : ℝ) : ℝ := N ^ 2 * Real.exp (√(π * N)) * Real.exp (-√(c * (N - d)))

noncomputable def h (N : ℝ) : ℝ := N ^ 2 * Real.exp (√(π * N)) * Real.exp (-(N - 1) / 2)

noncomputable def R (N : ℝ) : ℝ :=
  B1 * (Real.exp (π / 4) * g (3 * π / 2) 1 N + ξ * g (3 * π / 2) 2 N +
      (Real.exp ξ - 1 - ξ) * g (9 * π / 4) 3 N + Real.exp (Real.exp 1 * ξ) * h N) +
    B2 * g (9 * π / 4) 1 N + 2 * h N

lemma eqR (k : ℕ) :
    (((k + 1 : ℕ) : ℝ)) ^ 2 * Real.exp (√(π * ((k + 1 : ℕ) : ℝ))) * R06.restBound k =
      R ((k + 1 : ℕ) : ℝ) := by
  unfold R g h R06.restBound
  set N : ℝ := ((k + 1 : ℕ) : ℝ) with hN
  have hx : (k : ℝ) = N - 1 := by rw [hN]; push_cast; ring
  rw [hx, show 3 * π * (N - 1) / 2 = 3 * π / 2 * (N - 1) by ring,
    show 3 * π * (N - 1 - 1) / 2 = 3 * π / 2 * (N - 2) by ring,
    show 9 * π * (N - 1 - 2) / 4 = 9 * π / 4 * (N - 3) by ring,
    show 9 * π * (N - 1) / 4 = 9 * π / 4 * (N - 1) by ring]
  ring

lemma deriv_f (c d x : ℝ) (hx : 0 < x) (hxd : 0 < x - d) :
    HasDerivAt (fun N => 2 * Real.log N + √π * √N - √c * √(N - d))
      (2 * x⁻¹ + √π * (1 / (2 * √x)) - √c * (1 / (2 * √(x - d)))) x := by
  have h1 := (Real.hasDerivAt_log hx.ne').const_mul 2
  have h2 := (Real.hasDerivAt_sqrt hx.ne').const_mul √π
  have h3 := (((hasDerivAt_id x).sub_const d).sqrt hxd.ne').const_mul √c
  exact (h1.add h2).sub h3

lemma sqrt_ge_ten (N : ℝ) (hN : 100 ≤ N) : 10 ≤ √N := by
  rw [show (10 : ℝ) = √100 by rw [show (100 : ℝ) = 10 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]]
  exact Real.sqrt_le_sqrt hN

lemma dneg (c d x : ℝ) (hd1 : 1 ≤ d) (hd3 : d ≤ 3) (hc : √π + 0.398 ≤ √c) (hx : 100 ≤ x) :
    2 * x⁻¹ + √π * (1 / (2 * √x)) - √c * (1 / (2 * √(x - d))) ≤ 0 := by
  have hsp := Q05.sqrt_pi_gt
  set t := √x with ht
  set u := √(x - d) with hu
  have ht10 : 10 ≤ t := sqrt_ge_ten x hx
  have htt : t * t = x := Real.mul_self_sqrt (by linarith)
  have huu : u * u = x - d := Real.mul_self_sqrt (by linarith)
  have hu0 : 0 < u := Real.sqrt_pos.2 (by linarith)
  have htu : 2 * t * u ≤ 2 * t ^ 2 - d := by
    nlinarith [sq_nonneg d, mul_pos hu0 (by linarith : (0 : ℝ) < t)]
  have key : (2 * (2 * t) + √π * t ^ 2) * (2 * u) ≤ √c * (t ^ 2 * (2 * t)) := by
    have e1 : (2 * (2 * t) + √π * t ^ 2) * (2 * u) = (4 + √π * t) * (2 * t * u) := by ring
    have h4 : 0 ≤ 4 + √π * t := by nlinarith
    have s1 := mul_le_mul_of_nonneg_left htu h4
    have s2 : 0 ≤ t ^ 3 := by positivity
    have s3 : (4 + √π * t) * (2 * t ^ 2 - d) ≤ (√π + 0.398) * (t ^ 2 * (2 * t)) := by
      nlinarith [mul_le_mul_of_nonneg_left hd1 (by nlinarith : (0 : ℝ) ≤ √π * t)]
    have s4 : (√π + 0.398) * (t ^ 2 * (2 * t)) ≤ √c * (t ^ 2 * (2 * t)) :=
      mul_le_mul_of_nonneg_right hc (by positivity)
    linarith
  have hx' : x = t ^ 2 := by rw [← htt]; ring
  rw [hx']
  have : 2 * (t ^ 2)⁻¹ + √π * (1 / (2 * t)) - √c * (1 / (2 * u)) =
      ((2 * (2 * t) + √π * t ^ 2) * (2 * u) - √c * (t ^ 2 * (2 * t))) /
        (t ^ 2 * (2 * t) * (2 * u)) := by
    field_simp
  rw [this]
  apply div_nonpos_of_nonpos_of_nonneg (by linarith) (by positivity)


lemma anti (c d : ℝ) (hd1 : 1 ≤ d) (hd3 : d ≤ 3) (hc : √π + 0.398 ≤ √c) :
    AntitoneOn (fun N => 2 * Real.log N + √π * √N - √c * √(N - d)) (Set.Ici 100) := by
  have hD : ∀ x ∈ Set.Ici (100 : ℝ), HasDerivAt (fun N => 2 * Real.log N + √π * √N - √c * √(N - d))
      (2 * x⁻¹ + √π * (1 / (2 * √x)) - √c * (1 / (2 * √(x - d)))) x := fun x hx =>
    deriv_f c d x (by linarith [Set.mem_Ici.1 hx]) (by linarith [Set.mem_Ici.1 hx])
  apply antitoneOn_of_deriv_nonpos (convex_Ici 100)
  · exact fun x hx => (hD x hx).continuousAt.continuousWithinAt
  · rw [interior_Ici]
    exact fun x hx => (hD x (show (100 : ℝ) ≤ x from le_of_lt hx)).differentiableAt.differentiableWithinAt
  · rw [interior_Ici]
    intro x hx
    rw [(hD x (show (100 : ℝ) ≤ x from le_of_lt hx)).deriv]
    exact dneg c d x hd1 hd3 hc (le_of_lt hx)

lemma g_eq (c d N : ℝ) (hc : 0 ≤ c) (hN : 0 < N) :
    g c d N = Real.exp (2 * Real.log N + √π * √N - √c * √(N - d)) := by
  unfold g
  rw [Real.sqrt_mul Real.pi_pos.le, Real.sqrt_mul hc,
    show 2 * Real.log N + √π * √N - √c * √(N - d) =
      Real.log N + Real.log N + √π * √N + -(√c * √(N - d)) by ring,
    Real.exp_add, Real.exp_add, Real.exp_add, Real.exp_log hN]
  ring

lemma g_le (c d N : ℝ) (hd1 : 1 ≤ d) (hd3 : d ≤ 3) (hc : √π + 0.398 ≤ √c) (hN : 100 ≤ N) :
    g c d N ≤ g c d 100 := by
  have hc0 : 0 ≤ c := by
    have : 0 < √c := by linarith [Real.sqrt_nonneg π]
    exact (Real.sqrt_pos.1 this).le
  rw [g_eq c d N hc0 (by linarith), g_eq c d 100 hc0 (by norm_num)]
  exact Real.exp_le_exp.2 (anti c d hd1 hd3 hc (Set.mem_Ici.2 le_rfl) hN hN)

lemma deriv_l (x : ℝ) (hx : 0 < x) :
    HasDerivAt (fun N => 2 * Real.log N + √π * √N - (N - 1) / 2)
      (2 * x⁻¹ + √π * (1 / (2 * √x)) - 1 / 2) x := by
  have h1 := (Real.hasDerivAt_log hx.ne').const_mul 2
  have h2 := (Real.hasDerivAt_sqrt hx.ne').const_mul √π
  have h3 := ((hasDerivAt_id x).sub_const 1).div_const 2
  exact (h1.add h2).sub h3

lemma anti_l : AntitoneOn (fun N => 2 * Real.log N + √π * √N - (N - 1) / 2) (Set.Ici 100) := by
  have hD : ∀ x ∈ Set.Ici (100 : ℝ), HasDerivAt (fun N => 2 * Real.log N + √π * √N - (N - 1) / 2)
      (2 * x⁻¹ + √π * (1 / (2 * √x)) - 1 / 2) x := fun x hx =>
    deriv_l x (by linarith [Set.mem_Ici.1 hx])
  apply antitoneOn_of_deriv_nonpos (convex_Ici 100)
  · exact fun x hx => (hD x hx).continuousAt.continuousWithinAt
  · rw [interior_Ici]
    exact fun x hx => (hD x (show (100 : ℝ) ≤ x from le_of_lt hx)).differentiableAt.differentiableWithinAt
  · rw [interior_Ici]
    intro x hx
    rw [(hD x (show (100 : ℝ) ≤ x from le_of_lt hx)).deriv]
    have hx' : (100 : ℝ) ≤ x := le_of_lt hx
    have ht := sqrt_ge_ten x hx'
    have hsp := Q05.sqrt_pi_lt
    have h1 : 2 * x⁻¹ ≤ 0.02 := by
      rw [← div_eq_mul_inv, div_le_iff₀ (by linarith)]; linarith
    have h2 : √π * (1 / (2 * √x)) ≤ 0.09 := by
      rw [mul_one_div, div_le_iff₀ (by linarith)]; nlinarith
    linarith

lemma h_eq (N : ℝ) (hN : 0 < N) :
    h N = Real.exp (2 * Real.log N + √π * √N - (N - 1) / 2) := by
  unfold h
  rw [Real.sqrt_mul Real.pi_pos.le,
    show 2 * Real.log N + √π * √N - (N - 1) / 2 =
      Real.log N + Real.log N + √π * √N + (-(N - 1) / 2) by ring,
    Real.exp_add, Real.exp_add, Real.exp_add, Real.exp_log hN]
  ring

lemma h_le (N : ℝ) (hN : 100 ≤ N) : h N ≤ h 100 := by
  rw [h_eq N (by linarith), h_eq 100 (by norm_num)]
  exact Real.exp_le_exp.2 (anti_l (Set.mem_Ici.2 le_rfl) hN hN)

lemma c1 : √π + 0.398 ≤ √(3 * π / 2) := by
  have hsp := Q05.sqrt_pi_lt
  have hsg := Q05.sqrt_pi_gt
  have hsq : √π * √π = π := Real.mul_self_sqrt Real.pi_pos.le
  apply Real.le_sqrt_of_sq_le
  nlinarith

lemma c2 : √π + 0.398 ≤ √(9 * π / 4) := by
  have hsp := Q05.sqrt_pi_lt
  have hsg := Q05.sqrt_pi_gt
  have hsq : √π * √π = π := Real.mul_self_sqrt Real.pi_pos.le
  apply Real.le_sqrt_of_sq_le
  nlinarith

lemma R_le (N : ℝ) (hN : 100 ≤ N) : R N ≤ R 100 := by
  have g1 := g_le _ 1 N le_rfl (by norm_num) c1 hN
  have g2 := g_le _ 2 N (by norm_num) (by norm_num) c1 hN
  have g3 := g_le _ 3 N (by norm_num) le_rfl c2 hN
  have g4 := g_le _ 1 N le_rfl (by norm_num) c2 hN
  have h1 := h_le N hN
  have hB1 := Q07.PE.B1_pos
  have hB2 := Q07.PE.B2_pos
  have hxi := Q07.PE.xi_nonneg
  have hex := Q07.PE.expxi_nonneg
  have e1 : 0 ≤ Real.exp (π / 4) := (Real.exp_pos _).le
  have e2 : 0 ≤ Real.exp (Real.exp 1 * ξ) := (Real.exp_pos _).le
  unfold R
  have inner := add_le_add (add_le_add (add_le_add (mul_le_mul_of_nonneg_left g1 e1)
    (mul_le_mul_of_nonneg_left g2 hxi)) (mul_le_mul_of_nonneg_left g3 hex))
    (mul_le_mul_of_nonneg_left h1 e2)
  have := add_le_add (add_le_add (mul_le_mul_of_nonneg_left inner hB1.le)
    (mul_le_mul_of_nonneg_left g4 hB2.le)) (mul_le_mul_of_nonneg_left h1 (by norm_num : (0:ℝ) ≤ 2))
  exact this

end R06.P1

namespace R06.P1

open Real Q05

lemma exp_neg_le' (D E : ℝ) (hD : 0 ≤ D) (hE : 0 < E)
    (h : 1 / E ≤ ∑ i ∈ Finset.range 16, D ^ i / (i.factorial : ℝ)) : Real.exp (-D) ≤ E := by
  have h1 := Real.sum_le_exp_of_nonneg hD 16
  rw [Real.exp_neg, inv_le_comm₀ (Real.exp_pos _) hE, ← one_div]
  linarith

lemma sqrt_pi100 : √(π * 100) ≤ 17.72454 := by
  rw [Real.sqrt_le_left (by norm_num)]; nlinarith [Real.pi_lt_d20]

lemma g100 (c d L D E : ℝ) (hLc : L ^ 2 ≤ c * (100 - d)) (hD : 17.72454 - L ≤ -D)
    (hD0 : 0 ≤ D) (hE : 0 < E)
    (h : 1 / E ≤ ∑ i ∈ Finset.range 16, D ^ i / (i.factorial : ℝ)) :
    g c d 100 ≤ 10000 * E := by
  unfold g
  have hL := Real.le_sqrt_of_sq_le hLc
  have h1 := sqrt_pi100
  rw [mul_assoc, ← Real.exp_add]
  have h2 : Real.exp (√(π * 100) + -√(c * (100 - d))) ≤ E :=
    (Real.exp_le_exp.2 (by linarith)).trans (exp_neg_le' D E hD0 hE h)
  rw [show (100 : ℝ) ^ 2 = 10000 by norm_num]
  linarith

lemma h100 : h 100 ≤ 10000 * 2.2e-11 := by
  unfold h
  have h1 := sqrt_pi100
  rw [mul_assoc, ← Real.exp_add]
  have h2 : Real.exp (√(π * 100) + -(100 - 1) / 2) ≤ 2.2e-11 :=
    (Real.exp_le_exp.2 (by linarith)).trans (exp_neg_le' 31.77546 2.2e-11 (by norm_num)
      (by norm_num) (by norm_num [Finset.sum_range_succ, Nat.factorial]))
  rw [show (100 : ℝ) ^ 2 = 10000 by norm_num]
  linarith

end R06.P1

namespace R06.P1

open Real Q05

lemma consts : B1 < 70.1 ∧ B2 < 171 ∧ Real.exp (π / 4) < 2.1933 ∧ ξ < 4.3866 ∧
    Real.exp ξ < 80.37 ∧ Real.exp (Real.exp 1 * ξ) < 150860 := by
  unfold B1 B2 ξ
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
  exact ⟨hB1, hB2, he3, hxi, he4, he5⟩

lemma g_nonneg (c d N : ℝ) : 0 ≤ g c d N := by unfold g; positivity

lemma h_nonneg (N : ℝ) : 0 ≤ h N := by unfold h; positivity

lemma R100 : R 100 < 120000 := by
  have g1 := g100 (3 * π / 2) 1 21.59922 3.87468 0.020762 (by nlinarith [Real.pi_gt_d20])
    (by norm_num) (by norm_num) (by norm_num) (by norm_num [Finset.sum_range_succ, Nat.factorial])
  have g2 := g100 (3 * π / 2) 2 21.48986 3.76532 0.023161 (by nlinarith [Real.pi_gt_d20])
    (by norm_num) (by norm_num) (by norm_num) (by norm_num [Finset.sum_range_succ, Nat.factorial])
  have g3 := g100 (9 * π / 4) 3 26.18496 8.46042 0.0002146 (by nlinarith [Real.pi_gt_d20])
    (by norm_num) (by norm_num) (by norm_num) (by norm_num [Finset.sum_range_succ, Nat.factorial])
  have g4 := g100 (9 * π / 4) 1 26.45353 8.72899 0.0001647 (by nlinarith [Real.pi_gt_d20])
    (by norm_num) (by norm_num) (by norm_num) (by norm_num [Finset.sum_range_succ, Nat.factorial])
  have h1 := h100
  obtain ⟨hB1, hB2, he, hxi, hexi, heexi⟩ := consts
  have hxi0 := Q07.PE.xi_nonneg
  have hex0 := Q07.PE.expxi_nonneg
  have hxil : 1 + π / 4 ≤ Real.exp (π / 4) := by
    have := Real.add_one_le_exp (π / 4); linarith
  have hC : Real.exp ξ - 1 - ξ ≤ 76 := by
    unfold ξ at hexi ⊢; nlinarith [Real.pi_gt_three]
  have q1 := mul_le_mul he.le g1 (g_nonneg _ _ _) (by norm_num)
  have q2 := mul_le_mul hxi.le g2 (g_nonneg _ _ _) (by norm_num)
  have q3 := mul_le_mul hC g3 (g_nonneg _ _ _) (by norm_num)
  have q4 := mul_le_mul heexi.le h1 (h_nonneg _) (by norm_num)
  have q5 := mul_le_mul hB2.le g4 (g_nonneg _ _ _) (by norm_num)
  have i0 : 0 ≤ Real.exp (π / 4) * g (3 * π / 2) 1 100 + ξ * g (3 * π / 2) 2 100 +
      (Real.exp ξ - 1 - ξ) * g (9 * π / 4) 3 100 + Real.exp (Real.exp 1 * ξ) * h 100 := by
    have := g_nonneg (3 * π / 2) 1 100
    have := g_nonneg (3 * π / 2) 2 100
    have := g_nonneg (9 * π / 4) 3 100
    have := h_nonneg 100
    have := mul_nonneg (Real.exp_pos (π / 4)).le (g_nonneg (3 * π / 2) 1 100)
    have := mul_nonneg hxi0 (g_nonneg (3 * π / 2) 2 100)
    have := mul_nonneg hex0 (g_nonneg (9 * π / 4) 3 100)
    have := mul_nonneg (Real.exp_pos (Real.exp 1 * ξ)).le (h_nonneg 100)
    linarith
  have i1 : Real.exp (π / 4) * g (3 * π / 2) 1 100 + ξ * g (3 * π / 2) 2 100 +
      (Real.exp ξ - 1 - ξ) * g (9 * π / 4) 3 100 + Real.exp (Real.exp 1 * ξ) * h 100 ≤
      2.1933 * (10000 * 0.020762) + 4.3866 * (10000 * 0.023161) + 76 * (10000 * 0.0002146) +
        150860 * (10000 * 2.2e-11) := by linarith
  have q6 := mul_le_mul hB1.le i1 i0 (by norm_num)
  unfold R
  have := h_nonneg 100
  norm_num at q6 q5 h1 ⊢
  linarith

end R06.P1

theorem R06.proof_1 : R06.stmt_1 := by
  have hc : R06.crest = R06.P1.R 100 := by
    have := R06.P1.eqR 99
    rw [show (((99 + 1 : ℕ) : ℝ)) = 100 by norm_num] at this
    unfold R06.crest
    exact this
  refine ⟨fun k hk => ?_, hc ▸ R06.P1.R100⟩
  rw [R06.P1.eqR k, hc]
  apply R06.P1.R_le
  have : (100 : ℝ) ≤ ((k + 1 : ℕ) : ℝ) := by exact_mod_cast hk
  exact this

namespace R06.P3

open Real Q05

lemma sl : (1.77245385090551602 : ℝ) < √π := by
  rw [Real.lt_sqrt (by norm_num)]; nlinarith [Real.pi_gt_d20]

lemma sh : √π < 1.77245385090551603 := by
  rw [Real.sqrt_lt' (by norm_num)]; nlinarith [Real.pi_lt_d20]

end R06.P3

theorem R06.proof_2 : R06.stmt_2 := by
  obtain ⟨h1, h2⟩ := R06.proof_1
  refine ⟨fun k hk => ?_, h2⟩
  have h : |Q05.M k - Q05.Mmain k| ≤ R06.restBound k := Q07.proof_4 k hk
  set N : ℝ := ((k + 1 : ℕ) : ℝ)
  have hN : 0 < N := by positivity
  have hE : 0 < Real.exp (√(Real.pi * N)) := Real.exp_pos _
  have hA : 0 < N ^ 2 * Real.exp (√(Real.pi * N)) := by positivity
  have h3 : N ^ 2 * Real.exp (√(Real.pi * N)) * |Q05.M k - Q05.Mmain k| ≤ R06.crest :=
    (mul_le_mul_of_nonneg_left h hA.le).trans (h1 k hk)
  refine ⟨h3, ?_⟩
  have e : R06.crest / N ^ 2 * Real.exp (-√(Real.pi * N)) =
      R06.crest / (N ^ 2 * Real.exp (√(Real.pi * N))) := by
    rw [Real.exp_neg]; field_simp
  rw [e, le_div_iff₀ hA]
  linarith

theorem R06.proof_3 : R06.stmt_3 := by
  have hpl := Real.pi_gt_d20
  have hph := Real.pi_lt_d20
  have hsl := R06.P3.sl
  have hsh := R06.P3.sh
  have hs : 0 < √Real.pi := by linarith
  have hp2l := mul_lt_mul'' hpl hpl (by norm_num) (by norm_num)
  have hp2h := mul_lt_mul'' hph hph (by linarith) (by linarith)
  refine ⟨rfl, rfl, ?_, ?_, ?_⟩
  · show |(24 + 12 * Real.pi + Real.pi ^ 2) / (48 * √Real.pi) - 0.84121508777| < 5e-12
    rw [abs_lt]; constructor
    · rw [lt_sub_iff_add_lt, lt_div_iff₀ (by positivity)]; nlinarith
    · rw [sub_lt_iff_lt_add, div_lt_iff₀ (by positivity)]; nlinarith
  · show |(24 - 12 * Real.pi + Real.pi ^ 2) / (48 * √Real.pi) - (-0.04501183769)| < 5e-12
    rw [abs_lt]; constructor
    · rw [lt_sub_iff_add_lt, lt_div_iff₀ (by positivity)]; nlinarith
    · rw [sub_lt_iff_lt_add, div_lt_iff₀ (by positivity)]; nlinarith
  · intro N hN
    unfold Q05.X
    congr 1
    rw [mul_comm Complex.I, Complex.mul_re, Complex.exp_ofReal_mul_I_re,
      Complex.exp_ofReal_mul_I_im, Complex.add_re, Complex.add_im, Complex.one_re,
      Complex.one_im, Complex.div_ofReal_re, Complex.div_ofReal_im]
    simp only [Q05.varpi1]; ring

theorem R06.proof_4 : R06.stmt_4 := by
  intro k hk
  obtain ⟨h2, hc⟩ := R06.proof_2
  have ha := (h2 k hk).2
  have hb := Q07.proof_6 k hk
  set N : ℝ := ((k + 1 : ℕ) : ℝ)
  set e := Real.exp (-√(Real.pi * N))
  have hN : 0 < N := by positivity
  have he : 0 < e := Real.exp_pos _
  have h3 : |Q05.M k + Q05.X N| ≤ |Q05.M k - Q05.Mmain k| + |Q05.Mmain k + Q05.X N| := by
    rw [show Q05.M k + Q05.X N = (Q05.M k - Q05.Mmain k) + (Q05.Mmain k + Q05.X N) by ring]
    exact abs_add_le _ _
  have h4 : R06.crest / N ^ 2 * e + 71000 / N ^ 2 * e ≤ 191000 / N ^ 2 * e := by
    rw [← add_mul, ← add_div]
    apply mul_le_mul_of_nonneg_right _ he.le
    apply div_le_div_of_nonneg_right _ (by positivity)
    linarith
  linarith

theorem R06.proof_5 : R06.stmt_5 := by
  intro k hk
  have h := R06.proof_4 k hk
  unfold Q05.X at h
  refine h.trans (mul_le_mul_of_nonneg_right ?_ (Real.exp_pos _).le)
  exact div_le_div_of_nonneg_right (by norm_num) (by positivity)
