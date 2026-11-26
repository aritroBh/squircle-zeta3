import Mathlib
/-!
# X06 -- the area of the graph of `z = x^n y^n` (Ahmed's Open Problem 4)

Target id: X06.
Paper label: `thm:xnyn` (paper/ieee/main.tex, lines 1876-1891; definition of `A(n)` at
line 1876).
Detailed proof: paper/ieee/proofs/sec_g3_xnyn.tex (setting lines 5-10, theorem lines 13-24,
"Functions used" lines 28-36, Lemma `g3:lem:decomp` lines 43-48, Lemma `g3:lem:Phi` lines
64-66, Lemma `g3:lem:Eg` lines 84-86, Lemma `g3:lem:Tk` lines 97-99 with the Mellin
transform at line 102-106, Lemma `g3:lem:major` lines 117-119, the three pieces of `C(n)`
lines 131-139, Lemma `g3:lem:D1` lines 141-143, Lemma `g3:lem:D2` lines 153-155, assembly
lines 161-214) and the fragment paper/ieee/proofs/frag_xnyn_R3.tex (the hand identification
of `R_3`, now merged into the section at lines 174-212).
The printed decimals of `beta_1`, `beta_0` are target T11 (Wave1/T11_Decimals.lean); this
file does not import it and defines its own copies of the constants in namespace `X06`.

Quoted verbatim (main.tex 1876, the definition):

  The area $A(n)=\iint_{[0,1]^2}\sqrt{1+z_x^2+z_y^2}\,dx\,dy$ tends to $1$ ...

Quoted verbatim (main.tex 1880-1891):

  \begin{theorem}[Ahmed's Open Problem~4$^\dagger$]\label{thm:xnyn}
  Let $\ell=\log n$ and $\lambda=\log2$. For real $n\ge2$,
  \[
  A(n)=1+\frac{\sqrt2}{n}+\frac{R(n)}{n^2}+\frac{R_3(n)}{n^3}+O\Bigl(\frac{(1+\log n)^4}{n^4}\Bigr),
  \]
  where $R(n)=-\frac12\ell^2+\beta_1\ell+\beta_0$ and $R_3(n)=\frac13\ell^3+b_3\ell^2+c_3\ell+d_3$, with
  \begin{gather*}
  \beta_1=\frac{\lambda}{2}-1,\qquad \beta_0=-(1+\sqrt2)+\frac{\lambda}{2}-\frac{\pi^2}{12}-\frac{\lambda^2}{8},\\
  b_3=\frac{1-\lambda}{2},\qquad c_3=1-\frac{\lambda}{2}+\frac{\lambda^2}{4}+\frac{\pi^2}{6},\\
  d_3=1+\frac{3\sqrt2}{2}-\frac{\lambda}{2}+\frac{\lambda^2}{8}-\frac{\lambda^3}{24}+\frac{\pi^2}{12}-\frac{\pi^2\lambda}{12}-\frac{\zeta(3)}{2}.
  \end{gather*}
  \end{theorem}

Quoted verbatim (sec_g3_xnyn.tex 6-10, setting):

  $A(n)=\iint_{[0,1]^2}\sqrt{1+z_x^2+z_y^2}\dd x\dd y$ is the area of the graph of $z=x^ny^n$
  over $[0,1]^2$, where $z_x=nx^{n-1}y^n$ and $z_y=nx^ny^{n-1}$. ...
  Throughout, $e=1/n$ (not Euler's constant), $N=\sqrt2\,n$, $\Lg=\log N$, $\lambda=\log2$ as
  in the statement, and $\ell=\log n$, so $\Lg=\ell+\frac\lambda2$.

Quoted verbatim (sec_g3_xnyn.tex 26): "The $O$-constant is absolute (its value, $45$, is given
at the end of the proof)."

Quoted verbatim (sec_g3_xnyn.tex 29-36, functions used):

  $f(w)=1+w-\sqrt{1+w^2}$ and $g(w)=f(w)/w$ for $w>0$; $\psi(w)=wf'(w)=...$ and
  $\chi(w)=w\psi'(w)$. For a function $\phi$ on $(0,\infty)$ and $k\ge0$,
  $T_k(\phi)(N)=\int_0^N\log^k\frac Nw\,\phi(w)\dd w$. ...
  $E_g(s)=\frac{\Gamma(\frac{1+s}2)\,\Gamma(1-\frac s2)}{\sqrt\pi}=\sum_{j\ge0}\eta_js^j\quad(|s|<1)$
  and $\tilde m_k=\sum_{j\le k}\eta_j$ enter the polynomials
  $\mathrm{poly}_k=\frac{\Lg^{k+1}}{k+1}+\sum_{i=0}^k\frac{k!}{(k-i)!}\,\tilde m_{i+1}\Lg^{k-i}\quad(k\ge0).$

Quoted verbatim (Lemma `g3:lem:decomp`, lines 44-47):

  For $n>1$: $A(n)-1=\dfrac{\Phi(n)}{1+\frac e2}-e^2\,C(n)$, where
  $\Phi(n)=\int_0^1u^{n-1}\sqrt{1+u^2}\dd u$ and
  $C(n)=\int_0^\infty e^{-e\sigma}\int_0^\sigma f(W)\dd\delta\dd\sigma,\qquad W=N e^{-\sigma}K,\quad K=e^{e\sigma/2}\sqrt{\cosh(e\delta)}.$

Quoted verbatim (Lemma `g3:lem:Phi`, line 65):

  For real $n\ge2$ there is an absolute constant $c_\Phi<2.63$ (see the proof) with
  $\Bigl|\dfrac{\Phi(n)}{1+e/2}-\Bigl(\sqrt2\,e-\sqrt2\,e^2+\frac{5\sqrt2}4e^3\Bigr)\Bigr|\le c_\Phi e^4$.
  (The proof ends: "the lemma holds with $c_\Phi=\frac{3\sqrt2}2+\frac12<2.63$".)

Quoted verbatim (Lemma `g3:lem:Eg`, line 85):

  $E_g(0)=1$, $\eta_1=-\lambda$, $\eta_2=\frac{\lambda^2}2+\frac{\zeta(2)}2$,
  $\eta_3=-\frac{\lambda^3}6-\frac{\lambda\zeta(2)}2-\frac{\zeta(3)}4$.

Quoted verbatim (proof of `g3:lem:Tk`, line 102-106, Mellin transform):

  For $0<\operatorname{Re}z<1$ let $M_g(z)=\int_0^\infty w^{z-1}g(w)\dd w$ ... So
  $M_g(z)=\frac{E_g(z)}{z(1-z)}$.

Quoted verbatim (Lemma `g3:lem:Tk`, line 98):

  For $N\ge1$ and $k\ge0$: $T_k(g)(N)=\mathrm{poly}_k+(-1)^k\frac{k!}{2N}+(-1)^{k+1}\vartheta_k\frac{k!}{8\cdot3^{k+1}N^3}$
  with $\vartheta_k\in[0,1]$.

Quoted verbatim (Lemma `g3:lem:major`, line 118):

  (a) $0\le\psi(w)\le m_0(w):=\min\bigl(w,\frac1{2w}\bigr)$. (b) $|\chi(w)|\le\psi(w)$.
  (c) For $\sigma\ge0$ put $m(\sigma)=m_0(Ne^{-\sigma})=\min(Ne^{-\sigma},\frac{e^\sigma}{2N})$
  and $I_k=\int_0^\infty\sigma^km(\sigma)\dd\sigma$. Then $I_k\le(k!+\frac12)(\Lg+1)^k$ for
  $N\ge1$, and $T_k(g)(N)\le(k!+1)(\Lg+1)^{k+1}$.

Quoted verbatim (lines 131-139, the three pieces):

  Put $W_0=Ne^{-\sigma}$ and $\kappa(\sigma,\delta)=\log K=\frac{e\sigma}2+\frac12\log\cosh(e\delta)$.
  ... $C=C_0+\Delta_1+\Delta_2,\qquad C_0=\int_0^\infty\sigma f(W_0)\dd\sigma,$
  $\Delta_1=\int_0^\infty e^{-e\sigma}\!\!\int_0^\sigma\bigl[f(W_0e^\kappa)-f(W_0)\bigr]\dd\delta\dd\sigma,\qquad
  \Delta_2=\int_0^\infty(e^{-e\sigma}-1)\sigma f(W_0)\dd\sigma .$
  All three converge absolutely ... in particular $C_0=T_1(g)(N)$ ...

Quoted verbatim (Lemma `g3:lem:D1`, line 142):
  $\Delta_1=e\,T_1(g)(N)+r_1$, $|r_1|\le\frac{13}{12}e^2I_3\le7.1\,e^2(\Lg+1)^3$.

Quoted verbatim (Lemma `g3:lem:D2`, line 154):
  $\Delta_2=-e\,T_2(g)(N)+r_3$ with $0\le r_3\le\frac{e^2}2T_3(g)(N)\le3.5\,e^2(\Lg+1)^4$.

Quoted verbatim (assembly, lines 164-170, 172, 174, 214):

  C(n)=\mathrm{poly}_1-\frac{e}{2\sqrt2}+e\,(\mathrm{poly}_1-\mathrm{poly}_2)+\rho_C,\qquad|\rho_C|\le12\,e^2(\Lg+1)^4 .
  A(n)-1=\sqrt2e+\bigl(-\sqrt2-\mathrm{poly}_1\bigr)e^2+\Bigl(\frac{5\sqrt2}4+\frac1{2\sqrt2}-\mathrm{poly}_1+\mathrm{poly}_2\Bigr)e^3+\rho,\\
  |\rho|\le\bigl(c_\Phi+12(\Lg+1)^4\bigr)e^4 .
  ... the last bracket is $\frac{3\sqrt2}2+\mathrm{poly}_2-\mathrm{poly}_1$ ...
  It remains to identify $R(n)=-\sqrt2-\mathrm{poly}_1$ and
  $R_3(n)=\frac{3\sqrt2}2+\mathrm{poly}_2-\mathrm{poly}_1$ with the closed forms, substituting
  $\Lg=\ell+\frac\lambda2$.
  This proves Theorem~\ref{g3:thm:xnyn}: for real $n\ge2$ the error is at most
  $45(1+\log n)^4/n^4$.

## Formalization choices

* `zeta(3)` is `z3 := (riemannZeta 3).re` (a real number; `riemannZeta 3` is real, and the
  link to the series is Mathlib's `zeta_nat_eq_tsum_of_gt_one`). `zeta(2)` is written as
  `pi^2/6` (as the paper does in the final constants).
* `A n` is the Bochner integral over `Icc 0 1 ×ˢ Icc 0 1` of
  `sqrt (1 + (n x^(n-1) y^n)^2 + (n x^n y^(n-1))^2)` (real powers `rpow`, bases in `[0,1]`).
  For `n > 1` this integrand is continuous on the compact square, so the integral is genuine;
  `stmt_1` also asserts integrability explicitly.
* All other integrals are stated together with their integrability (or are interval integrals
  of functions that are integrable, which `stmt_4`, `stmt_5`, `stmt_7`, `stmt_8` assert) so that
  no claim is true by the `0`-for-non-integrable convention.
* `e = 1/n` is written `1/n` (or `n⁻¹`); `N = sqrt 2 * n` is `NN n`; `Lg = log N`.
* `psi w := w * deriv f w`, `chi w := w * deriv psi w` (`f` and `psi` are smooth on `ℝ`, so
  `deriv` is the true derivative).
* The Taylor coefficients `eta j` of `E_g` at `0` are defined as `iteratedDeriv j Eg 0 / j!`
  with `Eg` the real function `Gamma((1+s)/2) Gamma(1-s/2)/sqrt pi` (smooth near `0`; the
  `Gamma` arguments are near `1/2` and `1`). `mt k = sum_{j<=k} eta j`, and `poly k L` is the
  paper's `poly_k` with `Lg` replaced by the variable `L`; `k!/(k-i)!` is `Nat.descFactorial k i`.
* `T k N` is the interval integral `∫ w in 0..N, log(N/w)^k * g w` (the value of the integrand
  at the single point `w = 0` is irrelevant).
* Hypotheses: the paper's ranges are used (`n > 1` for the decomposition, as printed; the
  lemmas on `Delta_1`, `Delta_2` print no range, and we use `n > 1`, the range in which `C(n)`
  is defined; `n ≥ 2` wherever the paper says so). No hypothesis is strengthened.
* The `O`-claim is stated twice: with the supplement's explicit constant `45` (`stmt_15`) and,
  as the paper's theorem, with an existential constant (`stmt_16`, the last statement).
* `stmt_14` is the order-`n^-2` truncation with error `O((1+log n)^3/n^3)`, a consequence of the
  theorem (as allowed by the job), placed before the full statement.

## Statements
* `stmt_1`  Lemma decomp: `A(n) - 1 = Phi(n)/(1+e/2) - e^2 C(n)` for `n > 1` (with integrability).
* `stmt_2`  Lemma Phi: explicit expansion of `Phi(n)/(1+e/2)` with `c_Phi = 3 sqrt2/2 + 1/2`.
* `stmt_3`  Lemma Eg: the Taylor coefficients `eta_0..eta_3` of `E_g` (brings in `zeta(3)`).
* `stmt_4`  Mellin transform of `g`: `M_g(z) = E_g(z)/(z(1-z))` for `0<z<1`.
* `stmt_5`  Lemma Tk: the moment identity `T_k(g)(N) = poly_k + ...` for `N ≥ 1`, all `k`.
* `stmt_6`  Lemma major (a),(b): `0 ≤ psi ≤ min(w, 1/(2w))` and `|chi| ≤ psi`.
* `stmt_7`  Lemma major (c): bounds on `I_k` and `T_k(g)`.
* `stmt_8`  The split `C = C_0 + Delta_1 + Delta_2`, with `C_0 = T_1(g)(N)`.
* `stmt_9`  Lemma D1.
* `stmt_10` Lemma D2.
* `stmt_11` Expansion `(g3:eq:C)` of `C(n)` for `n ≥ 2`.
* `stmt_12` The expansion of `A(n) - 1` in `e` with coefficients in `poly_1`, `poly_2` (`n ≥ 2`).
* `stmt_13` Identification of `-sqrt2 - poly_1` and `3 sqrt2/2 + poly_2 - poly_1` with the closed
  forms `R`, `R_3` at `Lg = l + lambda/2`.
* `stmt_14` `A(n) = 1 + sqrt2/n + R(n)/n^2 + O((1+log n)^3/n^3)` for real `n ≥ 2`.
* `stmt_15` The theorem with the explicit constant `45` (supplement, line 214).
* `stmt_16` The theorem `thm:xnyn`: an absolute constant `C` with error `≤ C (1+log n)^4/n^4`
  for all real `n ≥ 2`.

NOT COVERED:
* The sentence (main.tex 1876) on the "short elementary proof" of the two-term expansion with
  `O(log n log log n)` (covered in a stronger form by `stmt_14`).
* The numerical sentences: `d_3 = 2.4722752293...`, `D(10^5) = -76.74`, "at $n=10^5$ the theorem
  gives $R(n)+R_3(n)/n=-76.7411$", and the high-precision fit (main.tex 1893-1902); the
  decimals of `beta_1`, `beta_0` are target T11.
* The intermediate facts inside proofs that the paper does not state as lemmas (the closed
  form of `rho(n)` in Lemma Phi, the constants `varpi_i = H^(i)(0)`, the tail bounds in Lemma Tk).
  (The Mellin-transform remark of main.tex 1900 is `stmt_4`.)
-/
-- BEGIN STATEMENTS
namespace X06

open MeasureTheory Set

/-- `λ = log 2` (`thm:xnyn`, line 1881). -/
noncomputable def lam : ℝ := Real.log 2

/-- `ζ(3)`, as the real part of Mathlib's `riemannZeta 3` (which is real). -/
noncomputable def z3 : ℝ := (riemannZeta 3).re

/-- The area integrand `sqrt(1 + z_x^2 + z_y^2)` with `z_x = n x^(n-1) y^n`,
`z_y = n x^n y^(n-1)` (sec_g3_xnyn.tex line 6). -/
noncomputable def areaIntegrand (n : ℝ) (p : ℝ × ℝ) : ℝ :=
  Real.sqrt (1 + (n * p.1 ^ (n - 1) * p.2 ^ n) ^ 2 + (n * p.1 ^ n * p.2 ^ (n - 1)) ^ 2)

/-- `A(n) = ∬_{[0,1]^2} sqrt(1+z_x^2+z_y^2) dx dy` (main.tex line 1876). -/
noncomputable def A (n : ℝ) : ℝ :=
  ∫ p in Icc (0 : ℝ) 1 ×ˢ Icc (0 : ℝ) 1, areaIntegrand n p

/-- `f(w) = 1 + w - sqrt(1+w^2)` (line 29). -/
noncomputable def f (w : ℝ) : ℝ := 1 + w - Real.sqrt (1 + w ^ 2)

/-- `g(w) = f(w)/w` (line 29). -/
noncomputable def g (w : ℝ) : ℝ := f w / w

/-- `psi(w) = w f'(w)` (line 29). -/
noncomputable def psi (w : ℝ) : ℝ := w * deriv f w

/-- `chi(w) = w psi'(w)` (line 29). -/
noncomputable def chi (w : ℝ) : ℝ := w * deriv psi w

/-- `Phi(n) = ∫_0^1 u^(n-1) sqrt(1+u^2) du` (Lemma decomp). -/
noncomputable def Phi (n : ℝ) : ℝ := ∫ u in (0 : ℝ)..1, u ^ (n - 1) * Real.sqrt (1 + u ^ 2)

/-- `N = sqrt 2 * n` (line 10). -/
noncomputable def NN (n : ℝ) : ℝ := Real.sqrt 2 * n

/-- `W = N e^{-σ} K`, `K = e^{eσ/2} sqrt(cosh(eδ))`, `e = 1/n` (Lemma decomp). -/
noncomputable def W (n σ δ : ℝ) : ℝ :=
  NN n * Real.exp (-σ) * (Real.exp (σ / n / 2) * Real.sqrt (Real.cosh (δ / n)))

/-- The inner integrand of `C(n)`: `σ ↦ e^{-eσ} ∫_0^σ f(W) dδ`. -/
noncomputable def Cinner (n σ : ℝ) : ℝ :=
  Real.exp (-σ / n) * ∫ δ in (0 : ℝ)..σ, f (W n σ δ)

/-- `C(n) = ∫_0^∞ e^{-eσ} ∫_0^σ f(W) dδ dσ` (Lemma decomp). -/
noncomputable def Cint (n : ℝ) : ℝ := ∫ σ in Ioi (0 : ℝ), Cinner n σ

/-- `E_g(s) = Γ((1+s)/2) Γ(1-s/2)/sqrt π` (line 31), as a real function. -/
noncomputable def Eg (s : ℝ) : ℝ :=
  Real.Gamma ((1 + s) / 2) * Real.Gamma (1 - s / 2) / Real.sqrt Real.pi

/-- `η_j`, the `j`-th Taylor coefficient of `E_g` at `0` (line 31). -/
noncomputable def eta (j : ℕ) : ℝ := iteratedDeriv j Eg 0 / (j.factorial : ℝ)

/-- `m̃_k = Σ_{j≤k} η_j` (line 33). -/
noncomputable def mt (k : ℕ) : ℝ := ∑ j ∈ Finset.range (k + 1), eta j

/-- `poly_k = L^{k+1}/(k+1) + Σ_{i=0}^k k!/(k-i)! m̃_{i+1} L^{k-i}` with `L` in place of
`Lg = log N` (line 35). -/
noncomputable def poly (k : ℕ) (L : ℝ) : ℝ :=
  L ^ (k + 1) / ((k : ℝ) + 1) +
    ∑ i ∈ Finset.range (k + 1), (k.descFactorial i : ℝ) * mt (i + 1) * L ^ (k - i)

/-- `T_k(g)(N) = ∫_0^N log^k(N/w) g(w) dw` (line 29). -/
noncomputable def T (k : ℕ) (N : ℝ) : ℝ := ∫ w in (0 : ℝ)..N, Real.log (N / w) ^ k * g w

/-- `M_g(z) = ∫_0^∞ w^{z-1} g(w) dw` (line 102). -/
noncomputable def Mg (z : ℝ) : ℝ := ∫ w in Ioi (0 : ℝ), w ^ (z - 1) * g w

/-- `W_0 = N e^{-σ}` (line 131). -/
noncomputable def W0 (n σ : ℝ) : ℝ := NN n * Real.exp (-σ)

/-- `κ(σ,δ) = eσ/2 + (1/2) log cosh(eδ)` (line 131). -/
noncomputable def kappa (n σ δ : ℝ) : ℝ := σ / n / 2 + Real.log (Real.cosh (δ / n)) / 2

/-- `C_0 = ∫_0^∞ σ f(W_0) dσ` (line 133). -/
noncomputable def C0 (n : ℝ) : ℝ := ∫ σ in Ioi (0 : ℝ), σ * f (W0 n σ)

/-- The integrand of `Δ_1` (line 136). -/
noncomputable def D1inner (n σ : ℝ) : ℝ :=
  Real.exp (-σ / n) * ∫ δ in (0 : ℝ)..σ, (f (W0 n σ * Real.exp (kappa n σ δ)) - f (W0 n σ))

/-- `Δ_1` (line 136). -/
noncomputable def Delta1 (n : ℝ) : ℝ := ∫ σ in Ioi (0 : ℝ), D1inner n σ

/-- `Δ_2 = ∫_0^∞ (e^{-eσ} - 1) σ f(W_0) dσ` (line 137). -/
noncomputable def Delta2 (n : ℝ) : ℝ :=
  ∫ σ in Ioi (0 : ℝ), (Real.exp (-σ / n) - 1) * σ * f (W0 n σ)

/-- The majorant `m(σ) = min(N e^{-σ}, e^σ/(2N))` (Lemma major (c)). -/
noncomputable def mmaj (N σ : ℝ) : ℝ := min (N * Real.exp (-σ)) (Real.exp σ / (2 * N))

/-- `I_k = ∫_0^∞ σ^k m(σ) dσ` (Lemma major (c)). -/
noncomputable def Ik (k : ℕ) (N : ℝ) : ℝ := ∫ σ in Ioi (0 : ℝ), σ ^ k * mmaj N σ

/-- `β_1 = λ/2 - 1` (line 1887). -/
noncomputable def beta1 : ℝ := lam / 2 - 1

/-- `β_0 = -(1+√2) + λ/2 - π²/12 - λ²/8` (line 1887). -/
noncomputable def beta0 : ℝ :=
  -(1 + Real.sqrt 2) + lam / 2 - Real.pi ^ 2 / 12 - lam ^ 2 / 8

/-- `b_3 = (1-λ)/2` (line 1888). -/
noncomputable def b3 : ℝ := (1 - lam) / 2

/-- `c_3 = 1 - λ/2 + λ²/4 + π²/6` (line 1888). -/
noncomputable def c3 : ℝ := 1 - lam / 2 + lam ^ 2 / 4 + Real.pi ^ 2 / 6

/-- `d_3 = 1 + 3√2/2 - λ/2 + λ²/8 - λ³/24 + π²/12 - π²λ/12 - ζ(3)/2` (line 1889). -/
noncomputable def d3 : ℝ :=
  1 + 3 * Real.sqrt 2 / 2 - lam / 2 + lam ^ 2 / 8 - lam ^ 3 / 24 + Real.pi ^ 2 / 12
    - Real.pi ^ 2 * lam / 12 - z3 / 2

/-- `R(n) = -ℓ²/2 + β_1 ℓ + β_0`, `ℓ = log n` (line 1885). -/
noncomputable def R (n : ℝ) : ℝ :=
  -(Real.log n) ^ 2 / 2 + beta1 * Real.log n + beta0

/-- `R_3(n) = ℓ³/3 + b_3 ℓ² + c_3 ℓ + d_3`, `ℓ = log n` (line 1885). -/
noncomputable def R3 (n : ℝ) : ℝ :=
  (Real.log n) ^ 3 / 3 + b3 * (Real.log n) ^ 2 + c3 * Real.log n + d3

/-- Lemma `g3:lem:decomp`: for `n > 1`, `A(n) - 1 = Φ(n)/(1+e/2) - e² C(n)`, with the area
integrand integrable on the square and the outer integrand of `C(n)` integrable on `(0,∞)`
("the integrals converge absolutely"). -/
def stmt_1 : Prop :=
  ∀ n : ℝ, 1 < n →
    IntegrableOn (areaIntegrand n) (Icc (0 : ℝ) 1 ×ˢ Icc (0 : ℝ) 1) ∧
    IntegrableOn (Cinner n) (Ioi 0) ∧
    A n - 1 = Phi n / (1 + 1 / n / 2) - (1 / n) ^ 2 * Cint n

/-- Lemma `g3:lem:Phi` with the proof's constant `c_Φ = 3√2/2 + 1/2`: for real `n ≥ 2`,
`|Φ(n)/(1+e/2) - (√2 e - √2 e² + 5√2/4 e³)| ≤ c_Φ e⁴`. -/
def stmt_2 : Prop :=
  ∀ n : ℝ, 2 ≤ n →
    |Phi n / (1 + 1 / n / 2) -
        (Real.sqrt 2 / n - Real.sqrt 2 / n ^ 2 + 5 * Real.sqrt 2 / 4 / n ^ 3)|
      ≤ (3 * Real.sqrt 2 / 2 + 1 / 2) / n ^ 4

/-- Lemma `g3:lem:Eg`: `E_g(0) = 1`, `η_1 = -λ`, `η_2 = λ²/2 + ζ(2)/2`,
`η_3 = -λ³/6 - λζ(2)/2 - ζ(3)/4` (with `ζ(2) = π²/6`). -/
def stmt_3 : Prop :=
  eta 0 = 1 ∧ eta 1 = -lam ∧ eta 2 = lam ^ 2 / 2 + Real.pi ^ 2 / 12 ∧
    eta 3 = -lam ^ 3 / 6 - lam * Real.pi ^ 2 / 12 - z3 / 4

/-- Proof of Lemma `g3:lem:Tk`, Mellin transform of `g`: for `0 < z < 1`,
`∫_0^∞ w^{z-1} g(w) dw` converges and equals `E_g(z)/(z(1-z))`. -/
def stmt_4 : Prop :=
  ∀ z : ℝ, 0 < z → z < 1 →
    IntegrableOn (fun w : ℝ => w ^ (z - 1) * g w) (Ioi 0) ∧ Mg z = Eg z / (z * (1 - z))

/-- Lemma `g3:lem:Tk`: for `N ≥ 1` and `k ≥ 0`, the integrand of `T_k(g)(N)` is integrable on
`(0,N]` and `T_k(g)(N) = poly_k + (-1)^k k!/(2N) + (-1)^{k+1} ϑ_k k!/(8·3^{k+1} N³)` for some
`ϑ_k ∈ [0,1]`, with `poly_k` evaluated at `Lg = log N`. -/
def stmt_5 : Prop :=
  ∀ (k : ℕ) (N : ℝ), 1 ≤ N →
    IntegrableOn (fun w : ℝ => Real.log (N / w) ^ k * g w) (Ioc 0 N) ∧
    ∃ θ : ℝ, θ ∈ Icc (0 : ℝ) 1 ∧
      T k N = poly k (Real.log N) + (-1 : ℝ) ^ k * (k.factorial : ℝ) / (2 * N)
        + (-1 : ℝ) ^ (k + 1) * θ * (k.factorial : ℝ) / (8 * 3 ^ (k + 1) * N ^ 3)

/-- Lemma `g3:lem:major` (a), (b): for `w > 0`, `0 ≤ ψ(w) ≤ min(w, 1/(2w))` and
`|χ(w)| ≤ ψ(w)`. -/
def stmt_6 : Prop :=
  ∀ w : ℝ, 0 < w → 0 ≤ psi w ∧ psi w ≤ min w (1 / (2 * w)) ∧ |chi w| ≤ psi w

/-- Lemma `g3:lem:major` (c): for `N ≥ 1` and `k ≥ 0`, `σ^k m(σ)` is integrable on `(0,∞)`,
`I_k ≤ (k! + 1/2)(Lg+1)^k` and `T_k(g)(N) ≤ (k!+1)(Lg+1)^{k+1}`, `Lg = log N`. -/
def stmt_7 : Prop :=
  ∀ (k : ℕ) (N : ℝ), 1 ≤ N →
    IntegrableOn (fun σ : ℝ => σ ^ k * mmaj N σ) (Ioi 0) ∧
    Ik k N ≤ ((k.factorial : ℝ) + 1 / 2) * (Real.log N + 1) ^ k ∧
    T k N ≤ ((k.factorial : ℝ) + 1) * (Real.log N + 1) ^ (k + 1)

/-- The three pieces of `C(n)` (lines 131-139): for `n > 1` the integrands of `C_0`, `Δ_1`,
`Δ_2` are integrable on `(0,∞)`, `C = C_0 + Δ_1 + Δ_2`, and `C_0 = T_1(g)(N)`. -/
def stmt_8 : Prop :=
  ∀ n : ℝ, 1 < n →
    IntegrableOn (fun σ : ℝ => σ * f (W0 n σ)) (Ioi 0) ∧
    IntegrableOn (D1inner n) (Ioi 0) ∧
    IntegrableOn (fun σ : ℝ => (Real.exp (-σ / n) - 1) * σ * f (W0 n σ)) (Ioi 0) ∧
    Cint n = C0 n + Delta1 n + Delta2 n ∧
    C0 n = T 1 (NN n)

/-- Lemma `g3:lem:D1`: `Δ_1 = e T_1(g)(N) + r_1` with
`|r_1| ≤ (13/12) e² I_3 ≤ 7.1 e² (Lg+1)³`. -/
def stmt_9 : Prop :=
  ∀ n : ℝ, 1 < n →
    |Delta1 n - (1 / n) * T 1 (NN n)| ≤ 13 / 12 * (1 / n) ^ 2 * Ik 3 (NN n) ∧
    13 / 12 * (1 / n) ^ 2 * Ik 3 (NN n) ≤ 7.1 * (1 / n) ^ 2 * (Real.log (NN n) + 1) ^ 3

/-- Lemma `g3:lem:D2`: `Δ_2 = -e T_2(g)(N) + r_3` with
`0 ≤ r_3 ≤ (e²/2) T_3(g)(N) ≤ 3.5 e² (Lg+1)⁴`. -/
def stmt_10 : Prop :=
  ∀ n : ℝ, 1 < n →
    0 ≤ Delta2 n + (1 / n) * T 2 (NN n) ∧
    Delta2 n + (1 / n) * T 2 (NN n) ≤ (1 / n) ^ 2 / 2 * T 3 (NN n) ∧
    (1 / n) ^ 2 / 2 * T 3 (NN n) ≤ 3.5 * (1 / n) ^ 2 * (Real.log (NN n) + 1) ^ 4

/-- Equation `(g3:eq:C)`: for `n ≥ 2`,
`C(n) = poly_1 - e/(2√2) + e(poly_1 - poly_2) + ρ_C`, `|ρ_C| ≤ 12 e² (Lg+1)⁴`. -/
def stmt_11 : Prop :=
  ∀ n : ℝ, 2 ≤ n →
    |Cint n - (poly 1 (Real.log (NN n)) - (1 / n) / (2 * Real.sqrt 2)
        + (1 / n) * (poly 1 (Real.log (NN n)) - poly 2 (Real.log (NN n))))|
      ≤ 12 * (1 / n) ^ 2 * (Real.log (NN n) + 1) ^ 4

/-- The assembled expansion (lines 168-172): for `n ≥ 2`,
`A(n) - 1 = √2 e + (-√2 - poly_1) e² + (3√2/2 + poly_2 - poly_1) e³ + ρ`,
`|ρ| ≤ (c_Φ + 12 (Lg+1)⁴) e⁴` with `c_Φ = 3√2/2 + 1/2`. -/
def stmt_12 : Prop :=
  ∀ n : ℝ, 2 ≤ n →
    |A n - 1 - (Real.sqrt 2 * (1 / n) + (-Real.sqrt 2 - poly 1 (Real.log (NN n))) * (1 / n) ^ 2
        + (3 * Real.sqrt 2 / 2 + poly 2 (Real.log (NN n)) - poly 1 (Real.log (NN n)))
          * (1 / n) ^ 3)|
      ≤ (3 * Real.sqrt 2 / 2 + 1 / 2 + 12 * (Real.log (NN n) + 1) ^ 4) * (1 / n) ^ 4

/-- Identification of the coefficients (lines 174-212): with `Lg = ℓ + λ/2`,
`-√2 - poly_1 = -ℓ²/2 + β_1 ℓ + β_0` and `3√2/2 + poly_2 - poly_1 = ℓ³/3 + b_3 ℓ² + c_3 ℓ + d_3`,
for every real `ℓ`. -/
def stmt_13 : Prop :=
  ∀ l : ℝ,
    -Real.sqrt 2 - poly 1 (l + lam / 2) = -l ^ 2 / 2 + beta1 * l + beta0 ∧
    3 * Real.sqrt 2 / 2 + poly 2 (l + lam / 2) - poly 1 (l + lam / 2)
      = l ^ 3 / 3 + b3 * l ^ 2 + c3 * l + d3

/-- Order-`n⁻²` truncation of `thm:xnyn`: there is `C` with
`|A(n) - 1 - √2/n - R(n)/n²| ≤ C (1+log n)³/n³` for all real `n ≥ 2`. -/
def stmt_14 : Prop :=
  ∃ C : ℝ, ∀ n : ℝ, 2 ≤ n →
    |A n - (1 + Real.sqrt 2 / n + R n / n ^ 2)| ≤ C * (1 + Real.log n) ^ 3 / n ^ 3

/-- `thm:xnyn` with the supplement's explicit constant (sec_g3_xnyn.tex line 214): for real
`n ≥ 2`, `|A(n) - (1 + √2/n + R(n)/n² + R_3(n)/n³)| ≤ 45 (1+log n)⁴/n⁴`. -/
def stmt_15 : Prop :=
  ∀ n : ℝ, 2 ≤ n →
    |A n - (1 + Real.sqrt 2 / n + R n / n ^ 2 + R3 n / n ^ 3)|
      ≤ 45 * (1 + Real.log n) ^ 4 / n ^ 4

/-- `thm:xnyn` (main.tex 1880-1891): for real `n ≥ 2`,
`A(n) = 1 + √2/n + R(n)/n² + R_3(n)/n³ + O((1+log n)⁴/n⁴)`, the `O` with an absolute constant. -/
def stmt_16 : Prop :=
  ∃ C : ℝ, ∀ n : ℝ, 2 ≤ n →
    |A n - (1 + Real.sqrt 2 / n + R n / n ^ 2 + R3 n / n ^ 3)|
      ≤ C * (1 + Real.log n) ^ 4 / n ^ 4

end X06
-- END STATEMENTS

open MeasureTheory Set

/-! All sixteen statements are proved below. `stmt_1` (exact decomposition) is at the end,
before the assembly; `stmt_3` (Taylor coefficients of `E_g`, via the log-Gamma series copied from
`X05`), `stmt_4` (Mellin transform) and `stmt_5` (moment identity) follow the earlier lemmas. -/

namespace X06
open Real

lemma int_mono (n : ℝ) (hn : 1 < n) (k : ℕ) :
    ∫ u in (0:ℝ)..1, u ^ (n - 1) * u ^ k = 1 / (n + k) := by
  have h : ∫ u in (0:ℝ)..1, u ^ (n - 1) * u ^ k = ∫ u in (0:ℝ)..1, u ^ (n - 1 + k) := by
    refine intervalIntegral.integral_congr (fun u hu => ?_)
    have hu0 : 0 ≤ u := by
      rw [Set.uIcc_of_le zero_le_one] at hu; exact hu.1
    rw [Real.rpow_add' hu0 (by have := (Nat.cast_nonneg k : (0:ℝ) ≤ k); intro h'; linarith),
      Real.rpow_natCast]
  rw [h, integral_rpow (Or.inl (by linarith [(Nat.cast_nonneg k : (0:ℝ) ≤ k)]))]
  have : n - 1 + (k:ℝ) + 1 = n + k := by ring
  rw [this]; simp [Real.zero_rpow (by positivity : (n + (k:ℝ)) ≠ 0)]

lemma cont_rp (n : ℝ) (hn : 1 ≤ n) : Continuous (fun u : ℝ => u ^ (n - 1)) :=
  Real.continuous_rpow_const (by linarith)

lemma ii_mono (n : ℝ) (hn : 1 ≤ n) (k : ℕ) (c : ℝ) :
    IntervalIntegrable (fun u : ℝ => c * (u ^ (n - 1) * u ^ k)) MeasureTheory.volume 0 1 :=
  (continuous_const.mul ((cont_rp n hn).mul (continuous_pow k))).intervalIntegrable 0 1

/-- the cubic Taylor polynomial of `sqrt(1+u^2)` at `u = 1`. -/
noncomputable def Pp (u : ℝ) : ℝ := Real.sqrt 2 * (11 + u + 5 * u ^ 2 - u ^ 3) / 16

lemma Pp_err (u : ℝ) (h0 : 0 ≤ u) (h1 : u ≤ 1) :
    0 ≤ Real.sqrt (1 + u ^ 2) - Pp u ∧
      Real.sqrt (1 + u ^ 2) - Pp u ≤ 3 / 32 * (1 - u) ^ 4 := by
  set r := Real.sqrt (1 + u ^ 2)
  have hr : r ^ 2 = 1 + u ^ 2 := Real.sq_sqrt (by positivity)
  have hr1 : 1 ≤ r := by rw [show (1:ℝ) = Real.sqrt 1 by simp]; exact Real.sqrt_le_sqrt (by nlinarith)
  have hs2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have hs0 : 0 < Real.sqrt 2 := by positivity
  have hP : 0 ≤ Pp u := by
    unfold Pp
    have : 0 ≤ 11 + u + 5 * u ^ 2 - u ^ 3 := by nlinarith
    positivity
  set t := 1 - u
  have key : (r - Pp u) * (r + Pp u) = t ^ 4 * (12 - 4 * t - t ^ 2) / 128 := by
    unfold Pp
    linear_combination hr - ((11 + u + 5 * u ^ 2 - u ^ 3) ^ 2 / 256) * hs2
  have ht0 : 0 ≤ t := by simp [t]; linarith
  have ht1 : t ≤ 1 := by simp [t]; linarith
  have hnum : 0 ≤ t ^ 4 * (12 - 4 * t - t ^ 2) / 128 := by
    have : 0 ≤ 12 - 4 * t - t ^ 2 := by nlinarith
    positivity
  have hpos : 0 ≤ r - Pp u := by
    by_contra hc
    push Not at hc
    nlinarith
  refine ⟨hpos, ?_⟩
  have : r - Pp u ≤ (r - Pp u) * (r + Pp u) := by nlinarith
  have h4 : 0 ≤ t ^ 4 := by positivity
  nlinarith

end X06

theorem X06.proof_2 : X06.stmt_2 := by
  intro n hn
  have hn1 : (1:ℝ) ≤ n := by linarith
  have hn0 : (0:ℝ) < n := by linarith
  set s := Real.sqrt 2 with hsdef
  have hs0 : 0 < s := by positivity
  -- Phi = main + E
  have hcr : Continuous (fun u : ℝ => u ^ (n - 1) * Real.sqrt (1 + u ^ 2)) :=
    (X06.cont_rp n hn1).mul (by fun_prop)
  have hcP : Continuous (fun u : ℝ => u ^ (n - 1) * X06.Pp u) :=
    (X06.cont_rp n hn1).mul (by unfold X06.Pp; fun_prop)
  have hcE : Continuous (fun u : ℝ => u ^ (n - 1) * (Real.sqrt (1 + u ^ 2) - X06.Pp u)) :=
    (X06.cont_rp n hn1).mul (by unfold X06.Pp; fun_prop)
  set E := ∫ u in (0:ℝ)..1, u ^ (n - 1) * (Real.sqrt (1 + u ^ 2) - X06.Pp u) with hE
  have hsplit : X06.Phi n = (∫ u in (0:ℝ)..1, u ^ (n - 1) * X06.Pp u) + E := by
    rw [hE, ← intervalIntegral.integral_add (hcP.intervalIntegrable 0 1)
      (hcE.intervalIntegrable 0 1)]
    unfold X06.Phi
    congr 1; funext u; ring
  have hM : (∫ u in (0:ℝ)..1, u ^ (n - 1) * X06.Pp u)
      = s * (11 / n + 1 / (n + 1) + 5 / (n + 2) - 1 / (n + 3)) / 16 := by
    have e : (fun u : ℝ => u ^ (n - 1) * X06.Pp u) = fun u =>
        ((s * 11 / 16) * (u ^ (n - 1) * u ^ 0) + (s / 16) * (u ^ (n - 1) * u ^ 1))
        + ((s * 5 / 16) * (u ^ (n - 1) * u ^ 2) + (-(s / 16)) * (u ^ (n - 1) * u ^ 3)) := by
      funext u; unfold X06.Pp; ring
    rw [e, intervalIntegral.integral_add ((X06.ii_mono n hn1 0 _).add (X06.ii_mono n hn1 1 _))
      ((X06.ii_mono n hn1 2 _).add (X06.ii_mono n hn1 3 _)),
      intervalIntegral.integral_add (X06.ii_mono n hn1 0 _) (X06.ii_mono n hn1 1 _),
      intervalIntegral.integral_add (X06.ii_mono n hn1 2 _) (X06.ii_mono n hn1 3 _)]
    simp only [intervalIntegral.integral_const_mul, X06.int_mono n (by linarith)]
    push_cast
    field_simp
    ring
  have hE0 : 0 ≤ E := by
    rw [hE]
    refine intervalIntegral.integral_nonneg zero_le_one (fun u hu => ?_)
    exact mul_nonneg (Real.rpow_nonneg hu.1 _) (X06.Pp_err u hu.1 hu.2).1
  have hE1 : E ≤ 9 / 4 / (n * (n + 1) * (n + 2) * (n + 3) * (n + 4)) := by
    have hb : E ≤ ∫ u in (0:ℝ)..1, u ^ (n - 1) * (3 / 32 * (1 - u) ^ 4) := by
      rw [hE]
      refine intervalIntegral.integral_mono_on zero_le_one (hcE.intervalIntegrable 0 1)
        (((X06.cont_rp n hn1).mul (by fun_prop)).intervalIntegrable 0 1) (fun u hu => ?_)
      exact mul_le_mul_of_nonneg_left (X06.Pp_err u hu.1 hu.2).2 (Real.rpow_nonneg hu.1 _)
    have hv : (∫ u in (0:ℝ)..1, u ^ (n - 1) * (3 / 32 * (1 - u) ^ 4))
        = 9 / 4 / (n * (n + 1) * (n + 2) * (n + 3) * (n + 4)) := by
      have e : (fun u : ℝ => u ^ (n - 1) * (3 / 32 * (1 - u) ^ 4)) = fun u =>
          (((3/32) * (u ^ (n - 1) * u ^ 0) + (-(3/8)) * (u ^ (n - 1) * u ^ 1))
          + ((9/16) * (u ^ (n - 1) * u ^ 2) + (-(3/8)) * (u ^ (n - 1) * u ^ 3)))
          + (3/32) * (u ^ (n - 1) * u ^ 4) := by
        funext u; ring
      rw [e, intervalIntegral.integral_add (((X06.ii_mono n hn1 0 _).add (X06.ii_mono n hn1 1 _)).add
          ((X06.ii_mono n hn1 2 _).add (X06.ii_mono n hn1 3 _))) (X06.ii_mono n hn1 4 _),
        intervalIntegral.integral_add ((X06.ii_mono n hn1 0 _).add (X06.ii_mono n hn1 1 _))
          ((X06.ii_mono n hn1 2 _).add (X06.ii_mono n hn1 3 _)),
        intervalIntegral.integral_add (X06.ii_mono n hn1 0 _) (X06.ii_mono n hn1 1 _),
        intervalIntegral.integral_add (X06.ii_mono n hn1 2 _) (X06.ii_mono n hn1 3 _)]
      simp only [intervalIntegral.integral_const_mul, X06.int_mono n (by linarith)]
      push_cast
      field_simp
      ring
    linarith
  rw [hsplit, hM]
  have hd : (0:ℝ) < 1 + 1 / n / 2 := by positivity
  have hid : (s * (11 / n + 1 / (n + 1) + 5 / (n + 2) - 1 / (n + 3)) / 16 + E) / (1 + 1 / n / 2)
      - (s / n - s / n ^ 2 + 5 * s / 4 / n ^ 3)
      = -(s * n * (12 * n ^ 3 + 72 * n ^ 2 + 91 * n + 30)
          / (4 * (n + 1) * (n + 2) * (n + 3) * (2 * n + 1))) / n ^ 4
        + E / (1 + 1 / n / 2) := by
    field_simp
    ring
  rw [hid]
  have hA : s * n * (12 * n ^ 3 + 72 * n ^ 2 + 91 * n + 30)
          / (4 * (n + 1) * (n + 2) * (n + 3) * (2 * n + 1)) ≤ 3 * s / 2 := by
    rw [div_le_iff₀ (by positivity)]
    have : n * (12 * n ^ 3 + 72 * n ^ 2 + 91 * n + 30)
        ≤ 6 * ((n + 1) * (n + 2) * (n + 3) * (2 * n + 1)) := by nlinarith
    nlinarith
  have hA0 : 0 ≤ s * n * (12 * n ^ 3 + 72 * n ^ 2 + 91 * n + 30)
          / (4 * (n + 1) * (n + 2) * (n + 3) * (2 * n + 1)) := by positivity
  have hB0 : 0 ≤ E / (1 + 1 / n / 2) := div_nonneg hE0 hd.le
  have hB : E / (1 + 1 / n / 2) ≤ 1 / 2 / n ^ 4 := by
    have h1 : E / (1 + 1 / n / 2) ≤ E := div_le_self hE0 (by
      have : 0 ≤ 1 / n / 2 := by positivity
      linarith)
    have h2 : 9 / 4 / (n * (n + 1) * (n + 2) * (n + 3) * (n + 4)) ≤ 1 / 2 / n ^ 4 := by
      rw [div_le_div_iff₀ (by positivity) (by positivity)]
      have : 6 * n ^ 3 ≤ (n + 1) * (n + 2) * (n + 3) * (n + 4) := by nlinarith
      nlinarith
    linarith
  have hn4 : 0 < n ^ 4 := by positivity
  rw [abs_le]
  constructor
  · have : -(s * n * (12 * n ^ 3 + 72 * n ^ 2 + 91 * n + 30)
          / (4 * (n + 1) * (n + 2) * (n + 3) * (2 * n + 1))) / n ^ 4 ≥ -(3 * s / 2) / n ^ 4 := by
      apply div_le_div_of_nonneg_right _ hn4.le; linarith
    have e2 : (3 * s / 2 + 1 / 2) / n ^ 4 = 3 * s / 2 / n ^ 4 + 1 / 2 / n ^ 4 := by ring
    have e3 : -(3 * s / 2) / n ^ 4 = -(3 * s / 2 / n ^ 4) := by ring
    have : 0 ≤ 1 / 2 / n ^ 4 := by positivity
    linarith
  · have : -(s * n * (12 * n ^ 3 + 72 * n ^ 2 + 91 * n + 30)
          / (4 * (n + 1) * (n + 2) * (n + 3) * (2 * n + 1))) / n ^ 4 ≤ 0 :=
      div_nonpos_of_nonpos_of_nonneg (by linarith) hn4.le
    have e2 : (3 * s / 2 + 1 / 2) / n ^ 4 = 3 * s / 2 / n ^ 4 + 1 / 2 / n ^ 4 := by ring
    have : 0 ≤ 3 * s / 2 / n ^ 4 := by positivity
    linarith

namespace X06
open Real

lemma hasDerivAt_f (w : ℝ) :
    HasDerivAt f (1 - w / Real.sqrt (1 + w ^ 2)) w := by
  have h1 : (1 + w ^ 2) ≠ 0 := by positivity
  have hs : 0 < Real.sqrt (1 + w ^ 2) := Real.sqrt_pos.2 (by positivity)
  have h := ((hasDerivAt_id w).const_add 1).sub
      (((hasDerivAt_pow 2 w).const_add 1).sqrt h1)
  convert h using 1
  · funext x; simp [f]
  · simp only [Nat.cast_ofNat]; field_simp; ring

lemma deriv_f_eq : deriv f = fun w => 1 - w / Real.sqrt (1 + w ^ 2) := by
  funext w; exact (hasDerivAt_f w).deriv

lemma psi_eq (w : ℝ) : psi w = w - w ^ 2 / Real.sqrt (1 + w ^ 2) := by
  unfold psi; rw [deriv_f_eq]; ring

lemma hasDerivAt_psi (w : ℝ) :
    HasDerivAt psi (1 - (2 * w + w ^ 3) / Real.sqrt (1 + w ^ 2) ^ 3) w := by
  have h1 : (1 + w ^ 2) ≠ 0 := by positivity
  have hs : 0 < Real.sqrt (1 + w ^ 2) := Real.sqrt_pos.2 (by positivity)
  have hss : Real.sqrt (1 + w ^ 2) ^ 2 = 1 + w ^ 2 := Real.sq_sqrt (by positivity)
  have hpsi : psi = fun w => w - w ^ 2 / Real.sqrt (1 + w ^ 2) := funext psi_eq
  rw [hpsi]
  have h := (hasDerivAt_id w).sub ((hasDerivAt_pow 2 w).div
    (((hasDerivAt_pow 2 w).const_add 1).sqrt h1) hs.ne')
  convert h using 1
  · funext x; simp
  simp only [Nat.cast_ofNat]
  field_simp
  rw [show Real.sqrt (1 + w ^ 2) ^ 3 = Real.sqrt (1 + w ^ 2) * Real.sqrt (1 + w ^ 2) ^ 2 by ring,
    hss]
  ring

lemma chi_eq (w : ℝ) : chi w = w - (2 * w ^ 2 + w ^ 4) / Real.sqrt (1 + w ^ 2) ^ 3 := by
  unfold chi; rw [(hasDerivAt_psi w).deriv]; ring

end X06

theorem X06.proof_6 : X06.stmt_6 := by
  intro w hw
  rw [X06.psi_eq, X06.chi_eq]
  set s := Real.sqrt (1 + w ^ 2) with hsdef
  have hs : 0 < s := Real.sqrt_pos.2 (by positivity)
  have hss : s ^ 2 = 1 + w ^ 2 := Real.sq_sqrt (by positivity)
  have hws : w < s := by nlinarith
  set T := w / s with hT
  have hT0 : 0 < T := div_pos hw hs
  have hT1 : T < 1 := (div_lt_one hs).2 hws
  have hpsi : w - w ^ 2 / s = w * (1 - T) := by rw [hT]; field_simp
  have hchi : w - (2 * w ^ 2 + w ^ 4) / s ^ 3 = w * (1 - T) * (1 - T - T ^ 2) := by
    rw [hT]; field_simp
    nlinarith [hss]
  rw [hpsi, hchi]
  refine ⟨by nlinarith, ?_, ?_⟩
  · refine le_min (by nlinarith) ?_
    rw [le_div_iff₀ (by positivity)]
    have : w * (1 - T) * (2 * w) = 2 * w ^ 2 / (s * (s + w)) := by
      rw [hT]; field_simp; nlinarith [hss]
    rw [this, div_le_one (by positivity)]
    nlinarith
  · rw [abs_mul, abs_of_nonneg (by nlinarith : 0 ≤ w * (1 - T))]
    have : |1 - T - T ^ 2| ≤ 1 := abs_le.2 ⟨by nlinarith, by nlinarith⟩
    nlinarith [this, (by nlinarith : 0 ≤ w * (1 - T))]

namespace X06
open Real MeasureTheory Set

lemma gammaMoment (i : ℕ) :
    IntegrableOn (fun t : ℝ => t ^ i * Real.exp (-t)) (Ioi 0) ∧
      ∫ t in Ioi (0:ℝ), t ^ i * Real.exp (-t) = (i.factorial : ℝ) := by
  have hk : (0 : ℝ) < ↑i + 1 := by positivity
  have key := integral_rpow_mul_exp_neg_mul_Ioi hk one_pos
  have e : ∀ t ∈ Ioi (0:ℝ), t ^ ((↑i + 1 : ℝ) - 1) * Real.exp (-(1 * t)) = t ^ i * Real.exp (-t) := by
    intro t ht; simp [add_sub_cancel_right, Real.rpow_natCast]
  rw [setIntegral_congr_fun measurableSet_Ioi e] at key
  rw [Real.Gamma_nat_eq_factorial] at key
  simp only [div_one, Real.one_rpow, one_mul] at key
  refine ⟨?_, key⟩
  exact .of_integral_ne_zero (by rw [key]; positivity)

lemma shift_Ioi (F : ℝ → ℝ) (L : ℝ) :
    (IntegrableOn F (Ioi L) ↔ IntegrableOn (fun t => F (t + L)) (Ioi 0)) ∧
      ∫ σ in Ioi L, F σ = ∫ t in Ioi (0:ℝ), F (t + L) := by
  have himg : (fun t : ℝ => t + L) '' Ioi 0 = Ioi L := by
    rw [Set.image_add_const_Ioi]; simp
  have hd : ∀ x ∈ Ioi (0:ℝ), HasDerivWithinAt (fun t : ℝ => t + L) 1 (Ioi 0) x :=
    fun x _ => ((hasDerivAt_id x).add_const L).hasDerivWithinAt
  have hinj : InjOn (fun t : ℝ => t + L) (Ioi 0) := fun x _ y _ h => by simpa using h
  constructor
  · rw [← himg, integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ioi hd hinj]
    simp
  · rw [← himg, integral_image_eq_integral_abs_deriv_smul measurableSet_Ioi hd hinj]
    simp

lemma tail_moment (k : ℕ) (L : ℝ) (hL : 0 ≤ L) :
    IntegrableOn (fun σ : ℝ => σ ^ k * Real.exp (-σ)) (Ioi L) ∧
      ∫ σ in Ioi L, σ ^ k * Real.exp (-σ) ≤ Real.exp (-L) * (k.factorial : ℝ) * (L + 1) ^ k := by
  refine ⟨(gammaMoment k).1.mono_set (Ioi_subset_Ioi hL), ?_⟩
  rw [(shift_Ioi _ L).2]
  have hexp : ∀ t : ℝ, (t + L) ^ k * Real.exp (-(t + L))
      = ∑ m ∈ Finset.range (k + 1),
          (Real.exp (-L) * L ^ (k - m) * (k.choose m : ℝ)) * (t ^ m * Real.exp (-t)) := by
    intro t
    rw [add_pow, Finset.sum_mul]
    refine Finset.sum_congr rfl (fun m _ => ?_)
    rw [show -(t + L) = -t + -L by ring, Real.exp_add]; ring
  simp_rw [hexp]
  rw [integral_finsetSum _ (fun m _ => ((gammaMoment m).1).const_mul _)]
  simp_rw [integral_const_mul, (gammaMoment _).2]
  have hb : ∀ m ∈ Finset.range (k + 1),
      Real.exp (-L) * L ^ (k - m) * (k.choose m : ℝ) * (m.factorial : ℝ)
        ≤ Real.exp (-L) * (k.factorial : ℝ) * (1 ^ m * L ^ (k - m) * (k.choose m : ℝ)) := by
    intro m hm
    have hmk : m ≤ k := Nat.lt_succ_iff.1 (Finset.mem_range.1 hm)
    have : (m.factorial : ℝ) ≤ k.factorial := by exact_mod_cast Nat.factorial_le hmk
    have h0 : 0 ≤ Real.exp (-L) * L ^ (k - m) * (k.choose m : ℝ) := by positivity
    calc Real.exp (-L) * L ^ (k - m) * (k.choose m : ℝ) * (m.factorial : ℝ)
        ≤ Real.exp (-L) * L ^ (k - m) * (k.choose m : ℝ) * (k.factorial : ℝ) :=
          mul_le_mul_of_nonneg_left this h0
      _ = _ := by ring
  refine (Finset.sum_le_sum hb).trans (le_of_eq ?_)
  rw [← Finset.mul_sum, ← add_pow, add_comm]

end X06

namespace X06
open Real MeasureTheory Set

lemma f_nonneg (w : ℝ) (hw : 0 ≤ w) : 0 ≤ f w := by
  unfold f
  have : Real.sqrt (1 + w ^ 2) ≤ 1 + w := by
    rw [Real.sqrt_le_left (by linarith)]; nlinarith
  linarith

lemma f_le_one (w : ℝ) : f w ≤ 1 := by
  unfold f
  have : w ≤ Real.sqrt (1 + w ^ 2) := by
    rcases le_or_gt w 0 with h | h
    · exact h.trans (Real.sqrt_nonneg _)
    · rw [Real.le_sqrt h.le (by positivity)]; nlinarith
  linarith

lemma f_le_self (w : ℝ) : f w ≤ w := by
  unfold f
  have : 1 ≤ Real.sqrt (1 + w ^ 2) := by
    rw [Real.le_sqrt (by norm_num) (by positivity)]; nlinarith
  linarith

lemma f_continuous : Continuous f := by unfold f; fun_prop

/-- substitution `w = N e^{-σ}` from `(0,N)` to `(0,∞)`. -/
lemma subst_exp (N : ℝ) (hN : 0 < N) (G : ℝ → ℝ) :
    (IntegrableOn G (Ioo 0 N) ↔
      IntegrableOn (fun σ => N * Real.exp (-σ) * G (N * Real.exp (-σ))) (Ioi 0)) ∧
      ∫ w in Ioo 0 N, G w = ∫ σ in Ioi (0:ℝ), N * Real.exp (-σ) * G (N * Real.exp (-σ)) := by
  have himg : (fun σ : ℝ => N * Real.exp (-σ)) '' Ioi 0 = Ioo 0 N := by
    ext w; constructor
    · rintro ⟨σ, hσ, rfl⟩
      refine ⟨by positivity, ?_⟩
      have : Real.exp (-σ) < 1 := Real.exp_lt_one_iff.2 (by simp at hσ ⊢; linarith)
      nlinarith
    · rintro ⟨h0, h1⟩
      refine ⟨Real.log (N / w), ?_, ?_⟩
      · exact Real.log_pos ((one_lt_div h0).2 h1)
      · simp only
        rw [Real.exp_neg, Real.exp_log (by positivity)]; field_simp
  have hd : ∀ x ∈ Ioi (0:ℝ), HasDerivWithinAt (fun σ : ℝ => N * Real.exp (-σ))
      (-(N * Real.exp (-x))) (Ioi 0) x := by
    intro x _
    have := ((Real.hasDerivAt_exp (-x)).comp x (hasDerivAt_neg x)).const_mul N
    convert this.hasDerivWithinAt using 1
    · funext y; simp [Function.comp]
    · ring
  have hinj : InjOn (fun σ : ℝ => N * Real.exp (-σ)) (Ioi 0) := by
    intro x _ y _ h
    simp only at h
    have := Real.exp_injective (mul_left_cancel₀ hN.ne' h)
    linarith
  have habs : ∀ x : ℝ, |-(N * Real.exp (-x))| = N * Real.exp (-x) := by
    intro x; rw [abs_neg, abs_of_pos (by positivity)]
  constructor
  · rw [← himg, integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ioi hd hinj]
    simp only [habs, smul_eq_mul]
  · rw [← himg, integral_image_eq_integral_abs_deriv_smul measurableSet_Ioi hd hinj]
    simp only [habs, smul_eq_mul]

lemma T_subst (k : ℕ) (N : ℝ) (hN : 0 < N) :
    T k N = ∫ σ in Ioi (0:ℝ), σ ^ k * f (N * Real.exp (-σ)) := by
  unfold T
  rw [intervalIntegral.integral_of_le hN.le, integral_Ioc_eq_integral_Ioo, (subst_exp N hN _).2]
  refine setIntegral_congr_fun measurableSet_Ioi (fun σ _ => ?_)
  simp only [g]
  have hw : 0 < N * Real.exp (-σ) := by positivity
  rw [show N / (N * Real.exp (-σ)) = Real.exp σ by
    rw [Real.exp_neg]; field_simp, Real.log_exp]
  field_simp

lemma Tint_subst (k : ℕ) (N : ℝ) (hN : 0 < N) :
    IntegrableOn (fun w : ℝ => Real.log (N / w) ^ k * g w) (Ioc 0 N) ↔
      IntegrableOn (fun σ : ℝ => σ ^ k * f (N * Real.exp (-σ))) (Ioi 0) := by
  rw [integrableOn_Ioc_iff_integrableOn_Ioo, (subst_exp N hN _).1]
  refine integrableOn_congr_fun (fun σ _ => ?_) measurableSet_Ioi
  simp only [g]
  have hw : 0 < N * Real.exp (-σ) := by positivity
  rw [show N / (N * Real.exp (-σ)) = Real.exp σ by
    rw [Real.exp_neg]; field_simp, Real.log_exp]
  field_simp

lemma fsub_integrable (k : ℕ) (N : ℝ) (hN : 0 < N) :
    IntegrableOn (fun σ : ℝ => σ ^ k * f (N * Real.exp (-σ))) (Ioi 0) := by
  refine (((gammaMoment k).1).const_mul N).mono' ?_ ?_
  · exact ((continuous_pow k).mul (f_continuous.comp (by fun_prop))).aestronglyMeasurable
  · filter_upwards [ae_restrict_mem measurableSet_Ioi] with σ hσ
    have hσ0 : (0:ℝ) ≤ σ := le_of_lt hσ
    rw [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg (pow_nonneg hσ0 k) (f_nonneg _ (by positivity)))]
    have := f_le_self (N * Real.exp (-σ))
    have h0 := f_nonneg (N * Real.exp (-σ)) (by positivity)
    calc σ ^ k * f (N * Real.exp (-σ)) ≤ σ ^ k * (N * Real.exp (-σ)) :=
          mul_le_mul_of_nonneg_left this (pow_nonneg hσ0 k)
      _ = N * (σ ^ k * Real.exp (-σ)) := by ring

end X06

namespace X06
open Real MeasureTheory Set

lemma split_Ioi (h : ℝ → ℝ) (L : ℝ) (hL : 0 ≤ L) (hi : IntegrableOn h (Ioi 0)) :
    ∫ σ in Ioi (0:ℝ), h σ = (∫ σ in Ioc 0 L, h σ) + ∫ σ in Ioi L, h σ := by
  rw [← Ioc_union_Ioi_eq_Ioi hL, setIntegral_union (Ioc_disjoint_Ioi le_rfl) measurableSet_Ioi
    (hi.mono_set Ioc_subset_Ioi_self) (hi.mono_set (Ioi_subset_Ioi hL))]

lemma mmaj_nonneg (N σ : ℝ) (hN : 0 < N) : 0 ≤ mmaj N σ := by
  unfold mmaj; exact le_min (by positivity) (by positivity)

lemma mmaj_continuous (N : ℝ) : Continuous (mmaj N) := by
  unfold mmaj; fun_prop

lemma mmaj_integrable (k : ℕ) (N : ℝ) (hN : 0 < N) :
    IntegrableOn (fun σ : ℝ => σ ^ k * mmaj N σ) (Ioi 0) := by
  refine (((gammaMoment k).1).const_mul N).mono' ?_ ?_
  · exact ((continuous_pow k).mul (mmaj_continuous N)).aestronglyMeasurable
  · filter_upwards [ae_restrict_mem measurableSet_Ioi] with σ hσ
    have hσ0 : (0:ℝ) ≤ σ := le_of_lt hσ
    rw [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg (pow_nonneg hσ0 k) (mmaj_nonneg N σ hN))]
    calc σ ^ k * mmaj N σ ≤ σ ^ k * (N * Real.exp (-σ)) :=
          mul_le_mul_of_nonneg_left (min_le_left _ _) (pow_nonneg hσ0 k)
      _ = N * (σ ^ k * Real.exp (-σ)) := by ring

/-- the tail `∫_{L}^∞ σ^k (N e^{-σ})` with `N = e^L` is at most `k! (L+1)^k`. -/
lemma tail_N (k : ℕ) (N : ℝ) (hN : 1 ≤ N) (h : ℝ → ℝ) (hi : IntegrableOn h (Ioi (Real.log N)))
    (hle : ∀ σ ∈ Ioi (Real.log N), h σ ≤ σ ^ k * (N * Real.exp (-σ))) :
    ∫ σ in Ioi (Real.log N), h σ ≤ (k.factorial : ℝ) * (Real.log N + 1) ^ k := by
  have hL : 0 ≤ Real.log N := Real.log_nonneg hN
  have hN0 : 0 < N := by linarith
  obtain ⟨hti, htb⟩ := tail_moment k (Real.log N) hL
  calc ∫ σ in Ioi (Real.log N), h σ ≤ ∫ σ in Ioi (Real.log N), N * (σ ^ k * Real.exp (-σ)) := by
        refine setIntegral_mono_on hi (hti.const_mul N) measurableSet_Ioi (fun σ hσ => ?_)
        rw [show N * (σ ^ k * Real.exp (-σ)) = σ ^ k * (N * Real.exp (-σ)) by ring]
        exact hle σ hσ
    _ = N * ∫ σ in Ioi (Real.log N), σ ^ k * Real.exp (-σ) := integral_const_mul _ _
    _ ≤ N * (Real.exp (-Real.log N) * (k.factorial : ℝ) * (Real.log N + 1) ^ k) :=
        mul_le_mul_of_nonneg_left htb hN0.le
    _ = (k.factorial : ℝ) * (Real.log N + 1) ^ k := by
        rw [Real.exp_neg, Real.exp_log hN0]; field_simp

end X06

theorem X06.proof_7 : X06.stmt_7 := by
  intro k N hN
  have hN0 : 0 < N := by linarith
  set L := Real.log N with hLdef
  have hL : 0 ≤ L := Real.log_nonneg hN
  have heL : Real.exp L = N := Real.exp_log hN0
  have hLk : L ^ k ≤ (L + 1) ^ k := pow_le_pow_left₀ hL (by linarith) k
  have hfac : (0:ℝ) ≤ k.factorial := by positivity
  refine ⟨X06.mmaj_integrable k N hN0, ?_, ?_⟩
  · -- I_k
    unfold X06.Ik
    rw [X06.split_Ioi _ L hL (X06.mmaj_integrable k N hN0)]
    have h1 : ∫ σ in Ioc 0 L, σ ^ k * X06.mmaj N σ ≤ L ^ k / 2 := by
      calc ∫ σ in Ioc 0 L, σ ^ k * X06.mmaj N σ
          ≤ ∫ σ in Ioc 0 L, L ^ k / (2 * N) * Real.exp σ := by
            refine setIntegral_mono_on
              (((continuous_pow k).mul (X06.mmaj_continuous N)).integrableOn_Ioc)
              ((continuous_const.mul Real.continuous_exp).integrableOn_Ioc)
              measurableSet_Ioc (fun σ hσ => ?_)
            have hσk : σ ^ k ≤ L ^ k := pow_le_pow_left₀ hσ.1.le hσ.2 k
            have hm : X06.mmaj N σ ≤ Real.exp σ / (2 * N) := min_le_right _ _
            calc σ ^ k * X06.mmaj N σ ≤ L ^ k * (Real.exp σ / (2 * N)) :=
                  mul_le_mul hσk hm (X06.mmaj_nonneg N σ hN0) (pow_nonneg hL k)
              _ = L ^ k / (2 * N) * Real.exp σ := by ring
        _ = L ^ k / (2 * N) * (Real.exp L - 1) := by
            rw [integral_const_mul, ← intervalIntegral.integral_of_le hL, integral_exp,
              Real.exp_zero]
        _ ≤ L ^ k / 2 := by
            rw [heL]
            have : L ^ k / (2 * N) * (N - 1) = L ^ k / 2 - L ^ k / (2 * N) := by
              field_simp
            rw [this]
            have : 0 ≤ L ^ k / (2 * N) := by positivity
            linarith
    have h2 := X06.tail_N k N hN (fun σ => σ ^ k * X06.mmaj N σ)
      ((X06.mmaj_integrable k N hN0).mono_set (Ioi_subset_Ioi hL))
      (fun σ hσ => mul_le_mul_of_nonneg_left (min_le_left _ _)
        (pow_nonneg (hL.trans hσ.le) k))
    nlinarith
  · -- T_k
    rw [X06.T_subst k N hN0, X06.split_Ioi _ L hL (X06.fsub_integrable k N hN0)]
    have h1 : ∫ σ in Ioc 0 L, σ ^ k * X06.f (N * Real.exp (-σ)) ≤ L ^ (k + 1) / (k + 1) := by
      calc ∫ σ in Ioc 0 L, σ ^ k * X06.f (N * Real.exp (-σ))
          ≤ ∫ σ in Ioc 0 L, σ ^ k := by
            refine setIntegral_mono_on
              (((continuous_pow k).mul (X06.f_continuous.comp (by fun_prop))).integrableOn_Ioc)
              ((continuous_pow k).integrableOn_Ioc)
              measurableSet_Ioc (fun σ hσ => ?_)
            have := X06.f_le_one (N * Real.exp (-σ))
            have h0 : 0 ≤ σ ^ k := pow_nonneg hσ.1.le k
            nlinarith
        _ = L ^ (k + 1) / (k + 1) := by
            rw [← intervalIntegral.integral_of_le hL, integral_pow]; simp
    have h2 := X06.tail_N k N hN (fun σ => σ ^ k * X06.f (N * Real.exp (-σ)))
      ((X06.fsub_integrable k N hN0).mono_set (Ioi_subset_Ioi hL))
      (fun σ hσ => mul_le_mul_of_nonneg_left (X06.f_le_self _)
        (pow_nonneg (hL.trans hσ.le) k))
    have h3 : L ^ (k + 1) / (k + 1) ≤ (L + 1) ^ (k + 1) := by
      have : L ^ (k + 1) ≤ (L + 1) ^ (k + 1) := pow_le_pow_left₀ hL (by linarith) (k + 1)
      have hk1 : (1:ℝ) ≤ k + 1 := by have := (Nat.cast_nonneg k : (0:ℝ) ≤ k); linarith
      calc L ^ (k + 1) / (k + 1) ≤ L ^ (k + 1) := div_le_self (by positivity) hk1
        _ ≤ _ := this
    have h4 : (k.factorial : ℝ) * (L + 1) ^ k ≤ (k.factorial : ℝ) * (L + 1) ^ (k + 1) := by
      apply mul_le_mul_of_nonneg_left _ hfac
      rw [pow_succ]
      have : 0 ≤ (L + 1) ^ k := by positivity
      nlinarith
    nlinarith

namespace X06
open Real

lemma NN_ge_one (n : ℝ) (hn : 1 < n) : 1 ≤ NN n := by
  unfold NN
  have : 1 ≤ Real.sqrt 2 := by rw [Real.le_sqrt (by norm_num) (by norm_num)]; norm_num
  nlinarith

lemma NN_pos (n : ℝ) (hn : 0 < n) : 0 < NN n := by unfold NN; positivity

lemma exp_neg_quad (x : ℝ) (hx : 0 ≤ x) :
    0 ≤ Real.exp (-x) - 1 + x ∧ Real.exp (-x) - 1 + x ≤ x ^ 2 / 2 := by
  constructor
  · have := Real.add_one_le_exp (-x); linarith
  · have h1 := Real.quadratic_le_exp_of_nonneg hx
    have hpos : 0 < 1 + x + x ^ 2 / 2 := by positivity
    have h2 : Real.exp (-x) ≤ 1 / (1 + x + x ^ 2 / 2) := by
      rw [Real.exp_neg, one_div]; exact inv_anti₀ hpos h1
    have h3 : 1 / (1 + x + x ^ 2 / 2) ≤ 1 - x + x ^ 2 / 2 := by
      rw [div_le_iff₀ hpos]; nlinarith [sq_nonneg (x ^ 2)]
    linarith

lemma W0_eq (n σ : ℝ) : W0 n σ = NN n * Real.exp (-σ) := rfl

lemma Delta2_integrable (n : ℝ) (hn : 1 < n) :
    IntegrableOn (fun σ : ℝ => (Real.exp (-σ / n) - 1) * σ * f (W0 n σ)) (Ioi 0) := by
  have hN := NN_pos n (by linarith)
  refine (fsub_integrable 1 (NN n) hN).mono' ?_ ?_
  · exact (((by fun_prop : Continuous fun σ : ℝ => (Real.exp (-σ / n) - 1) * σ)).mul
      (f_continuous.comp (by unfold W0; fun_prop))).aestronglyMeasurable
  · filter_upwards [ae_restrict_mem measurableSet_Ioi] with σ hσ
    have hσ0 : (0:ℝ) ≤ σ := le_of_lt hσ
    have hf0 := f_nonneg (W0 n σ) (by unfold W0; positivity)
    have he : Real.exp (-σ / n) ≤ 1 := Real.exp_le_one_iff.2 (by
      have : 0 ≤ σ / n := div_nonneg hσ0 (by linarith)
      rw [neg_div]; linarith)
    have he0 : 0 < Real.exp (-σ / n) := Real.exp_pos _
    rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_of_nonneg hσ0, abs_of_nonneg hf0, pow_one,
      ← W0_eq]
    have : |Real.exp (-σ / n) - 1| ≤ 1 := abs_le.2 ⟨by linarith, by linarith⟩
    have h0 : 0 ≤ σ * f (W0 n σ) := mul_nonneg hσ0 hf0
    nlinarith

end X06

theorem X06.proof_10 : X06.stmt_10 := by
  intro n hn
  have hn0 : 0 < n := by linarith
  have hN := X06.NN_pos n hn0
  have hN1 := X06.NN_ge_one n hn
  have key : X06.Delta2 n + (1 / n) * X06.T 2 (X06.NN n)
      = ∫ σ in Ioi (0:ℝ), (Real.exp (-(σ / n)) - 1 + σ / n) * (σ * X06.f (X06.W0 n σ)) := by
    rw [X06.T_subst 2 _ hN]
    unfold X06.Delta2
    rw [← integral_const_mul, ← integral_add (X06.Delta2_integrable n hn)
      ((X06.fsub_integrable 2 _ hN).const_mul _)]
    refine setIntegral_congr_fun measurableSet_Ioi (fun σ _ => ?_)
    simp only [X06.W0_eq, neg_div]
    ring
  have hT3 : X06.T 3 (X06.NN n)
      = ∫ σ in Ioi (0:ℝ), σ ^ 3 * X06.f (X06.W0 n σ) := X06.T_subst 3 _ hN
  have hint : IntegrableOn
      (fun σ : ℝ => (Real.exp (-(σ / n)) - 1 + σ / n) * (σ * X06.f (X06.W0 n σ))) (Ioi 0) := by
    have := (X06.Delta2_integrable n hn).add ((X06.fsub_integrable 2 _ hN).const_mul (1 / n))
    refine this.congr_fun (fun σ _ => ?_) measurableSet_Ioi
    simp only [Pi.add_apply, X06.W0_eq, neg_div]; ring
  refine ⟨?_, ?_, ?_⟩
  · rw [key]
    refine setIntegral_nonneg measurableSet_Ioi (fun σ hσ => ?_)
    have hσ0 : (0:ℝ) ≤ σ := le_of_lt hσ
    have := (X06.exp_neg_quad (σ / n) (div_nonneg hσ0 hn0.le)).1
    exact mul_nonneg this (mul_nonneg hσ0 (X06.f_nonneg _ (by unfold X06.W0; positivity)))
  · rw [key, hT3, ← integral_const_mul]
    refine setIntegral_mono_on hint ((X06.fsub_integrable 3 _ hN).const_mul _)
      measurableSet_Ioi (fun σ hσ => ?_)
    have hσ0 : (0:ℝ) ≤ σ := le_of_lt hσ
    have hq := (X06.exp_neg_quad (σ / n) (div_nonneg hσ0 hn0.le)).2
    have h0 : 0 ≤ σ * X06.f (X06.W0 n σ) :=
      mul_nonneg hσ0 (X06.f_nonneg _ (by unfold X06.W0; positivity))
    calc (Real.exp (-(σ / n)) - 1 + σ / n) * (σ * X06.f (X06.W0 n σ))
        ≤ (σ / n) ^ 2 / 2 * (σ * X06.f (X06.W0 n σ)) := mul_le_mul_of_nonneg_right hq h0
      _ = (1 / n) ^ 2 / 2 * (σ ^ 3 * X06.f (X06.W0 n σ)) := by ring
  · have h7 := (X06.proof_7 3 (X06.NN n) hN1).2.2
    have : ((Nat.factorial 3 : ℕ) : ℝ) = 6 := by norm_num [Nat.factorial]
    rw [this] at h7
    have hp : 0 ≤ (1 / n) ^ 2 / 2 := by positivity
    calc (1 / n) ^ 2 / 2 * X06.T 3 (X06.NN n)
        ≤ (1 / n) ^ 2 / 2 * ((6 + 1) * (Real.log (X06.NN n) + 1) ^ (3 + 1)) :=
          mul_le_mul_of_nonneg_left h7 hp
      _ = 3.5 * (1 / n) ^ 2 * (Real.log (X06.NN n) + 1) ^ 4 := by norm_num; ring

namespace X06
open Real

lemma W_eq (n σ δ : ℝ) : W n σ δ = W0 n σ * Real.exp (kappa n σ δ) := by
  unfold W W0 kappa
  have hc : 0 < Real.cosh (δ / n) := Real.cosh_pos _
  rw [Real.exp_add, Real.sqrt_eq_rpow, Real.rpow_def_of_pos hc]
  ring_nf

lemma W_continuous (n : ℝ) : Continuous (fun p : ℝ × ℝ => W n p.1 p.2) := by
  unfold W NN; fun_prop

lemma log_cosh_le (x : ℝ) (hx : 0 ≤ x) : Real.log (Real.cosh x) ≤ x := by
  have hc : Real.cosh x ≤ Real.exp x := by
    rw [Real.cosh_eq]
    have : Real.exp (-x) ≤ Real.exp x := Real.exp_le_exp.2 (by linarith)
    linarith
  calc Real.log (Real.cosh x) ≤ Real.log (Real.exp x) :=
        Real.log_le_log (Real.cosh_pos x) hc
    _ = x := Real.log_exp x

lemma log_cosh_nonneg (x : ℝ) : 0 ≤ Real.log (Real.cosh x) :=
  Real.log_nonneg (Real.one_le_cosh x)

lemma kappa_le (n σ δ : ℝ) (hn : 0 < n) (hδ : 0 ≤ δ) (hδσ : δ ≤ σ) :
    kappa n σ δ ≤ σ / n := by
  unfold kappa
  have := log_cosh_le (δ / n) (div_nonneg hδ hn.le)
  have : δ / n ≤ σ / n := div_le_div_of_nonneg_right hδσ hn.le
  linarith

lemma kappa_nonneg (n σ δ : ℝ) (hn : 0 < n) (hσ : 0 ≤ σ) : 0 ≤ kappa n σ δ := by
  unfold kappa
  have := log_cosh_nonneg (δ / n)
  have : 0 ≤ σ / n := div_nonneg hσ hn.le
  linarith

lemma D1inner_eq (n σ : ℝ) :
    D1inner n σ = Real.exp (-σ / n) * ∫ δ in (0:ℝ)..σ, (f (W n σ δ) - f (W0 n σ)) := by
  unfold D1inner
  simp only [W_eq]

lemma D1inner_continuous (n : ℝ) : Continuous (D1inner n) := by
  have : D1inner n = fun σ =>
      Real.exp (-σ / n) * ∫ δ in (0:ℝ)..σ, (f (W n σ δ) - f (W0 n σ)) := funext (D1inner_eq n)
  rw [this]
  refine (by fun_prop : Continuous fun σ : ℝ => Real.exp (-σ / n)).mul ?_
  refine intervalIntegral.continuous_parametric_intervalIntegral_of_continuous
    (f := fun σ δ => f (W n σ δ) - f (W0 n σ)) ?_ continuous_id
  have h1 : Continuous (fun p : ℝ × ℝ => f (W n p.1 p.2)) := f_continuous.comp (W_continuous n)
  have h2 : Continuous (fun p : ℝ × ℝ => f (W0 n p.1)) :=
    f_continuous.comp (by unfold W0 NN; fun_prop)
  exact h1.sub h2

lemma D1inner_bound (n σ : ℝ) (hn : 0 < n) (hσ : 0 ≤ σ) :
    |D1inner n σ| ≤ NN n * (σ * Real.exp (-σ)) := by
  unfold D1inner
  have hN := NN_pos n hn
  set M := W0 n σ * Real.exp (σ / n)
  have hM0 : 0 ≤ M := by unfold M W0; positivity
  have hb : ∀ δ ∈ Set.uIoc (0:ℝ) σ,
      ‖f (W0 n σ * Real.exp (kappa n σ δ)) - f (W0 n σ)‖ ≤ M := by
    intro δ hδ
    rw [Set.uIoc_of_le hσ] at hδ
    have hW0 : 0 ≤ W0 n σ := by unfold W0; positivity
    have ha0 := f_nonneg (W0 n σ * Real.exp (kappa n σ δ)) (by positivity)
    have hb0 := f_nonneg (W0 n σ) hW0
    have ha1 : f (W0 n σ * Real.exp (kappa n σ δ)) ≤ M :=
      (f_le_self _).trans (mul_le_mul_of_nonneg_left
        (Real.exp_le_exp.2 (kappa_le n σ δ hn hδ.1.le hδ.2)) hW0)
    have hb1 : f (W0 n σ) ≤ M := (f_le_self _).trans (le_mul_of_one_le_right hW0
      (Real.one_le_exp (div_nonneg hσ hn.le)))
    rw [Real.norm_eq_abs, abs_le]; constructor <;> linarith
  have hI := intervalIntegral.norm_integral_le_of_norm_le_const hb
  rw [Real.norm_eq_abs, sub_zero, abs_of_nonneg hσ] at hI
  rw [abs_mul, abs_of_pos (Real.exp_pos _)]
  calc Real.exp (-σ / n) * |∫ δ in (0:ℝ)..σ, (f (W0 n σ * Real.exp (kappa n σ δ)) - f (W0 n σ))|
      ≤ Real.exp (-σ / n) * (M * σ) := mul_le_mul_of_nonneg_left hI (Real.exp_pos _).le
    _ = NN n * (σ * Real.exp (-σ)) := by
      unfold M W0
      rw [show Real.exp (-σ / n) * (NN n * Real.exp (-σ) * Real.exp (σ / n) * σ)
          = NN n * (σ * Real.exp (-σ)) * (Real.exp (-σ / n) * Real.exp (σ / n)) by ring,
        ← Real.exp_add, show -σ / n + σ / n = 0 by ring, Real.exp_zero, mul_one]

lemma D1inner_integrable (n : ℝ) (hn : 0 < n) : IntegrableOn (D1inner n) (Ioi 0) := by
  refine (((gammaMoment 1).1).const_mul (NN n)).mono' (D1inner_continuous n).aestronglyMeasurable ?_
  filter_upwards [ae_restrict_mem measurableSet_Ioi] with σ hσ
  rw [Real.norm_eq_abs, pow_one]
  exact D1inner_bound n σ hn (le_of_lt hσ)

lemma Cinner_eq (n σ : ℝ) :
    Cinner n σ = D1inner n σ + (Real.exp (-σ / n) - 1) * σ * f (W0 n σ) + σ * f (W0 n σ) := by
  rw [D1inner_eq]
  unfold Cinner
  have hW : Continuous (fun δ => W n σ δ) := by unfold W; fun_prop
  have hc : Continuous (fun δ => f (W n σ δ)) := f_continuous.comp hW
  rw [intervalIntegral.integral_sub (hc.intervalIntegrable 0 σ) intervalIntegrable_const]
  simp only [intervalIntegral.integral_const, sub_zero, smul_eq_mul]
  ring

end X06

theorem X06.proof_8 : X06.stmt_8 := by
  intro n hn
  have hn0 : 0 < n := by linarith
  have hN := X06.NN_pos n hn0
  have hA : IntegrableOn (fun σ : ℝ => σ * X06.f (X06.W0 n σ)) (Ioi 0) := by
    simpa [X06.W0_eq] using X06.fsub_integrable 1 (X06.NN n) hN
  have hB := X06.D1inner_integrable n hn0
  have hC := X06.Delta2_integrable n hn
  refine ⟨hA, hB, hC, ?_, ?_⟩
  · unfold X06.Cint X06.C0 X06.Delta1 X06.Delta2
    simp_rw [X06.Cinner_eq]
    rw [integral_add (f := fun σ => X06.D1inner n σ + (Real.exp (-σ / n) - 1) * σ * X06.f (X06.W0 n σ))
      (hB.add hC) hA, integral_add hB hC]
    ring
  · unfold X06.C0
    rw [X06.T_subst 1 _ hN]
    simp [X06.W0_eq]

namespace X06
open Real Filter Topology

lemma psi_continuous : Continuous psi := by
  have : psi = fun w => w - w ^ 2 / Real.sqrt (1 + w ^ 2) := funext psi_eq
  rw [this]
  exact continuous_id.sub ((continuous_pow 2).div (by fun_prop)
    (fun w => (Real.sqrt_pos.2 (by positivity)).ne'))

lemma psi_nonneg_le (w : ℝ) (hw : 0 < w) : 0 ≤ psi w ∧ psi w ≤ min w (1 / (2 * w)) :=
  ⟨(X06.proof_6 w hw).1, (X06.proof_6 w hw).2.1⟩

/-- `d/dτ f(a e^τ) = ψ(a e^τ)`. -/
lemma hasDerivAt_f_exp (a τ : ℝ) :
    HasDerivAt (fun t => f (a * Real.exp t)) (psi (a * Real.exp τ)) τ := by
  have h1 := (Real.hasDerivAt_exp τ).const_mul a
  have := (hasDerivAt_f (a * Real.exp τ)).comp τ h1
  refine this.congr_deriv ?_
  rw [psi_eq]; ring

/-- `d/dτ ψ(a e^τ) = χ(a e^τ)`. -/
lemma hasDerivAt_psi_exp (a τ : ℝ) :
    HasDerivAt (fun t => psi (a * Real.exp t)) (chi (a * Real.exp τ)) τ := by
  have h1 := (Real.hasDerivAt_exp τ).const_mul a
  have := (hasDerivAt_psi (a * Real.exp τ)).comp τ h1
  refine this.congr_deriv ?_
  rw [chi_eq]; ring

/-- `d/dσ f(N e^{-σ}) = -ψ(N e^{-σ})`. -/
lemma hasDerivAt_f_Nexp (N σ : ℝ) :
    HasDerivAt (fun s => f (N * Real.exp (-s))) (-psi (N * Real.exp (-σ))) σ := by
  have h1 : HasDerivAt (fun s => -s) (-1) σ := hasDerivAt_neg σ
  have := (hasDerivAt_f_exp N (-σ)).comp σ h1
  convert this using 1
  · funext s; simp [Function.comp]
  · ring

lemma psi_sub_integrable (k : ℕ) (N : ℝ) (hN : 0 < N) :
    IntegrableOn (fun σ : ℝ => σ ^ k * psi (N * Real.exp (-σ))) (Ioi 0) := by
  refine (((gammaMoment k).1).const_mul N).mono' ?_ ?_
  · exact ((continuous_pow k).mul (psi_continuous.comp (by fun_prop))).aestronglyMeasurable
  · filter_upwards [ae_restrict_mem measurableSet_Ioi] with σ hσ
    have hσ0 : (0:ℝ) ≤ σ := le_of_lt hσ
    have hw : 0 < N * Real.exp (-σ) := by positivity
    obtain ⟨h0, h1⟩ := psi_nonneg_le _ hw
    rw [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg (pow_nonneg hσ0 k) h0)]
    calc σ ^ k * psi (N * Real.exp (-σ)) ≤ σ ^ k * (N * Real.exp (-σ)) :=
          mul_le_mul_of_nonneg_left (h1.trans (min_le_left _ _)) (pow_nonneg hσ0 k)
      _ = N * (σ ^ k * Real.exp (-σ)) := by ring

/-- integration by parts: `∫_0^∞ σ² ψ(N e^{-σ}) dσ = 2 T_1(g)(N)`. -/
lemma ibp_psi (N : ℝ) (hN : 0 < N) :
    ∫ σ in Ioi (0:ℝ), σ ^ 2 * psi (N * Real.exp (-σ)) = 2 * T 1 N := by
  have hu : ∀ x ∈ Ioi (0:ℝ), HasDerivAt (fun s : ℝ => s ^ 2) (2 * x) x := by
    intro x _; simpa using hasDerivAt_pow 2 x
  have hv : ∀ x ∈ Ioi (0:ℝ), HasDerivAt (fun s => f (N * Real.exp (-s)))
      (-psi (N * Real.exp (-x))) x := fun x _ => hasDerivAt_f_Nexp N x
  have huv' : IntegrableOn ((fun s : ℝ => s ^ 2) * fun x => -psi (N * Real.exp (-x))) (Ioi 0) := by
    have := (psi_sub_integrable 2 N hN).neg
    refine this.congr_fun (fun x _ => ?_) measurableSet_Ioi
    simp [Pi.mul_apply]
  have hu'v : IntegrableOn ((fun x : ℝ => 2 * x) * fun s => f (N * Real.exp (-s))) (Ioi 0) := by
    have := (fsub_integrable 1 N hN).const_mul 2
    refine IntegrableOn.congr_fun this (fun x _ => ?_) measurableSet_Ioi
    simp [Pi.mul_apply]; ring
  have hcont : Continuous ((fun s : ℝ => s ^ 2) * fun s => f (N * Real.exp (-s))) :=
    (continuous_pow 2).mul (f_continuous.comp (by fun_prop))
  have h0 : Tendsto ((fun s : ℝ => s ^ 2) * fun s => f (N * Real.exp (-s))) (𝓝[>] 0) (𝓝 0) := by
    have := (hcont.tendsto 0).mono_left (nhdsWithin_le_nhds (s := Ioi 0))
    simpa using this
  have hinf : Tendsto ((fun s : ℝ => s ^ 2) * fun s => f (N * Real.exp (-s))) atTop (𝓝 0) := by
    have hlim : Tendsto (fun s : ℝ => N * (s ^ 2 * Real.exp (-s))) atTop (𝓝 0) := by
      simpa using (tendsto_pow_mul_exp_neg_atTop_nhds_zero 2).const_mul N
    refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hlim ?_ ?_
    · filter_upwards [eventually_ge_atTop 0] with s hs
      simp only [Pi.mul_apply]
      exact mul_nonneg (by positivity) (f_nonneg _ (by positivity))
    · filter_upwards [eventually_ge_atTop 0] with s hs
      simp only [Pi.mul_apply]
      have := f_le_self (N * Real.exp (-s))
      have h2 : 0 ≤ s ^ 2 := by positivity
      nlinarith
  have := integral_Ioi_mul_deriv_eq_deriv_mul hu hv huv' hu'v h0 hinf
  rw [T_subst 1 N hN]
  have e1 : ∫ x in Ioi (0:ℝ), x ^ 2 * -psi (N * Real.exp (-x))
      = -∫ x in Ioi (0:ℝ), x ^ 2 * psi (N * Real.exp (-x)) := by
    rw [← integral_neg]; congr 1; funext x; ring
  have e2 : ∫ x in Ioi (0:ℝ), 2 * x * f (N * Real.exp (-x))
      = 2 * ∫ x in Ioi (0:ℝ), x ^ 1 * f (N * Real.exp (-x)) := by
    rw [← integral_const_mul]; congr 1; funext x; ring
  rw [e1, e2] at this
  linarith

end X06

namespace X06
open Real

lemma mmaj_eq (n σ : ℝ) : mmaj (NN n) σ = min (W0 n σ) (1 / (2 * W0 n σ)) := by
  unfold mmaj W0
  congr 1
  rw [Real.exp_neg]; field_simp

lemma chi_exp_bound (a s b : ℝ) (ha : 0 < a) (hs : 0 ≤ s) (hsb : s ≤ b) :
    |chi (a * Real.exp s)| ≤ Real.exp b * min a (1 / (2 * a)) := by
  have hx : 0 < a * Real.exp s := by positivity
  obtain ⟨_, h2, h3⟩ := X06.proof_6 _ hx
  have hes : Real.exp (-s) ≤ Real.exp s := Real.exp_le_exp.2 (by linarith)
  have h4 : min (a * Real.exp s) (1 / (2 * (a * Real.exp s))) ≤ Real.exp s * min a (1 / (2 * a)) := by
    rw [mul_min_of_nonneg _ _ (Real.exp_pos s).le]
    refine min_le_min (by rw [mul_comm]) ?_
    rw [show 1 / (2 * (a * Real.exp s)) = Real.exp (-s) * (1 / (2 * a)) by
      rw [Real.exp_neg]; field_simp]
    exact mul_le_mul_of_nonneg_right hes (by positivity)
  have h5 : Real.exp s * min a (1 / (2 * a)) ≤ Real.exp b * min a (1 / (2 * a)) :=
    mul_le_mul_of_nonneg_right (Real.exp_le_exp.2 hsb) (le_min ha.le (by positivity))
  linarith

lemma log_cosh_le_sq (x : ℝ) : Real.log (Real.cosh x) ≤ x ^ 2 / 2 := by
  calc Real.log (Real.cosh x) ≤ Real.log (Real.exp (x ^ 2 / 2)) :=
        Real.log_le_log (Real.cosh_pos x) (Real.cosh_le_exp_half_sq x)
    _ = x ^ 2 / 2 := Real.log_exp _

lemma kappa_continuous (n σ : ℝ) : Continuous (fun δ => kappa n σ δ) := by
  unfold kappa
  refine continuous_const.add (Continuous.div_const ?_ _)
  exact (Real.continuous_cosh.comp (by fun_prop)).log (fun x => (Real.cosh_pos _).ne')

/-- second-order Taylor bound for `τ ↦ f(a e^τ)` on `[0,b]`. -/
lemma taylor2 (a b : ℝ) (ha : 0 < a) (_hb : 0 ≤ b) :
    ∀ τ ∈ Icc 0 b, |f (a * Real.exp τ) - f a - psi a * τ|
      ≤ Real.exp b * min a (1 / (2 * a)) * τ ^ 2 / 2 := by
  set C := Real.exp b * min a (1 / (2 * a))
  -- step A
  have hA : ∀ τ ∈ Icc 0 b, ‖psi (a * Real.exp τ) - psi (a * Real.exp 0)‖ ≤ C * (τ - 0) := by
    refine norm_image_sub_le_of_norm_deriv_right_le_segment
      (fun x _ => (hasDerivAt_psi_exp a x).continuousAt.continuousWithinAt)
      (fun x _ => (hasDerivAt_psi_exp a x).hasDerivWithinAt) (fun x hx => ?_)
    exact chi_exp_bound a x b ha hx.1 hx.2.le
  simp only [Real.exp_zero, mul_one, sub_zero] at hA
  have hφ : ∀ τ, HasDerivAt (fun t => f (a * Real.exp t) - f a - psi a * t)
      (psi (a * Real.exp τ) - psi a) τ := by
    intro τ
    have := ((hasDerivAt_f_exp a τ).sub_const (f a)).sub ((hasDerivAt_id τ).const_mul (psi a))
    convert this using 1
    · funext t; simp
    · ring
  have hB : ∀ τ, HasDerivAt (fun t : ℝ => C * t ^ 2 / 2) (C * τ) τ := by
    intro τ
    have := ((hasDerivAt_pow 2 τ).const_mul C).div_const 2
    convert this using 1; push_cast; ring
  have key := image_norm_le_of_norm_deriv_right_le_deriv_boundary
    (f := fun t => f (a * Real.exp t) - f a - psi a * t)
    (fun x _ => (hφ x).continuousAt.continuousWithinAt)
    (fun x _ => (hφ x).hasDerivWithinAt) (B := fun t => C * t ^ 2 / 2)
    (by simp) hB (fun x hx => hA x (Ico_subset_Icc_self hx))
  intro τ hτ
  have := key hτ
  simpa [Real.norm_eq_abs] using this

lemma D1_pointwise (n σ : ℝ) (hn : 1 < n) (hσ : 0 ≤ σ) :
    |D1inner n σ - (1 / n) / 2 * (σ ^ 2 * psi (W0 n σ))|
      ≤ 13 / 12 * (1 / n) ^ 2 * (σ ^ 3 * mmaj (NN n) σ) := by
  have hn0 : 0 < n := by linarith
  set a := W0 n σ with ha_def
  have ha : 0 < a := by rw [ha_def]; unfold W0; exact mul_pos (NN_pos n hn0) (Real.exp_pos _)
  set b := σ / n with hb_def
  have hb : 0 ≤ b := div_nonneg hσ hn0.le
  set m := mmaj (NN n) σ with hm_def
  have hm : m = min a (1 / (2 * a)) := mmaj_eq n σ
  have hm0 : 0 ≤ m := by rw [hm]; exact le_min ha.le (by positivity)
  set C := Real.exp b * m with hC
  set P := psi a with hP
  obtain ⟨hP0, hP1⟩ := psi_nonneg_le a ha
  have hPm : P ≤ m := by rw [hm]; exact hP1
  -- remainder
  have hT := taylor2 a b ha hb
  rw [← hm] at hT
  have hκc := kappa_continuous n σ
  have hfc : Continuous (fun δ => f (a * Real.exp (kappa n σ δ))) :=
    f_continuous.comp (continuous_const.mul (Real.continuous_exp.comp hκc))
  have hR : ∀ δ ∈ Set.uIoc (0:ℝ) σ,
      ‖f (a * Real.exp (kappa n σ δ)) - f a - P * kappa n σ δ‖ ≤ C * b ^ 2 / 2 := by
    intro δ hδ
    rw [Set.uIoc_of_le hσ] at hδ
    have hk0 := kappa_nonneg n σ δ hn0 hσ
    have hk1 := kappa_le n σ δ hn0 hδ.1.le hδ.2
    have := hT (kappa n σ δ) ⟨hk0, hk1⟩
    rw [Real.norm_eq_abs]
    refine this.trans ?_
    have hC0 : 0 ≤ C := by positivity
    have : kappa n σ δ ^ 2 ≤ b ^ 2 := pow_le_pow_left₀ hk0 hk1 2
    have := mul_le_mul_of_nonneg_left this hC0
    linarith
  set Rint := ∫ δ in (0:ℝ)..σ, (f (a * Real.exp (kappa n σ δ)) - f a - P * kappa n σ δ)
  have hRb : |Rint| ≤ C * b ^ 2 / 2 * σ := by
    have := intervalIntegral.norm_integral_le_of_norm_le_const hR
    rwa [Real.norm_eq_abs, sub_zero, abs_of_nonneg hσ] at this
  set Kint := ∫ δ in (0:ℝ)..σ, kappa n σ δ
  have hK1 : σ ^ 2 / (2 * n) ≤ Kint := by
    have : ∫ δ in (0:ℝ)..σ, (σ / n / 2) ≤ Kint :=
      intervalIntegral.integral_mono_on hσ intervalIntegrable_const (hκc.intervalIntegrable 0 σ)
        (fun δ _ => by
          unfold kappa
          have := log_cosh_nonneg (δ / n)
          linarith)
    rw [intervalIntegral.integral_const, smul_eq_mul, sub_zero] at this
    have e : σ * (σ / n / 2) = σ ^ 2 / (2 * n) := by field_simp
    linarith
  have hK2 : Kint ≤ σ ^ 2 / (2 * n) + σ ^ 3 / (12 * n ^ 2) := by
    have hc2 : Continuous (fun δ : ℝ => σ / n / 2 + 1 / (4 * n ^ 2) * δ ^ 2) := by fun_prop
    have : Kint ≤ ∫ δ in (0:ℝ)..σ, (σ / n / 2 + 1 / (4 * n ^ 2) * δ ^ 2) :=
      intervalIntegral.integral_mono_on hσ (hκc.intervalIntegrable 0 σ) (hc2.intervalIntegrable 0 σ)
        (fun δ _ => by
          unfold kappa
          have := log_cosh_le_sq (δ / n)
          have e : (δ / n) ^ 2 / 2 / 2 = 1 / (4 * n ^ 2) * δ ^ 2 := by field_simp; ring
          linarith)
    rw [intervalIntegral.integral_add (f := fun _ => σ / n / 2)
      (g := fun δ => 1 / (4 * n ^ 2) * δ ^ 2) intervalIntegrable_const
      ((continuous_const.mul (continuous_pow 2)).intervalIntegrable 0 σ),
      intervalIntegral.integral_const, intervalIntegral.integral_const_mul, integral_pow] at this
    have e : (σ - 0) • (σ / n / 2) + 1 / (4 * n ^ 2) * ((σ ^ (2 + 1) - 0 ^ (2 + 1)) / (2 + 1))
        = σ ^ 2 / (2 * n) + σ ^ 3 / (12 * n ^ 2) := by
      simp only [smul_eq_mul]; field_simp; ring
    linarith
  -- decomposition
  have hdec : ∫ δ in (0:ℝ)..σ, (f (a * Real.exp (kappa n σ δ)) - f a) = Rint + P * Kint := by
    have hadd := intervalIntegral.integral_add (μ := MeasureTheory.volume)
      (f := fun δ => f (a * Real.exp (kappa n σ δ)) - f a - P * kappa n σ δ)
      (g := fun δ => P * kappa n σ δ)
      (((hfc.sub continuous_const).sub (continuous_const.mul hκc)).intervalIntegrable 0 σ)
      ((continuous_const.mul hκc).intervalIntegrable 0 σ)
    calc ∫ δ in (0:ℝ)..σ, (f (a * Real.exp (kappa n σ δ)) - f a)
        = ∫ δ in (0:ℝ)..σ, ((f (a * Real.exp (kappa n σ δ)) - f a - P * kappa n σ δ)
            + P * kappa n σ δ) := by congr 1; funext δ; ring
      _ = Rint + ∫ δ in (0:ℝ)..σ, P * kappa n σ δ := hadd
      _ = Rint + P * Kint := by rw [intervalIntegral.integral_const_mul]
  have hD : D1inner n σ = Real.exp (-σ / n) * (Rint + P * Kint) := by
    unfold D1inner; rw [← hdec]
  rw [hD]
  set E := Real.exp (-σ / n) with hE
  have hE0 : 0 < E := Real.exp_pos _
  have hE1 : E ≤ 1 := Real.exp_le_one_iff.2 (by rw [neg_div]; linarith)
  have hE2 : 1 - b ≤ E := by
    have h : -σ / n + 1 ≤ E := Real.add_one_le_exp (-σ / n)
    have e : -σ / n = -b := by rw [hb_def, neg_div]
    linarith
  have hEb : E * Real.exp b = 1 := by
    rw [hE, ← Real.exp_add, hb_def, neg_div, neg_add_cancel, Real.exp_zero]
  -- the two pieces
  have hX : |E * Rint| ≤ m * σ ^ 3 / (2 * n ^ 2) := by
    rw [abs_mul, abs_of_pos hE0]
    calc E * |Rint| ≤ E * (C * b ^ 2 / 2 * σ) := mul_le_mul_of_nonneg_left hRb hE0.le
      _ = (E * Real.exp b) * m * b ^ 2 / 2 * σ := by rw [hC]; ring
      _ = m * σ ^ 3 / (2 * n ^ 2) := by rw [hEb, hb_def]; field_simp
  have hY : |E * Kint - σ ^ 2 / (2 * n)| ≤ σ ^ 3 / (2 * n ^ 2) := by
    have hq : 0 ≤ σ ^ 2 / (2 * n) := by positivity
    have hK0 : 0 ≤ Kint := hq.trans hK1
    rw [abs_le]; constructor
    · have h1 : E * (σ ^ 2 / (2 * n)) ≤ E * Kint := mul_le_mul_of_nonneg_left hK1 hE0.le
      have h2 : (1 - b) * (σ ^ 2 / (2 * n)) ≤ E * (σ ^ 2 / (2 * n)) :=
        mul_le_mul_of_nonneg_right hE2 hq
      have e : (1 - b) * (σ ^ 2 / (2 * n)) = σ ^ 2 / (2 * n) - σ ^ 3 / (2 * n ^ 2) := by
        rw [hb_def]; field_simp
      linarith
    · have h1 : E * Kint ≤ Kint := by nlinarith
      have h3 : σ ^ 3 / (12 * n ^ 2) ≤ σ ^ 3 / (2 * n ^ 2) :=
        div_le_div_of_nonneg_left (by positivity) (by positivity) (by nlinarith)
      linarith
  have hsplit : E * (Rint + P * Kint) - 1 / n / 2 * (σ ^ 2 * P)
      = E * Rint + P * (E * Kint - σ ^ 2 / (2 * n)) := by field_simp; ring
  rw [hsplit]
  have hPY : |P * (E * Kint - σ ^ 2 / (2 * n))| ≤ m * (σ ^ 3 / (2 * n ^ 2)) := by
    rw [abs_mul, abs_of_nonneg hP0]
    exact mul_le_mul hPm hY (abs_nonneg _) hm0
  have htot := (abs_add_le _ _).trans (add_le_add hX hPY)
  have e : m * σ ^ 3 / (2 * n ^ 2) + m * (σ ^ 3 / (2 * n ^ 2)) = (1 / n) ^ 2 * (σ ^ 3 * m) := by
    field_simp; ring
  have h0 : 0 ≤ (1 / n) ^ 2 * (σ ^ 3 * m) := by positivity
  linarith

end X06

theorem X06.proof_9 : X06.stmt_9 := by
  intro n hn
  have hn0 : 0 < n := by linarith
  have hN := X06.NN_pos n hn0
  have hN1 := X06.NN_ge_one n hn
  have hI7 := X06.proof_7 3 (X06.NN n) hN1
  refine ⟨?_, ?_⟩
  · have hT1 : (1 / n) * X06.T 1 (X06.NN n)
        = ∫ σ in Ioi (0:ℝ), (1 / n) / 2 * (σ ^ 2 * X06.psi (X06.W0 n σ)) := by
      rw [integral_const_mul]
      have := X06.ibp_psi (X06.NN n) hN
      simp only [X06.W0_eq]
      rw [this]; ring
    have hint2 : IntegrableOn (fun σ : ℝ => (1 / n) / 2 * (σ ^ 2 * X06.psi (X06.W0 n σ))) (Ioi 0) :=
      (X06.psi_sub_integrable 2 (X06.NN n) hN).const_mul _
    have hdiff : X06.Delta1 n - (1 / n) * X06.T 1 (X06.NN n)
        = ∫ σ in Ioi (0:ℝ), (X06.D1inner n σ - (1 / n) / 2 * (σ ^ 2 * X06.psi (X06.W0 n σ))) := by
      rw [hT1, integral_sub (X06.D1inner_integrable n hn0) hint2]; rfl
    rw [hdiff]
    have hG : IntegrableOn (fun σ : ℝ => 13 / 12 * (1 / n) ^ 2 * (σ ^ 3 * X06.mmaj (X06.NN n) σ))
        (Ioi 0) := (X06.mmaj_integrable 3 (X06.NN n) hN).const_mul _
    have := norm_integral_le_of_norm_le hG (by
      filter_upwards [ae_restrict_mem measurableSet_Ioi] with σ hσ
      rw [Real.norm_eq_abs]
      exact X06.D1_pointwise n σ hn (le_of_lt hσ))
    rw [Real.norm_eq_abs, integral_const_mul] at this
    unfold X06.Ik
    linarith
  · have h7 := hI7.2.1
    have : ((Nat.factorial 3 : ℕ) : ℝ) = 6 := by norm_num [Nat.factorial]
    rw [this] at h7
    have hp : 0 ≤ (1 / n) ^ 2 := by positivity
    have hL : 0 ≤ (Real.log (X06.NN n) + 1) ^ 3 := by
      have := Real.log_nonneg hN1; positivity
    calc 13 / 12 * (1 / n) ^ 2 * X06.Ik 3 (X06.NN n)
        ≤ 13 / 12 * (1 / n) ^ 2 * ((6 + 1 / 2) * (Real.log (X06.NN n) + 1) ^ 3) :=
          mul_le_mul_of_nonneg_left h7 (by positivity)
      _ ≤ 7.1 * (1 / n) ^ 2 * (Real.log (X06.NN n) + 1) ^ 3 := by
          have : 13 / 12 * (1 / n) ^ 2 * ((6 + 1 / 2) * (Real.log (X06.NN n) + 1) ^ 3)
              = (169 / 24) * ((1 / n) ^ 2 * (Real.log (X06.NN n) + 1) ^ 3) := by ring
          rw [this]
          have h0 : 0 ≤ (1 / n) ^ 2 * (Real.log (X06.NN n) + 1) ^ 3 := by positivity
          nlinarith

namespace X06

lemma stmt11_of (h5 : stmt_5) : stmt_11 := by
  intro n hn
  have hn1 : 1 < n := by linarith
  have hn0 : 0 < n := by linarith
  have hN1 := NN_ge_one n hn1
  have hN0 := NN_pos n hn0
  obtain ⟨-, θ1, hθ1, hT1⟩ := h5 1 (NN n) hN1
  obtain ⟨-, θ2, hθ2, hT2⟩ := h5 2 (NN n) hN1
  obtain ⟨-, -, -, hC, hC0⟩ := proof_8 n hn1
  obtain ⟨hD1, hD1b⟩ := proof_9 n hn1
  obtain ⟨hD2a, hD2b, hD2c⟩ := proof_10 n hn1
  simp only [Nat.factorial, pow_one, Nat.succ_eq_add_one] at hT1 hT2
  norm_num at hT1 hT2
  have hL : 0 ≤ Real.log (NN n) := Real.log_nonneg hN1
  have he0 : 0 < 1 / n := by positivity
  have he1 : 1 / n ≤ 1 / 2 := by rw [div_le_div_iff₀ hn0 (by norm_num)]; linarith
  have hs0 : 0 < Real.sqrt 2 := by positivity
  have hs14 : (1.4:ℝ) ≤ Real.sqrt 2 := by
    rw [Real.le_sqrt (by norm_num) (by norm_num)]; norm_num
  -- `u = 1/N = e/√2`
  have hu : 1 / NN n = (1 / n) / Real.sqrt 2 := by unfold NN; field_simp
  have hu0 : 0 < 1 / NN n := by positivity
  have hue : 1 / NN n ≤ (1 / n) / 1.4 := by
    rw [hu]; exact div_le_div_of_nonneg_left he0.le (by norm_num) hs14
  have h2N : 1 / (2 * NN n) = 1 / NN n / 2 := by field_simp
  have h2N' : 2 / (2 * NN n) = 1 / NN n := by field_simp
  have h72 : θ1 / (72 * NN n ^ 3) = θ1 * (1 / NN n) ^ 3 / 72 := by field_simp
  have h216 : θ2 * 2 / (216 * NN n ^ 3) = θ2 * (1 / NN n) ^ 3 / 108 := by field_simp; ring
  rw [neg_div, h2N] at hT1
  rw [h2N', neg_div, h216] at hT2
  rw [h72] at hT1
  have hdiff : Cint n - (poly 1 (Real.log (NN n)) - (1 / n) / (2 * Real.sqrt 2)
        + (1 / n) * (poly 1 (Real.log (NN n)) - poly 2 (Real.log (NN n))))
      = ((1 / n) / (2 * Real.sqrt 2) - 1 / NN n / 2)
        + (θ1 * (1 / NN n) ^ 3 / 72 + (1 / n) * (θ1 * (1 / NN n) ^ 3 / 72)
            + (1 / n) * (θ2 * (1 / NN n) ^ 3 / 108))
        - 3 / 2 * (1 / n) * (1 / NN n)
        + (Delta1 n - (1 / n) * T 1 (NN n)) + (Delta2 n + (1 / n) * T 2 (NN n)) := by
    rw [hC, hC0]
    rw [hT1, hT2]
    ring
  have hzero : (1 / n) / (2 * Real.sqrt 2) - 1 / NN n / 2 = 0 := by rw [hu]; field_simp; ring
  rw [hdiff, hzero, zero_add]
  set e := 1 / n
  set u := 1 / NN n
  set L := Real.log (NN n)
  have he1' : e ≤ 1 := by linarith
  have hu3 : u ^ 3 ≤ e ^ 2 / 2 := by
    have h1 : u ≤ e := hue.trans (div_le_self he0.le (by norm_num))
    have h2 : u ^ 3 ≤ e ^ 3 := pow_le_pow_left₀ hu0.le h1 3
    have h3 : e ^ 3 ≤ e ^ 2 / 2 := by
      calc e ^ 3 = e * e ^ 2 := by ring
        _ ≤ 1 / 2 * e ^ 2 := mul_le_mul_of_nonneg_right he1 (sq_nonneg e)
        _ = e ^ 2 / 2 := by ring
    linarith
  have hu30 : 0 ≤ u ^ 3 := by positivity
  have hX0 : 0 ≤ θ1 * u ^ 3 / 72 + e * (θ1 * u ^ 3 / 72) + e * (θ2 * u ^ 3 / 108) := by
    have := hθ1.1; have := hθ2.1; positivity
  have hX1 : θ1 * u ^ 3 / 72 + e * (θ1 * u ^ 3 / 72) + e * (θ2 * u ^ 3 / 108) ≤ e ^ 2 / 50 := by
    have b1 : 0 ≤ θ1 * u ^ 3 := mul_nonneg hθ1.1 hu30
    have b2 : 0 ≤ θ2 * u ^ 3 := mul_nonneg hθ2.1 hu30
    have a1 : θ1 * u ^ 3 ≤ u ^ 3 := mul_le_of_le_one_left hu30 hθ1.2
    have a2 : θ2 * u ^ 3 ≤ u ^ 3 := mul_le_of_le_one_left hu30 hθ2.2
    have a3 : e * (θ1 * u ^ 3) ≤ θ1 * u ^ 3 := mul_le_of_le_one_left b1 he1'
    have a4 : e * (θ2 * u ^ 3) ≤ θ2 * u ^ 3 := mul_le_of_le_one_left b2 he1'
    have : θ1 * u ^ 3 / 72 + e * (θ1 * u ^ 3 / 72) + e * (θ2 * u ^ 3 / 108)
        = θ1 * u ^ 3 / 72 + e * (θ1 * u ^ 3) / 72 + e * (θ2 * u ^ 3) / 108 := by ring
    rw [this]
    linarith
  have hY0 : 0 ≤ 3 / 2 * e * u := by positivity
  have hY1 : 3 / 2 * e * u ≤ 1.08 * e ^ 2 := by
    have h := mul_le_mul_of_nonneg_left hue (by positivity : (0:ℝ) ≤ 3 / 2 * e)
    have e2 : 3 / 2 * e * (e / 1.4) = 15 / 14 * e ^ 2 := by ring
    have : 0 ≤ e ^ 2 := sq_nonneg e
    linarith
  have hr1 := abs_le.1 (hD1.trans hD1b)
  have hr3 : |Delta2 n + e * T 2 (NN n)| ≤ 3.5 * e ^ 2 * (L + 1) ^ 4 := by
    rw [abs_of_nonneg hD2a]; exact hD2b.trans hD2c
  have hr3' := abs_le.1 hr3
  have hL1 : 1 ≤ (L + 1) ^ 4 := one_le_pow₀ (by linarith)
  have hL34 : (L + 1) ^ 3 ≤ (L + 1) ^ 4 := pow_le_pow_right₀ (by linarith) (by norm_num)
  have he2 : 0 ≤ e ^ 2 := sq_nonneg e
  have hb1 : e ^ 2 * (L + 1) ^ 3 ≤ e ^ 2 * (L + 1) ^ 4 := mul_le_mul_of_nonneg_left hL34 he2
  have hb2 : e ^ 2 * 1 ≤ e ^ 2 * (L + 1) ^ 4 := mul_le_mul_of_nonneg_left hL1 he2
  have hr1a : |Delta1 n - e * T 1 (NN n)| ≤ 7.1 * (e ^ 2 * (L + 1) ^ 3) := by
    have := hD1.trans hD1b; linarith
  have hr1' := abs_le.1 hr1a
  clear hr1 hT1 hT2 hC hC0 hdiff hD1 hD1b hD2b hD2c hr3
  rw [abs_le]; constructor <;> linarith

end X06

namespace X06

lemma poly1_eq (L : ℝ) : poly 1 L = L ^ 2 / 2 + mt 1 * L + mt 2 := by
  simp [poly, Finset.sum_range_succ, Nat.descFactorial]; ring

lemma poly2_eq (L : ℝ) : poly 2 L = L ^ 3 / 3 + mt 1 * L ^ 2 + 2 * mt 2 * L + 2 * mt 3 := by
  simp [poly, Finset.sum_range_succ, Nat.descFactorial]; norm_num; ring

lemma stmt13_of_stmt3 (h : stmt_3) : stmt_13 := by
  obtain ⟨h0, h1, h2, h3⟩ := h
  have m1 : mt 1 = 1 - lam := by simp [mt, Finset.sum_range_succ, h0, h1]; ring
  have m2 : mt 2 = 1 - lam + (lam ^ 2 / 2 + Real.pi ^ 2 / 12) := by
    simp [mt, Finset.sum_range_succ, h0, h1, h2]; ring
  have m3 : mt 3 = 1 - lam + (lam ^ 2 / 2 + Real.pi ^ 2 / 12)
      + (-lam ^ 3 / 6 - lam * Real.pi ^ 2 / 12 - z3 / 4) := by
    simp [mt, Finset.sum_range_succ, h0, h1, h2, h3]; ring
  intro l
  rw [poly1_eq, poly2_eq, m1, m2, m3]
  unfold beta1 beta0 b3 c3 d3
  constructor <;> ring

end X06

namespace X06

lemma stmt12_of (h1 : stmt_1) (h11 : stmt_11) : stmt_12 := by
  intro n hn
  have hn0 : 0 < n := by linarith
  obtain ⟨-, -, hA⟩ := h1 n (by linarith)
  have h2 := proof_2 n hn
  have hC := h11 n hn
  set s := Real.sqrt 2 with hs
  have hs2 : s ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have hs0 : 0 < s := by positivity
  set p1 := poly 1 (Real.log (NN n))
  set p2 := poly 2 (Real.log (NN n))
  set L := Real.log (NN n)
  have hz : (1 / n) ^ 3 * (1 / (2 * s) - s / 4) = 0 := by
    have : 1 / (2 * s) - s / 4 = 0 := by
      field_simp; linarith
    rw [this, mul_zero]
  have key : A n - 1 - (s * (1 / n) + (-s - p1) * (1 / n) ^ 2 + (3 * s / 2 + p2 - p1) * (1 / n) ^ 3)
      = (Phi n / (1 + 1 / n / 2) - (s / n - s / n ^ 2 + 5 * s / 4 / n ^ 3))
        - (1 / n) ^ 2 * (Cint n - (p1 - (1 / n) / (2 * s) + (1 / n) * (p1 - p2)))
        + (1 / n) ^ 3 * (1 / (2 * s) - s / 4) := by
    rw [hA]; field_simp; ring
  rw [key, hz, add_zero]
  have hn4 : (1 / n) ^ 4 = 1 / n ^ 4 := by rw [div_pow, one_pow]
  have e1 : (3 * s / 2 + 1 / 2) / n ^ 4 = (3 * s / 2 + 1 / 2) * (1 / n) ^ 4 := by
    rw [hn4]; ring
  rw [e1] at h2
  have hp : 0 ≤ (1 / n) ^ 2 := by positivity
  have h3 : |(1 / n) ^ 2 * (Cint n - (p1 - (1 / n) / (2 * s) + (1 / n) * (p1 - p2)))|
      ≤ 12 * (L + 1) ^ 4 * (1 / n) ^ 4 := by
    rw [abs_mul, abs_of_nonneg hp]
    calc (1 / n) ^ 2 * |Cint n - (p1 - (1 / n) / (2 * s) + (1 / n) * (p1 - p2))|
        ≤ (1 / n) ^ 2 * (12 * (1 / n) ^ 2 * (L + 1) ^ 4) := mul_le_mul_of_nonneg_left hC hp
      _ = 12 * (L + 1) ^ 4 * (1 / n) ^ 4 := by ring
  calc |(Phi n / (1 + 1 / n / 2) - (s / n - s / n ^ 2 + 5 * s / 4 / n ^ 3))
        - (1 / n) ^ 2 * (Cint n - (p1 - (1 / n) / (2 * s) + (1 / n) * (p1 - p2)))|
      ≤ |Phi n / (1 + 1 / n / 2) - (s / n - s / n ^ 2 + 5 * s / 4 / n ^ 3)|
        + |(1 / n) ^ 2 * (Cint n - (p1 - (1 / n) / (2 * s) + (1 / n) * (p1 - p2)))| := abs_sub _ _
    _ ≤ (3 * s / 2 + 1 / 2) * (1 / n) ^ 4 + 12 * (L + 1) ^ 4 * (1 / n) ^ 4 := add_le_add h2 h3
    _ = (3 * s / 2 + 1 / 2 + 12 * (L + 1) ^ 4) * (1 / n) ^ 4 := by ring

lemma log_NN (n : ℝ) (hn : 0 < n) : Real.log (NN n) = Real.log n + lam / 2 := by
  unfold NN lam
  rw [Real.log_mul (by positivity) hn.ne', Real.log_sqrt (by norm_num)]; ring

lemma stmt15_of (h12 : stmt_12) (h13 : stmt_13) : stmt_15 := by
  intro n hn
  have hn0 : 0 < n := by linarith
  have h := h12 n hn
  obtain ⟨hR, hR3⟩ := h13 (Real.log n)
  rw [log_NN n hn0] at h
  rw [hR, hR3] at h
  set l := Real.log n
  have hl : 0 ≤ l := Real.log_nonneg (by linarith)
  have hlam : lam / 2 ≤ 0.35 := by
    unfold lam; have := Real.log_two_lt_d9; linarith
  have hlam0 : 0 ≤ lam := by unfold lam; exact Real.log_nonneg (by norm_num)
  have hs : Real.sqrt 2 < 1.415 := by
    rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have hR' : R n = -l ^ 2 / 2 + beta1 * l + beta0 := rfl
  have hR3' : R3 n = l ^ 3 / 3 + b3 * l ^ 2 + c3 * l + d3 := rfl
  have heq : A n - (1 + Real.sqrt 2 / n + R n / n ^ 2 + R3 n / n ^ 3)
      = A n - 1 - (Real.sqrt 2 * (1 / n) + (-l ^ 2 / 2 + beta1 * l + beta0) * (1 / n) ^ 2
        + (l ^ 3 / 3 + b3 * l ^ 2 + c3 * l + d3) * (1 / n) ^ 3) := by
    rw [hR', hR3']; field_simp; ring
  rw [heq]
  refine h.trans ?_
  have hb : l + lam / 2 + 1 ≤ 1.35 * (1 + l) := by linarith
  have hb4 : (l + lam / 2 + 1) ^ 4 ≤ (1.35 * (1 + l)) ^ 4 :=
    pow_le_pow_left₀ (by linarith) hb 4
  have h1l : 1 ≤ (1 + l) ^ 4 := one_le_pow₀ (by linarith)
  have hn4 : (1 / n) ^ 4 = 1 / n ^ 4 := by rw [div_pow, one_pow]
  rw [hn4, mul_one_div]
  apply div_le_div_of_nonneg_right _ (by positivity)
  have : (1.35 * (1 + l)) ^ 4 = 3.32150625 * (1 + l) ^ 4 := by ring
  nlinarith

lemma stmt16_of (h15 : stmt_15) : stmt_16 := ⟨45, h15⟩

lemma stmt14_of (h16 : stmt_16) : stmt_14 := by
  obtain ⟨C, hC⟩ := h16
  refine ⟨|C| + 1 / 3 + |b3| + |c3| + |d3|, fun n hn => ?_⟩
  have hn0 : 0 < n := by linarith
  have h := hC n hn
  set l := Real.log n
  have hl : 0 ≤ l := Real.log_nonneg (by linarith)
  have hln : 1 + l ≤ n := by have := Real.log_le_sub_one_of_pos hn0; linarith
  have h1l : 1 ≤ 1 + l := by linarith
  have hR3 : |R3 n| ≤ (1 / 3 + |b3| + |c3| + |d3|) * (1 + l) ^ 3 := by
    have hR3' : R3 n = l ^ 3 / 3 + b3 * l ^ 2 + c3 * l + d3 := rfl
    rw [hR3']
    have p3 : l ^ 3 ≤ (1 + l) ^ 3 := pow_le_pow_left₀ hl (by linarith) 3
    have p2 : l ^ 2 ≤ (1 + l) ^ 3 := by
      calc l ^ 2 ≤ (1 + l) ^ 2 := pow_le_pow_left₀ hl (by linarith) 2
        _ ≤ (1 + l) ^ 3 := pow_le_pow_right₀ h1l (by norm_num)
    have p1 : l ≤ (1 + l) ^ 3 := by
      calc l ≤ 1 + l := by linarith
        _ = (1 + l) ^ 1 := (pow_one _).symm
        _ ≤ (1 + l) ^ 3 := pow_le_pow_right₀ h1l (by norm_num)
    have p0 : 1 ≤ (1 + l) ^ 3 := one_le_pow₀ h1l
    calc |l ^ 3 / 3 + b3 * l ^ 2 + c3 * l + d3|
        ≤ |l ^ 3 / 3| + |b3 * l ^ 2| + |c3 * l| + |d3| := by
          have t1 := abs_add_le (l ^ 3 / 3 + b3 * l ^ 2 + c3 * l) d3
          have t2 := abs_add_le (l ^ 3 / 3 + b3 * l ^ 2) (c3 * l)
          have t3 := abs_add_le (l ^ 3 / 3) (b3 * l ^ 2)
          linarith
      _ ≤ 1 / 3 * (1 + l) ^ 3 + |b3| * (1 + l) ^ 3 + |c3| * (1 + l) ^ 3 + |d3| * (1 + l) ^ 3 := by
          rw [abs_div, abs_mul, abs_mul, abs_of_nonneg (pow_nonneg hl 3),
            abs_of_nonneg (pow_nonneg hl 2), abs_of_nonneg hl]
          have q1 : l ^ 3 / |(3:ℝ)| ≤ 1 / 3 * (1 + l) ^ 3 := by
            rw [abs_of_pos (by norm_num : (0:ℝ) < 3)]; linarith
          have q2 := mul_le_mul_of_nonneg_left p2 (abs_nonneg b3)
          have q3 := mul_le_mul_of_nonneg_left p1 (abs_nonneg c3)
          have q4 := mul_le_mul_of_nonneg_left p0 (abs_nonneg d3)
          linarith
      _ = (1 / 3 + |b3| + |c3| + |d3|) * (1 + l) ^ 3 := by ring
  have heq : A n - (1 + Real.sqrt 2 / n + R n / n ^ 2)
      = (A n - (1 + Real.sqrt 2 / n + R n / n ^ 2 + R3 n / n ^ 3)) + R3 n / n ^ 3 := by ring
  rw [heq]
  have hn3 : 0 < n ^ 3 := by positivity
  have t1 : C * (1 + l) ^ 4 / n ^ 4 ≤ |C| * (1 + l) ^ 3 / n ^ 3 := by
    have hq : (1 + l) / n ≤ 1 := (div_le_one hn0).2 hln
    have : C * (1 + l) ^ 4 / n ^ 4 = C * ((1 + l) ^ 3 / n ^ 3) * ((1 + l) / n) := by
      field_simp
    rw [this]
    have h0 : 0 ≤ (1 + l) ^ 3 / n ^ 3 := by positivity
    have h0' : 0 ≤ (1 + l) / n := by positivity
    calc C * ((1 + l) ^ 3 / n ^ 3) * ((1 + l) / n) ≤ |C| * ((1 + l) ^ 3 / n ^ 3) * ((1 + l) / n) := by
          apply mul_le_mul_of_nonneg_right _ h0'
          exact mul_le_mul_of_nonneg_right (le_abs_self C) h0
      _ ≤ |C| * ((1 + l) ^ 3 / n ^ 3) * 1 :=
          mul_le_mul_of_nonneg_left hq (mul_nonneg (abs_nonneg C) h0)
      _ = |C| * (1 + l) ^ 3 / n ^ 3 := by ring
  have t2 : |R3 n / n ^ 3| ≤ (1 / 3 + |b3| + |c3| + |d3|) * (1 + l) ^ 3 / n ^ 3 := by
    rw [abs_div, abs_of_pos hn3]
    exact div_le_div_of_nonneg_right hR3 hn3.le
  calc |(A n - (1 + Real.sqrt 2 / n + R n / n ^ 2 + R3 n / n ^ 3)) + R3 n / n ^ 3|
      ≤ |A n - (1 + Real.sqrt 2 / n + R n / n ^ 2 + R3 n / n ^ 3)| + |R3 n / n ^ 3| := abs_add_le _ _
    _ ≤ |C| * (1 + l) ^ 3 / n ^ 3 + (1 / 3 + |b3| + |c3| + |d3|) * (1 + l) ^ 3 / n ^ 3 :=
        add_le_add (h.trans t1) t2
    _ = (|C| + 1 / 3 + |b3| + |c3| + |d3|) * (1 + l) ^ 3 / n ^ 3 := by ring

end X06


/-! ### Mellin transform of `g` (`stmt_4`) -/

namespace X06
open Real MeasureTheory Set Filter Topology

lemma beta_real' {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    (∫ t in (0 : ℝ)..1, t ^ (a - 1) * (1 - t) ^ (b - 1)) * Real.Gamma (a + b) =
      Real.Gamma a * Real.Gamma b := by
  have key := Complex.Gamma_mul_Gamma_eq_betaIntegral (s := (a : ℂ)) (t := (b : ℂ))
    (by simpa using ha) (by simpa using hb)
  have hB : Complex.betaIntegral a b =
      ((∫ t in (0 : ℝ)..1, t ^ (a - 1) * (1 - t) ^ (b - 1) : ℝ) : ℂ) := by
    rw [Complex.betaIntegral, ← intervalIntegral.integral_ofReal]
    refine intervalIntegral.integral_congr_uIoo fun t ht => ?_
    rw [uIoo_of_le zero_le_one] at ht
    push_cast
    rw [Complex.ofReal_cpow ht.1.le, Complex.ofReal_cpow (by linarith [ht.2] : (0:ℝ) ≤ 1 - t)]
    push_cast
    ring_nf
  rw [hB, ← Complex.ofReal_add, Complex.Gamma_ofReal, Complex.Gamma_ofReal,
    Complex.Gamma_ofReal] at key
  have k2 := key.symm
  rw [mul_comm] at k2
  exact_mod_cast k2

/-- the substitution `t = w²/(1+w²)` from `(0,∞)` onto `(0,1)`. -/
lemma subst_beta (G : ℝ → ℝ) :
    (IntegrableOn G (Ioo 0 1) ↔
      IntegrableOn (fun w => 2 * w / (1 + w ^ 2) ^ 2 * G (w ^ 2 / (1 + w ^ 2))) (Ioi 0)) ∧
      ∫ t in Ioo 0 1, G t = ∫ w in Ioi (0:ℝ), 2 * w / (1 + w ^ 2) ^ 2 * G (w ^ 2 / (1 + w ^ 2)) := by
  have himg : (fun w : ℝ => w ^ 2 / (1 + w ^ 2)) '' Ioi 0 = Ioo 0 1 := by
    ext t; constructor
    · rintro ⟨w, hw, rfl⟩
      have hw' : (0:ℝ) < w := hw
      refine ⟨by positivity, ?_⟩
      rw [div_lt_one (by positivity)]; linarith
    · rintro ⟨h0, h1⟩
      refine ⟨Real.sqrt (t / (1 - t)), Real.sqrt_pos.2 (div_pos h0 (by linarith)), ?_⟩
      simp only
      rw [Real.sq_sqrt (div_pos h0 (by linarith)).le]
      have : (1 - t) ≠ 0 := by linarith
      field_simp
      ring
  have hd : ∀ x ∈ Ioi (0:ℝ), HasDerivWithinAt (fun w : ℝ => w ^ 2 / (1 + w ^ 2))
      (2 * x / (1 + x ^ 2) ^ 2) (Ioi 0) x := by
    intro x _
    have hq : (1 + x ^ 2) ≠ 0 := by positivity
    have := (hasDerivAt_pow 2 x).div ((hasDerivAt_pow 2 x).const_add 1) hq
    convert this.hasDerivWithinAt using 1
    simp only [Nat.cast_ofNat]; field_simp; ring
  have hinj : InjOn (fun w : ℝ => w ^ 2 / (1 + w ^ 2)) (Ioi 0) := by
    intro x hx y hy h
    simp only at h
    have hx' : (0:ℝ) < x := hx
    have hy' : (0:ℝ) < y := hy
    rw [div_eq_div_iff (by positivity) (by positivity)] at h
    have : x ^ 2 = y ^ 2 := by nlinarith
    nlinarith [sq_nonneg (x - y), sq_nonneg (x + y)]
  have habs : ∀ x ∈ Ioi (0:ℝ), |2 * x / (1 + x ^ 2) ^ 2| = 2 * x / (1 + x ^ 2) ^ 2 := by
    intro x hx; have : (0:ℝ) < x := hx; rw [abs_of_pos (by positivity)]
  constructor
  · rw [← himg, integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ioi hd hinj]
    refine integrableOn_congr_fun (fun x hx => ?_) measurableSet_Ioi
    simp only [habs x hx, smul_eq_mul]
  · rw [← himg, integral_image_eq_integral_abs_deriv_smul measurableSet_Ioi hd hinj]
    refine setIntegral_congr_fun measurableSet_Ioi (fun x hx => ?_)
    simp only [habs x hx, smul_eq_mul]

lemma beta_pointwise (z w : ℝ) (hw : 0 < w) :
    2 * w / (1 + w ^ 2) ^ 2 * ((w ^ 2 / (1 + w ^ 2)) ^ ((1 + z) / 2 - 1) *
      (1 - w ^ 2 / (1 + w ^ 2)) ^ ((1 - z / 2) - 1))
      = 2 * (w ^ z * (1 + w ^ 2) ^ (-(3 / 2 : ℝ))) := by
  have hq : 0 < 1 + w ^ 2 := by positivity
  have h1 : 1 - w ^ 2 / (1 + w ^ 2) = (1 + w ^ 2)⁻¹ := by field_simp; ring
  rw [h1]
  have hw2 : w ^ 2 = w ^ (2:ℝ) := by norm_cast
  generalize hqdef : 1 + w ^ 2 = q at *
  rw [Real.div_rpow (by positivity) hq.le, Real.inv_rpow hq.le, hw2, ← Real.rpow_mul hw.le]
  rw [Real.rpow_def_of_pos hw, Real.rpow_def_of_pos hw, Real.rpow_def_of_pos hq,
    Real.rpow_def_of_pos hq, Real.rpow_def_of_pos hq]
  have ew : w = Real.exp (Real.log w) := (Real.exp_log hw).symm
  have eq2 : q ^ 2 = Real.exp (Real.log q * 2) := by
    rw [Real.exp_mul, Real.exp_log hq]; norm_cast
  rw [eq2]
  conv_lhs => rw [show 2 * w = 2 * Real.exp (Real.log w) by rw [← ew]]
  rw [mul_div_assoc, ← Real.exp_sub, div_eq_mul_inv, ← Real.exp_neg,
    ← Real.exp_neg, mul_assoc, mul_assoc, ← Real.exp_add, ← Real.exp_add, ← Real.exp_add,
    ← Real.exp_add]
  congr 2
  ring

lemma K_integral (z : ℝ) (hz0 : 0 < z) (hz1 : z < 1) :
    IntegrableOn (fun w : ℝ => w ^ z * (1 + w ^ 2) ^ (-(3 / 2 : ℝ))) (Ioi 0) ∧
    ∫ w in Ioi (0:ℝ), w ^ z * (1 + w ^ 2) ^ (-(3 / 2 : ℝ)) = Eg z := by
  set a : ℝ := (1 + z) / 2
  set b : ℝ := 1 - z / 2
  have ha : 0 < a := by positivity
  have hb : 0 < b := by simp only [b]; linarith
  set G : ℝ → ℝ := fun t => t ^ (a - 1) * (1 - t) ^ (b - 1)
  have hab : a + b = 3 / 2 := by simp only [a, b]; ring
  have key := beta_real' ha hb
  rw [hab] at key
  have hG32 : Real.Gamma (3 / 2) = Real.sqrt Real.pi / 2 := by
    rw [show (3 / 2 : ℝ) = 1 / 2 + 1 by norm_num, Real.Gamma_add_one (by norm_num),
      Real.Gamma_one_half_eq]; ring
  rw [hG32] at key
  have hpi : 0 < Real.sqrt Real.pi := Real.sqrt_pos.2 Real.pi_pos
  -- integrability of G on (0,1)
  have hGint : IntegrableOn G (Ioo 0 1) := by
    have hconv := Complex.betaIntegral_convergent (u := (a : ℂ)) (v := (b : ℂ))
      (by simpa using ha) (by simpa using hb)
    have hI : IntervalIntegrable G volume 0 1 := by
      have h2 : IntervalIntegrable (fun t : ℝ => ((G t : ℝ) : ℂ)) volume 0 1 := by
        refine hconv.congr ?_
        intro t ht
        rw [uIoc_of_le zero_le_one] at ht
        simp only [G]
        push_cast
        rw [Complex.ofReal_cpow ht.1.le, Complex.ofReal_cpow (by linarith [ht.2] : (0:ℝ) ≤ 1 - t)]
        push_cast
        ring_nf
      exact ⟨by simpa using h2.1.re, by simp⟩
    exact (intervalIntegrable_iff_integrableOn_Ioo_of_le zero_le_one).1 hI
  obtain ⟨hiff, heq⟩ := subst_beta G
  have hpt : ∀ w ∈ Ioi (0:ℝ), 2 * w / (1 + w ^ 2) ^ 2 * G (w ^ 2 / (1 + w ^ 2))
      = 2 * (w ^ z * (1 + w ^ 2) ^ (-(3 / 2 : ℝ))) := fun w hw => beta_pointwise z w hw
  have hint2 := (hiff.1 hGint).congr_fun hpt measurableSet_Ioi
  refine ⟨?_, ?_⟩
  · have := hint2.const_mul (1 / 2)
    exact IntegrableOn.congr_fun (f := fun w => 1 / 2 * (2 * (w ^ z * (1 + w ^ 2) ^ (-(3 / 2 : ℝ)))))
      this (fun w _ => by ring) measurableSet_Ioi
  · rw [setIntegral_congr_fun measurableSet_Ioi hpt, integral_const_mul] at heq
    have h01 : ∫ t in (0:ℝ)..1, G t = ∫ t in Ioo 0 1, G t := by
      rw [intervalIntegral.integral_of_le zero_le_one, integral_Ioc_eq_integral_Ioo]
    rw [h01, heq] at key
    unfold Eg
    rw [eq_div_iff hpi.ne']
    linarith
/-! ### Mellin transform of `g` (stmt_4) -/

/-- `f'(w) = 1 - w/√(1+w²)`. -/
noncomputable def fp (w : ℝ) : ℝ := 1 - w / Real.sqrt (1 + w ^ 2)

lemma fp_continuous : Continuous fp := by
  unfold fp
  exact continuous_const.sub (continuous_id.div (by fun_prop)
    (fun w => (Real.sqrt_pos.2 (by positivity)).ne'))

lemma fp_bounds (w : ℝ) (hw : 0 < w) : 0 ≤ fp w ∧ fp w ≤ 1 ∧ fp w ≤ 1 / w ^ 2 := by
  unfold fp
  set s := Real.sqrt (1 + w ^ 2) with hsdef
  have hs : 0 < s := Real.sqrt_pos.2 (by positivity)
  have hss : s ^ 2 = 1 + w ^ 2 := Real.sq_sqrt (by positivity)
  have hws : w < s := by nlinarith
  refine ⟨?_, ?_, ?_⟩
  · rw [sub_nonneg, div_le_one hs]; exact hws.le
  · have : 0 ≤ w / s := by positivity
    linarith
  · have e : 1 - w / s = 1 / (s * (s + w)) := by
      field_simp; nlinarith
    rw [e]
    apply one_div_le_one_div_of_le (by positivity)
    nlinarith

lemma rpow32_eq (q : ℝ) (hq : 0 < q) : q ^ (-(3 / 2 : ℝ)) = 1 / (q * Real.sqrt q) := by
  rw [Real.rpow_neg hq.le, Real.sqrt_eq_rpow, ← Real.rpow_one_add' hq.le (by norm_num)]
  norm_num

lemma hasDerivAt_fp (w : ℝ) :
    HasDerivAt fp (-(1 + w ^ 2) ^ (-(3 / 2 : ℝ))) w := by
  have h1 : (1 + w ^ 2) ≠ 0 := by positivity
  have hq : 0 < 1 + w ^ 2 := by positivity
  have hs : 0 < Real.sqrt (1 + w ^ 2) := Real.sqrt_pos.2 hq
  have hss : Real.sqrt (1 + w ^ 2) ^ 2 = 1 + w ^ 2 := Real.sq_sqrt hq.le
  have h := (hasDerivAt_const w (1:ℝ)).sub ((hasDerivAt_id w).div
    (((hasDerivAt_pow 2 w).const_add 1).sqrt h1) hs.ne')
  convert h using 1
  · funext x; simp [fp]
  rw [rpow32_eq _ hq]
  simp only [Nat.cast_ofNat, id, pow_one, Nat.add_one_sub_one]
  field_simp
  rw [hss]; ring

lemma fp_le_one_div_sq (w : ℝ) (hw : 0 < w) (z : ℝ) : w ^ z * (1 / w ^ 2) = w ^ (z - 2) := by
  rw [Real.rpow_sub hw, show (2:ℝ) = ((2:ℕ):ℝ) by norm_num, Real.rpow_natCast]; ring

/-- an integrability criterion on `(0,∞)` from power bounds near `0` and near `∞`. -/
lemma int_pieces (h : ℝ → ℝ) (hc : ContinuousOn h (Ioi 0)) {p r : ℝ} (hp : -1 < p) (hr : r < -1)
    (h1 : ∀ w ∈ Ioc (0:ℝ) 1, |h w| ≤ w ^ p) (h2 : ∀ w ∈ Ioi (1:ℝ), |h w| ≤ w ^ r) :
    IntegrableOn h (Ioi 0) := by
  rw [← Ioc_union_Ioi_eq_Ioi zero_le_one]
  refine IntegrableOn.union ?_ ?_
  · have hi : IntegrableOn (fun x : ℝ => x ^ p) (Ioc 0 1) :=
      (intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one).1
        (intervalIntegral.intervalIntegrable_rpow' hp)
    refine hi.mono' ((hc.mono Ioc_subset_Ioi_self).aestronglyMeasurable measurableSet_Ioc) ?_
    filter_upwards [ae_restrict_mem measurableSet_Ioc] with w hw
    rw [Real.norm_eq_abs]; exact h1 w hw
  · have hi : IntegrableOn (fun x : ℝ => x ^ r) (Ioi 1) := integrableOn_Ioi_rpow_of_lt hr one_pos
    refine hi.mono' ((hc.mono (Ioi_subset_Ioi zero_le_one)).aestronglyMeasurable
      measurableSet_Ioi) ?_
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with w hw
    rw [Real.norm_eq_abs]; exact h2 w hw

lemma rpow_contOn (c : ℝ) : ContinuousOn (fun w : ℝ => w ^ c) (Ioi 0) :=
  fun w hw => (Real.continuousAt_rpow_const w c (Or.inl (ne_of_gt hw))).continuousWithinAt

lemma int_i (z : ℝ) (hz0 : 0 < z) (hz1 : z < 1) :
    IntegrableOn (fun w : ℝ => f w * w ^ (z - 2)) (Ioi 0) := by
  refine int_pieces _ (f_continuous.continuousOn.mul (rpow_contOn _)) (p := z - 1) (r := z - 2)
    (by linarith) (by linarith) ?_ ?_
  · intro w hw
    have hw0 : 0 < w := hw.1
    rw [abs_of_nonneg (mul_nonneg (f_nonneg w hw0.le) (by positivity))]
    calc f w * w ^ (z - 2) ≤ w * w ^ (z - 2) :=
          mul_le_mul_of_nonneg_right (f_le_self w) (by positivity)
      _ = w ^ (z - 1) := by
          rw [show z - 1 = (z - 2) + 1 by ring, Real.rpow_add_one hw0.ne']; ring
  · intro w hw
    have hw0 : (0:ℝ) < w := lt_trans one_pos hw
    rw [abs_of_nonneg (mul_nonneg (f_nonneg w hw0.le) (by positivity))]
    calc f w * w ^ (z - 2) ≤ 1 * w ^ (z - 2) :=
          mul_le_mul_of_nonneg_right (f_le_one w) (by positivity)
      _ = w ^ (z - 2) := one_mul _

lemma int_ii (z : ℝ) (hz0 : 0 < z) (hz1 : z < 1) :
    IntegrableOn (fun w : ℝ => fp w * w ^ (z - 1)) (Ioi 0) := by
  refine int_pieces _ (fp_continuous.continuousOn.mul (rpow_contOn _)) (p := z - 1) (r := z - 3)
    (by linarith) (by linarith) ?_ ?_
  · intro w hw
    have hw0 : 0 < w := hw.1
    obtain ⟨b0, b1, -⟩ := fp_bounds w hw0
    rw [abs_of_nonneg (mul_nonneg b0 (by positivity))]
    calc fp w * w ^ (z - 1) ≤ 1 * w ^ (z - 1) := mul_le_mul_of_nonneg_right b1 (by positivity)
      _ = w ^ (z - 1) := one_mul _
  · intro w hw
    have hw0 : (0:ℝ) < w := lt_trans one_pos hw
    obtain ⟨b0, -, b2⟩ := fp_bounds w hw0
    rw [abs_of_nonneg (mul_nonneg b0 (by positivity))]
    calc fp w * w ^ (z - 1) ≤ 1 / w ^ 2 * w ^ (z - 1) :=
          mul_le_mul_of_nonneg_right b2 (by positivity)
      _ = w ^ (z - 1) * (1 / w ^ 2) := by ring
      _ = w ^ (z - 3) := by rw [fp_le_one_div_sq w hw0]; ring_nf

lemma tendsto_rpow_zero_right (z : ℝ) (hz : 0 < z) :
    Tendsto (fun w : ℝ => w ^ z) (𝓝[>] 0) (𝓝 0) := by
  have := (Real.continuousAt_rpow_const 0 z (Or.inr hz.le)).tendsto
  rw [Real.zero_rpow hz.ne'] at this
  exact this.mono_left nhdsWithin_le_nhds

/-- first integration by parts. -/
lemma ibp1 (z : ℝ) (hz0 : 0 < z) (hz1 : z < 1) :
    ∫ w in Ioi (0:ℝ), f w * w ^ (z - 2) =
      1 / (1 - z) * ∫ w in Ioi (0:ℝ), fp w * w ^ (z - 1) := by
  have hz1' : z - 1 ≠ 0 := by linarith
  have hu : ∀ x ∈ Ioi (0:ℝ), HasDerivAt f (fp x) x := fun x _ => hasDerivAt_f x
  have hv : ∀ x ∈ Ioi (0:ℝ), HasDerivAt (fun w : ℝ => w ^ (z - 1) / (z - 1)) (x ^ (z - 2)) x := by
    intro x hx
    have := (Real.hasDerivAt_rpow_const (p := z - 1) (Or.inl (ne_of_gt hx))).div_const (z - 1)
    convert this using 1
    rw [show z - 1 - 1 = z - 2 by ring]; field_simp
  have huv' : IntegrableOn (f * fun x : ℝ => x ^ (z - 2)) (Ioi 0) := int_i z hz0 hz1
  have hu'v : IntegrableOn (fp * fun w : ℝ => w ^ (z - 1) / (z - 1)) (Ioi 0) := by
    have : IntegrableOn _ (Ioi (0:ℝ)) := (int_ii z hz0 hz1).div_const (z - 1)
    refine this.congr_fun (fun x _ => ?_) measurableSet_Ioi
    simp only [Pi.mul_apply]; ring
  have hbd : ∀ w : ℝ, 0 < w → |f w * (w ^ (z - 1) / (z - 1))| ≤ w ^ z / (1 - z) ∧
      |f w * (w ^ (z - 1) / (z - 1))| ≤ w ^ (z - 1) / (1 - z) := by
    intro w hw
    have hf0 := f_nonneg w hw.le
    rw [abs_mul, abs_div, abs_of_nonneg hf0, abs_of_pos (by positivity : 0 < w ^ (z - 1)),
      abs_of_neg (by linarith : z - 1 < 0), neg_sub]
    constructor
    · calc f w * (w ^ (z - 1) / (1 - z)) ≤ w * (w ^ (z - 1) / (1 - z)) :=
            mul_le_mul_of_nonneg_right (f_le_self w) (div_nonneg (by positivity) (by linarith))
        _ = w ^ z / (1 - z) := by
            rw [Real.rpow_sub_one hw.ne']; field_simp
    · calc f w * (w ^ (z - 1) / (1 - z)) ≤ 1 * (w ^ (z - 1) / (1 - z)) :=
            mul_le_mul_of_nonneg_right (f_le_one w) (div_nonneg (by positivity) (by linarith))
        _ = w ^ (z - 1) / (1 - z) := one_mul _
  have h0 : Tendsto (f * fun w : ℝ => w ^ (z - 1) / (z - 1)) (𝓝[>] 0) (𝓝 0) := by
    have hlim := (tendsto_rpow_zero_right z hz0).div_const (1 - z)
    rw [zero_div] at hlim
    refine squeeze_zero_norm' ?_ hlim
    filter_upwards [self_mem_nhdsWithin] with w hw
    rw [Real.norm_eq_abs]; exact (hbd w hw).1
  have hinf : Tendsto (f * fun w : ℝ => w ^ (z - 1) / (z - 1)) atTop (𝓝 0) := by
    have hlim := (tendsto_rpow_neg_atTop (y := 1 - z) (by linarith)).div_const (1 - z)
    rw [zero_div] at hlim
    refine squeeze_zero_norm' ?_ hlim
    filter_upwards [eventually_gt_atTop 0] with w hw
    rw [Real.norm_eq_abs, neg_sub, show (z - 1) = -(1 - z) by ring, ← show (z - 1) = -(1 - z) by ring]
    exact (hbd w hw).2
  have := integral_Ioi_mul_deriv_eq_deriv_mul hu hv huv' hu'v h0 hinf
  rw [this]
  have e : ∫ x in Ioi (0:ℝ), fp x * (x ^ (z - 1) / (z - 1))
      = (1 / (z - 1)) * ∫ x in Ioi (0:ℝ), fp x * x ^ (z - 1) := by
    rw [← integral_const_mul]; congr 1; funext x; field_simp
  rw [e, show (1 : ℝ) / (1 - z) = -(1 / (z - 1)) by field_simp; ring]
  ring

/-- second integration by parts. -/
lemma ibp2 (z : ℝ) (hz0 : 0 < z) (hz1 : z < 1) :
    ∫ w in Ioi (0:ℝ), fp w * w ^ (z - 1) =
      1 / z * ∫ w in Ioi (0:ℝ), w ^ z * (1 + w ^ 2) ^ (-(3 / 2 : ℝ)) := by
  have hu : ∀ x ∈ Ioi (0:ℝ), HasDerivAt fp (-(1 + x ^ 2) ^ (-(3 / 2 : ℝ))) x :=
    fun x _ => hasDerivAt_fp x
  have hv : ∀ x ∈ Ioi (0:ℝ), HasDerivAt (fun w : ℝ => w ^ z / z) (x ^ (z - 1)) x := by
    intro x hx
    have := (Real.hasDerivAt_rpow_const (p := z) (Or.inl (ne_of_gt hx))).div_const z
    convert this using 1
    field_simp
  have huv' : IntegrableOn (fp * fun x : ℝ => x ^ (z - 1)) (Ioi 0) := int_ii z hz0 hz1
  have hu'v : IntegrableOn ((fun x : ℝ => -(1 + x ^ 2) ^ (-(3 / 2 : ℝ))) *
      fun w : ℝ => w ^ z / z) (Ioi 0) := by
    have : IntegrableOn _ (Ioi (0:ℝ)) := ((K_integral z hz0 hz1).1.div_const z).neg
    refine this.congr_fun (fun x _ => ?_) measurableSet_Ioi
    simp only [Pi.mul_apply, Pi.neg_apply]; ring
  have hbd : ∀ w : ℝ, 0 < w → |fp w * (w ^ z / z)| ≤ w ^ z / z ∧
      |fp w * (w ^ z / z)| ≤ w ^ (z - 2) / z := by
    intro w hw
    obtain ⟨b0, b1, b2⟩ := fp_bounds w hw
    rw [abs_of_nonneg (mul_nonneg b0 (by positivity))]
    constructor
    · calc fp w * (w ^ z / z) ≤ 1 * (w ^ z / z) := mul_le_mul_of_nonneg_right b1 (by positivity)
        _ = w ^ z / z := one_mul _
    · calc fp w * (w ^ z / z) ≤ 1 / w ^ 2 * (w ^ z / z) :=
            mul_le_mul_of_nonneg_right b2 (by positivity)
        _ = (w ^ z * (1 / w ^ 2)) / z := by ring
        _ = w ^ (z - 2) / z := by rw [fp_le_one_div_sq w hw]
  have h0 : Tendsto (fp * fun w : ℝ => w ^ z / z) (𝓝[>] 0) (𝓝 0) := by
    have hlim := (tendsto_rpow_zero_right z hz0).div_const z
    rw [zero_div] at hlim
    refine squeeze_zero_norm' ?_ hlim
    filter_upwards [self_mem_nhdsWithin] with w hw
    rw [Real.norm_eq_abs]; exact (hbd w hw).1
  have hinf : Tendsto (fp * fun w : ℝ => w ^ z / z) atTop (𝓝 0) := by
    have hlim := (tendsto_rpow_neg_atTop (y := 2 - z) (by linarith)).div_const z
    rw [zero_div] at hlim
    refine squeeze_zero_norm' ?_ hlim
    filter_upwards [eventually_gt_atTop 0] with w hw
    rw [Real.norm_eq_abs]
    convert (hbd w hw).2 using 3; ring
  have := integral_Ioi_mul_deriv_eq_deriv_mul hu hv huv' hu'v h0 hinf
  rw [this]
  have e : ∫ x in Ioi (0:ℝ), -(1 + x ^ 2) ^ (-(3 / 2 : ℝ)) * (x ^ z / z)
      = -((1 / z) * ∫ x in Ioi (0:ℝ), x ^ z * (1 + x ^ 2) ^ (-(3 / 2 : ℝ))) := by
    rw [← integral_const_mul, ← integral_neg]; congr 1; funext x; field_simp
  rw [e]
  ring

lemma mellin_main : stmt_4 := by
  intro z hz0 hz1
  have hpt : ∀ w ∈ Ioi (0:ℝ), w ^ (z - 1) * g w = f w * w ^ (z - 2) := by
    intro w hw
    have hw0 : (0:ℝ) < w := hw
    unfold g
    rw [show z - 1 = (z - 2) + 1 by ring, Real.rpow_add_one hw0.ne']
    field_simp
  refine ⟨(int_i z hz0 hz1).congr_fun (fun w hw => (hpt w hw).symm) measurableSet_Ioi, ?_⟩
  unfold Mg
  rw [setIntegral_congr_fun measurableSet_Ioi hpt, ibp1 z hz0 hz1, ibp2 z hz0 hz1,
    (K_integral z hz0 hz1).2]
  have : 1 - z ≠ 0 := by linarith
  field_simp

end X06

theorem X06.proof_4 : X06.stmt_4 := X06.mellin_main

/-! ### Taylor coefficients of `E_g` (`stmt_3`) -/

namespace X06
open Real MeasureTheory Set Filter Topology Finset

/-- `ζ(n)` as a real number. -/
noncomputable def zr (n : ℕ) : ℝ := (riemannZeta n).re

lemma zr_hasSum {k : ℕ} (hk : 1 < k) :
    HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ k) (zr k) := by
  have hs : Summable (fun n : ℕ => 1 / (n : ℝ) ^ k) := Real.summable_one_div_nat_pow.mpr hk
  have hz : zr k = ∑' n : ℕ, 1 / (n : ℝ) ^ k := by
    unfold zr
    rw [zeta_nat_eq_tsum_of_gt_one hk]
    rw [show (fun n : ℕ => 1 / (n : ℂ) ^ k) = fun n : ℕ => ((1 / (n : ℝ) ^ k : ℝ) : ℂ) by
      ext n; push_cast; rfl]
    rw [← Complex.ofReal_tsum, Complex.ofReal_re]
  rw [hz]
  have := hs.hasSum
  rw [← hasSum_nat_add_iff' 1] at this
  simpa [zero_pow (show k ≠ 0 by omega)] using this

lemma lg_prod_shift (x : ℝ) (n : ℕ) :
    ∏ j ∈ range (n + 1), (1 + x + (j : ℝ)) =
      ((n + 1).factorial : ℝ) * ∏ k ∈ range (n + 1), (1 + x / ((k : ℝ) + 1)) := by
  rw [show (((n + 1).factorial : ℕ) : ℝ) = ∏ k ∈ range (n + 1), ((k : ℝ) + 1) by
    rw [← Finset.prod_range_add_one_eq_factorial]; push_cast; rfl, ← Finset.prod_mul_distrib]
  refine Finset.prod_congr rfl fun k _ => ?_
  have : (0 : ℝ) < (k : ℝ) + 1 := by positivity
  field_simp
  ring

lemma lg_fac_pos {x : ℝ} (hx : |x| < 1) (k : ℕ) : 0 < 1 + x / ((k : ℝ) + 1) := by
  have h1 : (1 : ℝ) ≤ (k : ℝ) + 1 := by have := k.cast_nonneg (α := ℝ); linarith
  have : |x / ((k : ℝ) + 1)| < 1 := by
    rw [abs_div, abs_of_pos (by linarith : (0:ℝ) < (k : ℝ) + 1), div_lt_one (by linarith)]
    linarith
  linarith [neg_abs_le (x / ((k : ℝ) + 1))]

/-- the partial products `R_m` -/
noncomputable def lgR (η u v : ℝ) (m : ℕ) : ℝ :=
  ∏ k ∈ range m, ((1 + η / ((k : ℝ) + 1)) / ((1 + u / ((k : ℝ) + 1)) * (1 + v / ((k : ℝ) + 1))))

lemma lg_gammaSeq_ratio {η u v : ℝ} (huv : u + v = η) (hu : |u| < 1) (hv : |v| < 1)
    (hη : |η| < 1) {n : ℕ} (hn : 1 ≤ n) :
    Real.GammaSeq (1 + u) n * Real.GammaSeq (1 + v) n / Real.GammaSeq (1 + η) n =
      (n : ℝ) / ((n : ℝ) + 1) * lgR η u v (n + 1) := by
  unfold Real.GammaSeq lgR
  rw [lg_prod_shift, lg_prod_shift, lg_prod_shift, Finset.prod_div_distrib, Finset.prod_mul_distrib]
  have hn0 : (0 : ℝ) < n := by exact_mod_cast hn
  have pu : 0 < ∏ k ∈ range (n + 1), (1 + u / ((k : ℝ) + 1)) :=
    Finset.prod_pos fun k _ => lg_fac_pos hu k
  have pv : 0 < ∏ k ∈ range (n + 1), (1 + v / ((k : ℝ) + 1)) :=
    Finset.prod_pos fun k _ => lg_fac_pos hv k
  have pη : 0 < ∏ k ∈ range (n + 1), (1 + η / ((k : ℝ) + 1)) :=
    Finset.prod_pos fun k _ => lg_fac_pos hη k
  have e : ∀ x : ℝ, (n : ℝ) ^ (1 + x) = n * (n : ℝ) ^ x := fun x => by
    rw [Real.rpow_add hn0, Real.rpow_one]
  have euv : (n : ℝ) ^ u * (n : ℝ) ^ v = (n : ℝ) ^ η := by rw [← Real.rpow_add hn0, huv]
  rw [e, e, e]
  have hpη : 0 < (n : ℝ) ^ η := Real.rpow_pos_of_pos hn0 _
  have hf : (0 : ℝ) < (n.factorial : ℝ) := by exact_mod_cast n.factorial_pos
  have hf1 : (0 : ℝ) < ((n + 1).factorial : ℝ) := by exact_mod_cast (n + 1).factorial_pos
  rw [Nat.factorial_succ] at hf1 ⊢
  push_cast at hf1 ⊢
  field_simp
  rw [← euv]

lemma lg_log_series_k {η u v : ℝ} (huv : u + v = η) (hu : |u| < 1) (hv : |v| < 1)
    (hη : |η| < 1) (k : ℕ) :
    HasSum (fun j : ℕ => (-1 : ℝ) ^ (j + 2) * (1 / ((k : ℝ) + 1) ^ (j + 2)) / ((j + 2 : ℕ) : ℝ) *
        (u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2)))
      (Real.log (1 + η / ((k : ℝ) + 1)) - Real.log (1 + u / ((k : ℝ) + 1))
        - Real.log (1 + v / ((k : ℝ) + 1))) := by
  have hk : (0 : ℝ) < (k : ℝ) + 1 := by positivity
  have sm : ∀ x : ℝ, |x| < 1 → |-(x / ((k : ℝ) + 1))| < 1 := by
    intro x hx
    have h1 : (1 : ℝ) ≤ (k : ℝ) + 1 := by have := k.cast_nonneg (α := ℝ); linarith
    rw [abs_neg, abs_div, abs_of_pos hk, div_lt_one hk]; linarith
  have Lu := Real.hasSum_pow_div_log_of_abs_lt_one (sm u hu)
  have Lv := Real.hasSum_pow_div_log_of_abs_lt_one (sm v hv)
  have Lη := Real.hasSum_pow_div_log_of_abs_lt_one (sm η hη)
  have L := (Lu.add Lv).sub Lη
  simp only [sub_neg_eq_add] at L
  rw [← hasSum_nat_add_iff' 1] at L
  simp only [Finset.range_one, Finset.sum_singleton] at L
  convert L using 1
  · ext j
    simp only [neg_pow (_ / _), div_pow]
    push_cast
    ring
  · rw [← huv]
    ring_nf

lemma lg_log_R {η u v : ℝ} (hu : |u| < 1) (hv : |v| < 1) (hη : |η| < 1) (m : ℕ) :
    Real.log (lgR η u v m) = ∑ k ∈ range m, (Real.log (1 + η / ((k : ℝ) + 1))
      - Real.log (1 + u / ((k : ℝ) + 1)) - Real.log (1 + v / ((k : ℝ) + 1))) := by
  unfold lgR
  rw [Real.log_prod]
  · refine Finset.sum_congr rfl fun k _ => ?_
    have a := lg_fac_pos hu k; have b := lg_fac_pos hv k; have c := lg_fac_pos hη k
    rw [Real.log_div c.ne' (by positivity), Real.log_mul a.ne' b.ne']
    ring
  · intro k _
    have a := lg_fac_pos hu k; have b := lg_fac_pos hv k; have c := lg_fac_pos hη k
    positivity

lemma lg_R_tendsto {η u v : ℝ} (huv : u + v = η) (hu : |u| < 1) (hv : |v| < 1)
    (hη : |η| < 1) (hη1 : -1 < η) :
    Tendsto (fun n : ℕ => lgR η u v (n + 1)) atTop
      (𝓝 (Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + η))) := by
  have hG : Real.Gamma (1 + η) ≠ 0 := (Real.Gamma_pos_of_pos (by linarith)).ne'
  have h1 := ((Real.GammaSeq_tendsto_Gamma (1 + u)).mul
    (Real.GammaSeq_tendsto_Gamma (1 + v))).div (Real.GammaSeq_tendsto_Gamma (1 + η)) hG
  have h2 : Tendsto (fun n : ℕ => ((n : ℝ) + 1) / n) atTop (𝓝 1) := by
    have := (tendsto_one_div_atTop_nhds_zero_nat).const_add (1 : ℝ)
    rw [add_zero] at this
    refine this.congr' ?_
    filter_upwards [eventually_ge_atTop 1] with n hn
    have : (0 : ℝ) < n := by exact_mod_cast hn
    field_simp
  have h3 := h1.mul h2
  rw [mul_one] at h3
  refine h3.congr' ?_
  filter_upwards [eventually_ge_atTop 1] with n hn
  have : (0 : ℝ) < n := by exact_mod_cast hn
  simp only [Pi.div_apply]
  rw [lg_gammaSeq_ratio huv hu hv hη hn]
  field_simp


/-- `log (Γ(1+u)Γ(1+v)/Γ(1+u+v)) = Σ_{n≥2} (-1)^n ζ(n)/n (u^n + v^n - (u+v)^n)` (adapted from
`X05.proof_3`). -/
lemma log_gamma_series {u v : ℝ} (hu : |u| < 3 / 4) (hv : |v| < 3 / 4) (hη : |u + v| < 3 / 4) :
    HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * zr (k + 2) / ((k + 2 : ℕ) : ℝ) *
        (u ^ (k + 2) + v ^ (k + 2) - (u + v) ^ (k + 2)))
      (Real.log (Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + (u + v)))) := by
  set η := u + v with hηdef
  have huv : u + v = η := rfl
  have hu' : |u| < 1 := by linarith
  have hv' : |v| < 1 := by linarith
  have hη' : |η| < 1 := by linarith
  have hη1 : -1 < η := by cases abs_lt.mp hη; linarith
  have hApos : 0 < Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + η) := by
    have := Real.Gamma_pos_of_pos (show 0 < 1 + u by cases abs_lt.mp hu; linarith)
    have := Real.Gamma_pos_of_pos (show 0 < 1 + v by cases abs_lt.mp hv; linarith)
    have := Real.Gamma_pos_of_pos (show 0 < 1 + η by linarith)
    positivity
  set f : ℕ × ℕ → ℝ := fun p => (-1 : ℝ) ^ (p.2 + 2) * (1 / ((p.1 : ℝ) + 1) ^ (p.2 + 2)) /
      ((p.2 + 2 : ℕ) : ℝ) * (u ^ (p.2 + 2) + v ^ (p.2 + 2) - η ^ (p.2 + 2)) with hfdef
  have hfs : Summable f := by
    have hg1 : Summable (fun k : ℕ => 1 / ((k : ℝ) + 1) ^ 2) := (zr_hasSum (by norm_num)).summable
    have hg2 : Summable (fun j : ℕ => 3 * (3 / 4 : ℝ) ^ j) :=
      (summable_geometric_of_lt_one (by norm_num) (by norm_num)).mul_left 3
    refine Summable.of_norm_bounded (hg1.mul_of_nonneg hg2 (fun k => by positivity)
      (fun j => by positivity)) ?_
    rintro ⟨k, j⟩
    simp only [hfdef, Real.norm_eq_abs]
    have hk1 : (1 : ℝ) ≤ (k : ℝ) + 1 := by have := k.cast_nonneg (α := ℝ); linarith
    have hj : (1 : ℝ) ≤ ((j + 2 : ℕ) : ℝ) := by push_cast; have := j.cast_nonneg (α := ℝ); linarith
    rw [abs_mul, abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
      abs_of_pos (by positivity : (0 : ℝ) < 1 / ((k : ℝ) + 1) ^ (j + 2)),
      abs_of_pos (by positivity : (0 : ℝ) < ((j + 2 : ℕ) : ℝ))]
    have b1 : 1 / ((k : ℝ) + 1) ^ (j + 2) ≤ 1 / ((k : ℝ) + 1) ^ 2 := by
      apply one_div_le_one_div_of_le (by positivity)
      exact pow_le_pow_right₀ hk1 (by omega)
    have b2 : |u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2)| ≤ 3 * (3 / 4 : ℝ) ^ j := by
      have pu : |u| ^ (j + 2) ≤ (3 / 4 : ℝ) ^ j :=
        (pow_le_pow_left₀ (abs_nonneg u) hu.le _).trans
          (pow_le_pow_of_le_one (by norm_num) (by norm_num) (by omega))
      have pv : |v| ^ (j + 2) ≤ (3 / 4 : ℝ) ^ j :=
        (pow_le_pow_left₀ (abs_nonneg v) hv.le _).trans
          (pow_le_pow_of_le_one (by norm_num) (by norm_num) (by omega))
      have pη : |η| ^ (j + 2) ≤ (3 / 4 : ℝ) ^ j :=
        (pow_le_pow_left₀ (abs_nonneg η) hη.le _).trans
          (pow_le_pow_of_le_one (by norm_num) (by norm_num) (by omega))
      calc |u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2)|
          ≤ |u ^ (j + 2)| + |v ^ (j + 2)| + |η ^ (j + 2)| := by
            have := abs_sub (u ^ (j + 2) + v ^ (j + 2)) (η ^ (j + 2))
            have := abs_add_le (u ^ (j + 2)) (v ^ (j + 2))
            linarith
        _ ≤ 3 * (3 / 4 : ℝ) ^ j := by rw [abs_pow, abs_pow, abs_pow]; linarith
    calc 1 / ((k : ℝ) + 1) ^ (j + 2) / ((j + 2 : ℕ) : ℝ) * |u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2)|
        ≤ 1 / ((k : ℝ) + 1) ^ 2 / 1 * (3 * (3 / 4 : ℝ) ^ j) := by
          gcongr
      _ = _ := by simp
  have hf := hfs.hasSum
  have hT := hf.prod_fiberwise (fun k => lg_log_series_k huv hu' hv' hη' k)
  have hf' : HasSum (f ∘ Equiv.prodComm ℕ ℕ) (∑' p, f p) := (Equiv.prodComm ℕ ℕ).hasSum_iff.mpr hf
  have hZ := hf'.prod_fiberwise (fun j => by
    have := ((zr_hasSum (show 1 < j + 2 by omega)).mul_left ((-1 : ℝ) ^ (j + 2))).div_const
      ((j + 2 : ℕ) : ℝ) |>.mul_right (u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2))
    exact this)
  have hlim1 := (tendsto_add_atTop_iff_nat 1).mpr hT.tendsto_sum_nat
  have hlim2 : Tendsto (fun n : ℕ => ∑ k ∈ range (n + 1), (Real.log (1 + η / ((k : ℝ) + 1))
      - Real.log (1 + u / ((k : ℝ) + 1)) - Real.log (1 + v / ((k : ℝ) + 1)))) atTop
      (𝓝 (Real.log (Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + η)))) := by
    have := ((Real.continuousAt_log hApos.ne').tendsto).comp
      (lg_R_tendsto huv hu' hv' hη' hη1)
    refine this.congr fun n => ?_
    simp only [Function.comp]
    exact lg_log_R hu' hv' hη' (n + 1)
  have heq := tendsto_nhds_unique hlim1 hlim2
  rw [heq] at hZ
  convert hZ using 1


lemma zr_two : zr 2 = π ^ 2 / 6 := by
  unfold zr
  rw [show ((2 : ℕ) : ℂ) = 2 by norm_num, riemannZeta_two]
  rw [show (π : ℂ) ^ 2 / 6 = ((π ^ 2 / 6 : ℝ) : ℂ) by push_cast; ring, Complex.ofReal_re]

lemma zr_three : zr 3 = z3 := by unfold zr z3; norm_num

lemma zr_bd {n : ℕ} (hn : 2 ≤ n) : |zr n| ≤ 2 := by
  have h := zr_hasSum (k := n) (by omega)
  have h2 := zr_hasSum (k := 2) (by omega)
  have hnn : 0 ≤ zr n := h.nonneg fun m => by positivity
  have hle : zr n ≤ zr 2 := hasSum_le (fun m => by
    have h1 : (1 : ℝ) ≤ (m : ℝ) + 1 := by have := m.cast_nonneg (α := ℝ); linarith
    exact one_div_le_one_div_of_le (by positivity) (pow_le_pow_right₀ h1 hn)) h h2
  rw [zr_two] at hle
  have := Real.pi_lt_d2
  have : π ^ 2 / 6 ≤ 2 := by nlinarith [Real.pi_pos]
  rw [abs_of_nonneg hnn]; linarith

/-- duplication form `E_g(s) = 2^{-s} Γ(1+s) Γ(1-s/2) / Γ(1+s/2)`. -/
lemma Eg_dup (s : ℝ) (hs : -2 < s) :
    Eg s = (2:ℝ) ^ (-s) * (Real.Gamma (1 + s) * Real.Gamma (1 - s / 2) / Real.Gamma (1 + s / 2)) := by
  have h := Real.Gamma_mul_Gamma_add_half ((1 + s) / 2)
  rw [show (1 + s) / 2 + 1 / 2 = 1 + s / 2 by ring, show 2 * ((1 + s) / 2) = 1 + s by ring,
    show 1 - (1 + s) = -s by ring] at h
  have hG : 0 < Real.Gamma (1 + s / 2) := Real.Gamma_pos_of_pos (by linarith)
  have hpi : 0 < Real.sqrt π := Real.sqrt_pos.2 Real.pi_pos
  have h' : Real.Gamma ((1 + s) / 2) =
      Real.Gamma (1 + s) * (2:ℝ) ^ (-s) * Real.sqrt π / Real.Gamma (1 + s / 2) := by
    rw [eq_div_iff hG.ne', h]
  unfold Eg
  rw [h']
  field_simp

lemma Eg_pos (s : ℝ) (hs : |s| < 1) : 0 < Eg s := by
  obtain ⟨h1, h2⟩ := abs_lt.mp hs
  unfold Eg
  have := Real.Gamma_pos_of_pos (show 0 < (1 + s) / 2 by linarith)
  have := Real.Gamma_pos_of_pos (show 0 < 1 - s / 2 by linarith)
  have : 0 < Real.sqrt π := Real.sqrt_pos.2 Real.pi_pos
  positivity

/-- the Taylor coefficients of `log E_g` at `0`. -/
noncomputable def ce : ℕ → ℝ
  | 0 => 0
  | 1 => -Real.log 2
  | (k + 2) => (-1 : ℝ) ^ (k + 2) * zr (k + 2) / ((k + 2 : ℕ) : ℝ) *
      (1 + (-1 / 2 : ℝ) ^ (k + 2) - (1 / 2 : ℝ) ^ (k + 2))

lemma logEg_hasSum (s : ℝ) (hs : |s| < 1 / 2) :
    HasSum (fun n => ce n * s ^ n) (Real.log (Eg s)) := by
  obtain ⟨hs1, hs2⟩ := abs_lt.mp hs
  have hu : |s| < 3 / 4 := by linarith
  have hv : |-s / 2| < 3 / 4 := by rw [abs_lt]; constructor <;> linarith
  have hη : |s + -s / 2| < 3 / 4 := by rw [abs_lt]; constructor <;> linarith
  have H := log_gamma_series hu hv hη
  rw [show s + -s / 2 = s / 2 by ring, show 1 + -s / 2 = 1 - s / 2 by ring] at H
  rw [← hasSum_nat_add_iff' 2]
  have hG1 := Real.Gamma_pos_of_pos (show 0 < 1 + s by linarith)
  have hG2 := Real.Gamma_pos_of_pos (show 0 < 1 - s / 2 by linarith)
  have hG3 := Real.Gamma_pos_of_pos (show 0 < 1 + s / 2 by linarith)
  have e : Real.log (Eg s) - ∑ i ∈ range 2, ce i * s ^ i =
      Real.log (Real.Gamma (1 + s) * Real.Gamma (1 - s / 2) / Real.Gamma (1 + s / 2)) := by
    rw [Eg_dup s (by linarith), Real.log_mul (by positivity) (by positivity),
      Real.log_rpow two_pos]
    simp [Finset.sum_range_succ, ce]
    ring
  rw [e]
  convert H using 1
  funext k
  simp only [ce]
  ring

/-- the formal power series of `log E_g`. -/
noncomputable def pS : FormalMultilinearSeries ℝ ℝ ℝ := FormalMultilinearSeries.ofScalars ℝ ce

lemma ce_bd (n : ℕ) : |ce n| ≤ 3 := by
  match n with
  | 0 => simp [ce]
  | 1 =>
    simp only [ce, abs_neg]
    rw [abs_of_pos (Real.log_pos (by norm_num))]
    have := Real.log_two_lt_d9; linarith
  | (k + 2) =>
    simp only [ce]
    have hz := zr_bd (n := k + 2) (by omega)
    have hk : (2:ℝ) ≤ ((k + 2 : ℕ) : ℝ) := by push_cast; have := k.cast_nonneg (α := ℝ); linarith
    have hf : |1 + (-1 / 2 : ℝ) ^ (k + 2) - (1 / 2 : ℝ) ^ (k + 2)| ≤ 3 := by
      have a1 : |(-1 / 2 : ℝ) ^ (k + 2)| ≤ 1 := by
        rw [abs_pow]; exact pow_le_one₀ (by positivity) (by norm_num [abs_div])
      have a2 : |(1 / 2 : ℝ) ^ (k + 2)| ≤ 1 := by
        rw [abs_pow]; exact pow_le_one₀ (by positivity) (by norm_num [abs_div])
      have := abs_sub (1 + (-1 / 2 : ℝ) ^ (k + 2)) ((1 / 2 : ℝ) ^ (k + 2))
      have := abs_add_le (1 : ℝ) ((-1 / 2 : ℝ) ^ (k + 2))
      rw [abs_one] at this
      linarith
    rw [abs_mul, abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
      abs_of_pos (by positivity : (0:ℝ) < ((k + 2 : ℕ) : ℝ))]
    have h1 : |zr (k + 2)| / ((k + 2 : ℕ) : ℝ) ≤ 1 := by
      rw [div_le_one (by positivity)]; linarith
    calc |zr (k + 2)| / ((k + 2 : ℕ) : ℝ) * |1 + (-1 / 2 : ℝ) ^ (k + 2) - (1 / 2 : ℝ) ^ (k + 2)|
        ≤ 1 * 3 := mul_le_mul h1 hf (abs_nonneg _) zero_le_one
      _ = 3 := by norm_num

lemma pS_radius : ((1 / 2 : NNReal) : ENNReal) ≤ pS.radius := by
  refine FormalMultilinearSeries.le_radius_of_bound _ 3 (fun n => ?_)
  rw [pS, FormalMultilinearSeries.ofScalars_norm, Real.norm_eq_abs]
  have : ((1 / 2 : NNReal) : ℝ) ^ n ≤ 1 := pow_le_one₀ (by positivity) (by norm_num)
  calc |ce n| * ((1 / 2 : NNReal) : ℝ) ^ n ≤ 3 * 1 :=
        mul_le_mul (ce_bd n) this (by positivity) (by norm_num)
    _ = 3 := by norm_num

lemma mem_eball_half (y : ℝ) (hy : |y| < 1 / 2) : y ∈ Metric.eball (0:ℝ) pS.radius := by
  refine lt_of_lt_of_le ?_ pS_radius
  rw [edist_dist, dist_zero_right, Real.norm_eq_abs]
  have : (1 / 2 : ENNReal) = ENNReal.ofReal (1 / 2) := by
    rw [ENNReal.ofReal_div_of_pos (by norm_num)]; simp
  rw [show ((1 / 2 : NNReal) : ENNReal) = ENNReal.ofReal (1 / 2) by
    rw [ENNReal.ofReal_div_of_pos (by norm_num)]; simp]
  exact (ENNReal.ofReal_lt_ofReal_iff (by norm_num)).2 hy

lemma pS_radius_pos : 0 < pS.radius :=
  lt_of_lt_of_le (by simp) pS_radius

lemma pS_hasFPS : HasFPowerSeriesOnBall pS.sum pS 0 pS.radius := pS.hasFPowerSeriesOnBall pS_radius_pos

lemma S_eq (y : ℝ) (hy : |y| < 1 / 2) : pS.sum y = Real.log (Eg y) := by
  have h := pS.hasSum (mem_eball_half y hy)
  simp only [pS, FormalMultilinearSeries.ofScalars_apply_eq, smul_eq_mul] at h
  exact h.unique (logEg_hasSum y hy)

lemma S_iteratedDeriv (n : ℕ) : iteratedDeriv n pS.sum 0 = (n.factorial : ℝ) * ce n := by
  have h1 := pS_hasFPS.analyticAt.hasFPowerSeriesAt
  have h2 := h1.eq_formalMultilinearSeries pS_hasFPS.hasFPowerSeriesAt
  have h2' : (FormalMultilinearSeries.ofScalars ℝ fun n => iteratedDeriv n pS.sum 0 /
      (n.factorial : ℝ)) = FormalMultilinearSeries.ofScalars ℝ ce := h2
  have h3 := FormalMultilinearSeries.ofScalars_series_injective (E := ℝ) (𝕜 := ℝ) h2'
  have h4 := congrFun h3 n
  rw [div_eq_iff (by positivity)] at h4
  rw [h4]; ring

lemma stmt3_main : stmt_3 := by
  obtain ⟨S, hS⟩ : ∃ S : ℝ → ℝ, S = pS.sum := ⟨_, rfl⟩
  set B := Metric.ball (0:ℝ) (1 / 2) with hBdef
  have hB : IsOpen B := Metric.isOpen_ball
  have h0B : (0:ℝ) ∈ B := Metric.mem_ball_self (by norm_num)
  have memB : ∀ y ∈ B, |y| < 1 / 2 := fun y hy => by
    simpa [hBdef, Metric.mem_ball, Real.dist_eq] using hy
  have hA : AnalyticOnNhd ℝ S B := fun y hy => hS ▸ pS_hasFPS.analyticAt_of_mem (mem_eball_half y (memB y hy))
  have hA1 := hA.deriv
  have hA2 := hA1.deriv
  have hEexp : ∀ y ∈ B, Eg y = Real.exp (S y) := by
    intro y hy
    rw [hS, S_eq y (memB y hy), Real.exp_log (Eg_pos y (by linarith [memB y hy]))]
  have hEd : ∀ y ∈ B, HasDerivAt Eg (Eg y * deriv S y) y := by
    intro y hy
    have hS := (hA y hy).differentiableAt.hasDerivAt
    have := hS.exp
    rw [← hEexp y hy] at this
    exact this.congr_of_eventuallyEq (Filter.eventually_of_mem (hB.mem_nhds hy) hEexp)
  have hS1d : ∀ y ∈ B, HasDerivAt (deriv S) (deriv (deriv S) y) y :=
    fun y hy => (hA1 y hy).differentiableAt.hasDerivAt
  have hS2d : ∀ y ∈ B, HasDerivAt (deriv (deriv S)) (deriv (deriv (deriv S)) y) y :=
    fun y hy => (hA2 y hy).differentiableAt.hasDerivAt
  have hD1 : ∀ y ∈ B, deriv Eg y = Eg y * deriv S y := fun y hy => (hEd y hy).deriv
  have hD2 : ∀ y ∈ B, deriv (deriv Eg) y =
      Eg y * deriv S y * deriv S y + Eg y * deriv (deriv S) y := by
    intro y hy
    have ev : deriv Eg =ᶠ[𝓝 y] fun x => Eg x * deriv S x :=
      Filter.eventually_of_mem (hB.mem_nhds hy) hD1
    rw [ev.deriv_eq]
    exact ((hEd y hy).mul (hS1d y hy)).deriv
  have hD3 : deriv (deriv (deriv Eg)) 0 =
      (Eg 0 * deriv S 0 * deriv S 0 + Eg 0 * deriv (deriv S) 0) * deriv S 0
        + Eg 0 * deriv S 0 * deriv (deriv S) 0
        + ((Eg 0 * deriv S 0) * deriv (deriv S) 0 + Eg 0 * deriv (deriv (deriv S)) 0) := by
    have ev : deriv (deriv Eg) =ᶠ[𝓝 0] fun x =>
        Eg x * deriv S x * deriv S x + Eg x * deriv (deriv S) x :=
      Filter.eventually_of_mem (hB.mem_nhds h0B) hD2
    rw [ev.deriv_eq]
    exact ((((hEd 0 h0B).mul (hS1d 0 h0B)).mul (hS1d 0 h0B)).add
      ((hEd 0 h0B).mul (hS2d 0 h0B))).deriv
  have hE0 : Eg 0 = 1 := by
    unfold Eg
    rw [show (1 + (0:ℝ)) / 2 = 1 / 2 by norm_num, show 1 - (0:ℝ) / 2 = 1 by norm_num,
      Real.Gamma_one_half_eq, Real.Gamma_one, mul_one, div_self (Real.sqrt_pos.2 Real.pi_pos).ne']
  have c1 : deriv S 0 = -lam := by
    have := S_iteratedDeriv 1
    rw [← hS] at this
    rw [iteratedDeriv_one] at this
    rw [this]; simp [ce, lam]
  have c2 : deriv (deriv S) 0 = π ^ 2 / 6 := by
    have := S_iteratedDeriv 2
    rw [← hS] at this
    rw [iteratedDeriv_succ, iteratedDeriv_one] at this
    rw [this]; simp [ce, zr_two]; norm_num; ring
  have c3 : deriv (deriv (deriv S)) 0 = -(3 / 2) * z3 := by
    have := S_iteratedDeriv 3
    rw [← hS] at this
    rw [iteratedDeriv_succ, iteratedDeriv_succ, iteratedDeriv_one] at this
    rw [this]; simp [ce, zr_three]; norm_num; ring
  refine ⟨?_, ?_, ?_, ?_⟩
  · simp [eta, hE0]
  · simp only [eta, iteratedDeriv_one, Nat.factorial_one, Nat.cast_one, div_one]
    rw [hD1 0 h0B, hE0, c1, one_mul]
  · simp only [eta]
    rw [iteratedDeriv_succ, iteratedDeriv_one, hD2 0 h0B, hE0, c1, c2]
    simp [Nat.factorial]; ring
  · simp only [eta]
    rw [iteratedDeriv_succ, iteratedDeriv_succ, iteratedDeriv_one, hD3, hE0, c1, c2, c3]
    simp [Nat.factorial]; ring

end X06

theorem X06.proof_3 : X06.stmt_3 := X06.stmt3_main

namespace X06
open Real MeasureTheory Set Filter Topology

lemma Eg_analyticAt0 : AnalyticAt ℝ Eg 0 := by
  have h : AnalyticAt ℝ (fun z => Real.exp (pS.sum z)) 0 := pS_hasFPS.analyticAt.rexp'
  refine h.congr ?_
  filter_upwards [Metric.ball_mem_nhds (0:ℝ) (show (0:ℝ) < 1 / 2 by norm_num)] with y hy
  have hy' : |y| < 1 / 2 := by simpa [Metric.mem_ball, Real.dist_eq] using hy
  rw [S_eq y hy', Real.exp_log (Eg_pos y (by linarith))]

/-! ### The moment identity (`stmt_5`) -/

/-- `h = g - 1_{w>1}/w`. -/
noncomputable def hh (w : ℝ) : ℝ := g w - (Ioi (1:ℝ)).indicator (fun w => 1 / w) w

/-- `r(w) = 1/w - g(w) = 1/(w(w+√(1+w²)))`. -/
noncomputable def rr (w : ℝ) : ℝ := 1 / (w * (w + Real.sqrt (1 + w ^ 2)))

lemma hh_small (w : ℝ) (hw1 : w ≤ 1) : hh w = g w := by
  unfold hh; rw [indicator_of_notMem (by simpa using hw1), sub_zero]

lemma hh_large (w : ℝ) (hw : 1 < w) : hh w = -rr w := by
  unfold hh rr g f
  rw [indicator_of_mem (show w ∈ Ioi (1:ℝ) from hw)]
  have hw0 : 0 < w := by linarith
  set s := Real.sqrt (1 + w ^ 2)
  have hs : 0 < s := Real.sqrt_pos.2 (by positivity)
  have hss : s ^ 2 = 1 + w ^ 2 := Real.sq_sqrt (by positivity)
  field_simp
  nlinarith

lemma g_bounds (w : ℝ) (hw : 0 < w) : 0 ≤ g w ∧ g w ≤ 1 := by
  unfold g
  exact ⟨div_nonneg (f_nonneg w hw.le) hw.le, (div_le_one hw).2 (f_le_self w)⟩

lemma rr_bounds (w : ℝ) (hw : 0 < w) :
    0 ≤ 1 / (2 * w ^ 2) - rr w ∧ 1 / (2 * w ^ 2) - rr w ≤ 1 / (8 * w ^ 4) := by
  unfold rr
  set s := Real.sqrt (1 + w ^ 2)
  have hs : 0 < s := Real.sqrt_pos.2 (by positivity)
  have hss : s ^ 2 = 1 + w ^ 2 := Real.sq_sqrt (by positivity)
  have hws : w < s := by nlinarith
  have e : 1 / (2 * w ^ 2) - 1 / (w * (w + s)) = 1 / (2 * w ^ 2 * (w + s) ^ 2) := by
    field_simp; nlinarith
  rw [e]
  refine ⟨by positivity, ?_⟩
  apply one_div_le_one_div_of_le (by positivity)
  have : (2 * w) ^ 2 ≤ (w + s) ^ 2 := by nlinarith
  nlinarith

lemma rr_nonneg (w : ℝ) (hw : 0 < w) : 0 ≤ rr w := by
  unfold rr; have : 0 ≤ Real.sqrt (1 + w ^ 2) := Real.sqrt_nonneg _; positivity

lemma hh_abs_small (w : ℝ) (hw0 : 0 < w) (hw1 : w ≤ 1) : |hh w| ≤ 1 := by
  rw [hh_small w hw1, abs_of_nonneg (g_bounds w hw0).1]; exact (g_bounds w hw0).2

lemma hh_abs_large (w : ℝ) (hw : 1 < w) : |hh w| ≤ 1 / (2 * w ^ 2) := by
  have hw0 : 0 < w := by linarith
  rw [hh_large w hw, abs_neg, abs_of_nonneg (rr_nonneg w hw0)]
  linarith [(rr_bounds w hw0).1]

lemma g_measurable : Measurable g := by
  unfold g; exact f_continuous.measurable.div measurable_id

lemma hh_measurable : Measurable hh := by
  unfold hh
  exact g_measurable.sub ((measurable_const.div measurable_id).indicator measurableSet_Ioi)

/-- integrability on `(0,∞)` from power bounds near `0` and `∞` (measurable version). -/
lemma int_pieces_m (h : ℝ → ℝ) (hm : Measurable h) {p r : ℝ} (hp : -1 < p) (hr : r < -1)
    (h1 : ∀ w ∈ Ioc (0:ℝ) 1, |h w| ≤ w ^ p) (h2 : ∀ w ∈ Ioi (1:ℝ), |h w| ≤ w ^ r) :
    IntegrableOn h (Ioi 0) := by
  rw [← Ioc_union_Ioi_eq_Ioi zero_le_one]
  refine IntegrableOn.union ?_ ?_
  · have hi : IntegrableOn (fun x : ℝ => x ^ p) (Ioc 0 1) :=
      (intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one).1
        (intervalIntegral.intervalIntegrable_rpow' hp)
    refine hi.mono' hm.aestronglyMeasurable ?_
    filter_upwards [ae_restrict_mem measurableSet_Ioc] with w hw
    rw [Real.norm_eq_abs]; exact h1 w hw
  · have hi : IntegrableOn (fun x : ℝ => x ^ r) (Ioi 1) := integrableOn_Ioi_rpow_of_lt hr one_pos
    refine hi.mono' hm.aestronglyMeasurable ?_
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with w hw
    rw [Real.norm_eq_abs]; exact h2 w hw

/-- the master majorant `e^{|log w|/2} |h(w)|`. -/
noncomputable def MM (w : ℝ) : ℝ := Real.exp (|Real.log w| / 2) * |hh w|

lemma MM_measurable : Measurable MM := by
  unfold MM
  exact (Real.measurable_exp.comp (Real.measurable_log.abs.div_const 2)).mul hh_measurable.abs

lemma MM_int : IntegrableOn MM (Ioi 0) := by
  refine int_pieces_m MM MM_measurable (p := -(1 / 2)) (r := -(3 / 2)) (by norm_num) (by norm_num)
    ?_ ?_
  · intro w hw
    have hw0 : 0 < w := hw.1
    have hl : Real.log w ≤ 0 := Real.log_nonpos hw0.le hw.2
    unfold MM
    rw [abs_of_nonneg (by positivity), abs_of_nonpos hl, Real.rpow_def_of_pos hw0]
    calc Real.exp (-Real.log w / 2) * |hh w| ≤ Real.exp (-Real.log w / 2) * 1 :=
          mul_le_mul_of_nonneg_left (hh_abs_small w hw0 hw.2) (by positivity)
      _ = Real.exp (Real.log w * -(1 / 2)) := by ring_nf
  · intro w hw
    have hw1 : (1:ℝ) < w := hw
    have hw0 : 0 < w := by linarith
    have hl : 0 ≤ Real.log w := Real.log_nonneg hw1.le
    unfold MM
    rw [abs_of_nonneg (by positivity), abs_of_nonneg hl]
    have e1 : Real.exp (Real.log w / 2) = w ^ (1 / 2 : ℝ) := by
      rw [Real.rpow_def_of_pos hw0]; ring_nf
    have e2 : w ^ (1 / 2 : ℝ) * (1 / (2 * w ^ 2)) = w ^ (-(3 / 2) : ℝ) / 2 := by
      rw [show (-(3 / 2) : ℝ) = 1 / 2 - 2 by norm_num, Real.rpow_sub hw0,
        show (2:ℝ) = ((2:ℕ):ℝ) by norm_num, Real.rpow_natCast]
      field_simp
    rw [e1]
    calc w ^ (1 / 2 : ℝ) * |hh w| ≤ w ^ (1 / 2 : ℝ) * (1 / (2 * w ^ 2)) :=
          mul_le_mul_of_nonneg_left (hh_abs_large w hw1) (by positivity)
      _ = w ^ (-(3 / 2) : ℝ) / 2 := e2
      _ ≤ w ^ (-(3 / 2) : ℝ) := by
          have : 0 ≤ w ^ (-(3 / 2) : ℝ) := by positivity
          linarith

lemma pow_log_le_MM (i : ℕ) (w : ℝ) :
    |Real.log w| ^ i / (i.factorial : ℝ) * |hh w| ≤ 2 ^ i * MM w := by
  unfold MM
  have h := Real.pow_div_factorial_le_exp (x := |Real.log w| / 2) (by positivity) i
  have e : |Real.log w| ^ i / (i.factorial : ℝ) = 2 ^ i * ((|Real.log w| / 2) ^ i / (i.factorial : ℝ)) := by
    rw [div_pow]; field_simp
  rw [e, mul_assoc]
  exact mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_right h (abs_nonneg _)) (by positivity)

lemma mom_int (i : ℕ) : IntegrableOn (fun w : ℝ => (-Real.log w) ^ i * hh w) (Ioi 0) := by
  refine ((MM_int.const_mul (2 ^ i * (i.factorial : ℝ)))).mono'
    ((Real.measurable_log.neg.pow_const i).mul hh_measurable).aestronglyMeasurable ?_
  refine Eventually.of_forall (fun w => ?_)
  rw [Real.norm_eq_abs, abs_mul, abs_pow, abs_neg]
  have h := pow_log_le_MM i w
  have hf : (0:ℝ) < i.factorial := by exact_mod_cast i.factorial_pos
  rw [div_mul_eq_mul_div, div_le_iff₀ hf] at h
  linarith


/-- `b_i = (1/i!) ∫_0^∞ (-log w)^i h(w) dw`. -/
noncomputable def bm (i : ℕ) : ℝ := ∫ w in Ioi (0:ℝ), (-Real.log w) ^ i / (i.factorial : ℝ) * hh w

/-- `H(s) = ∫_0^∞ w^{-s} h(w) dw`. -/
noncomputable def HH (s : ℝ) : ℝ := ∫ w in Ioi (0:ℝ), w ^ (-s) * hh w

lemma exp_hasSum (x : ℝ) : HasSum (fun n => x ^ n / (n.factorial : ℝ)) (Real.exp x) := by
  rw [Real.exp_eq_exp_ℝ]; exact NormedSpace.expSeries_div_hasSum_exp x

lemma HH_hasSum (s : ℝ) (hs : |s| ≤ 1 / 2) : HasSum (fun i => bm i * s ^ i) (HH s) := by
  have key := hasSum_integral_of_dominated_convergence (μ := volume.restrict (Ioi (0:ℝ)))
    (F := fun i w => (-Real.log w) ^ i / (i.factorial : ℝ) * hh w * s ^ i)
    (f := fun w => w ^ (-s) * hh w)
    (fun i w => (|Real.log w| / 2) ^ i / (i.factorial : ℝ) * |hh w|) ?_ ?_ ?_ ?_ ?_
  · convert key using 1
    · funext i
      unfold bm
      rw [integral_mul_const]
    · rfl
  · intro i
    exact (((Real.measurable_log.neg.pow_const i).div_const _).mul hh_measurable).mul_const _
      |>.aestronglyMeasurable
  · intro i
    refine Eventually.of_forall (fun w => ?_)
    rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_div, abs_pow, abs_pow, abs_neg,
      abs_of_pos (by exact_mod_cast i.factorial_pos : (0:ℝ) < i.factorial)]
    have h1 : |s| ^ i ≤ (1 / 2) ^ i := pow_le_pow_left₀ (abs_nonneg s) hs i
    have e : (|Real.log w| / 2) ^ i / (i.factorial : ℝ) * |hh w| =
        |Real.log w| ^ i / (i.factorial : ℝ) * |hh w| * (1 / 2) ^ i := by
      rw [div_pow, div_pow, one_pow]; ring
    rw [e]
    exact mul_le_mul_of_nonneg_left h1 (by positivity)
  · refine Eventually.of_forall (fun w => ?_)
    exact ((exp_hasSum (|Real.log w| / 2)).mul_right |hh w|).summable
  · have e : (fun w => ∑' n, (|Real.log w| / 2) ^ n / (n.factorial : ℝ) * |hh w|) = MM := by
      funext w
      unfold MM
      exact ((exp_hasSum (|Real.log w| / 2)).mul_right |hh w|).tsum_eq
    rw [e]; exact MM_int
  · filter_upwards [ae_restrict_mem measurableSet_Ioi] with w hw
    have hw0 : (0:ℝ) < w := hw
    have h := (exp_hasSum (-s * Real.log w)).mul_right (hh w)
    rw [show Real.exp (-s * Real.log w) = w ^ (-s) by
      rw [Real.rpow_def_of_pos hw0]; ring_nf] at h
    convert h using 1
    funext i
    rw [mul_pow]; ring


lemma Eg_reflect (s : ℝ) : Eg (1 - s) = Eg s := by
  unfold Eg
  rw [show (1 + (1 - s)) / 2 = 1 - s / 2 by ring, show 1 - (1 - s) / 2 = (1 + s) / 2 by ring,
    mul_comm (Real.Gamma (1 - s / 2))]

lemma ind_int (s : ℝ) (hs : 0 < s) :
    IntegrableOn (fun w : ℝ => w ^ (-s) * (Ioi (1:ℝ)).indicator (fun w => 1 / w) w) (Ioi 0) ∧
    ∫ w in Ioi (0:ℝ), w ^ (-s) * (Ioi (1:ℝ)).indicator (fun w => 1 / w) w = 1 / s := by
  have e : (fun w : ℝ => w ^ (-s) * (Ioi (1:ℝ)).indicator (fun w => 1 / w) w) =
      (Ioi (1:ℝ)).indicator (fun w => w ^ (-s - 1)) := by
    funext w
    by_cases hw : w ∈ Ioi (1:ℝ)
    · have hw0 : (0:ℝ) < w := lt_trans one_pos hw
      rw [indicator_of_mem hw, indicator_of_mem hw, Real.rpow_sub_one hw0.ne']; ring
    · rw [indicator_of_notMem hw, indicator_of_notMem hw, mul_zero]
  have hI : IntegrableOn (fun w : ℝ => w ^ (-s - 1)) (Ioi 1) :=
    integrableOn_Ioi_rpow_of_lt (by linarith) one_pos
  rw [e]
  refine ⟨(hI.integrable_indicator measurableSet_Ioi).integrableOn, ?_⟩
  rw [setIntegral_indicator measurableSet_Ioi,
    inter_eq_right.2 (Ioi_subset_Ioi zero_le_one), integral_Ioi_rpow_of_lt (by linarith) one_pos]
  rw [show -s - 1 + 1 = -s by ring, Real.one_rpow]
  field_simp

lemma HH_mellin (s : ℝ) (hs0 : 0 < s) (hs1 : s < 1) : HH s = Eg s / ((1 - s) * s) - 1 / s := by
  obtain ⟨hi, hM⟩ := mellin_main (1 - s) (by linarith) (by linarith)
  have e1 : (fun w : ℝ => w ^ ((1 - s) - 1) * g w) = fun w => w ^ (-s) * g w := by
    funext w; ring_nf
  rw [e1] at hi
  unfold Mg at hM
  rw [e1, Eg_reflect, show 1 - (1 - s) = s by ring] at hM
  obtain ⟨hi2, hI2⟩ := ind_int s hs0
  unfold HH hh
  simp_rw [mul_sub]
  rw [integral_sub hi hi2, hM, hI2]

/-- coefficients of `(1-s)(1+sH(s))`. -/
noncomputable def dd : ℕ → ℝ
  | 0 => 1
  | 1 => -1 + bm 0
  | (n + 2) => bm (n + 1) - bm n

lemma MM_nonneg (w : ℝ) : 0 ≤ MM w := by unfold MM; positivity

lemma bm_bd (i : ℕ) : |bm i| ≤ 2 ^ i * ∫ w in Ioi (0:ℝ), MM w := by
  unfold bm
  rw [← integral_const_mul, ← Real.norm_eq_abs]
  refine norm_integral_le_of_norm_le (MM_int.const_mul _) (Eventually.of_forall (fun w => ?_))
  rw [Real.norm_eq_abs, abs_mul, abs_div, abs_pow, abs_neg,
    abs_of_pos (by exact_mod_cast i.factorial_pos : (0:ℝ) < i.factorial)]
  exact pow_log_le_MM i w

lemma dd_bd (n : ℕ) : |dd n| * (1 / 4) ^ n ≤ 2 * (1 + ∫ w in Ioi (0:ℝ), MM w) := by
  set C := ∫ w in Ioi (0:ℝ), MM w
  have hC : 0 ≤ C := setIntegral_nonneg measurableSet_Ioi (fun w _ => MM_nonneg w)
  match n with
  | 0 => simp [dd]; linarith
  | 1 =>
    simp only [dd, pow_one]
    have := bm_bd 0
    simp only [pow_zero, one_mul] at this
    have h2 : |-1 + bm 0| ≤ 1 + C := by
      have := abs_add_le (-1 : ℝ) (bm 0); rw [abs_neg, abs_one] at this; linarith
    nlinarith [abs_nonneg (-1 + bm 0)]
  | (k + 2) =>
    simp only [dd]
    have h1 := bm_bd (k + 1)
    have h0 := bm_bd k
    have h2 : |bm (k + 1) - bm k| ≤ 2 ^ (k + 1) * C + 2 ^ k * C := by
      have := abs_sub (bm (k + 1)) (bm k); linarith
    have e : (2:ℝ) ^ (k + 1) * C + 2 ^ k * C ≤ 2 ^ (k + 2) * C := by
      have : (0:ℝ) ≤ 2 ^ k * C := by positivity
      rw [show (2:ℝ) ^ (k + 1) = 2 * 2 ^ k by ring, show (2:ℝ) ^ (k + 2) = 4 * 2 ^ k by ring]
      nlinarith
    have hq : (2:ℝ) ^ (k + 2) * (1 / 4) ^ (k + 2) ≤ 1 := by
      rw [← mul_pow]; exact pow_le_one₀ (by norm_num) (by norm_num)
    calc |bm (k + 1) - bm k| * (1 / 4) ^ (k + 2) ≤ 2 ^ (k + 2) * C * (1 / 4) ^ (k + 2) :=
          mul_le_mul_of_nonneg_right (h2.trans e) (by positivity)
      _ = C * (2 ^ (k + 2) * (1 / 4) ^ (k + 2)) := by ring
      _ ≤ C * 1 := mul_le_mul_of_nonneg_left hq hC
      _ ≤ 2 * (1 + C) := by linarith

noncomputable def pD : FormalMultilinearSeries ℝ ℝ ℝ := FormalMultilinearSeries.ofScalars ℝ dd

lemma pD_radius : ((1 / 4 : NNReal) : ENNReal) ≤ pD.radius := by
  refine FormalMultilinearSeries.le_radius_of_bound _ (2 * (1 + ∫ w in Ioi (0:ℝ), MM w))
    (fun n => ?_)
  rw [pD, FormalMultilinearSeries.ofScalars_norm, Real.norm_eq_abs]
  have := dd_bd n
  convert this using 2
  push_cast; ring

lemma mem_eball_quarter (y : ℝ) (hy : |y| < 1 / 4) : y ∈ Metric.eball (0:ℝ) pD.radius := by
  refine lt_of_lt_of_le ?_ pD_radius
  rw [edist_dist, dist_zero_right, Real.norm_eq_abs]
  rw [show ((1 / 4 : NNReal) : ENNReal) = ENNReal.ofReal (1 / 4) by
    rw [ENNReal.ofReal_div_of_pos (by norm_num)]; simp]
  exact (ENNReal.ofReal_lt_ofReal_iff (by norm_num)).2 hy

lemma pD_hasFPS : HasFPowerSeriesOnBall pD.sum pD 0 pD.radius :=
  pD.hasFPowerSeriesOnBall (lt_of_lt_of_le (by simp) pD_radius)

lemma Eg_series (s : ℝ) (hs0 : 0 < s) (hs1 : s < 1 / 4) :
    HasSum (fun n => dd n * s ^ n) (Eg s) := by
  have hH := HH_hasSum s (by rw [abs_of_pos hs0]; linarith)
  have hM := HH_mellin s hs0 (by linarith)
  have hE : Eg s = (1 - s) * (1 + s * HH s) := by
    rw [hM]; have : 1 - s ≠ 0 := by linarith
    field_simp; ring
  have A1 := hH.mul_left s
  rw [← hasSum_nat_add_iff' 1] at A1
  have A2 := hH.mul_left (s ^ 2)
  have A := A1.sub A2
  rw [← hasSum_nat_add_iff' 2]
  convert A using 1
  · funext n; simp only [dd]; ring
  · rw [hE]; simp [Finset.sum_range_succ, dd]; ring

lemma Eg_eq_D (s : ℝ) (hs0 : 0 < s) (hs1 : s < 1 / 4) : Eg s = pD.sum s := by
  have h := pD.hasSum (mem_eball_quarter s (by rw [abs_of_pos hs0]; exact hs1))
  simp only [pD, FormalMultilinearSeries.ofScalars_apply_eq, smul_eq_mul] at h
  exact (Eg_series s hs0 hs1).unique h

lemma eta_eq_dd (n : ℕ) : eta n = dd n := by
  have hfreq : ∃ᶠ z in 𝓝[≠] (0:ℝ), Eg z = pD.sum z := by
    have : ∀ᶠ z in 𝓝[>] (0:ℝ), Eg z = pD.sum z := by
      filter_upwards [Ioo_mem_nhdsGT (show (0:ℝ) < 1 / 4 by norm_num)] with z hz
      exact Eg_eq_D z hz.1 hz.2
    exact this.frequently.filter_mono (nhdsWithin_mono _ (fun x hx => ne_of_gt hx))
  have hev := (Eg_analyticAt0.frequently_eq_iff_eventually_eq pD_hasFPS.analyticAt).1 hfreq
  have hev' : Eg =ᶠ[𝓝 0] pD.sum := hev
  have hP : HasFPowerSeriesAt Eg pD 0 := pD_hasFPS.hasFPowerSeriesAt.congr hev'.symm
  have h1 := Eg_analyticAt0.hasFPowerSeriesAt
  have h2 := h1.eq_formalMultilinearSeries hP
  have h2' : (FormalMultilinearSeries.ofScalars ℝ fun n => iteratedDeriv n Eg 0 /
      (n.factorial : ℝ)) = FormalMultilinearSeries.ofScalars ℝ dd := h2
  have h3 := FormalMultilinearSeries.ofScalars_series_injective (E := ℝ) (𝕜 := ℝ) h2'
  exact congrFun h3 n

lemma mt_eq (i : ℕ) : mt (i + 1) = bm i := by
  induction i with
  | zero =>
    unfold mt
    simp [Finset.sum_range_succ, eta_eq_dd, dd]
  | succ i ih =>
    unfold mt at ih ⊢
    rw [Finset.sum_range_succ, ih, eta_eq_dd]
    simp only [dd]; ring


/-- substitution `w = N e^t` from `(0,∞)` onto `(N,∞)`. -/
lemma subst_Nexp (N : ℝ) (hN : 0 < N) (G : ℝ → ℝ) :
    (IntegrableOn G (Ioi N) ↔
      IntegrableOn (fun t => N * Real.exp t * G (N * Real.exp t)) (Ioi 0)) ∧
      ∫ w in Ioi N, G w = ∫ t in Ioi (0:ℝ), N * Real.exp t * G (N * Real.exp t) := by
  have himg : (fun t : ℝ => N * Real.exp t) '' Ioi 0 = Ioi N := by
    ext w; constructor
    · rintro ⟨t, ht, rfl⟩
      have : 1 < Real.exp t := Real.one_lt_exp_iff.2 ht
      show N < N * Real.exp t
      nlinarith
    · intro hw
      have hw' : N < w := hw
      refine ⟨Real.log (w / N), Real.log_pos ((one_lt_div hN).2 hw'), ?_⟩
      simp only
      rw [Real.exp_log (div_pos (by linarith) hN)]; field_simp
  have hd : ∀ x ∈ Ioi (0:ℝ), HasDerivWithinAt (fun t : ℝ => N * Real.exp t)
      (N * Real.exp x) (Ioi 0) x := fun x _ => ((Real.hasDerivAt_exp x).const_mul N).hasDerivWithinAt
  have hinj : InjOn (fun t : ℝ => N * Real.exp t) (Ioi 0) := by
    intro x _ y _ h
    exact Real.exp_injective (mul_left_cancel₀ hN.ne' h)
  have habs : ∀ x : ℝ, |N * Real.exp x| = N * Real.exp x := fun x => abs_of_pos (by positivity)
  constructor
  · rw [← himg, integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ioi hd hinj]
    simp only [habs, smul_eq_mul]
  · rw [← himg, integral_image_eq_integral_abs_deriv_smul measurableSet_Ioi hd hinj]
    simp only [habs, smul_eq_mul]

lemma moment3 (k : ℕ) :
    IntegrableOn (fun t : ℝ => t ^ k * Real.exp (-(3 * t))) (Ioi 0) ∧
      ∫ t in Ioi (0:ℝ), t ^ k * Real.exp (-(3 * t)) = (k.factorial : ℝ) / 3 ^ (k + 1) := by
  have hk : (0 : ℝ) < ↑k + 1 := by positivity
  have key := integral_rpow_mul_exp_neg_mul_Ioi hk (show (0:ℝ) < 3 by norm_num)
  have e : ∀ t ∈ Ioi (0:ℝ), t ^ ((↑k + 1 : ℝ) - 1) * Real.exp (-(3 * t)) =
      t ^ k * Real.exp (-(3 * t)) := by
    intro t ht; simp [add_sub_cancel_right, Real.rpow_natCast]
  rw [setIntegral_congr_fun measurableSet_Ioi e, Real.Gamma_nat_eq_factorial,
    show (↑k + 1 : ℝ) = ((k + 1 : ℕ) : ℝ) by push_cast; ring, Real.rpow_natCast] at key
  rw [key, div_pow, one_pow]
  refine ⟨?_, by ring⟩
  refine ((gammaMoment k).1).mono' ?_ ?_
  · exact ((continuous_pow k).mul (by fun_prop)).aestronglyMeasurable
  · filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
    have ht0 : (0:ℝ) < t := ht
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
    exact mul_le_mul_of_nonneg_left (Real.exp_le_exp.2 (by linarith)) (by positivity)

lemma rr_continuous_exp (N : ℝ) (hN : 0 < N) : Continuous (fun t : ℝ => rr (N * Real.exp t)) := by
  unfold rr
  have hc : Continuous (fun t : ℝ => N * Real.exp t) := continuous_const.mul Real.continuous_exp
  refine continuous_const.div (hc.mul (hc.add (Real.continuous_sqrt.comp
    (continuous_const.add (hc.pow 2))))) (fun t => ?_)
  have hp : 0 < N * Real.exp t := by positivity
  exact (mul_pos hp (add_pos_of_pos_of_nonneg hp (Real.sqrt_nonneg _))).ne'

lemma tail_rr (k : ℕ) (N : ℝ) (hN : 1 ≤ N) :
    IntegrableOn (fun w => Real.log (w / N) ^ k * rr w) (Ioi N) ∧
    ∃ J : ℝ, 0 ≤ J ∧ J ≤ (k.factorial : ℝ) / (8 * 3 ^ (k + 1) * N ^ 3) ∧
      ∫ w in Ioi N, Real.log (w / N) ^ k * rr w = (k.factorial : ℝ) / (2 * N) - J := by
  have hN0 : 0 < N := by linarith
  obtain ⟨hiff, heq⟩ := subst_Nexp N hN0 (fun w => Real.log (w / N) ^ k * rr w)
  set ρ : ℝ → ℝ := fun t => Real.exp (-t) / (2 * N) - N * Real.exp t * rr (N * Real.exp t)
  have hpt : ∀ t : ℝ, N * Real.exp t * (Real.log (N * Real.exp t / N) ^ k * rr (N * Real.exp t))
      = t ^ k * Real.exp (-t) / (2 * N) - t ^ k * ρ t := by
    intro t
    rw [show N * Real.exp t / N = Real.exp t by field_simp, Real.log_exp]
    simp only [ρ]; ring
  have hρ : ∀ t : ℝ, 0 ≤ ρ t ∧ ρ t ≤ Real.exp (-(3 * t)) / (8 * N ^ 3) := by
    intro t
    set w := N * Real.exp t with hw
    have hw0 : 0 < w := by positivity
    obtain ⟨b1, b2⟩ := rr_bounds w hw0
    have e1 : Real.exp (-t) / (2 * N) = w * (1 / (2 * w ^ 2)) := by
      rw [hw, Real.exp_neg]; field_simp
    have e2 : Real.exp (-(3 * t)) / (8 * N ^ 3) = w * (1 / (8 * w ^ 4)) := by
      rw [hw, show -(3 * t) = -t + -t + -t by ring, Real.exp_add, Real.exp_add, Real.exp_neg]
      field_simp
    simp only [ρ]
    rw [e1, e2]
    constructor
    · nlinarith
    · nlinarith
  have hρc : Continuous ρ := by
    simp only [ρ]
    exact (by fun_prop : Continuous fun t : ℝ => Real.exp (-t) / (2 * N)).sub
      ((by fun_prop : Continuous fun t : ℝ => N * Real.exp t).mul (rr_continuous_exp N hN0))
  have hρi : IntegrableOn (fun t => t ^ k * ρ t) (Ioi 0) := by
    refine (((moment3 k).1).div_const (8 * N ^ 3)).mono'
      ((continuous_pow k).mul hρc).aestronglyMeasurable ?_
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
    have ht0 : (0:ℝ) < t := ht
    rw [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg (by positivity) (hρ t).1), mul_div_assoc]
    exact mul_le_mul_of_nonneg_left (hρ t).2 (by positivity)
  have hm1 : IntegrableOn (fun t : ℝ => t ^ k * Real.exp (-t) / (2 * N)) (Ioi 0) :=
    ((gammaMoment k).1).div_const _
  have hint : IntegrableOn (fun t => N * Real.exp t * (Real.log (N * Real.exp t / N) ^ k *
      rr (N * Real.exp t))) (Ioi 0) :=
    (hm1.sub hρi).congr_fun (fun t _ => (hpt t).symm) measurableSet_Ioi
  refine ⟨hiff.2 hint, ∫ t in Ioi (0:ℝ), t ^ k * ρ t,
    setIntegral_nonneg measurableSet_Ioi (fun t ht => mul_nonneg (pow_nonneg (le_of_lt ht) k) (hρ t).1),
    ?_, ?_⟩
  · calc ∫ t in Ioi (0:ℝ), t ^ k * ρ t ≤ ∫ t in Ioi (0:ℝ), t ^ k * Real.exp (-(3 * t)) / (8 * N ^ 3) :=
          setIntegral_mono_on hρi (((moment3 k).1).div_const _) measurableSet_Ioi (fun t ht => by
            rw [mul_div_assoc]
            exact mul_le_mul_of_nonneg_left (hρ t).2 (pow_nonneg (le_of_lt ht) k))
      _ = (k.factorial : ℝ) / (8 * 3 ^ (k + 1) * N ^ 3) := by
          rw [integral_div, (moment3 k).2]; field_simp
  · rw [heq, setIntegral_congr_fun measurableSet_Ioi (fun t _ => hpt t), integral_sub hm1 hρi,
      integral_div, (gammaMoment k).2]

lemma ind_part (k : ℕ) (N : ℝ) (hN : 1 ≤ N) :
    ∫ w in Ioc 0 N, Real.log (N / w) ^ k * (Ioi (1:ℝ)).indicator (fun w => 1 / w) w
      = Real.log N ^ (k + 1) / ((k : ℝ) + 1) := by
  have hN0 : 0 < N := by linarith
  have e : (fun w => Real.log (N / w) ^ k * (Ioi (1:ℝ)).indicator (fun w => 1 / w) w) =
      (Ioi (1:ℝ)).indicator (fun w => Real.log (N / w) ^ k * (1 / w)) := by
    funext w; rw [indicator_mul_right]
  have hset : Ioc 0 N ∩ Ioi (1:ℝ) = Ioc 1 N := by
    ext w; simp only [mem_inter_iff, mem_Ioc, mem_Ioi]
    constructor
    · rintro ⟨⟨_, h2⟩, h3⟩; exact ⟨h3, h2⟩
    · rintro ⟨h1, h2⟩; exact ⟨⟨by linarith, h2⟩, h1⟩
  rw [e, setIntegral_indicator measurableSet_Ioi, hset, ← intervalIntegral.integral_of_le hN]
  have hderiv : ∀ x ∈ uIcc 1 N, HasDerivAt (fun w => -(Real.log (N / w) ^ (k + 1)) / ((k : ℝ) + 1))
      (Real.log (N / x) ^ k * (1 / x)) x := by
    intro x hx
    rw [uIcc_of_le hN] at hx
    have hx0 : 0 < x := by linarith [hx.1]
    have h1 : HasDerivAt (fun w => N / w) (-(N / x ^ 2)) x := by
      have := (hasDerivAt_inv hx0.ne').const_mul N
      convert this using 1
      · funext w; rw [div_eq_mul_inv]
      · ring
    have h2 := (h1.log (by positivity)).pow (k + 1)
    have h3 := (h2.neg).div_const ((k : ℝ) + 1)
    convert h3 using 1
    have hk : ((k : ℝ) + 1) ≠ 0 := by positivity
    push_cast
    field_simp
  have hcont : ContinuousOn (fun x => Real.log (N / x) ^ k * (1 / x)) (uIcc 1 N) := by
    intro x hx
    rw [uIcc_of_le hN] at hx
    have hx0 : 0 < x := by linarith [hx.1]
    apply ContinuousAt.continuousWithinAt
    exact (((continuousAt_const.div continuousAt_id hx0.ne').log (div_pos hN0 hx0).ne').pow k).mul
      (continuousAt_const.div continuousAt_id hx0.ne')
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hderiv (hcont.intervalIntegrable)]
  rw [div_self hN0.ne', Real.log_one, div_one]
  simp only [ne_eq, Nat.add_eq_zero_iff, one_ne_zero, and_false, not_false_eq_true, zero_pow,
    neg_zero, zero_div, zero_sub]
  ring

lemma mom_eq (m : ℕ) : ∫ w in Ioi (0:ℝ), (-Real.log w) ^ m * hh w = (m.factorial : ℝ) * bm m := by
  unfold bm
  rw [← integral_const_mul]
  refine setIntegral_congr_fun measurableSet_Ioi (fun w _ => ?_)
  have : (m.factorial : ℝ) ≠ 0 := by exact_mod_cast m.factorial_ne_zero
  field_simp

lemma stmt5_main : stmt_5 := by
  intro k N hN
  have hN0 : 0 < N := by linarith
  set L := Real.log N with hLdef
  have hintg := (Tint_subst k N hN0).2 (fsub_integrable k N hN0)
  refine ⟨hintg, ?_⟩
  -- binomial expansion
  have hexp : ∀ w ∈ Ioi (0:ℝ), Real.log (N / w) ^ k * hh w =
      ∑ m ∈ Finset.range (k + 1), (L ^ (k - m) * (k.choose m : ℝ)) * ((-Real.log w) ^ m * hh w) := by
    intro w hw
    have hw0 : (0:ℝ) < w := hw
    rw [Real.log_div hN0.ne' hw0.ne', show L - Real.log w = -Real.log w + L by rw [hLdef]; ring,
      add_pow, Finset.sum_mul]
    refine Finset.sum_congr rfl (fun m _ => ?_)
    ring
  have hI0' : IntegrableOn (fun w => ∑ m ∈ Finset.range (k + 1),
      (L ^ (k - m) * (k.choose m : ℝ)) * ((-Real.log w) ^ m * hh w)) (Ioi 0) :=
    integrable_finsetSum _ (fun m _ => (mom_int m).const_mul _)
  have hI0 : IntegrableOn (fun w => Real.log (N / w) ^ k * hh w) (Ioi 0) :=
    hI0'.congr_fun (fun w hw => (hexp w hw).symm) measurableSet_Ioi
  have hval0 : ∫ w in Ioi (0:ℝ), Real.log (N / w) ^ k * hh w =
      ∑ m ∈ Finset.range (k + 1), (L ^ (k - m) * (k.choose m : ℝ)) * ((m.factorial : ℝ) * bm m) := by
    rw [setIntegral_congr_fun measurableSet_Ioi hexp,
      integral_finsetSum _ (fun m _ => (mom_int m).const_mul _)]
    refine Finset.sum_congr rfl (fun m _ => ?_)
    rw [integral_const_mul, mom_eq]
  -- tail
  obtain ⟨-, J, hJ0, hJ1, hJ⟩ := tail_rr k N hN
  have htail : ∫ w in Ioi N, Real.log (N / w) ^ k * hh w =
      -((-1 : ℝ) ^ k * ((k.factorial : ℝ) / (2 * N) - J)) := by
    rw [← hJ, ← integral_const_mul, ← integral_neg]
    refine setIntegral_congr_fun measurableSet_Ioi (fun w hw => ?_)
    have hw1 : 1 < w := lt_of_le_of_lt hN hw
    have hw0 : 0 < w := by linarith
    rw [hh_large w hw1, show N / w = (w / N)⁻¹ by field_simp, Real.log_inv, neg_pow]
    ring
  have hsplit := split_Ioi _ N hN0.le hI0
  -- the indicator part
  have hind := ind_part k N hN
  have hhint : IntegrableOn (fun w => Real.log (N / w) ^ k * hh w) (Ioc 0 N) :=
    hI0.mono_set Ioc_subset_Ioi_self
  have hg_eq : ∀ w, Real.log (N / w) ^ k * g w = Real.log (N / w) ^ k * hh w +
      Real.log (N / w) ^ k * (Ioi (1:ℝ)).indicator (fun w => 1 / w) w := by
    intro w; unfold hh; ring
  have hindint : IntegrableOn (fun w => Real.log (N / w) ^ k *
      (Ioi (1:ℝ)).indicator (fun w => 1 / w) w) (Ioc 0 N) :=
    (hintg.sub hhint).congr_fun (fun w _ => by simp only [Pi.sub_apply]; rw [hg_eq w]; ring)
      measurableSet_Ioc
  have hT : T k N = (∫ w in Ioc 0 N, Real.log (N / w) ^ k * hh w) + L ^ (k + 1) / ((k : ℝ) + 1) := by
    unfold T
    rw [intervalIntegral.integral_of_le hN0.le, ← hind, ← integral_add hhint hindint]
    exact setIntegral_congr_fun measurableSet_Ioc (fun w _ => hg_eq w)
  -- theta
  have hB : 0 < (k.factorial : ℝ) / (8 * 3 ^ (k + 1) * N ^ 3) := by
    have : (0:ℝ) < k.factorial := by exact_mod_cast k.factorial_pos
    positivity
  refine ⟨J / ((k.factorial : ℝ) / (8 * 3 ^ (k + 1) * N ^ 3)),
    ⟨div_nonneg hJ0 hB.le, (div_le_one hB).2 hJ1⟩, ?_⟩
  rw [hT]
  have hc : ∫ w in Ioc 0 N, Real.log (N / w) ^ k * hh w =
      (∫ w in Ioi (0:ℝ), Real.log (N / w) ^ k * hh w) - ∫ w in Ioi N, Real.log (N / w) ^ k * hh w := by
    rw [hsplit]; ring
  rw [hc, hval0, htail]
  unfold poly
  have hsum : ∑ m ∈ Finset.range (k + 1), (L ^ (k - m) * (k.choose m : ℝ)) * ((m.factorial : ℝ) * bm m)
      = ∑ i ∈ Finset.range (k + 1), (k.descFactorial i : ℝ) * mt (i + 1) * L ^ (k - i) := by
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [mt_eq, Nat.descFactorial_eq_factorial_mul_choose]
    push_cast; ring
  rw [hsum]
  field_simp
  ring

end X06

theorem X06.proof_5 : X06.stmt_5 := X06.stmt5_main

/-! ### The decomposition (`stmt_1`): integrability parts -/

namespace X06
open Real MeasureTheory Set

lemma area_integrable (n : ℝ) (hn : 1 < n) :
    IntegrableOn (areaIntegrand n) (Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1) := by
  have hc1 : Continuous (fun x : ℝ => x ^ (n - 1)) := Real.continuous_rpow_const (by linarith)
  have hc2 : Continuous (fun x : ℝ => x ^ n) := Real.continuous_rpow_const (by linarith)
  have hc : Continuous (areaIntegrand n) := by
    unfold areaIntegrand
    exact Real.continuous_sqrt.comp ((continuous_const.add
      (((continuous_const.mul (hc1.comp continuous_fst)).mul (hc2.comp continuous_snd)).pow 2)).add
      (((continuous_const.mul (hc2.comp continuous_fst)).mul (hc1.comp continuous_snd)).pow 2))
  exact hc.continuousOn.integrableOn_compact (isCompact_Icc.prod isCompact_Icc)

lemma Cinner_integrable (n : ℝ) (hn : 1 < n) : IntegrableOn (Cinner n) (Ioi 0) := by
  obtain ⟨hA, hB, hC, -, -⟩ := X06.proof_8 n hn
  exact ((hB.add hC).add hA).congr_fun (fun σ _ => (Cinner_eq n σ).symm) measurableSet_Ioi

end X06

namespace X06
open Real MeasureTheory Set Filter Topology

/-! ### The exact decomposition (`stmt_1`) -/

/-- `z_x² + z_y²`. -/
noncomputable def Qs (n : ℝ) (p : ℝ × ℝ) : ℝ :=
  (n * p.1 ^ (n - 1) * p.2 ^ n) ^ 2 + (n * p.1 ^ n * p.2 ^ (n - 1)) ^ 2

/-- `G = √(z_x² + z_y²)`. -/
noncomputable def Gs (n : ℝ) (p : ℝ × ℝ) : ℝ := Real.sqrt (Qs n p)

lemma Qs_continuous (n : ℝ) (hn : 1 < n) : Continuous (Qs n) := by
  have hc1 : Continuous (fun x : ℝ => x ^ (n - 1)) := Real.continuous_rpow_const (by linarith)
  have hc2 : Continuous (fun x : ℝ => x ^ n) := Real.continuous_rpow_const (by linarith)
  unfold Qs
  exact (((continuous_const.mul (hc1.comp continuous_fst)).mul (hc2.comp continuous_snd)).pow 2).add
      (((continuous_const.mul (hc2.comp continuous_fst)).mul (hc1.comp continuous_snd)).pow 2)

lemma Gs_continuous (n : ℝ) (hn : 1 < n) : Continuous (Gs n) :=
  Real.continuous_sqrt.comp (Qs_continuous n hn)

lemma Gs_nonneg (n : ℝ) (p : ℝ × ℝ) : 0 ≤ Gs n p := Real.sqrt_nonneg _

lemma area_eq (n : ℝ) (p : ℝ × ℝ) : areaIntegrand n p = 1 + Gs n p - f (Gs n p) := by
  unfold areaIntegrand f Gs
  have hQ : 0 ≤ Qs n p := by unfold Qs; positivity
  rw [Real.sq_sqrt hQ]
  unfold Qs; ring_nf

lemma Gs_symm (n x y : ℝ) : Gs n (y, x) = Gs n (x, y) := by
  unfold Gs Qs; simp only; ring_nf

/-- `G(e^{-s/n}, e^{-t/n}) = W(s+t, s-t)`. -/
lemma Gs_exp (n : ℝ) (hn : 1 < n) (s t : ℝ) :
    Gs n (Real.exp (-s / n), Real.exp (-t / n)) = W n (s + t) (s - t) := by
  have hn0 : 0 < n := by linarith
  unfold Gs Qs
  simp only
  rw [← Real.exp_mul, ← Real.exp_mul, ← Real.exp_mul, ← Real.exp_mul]
  have hW : 0 ≤ W n (s + t) (s - t) := by
    unfold W NN; have := Real.sqrt_nonneg (Real.cosh ((s - t) / n)); positivity
  rw [Real.sqrt_eq_iff_eq_sq (by positivity) hW]
  unfold W NN
  have hne : n ≠ 0 := hn0.ne'
  have h1 : Real.exp (-s / n * (n - 1)) = Real.exp (-s) * Real.exp (s / n) := by
    rw [← Real.exp_add]; congr 1; field_simp; ring
  have h2 : Real.exp (-t / n * (n - 1)) = Real.exp (-t) * Real.exp (t / n) := by
    rw [← Real.exp_add]; congr 1; field_simp; ring
  have h3 : Real.exp (-s / n * n) = Real.exp (-s) := by congr 1; field_simp
  have h4 : Real.exp (-t / n * n) = Real.exp (-t) := by congr 1; field_simp
  have h5 : Real.exp (-(s + t)) = Real.exp (-s) * Real.exp (-t) := by
    rw [← Real.exp_add]; congr 1; ring
  have h6 : Real.exp ((s + t) / n / 2) ^ 2 = Real.exp (s / n) * Real.exp (t / n) := by
    rw [sq, ← Real.exp_add, ← Real.exp_add]; congr 1; ring
  have h7 : Real.cosh ((s - t) / n) =
      (Real.exp (s / n) / Real.exp (t / n) + Real.exp (t / n) / Real.exp (s / n)) / 2 := by
    rw [Real.cosh_eq, ← Real.exp_sub, ← Real.exp_sub]; congr 3 <;> ring_nf
  rw [h1, h2, h3, h4, h5]
  rw [show (Real.sqrt 2 * n * (Real.exp (-s) * Real.exp (-t)) * (Real.exp ((s + t) / n / 2) *
      Real.sqrt (Real.cosh ((s - t) / n)))) ^ 2 = Real.sqrt 2 ^ 2 * n ^ 2 *
      (Real.exp (-s) * Real.exp (-t)) ^ 2 * (Real.exp ((s + t) / n / 2) ^ 2 *
      Real.sqrt (Real.cosh ((s - t) / n)) ^ 2) by ring]
  rw [Real.sq_sqrt (by norm_num), h6, Real.sq_sqrt (Real.cosh_pos _).le, h7]
  have := Real.exp_pos (s / n)
  have := Real.exp_pos (t / n)
  field_simp


/-- substitution `y = e^{-t/n}` from `(0,∞)` onto `(0,1)`. -/
lemma subst_en (n : ℝ) (hn : 0 < n) (G : ℝ → ℝ) :
    ∫ y in Ioo 0 1, G y = ∫ t in Ioi (0:ℝ), 1 / n * Real.exp (-t / n) * G (Real.exp (-t / n)) := by
  have himg : (fun t : ℝ => Real.exp (-t / n)) '' Ioi 0 = Ioo 0 1 := by
    ext y; constructor
    · rintro ⟨t, ht, rfl⟩
      have ht' : (0:ℝ) < t := ht
      exact ⟨Real.exp_pos _, Real.exp_lt_one_iff.2 (by rw [neg_div]; exact neg_neg_of_pos (by positivity))⟩
    · rintro ⟨h0, h1⟩
      refine ⟨-n * Real.log y, ?_, ?_⟩
      · have := Real.log_neg h0 h1
        show 0 < -n * Real.log y
        nlinarith
      · simp only
        rw [show -(-n * Real.log y) / n = Real.log y by field_simp, Real.exp_log h0]
  have hd : ∀ x ∈ Ioi (0:ℝ), HasDerivWithinAt (fun t : ℝ => Real.exp (-t / n))
      (-(1 / n * Real.exp (-x / n))) (Ioi 0) x := by
    intro x _
    have h := (Real.hasDerivAt_exp (-x / n)).comp x ((hasDerivAt_neg x).div_const n)
    have h' : HasDerivAt (fun t : ℝ => Real.exp (-t / n)) (-(1 / n * Real.exp (-x / n))) x := by
      convert h using 1
      · funext t; simp only [Function.comp]
      · ring
    exact h'.hasDerivWithinAt
  have hinj : InjOn (fun t : ℝ => Real.exp (-t / n)) (Ioi 0) := by
    intro x _ y _ h
    have := Real.exp_injective h
    field_simp at this
    linarith
  rw [← himg, integral_image_eq_integral_abs_deriv_smul measurableSet_Ioi hd hinj]
  refine setIntegral_congr_fun measurableSet_Ioi (fun x _ => ?_)
  rw [abs_neg, abs_of_pos (by positivity), smul_eq_mul]

/-- swapping the order of integration over `0 < s < σ`. -/
lemma tri_swap (F : ℝ → ℝ → ℝ)
    (hF : Integrable (fun p : ℝ × ℝ => {p : ℝ × ℝ | 0 < p.1 ∧ p.1 < p.2}.indicator
      (fun p => F p.1 p.2) p) (volume.prod volume)) :
    ∫ s in Ioi (0:ℝ), ∫ σ in Ioi s, F s σ = ∫ σ in Ioi (0:ℝ), ∫ s in Ioo 0 σ, F s σ := by
  set R := {p : ℝ × ℝ | 0 < p.1 ∧ p.1 < p.2}
  have h := integral_integral_swap (μ := volume) (ν := volume)
    (f := fun s σ => R.indicator (fun p => F p.1 p.2) (s, σ)) hF
  have hL : ∀ s : ℝ, ∫ σ, R.indicator (fun p => F p.1 p.2) (s, σ) =
      (Ioi (0:ℝ)).indicator (fun s => ∫ σ in Ioi s, F s σ) s := by
    intro s
    by_cases hs : 0 < s
    · rw [indicator_of_mem (show s ∈ Ioi (0:ℝ) from hs), ← integral_indicator measurableSet_Ioi]
      congr 1; funext σ
      by_cases hσ : s < σ
      · rw [indicator_of_mem (show (s, σ) ∈ R from ⟨hs, hσ⟩), indicator_of_mem (show σ ∈ Ioi s from hσ)]
      · rw [indicator_of_notMem (show (s, σ) ∉ R from fun h => hσ (show 0 < s ∧ s < σ from h).2),
          indicator_of_notMem (show σ ∉ Ioi s from hσ)]
    · rw [indicator_of_notMem (show s ∉ Ioi (0:ℝ) from hs)]
      have : (fun σ => R.indicator (fun p => F p.1 p.2) (s, σ)) = fun _ => 0 := by
        funext σ; exact indicator_of_notMem (fun h => hs (show 0 < s ∧ s < σ from h).1) _
      rw [this, integral_zero]
  have hR : ∀ σ : ℝ, ∫ s, R.indicator (fun p => F p.1 p.2) (s, σ) =
      (Ioi (0:ℝ)).indicator (fun σ => ∫ s in Ioo 0 σ, F s σ) σ := by
    intro σ
    by_cases hσ : 0 < σ
    · rw [indicator_of_mem (show σ ∈ Ioi (0:ℝ) from hσ), ← integral_indicator measurableSet_Ioo]
      congr 1
    · rw [indicator_of_notMem (show σ ∉ Ioi (0:ℝ) from hσ)]
      have : (fun s => R.indicator (fun p => F p.1 p.2) (s, σ)) = fun _ => 0 := by
        funext s; exact indicator_of_notMem (fun h => hσ (lt_trans (show 0 < s ∧ s < σ from h).1 (show 0 < s ∧ s < σ from h).2)) _
      rw [this, integral_zero]
  simp only [hL, hR] at h
  rw [integral_indicator measurableSet_Ioi, integral_indicator measurableSet_Ioi] at h
  exact h


/-- the integrand `e^{-σ/n} f(W(σ, 2s-σ))` on `0 < s < σ`. -/
noncomputable def Kt (n s σ : ℝ) : ℝ := Real.exp (-σ / n) * f (W n σ (2 * s - σ))

lemma W_nonneg (n σ δ : ℝ) (hn : 0 < n) : 0 ≤ W n σ δ := by
  unfold W NN; have := Real.sqrt_nonneg (Real.cosh (δ / n)); positivity

lemma Kt_integrable (n : ℝ) (hn : 1 < n) :
    Integrable (fun p : ℝ × ℝ => {p : ℝ × ℝ | 0 < p.1 ∧ p.1 < p.2}.indicator
      (fun p => Kt n p.1 p.2) p) (volume.prod volume) := by
  have hn0 : 0 < n := by linarith
  set c := 1 / (2 * n) with hcdef
  have hc : 0 < c := by positivity
  have hI : Integrable ((Ioi (0:ℝ)).indicator (fun s => Real.exp (-c * s))) :=
    (exp_neg_integrableOn_Ioi 0 hc).integrable_indicator measurableSet_Ioi
  have hb := Integrable.mul_prod hI hI
  refine hb.mono' ?_ ?_
  · have hK : Continuous (fun p : ℝ × ℝ => Kt n p.1 p.2) := by
      unfold Kt; exact (by fun_prop : Continuous fun p : ℝ × ℝ => Real.exp (-p.2 / n)).mul
        (f_continuous.comp (by unfold W; fun_prop))
    have hR : MeasurableSet {p : ℝ × ℝ | 0 < p.1 ∧ p.1 < p.2} :=
      (measurableSet_lt (measurable_const (a := (0:ℝ))) measurable_fst).inter
        (measurableSet_lt measurable_fst measurable_snd)
    exact (hK.measurable.indicator hR).aestronglyMeasurable
  · refine Eventually.of_forall (fun p => ?_)
    by_cases hp : 0 < p.1 ∧ p.1 < p.2
    · rw [indicator_of_mem (show p ∈ {p : ℝ × ℝ | 0 < p.1 ∧ p.1 < p.2} from hp),
        indicator_of_mem (show p.1 ∈ Ioi (0:ℝ) from hp.1),
        indicator_of_mem (show p.2 ∈ Ioi (0:ℝ) from lt_trans hp.1 hp.2)]
      unfold Kt
      have h0 := f_nonneg (W n p.2 (2 * p.1 - p.2)) (W_nonneg _ _ _ hn0)
      have h1 := f_le_one (W n p.2 (2 * p.1 - p.2))
      rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
      calc Real.exp (-p.2 / n) * f (W n p.2 (2 * p.1 - p.2)) ≤ Real.exp (-p.2 / n) * 1 :=
            mul_le_mul_of_nonneg_left h1 (by positivity)
        _ = Real.exp (-c * p.2) * Real.exp (-c * p.2) := by
            rw [mul_one, ← Real.exp_add]; congr 1; rw [hcdef]; field_simp; ring
        _ ≤ Real.exp (-c * p.1) * Real.exp (-c * p.2) := by
            exact mul_le_mul_of_nonneg_right (Real.exp_le_exp.2 (by nlinarith [hp.2]))
              (Real.exp_pos _).le
    · rw [indicator_of_notMem (show p ∉ {p : ℝ × ℝ | 0 < p.1 ∧ p.1 < p.2} from hp), norm_zero]
      have : ∀ x : ℝ, 0 ≤ (Ioi (0:ℝ)).indicator (fun s => Real.exp (-c * s)) x := fun x =>
        indicator_nonneg (fun _ _ => (Real.exp_pos _).le) x
      exact mul_nonneg (this _) (this _)

lemma inner_Kt (n σ : ℝ) (hn : 1 < n) (hσ : 0 < σ) : ∫ s in Ioo 0 σ, Kt n s σ = Cinner n σ := by
  have hn0 : 0 < n := by linarith
  unfold Kt Cinner
  rw [integral_const_mul]
  congr 1
  rw [← integral_Ioc_eq_integral_Ioo, ← intervalIntegral.integral_of_le hσ.le]
  rw [intervalIntegral.integral_comp_mul_sub (f := fun δ => f (W n σ δ)) two_ne_zero σ]
  simp only [mul_zero, zero_sub, smul_eq_mul]
  rw [show 2 * σ - σ = σ by ring]
  have hev : ∀ δ, f (W n σ (-δ)) = f (W n σ δ) := by
    intro δ; unfold W; rw [neg_div, Real.cosh_neg]
  have hc : Continuous (fun δ => f (W n σ δ)) := f_continuous.comp (by unfold W; fun_prop)
  rw [← intervalIntegral.integral_add_adjacent_intervals (hc.intervalIntegrable (-σ) 0)
    (hc.intervalIntegrable 0 σ)]
  have h := intervalIntegral.integral_comp_neg (a := 0) (b := σ) (fun δ => f (W n σ δ))
  simp only [neg_zero, hev] at h
  rw [← h]
  ring

lemma partII (n : ℝ) (hn : 1 < n) :
    ∫ p in Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1, f (Gs n p) = (1 / n) ^ 2 * Cint n := by
  have hn0 : 0 < n := by linarith
  have hcont : Continuous (fun p => f (Gs n p)) := f_continuous.comp (Gs_continuous n hn)
  have hint : IntegrableOn (fun p => f (Gs n p)) (Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1)
      (volume.prod volume) :=
    hcont.continuousOn.integrableOn_compact (isCompact_Icc.prod isCompact_Icc)
  have hF := setIntegral_prod (fun p => f (Gs n p)) hint
  rw [show (∫ p in Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1, f (Gs n p)) =
      ∫ p in Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1, f (Gs n p) ∂(volume.prod volume) from rfl, hF]
  simp only [integral_Icc_eq_integral_Ioo]
  have step1 : ∀ x : ℝ, ∫ y in Ioo (0:ℝ) 1, f (Gs n (x, y)) = ∫ t in Ioi (0:ℝ),
      1 / n * Real.exp (-t / n) * f (Gs n (x, Real.exp (-t / n))) :=
    fun x => subst_en n hn0 (fun y => f (Gs n (x, y)))
  simp only [step1]
  rw [subst_en n hn0]
  have step2 : ∀ s : ℝ, 1 / n * Real.exp (-s / n) * ∫ t in Ioi (0:ℝ), 1 / n * Real.exp (-t / n) *
      f (Gs n (Real.exp (-s / n), Real.exp (-t / n))) = (1 / n) ^ 2 * ∫ σ in Ioi s, Kt n s σ := by
    intro s
    rw [(shift_Ioi (fun σ => Kt n s σ) s).2, ← integral_const_mul, ← integral_const_mul]
    congr 1; funext t
    rw [Gs_exp n hn]
    unfold Kt
    rw [show 2 * s - (t + s) = s - t by ring, show t + s = s + t by ring,
      show -(s + t) / n = -s / n + -t / n by ring, Real.exp_add]
    ring
  simp only [step2]
  rw [integral_const_mul, tri_swap _ (Kt_integrable n hn)]
  congr 1
  exact setIntegral_congr_fun measurableSet_Ioi (fun σ hσ => inner_Kt n σ hn hσ)


lemma Gs_scale (n x u : ℝ) (hn : 1 < n) (hx : 0 < x) (hu : 0 ≤ u) :
    Gs n (x, x * u) = n * x ^ (2 * n - 1) * (u ^ (n - 1) * Real.sqrt (1 + u ^ 2)) := by
  have hn0 : 0 < n := by linarith
  unfold Gs Qs
  simp only
  have hR : 0 ≤ n * x ^ (2 * n - 1) * (u ^ (n - 1) * Real.sqrt (1 + u ^ 2)) := by positivity
  rw [Real.sqrt_eq_iff_eq_sq (by positivity) hR]
  rw [Real.mul_rpow hx.le hu, Real.mul_rpow hx.le hu,
    show 2 * n - 1 = (n - 1) + n by ring, Real.rpow_add hx,
    show u ^ n = u ^ (n - 1) * u by
      rw [← Real.rpow_add_one' hu (by linarith : n - 1 + 1 ≠ 0)]; ring_nf]
  rw [show (n * (x ^ (n - 1) * x ^ n) * (u ^ (n - 1) * Real.sqrt (1 + u ^ 2))) ^ 2 =
      n ^ 2 * (x ^ (n - 1) * x ^ n) ^ 2 * (u ^ (n - 1)) ^ 2 * Real.sqrt (1 + u ^ 2) ^ 2 by ring,
    Real.sq_sqrt (by positivity)]
  ring

lemma Phi_eq (n : ℝ) : Phi n = ∫ u in (0:ℝ)..1, u ^ (n - 1) * Real.sqrt (1 + u ^ 2) := rfl

lemma inner_low (n x : ℝ) (hn : 1 < n) (hx : 0 < x) :
    ∫ y in Ioc 0 x, Gs n (x, y) = n * x ^ (2 * n) * Phi n := by
  rw [← intervalIntegral.integral_of_le hx.le]
  have h := intervalIntegral.integral_comp_mul_left (a := 0) (b := 1) (fun y => Gs n (x, y)) hx.ne'
  simp only [mul_zero, mul_one, smul_eq_mul] at h
  rw [eq_inv_mul_iff_mul_eq₀ hx.ne'] at h
  rw [← h]
  have e : ∫ u in (0:ℝ)..1, Gs n (x, x * u) =
      ∫ u in (0:ℝ)..1, n * x ^ (2 * n - 1) * (u ^ (n - 1) * Real.sqrt (1 + u ^ 2)) := by
    refine intervalIntegral.integral_congr (fun u hu => ?_)
    rw [uIcc_of_le zero_le_one] at hu
    exact Gs_scale n x u hn hx hu.1
  rw [e, intervalIntegral.integral_const_mul, ← Phi_eq]
  rw [show 2 * n = (2 * n - 1) + 1 by ring, Real.rpow_add_one hx.ne']
  ring

lemma Gs_bound (n : ℝ) (hn : 1 < n) : ∃ M : ℝ, 0 ≤ M ∧
    ∀ p ∈ Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1, |Gs n p| ≤ M := by
  obtain ⟨C, hC⟩ := (isCompact_Icc.prod isCompact_Icc).exists_bound_of_continuousOn
    (Gs_continuous n hn).continuousOn (f := Gs n) (s := Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1)
  refine ⟨max C 0, le_max_right _ _, fun p hp => ?_⟩
  have := hC p hp
  rw [Real.norm_eq_abs] at this
  exact this.trans (le_max_left _ _)

lemma partI_upper (n : ℝ) (hn : 1 < n) :
    ∫ x in Ioo (0:ℝ) 1, ∫ y in Ioo x 1, Gs n (x, y) = n * Phi n / (2 * n + 1) := by
  have hn0 : 0 < n := by linarith
  set F : ℝ → ℝ → ℝ := fun s σ => (Iio (1:ℝ)).indicator (fun σ => Gs n (s, σ)) σ with hFdef
  set R := {p : ℝ × ℝ | 0 < p.1 ∧ p.1 < p.2}
  have hRm : MeasurableSet R :=
    (measurableSet_lt (measurable_const (a := (0:ℝ))) measurable_fst).inter
      (measurableSet_lt measurable_fst measurable_snd)
  have hFint : Integrable (fun p : ℝ × ℝ => R.indicator (fun p => F p.1 p.2) p)
      (volume.prod volume) := by
    obtain ⟨M, hM0, hM⟩ := Gs_bound n hn
    have hI : Integrable ((Icc (0:ℝ) 1).indicator (fun _ => (1:ℝ))) :=
      (integrableOn_const (by simp)).integrable_indicator measurableSet_Icc
    have hb := (Integrable.mul_prod hI hI).const_mul M
    refine hb.mono' ?_ ?_
    · have e : (fun p : ℝ × ℝ => F p.1 p.2) = (Prod.snd ⁻¹' Iio (1:ℝ)).indicator (Gs n) := by
        funext p; simp only [hFdef, indicator, mem_preimage, mem_Iio]
      rw [e]
      exact (((Gs_continuous n hn).measurable.indicator
        (measurable_snd measurableSet_Iio)).indicator hRm).aestronglyMeasurable
    · refine Eventually.of_forall (fun p => ?_)
      by_cases hp : p ∈ R ∧ p.2 < 1
      · have hp1 : p.1 ∈ Icc (0:ℝ) 1 := ⟨(show 0 < p.1 ∧ p.1 < p.2 from hp.1).1.le,
          by linarith [(show 0 < p.1 ∧ p.1 < p.2 from hp.1).2, hp.2]⟩
        have hp2 : p.2 ∈ Icc (0:ℝ) 1 := ⟨by linarith [(show 0 < p.1 ∧ p.1 < p.2 from hp.1).1,
          (show 0 < p.1 ∧ p.1 < p.2 from hp.1).2], hp.2.le⟩
        rw [indicator_of_mem hp.1, hFdef]
        simp only
        rw [indicator_of_mem (show p.2 ∈ Iio (1:ℝ) from hp.2), indicator_of_mem hp1,
          indicator_of_mem hp2, Real.norm_eq_abs, mul_one, mul_one]
        exact hM (p.1, p.2) ⟨hp1, hp2⟩
      · have h0 : R.indicator (fun p => F p.1 p.2) p = 0 := by
          by_cases hR : p ∈ R
          · rw [indicator_of_mem hR, hFdef]
            simp only
            exact indicator_of_notMem (fun h => hp ⟨hR, h⟩) _
          · exact indicator_of_notMem hR _
        rw [h0, norm_zero]
        have : ∀ x : ℝ, 0 ≤ (Icc (0:ℝ) 1).indicator (fun _ => (1:ℝ)) x := fun x =>
          indicator_nonneg (fun _ _ => zero_le_one) x
        exact mul_nonneg hM0 (mul_nonneg (this _) (this _))
  have hsw := tri_swap F hFint
  -- left side
  have hL : ∀ s : ℝ, ∫ σ in Ioi s, F s σ = ∫ σ in Ioo s 1, Gs n (s, σ) := by
    intro s
    rw [hFdef]; simp only
    rw [setIntegral_indicator measurableSet_Iio, Ioi_inter_Iio]
  have hR' : ∀ σ : ℝ, ∫ s in Ioo 0 σ, F s σ = (Iio (1:ℝ)).indicator
      (fun σ => ∫ s in Ioo 0 σ, Gs n (s, σ)) σ := by
    intro σ
    by_cases hσ : σ < 1
    · rw [indicator_of_mem (show σ ∈ Iio (1:ℝ) from hσ)]
      refine setIntegral_congr_fun measurableSet_Ioo (fun s _ => ?_)
      rw [hFdef]; simp only
      rw [indicator_of_mem (show σ ∈ Iio (1:ℝ) from hσ)]
    · rw [indicator_of_notMem (show σ ∉ Iio (1:ℝ) from hσ)]
      have : ∀ s, F s σ = 0 := fun s => by
        rw [hFdef]; simp only; exact indicator_of_notMem (show σ ∉ Iio (1:ℝ) from hσ) _
      simp only [this, integral_zero]
  simp only [hL, hR'] at hsw
  have hLHS : ∫ s in Ioi (0:ℝ), ∫ σ in Ioo s 1, Gs n (s, σ) =
      ∫ x in Ioo (0:ℝ) 1, ∫ y in Ioo x 1, Gs n (x, y) := by
    refine setIntegral_eq_of_subset_of_forall_sdiff_eq_zero measurableSet_Ioi Ioo_subset_Ioi_self
      (fun s hs => ?_)
    have h1 : 1 ≤ s := by
      by_contra h; exact hs.2 ⟨hs.1, by linarith⟩
    rw [Ioo_eq_empty (by linarith), Measure.restrict_empty, integral_zero_measure]
  rw [hLHS] at hsw
  rw [hsw, setIntegral_indicator measurableSet_Iio, Ioi_inter_Iio]
  have hsym : ∀ σ ∈ Ioo (0:ℝ) 1, ∫ s in Ioo 0 σ, Gs n (s, σ) = n * σ ^ (2 * n) * Phi n := by
    intro σ hσ
    rw [← integral_Ioc_eq_integral_Ioo, ← inner_low n σ hn hσ.1]
    refine setIntegral_congr_fun measurableSet_Ioc (fun s _ => ?_)
    rw [Gs_symm]
  rw [setIntegral_congr_fun measurableSet_Ioo hsym]
  rw [← integral_Ioc_eq_integral_Ioo, ← intervalIntegral.integral_of_le zero_le_one]
  have e : (fun σ : ℝ => n * σ ^ (2 * n) * Phi n) = fun σ => (n * Phi n) * σ ^ (2 * n) := by
    funext σ; ring
  rw [e, intervalIntegral.integral_const_mul, integral_rpow (Or.inl (by linarith))]
  rw [Real.one_rpow, Real.zero_rpow (by linarith), sub_zero]
  field_simp


lemma int_A (n : ℝ) (hn : 1 < n) :
    ∫ x in Ioo (0:ℝ) 1, n * x ^ (2 * n) * Phi n = n * Phi n / (2 * n + 1) := by
  have hn0 : 0 < n := by linarith
  rw [← integral_Ioc_eq_integral_Ioo, ← intervalIntegral.integral_of_le zero_le_one]
  have e : (fun σ : ℝ => n * σ ^ (2 * n) * Phi n) = fun σ => (n * Phi n) * σ ^ (2 * n) := by
    funext σ; ring
  rw [e, intervalIntegral.integral_const_mul, integral_rpow (Or.inl (by linarith))]
  rw [Real.one_rpow, Real.zero_rpow (by linarith), sub_zero]
  field_simp

lemma partI (n : ℝ) (hn : 1 < n) :
    ∫ p in Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1, Gs n p = 2 * n * Phi n / (2 * n + 1) := by
  have hn0 : 0 < n := by linarith
  have hcont := Gs_continuous n hn
  have hint : IntegrableOn (Gs n) (Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1) (volume.prod volume) :=
    hcont.continuousOn.integrableOn_compact (isCompact_Icc.prod isCompact_Icc)
  have hF := setIntegral_prod (Gs n) hint
  rw [show (∫ p in Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1, Gs n p) =
      ∫ p in Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1, Gs n p ∂(volume.prod volume) from rfl, hF]
  simp only [integral_Icc_eq_integral_Ioo]
  -- integrability of `x ↦ ∫ y, G(x,y)`
  have hint' : IntegrableOn (Gs n) (Ioo (0:ℝ) 1 ×ˢ Ioo (0:ℝ) 1) (volume.prod volume) :=
    hint.mono_set (prod_mono Ioo_subset_Icc_self Ioo_subset_Icc_self)
  rw [IntegrableOn, ← Measure.prod_restrict] at hint'
  have hI : IntegrableOn (fun x => ∫ y in Ioo (0:ℝ) 1, Gs n (x, y)) (Ioo (0:ℝ) 1) :=
    hint'.integral_prod_left
  have hAc : Continuous (fun x : ℝ => n * x ^ (2 * n) * Phi n) :=
    (continuous_const.mul (Real.continuous_rpow_const (by linarith))).mul continuous_const
  have hA : IntegrableOn (fun x : ℝ => n * x ^ (2 * n) * Phi n) (Ioo (0:ℝ) 1) :=
    (hAc.integrableOn_Icc (a := 0) (b := 1)).mono_set Ioo_subset_Icc_self
  have hsplit : ∀ x ∈ Ioo (0:ℝ) 1, ∫ y in Ioo x 1, Gs n (x, y) =
      (∫ y in Ioo (0:ℝ) 1, Gs n (x, y)) - n * x ^ (2 * n) * Phi n := by
    intro x hx
    have hcy : Continuous (fun y => Gs n (x, y)) := hcont.comp (by fun_prop)
    rw [← inner_low n x hn hx.1, ← Ioc_union_Ioo_eq_Ioo hx.1.le hx.2,
      setIntegral_union (Set.disjoint_left.2 (fun y h1 h2 => by linarith [h1.2, h2.1]))
        measurableSet_Ioo ((hcy.integrableOn_Icc (a := 0) (b := 1)).mono_set
          (Ioc_subset_Icc_self.trans (Icc_subset_Icc le_rfl hx.2.le)))
        ((hcy.integrableOn_Icc (a := 0) (b := 1)).mono_set
          (Ioo_subset_Icc_self.trans (Icc_subset_Icc hx.1.le le_rfl)))]
    ring
  have hup := partI_upper n hn
  rw [setIntegral_congr_fun measurableSet_Ioo hsplit, integral_sub hI hA, int_A n hn] at hup
  have key := eq_add_of_sub_eq hup
  convert key using 1
  ring


lemma decomp (n : ℝ) (hn : 1 < n) :
    A n - 1 = Phi n / (1 + 1 / n / 2) - (1 / n) ^ 2 * Cint n := by
  have hn0 : 0 < n := by linarith
  have hK : IsCompact (Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1) := isCompact_Icc.prod isCompact_Icc
  have hG := Gs_continuous n hn
  have h1 : IntegrableOn (fun _ : ℝ × ℝ => (1:ℝ)) (Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1) :=
    continuous_const.continuousOn.integrableOn_compact hK
  have h2 : IntegrableOn (Gs n) (Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1) := hG.continuousOn.integrableOn_compact hK
  have h3 : IntegrableOn (fun p => f (Gs n p)) (Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1) :=
    (f_continuous.comp hG).continuousOn.integrableOn_compact hK
  have hvol : ∫ p in Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1, (1:ℝ) = 1 := by
    rw [setIntegral_const, smul_eq_mul, mul_one, Measure.real, Measure.volume_eq_prod,
      Measure.prod_prod, Real.volume_Icc]
    simp
  unfold A
  rw [setIntegral_congr_fun (measurableSet_Icc.prod measurableSet_Icc) (fun p _ => area_eq n p)]
  have h12 : IntegrableOn (fun p => 1 + Gs n p) (Icc (0:ℝ) 1 ×ˢ Icc (0:ℝ) 1) := h1.add h2
  rw [integral_sub h12 h3, integral_add (f := fun _ => (1:ℝ)) (g := fun p => Gs n p) h1 h2, hvol,
    partI n hn, partII n hn]
  field_simp
  ring

end X06

theorem X06.proof_1 : X06.stmt_1 := fun n hn =>
  ⟨X06.area_integrable n hn, X06.Cinner_integrable n hn, X06.decomp n hn⟩

/-! ### Assembly -/

theorem X06.proof_11 : X06.stmt_11 := X06.stmt11_of X06.proof_5

theorem X06.proof_12 : X06.stmt_12 := X06.stmt12_of X06.proof_1 X06.proof_11

theorem X06.proof_13 : X06.stmt_13 := X06.stmt13_of_stmt3 X06.proof_3

theorem X06.proof_15 : X06.stmt_15 := X06.stmt15_of X06.proof_12 X06.proof_13

theorem X06.proof_16 : X06.stmt_16 := X06.stmt16_of X06.proof_15

theorem X06.proof_14 : X06.stmt_14 := X06.stmt14_of X06.proof_16
