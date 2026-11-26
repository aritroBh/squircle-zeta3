import Mathlib
import Zeta3.Wave1.T10_PiP
import Zeta3.Wave4.X02_PiP
/-!
# Q02 — the expansion of `π_p` to every order, its divergence and the growth law of `N_k`
(paper label `thm:pip`)

Source: `paper/ieee/main.tex`, Theorem `thm:pip` ("Expansion of $\pi_p$"), lines 1336–1350,
and the complement printed after it, lines 1373–1379; detailed proof in
`paper/ieee/proofs/sec_pip.tex`: notation (lines 43–48), eq. `eq:pip:Ehat` (lines 70–73),
`W` (line 79), Lemma `lem:pip:W` (lines 83–103), Lemma `lem:pip:expansion` (lines 105–133),
Lemma `lem:pip:J` (lines 139–160), Section `sec:pip:N45` (lines 191–343) with
Lemma `lem:pip:C45` (lines 227–329, also `proofs/sec_pip_C45.tex`), and Section
`sec:pip:div` (lines 369–588): Lemmas `lem:pip:laplace`, `lem:pip:Tk`, `lem:pip:hankel`,
`lem:pip:TkT1`, `lem:pip:twoterm`, `prop:pip:growth`, `prop:pip:div`.

The parts of `thm:pip` already proved in `Zeta3.Wave4.X02_PiP` (`π_p = 4·2^{-1/p} + E(p)` and the
expansion through `N_3` with `O(p^{-5})`) are not restated; `X02.piP`, `X02.E`, `X02.W`, `X02.Q`,
`X02.calL`, `X02.z3`, `X02.N0`–`X02.N3` are reused. This file covers the rest of the theorem
(main.tex lines 1344–1349):

> The expansion continues to every order, $\pi_p/4\sim1+\sum_{k\ge0}N_kp^{-k-1}$, and this series
> diverges for every $p$. As $k\to\infty$, for odd $k$
> \[N_k=(-1)^{\frac{k-1}2}\,\frac{2\,(k-2)!}{\pi^{k-1}}\Bigl(1+\frac{\pi^2}{2(k-2)}+O\Bigl(\frac{\log^2k}{k^2}\Bigr)\Bigr),\]
> and for even $k$
> \[N_k=(-1)^{\frac{k-2}2}\,\frac{2\,(k-3)!}{\pi^{k-2}}\Bigl(\psi(k-2)+\gamma-\log\pi+O\Bigl(\frac{\log^2k}{k}\Bigr)\Bigr),\]
> where $\psi=\Gamma'/\Gamma$ and $\gamma$ is Euler's constant.

and the complement (main.tex lines 1373–1379):

> N_4&=-\frac{\log^52}{120}+\frac{\pi^2}{24}\Bigl(1-\log2+\frac{\log^22}{2}\Bigr)\\
> &\quad+\frac{\zeta(3)}{4}(1-\log2)+\frac{\pi^4}{2880},\\
> N_5&=\frac{\log^62}{720}+\frac{\pi^2}{24}\Bigl(1-\log2+\frac{\log^22}{2}-\frac{\log^32}{6}\Bigr)\\
> &\quad+\frac{\zeta(3)}{4}\Bigl(1-\log2+\frac{\log^22}{2}\Bigr)+\frac{\pi^4}{2880}(1-\log2)\\
> &\quad-\frac{7\pi^4}{1920}-\frac{\zeta(3)\pi^2}{48}+\frac{\zeta(5)}{8}.

Intermediate steps (sec_pip.tex, quoted):

> (line 46) J_{r,j}=\int_0^\infty s^j\mathcal L(s)^r\,ds\qquad(r\ge1,\ j\ge0).
>   These converge because $0<\mathcal L(s)\le\eu^{-s}$.
> (line 72) \hat E(h)=h\kappa\int_0^\infty \eu^{-h\kappa s}(1+\eu^{-\kappa s})^{-1-h}\bigl[(1+\eu^{-s})^h-1\bigr]ds,
>   \qquad \frac{\pi_p}4=2^{-h}+\hat E(h).
> (line 118) Write $\hat E=\frac{h^2}{1-h}Q(h)$ with $Q(h)=\int_0^\infty W(h,sh)\bigl[(1+\eu^{-s})^h-1\bigr]ds$.
> (lines 84–92, Lemma lem:pip:W) Write $W(x,y)=\sum_{i,j\ge0}w_{ij}x^iy^j$ and $W_m(s)=\sum_{i+j=m}w_{ij}s^j$.
>   ... for every $K\ge0$ there are $C_K,h_K>0$ with
>   \Bigl|W(h,sh)-\sum_{m\le K}h^mW_m(s)\Bigr|\le C_K\bigl(h(1+s)\bigr)^{K+1}\qquad(0<h\le h_K,\ 0\le s\le h^{-1/2}).
>   Moreover $0<W(h,sh)\le1$ for $0<h<1$ and $s\ge0$.
> (lines 106–113, Lemma lem:pip:expansion) Let $q_n=\sum_{m+r=n,\ r\ge1}\frac1{r!}\int_0^\infty W_m(s)\mathcal L(s)^r\,ds$
>   for $n\ge1$. For every $K\ge1$, as $h\to0^+$,
>   \hat E(h)=\frac{h^2}{1-h}\Bigl(\sum_{n=1}^Kq_nh^n+O(h^{K+1})\Bigr).
>   Consequently $\pi_p/4\sim1+\sum_{k\ge0}N_kh^{k+1}$ with
>   N_k=\frac{(-\log2)^{k+1}}{(k+1)!}+\sum_{n=1}^{k-1}q_n\qquad(k\ge0),
> (lines 140–143, Lemma lem:pip:J(i)) $J_{1,j}=j!\,\eta(j+2)=j!\,(1-2^{-j-1})\zeta(j+2)$. In particular
>   J_{1,0}=\frac{\pi^2}{12},\quad J_{1,1}=\frac{3\zeta(3)}4,\quad J_{1,2}=\frac{7\pi^4}{360},\quad J_{1,3}=\frac{45\zeta(5)}8 .
> (lines 203–204) W_2(s)=\frac{\log^22}4-\frac{\log2}4\,s,\qquad
>   W_3(s)=-\frac{\log^32}{12}+\frac{\log^22}8\,s-\frac3{16}s^2-\frac1{48}s^3 .
> (lines 230–231, Lemma lem:pip:C45) \frac1{12}J_{3,0}+\frac18J_{2,1}=\frac{\pi^4}{2880},\qquad
>   \frac1{48}\bigl(J_{4,0}+2J_{3,1}\bigr)=-\frac{\zeta(3)\pi^2}{48}+\frac{31\,\zeta(5)}{128}.
> (lines 378–381, Lemma lem:pip:laplace) For $0<h<1$,
>   \hat E(h)=h\sum_{n\ge1}\binom hn\int_0^\infty \eu^{-ny/h}H_n(y,h)\,dy,\qquad H_n(y,h)=\eu^{(n+1)y}(1+\eu^{y})^{-1-h}.
> (lines 433–436, Lemma lem:pip:Tk) Let $\mathcal A(y,h)=\eu^{2y}(1+\eu^y)^{-1-h}=H_1(y,h)$. Then
>   $N_k=\frac{(-\log2)^{k+1}}{(k+1)!}+T_k+U_k$, where
>   T_k=\sum_{m=0}^{k-2}m!\,[y^mh^{k-2-m}]\mathcal A(y,h),\qquad |U_k|\le C\bigl(2^k+k\,(k-2)!\,6^{-k}\bigr).
> (line 475 and lines 484, 489, Lemma lem:pip:hankel, case h = 0) $S_0^\pm(y,h)=(y_0-y)^{-1-h}$ for
>   $y_0=\pm i\pi$, and $R_1=\mathcal A-S_0^+-S_0^-$ ... |[y^m]R_1(\cdot,0)|\le C(2\pi)^{-m}.
> (lines 514–520, Lemma lem:pip:TkT1)
>   T_k^{(1)}=2\operatorname{Re}\sum_{m=0}^{k-2}[h^{k-2-m}]\Bigl(y_0^{-1-m}\eu^{-h\Log y_0}(1+h)_m\Bigr).
>   $T_k=T_k^{(1)}+O\bigl((k-2)!\,\pi^{-(k-1)}\,k^{-2}\log k\bigr)$ as $k\to\infty$.
> (lines 541–544, Lemma lem:pip:twoterm) Let $y_0=i\pi$, $X=\psi(k-2)+\gamma-\log\pi$. As $k\to\infty$,
>   N_k=2\operatorname{Re}\Bigl[(k-2)!\,y_0^{1-k}+(k-3)!\,y_0^{2-k}\bigl(X-\tfrac{i\pi}2\bigr)\Bigr]
>   +O\Bigl((k-2)!\,\pi^{-(k-1)}\frac{\log^2k}{k^2}\Bigr).
> (lines 582–584, Lemma prop:pip:div) For every real $p>0$ the series $\sum_{k\ge0}N_kp^{-k-1}$
>   diverges. The series $\sum N_kh^{k+1}$ has radius of convergence $0$.

Conventions.
* `ζ(3)` is `X02.z3 = (riemannZeta 3).re` and `ζ(5)` is `z5 = (riemannZeta 5).re`; more generally
  `ζ(j+2)` is `(riemannZeta (j+2)).re`. These are the real numbers `∑_{n≥1} n^{-s}`
  (`zeta_nat_eq_tsum_of_gt_one`).
* `h = 1/p`, and `\hat E(h) = E(p)/4` is `Ehat h = X02.E (1/h) / 4`.
* The Taylor coefficients `w_{ij}` of `W` at `(0,0)` are defined as mixed iterated derivatives
  `∂_x^i ∂_y^j W(0,0)/(i! j!)` (`iteratedDeriv`, which is local; `W` is smooth near `(0,0)`, since
  there `1 - x ≠ 0` and all `rpow` bases are positive). Likewise `[y^m h^i]𝓐` for `𝓐(y,h)`, which is
  smooth on all of `ℝ²`. `W_m`, `q_n` and `N_k` are then defined literally by the supplement's
  formulas (`N_k` by `eq:pip:Nk`, with the sum over `1 ≤ n ≤ k-1`, empty for `k = 0, 1`). So
  `stmt_8` (expansion) and `stmt_9`–`stmt_11` (values) have content: they concern this `N`.
* `ψ = Γ'/Γ` is `psi = logDeriv Real.Gamma` (real log-derivative, `deriv Γ x / Γ x`); at the
  points `k - 2 ≥ 1` used here `Γ` is differentiable and positive. `γ` is
  `Real.eulerMascheroniConstant`.
* The generalised binomial coefficient `\binom hn` is `(descPochhammer ℝ n).eval h / n!`, and
  `(a)_m` is `(ascPochhammer ℂ m).eval a`; `[h^i]` of an entire function `f : ℂ → ℂ` is
  `iteratedDeriv i f 0 / i!`. `Log` is `Complex.log` (principal branch, `Log(iπ) = log π + iπ/2`).
  Negative powers `y_0^{1-k}` are written `((I π)^(k-1))⁻¹`.
* Every `O(·)` is an explicit bound `∃ C, ∃ k₀ (or h₀ > 0, p₀ > 1), ∀ …, |error| ≤ C · (scale)`.
  The paper prints no constants. In the growth laws the paper's relative error
  `N_k = σ A (main + O(e_k))` with `A > 0` is written `|N_k - σ A·main| ≤ C·A·e_k`.
* All integrals over `(0,∞)` that appear have their integrability asserted (`stmt_1`, `stmt_6`,
  `stmt_12`), so no statement is true only because a non-integrable integral is `0`.
* Natural-number subtractions `k - 1`, `k - 2`, `k - 3`, `(k-1)/2`, `(k-2)/2` occur only for
  `k ≥ k₀` with `k₀` chosen by the prover; `(k:ℝ) - 2` is a real subtraction.
* Divergence (`Q06.stmt_6`) is stated for every real `p ≠ 0` (the supplement's lemma says `p > 0`,
  the theorem "every `p`"; for `p < 0` the terms have the same absolute value). "Diverges" is
  stated in the strong form "the terms do not tend to `0`", and radius `0` as "for every `x ≠ 0`
  the terms `N_k x^{k+1}` are unbounded".

NOT COVERED:
* Holomorphy of `W` on the bidisc `|x| < 1/2, |y| < π/2` (only its consequence, the bound of
  `stmt_4`); the bound on the degrees of `W_m` is implicit in the definition.
* The Euler sum of Lemma `lem:pip:J`(ii) and the values `J_{1,0}, J_{1,1}, J_{2,0}` (proved in
  `X02`), and the Beta-function identities inside the proof of Lemma `lem:pip:C45` (only its
  conclusion is stated), Remark `rem:pip:C45x`, and the machine reductions of `J_{3,0}, …, J_{3,1}`
  in `Li_4(1/2)`, `Li_5(1/2)`.
* Watson's lemma uniformly in `n` (Lemma `lem:pip:watson`, the functions `D_m` and formula
  `eq:pip:Nk2`); the local form of `𝓐` at `±iπ` with complex `h` (Lemma `lem:pip:local`) and the
  Hankel bound for complex `h` (Lemma `lem:pip:hankel`, only its `h = 0` case is `stmt_14`).
* Lemmas lem:pip:Tk, lem:pip:TkT1, lem:pip:twoterm, the growth law and the divergence were audited as stmt_13, stmt_15 to stmt_20 and moved, unproved, to Q06_PiPGrowth.lean (planner split, 2026-11-26).
* Remark `rem:pip:duality` (`p → 1` via `π_p = π_q`) and all numerical remarks (signs of period
  four for `k ≥ 9`, optimal truncation, sharpness of the error terms).
-/
-- BEGIN STATEMENTS
namespace Q02

/-- `ζ(5)` as a real number: the real part of Mathlib's `riemannZeta 5` (which is real). -/
noncomputable def z5 : ℝ := (riemannZeta 5).re

/-- `J_{r,j} = ∫_0^∞ s^j 𝓛(s)^r ds`, `𝓛(s) = log(1+e^{-s})` (sec_pip.tex line 46). -/
noncomputable def J (r j : ℕ) : ℝ :=
  ∫ s in Set.Ioi (0 : ℝ), s ^ j * X02.calL s ^ r

/-- The Taylor coefficient `w_{ij}` of `W(x,y)` at `(0,0)`:
`∂_x^i ∂_y^j W(0,0) / (i! j!)` (sec_pip.tex line 85). -/
noncomputable def wcoef (i j : ℕ) : ℝ :=
  iteratedDeriv i (fun x => iteratedDeriv j (fun y => X02.W x y) 0) 0 /
    ((i.factorial : ℝ) * (j.factorial : ℝ))

/-- `W_m(s) = Σ_{i+j=m} w_{ij} s^j` (sec_pip.tex line 86). -/
noncomputable def Wpoly (m : ℕ) (s : ℝ) : ℝ :=
  ∑ j ∈ Finset.range (m + 1), wcoef (m - j) j * s ^ j

/-- `q_n = Σ_{m+r=n, r≥1} (1/r!) ∫_0^∞ W_m(s) 𝓛(s)^r ds` (sec_pip.tex line 106). -/
noncomputable def q (n : ℕ) : ℝ :=
  ∑ r ∈ Finset.Icc 1 n,
    (1 / (r.factorial : ℝ)) * ∫ s in Set.Ioi (0 : ℝ), Wpoly (n - r) s * X02.calL s ^ r

/-- `N_k = (-log 2)^{k+1}/(k+1)! + Σ_{n=1}^{k-1} q_n` (sec_pip.tex eq. `eq:pip:Nk`, line 112). -/
noncomputable def N (k : ℕ) : ℝ :=
  (-Real.log 2) ^ (k + 1) / ((k + 1).factorial : ℝ) + ∑ n ∈ Finset.Ico 1 k, q n

/-- `\hat E(h) = E(p)/4` with `p = 1/h` (sec_pip.tex line 70). -/
noncomputable def Ehat (h : ℝ) : ℝ := X02.E (1 / h) / 4

/-- The generalised binomial coefficient `\binom hn = h(h-1)⋯(h-n+1)/n!`. -/
noncomputable def binom (h : ℝ) (n : ℕ) : ℝ :=
  (descPochhammer ℝ n).eval h / (n.factorial : ℝ)

/-- `H_n(y,h) = e^{(n+1)y} (1+e^y)^{-1-h}` (sec_pip.tex line 380). -/
noncomputable def Hn (n : ℕ) (y h : ℝ) : ℝ :=
  Real.exp (((n : ℝ) + 1) * y) * (1 + Real.exp y) ^ (-1 - h)

/-- `𝓐(y,h) = e^{2y} (1+e^y)^{-1-h} = H_1(y,h)` (sec_pip.tex line 433). -/
noncomputable def Acal (y h : ℝ) : ℝ :=
  Real.exp (2 * y) * (1 + Real.exp y) ^ (-1 - h)

/-- The Taylor coefficient `[y^m h^i] 𝓐(y,h)` at `(0,0)`: `∂_h^i ∂_y^m 𝓐(0,0)/(m! i!)`. -/
noncomputable def Acoef (m i : ℕ) : ℝ :=
  iteratedDeriv i (fun h => iteratedDeriv m (fun y => Acal y h) 0) 0 /
    ((m.factorial : ℝ) * (i.factorial : ℝ))

/-- `T_k = Σ_{m=0}^{k-2} m! [y^m h^{k-2-m}] 𝓐(y,h)` (sec_pip.tex line 435). -/
noncomputable def T (k : ℕ) : ℝ :=
  ∑ m ∈ Finset.range (k - 1), (m.factorial : ℝ) * Acoef m (k - 2 - m)

/-- `T_k^{(1)} = 2 Re Σ_{m=0}^{k-2} [h^{k-2-m}] (y_0^{-1-m} e^{-h Log y_0} (1+h)_m)`, `y_0 = iπ`
(sec_pip.tex line 514). -/
noncomputable def T1 (k : ℕ) : ℝ :=
  2 * (∑ m ∈ Finset.range (k - 1),
    iteratedDeriv (k - 2 - m)
        (fun h : ℂ => ((Complex.I * Real.pi) ^ (m + 1))⁻¹ *
          Complex.exp (-(h * Complex.log (Complex.I * Real.pi))) *
          (ascPochhammer ℂ m).eval (1 + h)) 0 /
      ((k - 2 - m).factorial : ℂ)).re

/-- The digamma function `ψ = Γ'/Γ`, as the real log-derivative of `Real.Gamma`. -/
noncomputable def psi : ℝ → ℝ := logDeriv Real.Gamma

/-- sec_pip.tex line 48 (and line 106): the integrals `J_{r,j}` (`r ≥ 1`) converge, and so do the
integrals `∫_0^∞ W_m 𝓛^r` defining `q_n`. -/
def stmt_1 : Prop :=
  ∀ r j m : ℕ, 1 ≤ r →
    MeasureTheory.IntegrableOn (fun s => s ^ j * X02.calL s ^ r) (Set.Ioi 0) ∧
    MeasureTheory.IntegrableOn (fun s => Wpoly m s * X02.calL s ^ r) (Set.Ioi 0)

/-- Lemma `lem:pip:J`(i): `J_{1,j} = j! (1 - 2^{-j-1}) ζ(j+2)` for every `j`, and in particular
`J_{1,2} = 7π⁴/360`, `J_{1,3} = 45ζ(5)/8`. -/
def stmt_2 : Prop :=
  (∀ j : ℕ, J 1 j = (j.factorial : ℝ) * (1 - 1 / 2 ^ (j + 1)) *
      (riemannZeta ((j : ℂ) + 2)).re) ∧
  J 1 2 = 7 * Real.pi ^ 4 / 360 ∧ J 1 3 = 45 * z5 / 8

/-- Lemma `lem:pip:C45`: `C_4 = J_{3,0}/12 + J_{2,1}/8 = π⁴/2880` and
`C_5 = (J_{4,0} + 2 J_{3,1})/48 = -ζ(3)π²/48 + 31ζ(5)/128`. -/
def stmt_3 : Prop :=
  J 3 0 / 12 + J 2 1 / 8 = Real.pi ^ 4 / 2880 ∧
  (J 4 0 + 2 * J 3 1) / 48 = -(X02.z3 * Real.pi ^ 2 / 48) + 31 * z5 / 128

/-- Lemma `lem:pip:W` for every `K`: `|W(h,sh) - Σ_{m≤K} h^m W_m(s)| ≤ C_K (h(1+s))^{K+1}` for
`0 < h ≤ h_K`, `0 ≤ s ≤ h^{-1/2}`; and `0 < W(h,sh) ≤ 1` for `0 < h < 1`, `s ≥ 0`. -/
def stmt_4 : Prop :=
  (∀ K : ℕ, ∃ C h₀ : ℝ, 0 < h₀ ∧ ∀ h s : ℝ, 0 < h → h ≤ h₀ → 0 ≤ s → s ≤ 1 / Real.sqrt h →
    |X02.W h (s * h) - ∑ m ∈ Finset.range (K + 1), h ^ m * Wpoly m s|
      ≤ C * (h * (1 + s)) ^ (K + 1)) ∧
  (∀ h s : ℝ, 0 < h → h < 1 → 0 ≤ s → 0 < X02.W h (s * h) ∧ X02.W h (s * h) ≤ 1)

/-- sec_pip.tex lines 178 and 203–204: `W_0 = 1/2`, `W_1(s) = s/4 - log2/2`,
`W_2(s) = log²2/4 - (log2/4) s`, `W_3(s) = -log³2/12 + (log²2/8) s - (3/16) s² - (1/48) s³`. -/
def stmt_5 : Prop :=
  ∀ s : ℝ,
    Wpoly 0 s = 1 / 2 ∧
    Wpoly 1 s = s / 4 - Real.log 2 / 2 ∧
    Wpoly 2 s = Real.log 2 ^ 2 / 4 - Real.log 2 / 4 * s ∧
    Wpoly 3 s = -(Real.log 2 ^ 3 / 12) + Real.log 2 ^ 2 / 8 * s - 3 / 16 * s ^ 2 - 1 / 48 * s ^ 3

/-- sec_pip.tex line 118 (from eq. `eq:pip:Ehat`): for `0 < h < 1` the integrand of `Q(h)` is
integrable and `\hat E(h) = h²/(1-h) · Q(h)`. -/
def stmt_6 : Prop :=
  ∀ h : ℝ, 0 < h → h < 1 →
    MeasureTheory.IntegrableOn (fun s => X02.W h (s * h) * ((1 + Real.exp (-s)) ^ h - 1))
      (Set.Ioi 0) ∧
    Ehat h = h ^ 2 / (1 - h) * X02.Q h

/-- Lemma `lem:pip:expansion`, the core: for every `K ≥ 1`,
`Q(h) = Σ_{n=1}^K q_n h^n + O(h^{K+1})` as `h → 0⁺`. -/
def stmt_7 : Prop :=
  ∀ K : ℕ, 1 ≤ K → ∃ C h₀ : ℝ, 0 < h₀ ∧ ∀ h : ℝ, 0 < h → h ≤ h₀ →
    |X02.Q h - ∑ n ∈ Finset.Icc 1 K, q n * h ^ n| ≤ C * h ^ (K + 1)

/-- thm:pip, "the expansion continues to every order, `π_p/4 ∼ 1 + Σ_{k≥0} N_k p^{-k-1}`":
for every `K` there are `C`, `p₀` with
`|π_p/4 - 1 - Σ_{k=0}^K N_k p^{-k-1}| ≤ C p^{-K-2}` for `p ≥ p₀`. -/
def stmt_8 : Prop :=
  ∀ K : ℕ, ∃ C p₀ : ℝ, 1 < p₀ ∧ ∀ p : ℝ, p₀ ≤ p →
    |X02.piP p / 4 - 1 - ∑ k ∈ Finset.range (K + 1), N k * (1 / p) ^ (k + 1)|
      ≤ C * (1 / p) ^ (K + 2)

/-- Lemma `prop:pip:first`: the `N_k` of `eq:pip:Nk` agree with the printed `N_0, …, N_3`. -/
def stmt_9 : Prop :=
  N 0 = X02.N0 ∧ N 1 = X02.N1 ∧ N 2 = X02.N2 ∧ N 3 = X02.N3

/-- The complement after thm:pip (main.tex lines 1374–1375): the printed value of `N_4`. -/
def stmt_10 : Prop :=
  N 4 = -(Real.log 2 ^ 5 / 120)
    + Real.pi ^ 2 / 24 * (1 - Real.log 2 + Real.log 2 ^ 2 / 2)
    + X02.z3 / 4 * (1 - Real.log 2) + Real.pi ^ 4 / 2880

/-- The complement after thm:pip (main.tex lines 1376–1378): the printed value of `N_5`. -/
def stmt_11 : Prop :=
  N 5 = Real.log 2 ^ 6 / 720
    + Real.pi ^ 2 / 24 * (1 - Real.log 2 + Real.log 2 ^ 2 / 2 - Real.log 2 ^ 3 / 6)
    + X02.z3 / 4 * (1 - Real.log 2 + Real.log 2 ^ 2 / 2)
    + Real.pi ^ 4 / 2880 * (1 - Real.log 2)
    - 7 * Real.pi ^ 4 / 1920 - X02.z3 * Real.pi ^ 2 / 48 + z5 / 8

/-- Lemma `lem:pip:laplace` (also main.tex line 1381): for `0 < h < 1` the integrals
`∫_0^∞ e^{-ny/h} H_n(y,h) dy` converge and
`\hat E(h) = h Σ_{n≥1} \binom hn ∫_0^∞ e^{-ny/h} H_n(y,h) dy`. -/
def stmt_12 : Prop :=
  ∀ h : ℝ, 0 < h → h < 1 →
    (∀ n : ℕ, 1 ≤ n →
      MeasureTheory.IntegrableOn (fun y => Real.exp (-((n : ℝ) * y / h)) * Hn n y h)
        (Set.Ioi 0)) ∧
    HasSum (fun n : ℕ => h * binom h (n + 1) *
        ∫ y in Set.Ioi (0 : ℝ), Real.exp (-(((n + 1 : ℕ) : ℝ) * y / h)) * Hn (n + 1) y h)
      (Ehat h)

/-- Lemma `lem:pip:hankel`, case `h = 0` (the simple poles at `±iπ`): with
`S_0^+ + S_0^- = 1/(iπ - y) + 1/(-iπ - y) = -2y/(y² + π²)`, the Taylor coefficients of
`R_1(y,0) = e^{2y}/(1+e^y) + 2y/(y²+π²)` at `0` satisfy `|[y^m]R_1(·,0)| ≤ C (2π)^{-m}`. -/
def stmt_14 : Prop :=
  ∃ C : ℝ, ∀ m : ℕ,
    |iteratedDeriv m (fun y : ℝ => Acal y 0 + 2 * y / (y ^ 2 + Real.pi ^ 2)) 0 /
        (m.factorial : ℝ)| ≤ C / (2 * Real.pi) ^ m

end Q02
-- END STATEMENTS

open MeasureTheory Set Real

namespace Q02

lemma intOn_pow_calL (r j : ℕ) (hr : 1 ≤ r) :
    IntegrableOn (fun s => s ^ j * X02.calL s ^ r) (Ioi 0) := by
  refine X02.intOn_of_bound ((continuous_pow j).mul (X02.continuous_calL.pow r))
    (X02.intOn_pow_exp j one_pos) (fun s hs => ?_)
  have h0 := X02.calL_nonneg s
  have h1 := X02.calL_le s
  have h2 : Real.exp (-s) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
  have hp : X02.calL s ^ r ≤ X02.calL s := pow_le_of_le_one h0 (h1.trans h2) (by omega)
  rw [abs_of_nonneg (by positivity), one_mul]
  exact mul_le_mul_of_nonneg_left (hp.trans h1) (pow_nonneg hs.le j)

lemma intOn_Wpoly_calL (m r : ℕ) (hr : 1 ≤ r) :
    IntegrableOn (fun s => Wpoly m s * X02.calL s ^ r) (Ioi 0) := by
  have : (fun s => Wpoly m s * X02.calL s ^ r) =
      fun s => ∑ j ∈ Finset.range (m + 1), wcoef (m - j) j * (s ^ j * X02.calL s ^ r) := by
    funext s; rw [Wpoly, Finset.sum_mul]; congr 1; funext j; ring
  rw [this]
  exact integrable_finsetSum _ (fun j _ => (intOn_pow_calL r j hr).const_mul _)

lemma hasSum_zeta_nat (k : ℕ) (hk : 2 ≤ k) :
    HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ k) (riemannZeta (k : ℂ)).re := by
  have h1 : riemannZeta (k : ℂ) = ((∑' n : ℕ, 1 / ((n : ℝ)) ^ k : ℝ) : ℂ) := by
    rw [zeta_nat_eq_tsum_of_gt_one (by omega), Complex.ofReal_tsum]
    push_cast; rfl
  have hs : Summable (fun n : ℕ => 1 / ((n : ℝ)) ^ k) :=
    (Real.summable_one_div_nat_pow (p := k)).2 (by omega)
  rw [h1, Complex.ofReal_re, hs.tsum_eq_zero_add]
  have hk0 : k ≠ 0 := by omega
  simp only [CharP.cast_eq_zero, ne_eq, hk0, not_false_eq_true, zero_pow, div_zero, zero_add]
  convert ((summable_nat_add_iff 1).2 hs).hasSum using 1
  funext n; push_cast; rfl

lemma J1_eq (j : ℕ) : J 1 j = (j.factorial : ℝ) * (1 - 1 / 2 ^ (j + 1)) *
    (riemannZeta ((j : ℂ) + 2)).re := by
  have h := X02.integral_series (fun n => (-1) ^ n) (fun n => by simp) j
  have e : J 1 j = ∫ s in Ioi (0 : ℝ), s ^ j * ∑' n : ℕ,
      (-1) ^ n * Real.exp (-(((n : ℝ) + 1) * s)) / ((n : ℝ) + 1) :=
    setIntegral_congr_fun measurableSet_Ioi (fun s hs => by rw [pow_one, X02.calL_series hs])
  have hz := X02.hasSum_alt (j + 2) (hasSum_zeta_nat (j + 2) (by omega))
  rw [e, h]
  have : (∑' n : ℕ, (-1 : ℝ) ^ n / ((n : ℝ) + 1) ^ (j + 2))
      = ∑' n : ℕ, (-1) ^ n * (1 / ((n : ℝ) + 1) ^ (j + 2)) := by
    congr 1; funext n; ring
  rw [this, hz.tsum_eq]
  push_cast
  rw [pow_succ, pow_succ]
  field_simp

end Q02

theorem Q02.proof_1 : Q02.stmt_1 := fun r j m hr =>
  ⟨Q02.intOn_pow_calL r j hr, Q02.intOn_Wpoly_calL m r hr⟩

theorem Q02.proof_2 : Q02.stmt_2 := by
  refine ⟨Q02.J1_eq, ?_, ?_⟩
  · rw [Q02.J1_eq]
    have : ((2 : ℕ) : ℂ) + 2 = 4 := by norm_num
    rw [this, riemannZeta_four]
    have : ((Real.pi : ℂ) ^ 4 / 90).re = Real.pi ^ 4 / 90 := by
      rw [show (Real.pi : ℂ) ^ 4 / 90 = ((Real.pi ^ 4 / 90 : ℝ) : ℂ) by push_cast; ring]
      exact Complex.ofReal_re _
    rw [this]; norm_num [Nat.factorial]; ring
  · rw [Q02.J1_eq]
    have : ((3 : ℕ) : ℂ) + 2 = 5 := by norm_num
    rw [this, ← Q02.z5]; norm_num [Nat.factorial]; ring

namespace Q02

/-- `tm x k c u = e^{-(x+k)u} (1+e^{-u})^{-c}`. -/
noncomputable def tm (x k c u : ℝ) : ℝ := Real.exp (-((x + k) * u)) * (1 + Real.exp (-u)) ^ (-c)

lemma tm_hasDerivAt (x k c u : ℝ) :
    HasDerivAt (tm x k c) (-(x + k) * tm x k c u + c * tm x (k + 1) (c + 1) u) u := by
  have hP : 0 < 1 + Real.exp (-u) := by positivity
  have h1 : HasDerivAt (fun u => Real.exp (-((x + k) * u)))
      (Real.exp (-((x + k) * u)) * (-((x + k) * 1))) u :=
    ((hasDerivAt_id u).const_mul (x + k)).neg.exp
  have h2 : HasDerivAt (fun u => 1 + Real.exp (-u)) (Real.exp (-u) * (-1)) u :=
    ((hasDerivAt_id u).neg.exp).const_add 1
  have h3 := h2.rpow_const (p := -c) (Or.inl hP.ne')
  have := h1.mul h3
  unfold tm
  convert this using 1
  rw [show -(c + 1) = -c - 1 by ring, show -((x + (k + 1)) * u) = -((x + k) * u) + -u by ring,
    Real.exp_add]
  ring

/-- `e x k c = d/du tm x k c`. -/
noncomputable def ed (x k c u : ℝ) : ℝ := -(x + k) * tm x k c u + c * tm x (k + 1) (c + 1) u

lemma ed_hasDerivAt (x k c u : ℝ) :
    HasDerivAt (ed x k c) (-(x + k) * ed x k c u + c * ed x (k + 1) (c + 1) u) u := by
  have := ((tm_hasDerivAt x k c u).const_mul (-(x + k))).add
    ((tm_hasDerivAt x (k + 1) (c + 1) u).const_mul c)
  convert this using 1
  · funext v; rfl
  · unfold ed; ring

noncomputable def d2 (x u : ℝ) : ℝ := -(x + 0) * ed x 0 (1 + x) u + (1 + x) * ed x (0 + 1) (1 + x + 1) u

noncomputable def d3 (x u : ℝ) : ℝ :=
  -(x + 0) * d2 x u + (1 + x) * (-(x + (0 + 1)) * ed x (0 + 1) (1 + x + 1) u +
    (1 + x + 1) * ed x (0 + 1 + 1) (1 + x + 1 + 1) u)

lemma d2_hasDerivAt (x u : ℝ) : HasDerivAt (d2 x) (d3 x u) u := by
  have := ((ed_hasDerivAt x 0 (1 + x) u).const_mul (-(x + 0))).add
    ((ed_hasDerivAt x (0 + 1) (1 + x + 1) u).const_mul (1 + x))
  convert this using 1
  · funext v; rfl
  · unfold d3 d2; ring

lemma deriv_comp_scale {f f' : ℝ → ℝ} (hf : ∀ u, HasDerivAt f (f' u) u) (a c : ℝ) :
    deriv (fun y => c * f (a * y)) = fun y => c * a * f' (a * y) := by
  funext y
  have := ((hf (a * y)).comp y ((hasDerivAt_id y).const_mul a)).const_mul c
  exact this.deriv.trans (by ring)

lemma W_eq (x : ℝ) : (fun y => X02.W x y) = fun y => 1 * tm x 0 (1 + x) (1 / (1 - x) * y) := by
  funext y; unfold X02.W tm
  rw [show -(1 + x) = -1 - x by ring, show x * y / (1 - x) = (x + 0) * (1 / (1 - x) * y) by ring,
    show y / (1 - x) = 1 / (1 - x) * y by ring, one_mul]

lemma iD1 (x : ℝ) : iteratedDeriv 1 (fun y => X02.W x y) =
    fun y => 1 * (1 / (1 - x)) * ed x 0 (1 + x) (1 / (1 - x) * y) := by
  rw [iteratedDeriv_one, W_eq, deriv_comp_scale (tm_hasDerivAt x 0 (1 + x))]; rfl

lemma iD2 (x : ℝ) : iteratedDeriv 2 (fun y => X02.W x y) =
    fun y => 1 * (1 / (1 - x)) * (1 / (1 - x)) * d2 x (1 / (1 - x) * y) := by
  rw [iteratedDeriv_succ, iD1, deriv_comp_scale (ed_hasDerivAt x 0 (1 + x))]; rfl

lemma iD3 (x : ℝ) : iteratedDeriv 3 (fun y => X02.W x y) =
    fun y => 1 * (1 / (1 - x)) * (1 / (1 - x)) * (1 / (1 - x)) * d3 x (1 / (1 - x) * y) := by
  rw [iteratedDeriv_succ, iD2, deriv_comp_scale (d2_hasDerivAt x)]

/-- `E(x) = e^{-x log 2}`. -/
lemma two_rpow (c : ℝ) : (2 : ℝ) ^ (-c) = Real.exp (-c * Real.log 2) := by
  rw [Real.rpow_def_of_pos (by norm_num)]; ring_nf

lemma tm_zero (x k c : ℝ) : tm x k c 0 = Real.exp (-c * Real.log 2) := by
  unfold tm; simp only [mul_zero, neg_zero, Real.exp_zero, one_mul]
  rw [show (1 : ℝ) + 1 = 2 by norm_num, two_rpow]

lemma expadd (c : ℝ) : Real.exp (-(c + 1) * Real.log 2) = Real.exp (-c * Real.log 2) / 2 := by
  rw [show -(c + 1) * Real.log 2 = -c * Real.log 2 + -Real.log 2 by ring, Real.exp_add,
    Real.exp_neg, Real.exp_log (by norm_num)]; ring

lemma expbase (x : ℝ) : Real.exp (-(1 + x) * Real.log 2) = (1 / 2) * Real.exp (-Real.log 2 * x) := by
  rw [show -(1 + x) * Real.log 2 = -Real.log 2 * x + -Real.log 2 by ring, Real.exp_add,
    Real.exp_neg, Real.exp_log (by norm_num)]; ring

lemma F0 : (fun x => iteratedDeriv 0 (fun y => X02.W x y) 0) =
    fun x => (1 / 2) * Real.exp (-Real.log 2 * x) := by
  funext x
  rw [iteratedDeriv_zero, show X02.W x 0 = 1 * tm x 0 (1 + x) (1 / (1 - x) * 0) from
    congrFun (W_eq x) 0, mul_zero, one_mul, tm_zero, expbase]

lemma F1 : (fun x => iteratedDeriv 1 (fun y => X02.W x y) 0) =ᶠ[nhds 0]
    fun x => (1 / 4) * Real.exp (-Real.log 2 * x) := by
  filter_upwards [Iio_mem_nhds (show (0 : ℝ) < 1 by norm_num)] with x hx
  have hx1 : 1 - x ≠ 0 := by have := mem_Iio.mp hx; linarith
  rw [iD1]; simp only [mul_zero]
  unfold ed; simp only [tm_zero, expadd, expbase]
  field_simp; ring

lemma F2 : (fun x => iteratedDeriv 2 (fun y => X02.W x y) 0) =
    fun x => x * ((1 / 8) * Real.exp (-Real.log 2 * x) * (x - 3) * (1 - x)⁻¹ ^ 2) := by
  funext x
  rw [iD2]; simp only [mul_zero]
  unfold d2 ed; simp only [tm_zero, expadd, expbase]
  ring

lemma w00 : wcoef 0 0 = 1 / 2 := by
  have := congrFun F0 0
  unfold wcoef; rw [iteratedDeriv_zero, this]; simp

lemma wi0 (i : ℕ) : wcoef i 0 = (1 / 2) * (-Real.log 2) ^ i / (i.factorial : ℝ) := by
  unfold wcoef
  rw [show (fun x => iteratedDeriv 0 (fun y => X02.W x y) 0) = _ from F0,
    iteratedDeriv_const_mul_field, iteratedDeriv_exp_const_mul]
  simp

lemma wi1 (i : ℕ) : wcoef i 1 = (1 / 4) * (-Real.log 2) ^ i / (i.factorial : ℝ) := by
  unfold wcoef
  rw [F1.iteratedDeriv_eq, iteratedDeriv_const_mul_field, iteratedDeriv_exp_const_mul]
  simp

lemma w02 : wcoef 0 2 = 0 := by
  have := congrFun F2 0
  unfold wcoef; rw [iteratedDeriv_zero, this]; simp

lemma w12 : wcoef 1 2 = -3 / 16 := by
  unfold wcoef
  rw [iteratedDeriv_one, show (fun x => iteratedDeriv 2 (fun y => X02.W x y) 0) = _ from F2]
  have hnum : DifferentiableAt ℝ (fun x : ℝ => (1 / 8) * Real.exp (-Real.log 2 * x) * (x - 3)) 0 := by
    fun_prop
  have hden : DifferentiableAt ℝ (fun x : ℝ => 1 - x) 0 := by fun_prop
  have hφ := hnum.mul ((hden.inv (by norm_num)).pow 2)
  have := (hasDerivAt_id (0 : ℝ)).mul hφ.hasDerivAt
  have e : deriv (fun x : ℝ => x * ((1 / 8) * Real.exp (-Real.log 2 * x) * (x - 3) *
      (1 - x)⁻¹ ^ 2)) 0 = _ := this.deriv
  rw [e]; simp; norm_num

lemma w03 : wcoef 0 3 = -1 / 48 := by
  unfold wcoef
  rw [iteratedDeriv_zero, iD3]; simp only [mul_zero]
  unfold d3 d2 ed; simp only [tm_zero, expadd, expbase]
  norm_num

end Q02

theorem Q02.proof_5 : Q02.stmt_5 := by
  intro s
  simp only [Q02.Wpoly, Finset.sum_range_succ, Finset.sum_range_zero, zero_add, Nat.sub_self,
    Nat.sub_zero, Q02.wi0, Q02.wi1, Q02.w02, Q02.w12, Q02.w03, Nat.reduceSub]
  norm_num [Nat.factorial]
  refine ⟨by ring, by ring, by ring⟩

namespace Q02

lemma q1_eq : q 1 = Real.pi ^ 2 / 24 := by
  unfold q
  simp only [Finset.Icc_self, Finset.sum_singleton, Nat.sub_self, Nat.factorial_one]
  have : (fun s => Wpoly 0 s * X02.calL s ^ 1) = fun s => (1 / 2) * X02.calL s := by
    funext s; rw [(proof_5 s).1, pow_one]
  rw [this, integral_const_mul, X02.J10]; push_cast; ring

lemma q2_eq : q 2 = X02.z3 / 4 - Real.pi ^ 2 * Real.log 2 / 24 := by
  unfold q
  rw [show Finset.Icc 1 2 = {1, 2} by rfl, Finset.sum_pair (by norm_num)]
  simp only [Nat.reduceSub, Nat.factorial_one, Nat.factorial_two]
  have e1 : (fun s => Wpoly 1 s * X02.calL s ^ 1) =
      fun s => (1 / 4) * (s * X02.calL s) + (-(Real.log 2 / 2)) * X02.calL s := by
    funext s; rw [(proof_5 s).2.1, pow_one]; ring
  have e2 : (fun s => Wpoly 0 s * X02.calL s ^ 2) = fun s => (1 / 2) * X02.calL s ^ 2 := by
    funext s; rw [(proof_5 s).1]
  rw [e1, e2, integral_add (X02.intOn_s_calL.const_mul _) (X02.intOn_calL.const_mul _),
    integral_const_mul, integral_const_mul, integral_const_mul, X02.J10, X02.J11, X02.J20]
  push_cast; ring

end Q02

theorem Q02.proof_9 : Q02.stmt_9 := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · simp [Q02.N, X02.N0]
  · simp [Q02.N, X02.N1]; ring
  · simp only [Q02.N, X02.N2, show Finset.Ico 1 2 = {1} by rfl, Finset.sum_singleton, Q02.q1_eq]
    norm_num [Nat.factorial]; ring
  · simp only [Q02.N, X02.N3, show Finset.Ico 1 3 = {1, 2} by rfl,
      Finset.sum_pair (show (1 : ℕ) ≠ 2 by norm_num), Q02.q1_eq, Q02.q2_eq]
    norm_num [Nat.factorial]; ring

namespace Q02

lemma EInt_eq {h : ℝ} (hh : 0 < h) (h1 : h < 1) (s : ℝ) :
    X02.EIntegrand (1 / h) s = X02.W h (s * h) * ((1 + Real.exp (-s)) ^ h - 1) := by
  unfold X02.EIntegrand X02.W
  have h1h : 1 - h ≠ 0 := by linarith
  have a1 : s / (1 / h * (1 / h - 1)) = h * (s * h) / (1 - h) := by
    field_simp
  have a2 : s / (1 / h - 1) = s * h / (1 - h) := by
    field_simp
  rw [a1, a2, one_div_one_div]

end Q02

theorem Q02.proof_6 : Q02.stmt_6 := by
  intro h hh h1
  have hp : 1 < 1 / h := by rw [lt_one_div one_pos hh]; simpa using h1
  refine ⟨?_, ?_⟩
  · have := (X02.proof_1 (1 / h) hp).2.1
    exact this.congr_fun (fun s _ => Q02.EInt_eq hh h1 s) measurableSet_Ioi
  · unfold Q02.Ehat
    rw [X02.E_eq_Q hp, one_div_one_div]
    have h1h : 1 - h ≠ 0 := by linarith
    field_simp

namespace Q02

lemma intW (m r : ℕ) (hr : 1 ≤ r) :
    ∫ s in Ioi (0 : ℝ), Wpoly m s * X02.calL s ^ r =
      ∑ j ∈ Finset.range (m + 1), wcoef (m - j) j * J r j := by
  have : (fun s => Wpoly m s * X02.calL s ^ r) =
      fun s => ∑ j ∈ Finset.range (m + 1), wcoef (m - j) j * (s ^ j * X02.calL s ^ r) := by
    funext s; rw [Wpoly, Finset.sum_mul]; congr 1; funext j; ring
  rw [this, integral_finsetSum _ (fun j _ => (intOn_pow_calL r j hr).const_mul _)]
  congr 1; funext j; rw [integral_const_mul]; rfl

lemma J10' : J 1 0 = Real.pi ^ 2 / 12 := by
  unfold J; simp only [pow_zero, one_mul, pow_one]; exact X02.J10
lemma J11' : J 1 1 = 3 * X02.z3 / 4 := by
  unfold J; simp only [pow_one]; exact X02.J11
lemma J20' : J 2 0 = X02.z3 / 4 := by
  unfold J; simp only [pow_zero, one_mul]; exact X02.J20

lemma q_eq (n : ℕ) : q n = ∑ r ∈ Finset.Icc 1 n, 1 / (r.factorial : ℝ) *
    ∑ j ∈ Finset.range (n - r + 1), wcoef (n - r - j) j * J r j := by
  unfold q
  refine Finset.sum_congr rfl (fun r hr => ?_)
  rw [intW _ _ (Finset.mem_Icc.mp hr).1]

lemma q3_eq : q 3 = Real.pi ^ 2 * Real.log 2 ^ 2 / 48 - X02.z3 * Real.log 2 / 4 +
    (J 3 0 / 12 + J 2 1 / 8) := by
  rw [q_eq, show Finset.Icc 1 3 = {1, 2, 3} from rfl]
  simp only [Finset.mem_insert, Finset.mem_singleton, not_false_eq_true, Finset.sum_insert,
    Finset.sum_singleton, or_self, Nat.reduceEqDiff, Nat.reduceSub, Nat.sub_self,
    Finset.sum_range_succ, Finset.sum_range_zero, zero_add, Nat.reduceAdd, wi0, wi1, w02, w12,
    w03, J10', J11', J20']
  norm_num [Nat.factorial]
  ring

lemma q4_eq : q 4 = -(Real.pi ^ 2 * Real.log 2 ^ 3 / 144) + X02.z3 * Real.log 2 ^ 2 / 8
    - 3 / 16 * J 1 2 - 1 / 48 * J 1 3 - Real.log 2 * (J 3 0 / 12 + J 2 1 / 8)
    + (J 4 0 + 2 * J 3 1) / 48 := by
  rw [q_eq, show Finset.Icc 1 4 = {1, 2, 3, 4} from rfl]
  simp only [Finset.mem_insert, Finset.mem_singleton, not_false_eq_true, Finset.sum_insert,
    Finset.sum_singleton, or_self, Nat.reduceEqDiff, Nat.reduceSub, Nat.sub_self,
    Finset.sum_range_succ, Finset.sum_range_zero, zero_add, Nat.reduceAdd, wi0, wi1, w02, w12,
    w03, J10', J11', J20']
  norm_num [Nat.factorial]
  ring

lemma N4_of (hC : J 3 0 / 12 + J 2 1 / 8 = Real.pi ^ 4 / 2880) :
    N 4 = -(Real.log 2 ^ 5 / 120)
      + Real.pi ^ 2 / 24 * (1 - Real.log 2 + Real.log 2 ^ 2 / 2)
      + X02.z3 / 4 * (1 - Real.log 2) + Real.pi ^ 4 / 2880 := by
  rw [N, show Finset.Ico 1 4 = {1, 2, 3} from rfl]
  simp only [Finset.mem_insert, Finset.mem_singleton, not_false_eq_true, Finset.sum_insert,
    Finset.sum_singleton, or_self, Nat.reduceEqDiff, q1_eq, q2_eq, q3_eq, hC]
  norm_num [Nat.factorial]
  ring

lemma N5_of (hC : J 3 0 / 12 + J 2 1 / 8 = Real.pi ^ 4 / 2880)
    (hD : (J 4 0 + 2 * J 3 1) / 48 = -(X02.z3 * Real.pi ^ 2 / 48) + 31 * z5 / 128) :
    N 5 = Real.log 2 ^ 6 / 720
    + Real.pi ^ 2 / 24 * (1 - Real.log 2 + Real.log 2 ^ 2 / 2 - Real.log 2 ^ 3 / 6)
    + X02.z3 / 4 * (1 - Real.log 2 + Real.log 2 ^ 2 / 2)
    + Real.pi ^ 4 / 2880 * (1 - Real.log 2)
    - 7 * Real.pi ^ 4 / 1920 - X02.z3 * Real.pi ^ 2 / 48 + z5 / 8 := by
  rw [N, show Finset.Ico 1 5 = {1, 2, 3, 4} from rfl]
  simp only [Finset.mem_insert, Finset.mem_singleton, not_false_eq_true, Finset.sum_insert,
    Finset.sum_singleton, or_self, Nat.reduceEqDiff, q1_eq, q2_eq, q3_eq, q4_eq, hC, hD,
    proof_2.2.1, proof_2.2.2]
  norm_num [Nat.factorial]
  ring

end Q02

open scoped ContDiff

namespace Q02

/-! ### Lemma `lem:pip:W` for every `K`: smoothness of `W` and Taylor's theorem twice -/

/-- the open set `x < 1` where `W` is smooth -/
def Uw : Set (ℝ × ℝ) := Iio 1 ×ˢ univ

lemma isOpen_Uw : IsOpen Uw := isOpen_Iio.prod isOpen_univ

lemma W_contDiffAt {p : ℝ × ℝ} (hp : p.1 < 1) :
    ContDiffAt ℝ ∞ (Function.uncurry X02.W) p := by
  have h1 : (1 - p.1) ≠ 0 := by linarith
  have e : Function.uncurry X02.W =ᶠ[nhds p] fun q : ℝ × ℝ =>
      Real.exp (-(q.1 * q.2 / (1 - q.1))) *
        Real.exp (Real.log (1 + Real.exp (-(q.2 / (1 - q.1)))) * (-1 - q.1)) := by
    filter_upwards with q
    simp only [Function.uncurry, X02.W]
    rw [Real.rpow_def_of_pos (by positivity)]
  refine ContDiffAt.congr_of_eventuallyEq ?_ e
  have hd : ContDiffAt ℝ ∞ (fun q : ℝ × ℝ => 1 - q.1) p := by fun_prop
  have hq1 : ContDiffAt ℝ ∞ (fun q : ℝ × ℝ => q.1 * q.2 / (1 - q.1)) p :=
    (contDiffAt_fst.mul contDiffAt_snd).div hd h1
  have hq2 : ContDiffAt ℝ ∞ (fun q : ℝ × ℝ => q.2 / (1 - q.1)) p := contDiffAt_snd.div hd h1
  have hl : ContDiffAt ℝ ∞ (fun q : ℝ × ℝ => Real.log (1 + Real.exp (-(q.2 / (1 - q.1))))) p :=
    (contDiffAt_const.add hq2.neg.exp).log (by positivity)
  exact hq1.neg.exp.mul ((hl.mul (contDiffAt_const.sub contDiffAt_fst)).exp)

/-- `G j (x, y) = ∂_y^j W(x, y)` -/
noncomputable def G (j : ℕ) (p : ℝ × ℝ) : ℝ := iteratedDeriv j (fun y => X02.W p.1 y) p.2

lemma G_contDiffOn (j : ℕ) : ContDiffOn ℝ ∞ (G j) Uw := by
  induction j with
  | zero =>
    intro p hp
    have : G 0 = Function.uncurry X02.W := by funext q; simp [G, Function.uncurry]
    rw [this]
    exact (W_contDiffAt hp.1).contDiffWithinAt
  | succ j ih =>
    have hf : ContDiffOn ℝ (∞ + 1) (Function.uncurry fun x y => G j (x, y)) (Iio 1 ×ˢ univ) := by
      have : (Function.uncurry fun x y => G j (x, y)) = G j := by funext q; rfl
      rw [this]; exact ih
    have h2 := hf.fderiv_right (u := univ) isOpen_univ
    have h3 := h2.clm_apply (contDiffOn_const (c := (1 : ℝ)))
    refine h3.congr (fun p _ => ?_)
    simp only [G, iteratedDeriv_succ]
    rfl

lemma contDiffOn_le_top {f : ℝ → ℝ} {S : Set ℝ} (hf : ContDiffOn ℝ ∞ f S) (n : ℕ) :
    ContDiffOn ℝ n f S := hf.of_le (by exact_mod_cast le_top)

lemma contDiffAt_le_top {f : ℝ → ℝ} {x : ℝ} (hf : ContDiffAt ℝ ∞ f x) (n : ℕ) :
    ContDiffAt ℝ n f x := hf.of_le (by exact_mod_cast le_top)

/-- Taylor's theorem at `0` on `[0, 1/2]` for a function smooth on `x < 1`. -/
lemma taylor_x {f : ℝ → ℝ} (hf : ContDiffOn ℝ ∞ f (Iio 1)) (n : ℕ) :
    ∃ C, ∀ x ∈ Icc (0 : ℝ) (1 / 2),
      |f x - ∑ k ∈ Finset.range (n + 1), iteratedDeriv k f 0 / (k.factorial : ℝ) * x ^ k| ≤
        C * x ^ (n + 1) := by
  have hsub : Icc (0 : ℝ) (1 / 2) ⊆ Iio 1 := fun x hx => by simp only [mem_Iio]; linarith [hx.2]
  have hcd : ContDiffOn ℝ ((n + 1 : ℕ) : ℕ∞ω) f (Icc 0 (1 / 2)) :=
    (contDiffOn_le_top hf (n + 1)).mono hsub
  obtain ⟨C, hC⟩ := exists_taylor_mean_remainder_bound (n := n) (by exact_mod_cast hcd)
  refine ⟨C, fun x hx => ?_⟩
  have h := hC x hx
  rw [taylor_within_apply, Real.norm_eq_abs, sub_zero] at h
  convert h using 3
  refine Finset.sum_congr rfl (fun k _ => ?_)
  have hat : ContDiffAt ℝ ∞ f 0 := hf.contDiffAt (Iio_mem_nhds (by norm_num))
  rw [iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_Icc (by norm_num))
    (contDiffAt_le_top hat k) (left_mem_Icc.2 (by norm_num)), smul_eq_mul]
  ring

lemma Wy_contDiffAt {h : ℝ} (hh : h < 1) (y : ℝ) : ContDiffAt ℝ ∞ (fun y => X02.W h y) y := by
  have := (W_contDiffAt (p := (h, y)) hh).comp y (contDiffAt_const.prodMk contDiffAt_id)
  exact this

/-- Taylor's theorem in `y`, uniformly for `h ∈ [0, 1/2]`. -/
lemma taylor_y (K : ℕ) : ∃ C, ∀ h ∈ Icc (0 : ℝ) (1 / 2), ∀ y ∈ Icc (0 : ℝ) 1,
    |X02.W h y - ∑ j ∈ Finset.range (K + 1), G j (h, 0) / (j.factorial : ℝ) * y ^ j| ≤
      C * y ^ (K + 1) := by
  have hS : Icc (0 : ℝ) (1 / 2) ×ˢ Icc (0 : ℝ) 1 ⊆ Uw := fun p hp => by
    simp only [Uw, mem_prod, mem_Iio, mem_univ, and_true]; linarith [hp.1.2]
  have hcont : ContinuousOn (G (K + 1)) (Icc (0 : ℝ) (1 / 2) ×ˢ Icc (0 : ℝ) 1) :=
    ((G_contDiffOn (K + 1)).continuousOn).mono hS
  obtain ⟨M, hM⟩ := (isCompact_Icc.prod isCompact_Icc).exists_bound_of_continuousOn hcont
  refine ⟨M / (K.factorial : ℝ), fun h hh y hy => ?_⟩
  have hh1 : h < 1 := by linarith [hh.2]
  have hcd : ContDiffOn ℝ ((K + 1 : ℕ) : ℕ∞ω) (fun y => X02.W h y) (Icc 0 1) :=
    fun y' _ => (contDiffAt_le_top (Wy_contDiffAt hh1 y') (K + 1)).contDiffWithinAt
  have hb : ∀ y' ∈ Icc (0 : ℝ) 1,
      ‖iteratedDerivWithin (K + 1) (fun y => X02.W h y) (Icc 0 1) y'‖ ≤ M := by
    intro y' hy'
    rw [iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_Icc (by norm_num))
      (contDiffAt_le_top (Wy_contDiffAt hh1 y') (K + 1)) hy']
    exact hM (h, y') ⟨hh, hy'⟩
  have T := taylor_mean_remainder_bound (by exact_mod_cast hcd) hy hb
  rw [taylor_within_apply, Real.norm_eq_abs, sub_zero] at T
  have e : ∑ k ∈ Finset.range (K + 1), ((k.factorial : ℝ)⁻¹ * y ^ k) •
        iteratedDerivWithin k (fun y => X02.W h y) (Icc 0 1) 0 =
      ∑ j ∈ Finset.range (K + 1), G j (h, 0) / (j.factorial : ℝ) * y ^ j := by
    refine Finset.sum_congr rfl (fun k _ => ?_)
    rw [iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_Icc (by norm_num))
      (contDiffAt_le_top (Wy_contDiffAt hh1 0) k) (left_mem_Icc.2 (by norm_num)), smul_eq_mul]
    simp only [G]; ring
  rw [e] at T
  calc _ ≤ M * y ^ (K + 1) / (K.factorial : ℝ) := T
    _ = _ := by ring

lemma aj_contDiffOn (j : ℕ) : ContDiffOn ℝ ∞ (fun x => G j (x, 0)) (Iio 1) := by
  intro x hx
  have h1 : ContDiffAt ℝ ∞ (G j) (x, 0) :=
    (G_contDiffOn j).contDiffAt (isOpen_Uw.mem_nhds (by simp [Uw]; exact hx))
  exact (h1.comp x (contDiffAt_id.prodMk contDiffAt_const)).contDiffWithinAt

lemma reindex (K : ℕ) (h s : ℝ) :
    ∑ m ∈ Finset.range (K + 1), h ^ m * Wpoly m s =
      ∑ j ∈ Finset.range (K + 1), (s * h) ^ j / (j.factorial : ℝ) *
        ∑ i ∈ Finset.range (K - j + 1),
          iteratedDeriv i (fun x => G j (x, 0)) 0 / (i.factorial : ℝ) * h ^ i := by
  have step1 : ∑ m ∈ Finset.range (K + 1), h ^ m * Wpoly m s =
      ∑ m ∈ Finset.range (K + 1), ∑ j ∈ Finset.range (m + 1),
        wcoef (m - j) j * h ^ (m - j) * (s * h) ^ j := by
    refine Finset.sum_congr rfl (fun m _ => ?_)
    rw [Wpoly, Finset.mul_sum]
    refine Finset.sum_congr rfl (fun j hj => ?_)
    have hjm : j ≤ m := by simp at hj; omega
    have : h ^ m = h ^ (m - j) * h ^ j := by rw [← pow_add, Nat.sub_add_cancel hjm]
    rw [this, mul_pow]; ring
  rw [step1, Finset.sum_comm' (t' := Finset.range (K + 1)) (s' := fun j => Finset.Ico j (K + 1))
    (by intro m j; simp only [Finset.mem_range, Finset.mem_Ico]; omega)]
  refine Finset.sum_congr rfl (fun j hj => ?_)
  have hjK : j ≤ K := by simp at hj; omega
  rw [Finset.sum_Ico_eq_sum_range, Finset.mul_sum,
    show K + 1 - j = K - j + 1 by omega]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [show j + i - j = i by omega, wcoef]
  simp only [G]
  field_simp

lemma mono_bd {h y : ℝ} (hh : 0 ≤ h) (hy : 0 ≤ y) {j K : ℕ} (hj : j ≤ K) :
    y ^ j * h ^ (K - j + 1) ≤ (h + y) ^ (K + 1) := by
  have e : (h + y) ^ (K + 1) = (h + y) ^ j * (h + y) ^ (K - j + 1) := by
    rw [← pow_add]; congr 1; omega
  rw [e]
  exact mul_le_mul (pow_le_pow_left₀ hy (by linarith) j) (pow_le_pow_left₀ hh (by linarith) _)
    (by positivity) (by positivity)

theorem W_taylor (K : ℕ) : ∃ C h₀ : ℝ, 0 < h₀ ∧ ∀ h s : ℝ, 0 < h → h ≤ h₀ → 0 ≤ s →
    s ≤ 1 / Real.sqrt h →
    |X02.W h (s * h) - ∑ m ∈ Finset.range (K + 1), h ^ m * Wpoly m s|
      ≤ C * (h * (1 + s)) ^ (K + 1) := by
  obtain ⟨C1, hC1⟩ := taylor_y K
  have hx := fun j : ℕ => taylor_x (aj_contDiffOn j) (K - j)
  choose Cx hCx using hx
  refine ⟨|C1| + ∑ j ∈ Finset.range (K + 1), |Cx j|, 1 / 4, by norm_num,
    fun h s hh hh4 hs hs1 => ?_⟩
  set y := s * h with hy
  have hsq : 0 < Real.sqrt h := Real.sqrt_pos.2 hh
  have hsq1 : Real.sqrt h ≤ 1 / 2 := by
    have : Real.sqrt (1 / 4) = 1 / 2 := by
      rw [show (1 / 4 : ℝ) = (1 / 2) ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]
    rw [← this]; exact Real.sqrt_le_sqrt hh4
  have hy0 : 0 ≤ y := mul_nonneg hs hh.le
  have hy1 : y ≤ 1 := by
    have : s * h ≤ 1 / Real.sqrt h * h := mul_le_mul_of_nonneg_right hs1 hh.le
    have e : 1 / Real.sqrt h * h = Real.sqrt h := by
      have := Real.sq_sqrt hh.le
      field_simp; linarith
    linarith
  have hhI : h ∈ Icc (0 : ℝ) (1 / 2) := ⟨hh.le, by linarith⟩
  have hyI : y ∈ Icc (0 : ℝ) 1 := ⟨hy0, hy1⟩
  have hhy : h * (1 + s) = h + y := by rw [hy]; ring
  rw [hhy, reindex]
  set A := ∑ j ∈ Finset.range (K + 1), G j (h, 0) / (j.factorial : ℝ) * y ^ j
  have key : X02.W h y - ∑ j ∈ Finset.range (K + 1), y ^ j / (j.factorial : ℝ) *
        ∑ i ∈ Finset.range (K - j + 1),
          iteratedDeriv i (fun x => G j (x, 0)) 0 / (i.factorial : ℝ) * h ^ i =
      (X02.W h y - A) + ∑ j ∈ Finset.range (K + 1), y ^ j / (j.factorial : ℝ) *
        (G j (h, 0) - ∑ i ∈ Finset.range (K - j + 1),
          iteratedDeriv i (fun x => G j (x, 0)) 0 / (i.factorial : ℝ) * h ^ i) := by
    simp only [A, mul_sub, Finset.sum_sub_distrib]
    have : ∑ j ∈ Finset.range (K + 1), G j (h, 0) / (j.factorial : ℝ) * y ^ j =
        ∑ j ∈ Finset.range (K + 1), y ^ j / (j.factorial : ℝ) * G j (h, 0) := by
      refine Finset.sum_congr rfl (fun j _ => ?_); ring
    rw [this]; ring
  rw [key]
  have b1 : |X02.W h y - A| ≤ |C1| * (h + y) ^ (K + 1) := by
    have := hC1 h hhI y hyI
    calc _ ≤ C1 * y ^ (K + 1) := this
      _ ≤ |C1| * y ^ (K + 1) := mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity)
      _ ≤ |C1| * (h + y) ^ (K + 1) :=
          mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hy0 (by linarith) _) (abs_nonneg _)
  have b2 : ∀ j ∈ Finset.range (K + 1), |y ^ j / (j.factorial : ℝ) *
        (G j (h, 0) - ∑ i ∈ Finset.range (K - j + 1),
          iteratedDeriv i (fun x => G j (x, 0)) 0 / (i.factorial : ℝ) * h ^ i)|
        ≤ |Cx j| * (h + y) ^ (K + 1) := by
    intro j hj
    have hjK : j ≤ K := by simp at hj; omega
    have hf : (1 : ℝ) ≤ (j.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos j
    have hyj : y ^ j / (j.factorial : ℝ) ≤ y ^ j := div_le_self (by positivity) hf
    have hT := hCx j h hhI
    rw [abs_mul, abs_of_nonneg (by positivity : (0 : ℝ) ≤ y ^ j / (j.factorial : ℝ))]
    calc y ^ j / (j.factorial : ℝ) * |_| ≤ y ^ j * (|Cx j| * h ^ (K - j + 1)) := by
          refine mul_le_mul hyj (hT.trans ?_) (abs_nonneg _) (by positivity)
          exact mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity)
      _ = |Cx j| * (y ^ j * h ^ (K - j + 1)) := by ring
      _ ≤ |Cx j| * (h + y) ^ (K + 1) :=
          mul_le_mul_of_nonneg_left (mono_bd hh.le hy0 hjK) (abs_nonneg _)
  calc _ ≤ |X02.W h y - A| + |∑ j ∈ Finset.range (K + 1), y ^ j / (j.factorial : ℝ) *
        (G j (h, 0) - ∑ i ∈ Finset.range (K - j + 1),
          iteratedDeriv i (fun x => G j (x, 0)) 0 / (i.factorial : ℝ) * h ^ i)| := abs_add_le _ _
    _ ≤ |C1| * (h + y) ^ (K + 1) + ∑ j ∈ Finset.range (K + 1), |Cx j| * (h + y) ^ (K + 1) :=
        add_le_add b1 ((Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum b2))
    _ = _ := by rw [← Finset.sum_mul]; ring

end Q02

/-- Lemma `lem:pip:W`: the Taylor bound of every order (Taylor's theorem in `y`, then in `x`,
for the smooth function `W` on `x < 1`), and `0 < W ≤ 1` (`X02.W_pos_le`). -/
theorem Q02.proof_4 : Q02.stmt_4 :=
  ⟨Q02.W_taylor, fun _ _ hh h1 hs => X02.W_pos_le hh.le h1 hs⟩


namespace Q02

lemma choose_eq_binom (a : ℝ) (n : ℕ) : Ring.choose a n = binom a n := by
  rw [Ring.choose_eq_smul, binom, smul_eq_mul, ← descPochhammer_map (algebraMap ℤ ℝ),
    Polynomial.eval_map, ← Polynomial.aeval_def, Polynomial.aeval_eq_smeval]
  ring

lemma binom_hasSum (a x : ℝ) (hx : |x| < 1) :
    HasSum (fun n => binom a n * x ^ n) ((1 + x) ^ a) := by
  have hb : x ∈ Metric.eball (0 : ℝ) 1 := by
    rw [Metric.mem_eball, edist_dist, ENNReal.ofReal_lt_one, Real.dist_eq, sub_zero]; exact hx
  have := (Real.one_add_rpow_hasFPowerSeriesOnBall_zero (a := a)).hasSum hb
  simp only [binomialSeries_apply, List.ofFn_const, List.prod_replicate, smul_eq_mul,
    zero_add] at this
  convert this using 2 with n
  rw [choose_eq_binom]

lemma desc_abs {h : ℝ} (hh : 0 < h) (h1 : h < 1) (n : ℕ) :
    |(descPochhammer ℝ (n + 1)).eval h| ≤ h * (n.factorial : ℝ) := by
  induction n with
  | zero => simp [abs_of_pos hh]
  | succ n ih =>
    rw [descPochhammer_succ_eval, abs_mul, Nat.factorial_succ]
    have : |h - ((n + 1 : ℕ) : ℝ)| ≤ (n : ℝ) + 1 := by
      rw [abs_le]; push_cast; constructor <;> linarith
    calc _ ≤ h * (n.factorial : ℝ) * ((n : ℝ) + 1) :=
          mul_le_mul ih this (abs_nonneg _) (by positivity)
      _ = _ := by push_cast; ring

lemma binom_abs {h : ℝ} (hh : 0 < h) (h1 : h < 1) (n : ℕ) :
    |binom h (n + 1)| ≤ h / ((n : ℝ) + 1) := by
  rw [binom, abs_div, abs_of_pos (by positivity : (0 : ℝ) < ((n + 1).factorial : ℝ)),
    Nat.factorial_succ, div_le_div_iff₀ (by positivity) (by positivity)]
  have := desc_abs hh h1 n
  push_cast
  nlinarith [show (0 : ℝ) ≤ (n.factorial : ℝ) by positivity]

lemma continuous_Wsh (h : ℝ) : Continuous (fun s => X02.W h (s * h)) := by
  unfold X02.W
  have hb : ∀ s : ℝ, 1 + Real.exp (-(s * h / (1 - h))) ≠ 0 := fun s => by positivity
  exact (by fun_prop : Continuous fun s : ℝ => Real.exp (-(h * (s * h) / (1 - h)))).mul
    ((by fun_prop : Continuous fun s : ℝ => 1 + Real.exp (-(s * h / (1 - h)))).rpow_const
      (fun s => Or.inl (hb s)))

lemma intOn_exp_mul {c : ℝ} (hc : 0 < c) : IntegrableOn (fun s => Real.exp (-(c * s))) (Ioi 0) := by
  have := X02.intOn_pow_exp 0 hc; simpa using this

lemma integral_exp_mul {c : ℝ} (hc : 0 < c) : ∫ s in Ioi (0 : ℝ), Real.exp (-(c * s)) = 1 / c := by
  have := X02.integral_pow_exp 0 hc; simpa using this

/-- The substitution `y = κ s`, `κ = h/(1-h)`. -/
lemma subst_point {h : ℝ} (hh : 0 < h) (h1 : h < 1) (m : ℕ) (s : ℝ) :
    X02.W h (s * h) * Real.exp (-((m : ℝ) * s)) =
      Real.exp (-((m : ℝ) * (h / (1 - h) * s) / h)) * Hn m (h / (1 - h) * s) h := by
  unfold X02.W Hn
  have h1h : 1 - h ≠ 0 := by linarith
  set y := h / (1 - h) * s with hy
  have e1 : s * h / (1 - h) = y := by rw [hy]; field_simp
  have e2 : h * (s * h) / (1 - h) = h * y := by rw [hy]; field_simp
  rw [e1, e2]
  have e3 : 1 + Real.exp (-y) = Real.exp (-y) * (1 + Real.exp y) := by
    rw [mul_add, ← Real.exp_add]; simp [add_comm]
  rw [e3, Real.mul_rpow (Real.exp_pos _).le (by positivity), ← Real.exp_mul]
  have e4 : s = y * (1 - h) / h := by rw [hy]; field_simp
  have key : Real.exp (-(h * y)) * Real.exp (-y * (-1 - h)) * Real.exp (-((m : ℝ) * s)) =
      Real.exp (-((m : ℝ) * y / h)) * Real.exp (((m : ℝ) + 1) * y) := by
    rw [← Real.exp_add, ← Real.exp_add, ← Real.exp_add]; congr 1; rw [e4]; field_simp; ring
  linear_combination ((1 + Real.exp y) ^ (-1 - h)) * key

lemma intOn_subst {h : ℝ} (hh : 0 < h) (h1 : h < 1) (m : ℕ) (hm : 1 ≤ m) :
    IntegrableOn (fun y => Real.exp (-((m : ℝ) * y / h)) * Hn m y h) (Ioi 0) := by
  have hc : 0 < (m : ℝ) * (1 / h - 1) := by
    have : (1 : ℝ) ≤ m := by exact_mod_cast hm
    have : 1 < 1 / h := by rw [lt_one_div one_pos hh]; simpa using h1
    nlinarith
  refine X02.intOn_of_bound ?_ (intOn_exp_mul hc) (fun y hy => ?_)
  · unfold Hn
    have hb : ∀ y : ℝ, 1 + Real.exp y ≠ 0 := fun y => by positivity
    exact (by fun_prop : Continuous fun y : ℝ => Real.exp (-((m : ℝ) * y / h))).mul
      ((by fun_prop : Continuous fun y : ℝ => Real.exp (((m : ℝ) + 1) * y)).mul
      ((by fun_prop : Continuous fun y : ℝ => 1 + Real.exp y).rpow_const (fun y => Or.inl (hb y))))
  · unfold Hn
    rw [abs_of_pos (by positivity)]
    have hle : (1 + Real.exp y) ^ (-1 - h) ≤ (Real.exp y) ^ (-1 - h) :=
      Real.rpow_le_rpow_of_nonpos (Real.exp_pos y) (by linarith [Real.exp_pos y]) (by linarith)
    rw [← Real.exp_mul] at hle
    calc _ ≤ Real.exp (-((m : ℝ) * y / h)) * (Real.exp (((m : ℝ) + 1) * y) *
          Real.exp (y * (-1 - h))) := by gcongr
      _ = Real.exp (-((m : ℝ) * (1 / h - 1) * y) - h * y) := by
          rw [← Real.exp_add, ← Real.exp_add]; congr 1; field_simp; ring
      _ ≤ _ := by
          apply Real.exp_le_exp.mpr
          nlinarith

lemma integral_subst {h : ℝ} (hh : 0 < h) (h1 : h < 1) (m : ℕ) :
    ∫ s in Ioi (0 : ℝ), X02.W h (s * h) * Real.exp (-((m : ℝ) * s)) =
      (1 - h) / h * ∫ y in Ioi (0 : ℝ), Real.exp (-((m : ℝ) * y / h)) * Hn m y h := by
  have hk : 0 < h / (1 - h) := div_pos hh (by linarith)
  rw [setIntegral_congr_fun measurableSet_Ioi (fun s _ => subst_point hh h1 m s),
    integral_comp_mul_left_Ioi (fun y => Real.exp (-((m : ℝ) * y / h)) * Hn m y h) 0 hk,
    mul_zero, smul_eq_mul, inv_div]

end Q02

theorem Q02.proof_12 : Q02.stmt_12 := by
  intro h hh h1
  refine ⟨fun n hn => Q02.intOn_subst hh h1 n hn, ?_⟩
  set F : ℕ → ℝ → ℝ := fun n s =>
    Q02.binom h (n + 1) * (X02.W h (s * h) * Real.exp (-(((n + 1 : ℕ) : ℝ) * s))) with hF
  have hpos : ∀ n : ℕ, (0 : ℝ) < ((n + 1 : ℕ) : ℝ) := fun n => by positivity
  have hWb : ∀ s : ℝ, 0 < s → |X02.W h (s * h)| ≤ 1 := fun s hs => by
    have := X02.W_pos_le hh.le h1 hs.le
    rw [abs_of_pos this.1]; exact this.2
  have hint : ∀ n, Integrable (F n) (volume.restrict (Ioi 0)) := fun n => by
    refine X02.intOn_of_bound (continuous_const.mul ((Q02.continuous_Wsh h).mul (by fun_prop)))
      ((Q02.intOn_exp_mul (hpos n)).const_mul |Q02.binom h (n + 1)|) (fun s hs => ?_)
    simp only [hF, abs_mul, abs_of_pos (Real.exp_pos _)]
    have := hWb s hs
    gcongr
    calc _ ≤ 1 * Real.exp (-(((n + 1 : ℕ) : ℝ) * s)) :=
          mul_le_mul_of_nonneg_right this (Real.exp_pos _).le
      _ = _ := one_mul _
  have hnorm : ∀ n, ∫ s in Ioi (0 : ℝ), ‖F n s‖ ≤ h / ((n : ℝ) + 1) ^ 2 := fun n => by
    calc ∫ s in Ioi (0 : ℝ), ‖F n s‖
        ≤ ∫ s in Ioi (0 : ℝ), |Q02.binom h (n + 1)| * Real.exp (-(((n + 1 : ℕ) : ℝ) * s)) := by
          refine setIntegral_mono_on (hint n).norm
            ((Q02.intOn_exp_mul (hpos n)).const_mul _) measurableSet_Ioi (fun s hs => ?_)
          simp only [hF, Real.norm_eq_abs, abs_mul, abs_of_pos (Real.exp_pos _)]
          have := hWb s hs
          gcongr
          calc _ ≤ 1 * Real.exp (-(((n + 1 : ℕ) : ℝ) * s)) :=
                mul_le_mul_of_nonneg_right this (Real.exp_pos _).le
            _ = _ := one_mul _
      _ = |Q02.binom h (n + 1)| * (1 / ((n : ℝ) + 1)) := by
          rw [integral_const_mul, Q02.integral_exp_mul (hpos n)]; push_cast; rfl
      _ ≤ h / ((n : ℝ) + 1) * (1 / ((n : ℝ) + 1)) :=
          mul_le_mul_of_nonneg_right (Q02.binom_abs hh h1 n) (by positivity)
      _ = _ := by rw [div_mul_div_comm, mul_one, sq]
  have hsumN : Summable fun n => ∫ s in Ioi (0 : ℝ), ‖F n s‖ := by
    refine Summable.of_nonneg_of_le (fun n => integral_nonneg (fun s => norm_nonneg _)) hnorm ?_
    exact ((X02.summable_inv_pow 2 le_rfl).mul_left h).congr (fun n => by rw [mul_one_div])
  have key := integral_tsum_of_summable_integral_norm hint hsumN
  have hsumI : Summable fun n => ∫ s in Ioi (0 : ℝ), F n s :=
    Summable.of_norm_bounded hsumN (fun n => norm_integral_le_integral_norm _)
  have hpt : ∀ s ∈ Ioi (0 : ℝ), ∑' n, F n s =
      X02.W h (s * h) * ((1 + Real.exp (-s)) ^ h - 1) := fun s hs => by
    have hx : |Real.exp (-s)| < 1 := by
      rw [abs_of_pos (Real.exp_pos _)]; exact Real.exp_lt_one_iff.mpr (by linarith [mem_Ioi.mp hs])
    have hb := Q02.binom_hasSum h _ hx
    have hb1 := (hasSum_nat_add_iff' 1).mpr hb
    simp only [Finset.range_one, Finset.sum_singleton, pow_zero, mul_one] at hb1
    have hb0 : Q02.binom h 0 = 1 := by simp [Q02.binom]
    rw [hb0] at hb1
    have := (hb1.mul_left (X02.W h (s * h))).tsum_eq
    rw [← this]
    congr 1; funext n; simp only [hF]
    rw [← Real.exp_nat_mul]; push_cast; ring_nf
  have hQ : X02.Q h = ∑' n, ∫ s in Ioi (0 : ℝ), F n s := by
    rw [key, X02.Q]; exact (setIntegral_congr_fun measurableSet_Ioi hpt).symm
  rw [(Q02.proof_6 h hh h1).2, hQ, ← (hsumI.hasSum.mul_left (h ^ 2 / (1 - h))).tsum_eq]
  convert (hsumI.hasSum.mul_left (h ^ 2 / (1 - h))) using 1
  · funext n
    simp only [hF]
    rw [integral_const_mul, Q02.integral_subst hh h1]
    have h1h : 1 - h ≠ 0 := by linarith
    field_simp
  · rw [(hsumI.hasSum.mul_left (h ^ 2 / (1 - h))).tsum_eq]

/-! ### Lemma `lem:pip:C45`: the Beta function split at `1` (proofs/sec_pip_C45.tex) -/

namespace Q02

noncomputable def phiI (a b s : ℝ) : ℝ := Real.exp (-(a * s)) * (Real.exp (-(b * X02.calL s)) - 1)
noncomputable def Phi (a b : ℝ) : ℝ := ∫ s in Ioi (0 : ℝ), phiI a b s
lemma continuous_phiI (a b : ℝ) : Continuous (phiI a b) := by
  unfold phiI
  exact (by fun_prop : Continuous fun s : ℝ => Real.exp (-(a * s))).mul
    ((Real.continuous_exp.comp ((continuous_const.mul X02.continuous_calL).neg)).sub
      continuous_const)
lemma intOn_phiI {a b : ℝ} (ha : 0 < a) (hb : 0 ≤ b) : IntegrableOn (phiI a b) (Ioi 0) := by
  refine X02.intOn_of_bound (continuous_phiI a b) (intOn_exp_mul ha) (fun s _ => ?_)
  unfold phiI
  have h0 := X02.calL_nonneg s
  have h1 : Real.exp (-(b * X02.calL s)) ≤ 1 := Real.exp_le_one_iff.mpr (by nlinarith)
  have h2 : 0 < Real.exp (-(b * X02.calL s)) := Real.exp_pos _
  rw [abs_mul, abs_of_pos (Real.exp_pos _), abs_of_nonpos (by linarith)]
  nlinarith [Real.exp_pos (-(a * s))]

/-! ### The `J` side: `Φ` to order five -/

/-- the index set `r ≥ 1`, `r + j ≤ 4` -/
def idx : Finset (ℕ × ℕ) := {(1, 0), (1, 1), (1, 2), (1, 3), (2, 0), (2, 1), (2, 2), (3, 0), (3, 1), (4, 0)}

lemma idx_r {p : ℕ × ℕ} (hp : p ∈ idx) : 1 ≤ p.1 := by
  simp only [idx, Finset.mem_insert, Finset.mem_singleton] at hp
  rcases hp with h | h | h | h | h | h | h | h | h | h <;> simp [h]

/-- the coefficient `(-a)^j (-b)^r / (j! r!)` -/
noncomputable def cc (a b : ℝ) (p : ℕ × ℕ) : ℝ :=
  (-a) ^ p.2 * (-b) ^ p.1 / ((p.2.factorial : ℝ) * (p.1.factorial : ℝ))

noncomputable def Ppt (a b s : ℝ) : ℝ := ∑ p ∈ idx, cc a b p * (s ^ p.2 * X02.calL s ^ p.1)

noncomputable def Pint (a b : ℝ) : ℝ := ∑ p ∈ idx, cc a b p * J p.1 p.2

lemma intOn_Ppt (a b : ℝ) : IntegrableOn (Ppt a b) (Ioi 0) := by
  unfold Ppt
  exact integrable_finsetSum _ (fun p hp => (intOn_pow_calL p.1 p.2 (idx_r hp)).const_mul _)

lemma integral_Ppt (a b : ℝ) : ∫ s in Ioi (0 : ℝ), Ppt a b s = Pint a b := by
  unfold Ppt Pint
  rw [integral_finsetSum _ (fun p hp => (intOn_pow_calL p.1 p.2 (idx_r hp)).const_mul _)]
  congr 1; funext p; rw [integral_const_mul]; rfl

lemma Ppt_eq (a b s : ℝ) : Ppt a b s =
    let u := a * s
    let w := b * X02.calL s
    (-w) * (1 - u + u ^ 2 / 2 - u ^ 3 / 6) + w ^ 2 / 2 * (1 - u + u ^ 2 / 2)
      + (-(w ^ 3 / 6)) * (1 - u) + w ^ 4 / 24 * 1 := by
  simp only [Ppt, idx, cc]
  simp [Finset.sum_insert, Nat.factorial]
  ring

lemma exp_neg_rem {u : ℝ} (hu : 0 ≤ u) (n : ℕ) :
    |Real.exp (-u) - ∑ m ∈ Finset.range n, (-u) ^ m / (m.factorial : ℝ)| ≤ u ^ n * Real.exp u := by
  have h := Complex.norm_exp_sub_sum_le_norm_mul_exp ((-u : ℝ) : ℂ) n
  have e : Complex.exp ((-u : ℝ) : ℂ) - ∑ m ∈ Finset.range n, ((-u : ℝ) : ℂ) ^ m / (m.factorial : ℂ) =
      ((Real.exp (-u) - ∑ m ∈ Finset.range n, (-u) ^ m / (m.factorial : ℝ) : ℝ) : ℂ) := by
    push_cast; rfl
  rw [e, Complex.norm_real, Real.norm_eq_abs, Complex.norm_real, Real.norm_eq_abs,
    abs_neg, abs_of_nonneg hu] at h
  exact h

/-- the majorant `64 e^{-3s/4} (1 + s + s² + s³ + s⁴)` -/
noncomputable def gmaj (s : ℝ) : ℝ := 64 * Real.exp (-(3 / 4 * s)) * (1 + s + s ^ 2 + s ^ 3 + s ^ 4)

lemma intOn_gmaj : IntegrableOn gmaj (Ioi 0) := by
  have h : ∀ j : ℕ, IntegrableOn (fun s : ℝ => s ^ j * Real.exp (-(3 / 4 * s))) (Ioi 0) :=
    fun j => X02.intOn_pow_exp j (by norm_num)
  have h5 := ((((h 0).add (h 1)).add (h 2)).add (h 3)).add (h 4)
  have h64 : IntegrableOn (fun s : ℝ => 64 * (s ^ 0 * Real.exp (-(3 / 4 * s)) +
      s ^ 1 * Real.exp (-(3 / 4 * s)) + s ^ 2 * Real.exp (-(3 / 4 * s)) +
      s ^ 3 * Real.exp (-(3 / 4 * s)) + s ^ 4 * Real.exp (-(3 / 4 * s)))) (Ioi 0) :=
    h5.const_mul 64
  refine h64.congr_fun (fun s _ => ?_) measurableSet_Ioi
  simp only [gmaj, pow_zero, pow_one, one_mul]; ring

lemma pt_bound {x a b s : ℝ} (hx : x ≤ 1 / 4) (ha0 : 0 ≤ a) (hax : a ≤ x)
    (hb0 : 0 ≤ b) (hbx : b ≤ 2 * x) (hs : 0 ≤ s) :
    |phiI a b s - Ppt a b s| ≤ x ^ 5 * gmaj s := by
  set L := X02.calL s with hL
  have hL0 : 0 ≤ L := X02.calL_nonneg s
  have hLe : L ≤ Real.exp (-s) := X02.calL_le s
  have he1 : Real.exp (-s) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
  have hL1 : L ≤ 1 := hLe.trans he1
  have hx0 : 0 ≤ x := ha0.trans hax
  set u := a * s with hu
  set w := b * L with hw
  have hu0 : 0 ≤ u := mul_nonneg ha0 hs
  have hw0 : 0 ≤ w := mul_nonneg hb0 hL0
  have hux : u ≤ x * s := mul_le_mul_of_nonneg_right hax hs
  have hwx : w ≤ 2 * x * L := mul_le_mul_of_nonneg_right hbx hL0
  have hw1 : w ≤ 1 / 2 := by nlinarith
  have hew : Real.exp w ≤ 2 := by
    have := Real.exp_le_exp.mpr hw1
    have h2 : Real.exp (1 / 2 : ℝ) ≤ 2 := by
      have := Real.exp_one_lt_d9
      have h3 : Real.exp (1 / 2 : ℝ) ^ 2 = Real.exp 1 := by
        rw [← Real.exp_nat_mul]; norm_num
      nlinarith [Real.exp_pos (1 / 2 : ℝ)]
    linarith
  have heu : Real.exp u ≤ Real.exp (s / 4) := Real.exp_le_exp.mpr (by nlinarith)
  have eid : phiI a b s - Ppt a b s =
      Real.exp (-u) * (Real.exp (-w) - ∑ m ∈ Finset.range 5, (-w) ^ m / (m.factorial : ℝ))
      + (-w) * (Real.exp (-u) - ∑ m ∈ Finset.range 4, (-u) ^ m / (m.factorial : ℝ))
      + w ^ 2 / 2 * (Real.exp (-u) - ∑ m ∈ Finset.range 3, (-u) ^ m / (m.factorial : ℝ))
      + (-(w ^ 3 / 6)) * (Real.exp (-u) - ∑ m ∈ Finset.range 2, (-u) ^ m / (m.factorial : ℝ))
      + w ^ 4 / 24 * (Real.exp (-u) - ∑ m ∈ Finset.range 1, (-u) ^ m / (m.factorial : ℝ)) := by
    rw [Ppt_eq]
    simp only [phiI, Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial]
    rw [← hL, ← hu, ← hw]
    push_cast
    ring
  have r5 := exp_neg_rem hw0 5
  have r4 := exp_neg_rem hu0 4
  have r3 := exp_neg_rem hu0 3
  have r2 := exp_neg_rem hu0 2
  have r1 := exp_neg_rem hu0 1
  have heu1 : Real.exp (-u) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
  -- crude bounds of each piece
  have hLr : ∀ r : ℕ, 1 ≤ r → L ^ r ≤ Real.exp (-s) := fun r hr =>
    (pow_le_of_le_one hL0 hL1 (by omega)).trans hLe
  have hes : Real.exp (-s) * Real.exp (s / 4) = Real.exp (-(3 / 4 * s)) := by
    rw [← Real.exp_add]; congr 1; ring
  have hwr : ∀ r : ℕ, 1 ≤ r → w ^ r ≤ 2 ^ r * x ^ r * Real.exp (-s) := by
    intro r hr
    calc w ^ r ≤ (2 * x * L) ^ r := pow_le_pow_left₀ hw0 hwx r
      _ = 2 ^ r * x ^ r * L ^ r := by ring
      _ ≤ 2 ^ r * x ^ r * Real.exp (-s) := by gcongr; exact hLr r hr
  have hur : ∀ n : ℕ, u ^ n * Real.exp u ≤ x ^ n * s ^ n * Real.exp (s / 4) := by
    intro n
    calc u ^ n * Real.exp u ≤ (x * s) ^ n * Real.exp (s / 4) := by
          gcongr
      _ = _ := by ring
  have E : Real.exp (-(3 / 4 * s)) ≥ 0 := (Real.exp_pos _).le
  have T0 : |Real.exp (-u) * (Real.exp (-w) - ∑ m ∈ Finset.range 5, (-w) ^ m / (m.factorial : ℝ))|
      ≤ 64 * x ^ 5 * Real.exp (-(3 / 4 * s)) := by
    rw [abs_mul, abs_of_pos (Real.exp_pos _)]
    have h1 : w ^ 5 * Real.exp w ≤ 2 ^ 5 * x ^ 5 * Real.exp (-s) * 2 :=
      mul_le_mul (hwr 5 (by norm_num)) hew (Real.exp_pos _).le (by positivity)
    have h2 : Real.exp (-s) ≤ Real.exp (-(3 / 4 * s)) := Real.exp_le_exp.mpr (by linarith)
    have h3 : Real.exp (-u) * |Real.exp (-w) - ∑ m ∈ Finset.range 5, (-w) ^ m / (m.factorial : ℝ)|
        ≤ 1 * (w ^ 5 * Real.exp w) := mul_le_mul heu1 r5 (abs_nonneg _) zero_le_one
    have : (0 : ℝ) ≤ x ^ 5 := pow_nonneg hx0 5
    nlinarith
  have Tr : ∀ (r n : ℕ) (c : ℝ), 1 ≤ r → r ≤ 4 → r + n = 5 → 0 < c → c ≤ 1 →
      |c * w ^ r| * |Real.exp (-u) - ∑ m ∈ Finset.range n, (-u) ^ m / (m.factorial : ℝ)|
        ≤ 16 * x ^ 5 * s ^ n * Real.exp (-(3 / 4 * s)) := by
    intro r n c hr hr4 hrn hc hc1
    rw [abs_mul, abs_of_pos hc, abs_of_nonneg (pow_nonneg hw0 r)]
    have hrem := exp_neg_rem hu0 n
    have hcw : c * w ^ r ≤ 2 ^ r * x ^ r * Real.exp (-s) :=
      (mul_le_of_le_one_left (pow_nonneg hw0 r) hc1).trans (hwr r hr)
    have h2r : (2 : ℝ) ^ r ≤ 16 := by
      calc (2 : ℝ) ^ r ≤ 2 ^ 4 := pow_le_pow_right₀ (by norm_num) hr4
        _ = 16 := by norm_num
    calc c * w ^ r * |Real.exp (-u) - ∑ m ∈ Finset.range n, (-u) ^ m / (m.factorial : ℝ)|
        ≤ (2 ^ r * x ^ r * Real.exp (-s)) * (x ^ n * s ^ n * Real.exp (s / 4)) :=
          mul_le_mul hcw (hrem.trans (hur n)) (abs_nonneg _) (by positivity)
      _ = 2 ^ r * (x ^ (r + n)) * s ^ n * (Real.exp (-s) * Real.exp (s / 4)) := by ring
      _ ≤ 16 * x ^ 5 * s ^ n * Real.exp (-(3 / 4 * s)) := by
          rw [hrn, hes]; gcongr
  have T1 := Tr 1 4 1 le_rfl (by norm_num) rfl one_pos le_rfl
  have T2 := Tr 2 3 (1 / 2) (by norm_num) (by norm_num) rfl (by norm_num) (by norm_num)
  have T3 := Tr 3 2 (1 / 6) (by norm_num) (by norm_num) rfl (by norm_num) (by norm_num)
  have T4 := Tr 4 1 (1 / 24) (by norm_num) (by norm_num) rfl (by norm_num) (by norm_num)
  rw [eid]
  have a1 : |(-w) * (Real.exp (-u) - ∑ m ∈ Finset.range 4, (-u) ^ m / (m.factorial : ℝ))| =
      |1 * w ^ 1| * |Real.exp (-u) - ∑ m ∈ Finset.range 4, (-u) ^ m / (m.factorial : ℝ)| := by
    rw [abs_mul, abs_neg, one_mul, pow_one]
  have a2 : |w ^ 2 / 2 * (Real.exp (-u) - ∑ m ∈ Finset.range 3, (-u) ^ m / (m.factorial : ℝ))| =
      |1 / 2 * w ^ 2| * |Real.exp (-u) - ∑ m ∈ Finset.range 3, (-u) ^ m / (m.factorial : ℝ)| := by
    rw [abs_mul]; ring_nf
  have a3 : |(-(w ^ 3 / 6)) * (Real.exp (-u) - ∑ m ∈ Finset.range 2, (-u) ^ m / (m.factorial : ℝ))| =
      |1 / 6 * w ^ 3| * |Real.exp (-u) - ∑ m ∈ Finset.range 2, (-u) ^ m / (m.factorial : ℝ)| := by
    rw [abs_mul, abs_neg]; ring_nf
  have a4 : |w ^ 4 / 24 * (Real.exp (-u) - ∑ m ∈ Finset.range 1, (-u) ^ m / (m.factorial : ℝ))| =
      |1 / 24 * w ^ 4| * |Real.exp (-u) - ∑ m ∈ Finset.range 1, (-u) ^ m / (m.factorial : ℝ)| := by
    rw [abs_mul]; ring_nf
  set A0 := Real.exp (-u) * (Real.exp (-w) - ∑ m ∈ Finset.range 5, (-w) ^ m / (m.factorial : ℝ))
  set A1 := (-w) * (Real.exp (-u) - ∑ m ∈ Finset.range 4, (-u) ^ m / (m.factorial : ℝ))
  set A2 := w ^ 2 / 2 * (Real.exp (-u) - ∑ m ∈ Finset.range 3, (-u) ^ m / (m.factorial : ℝ))
  set A3 := (-(w ^ 3 / 6)) * (Real.exp (-u) - ∑ m ∈ Finset.range 2, (-u) ^ m / (m.factorial : ℝ))
  set A4 := w ^ 4 / 24 * (Real.exp (-u) - ∑ m ∈ Finset.range 1, (-u) ^ m / (m.factorial : ℝ))
  have tri : |A0 + A1 + A2 + A3 + A4| ≤ |A0| + |A1| + |A2| + |A3| + |A4| := by
    have := abs_add_le (A0 + A1 + A2 + A3) A4
    have := abs_add_le (A0 + A1 + A2) A3
    have := abs_add_le (A0 + A1) A2
    have := abs_add_le A0 A1
    linarith
  have hs1 : s ^ 4 * Real.exp (-(3 / 4 * s)) ≥ 0 := by positivity
  unfold gmaj
  have : (0 : ℝ) ≤ x ^ 5 := pow_nonneg hx0 5
  have key : 64 * x ^ 5 * Real.exp (-(3 / 4 * s)) + 16 * x ^ 5 * s ^ 4 * Real.exp (-(3 / 4 * s))
      + 16 * x ^ 5 * s ^ 3 * Real.exp (-(3 / 4 * s)) + 16 * x ^ 5 * s ^ 2 * Real.exp (-(3 / 4 * s))
      + 16 * x ^ 5 * s ^ 1 * Real.exp (-(3 / 4 * s))
      ≤ x ^ 5 * (64 * Real.exp (-(3 / 4 * s)) * (1 + s + s ^ 2 + s ^ 3 + s ^ 4)) := by
    have e : x ^ 5 * (64 * Real.exp (-(3 / 4 * s)) * (1 + s + s ^ 2 + s ^ 3 + s ^ 4)) -
        (64 * x ^ 5 * Real.exp (-(3 / 4 * s)) + 16 * x ^ 5 * s ^ 4 * Real.exp (-(3 / 4 * s))
        + 16 * x ^ 5 * s ^ 3 * Real.exp (-(3 / 4 * s)) + 16 * x ^ 5 * s ^ 2 * Real.exp (-(3 / 4 * s))
        + 16 * x ^ 5 * s ^ 1 * Real.exp (-(3 / 4 * s))) =
        x ^ 5 * Real.exp (-(3 / 4 * s)) * (48 * s + 48 * s ^ 2 + 48 * s ^ 3 + 48 * s ^ 4) := by ring
    have : 0 ≤ x ^ 5 * Real.exp (-(3 / 4 * s)) * (48 * s + 48 * s ^ 2 + 48 * s ^ 3 + 48 * s ^ 4) := by
      positivity
    linarith
  have T1' : |A1| ≤ 16 * x ^ 5 * s ^ 4 * Real.exp (-(3 / 4 * s)) := by rw [a1]; exact T1
  have T2' : |A2| ≤ 16 * x ^ 5 * s ^ 3 * Real.exp (-(3 / 4 * s)) := by rw [a2]; exact T2
  have T3' : |A3| ≤ 16 * x ^ 5 * s ^ 2 * Real.exp (-(3 / 4 * s)) := by rw [a3]; exact T3
  have T4' : |A4| ≤ 16 * x ^ 5 * s ^ 1 * Real.exp (-(3 / 4 * s)) := by rw [a4]; exact T4
  exact tri.trans ((add_le_add (add_le_add (add_le_add (add_le_add T0 T1') T2') T3') T4').trans key)

lemma phi_approx {x a b : ℝ} (hx : x ≤ 1 / 4) (ha0 : 0 < a) (hax : a ≤ x)
    (hb0 : 0 ≤ b) (hbx : b ≤ 2 * x) :
    |Phi a b - Pint a b| ≤ x ^ 5 * ∫ s in Ioi (0 : ℝ), gmaj s := by
  rw [← integral_Ppt, Phi, ← integral_sub (intOn_phiI ha0 hb0) (intOn_Ppt a b),
    ← integral_const_mul]
  rw [← Real.norm_eq_abs]
  refine norm_integral_le_of_norm_le (intOn_gmaj.const_mul _) ?_
  refine (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall (fun s hs => ?_))
  rw [Real.norm_eq_abs]
  exact pt_bound hx ha0.le hax hb0 hbx (le_of_lt hs)


/-- the Beta integrand `t^{a-1}(1-t)^{b-1}` -/
noncomputable def bI (a b t : ℝ) : ℝ := t ^ (a - 1) * (1 - t) ^ (b - 1)

lemma subst_pt (a b s : ℝ) :
    |-(Real.exp s) / (1 + Real.exp s) ^ 2| * bI a b (1 + Real.exp s)⁻¹ =
      Real.exp (-(a * s)) * Real.exp (-((a + b) * X02.calL s)) := by
  set E := Real.exp s with hE_def
  have hE : 0 < E := Real.exp_pos _
  have h1E : 0 < 1 + E := by linarith
  set ℓ := Real.log (1 + E) with hℓ_def
  have hℓ : Real.exp ℓ = 1 + E := Real.exp_log h1E
  have ht : 0 < (1 + E)⁻¹ := inv_pos.2 h1E
  have hlt : Real.log (1 + E)⁻¹ = -ℓ := Real.log_inv _
  have h1t : 1 - (1 + E)⁻¹ = E * (1 + E)⁻¹ := by field_simp; ring
  have h1t0 : 0 < E * (1 + E)⁻¹ := mul_pos hE ht
  have hl1t : Real.log (E * (1 + E)⁻¹) = s - ℓ := by
    rw [Real.log_mul hE.ne' ht.ne', hlt, hE_def, Real.log_exp]; ring
  have f1 : ((1 + E)⁻¹) ^ (a - 1) = Real.exp (-ℓ * (a - 1)) := by
    rw [Real.rpow_def_of_pos ht, hlt]
  have f2 : (1 - (1 + E)⁻¹) ^ (b - 1) = Real.exp ((s - ℓ) * (b - 1)) := by
    rw [h1t, Real.rpow_def_of_pos h1t0, hl1t]
  have f3 : |-E / (1 + E) ^ 2| = Real.exp (s - 2 * ℓ) := by
    rw [abs_div, abs_neg, abs_of_pos hE, abs_of_pos (by positivity), Real.exp_sub, ← hE_def,
      show 2 * ℓ = ℓ + ℓ by ring, Real.exp_add, hℓ]
    ring
  have hL : X02.calL s = ℓ - s := by rw [hℓ_def, hE_def, X02.phi_eq]; ring
  unfold bI
  rw [f1, f2, f3, hL, ← Real.exp_add, ← Real.exp_add, ← Real.exp_add]
  congr 1; ring

lemma subst_half {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    IntegrableOn (bI a b) (Ioo 0 (1 / 2)) ∧
      ∫ t in Ioo (0 : ℝ) (1 / 2), bI a b t = 1 / a + Phi a (a + b) := by
  set σ : ℝ → ℝ := fun s => (1 + Real.exp s)⁻¹ with hσ
  set σ' : ℝ → ℝ := fun s => -(Real.exp s) / (1 + Real.exp s) ^ 2 with hσ'
  have hder : ∀ s, HasDerivAt σ (σ' s) s := fun s => by
    have h1 := ((Real.hasDerivAt_exp s).const_add 1)
    have := h1.inv (by positivity)
    convert this using 1
  have hinj : InjOn σ (Ioi 0) := by
    intro u _ w _ huw
    simp only [hσ] at huw
    have := inv_injective huw
    have h2 : Real.exp u = Real.exp w := by linarith
    exact Real.exp_injective h2
  have himg : σ '' Ioi 0 = Ioo 0 (1 / 2) := by
    ext t; constructor
    · rintro ⟨s, hs, rfl⟩
      have h1 : 1 < Real.exp s := Real.one_lt_exp_iff.mpr hs
      simp only [hσ]
      refine ⟨by positivity, ?_⟩
      rw [inv_lt_comm₀ (by positivity) (by norm_num)]; linarith
    · rintro ⟨h0, h1⟩
      have hq : 1 < 1 / t - 1 := by
        rw [lt_sub_iff_add_lt, lt_div_iff₀ h0]; linarith
      refine ⟨Real.log (1 / t - 1), Real.log_pos hq, ?_⟩
      simp only [hσ]
      rw [Real.exp_log (by linarith)]
      field_simp; ring
  have hpt : ∀ s, |σ' s| • bI a b (σ s) =
      Real.exp (-(a * s)) * Real.exp (-((a + b) * X02.calL s)) := fun s => by
    simp only [hσ, hσ', smul_eq_mul]; exact subst_pt a b s
  have hsplit : ∀ s, Real.exp (-(a * s)) * Real.exp (-((a + b) * X02.calL s)) =
      Real.exp (-(a * s)) + phiI a (a + b) s := fun s => by unfold phiI; ring
  have hR : IntegrableOn (fun s => Real.exp (-(a * s)) * Real.exp (-((a + b) * X02.calL s)))
      (Ioi 0) := by
    simp_rw [hsplit]; exact (intOn_exp_mul ha).add (intOn_phiI ha (by linarith))
  have hII := (integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ioi
    (fun x _ => (hder x).hasDerivWithinAt) hinj (bI a b))
  have hI := integral_image_eq_integral_abs_deriv_smul measurableSet_Ioi
    (fun x _ => (hder x).hasDerivWithinAt) hinj (bI a b)
  simp_rw [hpt] at hII hI
  rw [himg] at hI hII
  refine ⟨hII.2 hR, ?_⟩
  rw [hI]
  simp_rw [hsplit]
  rw [integral_add (intOn_exp_mul ha) (intOn_phiI ha (by linarith)), integral_exp_mul ha]
  rfl

lemma bI_symm (a b t : ℝ) : bI a b (1 - t) = bI b a t := by
  unfold bI; rw [sub_sub_cancel]; ring

lemma beta_split {x y : ℝ} (hx : 0 < x) (hy : 0 < y) :
    ∫ t in (0 : ℝ)..1, t ^ (x - 1) * (1 - t) ^ (y - 1) =
      1 / x + 1 / y + Phi x (x + y) + Phi y (y + x) := by
  obtain ⟨i1, e1⟩ := subst_half hx hy
  obtain ⟨i2, e2⟩ := subst_half hy hx
  have ii1 : IntervalIntegrable (bI x y) volume 0 (1 / 2) := by
    rw [intervalIntegrable_iff_integrableOn_Ioo_of_le (by norm_num)]; exact i1
  have ii2' : IntervalIntegrable (bI y x) volume 0 (1 / 2) := by
    rw [intervalIntegrable_iff_integrableOn_Ioo_of_le (by norm_num)]; exact i2
  have hsym : ∫ t in (1 / 2 : ℝ)..1, bI x y t = ∫ t in (0 : ℝ)..(1 / 2), bI y x t := by
    have := intervalIntegral.integral_comp_sub_left (a := (0 : ℝ)) (b := 1 / 2) (bI x y) 1
    simp only [bI_symm] at this
    rw [this]; norm_num
  have ii2 : IntervalIntegrable (bI x y) volume (1 / 2) 1 := by
    have := ii2'.comp_sub_left 1
    simp only [bI_symm] at this
    norm_num at this
    exact this.symm
  have hJ : ∫ t in (0 : ℝ)..1, t ^ (x - 1) * (1 - t) ^ (y - 1) = ∫ t in (0 : ℝ)..1, bI x y t := rfl
  rw [hJ, ← intervalIntegral.integral_add_adjacent_intervals ii1 ii2, hsym,
    intervalIntegral.integral_of_le (by norm_num), intervalIntegral.integral_of_le (by norm_num),
    integral_Ioc_eq_integral_Ioo, integral_Ioc_eq_integral_Ioo, e1, e2]
  ring


end Q02

open Filter Topology Finset

namespace Q02

/-- `ζ(n)` as a real number (helper, used for the Gamma side of `lem:pip:C45`). -/
noncomputable def zz (n : ℕ) : ℝ := (riemannZeta n).re

lemma beta_real {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
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


lemma z_hasSum {k : ℕ} (hk : 1 < k) :
    HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ k) (zz k) := by
  have hs : Summable (fun n : ℕ => 1 / (n : ℝ) ^ k) := Real.summable_one_div_nat_pow.mpr hk
  have hz : zz k = ∑' n : ℕ, 1 / (n : ℝ) ^ k := by
    unfold zz
    rw [zeta_nat_eq_tsum_of_gt_one hk]
    rw [show (fun n : ℕ => 1 / (n : ℂ) ^ k) = fun n : ℕ => ((1 / (n : ℝ) ^ k : ℝ) : ℂ) by
      ext n; push_cast; rfl]
    rw [← Complex.ofReal_tsum, Complex.ofReal_re]
  rw [hz]
  have := hs.hasSum
  rw [← hasSum_nat_add_iff' 1] at this
  simpa [zero_pow (show k ≠ 0 by omega)] using this

lemma prod_shift (x : ℝ) (n : ℕ) :
    ∏ j ∈ range (n + 1), (1 + x + (j : ℝ)) =
      ((n + 1).factorial : ℝ) * ∏ k ∈ range (n + 1), (1 + x / ((k : ℝ) + 1)) := by
  rw [show (((n + 1).factorial : ℕ) : ℝ) = ∏ k ∈ range (n + 1), ((k : ℝ) + 1) by
    rw [← Finset.prod_range_add_one_eq_factorial]; push_cast; rfl, ← Finset.prod_mul_distrib]
  refine Finset.prod_congr rfl fun k _ => ?_
  have : (0 : ℝ) < (k : ℝ) + 1 := by positivity
  field_simp
  ring

lemma fac_pos {x : ℝ} (hx : |x| < 1) (k : ℕ) : 0 < 1 + x / ((k : ℝ) + 1) := by
  have h1 : (1 : ℝ) ≤ (k : ℝ) + 1 := by have := k.cast_nonneg (α := ℝ); linarith
  have : |x / ((k : ℝ) + 1)| < 1 := by
    rw [abs_div, abs_of_pos (by linarith : (0:ℝ) < (k : ℝ) + 1), div_lt_one (by linarith)]
    linarith
  linarith [neg_abs_le (x / ((k : ℝ) + 1))]

/-- the partial products `R_m` -/
noncomputable def R (η u v : ℝ) (m : ℕ) : ℝ :=
  ∏ k ∈ range m, ((1 + η / ((k : ℝ) + 1)) / ((1 + u / ((k : ℝ) + 1)) * (1 + v / ((k : ℝ) + 1))))

lemma gammaSeq_ratio {η u v : ℝ} (huv : u + v = η) (hu : |u| < 1) (hv : |v| < 1)
    (hη : |η| < 1) {n : ℕ} (hn : 1 ≤ n) :
    Real.GammaSeq (1 + u) n * Real.GammaSeq (1 + v) n / Real.GammaSeq (1 + η) n =
      (n : ℝ) / ((n : ℝ) + 1) * R η u v (n + 1) := by
  unfold Real.GammaSeq R
  rw [prod_shift, prod_shift, prod_shift, Finset.prod_div_distrib, Finset.prod_mul_distrib]
  have hn0 : (0 : ℝ) < n := by exact_mod_cast hn
  have pu : 0 < ∏ k ∈ range (n + 1), (1 + u / ((k : ℝ) + 1)) :=
    Finset.prod_pos fun k _ => fac_pos hu k
  have pv : 0 < ∏ k ∈ range (n + 1), (1 + v / ((k : ℝ) + 1)) :=
    Finset.prod_pos fun k _ => fac_pos hv k
  have pη : 0 < ∏ k ∈ range (n + 1), (1 + η / ((k : ℝ) + 1)) :=
    Finset.prod_pos fun k _ => fac_pos hη k
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

lemma log_series_k {η u v : ℝ} (huv : u + v = η) (hu : |u| < 1) (hv : |v| < 1)
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

lemma log_R {η u v : ℝ} (hu : |u| < 1) (hv : |v| < 1) (hη : |η| < 1) (m : ℕ) :
    Real.log (R η u v m) = ∑ k ∈ range m, (Real.log (1 + η / ((k : ℝ) + 1))
      - Real.log (1 + u / ((k : ℝ) + 1)) - Real.log (1 + v / ((k : ℝ) + 1))) := by
  unfold R
  rw [Real.log_prod]
  · refine Finset.sum_congr rfl fun k _ => ?_
    have a := fac_pos hu k; have b := fac_pos hv k; have c := fac_pos hη k
    rw [Real.log_div c.ne' (by positivity), Real.log_mul a.ne' b.ne']
    ring
  · intro k _
    have a := fac_pos hu k; have b := fac_pos hv k; have c := fac_pos hη k
    positivity

lemma R_tendsto {η u v : ℝ} (huv : u + v = η) (hu : |u| < 1) (hv : |v| < 1)
    (hη : |η| < 1) (hη1 : -1 < η) :
    Tendsto (fun n : ℕ => R η u v (n + 1)) atTop
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
  rw [gammaSeq_ratio huv hu hv hη hn]
  field_simp


lemma abs_mono_bd {a b x B : ℝ} (ha : |a| ≤ B) (hb : |b| ≤ B) (hx : |x| ≤ 1) (m : ℕ) :
    |a * b * x ^ m| ≤ B ^ 2 := by
  rw [abs_mul, abs_mul, abs_pow]
  have h1 : |x| ^ m ≤ 1 := pow_le_one₀ (abs_nonneg _) hx
  have h2 : |a| * |b| ≤ B * B := mul_le_mul ha hb (abs_nonneg _) ((abs_nonneg _).trans ha)
  calc |a| * |b| * |x| ^ m ≤ (B * B) * 1 :=
        mul_le_mul h2 h1 (by positivity) (by nlinarith [abs_nonneg a])
    _ = B ^ 2 := by ring

lemma abs_coef_bd {a x B : ℝ} (ha : |a| ≤ B) (hx : |x| ≤ 1) (m : ℕ) :
    |a * x ^ (m + 2)| ≤ B * x ^ 2 := by
  rw [abs_mul, abs_pow]
  have : |x| ^ (m + 2) ≤ x ^ 2 := by
    rw [pow_add, ← sq_abs x, mul_comm]
    exact mul_le_of_le_one_right (by positivity) (pow_le_one₀ (abs_nonneg _) hx)
  exact mul_le_mul ha this (by positivity) ((abs_nonneg _).trans ha)

lemma exp_trunc {L a0 a1 a2 a3 x B K : ℝ} (hB : 0 ≤ B) (hK : 0 ≤ K)
    (ha0 : |a0| ≤ B) (ha1 : |a1| ≤ B) (ha2 : |a2| ≤ B) (ha3 : |a3| ≤ B)
    (hx : |x| ≤ 1) (hsmall : (4 * B + K) * x ^ 2 ≤ 1)
    (hL : |L - (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5)| ≤ K * x ^ 6) :
    |Real.exp L - (1 + a0 * x ^ 2 + a1 * x ^ 3 + (a2 + a0 ^ 2 / 2) * x ^ 4
        + (a3 + a0 * a1) * x ^ 5)| ≤ (K + K * (8 * B + K) + 7 * B ^ 2 + (4 * B + K) ^ 3) * x ^ 6 := by
  obtain ⟨P, hPdef⟩ : ∃ P, P = a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5 := ⟨_, rfl⟩
  rw [← hPdef] at hL
  have hx2 : 0 ≤ x ^ 2 := sq_nonneg x
  have hx21 : x ^ 2 ≤ 1 := by rw [← sq_abs]; exact pow_le_one₀ (abs_nonneg _) hx
  have hx6 : x ^ 6 ≤ x ^ 2 := by
    rw [show x ^ 6 = x ^ 2 * (x ^ 2) ^ 2 by ring]
    exact mul_le_of_le_one_right hx2 (pow_le_one₀ hx2 hx21)
  have hP : |P| ≤ 4 * B * x ^ 2 := by
    have e0 := abs_coef_bd ha0 hx 0
    have e1 := abs_coef_bd ha1 hx 1
    have e2 := abs_coef_bd ha2 hx 2
    have e3 := abs_coef_bd ha3 hx 3
    simp only [zero_add] at e0 e1 e2 e3
    have t1 := abs_add_le (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4) (a3 * x ^ 5)
    have t2 := abs_add_le (a0 * x ^ 2 + a1 * x ^ 3) (a2 * x ^ 4)
    have t3 := abs_add_le (a0 * x ^ 2) (a1 * x ^ 3)
    rw [hPdef]
    linarith
  have hLb : |L| ≤ (4 * B + K) * x ^ 2 := by
    have : |L| ≤ |P| + |L - P| := by
      have := abs_add_le P (L - P); rwa [add_sub_cancel] at this
    nlinarith
  have hL1 : |L| ≤ 1 := hLb.trans hsmall
  have he := Real.exp_bound hL1 (n := 3) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
  norm_num at he
  have hL3 : |L| ^ 3 ≤ (4 * B + K) ^ 3 * x ^ 6 := by
    calc |L| ^ 3 ≤ ((4 * B + K) * x ^ 2) ^ 3 := pow_le_pow_left₀ (abs_nonneg _) hLb 3
      _ = (4 * B + K) ^ 3 * x ^ 6 := by ring
  have hE1 : |Real.exp L - (1 + L + L ^ 2 / 2)| ≤ (4 * B + K) ^ 3 * x ^ 6 := by
    have : |L| ^ 3 * (2 / 9) ≤ |L| ^ 3 := by nlinarith [pow_nonneg (abs_nonneg L) 3]
    linarith
  have hE2 : |(L ^ 2 - P ^ 2) / 2| ≤ K * (8 * B + K) * x ^ 6 := by
    have hLP : |L + P| ≤ (8 * B + K) * x ^ 2 := by
      calc |L + P| ≤ |L| + |P| := abs_add_le _ _
        _ ≤ (8 * B + K) * x ^ 2 := by linarith
    rw [show (L ^ 2 - P ^ 2) / 2 = (L - P) * (L + P) / 2 by ring, abs_div, abs_mul,
      abs_of_pos (by norm_num : (0:ℝ) < 2)]
    have := mul_le_mul hL hLP (abs_nonneg _) (by positivity)
    have h8 : 0 ≤ K * (8 * B + K) * x ^ 6 := by positivity
    have : K * x ^ 6 * ((8 * B + K) * x ^ 2) ≤ K * (8 * B + K) * x ^ 6 := by
      have := mul_le_mul_of_nonneg_left hx21 h8
      nlinarith
    linarith
  obtain ⟨W, hW⟩ : ∃ W : ℝ, W = a1 ^ 2 / 2 + a0 * a2 + (a0 * a3 + a1 * a2) * x + (a2 ^ 2 / 2 + a1 * a3) * x ^ 2
    + a2 * a3 * x ^ 3 + a3 ^ 2 / 2 * x ^ 4 := ⟨_, rfl⟩
  have hWb : |W| ≤ 7 * B ^ 2 := by
    have b1 := abs_mono_bd ha1 ha1 hx 0
    have b2 := abs_mono_bd ha0 ha2 hx 0
    have b3 := abs_mono_bd ha0 ha3 hx 1
    have b4 := abs_mono_bd ha1 ha2 hx 1
    have b5 := abs_mono_bd ha2 ha2 hx 2
    have b6 := abs_mono_bd ha1 ha3 hx 2
    have b7 := abs_mono_bd ha2 ha3 hx 3
    have b8 := abs_mono_bd ha3 ha3 hx 4
    have eW : W = (a1 * a1 * x ^ 0) / 2 + a0 * a2 * x ^ 0 + a0 * a3 * x ^ 1 + a1 * a2 * x ^ 1
        + (a2 * a2 * x ^ 2) / 2 + a1 * a3 * x ^ 2 + a2 * a3 * x ^ 3 + (a3 * a3 * x ^ 4) / 2 := by
      rw [hW]; ring
    rw [eW]
    rw [abs_le] at b1 b2 b3 b4 b5 b6 b7 b8 ⊢
    constructor <;> linarith [b1.1, b1.2, b2.1, b2.2, b3.1, b3.2, b4.1, b4.2, b5.1, b5.2,
      b6.1, b6.2, b7.1, b7.2, b8.1, b8.2]
  have hE3 : |x ^ 6 * W| ≤ 7 * B ^ 2 * x ^ 6 := by
    rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 6), mul_comm]
    exact mul_le_mul_of_nonneg_right hWb (by positivity)
  have key : Real.exp L - (1 + a0 * x ^ 2 + a1 * x ^ 3 + (a2 + a0 ^ 2 / 2) * x ^ 4
        + (a3 + a0 * a1) * x ^ 5) =
      (Real.exp L - (1 + L + L ^ 2 / 2)) + (L - P) + (L ^ 2 - P ^ 2) / 2 + x ^ 6 * W := by
    rw [hW, hPdef]; ring
  rw [key]
  have t1 := abs_add_le ((Real.exp L - (1 + L + L ^ 2 / 2)) + (L - P) + (L ^ 2 - P ^ 2) / 2)
    (x ^ 6 * W)
  have t2 := abs_add_le ((Real.exp L - (1 + L + L ^ 2 / 2)) + (L - P)) ((L ^ 2 - P ^ 2) / 2)
  have t3 := abs_add_le (Real.exp L - (1 + L + L ^ 2 / 2)) (L - P)
  have hk : 0 ≤ K * (8 * B + K) * x ^ 6 := by positivity
  linarith



lemma z_two : zz 2 = π ^ 2 / 6 := by
  unfold zz
  rw [show ((2 : ℕ) : ℂ) = 2 by norm_num, riemannZeta_two]
  rw [show (π : ℂ) ^ 2 / 6 = ((π ^ 2 / 6 : ℝ) : ℂ) by push_cast; ring, Complex.ofReal_re]

lemma z_four : zz 4 = π ^ 4 / 90 := by
  unfold zz
  rw [show ((4 : ℕ) : ℂ) = 4 by norm_num, riemannZeta_four]
  rw [show (π : ℂ) ^ 4 / 90 = ((π ^ 4 / 90 : ℝ) : ℂ) by push_cast; ring, Complex.ofReal_re]

lemma z_bd {n : ℕ} (hn : 2 ≤ n) : |zz n| ≤ 2 := by
  have h := z_hasSum (k := n) (by omega)
  have h2 := z_hasSum (k := 2) (by omega)
  have hnn : 0 ≤ zz n := h.nonneg fun m => by positivity
  have hle : zz n ≤ zz 2 := hasSum_le (fun m => by
    have h1 : (1 : ℝ) ≤ (m : ℝ) + 1 := by have := m.cast_nonneg (α := ℝ); linarith
    exact one_div_le_one_div_of_le (by positivity) (pow_le_pow_right₀ h1 hn)) h h2
  rw [z_two] at hle
  have := Real.pi_lt_d2
  have : π ^ 2 / 6 ≤ 2 := by nlinarith [Real.pi_pos]
  rw [abs_of_nonneg hnn]; linarith

/-- coefficient `a_k` of the one-variable series -/
noncomputable def cf (α β γ : ℝ) (k : ℕ) : ℝ :=
  (-1 : ℝ) ^ (k + 2) * zz (k + 2) / ((k + 2 : ℕ) : ℝ) * (α ^ (k + 2) + β ^ (k + 2) - γ ^ (k + 2))

lemma cf_bd {α β γ : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2) (k : ℕ) :
    |cf α β γ k| ≤ 3 * (3 / 2) ^ (k + 2) := by
  unfold cf
  have hz := z_bd (n := k + 2) (by omega)
  have hn : (2 : ℝ) ≤ ((k + 2 : ℕ) : ℝ) := by push_cast; have := k.cast_nonneg (α := ℝ); linarith
  have hb : |α ^ (k + 2) + β ^ (k + 2) - γ ^ (k + 2)| ≤ 3 * (3 / 2) ^ (k + 2) := by
    have a := pow_le_pow_left₀ (abs_nonneg α) hα (k + 2)
    have b := pow_le_pow_left₀ (abs_nonneg β) hβ (k + 2)
    have c := pow_le_pow_left₀ (abs_nonneg γ) hγ (k + 2)
    rw [← abs_pow] at a b c
    have t1 := abs_sub (α ^ (k + 2) + β ^ (k + 2)) (γ ^ (k + 2))
    have t2 := abs_add_le (α ^ (k + 2)) (β ^ (k + 2))
    linarith
  rw [abs_mul, abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
    abs_of_pos (by linarith : (0:ℝ) < ((k + 2 : ℕ) : ℝ))]
  have : |zz (k + 2)| / ((k + 2 : ℕ) : ℝ) ≤ 1 := by
    rw [div_le_one (by linarith)]; linarith
  calc |zz (k + 2)| / ((k + 2 : ℕ) : ℝ) * |α ^ (k + 2) + β ^ (k + 2) - γ ^ (k + 2)|
      ≤ 1 * (3 * (3 / 2) ^ (k + 2)) := mul_le_mul this hb (abs_nonneg _) (by norm_num)
    _ = _ := one_mul _

lemma series_trunc {α β γ x L : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2)
    (hx : |x| ≤ 1 / 3)
    (hL : HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * zz (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) L) :
    |L - (cf α β γ 0 * x ^ 2 + cf α β γ 1 * x ^ 3 + cf α β γ 2 * x ^ 4 + cf α β γ 3 * x ^ 5)|
      ≤ 69 * x ^ 6 := by
  have hf : (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * zz (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) =
      fun k => cf α β γ k * x ^ (k + 2) := by
    ext k; unfold cf; rw [mul_pow, mul_pow, mul_pow]; ring
  rw [hf] at hL
  set r := 3 / 2 * |x| with hr
  have hr0 : 0 ≤ r := by positivity
  have hr1 : r ≤ 1 / 2 := by rw [hr]; linarith
  have hbd : ∀ k, |cf α β γ k * x ^ (k + 2)| ≤ 3 * r ^ (k + 2) := by
    intro k
    rw [abs_mul, abs_pow, hr, mul_pow]
    have := cf_bd hα hβ hγ k
    calc |cf α β γ k| * |x| ^ (k + 2) ≤ 3 * (3 / 2) ^ (k + 2) * |x| ^ (k + 2) :=
          mul_le_mul_of_nonneg_right this (by positivity)
      _ = _ := by ring
  have hg := (hasSum_nat_add_iff' 4).mpr hL
  have hgeo : HasSum (fun k : ℕ => 3 * r ^ 6 * r ^ k) (3 * r ^ 6 * (1 - r)⁻¹) :=
    (hasSum_geometric_of_lt_one hr0 (by linarith)).mul_left _
  have hb1 : ∀ k : ℕ, cf α β γ (k + 4) * x ^ (k + 4 + 2) ≤ 3 * r ^ 6 * r ^ k := by
    intro k
    have := (abs_le.mp (hbd (k + 4))).2
    calc _ ≤ 3 * r ^ (k + 4 + 2) := this
      _ = _ := by ring
  have hb2 : ∀ k : ℕ, -(cf α β γ (k + 4) * x ^ (k + 4 + 2)) ≤ 3 * r ^ 6 * r ^ k := by
    intro k
    have := (abs_le.mp (hbd (k + 4))).1
    calc _ ≤ 3 * r ^ (k + 4 + 2) := by linarith
      _ = _ := by ring
  have u1 := hasSum_le hb1 hg hgeo
  have u2 := hasSum_le hb2 hg.neg hgeo
  have hinv : 3 * r ^ 6 * (1 - r)⁻¹ ≤ 6 * r ^ 6 := by
    rw [mul_inv_le_iff₀ (by linarith)]
    nlinarith [pow_nonneg hr0 6]
  have hr6 : 6 * r ^ 6 ≤ 69 * x ^ 6 := by
    rw [hr, mul_pow, pow_abs, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 6)]
    nlinarith [pow_nonneg (sq_nonneg x) 3]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add] at u1 u2
  rw [abs_le]
  constructor <;> norm_num at u1 u2 ⊢ <;> linarith

lemma A_expand {α β γ x L : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2)
    (hx : |x| ≤ 1 / 20)
    (hL : HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * zz (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) L) :
    |Real.exp L - (1 + cf α β γ 0 * x ^ 2 + cf α β γ 1 * x ^ 3
        + (cf α β γ 2 + cf α β γ 0 ^ 2 / 2) * x ^ 4
        + (cf α β γ 3 + cf α β γ 0 * cf α β γ 1) * x ^ 5)|
      ≤ (69 + 69 * (8 * 23 + 69) + 7 * 23 ^ 2 + (4 * 23 + 69) ^ 3) * x ^ 6 := by
  have h := series_trunc hα hβ hγ (by linarith) hL
  have hB : ∀ k ≤ 3, |cf α β γ k| ≤ 23 := by
    intro k hk
    have := cf_bd hα hβ hγ k
    have : (3 / 2 : ℝ) ^ (k + 2) ≤ (3 / 2) ^ 5 := pow_le_pow_right₀ (by norm_num) (by omega)
    norm_num at this ⊢; linarith
  have hx2 : x ^ 2 ≤ 1 / 400 := by
    have := pow_le_pow_left₀ (abs_nonneg x) hx 2
    rw [sq_abs] at this; norm_num at this ⊢; linarith
  exact exp_trunc (by norm_num) (by norm_num) (hB 0 (by norm_num)) (hB 1 (by norm_num))
    (hB 2 (by norm_num)) (hB 3 (by norm_num)) (by linarith) (by nlinarith) h



lemma logG_hasSum {u v : ℝ} (hu0 : 0 < u) (hv0 : 0 < v) (huv1 : u + v < 1 / 2) :
    0 < Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + (u + v)) ∧
    HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * zz (k + 2) / ((k + 2 : ℕ) : ℝ) *
        (u ^ (k + 2) + v ^ (k + 2) - (u + v) ^ (k + 2)))
      (Real.log (Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + (u + v)))) := by
  set η := u + v with hη_def
  have huv : u + v = η := rfl
  have hu : |u| < 3 / 4 := by rw [abs_lt]; constructor <;> linarith
  have hv : |v| < 3 / 4 := by rw [abs_lt]; constructor <;> linarith
  have hu' : |u| < 1 := by linarith
  have hv' : |v| < 1 := by linarith
  have hηa : |η| < 1 / 2 := by rw [abs_lt]; constructor <;> linarith
  have hη' : |η| < 1 := by linarith
  have hη1 : -1 < η := by linarith
  have hApos : 0 < Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + η) := by
    have := Real.Gamma_pos_of_pos (show 0 < 1 + u by linarith)
    have := Real.Gamma_pos_of_pos (show 0 < 1 + v by linarith)
    have := Real.Gamma_pos_of_pos (show 0 < 1 + η by linarith)
    positivity
  refine ⟨hApos, ?_⟩
  set f : ℕ × ℕ → ℝ := fun p => (-1 : ℝ) ^ (p.2 + 2) * (1 / ((p.1 : ℝ) + 1) ^ (p.2 + 2)) /
      ((p.2 + 2 : ℕ) : ℝ) * (u ^ (p.2 + 2) + v ^ (p.2 + 2) - η ^ (p.2 + 2)) with hfdef
  have hfs : Summable f := by
    have hg1 : Summable (fun k : ℕ => 1 / ((k : ℝ) + 1) ^ 2) := (z_hasSum (by norm_num)).summable
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
        (pow_le_pow_left₀ (abs_nonneg η) (show |η| ≤ 3 / 4 by linarith) _).trans
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
  have hT := hf.prod_fiberwise (fun k => log_series_k huv hu' hv' hη' k)
  have hf' : HasSum (f ∘ Equiv.prodComm ℕ ℕ) (∑' p, f p) := (Equiv.prodComm ℕ ℕ).hasSum_iff.mpr hf
  have hZ := hf'.prod_fiberwise (fun j => by
    have := ((z_hasSum (show 1 < j + 2 by omega)).mul_left ((-1 : ℝ) ^ (j + 2))).div_const
      ((j + 2 : ℕ) : ℝ) |>.mul_right (u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2))
    exact this)
  have hlim1 := (tendsto_add_atTop_iff_nat 1).mpr hT.tendsto_sum_nat
  have hlim2 : Tendsto (fun n : ℕ => ∑ k ∈ range (n + 1), (Real.log (1 + η / ((k : ℝ) + 1))
      - Real.log (1 + u / ((k : ℝ) + 1)) - Real.log (1 + v / ((k : ℝ) + 1)))) atTop
      (𝓝 (Real.log (Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + η)))) := by
    have := ((Real.continuousAt_log hApos.ne').tendsto).comp
      (R_tendsto huv hu' hv' hη' hη1)
    refine this.congr fun n => ?_
    simp only [Function.comp]
    exact log_R hu' hv' hη' (n + 1)
  have heq := tendsto_nhds_unique hlim1 hlim2
  rw [heq] at hZ
  convert hZ using 1

/-- the Gamma side along the ray `y = t x` -/
lemma Gside {t x : ℝ} (ht0 : 0 < t) (ht : t ≤ 1 / 2) (hx0 : 0 < x) (hx : x ≤ 1 / 20) :
    |Real.Gamma (1 + x) * Real.Gamma (1 + t * x) / Real.Gamma (1 + (x + t * x)) -
      (1 + cf 1 t (1 + t) 0 * x ^ 2 + cf 1 t (1 + t) 1 * x ^ 3
        + (cf 1 t (1 + t) 2 + cf 1 t (1 + t) 0 ^ 2 / 2) * x ^ 4
        + (cf 1 t (1 + t) 3 + cf 1 t (1 + t) 0 * cf 1 t (1 + t) 1) * x ^ 5)|
      ≤ (69 + 69 * (8 * 23 + 69) + 7 * 23 ^ 2 + (4 * 23 + 69) ^ 3) * x ^ 6 := by
  have htx : 0 < t * x := by positivity
  obtain ⟨hpos, hL⟩ := logG_hasSum hx0 htx (by nlinarith)
  have hL' : HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * zz (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((1 * x) ^ (k + 2) + (t * x) ^ (k + 2) - ((1 + t) * x) ^ (k + 2)))
      (Real.log (Real.Gamma (1 + x) * Real.Gamma (1 + t * x) / Real.Gamma (1 + (x + t * x)))) := by
    convert hL using 3; ring
  have h := A_expand (α := 1) (β := t) (γ := 1 + t) (by rw [abs_one]; norm_num)
    (by rw [abs_of_pos ht0]; linarith) (by rw [abs_of_pos (by linarith)]; linarith)
    (by rw [abs_of_pos hx0]; exact hx) hL'
  rwa [Real.exp_log hpos] at h

/-- `B(x,y) = (1/x + 1/y) Γ(1+x)Γ(1+y)/Γ(1+x+y)` -/
lemma beta_gamma {x y : ℝ} (hx : 0 < x) (hy : 0 < y) :
    Real.Gamma x * Real.Gamma y / Real.Gamma (x + y) =
      (1 / x + 1 / y) * (Real.Gamma (1 + x) * Real.Gamma (1 + y) / Real.Gamma (1 + (x + y))) := by
  rw [add_comm 1 x, add_comm 1 y, add_comm 1 (x + y), Real.Gamma_add_one hx.ne',
    Real.Gamma_add_one hy.ne', Real.Gamma_add_one (by linarith)]
  have := Real.Gamma_pos_of_pos (show 0 < x + y by linarith)
  field_simp
  ring


end Q02

namespace Q02

lemma lin_zero {a M δ : ℝ} (hδ : 0 < δ) (h : ∀ x : ℝ, 0 < x → x ≤ δ → |a| ≤ M * x) : a = 0 := by
  by_contra ha
  have ha' : 0 < |a| := abs_pos.mpr ha
  set x := min δ (|a| / (2 * (|M| + 1))) with hx
  have hx0 : 0 < x := lt_min hδ (by positivity)
  have hx1 : x ≤ |a| / (2 * (|M| + 1)) := min_le_right _ _
  have h1 := h x hx0 (min_le_left _ _)
  have h2 : M * x ≤ |M| * x := mul_le_mul_of_nonneg_right (le_abs_self M) hx0.le
  have h3 : |M| * x ≤ |M| * (|a| / (2 * (|M| + 1))) := mul_le_mul_of_nonneg_left hx1 (abs_nonneg M)
  have h4 : |M| * (|a| / (2 * (|M| + 1))) < |a| := by
    rw [mul_div_assoc', div_lt_iff₀ (by positivity)]
    nlinarith [abs_nonneg M]
  linarith

lemma pow_bd {e x : ℝ} (hx0 : 0 < x) (hx1 : x ≤ 1) {n m : ℕ} (hnm : m ≤ n) :
    |e * x ^ n| ≤ |e| * x ^ m := by
  rw [abs_mul, abs_of_pos (pow_pos hx0 n)]
  exact mul_le_mul_of_nonneg_left (pow_le_pow_of_le_one hx0.le hx1 hnm) (abs_nonneg e)

lemma div_pow_bd {e M x : ℝ} {n : ℕ} (hx0 : 0 < x) (h : |e| * x ^ n ≤ M * x ^ (n + 1)) :
    |e| ≤ M * x := by
  have hp : 0 < x ^ n := pow_pos hx0 n
  have : |e| * x ^ n ≤ (M * x) * x ^ n := by rw [pow_succ] at h; linarith
  exact le_of_mul_le_mul_right this hp

lemma poly_zero {e1 e2 e3 e4 K δ : ℝ} (hδ : 0 < δ) (hδ1 : δ ≤ 1)
    (h : ∀ x : ℝ, 0 < x → x ≤ δ → |e1 * x + e2 * x ^ 2 + e3 * x ^ 3 + e4 * x ^ 4| ≤ K * x ^ 5) :
    e1 = 0 ∧ e2 = 0 ∧ e3 = 0 ∧ e4 = 0 := by
  have tri : ∀ a b : ℝ, |a| ≤ |a + b| + |b| := fun a b => by
    have := abs_sub (a + b) b; rwa [add_sub_cancel_right] at this
  have hK : ∀ x : ℝ, 0 < x → K * x ^ 5 ≤ |K| * x ^ 5 := fun x hx =>
    mul_le_mul_of_nonneg_right (le_abs_self K) (pow_pos hx 5).le
  have h1 : e1 = 0 := by
    refine lin_zero hδ (M := |K| + |e2| + |e3| + |e4|) (fun x hx hxd => ?_)
    have hx1 : x ≤ 1 := hxd.trans hδ1
    have H := h x hx hxd
    have t := tri (e1 * x ^ 1) (e2 * x ^ 2 + e3 * x ^ 3 + e4 * x ^ 4)
    have b2 := pow_bd (e := e2) hx hx1 (show 2 ≤ 2 by norm_num)
    have b3 := pow_bd (e := e3) hx hx1 (show 2 ≤ 3 by norm_num)
    have b4 := pow_bd (e := e4) hx hx1 (show 2 ≤ 4 by norm_num)
    have b5 := pow_bd (e := |K|) hx hx1 (show 2 ≤ 5 by norm_num)
    rw [abs_abs, abs_of_nonneg (by have := pow_pos hx 5; positivity)] at b5
    have tt := abs_add_le (e2 * x ^ 2 + e3 * x ^ 3) (e4 * x ^ 4)
    have tt2 := abs_add_le (e2 * x ^ 2) (e3 * x ^ 3)
    have hx5 := hK x hx
    rw [show e1 * x ^ 1 + (e2 * x ^ 2 + e3 * x ^ 3 + e4 * x ^ 4) =
      e1 * x + e2 * x ^ 2 + e3 * x ^ 3 + e4 * x ^ 4 by ring, abs_mul, pow_one,
      abs_of_pos hx] at t
    apply div_pow_bd (n := 1) hx
    have : (|K| + |e2| + |e3| + |e4|) * x ^ (1 + 1) =
        |K| * x ^ 2 + |e2| * x ^ 2 + |e3| * x ^ 2 + |e4| * x ^ 2 := by ring
    rw [this, pow_one]
    linarith
  subst h1
  have h2 : e2 = 0 := by
    refine lin_zero hδ (M := |K| + |e3| + |e4|) (fun x hx hxd => ?_)
    have hx1 : x ≤ 1 := hxd.trans hδ1
    have H := h x hx hxd
    have t := tri (e2 * x ^ 2) (e3 * x ^ 3 + e4 * x ^ 4)
    have b3 := pow_bd (e := e3) hx hx1 (show 3 ≤ 3 by norm_num)
    have b4 := pow_bd (e := e4) hx hx1 (show 3 ≤ 4 by norm_num)
    have b5 := pow_bd (e := |K|) hx hx1 (show 3 ≤ 5 by norm_num)
    rw [abs_abs, abs_of_nonneg (by have := pow_pos hx 5; positivity)] at b5
    have tt := abs_add_le (e3 * x ^ 3) (e4 * x ^ 4)
    have hx5 := hK x hx
    rw [show e2 * x ^ 2 + (e3 * x ^ 3 + e4 * x ^ 4) =
      0 * x + e2 * x ^ 2 + e3 * x ^ 3 + e4 * x ^ 4 by ring, abs_mul,
      abs_of_pos (pow_pos hx 2)] at t
    apply div_pow_bd (n := 2) hx
    have : (|K| + |e3| + |e4|) * x ^ (2 + 1) = |K| * x ^ 3 + |e3| * x ^ 3 + |e4| * x ^ 3 := by ring
    rw [this]
    linarith
  subst h2
  have h3 : e3 = 0 := by
    refine lin_zero hδ (M := |K| + |e4|) (fun x hx hxd => ?_)
    have hx1 : x ≤ 1 := hxd.trans hδ1
    have H := h x hx hxd
    have t := tri (e3 * x ^ 3) (e4 * x ^ 4)
    have b4 := pow_bd (e := e4) hx hx1 (show 4 ≤ 4 by norm_num)
    have b5 := pow_bd (e := |K|) hx hx1 (show 4 ≤ 5 by norm_num)
    rw [abs_abs, abs_of_nonneg (by have := pow_pos hx 5; positivity)] at b5
    have hx5 := hK x hx
    rw [show e3 * x ^ 3 + e4 * x ^ 4 =
      0 * x + 0 * x ^ 2 + e3 * x ^ 3 + e4 * x ^ 4 by ring, abs_mul,
      abs_of_pos (pow_pos hx 3)] at t
    apply div_pow_bd (n := 3) hx
    have : (|K| + |e4|) * x ^ (3 + 1) = |K| * x ^ 4 + |e4| * x ^ 4 := by ring
    rw [this]
    linarith
  subst h3
  have h4 : e4 = 0 := by
    refine lin_zero hδ (M := |K|) (fun x hx hxd => ?_)
    have hx1 : x ≤ 1 := hxd.trans hδ1
    have H := h x hx hxd
    have hx5 := hK x hx
    rw [show 0 * x + 0 * x ^ 2 + 0 * x ^ 3 + e4 * x ^ 4 = e4 * x ^ 4 by ring, abs_mul,
      abs_of_pos (pow_pos hx 4)] at H
    exact div_pow_bd (n := 4) hx (by rw [show 4 + 1 = 5 from rfl]; linarith)
  exact ⟨rfl, rfl, rfl, h4⟩

/-- the coefficients `e_1, …, e_4` of the comparison along `y = t x` (`c = 1 + t`) -/
noncomputable def ecoef (t : ℝ) (d : ℕ) : ℝ :=
  let c := 1 + t
  let g0 := cf 1 t c 0
  let g1 := cf 1 t c 1
  let g2 := cf 1 t c 2 + g0 ^ 2 / 2
  let g3 := cf 1 t c 3 + g0 * g1
  match d with
  | 1 => -2 * c * J 1 0 - (1 + 1 / t) * g0
  | 2 => c * (1 + t) * J 1 1 + c ^ 2 * J 2 0 - (1 + 1 / t) * g1
  | 3 => -(c * (1 + t ^ 2) / 2 * J 1 2 + c ^ 2 * (1 + t) / 2 * J 2 1 + c ^ 3 / 3 * J 3 0)
      - (1 + 1 / t) * g2
  | _ => c * (1 + t ^ 3) / 6 * J 1 3 + c ^ 2 * (1 + t ^ 2) / 4 * J 2 2
      + c ^ 3 * (1 + t) / 6 * J 3 1 + c ^ 4 / 12 * J 4 0 - (1 + 1 / t) * g3

lemma Pint_expand (t x : ℝ) :
    Pint x (x + t * x) + Pint (t * x) (t * x + x) -
      (1 + 1 / t) * (cf 1 t (1 + t) 0 * x + cf 1 t (1 + t) 1 * x ^ 2
        + (cf 1 t (1 + t) 2 + cf 1 t (1 + t) 0 ^ 2 / 2) * x ^ 3
        + (cf 1 t (1 + t) 3 + cf 1 t (1 + t) 0 * cf 1 t (1 + t) 1) * x ^ 4) =
    ecoef t 1 * x + ecoef t 2 * x ^ 2 + ecoef t 3 * x ^ 3 + ecoef t 4 * x ^ 4 := by
  simp only [Pint, idx, cc, ecoef]
  simp [Finset.sum_insert, Nat.factorial]
  ring

lemma ecoef_zero {t : ℝ} (ht0 : 0 < t) (ht : t ≤ 1 / 2) :
    ecoef t 1 = 0 ∧ ecoef t 2 = 0 ∧ ecoef t 3 = 0 ∧ ecoef t 4 = 0 := by
  set KJ := ∫ s in Ioi (0 : ℝ), gmaj s
  set C1 : ℝ := 69 + 69 * (8 * 23 + 69) + 7 * 23 ^ 2 + (4 * 23 + 69) ^ 3
  refine poly_zero (δ := 1 / 20) (K := 2 * KJ + (1 + 1 / t) * C1) (by norm_num) (by norm_num)
    (fun x hx hxd => ?_)
  have htx : 0 < t * x := by positivity
  -- the Beta identity
  have hb := beta_real hx htx
  rw [beta_split hx htx, add_comm (t * x) x] at hb
  have hG0 : 0 < Real.Gamma (x + t * x) := Real.Gamma_pos_of_pos (by linarith)
  have hB : Real.Gamma x * Real.Gamma (t * x) / Real.Gamma (x + t * x) =
      1 / x + 1 / (t * x) + Phi x (x + t * x) + Phi (t * x) (t * x + x) := by
    rw [div_eq_iff hG0.ne', ← hb]; ring_nf
  rw [beta_gamma hx htx] at hB
  set G := Real.Gamma (1 + x) * Real.Gamma (1 + t * x) / Real.Gamma (1 + (x + t * x))
  have hPhi : Phi x (x + t * x) + Phi (t * x) (t * x + x) = (1 + 1 / t) / x * (G - 1) := by
    have : (1 / x + 1 / (t * x)) * G = 1 / x + 1 / (t * x) + Phi x (x + t * x)
        + Phi (t * x) (t * x + x) := hB
    field_simp at this ⊢
    linarith
  have hGs := Gside ht0 ht hx hxd
  have hx4 : x ≤ 1 / 4 := by linarith
  have p1 := phi_approx (a := x) (b := x + t * x) hx4 hx le_rfl (by positivity) (by nlinarith)
  have p2 := phi_approx (a := t * x) (b := t * x + x) hx4 htx (by nlinarith) (by positivity)
    (by nlinarith)
  rw [← Pint_expand t x]
  set P := cf 1 t (1 + t) 0 * x + cf 1 t (1 + t) 1 * x ^ 2
        + (cf 1 t (1 + t) 2 + cf 1 t (1 + t) 0 ^ 2 / 2) * x ^ 3
        + (cf 1 t (1 + t) 3 + cf 1 t (1 + t) 0 * cf 1 t (1 + t) 1) * x ^ 4 with hP
  have hGP : (1 + 1 / t) / x * (G - 1) - (1 + 1 / t) * P =
      (1 + 1 / t) / x * (G - (1 + cf 1 t (1 + t) 0 * x ^ 2 + cf 1 t (1 + t) 1 * x ^ 3
        + (cf 1 t (1 + t) 2 + cf 1 t (1 + t) 0 ^ 2 / 2) * x ^ 4
        + (cf 1 t (1 + t) 3 + cf 1 t (1 + t) 0 * cf 1 t (1 + t) 1) * x ^ 5)) := by
    rw [hP]; field_simp; ring
  have hc : 0 < (1 + 1 / t) / x := by positivity
  have q1 : |(1 + 1 / t) / x * (G - 1) - (1 + 1 / t) * P| ≤ (1 + 1 / t) * C1 * x ^ 5 := by
    rw [hGP, abs_mul, abs_of_pos hc]
    calc (1 + 1 / t) / x * |_| ≤ (1 + 1 / t) / x * (C1 * x ^ 6) :=
          mul_le_mul_of_nonneg_left hGs hc.le
      _ = (1 + 1 / t) * C1 * x ^ 5 := by field_simp
  have split : Pint x (x + t * x) + Pint (t * x) (t * x + x) - (1 + 1 / t) * P =
      -(Phi x (x + t * x) - Pint x (x + t * x)) - (Phi (t * x) (t * x + x) - Pint (t * x) (t * x + x))
      + ((1 + 1 / t) / x * (G - 1) - (1 + 1 / t) * P) := by
    rw [← hPhi]; ring
  rw [split]
  have t1 := abs_add_le (-(Phi x (x + t * x) - Pint x (x + t * x)) -
    (Phi (t * x) (t * x + x) - Pint (t * x) (t * x + x)))
    ((1 + 1 / t) / x * (G - 1) - (1 + 1 / t) * P)
  have t2 := abs_sub (-(Phi x (x + t * x) - Pint x (x + t * x)))
    (Phi (t * x) (t * x + x) - Pint (t * x) (t * x + x))
  rw [abs_neg] at t2
  have : (2 * KJ + (1 + 1 / t) * C1) * x ^ 5 = x ^ 5 * KJ + x ^ 5 * KJ + (1 + 1 / t) * C1 * x ^ 5 := by
    ring
  rw [this]
  linarith

lemma zz_three : zz 3 = X02.z3 := by unfold zz X02.z3; norm_num
lemma zz_five : zz 5 = z5 := by unfold zz z5; norm_num

lemma C45 (hJ12 : J 1 2 = 7 * Real.pi ^ 4 / 360) (hJ13 : J 1 3 = 45 * z5 / 8) :
    J 3 0 / 12 + J 2 1 / 8 = Real.pi ^ 4 / 2880 ∧
    (J 4 0 + 2 * J 3 1) / 48 = -(X02.z3 * Real.pi ^ 2 / 48) + 31 * z5 / 128 := by
  obtain ⟨-, -, a3, a4⟩ := ecoef_zero (t := 1 / 2) (by norm_num) le_rfl
  obtain ⟨-, -, -, b4⟩ := ecoef_zero (t := 1 / 4) (by norm_num) (by norm_num)
  simp only [ecoef, cf] at a3 a4 b4
  norm_num [z_two, z_four, zz_three, zz_five] at a3 a4 b4
  constructor
  · rw [hJ12] at a3
    linear_combination (-(2 : ℝ) / 27) * a3
  · rw [hJ13] at a4 b4
    linear_combination (1 / 48 : ℝ) * ((272 / 21) * a4 - (768 / 35) * b4)

end Q02

/-- Lemma `lem:pip:C45`, from the Beta function split at `1`, `J_{1,2}` and `J_{1,3}`. -/
theorem Q02.proof_3 : Q02.stmt_3 := Q02.C45 Q02.proof_2.2.1 Q02.proof_2.2.2

/-- `N_4` from `C_4` (Lemma `lem:pip:C45`, first identity). -/
theorem Q02.proof_10 : Q02.stmt_10 := Q02.N4_of Q02.proof_3.1

/-- `N_5` from `C_4`, `C_5` (Lemma `lem:pip:C45`) and `J_{1,2}`, `J_{1,3}`. -/
theorem Q02.proof_11 : Q02.stmt_11 := Q02.N5_of Q02.proof_3.1 Q02.proof_3.2

namespace Q02


/-! ### Lemma `lem:pip:expansion`: `Q(h) = Σ_{n ≤ K} q_n h^n + O(h^{K+1})` -/

lemma intOn_onepow (N : ℕ) : IntegrableOn (fun s => (1 + s) ^ N * Real.exp (-s)) (Ioi 0) := by
  have : (fun s : ℝ => (1 + s) ^ N * Real.exp (-s)) = fun s =>
      ∑ k ∈ Finset.range (N + 1), (N.choose k : ℝ) * (s ^ k * Real.exp (-(1 * s))) := by
    funext s
    rw [add_comm, add_pow, Finset.sum_mul]
    refine Finset.sum_congr rfl (fun k _ => ?_)
    rw [one_pow, one_mul]; ring
  rw [this]
  exact integrable_finsetSum _ (fun k _ => (X02.intOn_pow_exp k one_pos).const_mul _)

lemma Wpoly_bd (m : ℕ) : ∃ D : ℝ, 0 ≤ D ∧ ∀ s : ℝ, 0 ≤ s → |Wpoly m s| ≤ D * (1 + s) ^ m := by
  refine ⟨∑ j ∈ Finset.range (m + 1), |wcoef (m - j) j|, Finset.sum_nonneg (fun _ _ => abs_nonneg _),
    fun s hs => ?_⟩
  rw [Wpoly, Finset.sum_mul]
  refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum (fun j hj => ?_))
  have hjm : j ≤ m := by simp at hj; omega
  rw [abs_mul, abs_of_nonneg (pow_nonneg hs j)]
  exact mul_le_mul_of_nonneg_left
    ((pow_le_pow_left₀ hs (by linarith) j).trans (pow_le_pow_right₀ (by linarith) hjm))
    (abs_nonneg _)

lemma P_bd (K : ℕ) : ∃ D : ℝ, 0 ≤ D ∧ ∀ h s : ℝ, 0 ≤ h → h ≤ 1 → 0 ≤ s →
    |∑ m ∈ Finset.range (K + 1), h ^ m * Wpoly m s| ≤ D * (1 + s) ^ K := by
  have := fun m : ℕ => Wpoly_bd m
  choose D hD0 hD using this
  refine ⟨∑ m ∈ Finset.range (K + 1), D m, Finset.sum_nonneg (fun m _ => hD0 m),
    fun h s hh h1 hs => ?_⟩
  rw [Finset.sum_mul]
  refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum (fun m hm => ?_))
  have hmK : m ≤ K := by simp at hm; omega
  rw [abs_mul, abs_of_nonneg (pow_nonneg hh m)]
  calc h ^ m * |Wpoly m s| ≤ 1 * (D m * (1 + s) ^ m) :=
        mul_le_mul (pow_le_one₀ hh h1) (hD m s hs) (abs_nonneg _) zero_le_one
    _ ≤ D m * (1 + s) ^ K := by
        rw [one_mul]
        exact mul_le_mul_of_nonneg_left (pow_le_pow_right₀ (by linarith) hmK) (hD0 m)

/-- the Taylor bound of `W(h,sh)` for all `s ≥ 0` (polynomial weight) -/
lemma W_global_K (K : ℕ) : ∃ C h₀ : ℝ, 0 < h₀ ∧ h₀ ≤ 1 / 4 ∧ 0 ≤ C ∧ ∀ h s : ℝ, 0 < h → h ≤ h₀ →
    0 ≤ s → |X02.W h (s * h) - ∑ m ∈ Finset.range (K + 1), h ^ m * Wpoly m s|
      ≤ C * h ^ (K + 1) * (1 + s) ^ (3 * K + 3) := by
  obtain ⟨C0, h0, h0pos, hT⟩ := W_taylor K
  obtain ⟨D, hD0, hD⟩ := P_bd K
  refine ⟨|C0| + 1 + D, min h0 (1 / 4), lt_min h0pos (by norm_num), min_le_right _ _,
    by positivity, fun h s hh hh0 hs => ?_⟩
  have hh1 : h ≤ h0 := hh0.trans (min_le_left _ _)
  have hh4 : h ≤ 1 / 4 := hh0.trans (min_le_right _ _)
  have hs1 : 1 ≤ 1 + s := by linarith
  have hpow : (1 + s) ^ (K + 1) ≤ (1 + s) ^ (3 * K + 3) := pow_le_pow_right₀ hs1 (by omega)
  have hpowK : (1 + s) ^ K ≤ (1 + s) ^ K := le_rfl
  rcases le_or_gt s (1 / Real.sqrt h) with hsl | hsl
  · have := hT h s hh hh1 hs hsl
    rw [mul_pow] at this
    calc _ ≤ C0 * (h ^ (K + 1) * (1 + s) ^ (K + 1)) := this
      _ ≤ |C0| * (h ^ (K + 1) * (1 + s) ^ (3 * K + 3)) :=
          mul_le_mul (le_abs_self _) (mul_le_mul_of_nonneg_left hpow (by positivity))
            (by positivity) (abs_nonneg _)
      _ ≤ (|C0| + 1 + D) * h ^ (K + 1) * (1 + s) ^ (3 * K + 3) := by
          have : 0 ≤ h ^ (K + 1) * (1 + s) ^ (3 * K + 3) := by positivity
          nlinarith
  · -- the tail `s > h^{-1/2}`: `1 < s² h`
    have hsq : 0 < Real.sqrt h := Real.sqrt_pos.2 hh
    have h1 : 1 < s * Real.sqrt h := by rwa [div_lt_iff₀ hsq] at hsl
    have h2 : 1 ≤ (1 + s) ^ 2 * h := by
      have : (s * Real.sqrt h) ^ 2 = s ^ 2 * h := by rw [mul_pow, Real.sq_sqrt hh.le]
      nlinarith
    have h3 : 1 ≤ (1 + s) ^ (2 * K + 2) * h ^ (K + 1) := by
      have := one_le_pow₀ (n := K + 1) h2
      calc 1 ≤ ((1 + s) ^ 2 * h) ^ (K + 1) := this
        _ = _ := by rw [mul_pow, ← pow_mul]; ring_nf
    obtain ⟨w0, w1⟩ := X02.W_pos_le hh.le (by linarith) hs
    have hP := hD h s hh.le (by linarith) hs
    have hW : |X02.W h (s * h) - ∑ m ∈ Finset.range (K + 1), h ^ m * Wpoly m s|
        ≤ (1 + D) * (1 + s) ^ K := by
      have := abs_sub (X02.W h (s * h)) (∑ m ∈ Finset.range (K + 1), h ^ m * Wpoly m s)
      rw [abs_of_pos w0] at this
      have : (1 : ℝ) ≤ (1 + s) ^ K := one_le_pow₀ hs1
      nlinarith
    have hK3 : (1 + s) ^ (3 * K + 3) = (1 + s) ^ K * (1 + s) ^ (2 * K + 2) * (1 + s) := by
      rw [← pow_add, ← pow_succ]; congr 1; omega
    calc _ ≤ (1 + D) * (1 + s) ^ K := hW
      _ ≤ (1 + D) * (1 + s) ^ K * ((1 + s) ^ (2 * K + 2) * h ^ (K + 1)) := by
          have : 0 ≤ (1 + D) * (1 + s) ^ K := by positivity
          nlinarith
      _ ≤ (|C0| + 1 + D) * h ^ (K + 1) * (1 + s) ^ (3 * K + 3) := by
          rw [hK3]
          have : 0 ≤ (1 + s) ^ K * (1 + s) ^ (2 * K + 2) * h ^ (K + 1) := by positivity
          have hC : 1 + D ≤ |C0| + 1 + D := by linarith [abs_nonneg C0]
          calc (1 + D) * (1 + s) ^ K * ((1 + s) ^ (2 * K + 2) * h ^ (K + 1))
              = (1 + D) * ((1 + s) ^ K * (1 + s) ^ (2 * K + 2) * h ^ (K + 1)) * 1 := by ring
            _ ≤ (|C0| + 1 + D) * ((1 + s) ^ K * (1 + s) ^ (2 * K + 2) * h ^ (K + 1)) * (1 + s) :=
                mul_le_mul (mul_le_mul_of_nonneg_right hC this) hs1 zero_le_one (by positivity)
            _ = _ := by ring

/-- the product of the two truncations, minus its part of degree `≤ K` -/
lemma prod_trunc (K : ℕ) {h : ℝ} (hh : 0 ≤ h) (h1 : h ≤ 1) (a b : ℕ → ℝ) :
    |(∑ m ∈ Finset.range (K + 1), h ^ m * a m) * (∑ r ∈ Finset.Icc 1 K, h ^ r * b r) -
      ∑ n ∈ Finset.Icc 1 K, h ^ n * ∑ r ∈ Finset.Icc 1 n, b r * a (n - r)| ≤
      h ^ (K + 1) * ∑ r ∈ Finset.Icc 1 K, ∑ m ∈ Finset.range (K + 1), |a m| * |b r| := by
  have e1 : (∑ m ∈ Finset.range (K + 1), h ^ m * a m) * (∑ r ∈ Finset.Icc 1 K, h ^ r * b r) =
      ∑ r ∈ Finset.Icc 1 K, ∑ m ∈ Finset.range (K + 1), h ^ (m + r) * a m * b r := by
    rw [Finset.sum_mul_sum, Finset.sum_comm]
    refine Finset.sum_congr rfl (fun r _ => Finset.sum_congr rfl (fun m _ => ?_))
    rw [pow_add]; ring
  have e2 : ∑ n ∈ Finset.Icc 1 K, h ^ n * ∑ r ∈ Finset.Icc 1 n, b r * a (n - r) =
      ∑ r ∈ Finset.Icc 1 K, ∑ m ∈ Finset.range (K + 1 - r), h ^ (m + r) * a m * b r := by
    simp_rw [Finset.mul_sum]
    rw [Finset.sum_comm' (t' := Finset.Icc 1 K) (s' := fun r => Finset.Ico r (K + 1))
      (by intro n r; simp only [Finset.mem_Icc, Finset.mem_Ico]; omega)]
    refine Finset.sum_congr rfl (fun r hr => ?_)
    rw [Finset.sum_Ico_eq_sum_range]
    refine Finset.sum_congr rfl (fun m _ => ?_)
    rw [show r + m - r = m by omega, show r + m = m + r by omega]; ring
  have e3 : ∀ r ∈ Finset.Icc 1 K, ∑ m ∈ Finset.range (K + 1), h ^ (m + r) * a m * b r =
      ∑ m ∈ Finset.range (K + 1 - r), h ^ (m + r) * a m * b r +
        ∑ m ∈ Finset.Ico (K + 1 - r) (K + 1), h ^ (m + r) * a m * b r := by
    intro r _
    rw [Finset.sum_range_add_sum_Ico _ (by omega)]
  rw [e1, e2, Finset.sum_congr rfl e3, Finset.sum_add_distrib, add_sub_cancel_left,
    Finset.mul_sum]
  refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum (fun r hr => ?_))
  refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
  rw [Finset.mul_sum]
  calc ∑ m ∈ Finset.Ico (K + 1 - r) (K + 1), |h ^ (m + r) * a m * b r|
      ≤ ∑ m ∈ Finset.Ico (K + 1 - r) (K + 1), h ^ (K + 1) * (|a m| * |b r|) := by
        refine Finset.sum_le_sum (fun m hm => ?_)
        have hmr : K + 1 ≤ m + r := by simp at hm hr; omega
        rw [abs_mul, abs_mul, abs_of_nonneg (pow_nonneg hh _)]
        have : h ^ (m + r) ≤ h ^ (K + 1) := pow_le_pow_of_le_one hh h1 hmr
        calc h ^ (m + r) * |a m| * |b r| = h ^ (m + r) * (|a m| * |b r|) := by ring
          _ ≤ h ^ (K + 1) * (|a m| * |b r|) := mul_le_mul_of_nonneg_right this (by positivity)
    _ ≤ ∑ m ∈ Finset.range (K + 1), h ^ (K + 1) * (|a m| * |b r|) :=
        Finset.sum_le_sum_of_subset_of_nonneg (fun m hm => by simp at hm ⊢; omega)
          (fun m _ _ => by positivity)


lemma exp_rem_pos {x : ℝ} (hx : 0 ≤ x) (n : ℕ) :
    |Real.exp x - ∑ r ∈ Finset.range n, x ^ r / (r.factorial : ℝ)| ≤ x ^ n * Real.exp x := by
  have h := Complex.norm_exp_sub_sum_le_norm_mul_exp (x : ℂ) n
  have e : Complex.exp (x : ℂ) - ∑ m ∈ Finset.range n, (x : ℂ) ^ m / (m.factorial : ℂ) =
      ((Real.exp x - ∑ m ∈ Finset.range n, x ^ m / (m.factorial : ℝ) : ℝ) : ℂ) := by
    push_cast; rfl
  rw [e, Complex.norm_real, Real.norm_eq_abs, Complex.norm_real, Real.norm_eq_abs,
    abs_of_nonneg hx] at h
  exact h

lemma range_split (K : ℕ) (f : ℕ → ℝ) :
    ∑ r ∈ Finset.range (K + 1), f r = f 0 + ∑ r ∈ Finset.Icc 1 K, f r := by
  rw [Finset.range_eq_Ico, Finset.sum_eq_sum_Ico_succ_bot (by omega)]
  congr 1

/-- the truncated integrand whose integral is `Σ_{n ≤ K} q_n h^n` -/
noncomputable def Tq (K : ℕ) (h s : ℝ) : ℝ :=
  ∑ n ∈ Finset.Icc 1 K, h ^ n *
    ∑ r ∈ Finset.Icc 1 n, X02.calL s ^ r / (r.factorial : ℝ) * Wpoly (n - r) s

lemma intOn_Tq (K : ℕ) (h : ℝ) : IntegrableOn (Tq K h) (Ioi 0) := by
  unfold Tq
  refine integrable_finsetSum _ (fun n _ => Integrable.const_mul ?_ _)
  refine integrable_finsetSum _ (fun r hr => ?_)
  have hr1 : 1 ≤ r := (Finset.mem_Icc.mp hr).1
  have : IntegrableOn (fun s => 1 / (r.factorial : ℝ) * (Wpoly (n - r) s * X02.calL s ^ r))
      (Ioi 0) := (intOn_Wpoly_calL (n - r) r hr1).const_mul (1 / (r.factorial : ℝ))
  refine this.congr_fun (fun s _ => ?_) measurableSet_Ioi
  ring

lemma integral_Tq (K : ℕ) (h : ℝ) :
    ∫ s in Ioi (0 : ℝ), Tq K h s = ∑ n ∈ Finset.Icc 1 K, q n * h ^ n := by
  have hI : ∀ n r : ℕ, r ∈ Finset.Icc 1 n → IntegrableOn
      (fun s => X02.calL s ^ r / (r.factorial : ℝ) * Wpoly (n - r) s) (Ioi 0) := by
    intro n r hr
    have hr1 : 1 ≤ r := (Finset.mem_Icc.mp hr).1
    have : IntegrableOn (fun s => 1 / (r.factorial : ℝ) * (Wpoly (n - r) s * X02.calL s ^ r))
        (Ioi 0) := (intOn_Wpoly_calL (n - r) r hr1).const_mul (1 / (r.factorial : ℝ))
    refine this.congr_fun (fun s _ => ?_) measurableSet_Ioi
    ring
  have hIn : ∀ n ∈ Finset.Icc 1 K, IntegrableOn (fun s => h ^ n *
      ∑ r ∈ Finset.Icc 1 n, X02.calL s ^ r / (r.factorial : ℝ) * Wpoly (n - r) s) (Ioi 0) :=
    fun n _ => (integrable_finsetSum _ (fun r hr => hI n r hr)).const_mul _
  have e1 : ∫ s in Ioi (0 : ℝ), Tq K h s = ∑ n ∈ Finset.Icc 1 K, ∫ s in Ioi (0 : ℝ), h ^ n *
      ∑ r ∈ Finset.Icc 1 n, X02.calL s ^ r / (r.factorial : ℝ) * Wpoly (n - r) s :=
    integral_finsetSum _ hIn
  rw [e1]
  refine Finset.sum_congr rfl (fun n _ => ?_)
  have e2 : ∫ s in Ioi (0 : ℝ), h ^ n *
      ∑ r ∈ Finset.Icc 1 n, X02.calL s ^ r / (r.factorial : ℝ) * Wpoly (n - r) s =
      h ^ n * ∑ r ∈ Finset.Icc 1 n, ∫ s in Ioi (0 : ℝ),
        X02.calL s ^ r / (r.factorial : ℝ) * Wpoly (n - r) s := by
    rw [integral_const_mul, integral_finsetSum _ (fun r hr => hI n r hr)]
  rw [e2, mul_comm (q n)]
  refine congrArg (fun z => h ^ n * z) ?_
  unfold q
  refine Finset.sum_congr rfl (fun r _ => ?_)
  rw [← integral_const_mul]
  refine setIntegral_congr_fun measurableSet_Ioi (fun s _ => ?_)
  ring

theorem Q_expand (K : ℕ)
    (hI : ∀ h : ℝ, 0 < h → h < 1 →
      IntegrableOn (fun s => X02.W h (s * h) * ((1 + Real.exp (-s)) ^ h - 1)) (Ioi 0)) :
    ∃ C h₀ : ℝ, 0 < h₀ ∧ ∀ h : ℝ, 0 < h → h ≤ h₀ →
      |X02.Q h - ∑ n ∈ Finset.Icc 1 K, q n * h ^ n| ≤ C * h ^ (K + 1) := by
  obtain ⟨C, h0, h0pos, h04, hC0, hW⟩ := W_global_K K
  obtain ⟨D, hD0, hD⟩ := P_bd K
  have hWp := fun m : ℕ => Wpoly_bd m
  choose Dm hDm0 hDm using hWp
  set Ds := ∑ m ∈ Finset.range (K + 1), Dm m with hDs
  have hDs0 : 0 ≤ Ds := Finset.sum_nonneg (fun m _ => hDm0 m)
  set M := C + 2 * D + K * Ds with hM
  have hM0 : 0 ≤ M := by positivity
  set Nn := 3 * K + 3
  refine ⟨M * ∫ s in Ioi (0 : ℝ), (1 + s) ^ Nn * Real.exp (-s), h0, h0pos, fun h hh hh0 => ?_⟩
  have hh1 : h ≤ 1 := by linarith
  have hpt : ∀ s : ℝ, 0 < s →
      |X02.W h (s * h) * ((1 + Real.exp (-s)) ^ h - 1) - Tq K h s| ≤
        h ^ (K + 1) * (M * ((1 + s) ^ Nn * Real.exp (-s))) := by
    intro s hs
    set L := X02.calL s with hLdef
    have hL0 : 0 ≤ L := X02.calL_nonneg s
    have hLe : L ≤ Real.exp (-s) := X02.calL_le s
    have he1 : Real.exp (-s) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
    have hL1 : L ≤ 1 := hLe.trans he1
    have hs1 : 1 ≤ 1 + s := by linarith
    have hE : (1 + Real.exp (-s)) ^ h = Real.exp (h * L) := by
      rw [Real.rpow_def_of_pos (by positivity), mul_comm]; rfl
    have hE1 : Real.exp (h * L) ≤ 1 + Real.exp (-s) := by
      rw [← hE]
      calc (1 + Real.exp (-s)) ^ h ≤ (1 + Real.exp (-s)) ^ (1 : ℝ) :=
            Real.rpow_le_rpow_of_exponent_le (by linarith [Real.exp_pos (-s)]) hh1
        _ = _ := Real.rpow_one _
    have hE0 : 1 ≤ Real.exp (h * L) := Real.one_le_exp (by positivity)
    rw [hE]
    set P := ∑ m ∈ Finset.range (K + 1), h ^ m * Wpoly m s with hP
    set S := ∑ r ∈ Finset.Icc 1 K, h ^ r * (L ^ r / (r.factorial : ℝ)) with hS
    have hTq : Tq K h s = ∑ n ∈ Finset.Icc 1 K, h ^ n *
        ∑ r ∈ Finset.Icc 1 n, L ^ r / (r.factorial : ℝ) * Wpoly (n - r) s := rfl
    have hid : X02.W h (s * h) * (Real.exp (h * L) - 1) - Tq K h s =
        (X02.W h (s * h) - P) * (Real.exp (h * L) - 1) + P * (Real.exp (h * L) - 1 - S)
          + (P * S - Tq K h s) := by ring
    rw [hid]
    -- the three pieces
    have b1 : |(X02.W h (s * h) - P) * (Real.exp (h * L) - 1)| ≤
        C * h ^ (K + 1) * (1 + s) ^ Nn * Real.exp (-s) := by
      rw [abs_mul, abs_of_nonneg (by linarith : 0 ≤ Real.exp (h * L) - 1)]
      have := hW h s hh hh0 hs.le
      calc _ ≤ (C * h ^ (K + 1) * (1 + s) ^ Nn) * Real.exp (-s) :=
            mul_le_mul this (by linarith) (by linarith) (by positivity)
        _ = _ := by ring
    have hς : |Real.exp (h * L) - 1 - S| ≤ 2 * h ^ (K + 1) * Real.exp (-s) := by
      have hr := exp_rem_pos (x := h * L) (by positivity) (K + 1)
      rw [range_split] at hr
      have e : Real.exp (h * L) - ((h * L) ^ 0 / ((0 : ℕ).factorial : ℝ) +
          ∑ r ∈ Finset.Icc 1 K, (h * L) ^ r / (r.factorial : ℝ)) = Real.exp (h * L) - 1 - S := by
        rw [hS]; simp only [pow_zero, Nat.factorial_zero, Nat.cast_one, div_one]
        have : ∑ r ∈ Finset.Icc 1 K, (h * L) ^ r / (r.factorial : ℝ) =
            ∑ r ∈ Finset.Icc 1 K, h ^ r * (L ^ r / (r.factorial : ℝ)) :=
          Finset.sum_congr rfl (fun r _ => by rw [mul_pow]; ring)
        rw [this]; ring
      rw [e] at hr
      have hLK : L ^ (K + 1) ≤ Real.exp (-s) :=
        (pow_le_of_le_one hL0 hL1 (by omega)).trans hLe
      calc _ ≤ (h * L) ^ (K + 1) * Real.exp (h * L) := hr
        _ = h ^ (K + 1) * L ^ (K + 1) * Real.exp (h * L) := by ring
        _ ≤ h ^ (K + 1) * Real.exp (-s) * 2 := by
            refine mul_le_mul (mul_le_mul_of_nonneg_left hLK (by positivity)) (by linarith)
              (by positivity) (by positivity)
        _ = _ := by ring
    have hPb := hD h s hh.le hh1 hs.le
    have b2 : |P * (Real.exp (h * L) - 1 - S)| ≤
        2 * D * h ^ (K + 1) * (1 + s) ^ Nn * Real.exp (-s) := by
      rw [abs_mul]
      have hKN : (1 + s) ^ K ≤ (1 + s) ^ Nn := pow_le_pow_right₀ hs1 (by omega)
      calc _ ≤ (D * (1 + s) ^ K) * (2 * h ^ (K + 1) * Real.exp (-s)) :=
            mul_le_mul hPb hς (abs_nonneg _) (by positivity)
        _ ≤ (D * (1 + s) ^ Nn) * (2 * h ^ (K + 1) * Real.exp (-s)) := by gcongr
        _ = _ := by ring
    have b3 : |P * S - Tq K h s| ≤ K * Ds * h ^ (K + 1) * (1 + s) ^ Nn * Real.exp (-s) := by
      have := prod_trunc K hh.le hh1 (fun m => Wpoly m s) (fun r => L ^ r / (r.factorial : ℝ))
      rw [hP, hS, hTq]
      refine this.trans ?_
      have hab : ∀ r ∈ Finset.Icc 1 K, ∀ m ∈ Finset.range (K + 1),
          |Wpoly m s| * |L ^ r / (r.factorial : ℝ)| ≤ Dm m * ((1 + s) ^ Nn * Real.exp (-s)) := by
        intro r hr m hm
        have hr1 : 1 ≤ r := (Finset.mem_Icc.mp hr).1
        have hmK : m ≤ K := by simp at hm; omega
        have hf : (1 : ℝ) ≤ (r.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos r
        have hb : |L ^ r / (r.factorial : ℝ)| ≤ Real.exp (-s) := by
          rw [abs_of_nonneg (by positivity)]
          exact (div_le_self (by positivity) hf).trans
            ((pow_le_of_le_one hL0 hL1 (by omega)).trans hLe)
        have ha : |Wpoly m s| ≤ Dm m * (1 + s) ^ Nn :=
          (hDm m s hs.le).trans (mul_le_mul_of_nonneg_left
            (pow_le_pow_right₀ hs1 (by omega)) (hDm0 m))
        calc _ ≤ (Dm m * (1 + s) ^ Nn) * Real.exp (-s) :=
              mul_le_mul ha hb (abs_nonneg _) (by have := hDm0 m; positivity)
          _ = _ := by ring
      have hsum : ∑ r ∈ Finset.Icc 1 K, ∑ m ∈ Finset.range (K + 1),
          |Wpoly m s| * |L ^ r / (r.factorial : ℝ)| ≤ K * Ds * ((1 + s) ^ Nn * Real.exp (-s)) := by
        calc _ ≤ ∑ r ∈ Finset.Icc 1 K, ∑ m ∈ Finset.range (K + 1),
              Dm m * ((1 + s) ^ Nn * Real.exp (-s)) :=
              Finset.sum_le_sum (fun r hr => Finset.sum_le_sum (fun m hm => hab r hr m hm))
          _ = K * Ds * ((1 + s) ^ Nn * Real.exp (-s)) := by
              rw [Finset.sum_const, Nat.card_Icc, nsmul_eq_mul, ← Finset.sum_mul, hDs]
              push_cast; ring
      calc h ^ (K + 1) * _ ≤ h ^ (K + 1) * (K * Ds * ((1 + s) ^ Nn * Real.exp (-s))) :=
            mul_le_mul_of_nonneg_left hsum (by positivity)
        _ = _ := by ring
    have t1 := abs_add_le ((X02.W h (s * h) - P) * (Real.exp (h * L) - 1) +
      P * (Real.exp (h * L) - 1 - S)) (P * S - Tq K h s)
    have t2 := abs_add_le ((X02.W h (s * h) - P) * (Real.exp (h * L) - 1))
      (P * (Real.exp (h * L) - 1 - S))
    have : h ^ (K + 1) * (M * ((1 + s) ^ Nn * Real.exp (-s))) =
        C * h ^ (K + 1) * (1 + s) ^ Nn * Real.exp (-s)
        + 2 * D * h ^ (K + 1) * (1 + s) ^ Nn * Real.exp (-s)
        + K * Ds * h ^ (K + 1) * (1 + s) ^ Nn * Real.exp (-s) := by rw [hM]; ring
    rw [this]
    linarith
  -- integrate
  have hQ : X02.Q h - ∑ n ∈ Finset.Icc 1 K, q n * h ^ n =
      ∫ s in Ioi (0 : ℝ), (X02.W h (s * h) * ((1 + Real.exp (-s)) ^ h - 1) - Tq K h s) := by
    rw [integral_sub (hI h hh (by linarith)) (intOn_Tq K h), integral_Tq]; rfl
  rw [hQ, ← Real.norm_eq_abs]
  have hB : IntegrableOn (fun s => h ^ (K + 1) * (M * ((1 + s) ^ Nn * Real.exp (-s)))) (Ioi 0) :=
    ((intOn_onepow Nn).const_mul M).const_mul _
  refine (norm_integral_le_of_norm_le hB ((ae_restrict_iff' measurableSet_Ioi).2
    (Filter.Eventually.of_forall (fun s hs => by rw [Real.norm_eq_abs]; exact hpt s hs)))).trans
    (le_of_eq ?_)
  rw [integral_const_mul, integral_const_mul]; ring



/-! ### `π_p/4 ∼ 1 + Σ N_k p^{-k-1}` to every order -/

lemma geom_id (K : ℕ) (h : ℝ) :
    (1 - h) * ∑ k ∈ Finset.range (K + 1), h ^ (k + 1) * ∑ n ∈ Finset.Ico 1 k, q n =
      h ^ 2 * (∑ n ∈ Finset.Icc 1 K, q n * h ^ n - h ^ K * ∑ n ∈ Finset.Icc 1 K, q n) := by
  induction K with
  | zero => simp
  | succ K ih =>
    rw [Finset.sum_range_succ, mul_add, ih, Finset.sum_Icc_succ_top (by omega),
      Finset.sum_Icc_succ_top (by omega)]
    have e : Finset.Ico 1 (K + 1) = Finset.Icc 1 K := by ext n; simp
    rw [e]
    ring

lemma two_rpow_taylor {h : ℝ} (hh : 0 ≤ h) (h1 : h ≤ 1) (K : ℕ) :
    |(2 : ℝ) ^ (-h) - 1 - ∑ k ∈ Finset.range (K + 1),
        (-Real.log 2) ^ (k + 1) / ((k + 1).factorial : ℝ) * h ^ (k + 1)| ≤ 2 * h ^ (K + 2) := by
  have hl0 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hl1 : Real.log 2 < 1 := by have := Real.log_two_lt_d9; linarith
  have h2 : (2 : ℝ) ^ (-h) = Real.exp (-(h * Real.log 2)) := by
    rw [Real.rpow_def_of_pos (by norm_num)]; congr 1; ring
  have hr := Complex.norm_exp_sub_sum_le_norm_mul_exp ((-(h * Real.log 2) : ℝ) : ℂ) (K + 2)
  have e : Complex.exp ((-(h * Real.log 2) : ℝ) : ℂ) - ∑ m ∈ Finset.range (K + 2),
      ((-(h * Real.log 2) : ℝ) : ℂ) ^ m / (m.factorial : ℂ) =
      ((Real.exp (-(h * Real.log 2)) - ∑ m ∈ Finset.range (K + 2),
        (-(h * Real.log 2)) ^ m / (m.factorial : ℝ) : ℝ) : ℂ) := by push_cast; rfl
  rw [e, Complex.norm_real, Real.norm_eq_abs, Complex.norm_real, Real.norm_eq_abs, abs_neg,
    abs_of_nonneg (mul_nonneg hh hl0.le : 0 ≤ h * Real.log 2)] at hr
  have e2 : ∑ m ∈ Finset.range (K + 2), (-(h * Real.log 2)) ^ m / (m.factorial : ℝ) =
      1 + ∑ k ∈ Finset.range (K + 1),
        (-Real.log 2) ^ (k + 1) / ((k + 1).factorial : ℝ) * h ^ (k + 1) := by
    rw [Finset.sum_range_succ']
    simp only [pow_zero, Nat.factorial_zero, Nat.cast_one, div_one]
    rw [add_comm]; congr 1
    refine Finset.sum_congr rfl (fun k _ => ?_)
    rw [neg_mul_eq_mul_neg, mul_pow]; ring
  rw [e2] at hr
  rw [h2, show Real.exp (-(h * Real.log 2)) - 1 - ∑ k ∈ Finset.range (K + 1),
      (-Real.log 2) ^ (k + 1) / ((k + 1).factorial : ℝ) * h ^ (k + 1) =
      Real.exp (-(h * Real.log 2)) - (1 + ∑ k ∈ Finset.range (K + 1),
      (-Real.log 2) ^ (k + 1) / ((k + 1).factorial : ℝ) * h ^ (k + 1)) by ring]
  refine hr.trans ?_
  have ha : h * Real.log 2 ≤ h := by nlinarith
  have hb : Real.exp (h * Real.log 2) ≤ 2 := by
    calc Real.exp (h * Real.log 2) ≤ Real.exp (Real.log 2) :=
          Real.exp_le_exp.mpr (by nlinarith)
      _ = 2 := Real.exp_log (by norm_num)
  calc (h * Real.log 2) ^ (K + 2) * Real.exp (h * Real.log 2) ≤ h ^ (K + 2) * 2 :=
        mul_le_mul (pow_le_pow_left₀ (by positivity) ha _) hb (by positivity) (by positivity)
    _ = _ := by ring

theorem piP_all_orders
    (hI : ∀ h : ℝ, 0 < h → h < 1 →
      IntegrableOn (fun s => X02.W h (s * h) * ((1 + Real.exp (-s)) ^ h - 1)) (Ioi 0))
    (hE : ∀ h : ℝ, 0 < h → h < 1 → Ehat h = h ^ 2 / (1 - h) * X02.Q h)
    (hsplit : ∀ p : ℝ, 1 < p → X02.piP p = 4 * (2 : ℝ) ^ (-(1 / p)) + X02.E p) (K : ℕ) :
    ∃ C p₀ : ℝ, 1 < p₀ ∧ ∀ p : ℝ, p₀ ≤ p →
      |X02.piP p / 4 - 1 - ∑ k ∈ Finset.range (K + 1), N k * (1 / p) ^ (k + 1)|
        ≤ C * (1 / p) ^ (K + 2) := by
  obtain ⟨C, h0, h0pos, hQ⟩ := Q_expand (K + 1) hI
  set S := ∑ n ∈ Finset.Icc 1 K, q n
  refine ⟨2 + 2 * (|C| + |q (K + 1)| + |S|), max 2 (1 / h0), lt_of_lt_of_le (by norm_num)
    (le_max_left _ _), fun p hp => ?_⟩
  have hp2 : 2 ≤ p := (le_max_left _ _).trans hp
  have hp0 : 0 < p := by linarith
  set h := 1 / p with hhdef
  have hh : 0 < h := by positivity
  have hh2 : h ≤ 1 / 2 := by
    rw [hhdef, div_le_div_iff₀ hp0 (by norm_num)]; linarith
  have hh0 : h ≤ h0 := by
    have : 1 / h0 ≤ p := (le_max_right _ _).trans hp
    rw [hhdef, div_le_iff₀ hp0]; rw [div_le_iff₀ h0pos] at this; linarith
  have hpi : X02.piP p / 4 = (2 : ℝ) ^ (-h) + h ^ 2 / (1 - h) * X02.Q h := by
    rw [hsplit p (by linarith), ← hE h hh (by linarith)]
    unfold Ehat; rw [hhdef, one_div_one_div]; ring
  have hN : ∑ k ∈ Finset.range (K + 1), N k * h ^ (k + 1) =
      ∑ k ∈ Finset.range (K + 1), (-Real.log 2) ^ (k + 1) / ((k + 1).factorial : ℝ) * h ^ (k + 1)
      + ∑ k ∈ Finset.range (K + 1), h ^ (k + 1) * ∑ n ∈ Finset.Ico 1 k, q n := by
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl (fun k _ => ?_)
    unfold N; ring
  have h1h : 0 < 1 - h := by linarith
  have hG : ∑ k ∈ Finset.range (K + 1), h ^ (k + 1) * ∑ n ∈ Finset.Ico 1 k, q n =
      h ^ 2 / (1 - h) * (∑ n ∈ Finset.Icc 1 K, q n * h ^ n - h ^ K * S) := by
    rw [div_mul_eq_mul_div, eq_div_iff h1h.ne', mul_comm]
    exact geom_id K h
  rw [hpi, hN, hG]
  set Tt := ∑ n ∈ Finset.Icc 1 K, q n * h ^ n with hT
  set A := ∑ k ∈ Finset.range (K + 1),
    (-Real.log 2) ^ (k + 1) / ((k + 1).factorial : ℝ) * h ^ (k + 1) with hA
  have hQK := hQ h hh hh0
  rw [Finset.sum_Icc_succ_top (by omega), ← hT] at hQK
  have hh1 : h ≤ 1 := by linarith
  have hhK : h ^ (K + 1) ≤ h ^ K := pow_le_pow_of_le_one hh.le hh1 (by omega)
  have hhK2 : h ^ (K + 2) ≤ h ^ K := pow_le_pow_of_le_one hh.le hh1 (by omega)
  have hB : |X02.Q h - Tt + h ^ K * S| ≤ (|C| + |q (K + 1)| + |S|) * h ^ K := by
    have e : X02.Q h - Tt + h ^ K * S = (X02.Q h - (Tt + q (K + 1) * h ^ (K + 1)))
        + q (K + 1) * h ^ (K + 1) + h ^ K * S := by ring
    rw [e]
    have t1 := abs_add_le (X02.Q h - (Tt + q (K + 1) * h ^ (K + 1)) + q (K + 1) * h ^ (K + 1))
      (h ^ K * S)
    have t2 := abs_add_le (X02.Q h - (Tt + q (K + 1) * h ^ (K + 1))) (q (K + 1) * h ^ (K + 1))
    have a1 : C * h ^ (K + 1 + 1) ≤ |C| * h ^ K :=
      mul_le_mul (le_abs_self C) hhK2 (by positivity) (abs_nonneg _)
    have a2 : |q (K + 1) * h ^ (K + 1)| ≤ |q (K + 1)| * h ^ K := by
      rw [abs_mul, abs_of_nonneg (pow_nonneg hh.le (K + 1))]
      exact mul_le_mul_of_nonneg_left hhK (abs_nonneg _)
    have a3 : |h ^ K * S| = |S| * h ^ K := by
      rw [abs_mul, abs_of_nonneg (pow_nonneg hh.le K), mul_comm]
    linarith
  have hc : h ^ 2 / (1 - h) ≤ 2 * h ^ 2 := by
    rw [div_le_iff₀ h1h]; nlinarith
  have e : (2 : ℝ) ^ (-h) + h ^ 2 / (1 - h) * X02.Q h - 1 - (A + h ^ 2 / (1 - h) * (Tt - h ^ K * S)) =
      ((2 : ℝ) ^ (-h) - 1 - A) + h ^ 2 / (1 - h) * (X02.Q h - Tt + h ^ K * S) := by ring
  rw [e]
  have t := abs_add_le ((2 : ℝ) ^ (-h) - 1 - A) (h ^ 2 / (1 - h) * (X02.Q h - Tt + h ^ K * S))
  have r1 := two_rpow_taylor hh.le hh1 K
  have r2 : |h ^ 2 / (1 - h) * (X02.Q h - Tt + h ^ K * S)| ≤
      2 * (|C| + |q (K + 1)| + |S|) * h ^ (K + 2) := by
    rw [abs_mul, abs_of_nonneg (by positivity)]
    calc h ^ 2 / (1 - h) * |X02.Q h - Tt + h ^ K * S| ≤
        (2 * h ^ 2) * ((|C| + |q (K + 1)| + |S|) * h ^ K) :=
          mul_le_mul hc hB (abs_nonneg _) (by positivity)
      _ = _ := by rw [show K + 2 = 2 + K by ring, pow_add]; ring
  calc _ ≤ _ := t
    _ ≤ 2 * h ^ (K + 2) + 2 * (|C| + |q (K + 1)| + |S|) * h ^ (K + 2) := add_le_add r1 r2
    _ = _ := by ring

end Q02

/-- Lemma `lem:pip:expansion`, the core: `Q(h) = Σ_{n ≤ K} q_n h^n + O(h^{K+1})`. -/
theorem Q02.proof_7 : Q02.stmt_7 := fun K _ =>
  Q02.Q_expand K (fun h hh h1 => (Q02.proof_6 h hh h1).1)

/-- thm:pip: `π_p/4 ∼ 1 + Σ N_k p^{-k-1}` to every order. -/
theorem Q02.proof_8 : Q02.stmt_8 :=
  Q02.piP_all_orders (fun h hh h1 => (Q02.proof_6 h hh h1).1) (fun h hh h1 => (Q02.proof_6 h hh h1).2)
    (fun p hp => (X02.proof_1 p hp).2.2)

namespace Q02

/-! ### Lemma `lem:pip:hankel`, case `h = 0`: Cauchy's estimate on `|y| = 2π` -/

/-- `g(z) = e^{2z}/(1+e^z) + 2z/(z²+π²)`, with removable singularities at `±iπ`. -/
noncomputable def gR (z : ℂ) : ℂ :=
  Complex.exp (2 * z) / (1 + Complex.exp z) + 2 * z / (z ^ 2 + (Real.pi : ℂ) ^ 2)

lemma exp_ipi : Complex.exp (Complex.I * Real.pi) = -1 := by
  rw [mul_comm]; exact Complex.exp_pi_mul_I

lemma exp_neg_ipi : Complex.exp (-(Complex.I * Real.pi)) = -1 := by
  rw [Complex.exp_neg, exp_ipi]; norm_num

lemma sq_ipi (c : ℂ) (hc : c = Complex.I * Real.pi ∨ c = -(Complex.I * Real.pi)) :
    c ^ 2 = -(Real.pi : ℂ) ^ 2 := by
  rcases hc with rfl | rfl <;> ring_nf <;> simp [Complex.I_sq]

lemma one_add_exp_ne {z : ℂ} (hz : ‖z‖ < 3 * Real.pi) (h1 : z ≠ Complex.I * Real.pi)
    (h2 : z ≠ -(Complex.I * Real.pi)) : 1 + Complex.exp z ≠ 0 := by
  intro h
  have he : Complex.exp (z - Complex.I * Real.pi) = 1 := by
    rw [Complex.exp_sub, exp_ipi]
    have : Complex.exp z = -1 := by linear_combination h
    rw [this]; norm_num
  obtain ⟨n, hn⟩ := Complex.exp_eq_one_iff.mp he
  have hz' : z = ((2 * n + 1 : ℤ) : ℂ) * (Real.pi * Complex.I) := by
    rw [sub_eq_iff_eq_add] at hn; rw [hn]; push_cast; ring
  have hnorm : ‖z‖ = |((2 * n + 1 : ℤ) : ℝ)| * Real.pi := by
    rw [hz', norm_mul, norm_mul, Complex.norm_I, mul_one, Complex.norm_real,
      Real.norm_eq_abs, abs_of_pos Real.pi_pos]
    rw [show (((2 * n + 1 : ℤ)) : ℂ) = (((2 * n + 1 : ℤ) : ℝ) : ℂ) by push_cast; rfl,
      Complex.norm_real, Real.norm_eq_abs]
  rw [hnorm] at hz
  have hlt : |((2 * n + 1 : ℤ) : ℝ)| < 3 := by
    have := Real.pi_pos
    nlinarith [abs_nonneg ((2 * n + 1 : ℤ) : ℝ)]
  have hlt' : |2 * n + 1| < 3 := by exact_mod_cast hlt
  rw [abs_lt] at hlt'
  have : n = 0 ∨ n = -1 := by omega
  rcases this with rfl | rfl
  · apply h1; rw [hz']; push_cast; ring
  · apply h2; rw [hz']; push_cast; ring

lemma sq_add_ne {z : ℂ} (h1 : z ≠ Complex.I * Real.pi) (h2 : z ≠ -(Complex.I * Real.pi)) :
    z ^ 2 + (Real.pi : ℂ) ^ 2 ≠ 0 := by
  have e : z ^ 2 + (Real.pi : ℂ) ^ 2 = (z - Complex.I * Real.pi) * (z + Complex.I * Real.pi) := by
    ring_nf; rw [Complex.I_sq]; ring
  rw [e]
  exact mul_ne_zero (sub_ne_zero.mpr h1) (fun h => h2 (by linear_combination h))

/-- the punctured disc -/
noncomputable def Dp : Set ℂ := Metric.ball 0 (3 * Real.pi) \ {Complex.I * Real.pi, -(Complex.I * Real.pi)}

lemma gR_diffOn : DifferentiableOn ℂ gR Dp := by
  intro z hz
  simp only [Dp, Set.mem_sdiff, Metric.mem_ball, dist_zero_right, mem_insert_iff,
    mem_singleton_iff, not_or] at hz
  obtain ⟨hz, h1, h2⟩ := hz
  apply DifferentiableAt.differentiableWithinAt
  unfold gR
  exact ((Complex.differentiableAt_exp.comp z (by fun_prop)).div
    (by fun_prop) (one_add_exp_ne hz h1 h2)).add
    ((by fun_prop : DifferentiableAt ℂ (fun z : ℂ => 2 * z) z).div (by fun_prop) (sq_add_ne h1 h2))

/-- the pole of `e^{2z}/(1+e^z)` at `c = ±iπ` is cancelled -/
lemma gR_tendsto (c : ℂ) (hc : c = Complex.I * Real.pi ∨ c = -(Complex.I * Real.pi)) :
    Tendsto (fun z => (z - c) * gR z) (𝓝[≠] c) (𝓝 0) := by
  have hec : Complex.exp c = -1 := by rcases hc with rfl | rfl; exacts [exp_ipi, exp_neg_ipi]
  have hc2 := sq_ipi c hc
  have hc0 : c ≠ 0 := by
    rcases hc with rfl | rfl <;> simp [Complex.I_ne_zero, Real.pi_ne_zero]
  -- the slope of `1 + e^z` at `c`
  have hd : HasDerivAt (fun z => 1 + Complex.exp z) (Complex.exp c) c :=
    (Complex.hasDerivAt_exp c).const_add 1
  have hs := hasDerivAt_iff_tendsto_slope.mp hd
  rw [hec] at hs
  have hs' : Tendsto (fun z => (z - c) / (1 + Complex.exp z)) (𝓝[≠] c) (𝓝 (-1)⁻¹) := by
    refine (hs.inv₀ (by norm_num)).congr' ?_
    filter_upwards [self_mem_nhdsWithin] with z hz
    simp only [slope_def_field, hec]
    rw [show (1 : ℂ) + -1 = 0 by ring, sub_zero, inv_div]
  have h2 : Tendsto (fun z => Complex.exp (2 * z)) (𝓝[≠] c) (𝓝 (Complex.exp (2 * c))) :=
    ((Complex.continuous_exp.comp (continuous_const.mul continuous_id)).tendsto c).mono_left
      nhdsWithin_le_nhds
  have h3 : Tendsto (fun z : ℂ => 2 * z / (z + c)) (𝓝[≠] c) (𝓝 (2 * c / (c + c))) :=
    (((continuous_const.mul continuous_id).tendsto c).div
      ((continuous_id.add continuous_const).tendsto c) (by
        intro h; simp only [Pi.add_apply, id] at h; apply hc0; linear_combination h / 2)).mono_left
        nhdsWithin_le_nhds
  have hlim := (h2.mul hs').add h3
  have hv : Complex.exp (2 * c) * (-1)⁻¹ + 2 * c / (c + c) = 0 := by
    rw [show 2 * c = c + c by ring, Complex.exp_add, hec, div_self (by
      intro h; apply hc0; linear_combination h / 2)]
    norm_num
  rw [hv] at hlim
  refine hlim.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with z hz
  have hzc : z - c ≠ 0 := sub_ne_zero.mpr hz
  have e : z ^ 2 + (Real.pi : ℂ) ^ 2 = (z - c) * (z + c) := by
    ring_nf; rw [hc2]; ring
  unfold gR
  rw [e]
  by_cases h1 : 1 + Complex.exp z = 0
  · -- at a zero of `1 + e^z` both sides use Lean's `x / 0 = 0`
    rw [h1]; simp
    by_cases h4 : z + c = 0
    · simp [h4]
    · field_simp
  · by_cases h4 : z + c = 0
    · simp [h4]; ring
    · field_simp


lemma littleO_of_tendsto (f : ℂ → ℂ) (c : ℂ)
    (hT : Tendsto (fun z => (z - c) * f z) (𝓝[≠] c) (𝓝 0)) :
    (fun z => f z - f c) =o[𝓝[≠] c] fun z => (z - c)⁻¹ := by
  rw [Asymptotics.isLittleO_iff_tendsto']
  · have h2 : Tendsto (fun z => (z - c) * f c) (𝓝[≠] c) (𝓝 0) := by
      have : Tendsto (fun z : ℂ => (z - c) * f c) (𝓝 c) (𝓝 ((c - c) * f c)) :=
        ((continuous_id.sub continuous_const).mul continuous_const).tendsto c
      rw [sub_self, zero_mul] at this
      exact this.mono_left nhdsWithin_le_nhds
    have := hT.sub h2
    rw [sub_zero] at this
    refine this.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with z hz
    have hzc : z - c ≠ 0 := sub_ne_zero.mpr hz
    field_simp
  · filter_upwards [self_mem_nhdsWithin] with z hz
    intro h
    exact absurd h (inv_ne_zero (sub_ne_zero.mpr hz))

/-- the first removable singularity, at `iπ` -/
noncomputable def g1 : ℂ → ℂ :=
  Function.update gR (Complex.I * Real.pi) (limUnder (𝓝[≠] (Complex.I * Real.pi)) gR)

/-- the second removable singularity, at `-iπ` -/
noncomputable def g2 : ℂ → ℂ :=
  Function.update g1 (-(Complex.I * Real.pi)) (limUnder (𝓝[≠] (-(Complex.I * Real.pi))) g1)

lemma ipi_ne : Complex.I * Real.pi ≠ -(Complex.I * Real.pi) := by
  intro h
  have : Complex.I * Real.pi = 0 := by linear_combination h / 2
  simp [Complex.I_ne_zero, Real.pi_ne_zero] at this

lemma norm_ipi : ‖Complex.I * (Real.pi : ℂ)‖ = Real.pi := by
  rw [norm_mul, Complex.norm_I, one_mul, Complex.norm_real, Real.norm_eq_abs,
    abs_of_pos Real.pi_pos]

lemma g1_diffOn : DifferentiableOn ℂ g1 (Metric.ball 0 (3 * Real.pi) \ {-(Complex.I * Real.pi)}) := by
  refine Complex.differentiableOn_update_limUnder_of_isLittleO ?_ ?_ ?_
  · refine (Metric.isOpen_ball.sdiff isClosed_singleton).mem_nhds ⟨?_, ?_⟩
    · rw [Metric.mem_ball, dist_zero_right, norm_ipi]; linarith [Real.pi_pos]
    · exact ipi_ne
  · refine gR_diffOn.mono (fun z hz => ?_)
    simp only [Dp, Set.mem_sdiff, Set.mem_insert_iff, Set.mem_singleton_iff, not_or] at hz ⊢
    exact ⟨hz.1.1, hz.2, hz.1.2⟩
  · exact littleO_of_tendsto gR _ (gR_tendsto _ (Or.inl rfl))

lemma g2_diffOn : DifferentiableOn ℂ g2 (Metric.ball 0 (3 * Real.pi)) := by
  refine Complex.differentiableOn_update_limUnder_of_isLittleO ?_ g1_diffOn ?_
  · refine Metric.isOpen_ball.mem_nhds ?_
    rw [Metric.mem_ball, dist_zero_right, norm_neg, norm_ipi]; linarith [Real.pi_pos]
  · refine littleO_of_tendsto g1 _ ?_
    have hT := gR_tendsto (-(Complex.I * Real.pi)) (Or.inr rfl)
    refine hT.congr' ?_
    have hne : ∀ᶠ z in 𝓝[≠] (-(Complex.I * Real.pi)), z ≠ Complex.I * Real.pi :=
      (isOpen_ne.mem_nhds ipi_ne.symm |> nhdsWithin_le_nhds)
    filter_upwards [hne] with z hz
    rw [g1, Function.update_of_ne hz]

lemma g2_real (y : ℝ) : g2 (y : ℂ) = gR (y : ℂ) := by
  have h1 : (y : ℂ) ≠ Complex.I * Real.pi := by
    intro h; have := congrArg Complex.im h; simp at this; exact Real.pi_ne_zero this.symm
  have h2 : (y : ℂ) ≠ -(Complex.I * Real.pi) := by
    intro h; have := congrArg Complex.im h; simp at this
  rw [g2, Function.update_of_ne h2, g1, Function.update_of_ne h1]

lemma gR_real (y : ℝ) :
    gR (y : ℂ) = ((Acal y 0 + 2 * y / (y ^ 2 + Real.pi ^ 2) : ℝ) : ℂ) := by
  unfold gR Acal
  rw [sub_zero, Real.rpow_neg_one]
  push_cast
  ring

/-- the real Taylor coefficients are the real parts of the complex ones -/
lemma iteratedDeriv_real (k : ℕ) : ∀ y : ℝ, |y| < 3 * Real.pi →
    iteratedDeriv k (fun y : ℝ => Acal y 0 + 2 * y / (y ^ 2 + Real.pi ^ 2)) y =
      (iteratedDeriv k g2 (y : ℂ)).re := by
  have han : AnalyticOnNhd ℂ g2 (Metric.ball 0 (3 * Real.pi)) :=
    g2_diffOn.analyticOnNhd Metric.isOpen_ball
  induction k with
  | zero =>
    intro y _
    simp only [iteratedDeriv_zero]
    rw [g2_real, gR_real, Complex.ofReal_re]
  | succ k ih =>
    intro y hy
    have hopen : IsOpen {t : ℝ | |t| < 3 * Real.pi} :=
      isOpen_lt continuous_abs continuous_const
    have heq : iteratedDeriv k (fun y : ℝ => Acal y 0 + 2 * y / (y ^ 2 + Real.pi ^ 2)) =ᶠ[𝓝 y]
        fun t : ℝ => (iteratedDeriv k g2 (t : ℂ)).re := by
      filter_upwards [hopen.mem_nhds hy] with t ht using ih t ht
    rw [iteratedDeriv_succ, heq.deriv_eq]
    have hmem : (y : ℂ) ∈ Metric.ball (0 : ℂ) (3 * Real.pi) := by
      rw [Metric.mem_ball, dist_zero_right, Complex.norm_real, Real.norm_eq_abs]; exact hy
    have hdk : DifferentiableAt ℂ (iteratedDeriv k g2) (y : ℂ) := by
      rw [iteratedDeriv_eq_iterate]
      exact ((han.iterated_deriv k) (y : ℂ) hmem).differentiableAt
    have := hdk.hasDerivAt.real_of_complex
    rw [this.deriv, iteratedDeriv_succ]

theorem hankel0 : ∃ C : ℝ, ∀ m : ℕ,
    |iteratedDeriv m (fun y : ℝ => Acal y 0 + 2 * y / (y ^ 2 + Real.pi ^ 2)) 0 /
        (m.factorial : ℝ)| ≤ C / (2 * Real.pi) ^ m := by
  have hpi := Real.pi_pos
  have hsub : Metric.closedBall (0 : ℂ) (2 * Real.pi) ⊆ Metric.ball 0 (3 * Real.pi) :=
    Metric.closedBall_subset_ball (by linarith)
  have hdc : DiffContOnCl ℂ g2 (Metric.ball 0 (2 * Real.pi)) := g2_diffOn.diffContOnCl_ball hsub
  have hcont : ContinuousOn g2 (Metric.sphere (0 : ℂ) (2 * Real.pi)) :=
    (g2_diffOn.continuousOn).mono (Metric.sphere_subset_closedBall.trans hsub)
  obtain ⟨C, hC⟩ := (isCompact_sphere (0 : ℂ) (2 * Real.pi)).exists_bound_of_continuousOn hcont
  refine ⟨C, fun m => ?_⟩
  have hcau := Complex.norm_iteratedDeriv_le_of_forall_mem_sphere_norm_le m (by positivity) hdc hC
  rw [iteratedDeriv_real m 0 (by simp; positivity)]
  have hf : (0 : ℝ) < (m.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos m
  rw [abs_div, abs_of_pos hf, div_le_iff₀ hf]
  have h0 : ((0 : ℝ) : ℂ) = 0 := by simp
  rw [h0]
  calc |(iteratedDeriv m g2 0).re| ≤ ‖iteratedDeriv m g2 0‖ := Complex.abs_re_le_norm _
    _ ≤ m.factorial * C / (2 * Real.pi) ^ m := hcau
    _ = _ := by ring

end Q02

/-- Lemma `lem:pip:hankel`, case `h = 0`: Cauchy's estimate on `|y| = 2π` for the function
`e^{2y}/(1+e^y) + 2y/(y²+π²)`, holomorphic on `|y| < 3π` (its poles at `±iπ` cancel). -/
theorem Q02.proof_14 : Q02.stmt_14 := Q02.hankel0
