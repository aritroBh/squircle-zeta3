import Mathlib
import Zeta3.Wave2.D04_MalphaForms
import Zeta3.Wave4.X05_DiamondA4A5
import Zeta3.Wave7.V05_FlatAngle
/-!
# Q04 -- expanding the corner coefficients `M_k(α)` about the flat angle, every order `k`

Target id: Q04.
Paper labels: `thm:flat` (paper/ieee/main.tex, lines 1501-1511) and the first rung of
`thm:zladder` (main.tex, lines 1591-1593).
Supplement: paper/ieee/proofs/sec_g6_flat.tex (`g6:thm:flat` 13-20; "What the statement means"
22-23; Lemma `g6:lem:series` 36-48; Lemma `g6:lem:En` 82-88 with proof 90-97; Section
`g6:sec:flat:beta` 99-104 with `g6:eq:iso`; Lemma `g6:lem:gamma` 110-112 with proof 114-120;
Definition `g6:def:formal` 122-127; Lemma `g6:prop:formula` 129-131; proof of the theorem
180-186, including the coefficient bound) and paper/ieee/proofs/sec_g6_ladder.tex (meaning of
`[z]M` 8-12, `g6:thm:zladder` 15-17, precise form `g6:thm:ladderprecise` 29-35, the rung for
`ζ(3)` line 109).

Wave 7 file `Zeta3.Wave7.V05_FlatAngle` proved `thm:flat` for `k ≤ 3` (with a formal coefficient
in two variables `x2, x3`), and its `V05.stmt_4`, `V05.stmt_5` (Lemma `g6:lem:series`(c)) hold
for every `k`. They are imported, not restated: `M_k(α)` is `V05.Mk k α`
`= iteratedDeriv k (D04.Mint · α) 0 / k!` (the paper's Taylor coefficient `M^{(k)}(0;α)/k!` of
the corner integral `M(u;α) = D04.Mint u α`), `c = D04.c α = cos(α/2)`, the kernel
`E_n(u,y)` is `V05.E n u y`, `H^{(r)}_N` is `V05.Hr r N` (real) and `V05.HrQ r N` (rational),
`ρ_j(n)` is `V05.rhoQ j n`, `B(p,q) = (p-1)!(q-1)!/(p+q-1)!` is `V05.Bet p q`, and
`Φ(c) = 1/2 - (1-c²)/(2c) artanh c` is `V05.Phi`.

Quoted verbatim (main.tex 1501-1511, `thm:flat`):

  \begin{theorem}[Expanding about the flat angle$^\dagger$]\label{thm:flat}
  For every $k\ge0$ and $0<\alpha<\pi$, $M_k(\alpha)=\sum_{n\ge1}\mu_{k,n}c^{2n}$, where each
  $\mu_{k,n}$ is a polynomial with rational coefficients in $\zeta(2),\dots,\zeta(k)$ of weight at
  most $k$. In particular, with $r_n,r_n',q_n'\in\mathbb Q$,
  \[
  \mu_{2,n}=r_n+\frac{\zeta(2)}{4n^2-1},\qquad
  \mu_{3,n}=r_n'+q_n'\,\zeta(2)-\frac{\zeta(3)}{4n^2-1}.
  \]
  So, with $\Phi(c)=\sum_{n\ge1}\frac{c^{2n}}{4n^2-1}=\frac12-\frac{1-c^2}{2c}\operatorname{artanh}c$,
  the functions $M_2-\zeta(2)\Phi$ and $M_3+\zeta(3)\Phi-\zeta(2)\Psi$ are power series in $c^2$
  with rational coefficients, for some $\Psi\in\mathbb Q[[c^2]]$.
  \end{theorem}

Quoted verbatim (main.tex 1591-1593, `thm:zladder`, first sentence formalized):

  \begin{theorem}[The zeta ladder$^\dagger$]\label{thm:zladder}
  For every $k\ge0$ and $0<\alpha<\pi$, $[\zeta(3)]M_{k+1}(\alpha)=-[\zeta(2)]M_k(\alpha)$. More
  generally, ...

Quoted verbatim (main.tex 1589): "As in Theorem~\ref{thm:flat}, we treat the zeta values as
symbols, and we write $[z]M$ for the coefficient of the product $z$ of zeta values in $M$."

Quoted verbatim (supplement, flat, line 23): "``A polynomial in $\zeta(2),\dots,\zeta(k)$ of weight
at most $k$'' means an expression $P(\zeta(2),\dots,\zeta(k))$ with $P\in\QQ[\gvs2,\dots,\gvs k]$
in which every monomial has weight at most $k$ ($\gvs w$ has weight $w$)."

Quoted verbatim (supplement, ladder, line 10):
  [z]M_k(\alpha):=\sum_{n\ge1}\bigl([\gvs z]P_{k,n}\bigr)c^{2n}\ ,\qquad [z]M_k:=0\ \text{ for }k<0 .

Quoted verbatim (supplement, Lemma `g6:lem:En`, line 84): "For $n\ge1$ and $k\ge0$ there is a
polynomial $E_{n,k}(y,a,b)=\sum e^{(k)}_{j,i,m}(n)\,y^ja^ib^m$, with coefficients
$e^{(k)}_{j,i,m}\in\QQ[n]$, such that $\partial_u^kE_n(0,y)/k!=E_{n,k}(y,\log y,\log(1-y))$. Here
$j\le k$ and $i+m=k$ in every monomial." (proof, line 91: "$\log f=\sum_{m\ge1}(-1)^{m-1}y^m
(\eu^x-1)^m/m$ ... $\log E_n=uA+(1-2n)\log f$ is a formal series in $u$ ... Exponentiating gives
the same for $E_n$".)

Quoted verbatim (supplement, Section `g6:sec:flat:beta`, line 102):
  \Delta_\star=\sum_{r\ge1}\frac{(-1)^{r-1}}r\Bigl[x^rh^{(r)}_1+x'^rh^{(r)}_2-(x+x')^rh^{(r)}_3\Bigr]
  +\sum_{w\ge2}\frac{(-1)^{w+1}}w\,\gvs w\Bigl[(x+x')^w-x^w-x'^w\Bigr],
"and define $\mathfrak D_{i,m}$ by $\exp(\Delta_\star)=\sum_{i,m\ge0}\mathfrak D_{i,m}\,x^ix'^m/
(i!\,m!)$." ... "let $D^\sharp_{i,m}(p,q)\in\QQ[\gvs2,\gvs3,\dots]$ be $\mathfrak D_{i,m}$ with
$h^{(r)}_1,h^{(r)}_2,h^{(r)}_3$ replaced by the harmonic numbers $H^{(r)}_{p-1}$, $H^{(r)}_{q-1}$,
$H^{(r)}_{p+q-1}$".

Quoted verbatim (supplement, Lemma `g6:lem:gamma`, line 111): "For integers $p,q\ge1$ and
$i,m\ge0$, $D_{i,m}(p,q):=\partial_p^i\partial_q^m\gBeta(p,q)/\gBeta(p,q)$ is the value of
$D^\sharp_{i,m}(p,q)$ at $\gvs w=\zeta(w)$. In particular $D^\sharp_{i,m}$ has rational
coefficients, only $\gvs2,\dots,\gvs{i+m}$ occur in it, and its monomials have weight at most
$i+m$." (proof, line 116: "$\log\Gamma(1+x)=-\gamma_{\mathrm E}x+\sum_{w\ge2}(-1)^w\zeta(w)x^w/w$
for $|x|<1$ ... it follows from $\psi(1)=-\gamma_{\mathrm E}$ and $\psi^{(m)}(1)=(-1)^{m+1}m!\,
\zeta(m+1)$".)

Quoted verbatim (supplement, Definition `g6:def:formal` and Lemma `g6:prop:formula`, 123-130):
  P_{k,n}=-\frac1{n(2n-1)}\sum_{j=0}^k\ \sum_{i+m=k}e^{(k)}_{j,i,m}(n)\,\rho_j(n)\,
  D^\sharp_{i,m}(n+j,n)\ \in\QQ[\gvs2,\dots,\gvs k].
  "For $0<\alpha<\pi$ and $k\ge0$, $M_k(\alpha)=\sum_{n\ge1}\mu_{k,n}c^{2n}$ with
  $\mu_{k,n}=P_{k,n}(\zeta(2),\dots,\zeta(k))$, absolutely convergent. Every monomial of $P_{k,n}$
  has weight at most $k$."

Quoted verbatim (supplement, proof of `g6:thm:flat`, line 185): "Hence every coefficient of
$P_{k,n}$ (in particular $r_n,r'_n,q'_n$) is at most $C_k\,n^{k}(1+\log(2n+k))^{k}$ in absolute
value".

## Formalization choices

* Real zeta values: `ζ(w)` is `X05.z w = (riemannZeta w).re` (imported from
  `Zeta3.Wave4.X05_DiamondA4A5`; `riemannZeta w` is real for real `w > 1`, and
  `zeta_nat_eq_tsum_of_gt_one` identifies it with `Σ_{n≥1} n^{-w}`; `X05.z 2 = π²/6`). Note that
  `V05.z2 = π²/6 = X05.z 2` and `V05.z3 = X05.z 3`.
* Formal zeta values: `MvPolynomial ℕ ℚ`, the variable `X w` standing for `𝔷_w` (only `w ≥ 2`
  occurs, which is part of the weight predicate). The specialisation `𝔷_w ↦ ζ(w)` is
  `MvPolynomial.aeval X05.z`. "Polynomial with rational coefficients in `ζ(2),…,ζ(k)` of weight at
  most `k`" is the predicate `WeightLE k P`: every monomial `∏ 𝔷_w^{e_w}` of the support has all
  `w ≥ 2` and `Σ e_w w ≤ k` (so only `𝔷_2,…,𝔷_k` occur). "The coefficient `[𝔷_z]P`" is
  `MvPolynomial` coefficient extraction `P.coeff s`; `[𝔷_2]` is `coeff (single 2 1)`, `[𝔷_3]` is `coeff (single 3 1)`.
* `E_{n,k}` (Lemma `g6:lem:En`) is constructed exactly as in its proof, as a truncated formal
  computation in `ℚ[y,a,b][u]` (variables `y = X 0`, `a = X 1`, `b = X 2` of
  `MvPolynomial (Fin 3) ℚ`): `x = uℓ`, `ℓ = b - a`, `e^{x} - 1 ↦ Σ_{r=1}^K x^r/r!`,
  `log(1+w) ↦ Σ_{m=1}^K (-1)^{m-1} w^m/m`, `L = uA + (1-2n) log f`, `exp L ↦ Σ_{r=0}^K L^r/r!`,
  and `E_{n,k}` is the `u^k` coefficient with `K = k`. Since `L` has no constant term in `u`, the
  truncation does not change any coefficient of `u^k`, `k ≤ K`. `e^{(k)}_{j,i,m}(n)` is the
  coefficient of `y^j a^i b^m` (`eF k n j i m`).
* `D^♯_{i,m}(p,q)` is constructed the same way: `Δ_⋆` truncated at degree `K = i+m` in
  `x = X 0`, `x' = X 1` of `MvPolynomial (Fin 2) (MvPolynomial ℕ ℚ)`, with `h^{(r)}` specialised
  to the rational harmonic numbers, `exp ↦ Σ_{r=0}^K Δ^r/r!` (exact in degrees `≤ K`, as `Δ_⋆`
  has no constant term), and `D^♯_{i,m} = i! m! [x^i x'^m]`.
* `P_{k,n}` (`Pf k n`) is Definition `g6:def:formal` verbatim from these pieces. For `k ≤ 3` it
  agrees with `V05.P` (not stated; a consistency fact, not a claim of the paper).
* `D_{i,m}(p,q) = ∂_p^i∂_q^m B(p,q)/B(p,q)` is used, as in the supplement (line 115), through
  `∂_p^i∂_q^m B(p,q) = ∫_0^1 y^{p-1}(1-y)^{q-1} log^i y log^m(1-y) dy`; the statement of
  `g6:lem:gamma` is given in this integral form (with integrability stated, so the interval
  integral is not a junk value).
* The polygamma function `ψ^{(j)} = d^{j+1}/dx^{j+1} log Γ` is written with `iteratedDeriv` of
  `x ↦ Real.log (Real.Gamma x)`, which is smooth on `(0,∞)` (`Γ > 0` there), so `iteratedDeriv`
  at an integer `m ≥ 1` is the true derivative. Euler's constant is
  `Real.eulerMascheroniConstant`.
* Sums over `n ≥ 1` are over `n : ℕ` with index `n + 1`. `HasSum` is used throughout (never
  `tsum`). Natural subtraction only appears as `k - i` with `i ≤ k`, `p - 1`, `q - 1`, `m - 1`
  with `p, q, m ≥ 1`.
* The fallback statements for `k = 4, 5` (stmt_7, stmt_8) use real numbers and existential
  rational coefficients chosen before the angle; they do not depend on the formal
  construction. For `k = 4` the monomials of weight `≤ 4` are `1, 𝔷2, 𝔷3, 𝔷4, 𝔷2²`; for `k = 5`
  also `𝔷5, 𝔷2𝔷3`. (They are the cases `k = 4, 5` of stmt_10; since `ζ(4) = (2/5)ζ(2)²`, the
  real-number form is weaker than the formal one.)
* The final statement (stmt_10) is `thm:flat` with an existential family of formal polynomials
  `μ_{k,n}` independent of `α` (the paper's "where each `μ_{k,n}` is a polynomial ..."), the
  formal shapes of `μ_{2,n}`, `μ_{3,n}` (this is the meaning the supplement gives to "the
  coefficient of `ζ(2)`, `ζ(3)`", line 23), and the two power-series consequences with
  `r_n = [1]μ_{2,n}`, `r'_n = [1]μ_{3,n}`, `q'_n = [𝔷2]μ_{3,n}`, `Ψ = Σ q'_n c^{2n}`.

## Statements

* stmt_1: `log Γ(1+x) = -γx + Σ_{w≥2} (-1)^w ζ(w) x^w / w` for `|x| < 1` (DLMF 5.7, used in the
  proof of `g6:lem:gamma`).
* stmt_2: polygamma at integers: `ψ(m) = -γ + H_{m-1}` and, for `j ≥ 1`,
  `ψ^{(j)}(m) = (-1)^{j+1} j! (ζ(j+1) - H^{(j+1)}_{m-1})`, `m ≥ 1`.
* stmt_3: `g6:lem:En` for every `k`: each `e^{(k)}_{j,i,m}(n)` is a rational polynomial in `n`,
  vanishes unless `j ≤ k` and `i + m = k`, and `∂_u^k E_n(0,y)/k! = E_{n,k}(y, log y, log(1-y))`.
* stmt_4: `g6:lem:gamma` for every `(i,m)`: for integers `p, q ≥ 1`,
  `∫_0^1 y^{p-1}(1-y)^{q-1} log^i y log^m(1-y) dy = B(p,q) D^♯_{i,m}(p,q)|_{𝔷=ζ}`, and
  `D^♯_{i,m}(p,q)` has weight `≤ i + m` (only `𝔷_2, …, 𝔷_{i+m}` occur).
* stmt_5: `g6:prop:formula` for every `k`: `Σ_{n≥1} P_{k,n}(ζ) c^{2n}` converges absolutely to
  `M_k(α)`, and every `P_{k,n}` has weight `≤ k`.
* stmt_6: the coefficient bound: every coefficient of `P_{k,n}` is at most
  `C_k n^k (1 + log(2n+k))^k` in absolute value.
* stmt_7: `thm:flat` for `k = 4` (real form): rationals independent of `α` with
  `M_4(α) = Σ (a_n + b_n ζ(2) + c_n ζ(3) + d_n ζ(4) + e_n ζ(2)²) c^{2n}`.
* stmt_8: `thm:flat` for `k = 5` (real form), monomials `1, ζ2, ζ3, ζ4, ζ5, ζ2², ζ2ζ3`.
* stmt_9: first rung of `thm:zladder`, formal sense: `[𝔷3]P_{k+1,n} = -[𝔷2]P_{k,n}` for every
  `k ≥ 0`, `n ≥ 1`; and the series `[ζ(2)]M_k(α) = Σ [𝔷2]P_{k,n} c^{2n}` and
  `[ζ(3)]M_{k+1}(α) = Σ [𝔷3]P_{k+1,n} c^{2n}` converge, to opposite values.
* stmt_10: `thm:flat` in full (see the last bullet of the formalization choices).

## Numerical check before freezing

A Python/SymPy script mirroring these definitions exactly (the same truncations, exact
rationals) gave: the weight predicate holds for `P_{k,n}`, `k ≤ 6`, `n ≤ 4`; every monomial of
`E_{n,k}` has `j ≤ k`, `i + m = k`; `P_{3,1} = 35/72 - 𝔷2/36 - 𝔷3/3` (as in the supplement and
V05); `[𝔷3]P_{k+1,n} = -[𝔷2]P_{k,n}` exactly for `k ≤ 5`, `n ≤ 4`; `∂_u^k E_n(0,y)/k!`
(mpmath numerical differentiation, `y = 0.3`, `n = 1, 3`, `k = 4, 5`) agrees with
`E_{n,k}(y, log y, log(1-y))`; and `-(a_n/2)∫_0^1 [y(1-y)]^{n-1}E_{n,k} dy` agrees with
`P_{k,n}(ζ)` for `k = 4, 5`, `n = 1, 2` (differences `< 1e-30`). Examples:
`P_{2,2} = -103/864 + 𝔷2/15`, `P_{4,1} = -571/1080 + 103/432 𝔷2 + 𝔷3/36 + 4/15 𝔷4 - 𝔷2²/15`,
`([𝔷3]P_{4,2}, [𝔷2]P_{3,2}) = (-43/900, 43/900)` (as in the ladder supplement). `M_4(π/2)` and
`M_5(π/2)` from the defining integral (Cauchy integral in `u`, radius 0.2) agree with
`Σ_{n<45} P_{k,n}(ζ) c^{2n}` to 15 digits (`M_4(π/2) = 0.00602799311260912`,
`M_5(π/2) = 0.00316164871041226`).

## NOT COVERED

* The general rungs of `thm:zladder` (every product `z` of zeta values; the `ζ(5)` example;
  the precise form `g6:thm:ladderprecise` with `γ_e(z)`), and Lemmas `g6:lem:op`, `g6:lem:DE`.
* `g6:eq:iso` in its own form (weight exactly `i + m` with formal harmonic variables `h^{(r)}`
  given weight `r`): only its consequence, weight `≤ i + m` in the zeta variables, is stated.
* `D_{i,m} = ∂_p^i∂_q^m B/B` as derivatives of the Gamma quotient in real `(p,q)`: only the
  integral form (which is what the proof of `g6:prop:formula` uses) is stated.
* The remark that the split into zeta parts of the real number `μ_{3,n}` would need the linear
  independence of `1, π², ζ(3)` (open); everything is stated for the formal coefficient.
* `Ψ ∈ ℚ[[c²]]` is recorded through its rational coefficients `q'_n` and the convergence of
  `Σ q'_n c^{2n}`.
-/
-- BEGIN STATEMENTS
namespace Q04

open Real MeasureTheory Set Filter

/-- Truncated exponential `Σ_{r=0}^K L^r / r!` in a `ℚ`-algebra. -/
noncomputable def expT {A : Type*} [CommRing A] [Algebra ℚ A] (K : ℕ) (L : A) : A :=
  ∑ r ∈ Finset.range (K + 1), ((r.factorial : ℚ)⁻¹) • L ^ r

/-- Truncated logarithm `log(1+w) ↦ Σ_{m=1}^K (-1)^{m-1} w^m / m` in a `ℚ`-algebra. -/
noncomputable def logT {A : Type*} [CommRing A] [Algebra ℚ A] (K : ℕ) (w : A) : A :=
  ∑ m ∈ Finset.Icc 1 K, ((-1 : ℚ) ^ (m - 1) / (m : ℚ)) • w ^ m

/-- The kernel `E_n = exp(uA + (1-2n) log(1 + y(e^{uℓ}-1)))`, `ℓ = b - a`,
`A = -(n-1)a + nb`, as a polynomial in `u` over `ℚ[y,a,b]` (`y = X 0`, `a = X 1`, `b = X 2`),
computed modulo terms that do not affect the coefficients of `u^0, …, u^K`. -/
noncomputable def Eser (n K : ℕ) : Polynomial (MvPolynomial (Fin 3) ℚ) :=
  let yv : MvPolynomial (Fin 3) ℚ := MvPolynomial.X 0
  let av : MvPolynomial (Fin 3) ℚ := MvPolynomial.X 1
  let bv : MvPolynomial (Fin 3) ℚ := MvPolynomial.X 2
  let ℓ := bv - av
  let A := (-((n : ℚ) - 1)) • av + (n : ℚ) • bv
  let w := Polynomial.C yv *
    ∑ r ∈ Finset.Icc 1 K, ((r.factorial : ℚ)⁻¹) • (Polynomial.C ℓ * Polynomial.X) ^ r
  expT K (Polynomial.C A * Polynomial.X + (1 - 2 * (n : ℚ)) • logT K w)

/-- `E_{n,k} ∈ ℚ[y,a,b]`, the coefficient of `u^k` in `E_n` (Lemma `g6:lem:En`). -/
noncomputable def Ecoef (k n : ℕ) : MvPolynomial (Fin 3) ℚ := (Eser n k).coeff k

/-- `e^{(k)}_{j,i,m}(n)`, the coefficient of `y^j a^i b^m` in `E_{n,k}`. -/
noncomputable def eF (k n j i m : ℕ) : ℚ :=
  (Ecoef k n).coeff (Finsupp.single 0 j + Finsupp.single 1 i + Finsupp.single 2 m)

/-- `Δ_⋆` truncated at degree `K`, with `h^{(r)}_1, h^{(r)}_2, h^{(r)}_3` replaced by
`H^{(r)}_{p-1}, H^{(r)}_{q-1}, H^{(r)}_{p+q-1}`; `x = X 0`, `x' = X 1`, coefficients in the
formal zeta values `ℚ[𝔷_2, 𝔷_3, …]` (`𝔷_w = X w`). -/
noncomputable def DeltaT (K p q : ℕ) : MvPolynomial (Fin 2) (MvPolynomial ℕ ℚ) :=
  let x : MvPolynomial (Fin 2) (MvPolynomial ℕ ℚ) := MvPolynomial.X 0
  let x' : MvPolynomial (Fin 2) (MvPolynomial ℕ ℚ) := MvPolynomial.X 1
  (∑ r ∈ Finset.Icc 1 K, ((-1 : ℚ) ^ (r - 1) / (r : ℚ)) •
      (V05.HrQ r (p - 1) • x ^ r + V05.HrQ r (q - 1) • x' ^ r
        - V05.HrQ r (p + q - 1) • (x + x') ^ r)) +
  ∑ w ∈ Finset.Icc 2 K, ((-1 : ℚ) ^ (w + 1) / (w : ℚ)) •
      (MvPolynomial.C (MvPolynomial.X w) * ((x + x') ^ w - x ^ w - x' ^ w))

/-- `D^♯_{i,m}(p,q) = i! m! [x^i x'^m] exp(Δ_⋆) ∈ ℚ[𝔷_2, 𝔷_3, …]`. -/
noncomputable def DsharpF (i m p q : ℕ) : MvPolynomial ℕ ℚ :=
  ((i.factorial * m.factorial : ℕ) : ℚ) •
    (expT (i + m) (DeltaT (i + m) p q)).coeff (Finsupp.single 0 i + Finsupp.single 1 m)

/-- The formal coefficient `P_{k,n} ∈ ℚ[𝔷_2, …, 𝔷_k]` (Definition `g6:def:formal`). -/
noncomputable def Pf (k n : ℕ) : MvPolynomial ℕ ℚ :=
  (-(1 / ((n : ℚ) * (2 * (n : ℚ) - 1)))) •
    ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
      (eF k n j i (k - i) * V05.rhoQ j n) • DsharpF i (k - i) (n + j) n

/-- "A polynomial with rational coefficients in `𝔷_2, …, 𝔷_k` of weight at most `k`": every
monomial `∏ 𝔷_w^{e_w}` has all `w ≥ 2` and weight `Σ e_w w ≤ k`. -/
def WeightLE (k : ℕ) (P : MvPolynomial ℕ ℚ) : Prop :=
  ∀ s ∈ P.support, (∀ w ∈ s.support, 2 ≤ w) ∧ (s.sum fun w e => e * w) ≤ k

/-- DLMF 5.7 (proof of `g6:lem:gamma`): for `|x| < 1`,
`log Γ(1+x) = -γ x + Σ_{w≥2} (-1)^w ζ(w) x^w / w`. -/
def stmt_1 : Prop :=
  ∀ x : ℝ, |x| < 1 →
    HasSum (fun w : ℕ => (-1 : ℝ) ^ (w + 2) * X05.z (w + 2) * x ^ (w + 2) / ((w : ℝ) + 2))
      (Real.log (Real.Gamma (1 + x)) + Real.eulerMascheroniConstant * x)

/-- Polygamma values at integers `m ≥ 1`: `ψ(m) = -γ + H_{m-1}` and, for `j ≥ 1`,
`ψ^{(j)}(m) = (-1)^{j+1} j! (ζ(j+1) - H^{(j+1)}_{m-1})`, with `ψ^{(j)} = (log Γ)^{(j+1)}`. -/
def stmt_2 : Prop :=
  ∀ m : ℕ, 1 ≤ m →
    deriv (fun x => Real.log (Real.Gamma x)) m = -Real.eulerMascheroniConstant + V05.Hr 1 (m - 1) ∧
    ∀ j : ℕ, 1 ≤ j →
      iteratedDeriv (j + 1) (fun x => Real.log (Real.Gamma x)) m =
        (-1 : ℝ) ^ (j + 1) * (j.factorial : ℝ) * (X05.z (j + 1) - V05.Hr (j + 1) (m - 1))

/-- `g6:lem:En` for every `k`: `e^{(k)}_{j,i,m}(n) ∈ ℚ[n]`, it vanishes unless `j ≤ k` and
`i + m = k`, and `∂_u^k E_n(0,y)/k! = E_{n,k}(y, log y, log(1-y))` for `n ≥ 1`, `0 < y < 1`. -/
def stmt_3 : Prop :=
  ∀ k : ℕ,
    (∀ j i m : ℕ, ∃ Q : Polynomial ℚ, ∀ n : ℕ, eF k n j i m = Q.eval (n : ℚ)) ∧
    (∀ n j i m : ℕ, eF k n j i m ≠ 0 → j ≤ k ∧ i + m = k) ∧
    ∀ n : ℕ, 1 ≤ n → ∀ y : ℝ, 0 < y → y < 1 →
      iteratedDeriv k (fun u => V05.E n u y) 0 / (k.factorial : ℝ) =
        ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
          (eF k n j i (k - i) : ℝ) * y ^ j * Real.log y ^ i * Real.log (1 - y) ^ (k - i)

/-- `g6:lem:gamma` for every `(i,m)`: for integers `p, q ≥ 1`,
`∂_p^i∂_q^m B(p,q) = ∫_0^1 y^{p-1}(1-y)^{q-1} log^i y log^m(1-y) dy = B(p,q) D^♯_{i,m}(p,q)` at
`𝔷_w = ζ(w)`, and `D^♯_{i,m}(p,q)` is a rational polynomial in `𝔷_2, …, 𝔷_{i+m}` of weight
`≤ i + m`. -/
def stmt_4 : Prop :=
  ∀ p q : ℕ, 1 ≤ p → 1 ≤ q → ∀ i m : ℕ,
    IntervalIntegrable (fun y : ℝ =>
      y ^ (p - 1) * (1 - y) ^ (q - 1) * Real.log y ^ i * Real.log (1 - y) ^ m) volume 0 1 ∧
    (∫ y in (0 : ℝ)..1, y ^ (p - 1) * (1 - y) ^ (q - 1) * Real.log y ^ i * Real.log (1 - y) ^ m)
      = V05.Bet p q * MvPolynomial.aeval X05.z (DsharpF i m p q) ∧
    WeightLE (i + m) (DsharpF i m p q)

/-- `g6:prop:formula` for every `k`: for `0 < α < π`, `M_k(α) = Σ_{n≥1} P_{k,n}(ζ(2),…,ζ(k)) c^{2n}`,
absolutely convergent; every monomial of `P_{k,n}` has weight `≤ k`. -/
def stmt_5 : Prop :=
  (∀ k n : ℕ, WeightLE k (Pf k n)) ∧
  ∀ α : ℝ, 0 < α → α < π → ∀ k : ℕ,
    Summable (fun n : ℕ => |MvPolynomial.aeval X05.z (Pf k (n + 1)) * D04.c α ^ (2 * (n + 1))|) ∧
    HasSum (fun n : ℕ => MvPolynomial.aeval X05.z (Pf k (n + 1)) * D04.c α ^ (2 * (n + 1)))
      (V05.Mk k α)

/-- Proof of `g6:thm:flat` (line 185): every coefficient of `P_{k,n}` is at most
`C_k n^k (1 + log(2n+k))^k` in absolute value. -/
def stmt_6 : Prop :=
  ∀ k : ℕ, ∃ C : ℝ, ∀ n : ℕ, 1 ≤ n → ∀ s : ℕ →₀ ℕ,
    |(((Pf k n).coeff s : ℚ) : ℝ)| ≤
      C * (n : ℝ) ^ k * (1 + Real.log (2 * (n : ℝ) + k)) ^ k

/-- `thm:flat` for `k = 4` (real form): there are rationals, independent of the angle, with
`M_4(α) = Σ_{n≥1} (a_n + b_n ζ(2) + c_n ζ(3) + d_n ζ(4) + e_n ζ(2)²) c^{2n}` for `0 < α < π`
(the monomials of weight `≤ 4`). -/
def stmt_7 : Prop :=
  ∃ a b c d e : ℕ → ℚ, ∀ α : ℝ, 0 < α → α < π →
    HasSum (fun n : ℕ => ((a (n + 1) : ℝ) + (b (n + 1) : ℝ) * X05.z 2 + (c (n + 1) : ℝ) * X05.z 3
        + (d (n + 1) : ℝ) * X05.z 4 + (e (n + 1) : ℝ) * X05.z 2 ^ 2) * D04.c α ^ (2 * (n + 1)))
      (V05.Mk 4 α)

/-- `thm:flat` for `k = 5` (real form): rationals, independent of the angle, with
`M_5(α) = Σ_{n≥1} (a_n + b_n ζ(2) + c_n ζ(3) + d_n ζ(4) + e_n ζ(5) + f_n ζ(2)² + g_n ζ(2)ζ(3))
c^{2n}` for `0 < α < π` (the monomials of weight `≤ 5`). -/
def stmt_8 : Prop :=
  ∃ a b c d e f g : ℕ → ℚ, ∀ α : ℝ, 0 < α → α < π →
    HasSum (fun n : ℕ => ((a (n + 1) : ℝ) + (b (n + 1) : ℝ) * X05.z 2 + (c (n + 1) : ℝ) * X05.z 3
        + (d (n + 1) : ℝ) * X05.z 4 + (e (n + 1) : ℝ) * X05.z 5 + (f (n + 1) : ℝ) * X05.z 2 ^ 2
        + (g (n + 1) : ℝ) * (X05.z 2 * X05.z 3)) * D04.c α ^ (2 * (n + 1)))
      (V05.Mk 5 α)

/-- First rung of `thm:zladder`, in the supplement's formal sense: `[𝔷3]P_{k+1,n} = -[𝔷2]P_{k,n}`
for every `k ≥ 0` and `n ≥ 1`; and for `0 < α < π` the series
`[ζ(2)]M_k(α) = Σ_n [𝔷2]P_{k,n} c^{2n}` converges and `[ζ(3)]M_{k+1}(α) = Σ_n [𝔷3]P_{k+1,n} c^{2n}`
converges to its negative. -/
def stmt_9 : Prop :=
  (∀ k n : ℕ, 1 ≤ n →
    (Pf (k + 1) n).coeff (Finsupp.single 3 1) =
      -(Pf k n).coeff (Finsupp.single 2 1)) ∧
  ∀ k : ℕ, ∀ α : ℝ, 0 < α → α < π → ∃ S : ℝ,
    HasSum (fun n : ℕ => (((Pf k (n + 1)).coeff (Finsupp.single 2 1) : ℚ) : ℝ)
      * D04.c α ^ (2 * (n + 1))) S ∧
    HasSum (fun n : ℕ => (((Pf (k + 1) (n + 1)).coeff (Finsupp.single 3 1) : ℚ) : ℝ)
      * D04.c α ^ (2 * (n + 1))) (-S)

/-- `thm:flat` in full: there are formal polynomials `μ_{k,n} ∈ ℚ[𝔷_2, …, 𝔷_k]` of weight `≤ k`,
independent of the angle, with `M_k(α) = Σ_{n≥1} μ_{k,n}(ζ) c^{2n}` for every `k ≥ 0` and
`0 < α < π`; `μ_{2,n} = r_n + 𝔷2/(4n²-1)` and `μ_{3,n} = r'_n + q'_n 𝔷2 - 𝔷3/(4n²-1)` with
`r_n, r'_n, q'_n ∈ ℚ`; and with `Ψ = Σ q'_n c^{2n}`, `M_2 - ζ(2)Φ = Σ r_n c^{2n}` and
`M_3 + ζ(3)Φ - ζ(2)Ψ = Σ r'_n c^{2n}`. -/
def stmt_10 : Prop :=
  ∃ μ : ℕ → ℕ → MvPolynomial ℕ ℚ,
    (∀ k n : ℕ, WeightLE k (μ k n)) ∧
    (∀ n : ℕ, 1 ≤ n → ∃ r r' q' : ℚ,
      μ 2 n = MvPolynomial.C r + MvPolynomial.C (1 / (4 * (n : ℚ) ^ 2 - 1)) * MvPolynomial.X 2 ∧
      μ 3 n = MvPolynomial.C r' + MvPolynomial.C q' * MvPolynomial.X 2
        - MvPolynomial.C (1 / (4 * (n : ℚ) ^ 2 - 1)) * MvPolynomial.X 3) ∧
    ∀ α : ℝ, 0 < α → α < π →
      (∀ k : ℕ, HasSum (fun n : ℕ => MvPolynomial.aeval X05.z (μ k (n + 1)) * D04.c α ^ (2 * (n + 1)))
        (V05.Mk k α)) ∧
      HasSum (fun n : ℕ => (((μ 2 (n + 1)).coeff 0 : ℚ) : ℝ) * D04.c α ^ (2 * (n + 1)))
        (V05.Mk 2 α - X05.z 2 * V05.Phi (D04.c α)) ∧
      ∃ Ψ : ℝ,
        HasSum (fun n : ℕ => (((μ 3 (n + 1)).coeff (Finsupp.single 2 1) : ℚ) : ℝ)
          * D04.c α ^ (2 * (n + 1))) Ψ ∧
        HasSum (fun n : ℕ => (((μ 3 (n + 1)).coeff 0 : ℚ) : ℝ) * D04.c α ^ (2 * (n + 1)))
          (V05.Mk 3 α + X05.z 3 * V05.Phi (D04.c α) - X05.z 2 * Ψ)

end Q04
-- END STATEMENTS

open Real Filter Topology

namespace Q04

lemma sum_log_shift {x : ℝ} (hx : |x| < 1) (n : ℕ) :
    ∑ m ∈ Finset.range (n + 1), Real.log (1 + x + m) =
      Real.log ((n + 1).factorial : ℝ) +
        ∑ m ∈ Finset.range (n + 1), Real.log (1 + x / ((m : ℝ) + 1)) := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [Finset.sum_range_succ, ih,
      Finset.sum_range_succ (fun m : ℕ => Real.log (1 + x / ((m : ℝ) + 1))) (n + 1),
      Nat.factorial_succ (n + 1)]
    have h1 : (0 : ℝ) < (n : ℝ) + 1 + 1 := by positivity
    have h2 := X05.fac_pos hx (n + 1)
    have hf : (0 : ℝ) < ((n + 1).factorial : ℝ) := by exact_mod_cast (n + 1).factorial_pos
    push_cast at h2 ⊢
    rw [Real.log_mul h1.ne' hf.ne',
      show 1 + x + ((n : ℝ) + 1) = ((n : ℝ) + 1 + 1) * (1 + x / ((n : ℝ) + 1 + 1)) by
        field_simp; ring, Real.log_mul h1.ne' h2.ne']
    ring

lemma T_hasSum {x : ℝ} (hx : |x| < 1) (k : ℕ) :
    HasSum (fun w : ℕ => (-1 : ℝ) ^ (w + 2) * (1 / ((k : ℝ) + 1) ^ (w + 2)) * x ^ (w + 2)
        / ((w : ℝ) + 2))
      (x / ((k : ℝ) + 1) - Real.log (1 + x / ((k : ℝ) + 1))) := by
  have hk : (0 : ℝ) < (k : ℝ) + 1 := by positivity
  have sm : |-(x / ((k : ℝ) + 1))| < 1 := by
    have h1 : (1 : ℝ) ≤ (k : ℝ) + 1 := by have := k.cast_nonneg (α := ℝ); linarith
    rw [abs_neg, abs_div, abs_of_pos hk, div_lt_one hk]; linarith
  have L := Real.hasSum_pow_div_log_of_abs_lt_one sm
  rw [← hasSum_nat_add_iff' 1] at L
  simp only [Finset.range_one, Finset.sum_singleton, sub_neg_eq_add] at L
  convert L using 1
  · ext w
    simp only [neg_pow (_ / _), div_pow]
    push_cast
    ring
  · simp
    ring

end Q04

theorem Q04.proof_1 : Q04.stmt_1 := by
  intro x hx
  set f : ℕ × ℕ → ℝ := fun p => (-1 : ℝ) ^ (p.2 + 2) * (1 / ((p.1 : ℝ) + 1) ^ (p.2 + 2)) *
      x ^ (p.2 + 2) / ((p.2 : ℝ) + 2) with hfdef
  have hfs : Summable f := by
    have hg1 : Summable (fun k : ℕ => 1 / ((k : ℝ) + 1) ^ 2) :=
      (X05.z_hasSum (by norm_num)).summable
    have hg2 : Summable (fun j : ℕ => |x| ^ j) :=
      summable_geometric_of_lt_one (abs_nonneg x) hx
    refine Summable.of_norm_bounded (hg1.mul_of_nonneg hg2 (fun k => by positivity)
      (fun j => by positivity)) ?_
    rintro ⟨k, j⟩
    simp only [hfdef, Real.norm_eq_abs]
    have hk1 : (1 : ℝ) ≤ (k : ℝ) + 1 := by have := k.cast_nonneg (α := ℝ); linarith
    have hj : (1 : ℝ) ≤ (j : ℝ) + 2 := by have := j.cast_nonneg (α := ℝ); linarith
    rw [abs_div, abs_mul, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
      abs_of_pos (by positivity : (0 : ℝ) < 1 / ((k : ℝ) + 1) ^ (j + 2)),
      abs_of_pos (by positivity : (0 : ℝ) < (j : ℝ) + 2), abs_pow]
    have b1 : 1 / ((k : ℝ) + 1) ^ (j + 2) ≤ 1 / ((k : ℝ) + 1) ^ 2 := by
      apply one_div_le_one_div_of_le (by positivity)
      exact pow_le_pow_right₀ hk1 (by omega)
    have b2 : |x| ^ (j + 2) ≤ |x| ^ j :=
      pow_le_pow_of_le_one (abs_nonneg x) hx.le (by omega)
    calc 1 / ((k : ℝ) + 1) ^ (j + 2) * |x| ^ (j + 2) / ((j : ℝ) + 2)
        ≤ 1 / ((k : ℝ) + 1) ^ (j + 2) * |x| ^ (j + 2) / 1 := by
          gcongr
      _ ≤ 1 / ((k : ℝ) + 1) ^ 2 * |x| ^ j := by
          rw [div_one]; gcongr
  have hf := hfs.hasSum
  have hT := hf.prod_fiberwise (fun k => Q04.T_hasSum hx k)
  have hf' : HasSum (f ∘ Equiv.prodComm ℕ ℕ) (∑' p, f p) :=
    (Equiv.prodComm ℕ ℕ).hasSum_iff.mpr hf
  have hZ := hf'.prod_fiberwise (fun j => by
    have := (((X05.z_hasSum (show 1 < j + 2 by omega)).mul_left ((-1 : ℝ) ^ (j + 2))).mul_right
      (x ^ (j + 2))).div_const ((j : ℝ) + 2)
    exact this)
  -- identify the sum
  have hx0 : 0 < 1 + x := by cases abs_lt.mp hx; linarith
  have hlim1 := (tendsto_add_atTop_iff_nat 1).mpr hT.tendsto_sum_nat
  have hH : Tendsto (fun n : ℕ => ((harmonic (n + 1) : ℚ) : ℝ) - Real.log ((n : ℝ) + 1)) atTop
      (𝓝 Real.eulerMascheroniConstant) := by
    have := (tendsto_add_atTop_iff_nat 1).mpr Real.tendsto_harmonic_sub_log
    simpa using this
  have hL : Tendsto (fun n : ℕ => Real.log (n : ℝ) - Real.log ((n : ℝ) + 1)) atTop (𝓝 0) := by
    have := (Real.tendsto_log_comp_add_sub_log 1).comp tendsto_natCast_atTop_atTop
    simpa using this.neg
  have hlim := ((hH.const_mul (-x)).add (hL.const_mul (1 + x))).add hlim1
  have hG := Real.BohrMollerup.tendsto_log_gamma hx0
  have heq : ∀ n : ℕ, 1 ≤ n → Real.BohrMollerup.logGammaSeq (1 + x) n =
      -x * (((harmonic (n + 1) : ℚ) : ℝ) - Real.log ((n : ℝ) + 1)) +
        (1 + x) * (Real.log (n : ℝ) - Real.log ((n : ℝ) + 1)) +
        ∑ k ∈ Finset.range (n + 1), (x / ((k : ℝ) + 1) - Real.log (1 + x / ((k : ℝ) + 1))) := by
    intro n hn
    unfold Real.BohrMollerup.logGammaSeq
    rw [Q04.sum_log_shift hx n, Nat.factorial_succ]
    have hf : (0 : ℝ) < (n.factorial : ℝ) := by exact_mod_cast n.factorial_pos
    have h1 : (0 : ℝ) < (n : ℝ) + 1 := by positivity
    push_cast
    have hh : ∑ k ∈ Finset.range (n + 1), x / ((k : ℝ) + 1) =
        x * ((harmonic (n + 1) : ℚ) : ℝ) := by
      simp only [harmonic]
      push_cast
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl fun k _ => ?_
      ring
    rw [Real.log_mul h1.ne' hf.ne', Finset.sum_sub_distrib, hh]
    ring
  have hG' : Tendsto (fun n : ℕ => -x * (((harmonic (n + 1) : ℚ) : ℝ) - Real.log ((n : ℝ) + 1)) +
        (1 + x) * (Real.log (n : ℝ) - Real.log ((n : ℝ) + 1)) +
        ∑ k ∈ Finset.range (n + 1), (x / ((k : ℝ) + 1) - Real.log (1 + x / ((k : ℝ) + 1))))
      atTop (𝓝 (Real.log (Real.Gamma (1 + x)))) := by
    refine hG.congr' ?_
    filter_upwards [eventually_ge_atTop 1] with n hn
    exact heq n hn
  have hu := tendsto_nhds_unique hlim hG'
  have hS : ∑' p, f p = Real.log (Real.Gamma (1 + x)) + Real.eulerMascheroniConstant * x := by
    rw [← hu]; ring
  rw [hS] at hZ
  convert hZ using 1

namespace Q04

/-- the Taylor coefficients of `log Γ(1+y)` at `0` -/
noncomputable def aL : ℕ → ℝ
  | 0 => 0
  | 1 => -Real.eulerMascheroniConstant
  | (w + 2) => (-1 : ℝ) ^ (w + 2) * X05.z (w + 2) / ((w : ℝ) + 2)

lemma F_hasSum {y : ℝ} (hy : |y| < 1) :
    HasSum (fun n => aL n * y ^ n) (Real.log (Real.Gamma (1 + y))) := by
  rw [← hasSum_nat_add_iff' 2]
  have := Q04.proof_1 y hy
  convert this using 1
  · ext w; simp only [aL]; ring
  · simp [Finset.sum_range_succ, aL]

lemma aL_bd (n : ℕ) : |aL n| ≤ Real.eulerMascheroniConstant + 2 := by
  have hγ : 0 ≤ Real.eulerMascheroniConstant := by
    linarith [Real.one_half_lt_eulerMascheroniConstant]
  match n with
  | 0 => simp only [aL, abs_zero]; positivity
  | 1 => simp only [aL, abs_neg, abs_of_nonneg hγ]; linarith
  | (w + 2) =>
    simp only [aL]
    have h1 : (1 : ℝ) ≤ (w : ℝ) + 2 := by have := w.cast_nonneg (α := ℝ); linarith
    rw [abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
      abs_of_pos (by positivity : (0 : ℝ) < (w : ℝ) + 2)]
    have := X05.z_bd (show 2 ≤ w + 2 by omega)
    calc |X05.z (w + 2)| / ((w : ℝ) + 2) ≤ |X05.z (w + 2)| := div_le_self (abs_nonneg _) h1
      _ ≤ _ := by linarith

lemma F_series : HasFPowerSeriesOnBall (fun y => Real.log (Real.Gamma (1 + y)))
    (FormalMultilinearSeries.ofScalars ℝ aL) 0 1 where
  r_le := by
    refine ENNReal.le_of_forall_nnreal_lt fun r hr => ?_
    refine FormalMultilinearSeries.le_radius_of_bound _ (Real.eulerMascheroniConstant + 2)
      fun n => ?_
    rw [FormalMultilinearSeries.ofScalars_norm]
    have hr1 : (r : ℝ) < 1 := by exact_mod_cast hr
    calc ‖aL n‖ * (r : ℝ) ^ n ≤ (Real.eulerMascheroniConstant + 2) * 1 :=
          mul_le_mul (by simpa [Real.norm_eq_abs] using aL_bd n) (pow_le_one₀ r.2 hr1.le)
            (by positivity) (by linarith [Real.one_half_lt_eulerMascheroniConstant])
      _ = _ := mul_one _
  r_pos := one_pos
  hasSum := by
    intro y hy
    have hy' : |y| < 1 := by
      simpa [Metric.mem_eball, edist_dist, Real.dist_eq] using hy
    simp only [FormalMultilinearSeries.ofScalars_apply_eq, smul_eq_mul, zero_add]
    exact F_hasSum hy'

/-- `g = log ∘ Γ` -/
lemma g_cd1 (k : ℕ) : ContDiffAt ℝ k (fun x => Real.log (Real.Gamma x)) 1 := by
  have hF : ContDiffAt ℝ k (fun y => Real.log (Real.Gamma (1 + y))) ((fun x : ℝ => x - 1) 1) := by
    simp only [sub_self]; exact F_series.analyticAt.contDiffAt
  have hs : ContDiffAt ℝ k (fun x : ℝ => x - 1) 1 := contDiffAt_id.sub contDiffAt_const
  have := ContDiffAt.comp (g := fun y => Real.log (Real.Gamma (1 + y))) (f := fun x : ℝ => x - 1) (1 : ℝ) hF hs
  convert this using 1
  ext x; simp

lemma g_iter1 (k : ℕ) :
    iteratedDeriv k (fun x => Real.log (Real.Gamma x)) 1 = (k.factorial : ℝ) * aL k := by
  have h1 := F_series.factorial_smul 1 k
  simp only [FormalMultilinearSeries.ofScalars_apply_eq, smul_eq_mul, one_pow, mul_one,
    nsmul_eq_mul] at h1
  rw [h1, ← iteratedDeriv_eq_iteratedFDeriv,
    congrFun (iteratedDeriv_comp_const_add (n := k) (f := fun x => Real.log (Real.Gamma x))
      (s := 1)) 0]
  simp

lemma g_feq {x : ℝ} (hx : 0 < x) :
    Real.log (Real.Gamma (x + 1)) = Real.log (Real.Gamma x) + Real.log x := by
  rw [Real.Gamma_add_one hx.ne', Real.log_mul hx.ne' (Real.Gamma_pos_of_pos hx).ne']
  ring

lemma g_cd (k n : ℕ) : ContDiffAt ℝ k (fun x => Real.log (Real.Gamma x)) ((n : ℝ) + 1) := by
  induction n with
  | zero => simpa using g_cd1 k
  | succ n ih =>
    have hev : (fun x => Real.log (Real.Gamma x)) =ᶠ[𝓝 (((n + 1 : ℕ) : ℝ) + 1)]
        (fun x => Real.log (Real.Gamma (x - 1)) + Real.log (x - 1)) := by
      have : ∀ᶠ x in 𝓝 (((n + 1 : ℕ) : ℝ) + 1), (1 : ℝ) < x :=
        eventually_gt_nhds (by push_cast; linarith [n.cast_nonneg (α := ℝ)])
      filter_upwards [this] with x hx
      have := g_feq (show 0 < x - 1 by linarith)
      rw [sub_add_cancel] at this
      exact this
    refine ContDiffAt.congr_of_eventuallyEq ?_ hev
    have hp : ((n + 1 : ℕ) : ℝ) + 1 - 1 = (n : ℝ) + 1 := by push_cast; ring
    refine ContDiffAt.add ?_ ?_
    · have hs : ContDiffAt ℝ k (fun x : ℝ => x - 1) (((n + 1 : ℕ) : ℝ) + 1) :=
        contDiffAt_id.sub contDiffAt_const
      have := ContDiffAt.comp (g := fun x => Real.log (Real.Gamma x)) (f := fun x : ℝ => x - 1) _
        (show ContDiffAt ℝ k (fun x => Real.log (Real.Gamma x))
          ((fun x : ℝ => x - 1) (((n + 1 : ℕ) : ℝ) + 1)) by simp only []; rw [hp]; exact ih) hs
      exact this
    · have hs : ContDiffAt ℝ k (fun x : ℝ => x - 1) (((n + 1 : ℕ) : ℝ) + 1) :=
        contDiffAt_id.sub contDiffAt_const
      have h0 : (0 : ℝ) < ((n + 1 : ℕ) : ℝ) + 1 - 1 := by rw [hp]; positivity
      exact ContDiffAt.comp (g := Real.log) (f := fun x : ℝ => x - 1) _
        (show ContDiffAt ℝ k Real.log ((fun x : ℝ => x - 1) (((n + 1 : ℕ) : ℝ) + 1)) from
          Real.contDiffAt_log.mpr h0.ne') hs

lemma log_iter {k : ℕ} {x : ℝ} (_hx : 0 < x) :
    iteratedDeriv (k + 1) Real.log x = (-1 : ℝ) ^ k * (k.factorial : ℝ) / x ^ (k + 1) := by
  rw [iteratedDeriv_succ', Real.deriv_log', iteratedDeriv_eq_iterate, iter_deriv_inv]
  rw [show (-1 - (k : ℤ)) = -((k + 1 : ℕ) : ℤ) by push_cast; ring, zpow_neg, zpow_natCast]
  ring

lemma g_step (k n : ℕ) :
    iteratedDeriv (k + 1) (fun x => Real.log (Real.Gamma x)) ((n : ℝ) + 1 + 1) =
      iteratedDeriv (k + 1) (fun x => Real.log (Real.Gamma x)) ((n : ℝ) + 1) +
        (-1 : ℝ) ^ k * (k.factorial : ℝ) / ((n : ℝ) + 1) ^ (k + 1) := by
  have hpos : (0 : ℝ) < (n : ℝ) + 1 := by positivity
  have h1 : iteratedDeriv (k + 1) (fun x => Real.log (Real.Gamma (x + 1))) ((n : ℝ) + 1) =
      iteratedDeriv (k + 1) (fun x => Real.log (Real.Gamma x)) ((n : ℝ) + 1 + 1) := by
    rw [congrFun (iteratedDeriv_comp_add_const (n := k + 1)
      (f := fun x => Real.log (Real.Gamma x)) (s := 1)) _]
  have hev : (fun x => Real.log (Real.Gamma (x + 1))) =ᶠ[𝓝 ((n : ℝ) + 1)]
      (fun x => Real.log (Real.Gamma x)) + Real.log := by
    filter_upwards [eventually_gt_nhds hpos] with x hx
    exact g_feq hx
  rw [← h1, hev.iteratedDeriv_eq, iteratedDeriv_add (g_cd _ n)
    (Real.contDiffAt_log.mpr hpos.ne'), log_iter hpos]

end Q04

theorem Q04.proof_2 : Q04.stmt_2 := by
  intro m hm
  obtain ⟨n, rfl⟩ : ∃ n, m = n + 1 := ⟨m - 1, by omega⟩
  simp only [Nat.add_sub_cancel]
  push_cast
  clear hm
  refine ⟨?_, ?_⟩
  · rw [← iteratedDeriv_one]
    induction n with
    | zero => simp [Q04.g_iter1, Q04.aL, V05.Hr]
    | succ n ih =>
      push_cast
      rw [Q04.g_step 0 n, ih, V05.Hr_succ]
      simp only [pow_zero, Nat.factorial_zero, Nat.cast_one, mul_one, pow_one, one_div]
      ring
  · intro j hj
    obtain ⟨i, rfl⟩ : ∃ i, j = i + 1 := ⟨j - 1, by omega⟩
    clear hj
    induction n with
    | zero =>
      simp only [Nat.cast_zero, zero_add, Q04.g_iter1, Q04.aL, V05.Hr, Finset.range_zero,
        Finset.sum_empty, sub_zero]
      rw [Nat.factorial_succ (i + 1)]
      push_cast
      field_simp
      ring
    | succ n ih =>
      push_cast
      rw [Q04.g_step (i + 1) n, ih, V05.Hr_succ]
      ring

namespace Q04
open Polynomial

lemma expT_map {A B : Type*} [CommRing A] [Algebra ℚ A] [CommRing B] [Algebra ℚ B]
    (f : A →+* B) (K : ℕ) (L : A) : f (expT K L) = expT K (f L) := by
  simp [expT, map_sum, map_rat_smul, map_pow]

lemma logT_map {A B : Type*} [CommRing A] [Algebra ℚ A] [CommRing B] [Algebra ℚ B]
    (f : A →+* B) (K : ℕ) (L : A) : f (logT K L) = logT K (f L) := by
  simp [logT, map_sum, map_rat_smul, map_pow]

/-- `Eser` with `n` replaced by the polynomial variable. -/
noncomputable def EserP (K : ℕ) : Polynomial (MvPolynomial (Fin 3) (Polynomial ℚ)) :=
  let yv : MvPolynomial (Fin 3) (Polynomial ℚ) := MvPolynomial.X 0
  let av : MvPolynomial (Fin 3) (Polynomial ℚ) := MvPolynomial.X 1
  let bv : MvPolynomial (Fin 3) (Polynomial ℚ) := MvPolynomial.X 2
  let nv : MvPolynomial (Fin 3) (Polynomial ℚ) := MvPolynomial.C Polynomial.X
  let ℓ := bv - av
  let A := (-(nv - 1)) * av + nv * bv
  let w := Polynomial.C yv *
    ∑ r ∈ Finset.Icc 1 K, ((r.factorial : ℚ)⁻¹) • (Polynomial.C ℓ * Polynomial.X) ^ r
  expT K (Polynomial.C A * Polynomial.X + Polynomial.C (1 - 2 * nv) * logT K w)

noncomputable def evN (n : ℕ) :
    Polynomial (MvPolynomial (Fin 3) (Polynomial ℚ)) →+* Polynomial (MvPolynomial (Fin 3) ℚ) :=
  Polynomial.mapRingHom (MvPolynomial.map (Polynomial.evalRingHom (n : ℚ)))

lemma evN_EserP (n K : ℕ) : evN n (EserP K) = Eser n K := by
  unfold EserP Eser
  dsimp only
  refine (expT_map (A := Polynomial (MvPolynomial (Fin 3) (Polynomial ℚ)))
    (B := Polynomial (MvPolynomial (Fin 3) ℚ)) (evN n) K _).trans ?_
  congr 1
  rw [map_add, map_mul]
  refine congrArg₂ _ ?_ ?_
  · simp [evN, Algebra.smul_def]
  · rw [map_mul, logT_map (A := Polynomial (MvPolynomial (Fin 3) (Polynomial ℚ)))
      (B := Polynomial (MvPolynomial (Fin 3) ℚ))]
    simp [evN, Algebra.smul_def, Polynomial.map_sum, Polynomial.C_ofNat]
    left; left; rfl

lemma eF_poly (k j i m : ℕ) : ∃ Q : Polynomial ℚ, ∀ n : ℕ, eF k n j i m = Q.eval (n : ℚ) := by
  refine ⟨((EserP k).coeff k).coeff
    (Finsupp.single 0 j + Finsupp.single 1 i + Finsupp.single 2 m), fun n => ?_⟩
  rw [eF, Ecoef, ← evN_EserP]
  simp [evN, Polynomial.coeff_map, MvPolynomial.coeff_map]

/-- Monomials `y^{s_0} a^{s_1} b^{s_2}` with `s_1 + s_2 = d`, `s_0 ≤ e`. -/
def Wt (d e : ℕ) (P : MvPolynomial (Fin 3) ℚ) : Prop :=
  ∀ s : Fin 3 →₀ ℕ, P.coeff s ≠ 0 → s 1 + s 2 = d ∧ s 0 ≤ e

lemma Wt_zero (d e : ℕ) : Wt d e 0 := by intro s h; simp at h

lemma Wt_add {d e : ℕ} {P Q : MvPolynomial (Fin 3) ℚ} (hP : Wt d e P) (hQ : Wt d e Q) :
    Wt d e (P + Q) := by
  intro s h
  rw [MvPolynomial.coeff_add] at h
  by_cases h1 : P.coeff s = 0
  · rw [h1, zero_add] at h; exact hQ s h
  · exact hP s h1

lemma Wt_smul {d e : ℕ} (c : ℚ) {P : MvPolynomial (Fin 3) ℚ} (hP : Wt d e P) :
    Wt d e (c • P) := by
  intro s h
  rw [MvPolynomial.coeff_smul] at h
  exact hP s (right_ne_zero_of_mul h)

lemma Wt_sum {ι : Type*} (S : Finset ι) {d e : ℕ} (f : ι → MvPolynomial (Fin 3) ℚ)
    (h : ∀ i ∈ S, Wt d e (f i)) : Wt d e (∑ i ∈ S, f i) := by
  classical
  induction S using Finset.induction_on with
  | empty => simpa using Wt_zero d e
  | insert a S ha ih =>
    rw [Finset.sum_insert ha]
    exact Wt_add (h a (Finset.mem_insert_self _ _))
      (ih fun i hi => h i (Finset.mem_insert_of_mem hi))

lemma Wt_mul {d e d' e' : ℕ} {P Q : MvPolynomial (Fin 3) ℚ} (hP : Wt d e P) (hQ : Wt d' e' Q) :
    Wt (d + d') (e + e') (P * Q) := by
  intro s h
  rw [MvPolynomial.coeff_mul] at h
  obtain ⟨⟨a, b⟩, hab, hne⟩ := Finset.exists_ne_zero_of_sum_ne_zero h
  rw [Finset.HasAntidiagonal.mem_antidiagonal] at hab
  obtain ⟨h1, h2⟩ := hP a (left_ne_zero_of_mul hne)
  obtain ⟨h3, h4⟩ := hQ b (right_ne_zero_of_mul hne)
  subst hab
  simp only [Finsupp.add_apply]
  constructor <;> omega

lemma Wt_mono {d e e' : ℕ} {P : MvPolynomial (Fin 3) ℚ} (hP : Wt d e P) (he : e ≤ e') :
    Wt d e' P := fun s h => ⟨(hP s h).1, (hP s h).2.trans he⟩

lemma Wt_one : Wt 0 0 (1 : MvPolynomial (Fin 3) ℚ) := by
  intro s h
  rw [MvPolynomial.coeff_one] at h
  split_ifs at h with h0
  · subst h0; simp
  · exact absurd rfl h

lemma Wt_X (i : Fin 3) : Wt (if i = 0 then 0 else 1) (if i = 0 then 1 else 0)
    (MvPolynomial.X i : MvPolynomial (Fin 3) ℚ) := by
  intro s h
  rw [MvPolynomial.coeff_X'] at h
  split_ifs at h with h0
  · subst h0; fin_cases i <;> simp
  · exact absurd rfl h

lemma Wt_pow {d e : ℕ} {P : MvPolynomial (Fin 3) ℚ} (hP : Wt d e P) (r : ℕ) :
    Wt (r * d) (r * e) (P ^ r) := by
  induction r with
  | zero => simpa using Wt_one
  | succ r ih =>
    rw [pow_succ, show (r + 1) * d = r * d + d by ring, show (r + 1) * e = r * e + e by ring]
    exact Wt_mul ih hP

/-- Every `u^k` coefficient has `(a,b)`-degree `k` and `y`-degree `≤ k`. -/
def GoodP (F : Polynomial (MvPolynomial (Fin 3) ℚ)) : Prop := ∀ k, Wt k k (F.coeff k)

lemma GoodP_add {F G : Polynomial (MvPolynomial (Fin 3) ℚ)} (hF : GoodP F) (hG : GoodP G) :
    GoodP (F + G) := fun k => by rw [Polynomial.coeff_add]; exact Wt_add (hF k) (hG k)

lemma GoodP_smul (c : ℚ) {F : Polynomial (MvPolynomial (Fin 3) ℚ)} (hF : GoodP F) :
    GoodP (c • F) := fun k => by
  rw [Polynomial.coeff_smul]; exact Wt_smul c (hF k)

lemma GoodP_mul {F G : Polynomial (MvPolynomial (Fin 3) ℚ)} (hF : GoodP F) (hG : GoodP G) :
    GoodP (F * G) := fun k => by
  rw [Polynomial.coeff_mul]
  refine Wt_sum _ _ fun x hx => ?_
  rw [Finset.HasAntidiagonal.mem_antidiagonal] at hx
  have := Wt_mul (hF x.1) (hG x.2)
  rwa [hx] at this

lemma GoodP_one : GoodP (1 : Polynomial (MvPolynomial (Fin 3) ℚ)) := fun k => by
  rw [Polynomial.coeff_one]
  split_ifs with h
  · subst h; exact Wt_one
  · exact Wt_zero _ _

lemma GoodP_pow {F : Polynomial (MvPolynomial (Fin 3) ℚ)} (hF : GoodP F) (r : ℕ) :
    GoodP (F ^ r) := by
  induction r with
  | zero => simpa using GoodP_one
  | succ r ih => rw [pow_succ]; exact GoodP_mul ih hF

lemma GoodP_sum {ι : Type*} (S : Finset ι) (f : ι → Polynomial (MvPolynomial (Fin 3) ℚ))
    (h : ∀ i ∈ S, GoodP (f i)) : GoodP (∑ i ∈ S, f i) := fun k => by
  rw [Polynomial.finsetSum_coeff]
  exact Wt_sum _ _ fun i hi => h i hi k

lemma GoodP_CX {P : MvPolynomial (Fin 3) ℚ} (hP : Wt 1 1 P) :
    GoodP (Polynomial.C P * Polynomial.X) := fun k => by
  rw [Polynomial.coeff_C_mul_X]
  split_ifs with h
  · subst h; exact hP
  · exact Wt_zero _ _

lemma GoodP_expT {F : Polynomial (MvPolynomial (Fin 3) ℚ)} (hF : GoodP F) (K : ℕ) :
    GoodP (expT K F) :=
  GoodP_sum _ _ fun r _ => GoodP_smul _ (GoodP_pow hF r)

lemma GoodP_logT {F : Polynomial (MvPolynomial (Fin 3) ℚ)} (hF : GoodP F) (K : ℕ) :
    GoodP (logT K F) :=
  GoodP_sum _ _ fun r _ => GoodP_smul _ (GoodP_pow hF r)

lemma Wt_ell : Wt 1 0 (MvPolynomial.X 2 - MvPolynomial.X 1 : MvPolynomial (Fin 3) ℚ) := by
  rw [sub_eq_add_neg, ← neg_one_smul ℚ (MvPolynomial.X 1)]
  exact Wt_add (by simpa using Wt_X 2) (Wt_smul _ (by simpa using Wt_X 1))

lemma GoodP_w (K : ℕ) : GoodP (Polynomial.C (MvPolynomial.X 0 : MvPolynomial (Fin 3) ℚ) *
    ∑ r ∈ Finset.Icc 1 K, ((r.factorial : ℚ)⁻¹) •
      (Polynomial.C (MvPolynomial.X 2 - MvPolynomial.X 1) * Polynomial.X) ^ r) := by
  rw [Finset.mul_sum]
  refine GoodP_sum _ _ fun r hr => ?_
  rw [mul_smul_comm]
  refine GoodP_smul _ fun k => ?_
  rw [mul_pow, ← Polynomial.C_pow, ← mul_assoc, ← Polynomial.C_mul, Polynomial.coeff_C_mul_X_pow]
  split_ifs with h
  · subst h
    have hr1 : 1 ≤ k := (Finset.mem_Icc.mp hr).1
    have := Wt_mul (by simpa using Wt_X 0) (Wt_pow Wt_ell k)
    simp only [zero_add, mul_one, mul_zero] at this
    exact Wt_mono this hr1
  · exact Wt_zero _ _

lemma GoodP_Eser (n K : ℕ) : GoodP (Eser n K) := by
  unfold Eser
  dsimp only
  refine GoodP_expT (GoodP_add (GoodP_CX ?_) (GoodP_smul _ (GoodP_logT (GoodP_w K) K))) K
  exact Wt_add (Wt_smul _ (Wt_mono (by simpa using Wt_X 1) zero_le_one))
    (Wt_smul _ (Wt_mono (by simpa using Wt_X 2) zero_le_one))

lemma eF_supp (k n j i m : ℕ) (h : eF k n j i m ≠ 0) : j ≤ k ∧ i + m = k := by
  have := GoodP_Eser n k k _ h
  simp at this
  omega

end Q04


namespace Q04
open Asymptotics Filter Topology

lemma poly_littleO_zero : ∀ (K : ℕ) (d : ℕ → ℝ),
    (fun u : ℝ => ∑ k ∈ Finset.range (K + 1), d k * u ^ k) =o[𝓝[≠] 0] (fun u => u ^ K) →
    ∀ k ≤ K, d k = 0 := by
  intro K
  induction K with
  | zero =>
    intro d h k hk
    have hk0 : k = 0 := by omega
    subst hk0
    simp only [zero_add, Finset.range_one, Finset.sum_singleton, pow_zero, mul_one] at h
    have t := (isLittleO_one_iff ℝ).mp h
    exact tendsto_nhds_unique tendsto_const_nhds t
  | succ K ih =>
    intro d h
    have hd0 : d 0 = 0 := by
      have t1 : Tendsto (fun u : ℝ => ∑ k ∈ Finset.range (K + 1 + 1), d k * u ^ k) (𝓝[≠] 0)
          (𝓝 0) := by
        refine h.trans_tendsto ?_
        have : Tendsto (fun u : ℝ => u ^ (K + 1)) (𝓝 0) (𝓝 0) := by
          have := ((continuous_pow (K + 1)).tendsto (0 : ℝ)); simpa using this
        exact this.mono_left nhdsWithin_le_nhds
      have t2 : Tendsto (fun u : ℝ => ∑ k ∈ Finset.range (K + 1 + 1), d k * u ^ k) (𝓝[≠] 0)
          (𝓝 (d 0)) := by
        have hc : Continuous (fun u : ℝ => ∑ k ∈ Finset.range (K + 1 + 1), d k * u ^ k) := by
          fun_prop
        have := hc.tendsto 0
        rw [Finset.sum_eq_single 0 (fun b _ hb => by simp [zero_pow hb]) (by simp)] at this
        simpa using this.mono_left nhdsWithin_le_nhds
      exact tendsto_nhds_unique t2 t1
    have h2 : (fun u : ℝ => ∑ k ∈ Finset.range (K + 1), d (k + 1) * u ^ k) =o[𝓝[≠] 0]
        (fun u => u ^ K) := by
      have hm := (isBigO_refl (fun u : ℝ => u⁻¹) (𝓝[≠] 0)).mul_isLittleO h
      refine (hm.congr' ?_ ?_)
      · filter_upwards [self_mem_nhdsWithin] with u hu
        rw [Finset.sum_range_succ', hd0, zero_mul, add_zero, Finset.mul_sum]
        refine Finset.sum_congr rfl fun k _ => ?_
        have hu' : u ≠ 0 := hu
        rw [pow_succ]; field_simp
      · filter_upwards [self_mem_nhdsWithin] with u hu
        have hu' : u ≠ 0 := hu
        rw [pow_succ]; field_simp
    intro k hk
    rcases k with _ | k
    · exact hd0
    · exact ih (fun k => d (k + 1)) h2 k (by omega)

lemma powdiff {l : Filter ℝ} {a b g : ℝ → ℝ} (h : (fun u => a u - b u) =O[l] g)
    (ha : a =O[l] (fun _ => (1 : ℝ))) (hb : b =O[l] (fun _ => (1 : ℝ))) (m : ℕ) :
    (fun u => a u ^ m - b u ^ m) =O[l] g := by
  induction m with
  | zero => simpa using isBigO_zero g l
  | succ m ih =>
    have e : (fun u => a u ^ (m + 1) - b u ^ (m + 1)) =
        fun u => a u * (a u ^ m - b u ^ m) + b u ^ m * (a u - b u) := by
      ext u; ring
    rw [e]
    have h1 := ha.mul ih
    have h2 := (hb.pow m).mul h
    simp only [one_pow, one_mul] at h1 h2
    exact h1.add h2

lemma expT_real (K : ℕ) (z : ℝ) :
    expT K z = ∑ r ∈ Finset.range (K + 1), z ^ r / (r.factorial : ℝ) := by
  unfold expT
  refine Finset.sum_congr rfl fun r _ => ?_
  rw [Rat.smul_def]; push_cast; ring

lemma logT_real (K : ℕ) (w : ℝ) :
    logT K w = ∑ i ∈ Finset.range K, (-1) ^ i * w ^ (i + 1) / ((i : ℝ) + 1) := by
  unfold logT
  induction K with
  | zero => simp
  | succ K ih =>
    rw [Finset.sum_Icc_succ_top (by omega), ih, Finset.sum_range_succ, Rat.smul_def,
      Nat.add_sub_cancel]
    push_cast; ring

lemma expT_diff {l : Filter ℝ} {a b g : ℝ → ℝ} (h : (fun u => a u - b u) =O[l] g)
    (ha : a =O[l] (fun _ => (1 : ℝ))) (hb : b =O[l] (fun _ => (1 : ℝ))) (K : ℕ) :
    (fun u => expT K (a u) - expT K (b u)) =O[l] g := by
  simp only [expT_real, ← Finset.sum_sub_distrib]
  have := IsBigO.sum (s := Finset.range (K + 1))
    (A := fun r u => a u ^ r / (r.factorial : ℝ) - b u ^ r / (r.factorial : ℝ)) fun r _ =>
      ((powdiff h ha hb r).const_mul_left ((r.factorial : ℝ))⁻¹).congr_left fun u => by ring
  refine this.congr_left fun u => ?_
  simp [Finset.sum_apply]

lemma logT_diff {l : Filter ℝ} {a b g : ℝ → ℝ} (h : (fun u => a u - b u) =O[l] g)
    (ha : a =O[l] (fun _ => (1 : ℝ))) (hb : b =O[l] (fun _ => (1 : ℝ))) (K : ℕ) :
    (fun u => logT K (a u) - logT K (b u)) =O[l] g := by
  simp only [logT_real, ← Finset.sum_sub_distrib]
  have := IsBigO.sum (s := Finset.range K)
    (A := fun (i : ℕ) (u : ℝ) => (-1 : ℝ) ^ i * a u ^ (i + 1) / ((i : ℝ) + 1) -
      (-1) ^ i * b u ^ (i + 1) / ((i : ℝ) + 1)) fun i _ =>
      ((powdiff h ha hb (i + 1)).const_mul_left ((-1) ^ i / ((i : ℝ) + 1))).congr_left
        fun u => by ring
  refine this.congr_left fun u => ?_
  simp [Finset.sum_apply]

lemma exp_sub_expT (K : ℕ) :
    (fun z : ℝ => Real.exp z - expT K z) =O[𝓝 0] (fun z => z ^ (K + 1)) := by
  refine IsBigO.of_bound (((K + 1).succ : ℝ) / ((K + 1).factorial * (K + 1))) ?_
  filter_upwards [Metric.ball_mem_nhds (0 : ℝ) one_pos] with z hz
  rw [Metric.mem_ball, dist_zero_right] at hz
  have := Real.exp_bound (x := z) (n := K + 1) (Real.norm_eq_abs z ▸ hz.le) (by omega)
  rw [expT_real, Real.norm_eq_abs, Real.norm_eq_abs, abs_pow]
  convert this using 1
  push_cast
  ring

lemma log_sub_logT (K : ℕ) :
    (fun w : ℝ => Real.log (1 + w) - logT K w) =O[𝓝 0] (fun w => w ^ (K + 1)) := by
  refine IsBigO.of_bound 2 ?_
  filter_upwards [Metric.ball_mem_nhds (0 : ℝ) (by norm_num : (0 : ℝ) < 1 / 2)] with w hw
  rw [Metric.mem_ball, dist_zero_right, Real.norm_eq_abs] at hw
  have hw' : |-w| < 1 := by rw [abs_neg]; linarith
  have := Real.abs_log_sub_add_sum_range_le hw' K
  rw [sub_neg_eq_add] at this
  rw [logT_real, Real.norm_eq_abs, Real.norm_eq_abs, abs_pow]
  have e : Real.log (1 + w) - ∑ i ∈ Finset.range K, (-1) ^ i * w ^ (i + 1) / ((i : ℝ) + 1) =
      (∑ i ∈ Finset.range K, (-w) ^ (i + 1) / ((i : ℝ) + 1)) + Real.log (1 + w) := by
    rw [sub_eq_neg_add, ← Finset.sum_neg_distrib]
    congr 1
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [neg_pow w]; ring
  rw [e]
  refine this.trans ?_
  rw [abs_neg]
  have h1 : 1 / 2 ≤ 1 - |w| := by linarith
  have h0 : 0 ≤ |w| ^ (K + 1) := by positivity
  rw [div_le_iff₀ (by linarith)]
  nlinarith

lemma log_isBigO :
    (fun w : ℝ => Real.log (1 + w)) =O[𝓝 0] (fun w => w) := by
  have h0 := log_sub_logT 0
  simpa [logT_real] using h0

lemma sum_Icc_one {M : Type*} [AddCommMonoid M] (g : ℕ → M) (K : ℕ) :
    ∑ x ∈ Finset.Icc 1 K, g x = ∑ k ∈ Finset.range K, g (k + 1) := by
  induction K with
  | zero => simp
  | succ K ih => rw [Finset.sum_Icc_succ_top (by omega), ih, Finset.sum_range_succ]

/-- Evaluation `y ↦ y, a ↦ log y, b ↦ log(1-y)` on coefficients. -/
noncomputable def phiY (y : ℝ) : MvPolynomial (Fin 3) ℚ →ₐ[ℚ] ℝ :=
  MvPolynomial.aeval (![y, Real.log y, Real.log (1 - y)] : Fin 3 → ℝ)

lemma ev_Eser (n K : ℕ) (y u : ℝ) :
    Polynomial.eval₂ (phiY y).toRingHom u (Eser n K) =
      expT K (u * (-((n : ℝ) - 1) * Real.log y + (n : ℝ) * Real.log (1 - y)) +
        (1 - 2 * (n : ℝ)) * logT K (y * (expT K (u * (Real.log (1 - y) - Real.log y)) - 1))) := by
  have hm := expT_map (A := Polynomial (MvPolynomial (Fin 3) ℚ)) (B := ℝ)
    (Polynomial.eval₂RingHom (phiY y).toRingHom u)
  have hl := logT_map (A := Polynomial (MvPolynomial (Fin 3) ℚ)) (B := ℝ)
    (Polynomial.eval₂RingHom (phiY y).toRingHom u)
  simp only [Polynomial.coe_eval₂RingHom] at hm hl
  unfold Eser
  dsimp only
  rw [hm]
  congr 1
  rw [Polynomial.eval₂_add, Polynomial.eval₂_mul, Polynomial.eval₂_C, Polynomial.eval₂_X]
  have hs := map_rat_smul (Polynomial.eval₂RingHom (phiY y).toRingHom u)
  simp only [Polynomial.coe_eval₂RingHom] at hs
  rw [hs, hl]
  congr 1
  · simp [phiY, Algebra.smul_def]; ring
  · rw [Rat.smul_def]
    push_cast
    congr 2
    rw [Polynomial.eval₂_mul, Polynomial.eval₂_finsetSum]
    simp only [hs]
    simp [phiY, expT_real, Finset.sum_range_succ', mul_pow, Polynomial.eval₂_pow]
    left
    rw [sum_Icc_one]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [Rat.smul_def]; push_cast; ring

lemma mulpow (c : ℝ) (m : ℕ) : (fun u : ℝ => (u * c) ^ m) =O[𝓝 0] (fun u => u ^ m) :=
  IsBigO.of_bound (|c| ^ m) (Eventually.of_forall fun u => by
    rw [Real.norm_eq_abs, Real.norm_eq_abs, mul_pow, abs_mul, abs_pow, abs_pow, mul_comm])

lemma chain (A ℓ s y : ℝ) (K : ℕ) :
    (fun u : ℝ => Real.exp (u * A + s * Real.log (1 + y * (Real.exp (u * ℓ) - 1))) -
      expT K (u * A + s * logT K (y * (expT K (u * ℓ) - 1)))) =O[𝓝 0] (fun u => u ^ (K + 1)) := by
  have hlin : ∀ c : ℝ, Tendsto (fun u : ℝ => u * c) (𝓝 0) (𝓝 0) := fun c => by
    simpa using (tendsto_id (x := 𝓝 (0 : ℝ))).mul_const c
  have hu1 : (fun u : ℝ => u ^ (K + 1)) =O[𝓝 0] (fun u => u) := by
    refine IsBigO.of_bound 1 ?_
    filter_upwards [Metric.ball_mem_nhds (0 : ℝ) one_pos] with u hu
    rw [Metric.mem_ball, dist_zero_right] at hu
    rw [one_mul, norm_pow]
    exact pow_le_of_le_one (norm_nonneg _) hu.le (by omega)
  have huO1 : (fun u : ℝ => u) =O[𝓝 0] (fun _ => (1 : ℝ)) :=
    (tendsto_id (x := 𝓝 (0 : ℝ))).isBigO_one ℝ
  have hcA : (fun u : ℝ => u * A) =O[𝓝 0] (fun u => u) := by
    simpa using mulpow A 1
  have E1 : (fun u : ℝ => Real.exp (u * ℓ) - expT K (u * ℓ)) =O[𝓝 0] (fun u => u ^ (K + 1)) :=
    ((exp_sub_expT K).comp_tendsto (hlin ℓ)).trans (mulpow ℓ (K + 1))
  have E0 : (fun u : ℝ => Real.exp (u * ℓ) - 1) =O[𝓝 0] (fun u => u) := by
    have := ((exp_sub_expT 0).comp_tendsto (hlin ℓ)).trans (mulpow ℓ 1)
    simp only [expT_real, Function.comp_def] at this
    simpa using this
  have hwW : (fun u : ℝ => y * (Real.exp (u * ℓ) - 1) - y * (expT K (u * ℓ) - 1)) =O[𝓝 0]
      (fun u => u ^ (K + 1)) :=
    (E1.const_mul_left y).congr_left fun u => by ring
  have hw : (fun u : ℝ => y * (Real.exp (u * ℓ) - 1)) =O[𝓝 0] (fun u => u) := E0.const_mul_left y
  have hW : (fun u : ℝ => y * (expT K (u * ℓ) - 1)) =O[𝓝 0] (fun u => u) :=
    (hw.sub (hwW.trans hu1)).congr_left fun u => by ring
  have hwt : Tendsto (fun u : ℝ => y * (Real.exp (u * ℓ) - 1)) (𝓝 0) (𝓝 0) :=
    hw.trans_tendsto tendsto_id
  have hlog1 : (fun u : ℝ => Real.log (1 + y * (Real.exp (u * ℓ) - 1)) -
      logT K (y * (Real.exp (u * ℓ) - 1))) =O[𝓝 0] (fun u => u ^ (K + 1)) :=
    ((log_sub_logT K).comp_tendsto hwt).trans (hw.pow (K + 1))
  have hlog2 := logT_diff hwW (hw.trans huO1) (hW.trans huO1) K
  have hL : (fun u : ℝ => (u * A + s * Real.log (1 + y * (Real.exp (u * ℓ) - 1))) -
      (u * A + s * logT K (y * (expT K (u * ℓ) - 1)))) =O[𝓝 0] (fun u => u ^ (K + 1)) :=
    ((hlog1.add hlog2).const_mul_left s).congr_left fun u => by ring
  have hLu : (fun u : ℝ => u * A + s * Real.log (1 + y * (Real.exp (u * ℓ) - 1))) =O[𝓝 0]
      (fun u => u) :=
    hcA.add (((log_isBigO.comp_tendsto hwt).trans hw).const_mul_left s)
  have hLp : (fun u : ℝ => u * A + s * logT K (y * (expT K (u * ℓ) - 1))) =O[𝓝 0]
      (fun u => u) :=
    (hLu.sub (hL.trans hu1)).congr_left fun u => by ring
  have hLt := hLu.trans_tendsto tendsto_id
  have hE1 := ((exp_sub_expT K).comp_tendsto hLt).trans (hLu.pow (K + 1))
  have hE2 := expT_diff hL (hLu.trans huO1) (hLp.trans huO1) K
  exact (hE1.add hE2).congr_left fun u => by simp only [Function.comp]; ring

lemma E_eq (n : ℕ) {y : ℝ} (hy0 : 0 < y) (hy1 : y < 1) (u : ℝ) :
    V05.E n u y = Real.exp (u * (-((n : ℝ) - 1) * Real.log y + (n : ℝ) * Real.log (1 - y)) +
      (1 - 2 * (n : ℝ)) * Real.log (1 + y * (Real.exp (u * (Real.log (1 - y) - Real.log y)) - 1))) := by
  have hf : 0 < V05.f u y := by
    unfold V05.f; have := Real.exp_pos (u * (Real.log (1 - y) - Real.log y)); nlinarith
  unfold V05.E
  rw [Real.rpow_def_of_pos hf, ← Real.exp_add]
  unfold V05.f
  congr 1
  rw [show 1 - y + y * Real.exp (u * (Real.log (1 - y) - Real.log y)) =
    1 + y * (Real.exp (u * (Real.log (1 - y) - Real.log y)) - 1) by ring]
  ring

lemma E_contDiff (n : ℕ) {y : ℝ} (hy0 : 0 < y) (hy1 : y < 1) (k : ℕ) :
    ContDiff ℝ k (fun u => V05.E n u y) := by
  have e : (fun u => V05.E n u y) = fun u => Real.exp (u * (-((n : ℝ) - 1) * Real.log y +
      (n : ℝ) * Real.log (1 - y)) + (1 - 2 * (n : ℝ)) *
        Real.log (1 + y * (Real.exp (u * (Real.log (1 - y) - Real.log y)) - 1))) :=
    funext (E_eq n hy0 hy1)
  rw [e]
  refine Real.contDiff_exp.comp (ContDiff.add (by fun_prop) (ContDiff.mul contDiff_const ?_))
  refine ContDiff.log (by fun_prop) fun u => ?_
  have := Real.exp_pos (u * (Real.log (1 - y) - Real.log y)); nlinarith

lemma tail_isBigO (P : Polynomial (MvPolynomial (Fin 3) ℚ)) (φ : MvPolynomial (Fin 3) ℚ →+* ℝ)
    (K : ℕ) :
    (fun u : ℝ => Polynomial.eval₂ φ u P -
      ∑ k ∈ Finset.range (K + 1), φ (P.coeff k) * u ^ k) =O[𝓝 0] (fun u => u ^ (K + 1)) := by
  set N := P.natDegree + K + 2
  have hN : P.natDegree < N := by omega
  have hKN : K + 1 ≤ N := by omega
  have hg : Continuous (fun u : ℝ => ∑ k ∈ Finset.Ico (K + 1) N, φ (P.coeff k) * u ^ (k - (K + 1))) := by
    fun_prop
  have hb := (hg.tendsto 0).isBigO_one ℝ
  have hm := (isBigO_refl (fun u : ℝ => u ^ (K + 1)) (𝓝 0)).mul hb
  simp only [mul_one] at hm
  refine hm.congr_left fun u => ?_
  rw [Polynomial.eval₂_eq_sum_range' φ hN, ← Finset.sum_range_add_sum_Ico _ hKN, add_sub_cancel_left,
    Finset.mul_sum]
  refine Finset.sum_congr rfl fun k hk => ?_
  have := (Finset.mem_Ico.mp hk).1
  rw [mul_left_comm, ← pow_add, Nat.add_sub_cancel' this]

lemma taylor_coeff_E (n k : ℕ) {y : ℝ} (hy0 : 0 < y) (hy1 : y < 1) :
    iteratedDeriv k (fun u => V05.E n u y) 0 / (k.factorial : ℝ) = phiY y (Ecoef k n) := by
  set f := fun u => V05.E n u y with hfdef
  have hT := taylor_isLittleO_univ (x₀ := 0) (E_contDiff n hy0 hy1 k)
  simp only [taylor_within_apply, iteratedDerivWithin_univ, sub_zero] at hT
  have hA := chain (-((n : ℝ) - 1) * Real.log y + (n : ℝ) * Real.log (1 - y))
    (Real.log (1 - y) - Real.log y) (1 - 2 * (n : ℝ)) y k
  have hB := tail_isBigO (Eser n k) (phiY y).toRingHom k
  have hC : (fun u => f u - ∑ j ∈ Finset.range (k + 1), (phiY y) ((Eser n k).coeff j) * u ^ j)
      =o[𝓝 0] (fun u : ℝ => u ^ k) := by
    refine ((hA.add hB).trans_isLittleO (isLittleO_pow_pow (Nat.lt_succ_self k))).congr_left ?_
    intro u
    simp only [hfdef, E_eq n hy0 hy1 u, ev_Eser]
    simp
  have hD := (hT.sub hC).mono (nhdsWithin_le_nhds (s := {(0 : ℝ)}ᶜ))
  have := poly_littleO_zero k (fun j => (phiY y) ((Eser n k).coeff j) -
    iteratedDeriv j f 0 / (j.factorial : ℝ)) (hD.congr_left fun u => by
      simp only [smul_eq_mul, ← Finset.sum_sub_distrib]
      rw [sub_sub_sub_cancel_left, ← Finset.sum_sub_distrib]
      refine Finset.sum_congr rfl fun j _ => by ring) k le_rfl
  rw [Ecoef]
  linarith

lemma phiY_Ecoef (k n : ℕ) (y : ℝ) :
    phiY y (Ecoef k n) = ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
      (eF k n j i (k - i) : ℝ) * y ^ j * Real.log y ^ i * Real.log (1 - y) ^ (k - i) := by
  set ι : ℕ × ℕ → (Fin 3 →₀ ℕ) := fun p =>
    Finsupp.single 0 p.1 + Finsupp.single 1 p.2 + Finsupp.single 2 (k - p.2) with hι
  set F : (Fin 3 →₀ ℕ) → ℝ := fun s => ((Ecoef k n).coeff s : ℝ) * y ^ s 0 *
    Real.log y ^ s 1 * Real.log (1 - y) ^ s 2 with hF
  have hinj : Set.InjOn ι ↑(Finset.range (k + 1) ×ˢ Finset.range (k + 1)) := by
    intro p _ q _ h
    have h0 := congrArg (fun s => s 0) h
    have h1 := congrArg (fun s => s 1) h
    simp [hι] at h0 h1
    exact Prod.ext h0 h1
  have hR : ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
      (eF k n j i (k - i) : ℝ) * y ^ j * Real.log y ^ i * Real.log (1 - y) ^ (k - i) =
      ∑ s ∈ (Finset.range (k + 1) ×ˢ Finset.range (k + 1)).image ι, F s := by
    rw [Finset.sum_image hinj, Finset.sum_product]
    refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun i _ => ?_
    simp [hF, hι, eF]
  rw [hR, phiY, MvPolynomial.aeval_eq_eval₂Hom, MvPolynomial.coe_eval₂Hom,
    MvPolynomial.eval₂_eq']
  simp only [Fin.prod_univ_three, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.head_cons, Matrix.tail_cons, eq_ratCast]
  refine (Finset.sum_subset (s₂ := (Finset.range (k + 1) ×ˢ Finset.range (k + 1)).image ι)
    ?_ ?_).trans ?_
  · intro s hs
    have hs' := GoodP_Eser n k k s (MvPolynomial.mem_support_iff.mp hs)
    rw [Finset.mem_image]
    refine ⟨(s 0, s 1), ?_, ?_⟩
    · simp only [Finset.mem_product, Finset.mem_range]; omega
    · ext t
      fin_cases t <;> simp [hι] <;> omega
  · intro s _ hs
    rw [MvPolynomial.notMem_support_iff.mp hs]; simp
  · refine Finset.sum_congr rfl fun s _ => ?_
    simp [hF]; ring

theorem stmt3_c (k n : ℕ) {y : ℝ} (hy0 : 0 < y) (hy1 : y < 1) :
    iteratedDeriv k (fun u => V05.E n u y) 0 / (k.factorial : ℝ) =
      ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
        (eF k n j i (k - i) : ℝ) * y ^ j * Real.log y ^ i * Real.log (1 - y) ^ (k - i) := by
  rw [taylor_coeff_E n k hy0 hy1, phiY_Ecoef]

end Q04

theorem Q04.proof_3 : Q04.stmt_3 := by
  intro k
  refine ⟨Q04.eF_poly k, fun n j i m h => Q04.eF_supp k n j i m h, ?_⟩
  intro n _ y hy0 hy1
  exact Q04.stmt3_c k n hy0 hy1
namespace Q04

/-- weight of a monomial in the zeta variables -/
lemma WeightLE_iff' (d : ℕ) (P : MvPolynomial ℕ ℚ) :
    WeightLE d P ↔ ∀ s, P.coeff s ≠ 0 → (∀ w ∈ s.support, 2 ≤ w) ∧ (s.sum fun w e => e * w) ≤ d :=
  ⟨fun h s hs => h s (MvPolynomial.mem_support_iff.mpr hs),
    fun h s hs => h s (MvPolynomial.mem_support_iff.mp hs)⟩

lemma WL_zero (d : ℕ) : WeightLE d 0 := by intro s h; simp at h

lemma WL_add {d : ℕ} {P Q : MvPolynomial ℕ ℚ} (hP : WeightLE d P) (hQ : WeightLE d Q) :
    WeightLE d (P + Q) := by
  rw [WeightLE_iff'] at *
  intro s h
  rw [MvPolynomial.coeff_add] at h
  by_cases h1 : P.coeff s = 0
  · rw [h1, zero_add] at h; exact hQ s h
  · exact hP s h1

lemma WL_smul {d : ℕ} (c : ℚ) {P : MvPolynomial ℕ ℚ} (hP : WeightLE d P) :
    WeightLE d (c • P) := by
  rw [WeightLE_iff'] at *
  intro s h
  rw [MvPolynomial.coeff_smul] at h
  exact hP s (right_ne_zero_of_mul h)

lemma WL_neg {d : ℕ} {P : MvPolynomial ℕ ℚ} (hP : WeightLE d P) : WeightLE d (-P) := by
  rw [← neg_one_smul ℚ P]; exact WL_smul _ hP

lemma WL_mono {d d' : ℕ} {P : MvPolynomial ℕ ℚ} (hP : WeightLE d P) (h : d ≤ d') :
    WeightLE d' P := fun s hs => ⟨(hP s hs).1, (hP s hs).2.trans h⟩

lemma WL_mul {d d' : ℕ} {P Q : MvPolynomial ℕ ℚ} (hP : WeightLE d P) (hQ : WeightLE d' Q) :
    WeightLE (d + d') (P * Q) := by
  rw [WeightLE_iff'] at *
  intro s h
  rw [MvPolynomial.coeff_mul] at h
  obtain ⟨⟨a, b⟩, hab, hne⟩ := Finset.exists_ne_zero_of_sum_ne_zero h
  rw [Finset.HasAntidiagonal.mem_antidiagonal] at hab
  dsimp only at hab hne
  obtain ⟨h1, h2⟩ := hP a (left_ne_zero_of_mul hne)
  obtain ⟨h3, h4⟩ := hQ b (right_ne_zero_of_mul hne)
  subst hab
  refine ⟨fun w hw => ?_, ?_⟩
  · rcases Finset.mem_union.mp (Finsupp.support_add hw) with hw | hw
    · exact h1 w hw
    · exact h3 w hw
  · rw [Finsupp.sum_add_index' (fun _ => by simp) (fun _ _ _ => by ring)]
    omega

lemma WL_C (q : ℚ) : WeightLE 0 (MvPolynomial.C q : MvPolynomial ℕ ℚ) := by
  rw [WeightLE_iff']
  intro s h
  rw [MvPolynomial.coeff_C] at h
  split_ifs at h with h0
  · subst h0; simp
  · exact absurd rfl h

lemma WL_X {w : ℕ} (hw : 2 ≤ w) : WeightLE w (MvPolynomial.X w : MvPolynomial ℕ ℚ) := by
  rw [WeightLE_iff']
  intro s h
  rw [MvPolynomial.coeff_X'] at h
  split_ifs at h with h0
  · subst h0
    refine ⟨fun v hv => ?_, ?_⟩
    · rw [Finsupp.support_single_ne_zero _ one_ne_zero, Finset.mem_singleton] at hv
      omega
    · simp [Finsupp.sum_single_index]
  · exact absurd rfl h

abbrev S2 := MvPolynomial (Fin 2) (MvPolynomial ℕ ℚ)

/-- homogeneous of degree `d` in `x, x'`, coefficients of weight `0` -/
def Hm (d : ℕ) (G : S2) : Prop :=
  ∀ t : Fin 2 →₀ ℕ, WeightLE 0 (G.coeff t) ∧ (G.coeff t ≠ 0 → t 0 + t 1 = d)

lemma Hm_add {d : ℕ} {F G : S2} (hF : Hm d F) (hG : Hm d G) : Hm d (F + G) := fun t => by
  rw [MvPolynomial.coeff_add]
  refine ⟨WL_add (hF t).1 (hG t).1, fun h => ?_⟩
  by_cases h1 : F.coeff t = 0
  · rw [h1, zero_add] at h; exact (hG t).2 h
  · exact (hF t).2 h1

lemma Hm_smul {d : ℕ} (c : ℚ) {F : S2} (hF : Hm d F) : Hm d (c • F) := fun t => by
  rw [MvPolynomial.coeff_smul]
  exact ⟨WL_smul c (hF t).1, fun h => (hF t).2 fun h0 => h (by rw [h0, smul_zero])⟩

lemma Hm_sub {d : ℕ} {F G : S2} (hF : Hm d F) (hG : Hm d G) : Hm d (F - G) := by
  rw [sub_eq_add_neg, ← neg_one_smul ℚ G]; exact Hm_add hF (Hm_smul _ hG)

lemma Hm_mul {d d' : ℕ} {F G : S2} (hF : Hm d F) (hG : Hm d' G) : Hm (d + d') (F * G) :=
  fun t => by
  rw [MvPolynomial.coeff_mul]
  constructor
  · refine Finset.sum_induction _ (WeightLE 0) (fun a b ha hb => WL_add ha hb) (WL_zero 0) ?_
    intro x _
    simpa using WL_mul (hF x.1).1 (hG x.2).1
  · intro h
    obtain ⟨⟨a, b⟩, hab, hne⟩ := Finset.exists_ne_zero_of_sum_ne_zero h
    rw [Finset.HasAntidiagonal.mem_antidiagonal] at hab
    have h1 := (hF a).2 (left_ne_zero_of_mul hne)
    have h2 := (hG b).2 (right_ne_zero_of_mul hne)
    subst hab
    simp only [Finsupp.add_apply]
    omega

lemma Hm_one : Hm 0 (1 : S2) := fun t => by
  rw [MvPolynomial.coeff_one]
  split_ifs with h
  · subst h; exact ⟨by simpa using WL_C 1, fun _ => rfl⟩
  · exact ⟨WL_zero 0, fun h' => absurd rfl h'⟩

lemma Hm_pow {d : ℕ} {F : S2} (hF : Hm d F) (r : ℕ) : Hm (r * d) (F ^ r) := by
  induction r with
  | zero => simpa using Hm_one
  | succ r ih => rw [pow_succ, show (r + 1) * d = r * d + d by ring]; exact Hm_mul ih hF

lemma Hm_X (i : Fin 2) : Hm 1 (MvPolynomial.X i : S2) := fun t => by
  rw [MvPolynomial.coeff_X']
  split_ifs with h
  · subst h; exact ⟨by simpa using WL_C 1, fun _ => by fin_cases i <;> simp⟩
  · exact ⟨WL_zero 0, fun h' => absurd rfl h'⟩

/-- every `x^{t_0} x'^{t_1}` coefficient has zeta weight `≤ t_0 + t_1` -/
def GoodD (F : S2) : Prop := ∀ t : Fin 2 →₀ ℕ, WeightLE (t 0 + t 1) (F.coeff t)

lemma GoodD_add {F G : S2} (hF : GoodD F) (hG : GoodD G) : GoodD (F + G) := fun t => by
  rw [MvPolynomial.coeff_add]; exact WL_add (hF t) (hG t)

lemma GoodD_smul (c : ℚ) {F : S2} (hF : GoodD F) : GoodD (c • F) := fun t => by
  rw [MvPolynomial.coeff_smul]; exact WL_smul c (hF t)

lemma GoodD_mul {F G : S2} (hF : GoodD F) (hG : GoodD G) : GoodD (F * G) := fun t => by
  rw [MvPolynomial.coeff_mul]
  refine Finset.sum_induction _ (WeightLE (t 0 + t 1)) (fun a b ha hb => WL_add ha hb)
    (WL_zero _) ?_
  intro x hx
  rw [Finset.HasAntidiagonal.mem_antidiagonal] at hx
  have := WL_mul (hF x.1) (hG x.2)
  subst hx
  simp only [Finsupp.add_apply] at this ⊢
  exact WL_mono this (by omega)

lemma GoodD_one : GoodD (1 : S2) := fun t => by
  rw [MvPolynomial.coeff_one]
  split_ifs with h
  · subst h; simpa using WL_C 1
  · exact WL_zero _

lemma GoodD_pow {F : S2} (hF : GoodD F) (r : ℕ) : GoodD (F ^ r) := by
  induction r with
  | zero => simpa using GoodD_one
  | succ r ih => rw [pow_succ]; exact GoodD_mul ih hF

lemma GoodD_sum {ι : Type*} (S : Finset ι) (f : ι → S2) (h : ∀ i ∈ S, GoodD (f i)) :
    GoodD (∑ i ∈ S, f i) := fun t => by
  rw [MvPolynomial.coeff_sum]
  exact Finset.sum_induction _ _ (fun a b ha hb => WL_add ha hb) (WL_zero _) fun i hi => h i hi t

lemma GoodD_of_Hm {d : ℕ} {F : S2} (hF : Hm d F) : GoodD F := fun t =>
  WL_mono (hF t).1 (Nat.zero_le _)

lemma GoodD_CX {w : ℕ} (hw : 2 ≤ w) {G : S2} (hG : Hm w G) :
    GoodD (MvPolynomial.C (MvPolynomial.X w) * G) := fun t => by
  rw [MvPolynomial.coeff_C_mul]
  by_cases h : G.coeff t = 0
  · rw [h, mul_zero]; exact WL_zero _
  · rw [(hG t).2 h]
    simpa using WL_mul (WL_X hw) (hG t).1

lemma GoodD_expT {F : S2} (hF : GoodD F) (K : ℕ) : GoodD (expT K F) :=
  GoodD_sum _ _ fun r _ => GoodD_smul _ (GoodD_pow hF r)

lemma GoodD_DeltaT (K p q : ℕ) : GoodD (DeltaT K p q) := by
  unfold DeltaT
  dsimp only
  refine GoodD_add (GoodD_sum _ _ fun r _ => GoodD_smul _ ?_)
    (GoodD_sum _ _ fun w hw => GoodD_smul _ ?_)
  · have hx := Hm_X 0
    have hx' := Hm_X 1
    exact GoodD_of_Hm (Hm_sub (Hm_add (Hm_smul _ (Hm_pow hx r)) (Hm_smul _ (Hm_pow hx' r)))
      (Hm_smul _ (Hm_pow (Hm_add hx hx') r)))
  · have h := Hm_sub (Hm_sub (Hm_pow (Hm_add (Hm_X 0) (Hm_X 1)) w) (Hm_pow (Hm_X 0) w))
      (Hm_pow (Hm_X 1) w)
    simp only [mul_one] at h
    exact GoodD_CX (Finset.mem_Icc.mp hw).1 h

lemma DsharpF_weight (i m p q : ℕ) : WeightLE (i + m) (DsharpF i m p q) := by
  unfold DsharpF
  refine WL_smul _ ?_
  have := GoodD_expT (GoodD_DeltaT (i + m) p q) (i + m) (Finsupp.single 0 i + Finsupp.single 1 m)
  simpa using this

lemma WL_sum {ι : Type*} (S : Finset ι) {d : ℕ} (f : ι → MvPolynomial ℕ ℚ)
    (h : ∀ i ∈ S, WeightLE d (f i)) : WeightLE d (∑ i ∈ S, f i) :=
  Finset.sum_induction _ _ (fun a b ha hb => WL_add ha hb) (WL_zero _) h

lemma Pf_weight (k n : ℕ) : WeightLE k (Pf k n) := by
  unfold Pf
  refine WL_smul _ (WL_sum _ _ fun j _ => WL_sum _ _ fun i hi => WL_smul _ ?_)
  have := DsharpF_weight i (k - i) (n + j) n
  rwa [Nat.add_sub_cancel' (Nat.lt_succ_iff.mp (Finset.mem_range.mp hi))] at this

end Q04

namespace Q04
open Asymptotics Filter Topology MeasureTheory Set

lemma exp_hasSum' (z : ℝ) : HasSum (fun k : ℕ => z ^ k / (k.factorial : ℝ)) (Real.exp z) := by
  rw [Real.exp_eq_exp_ℝ]; exact NormedSpace.expSeries_div_hasSum_exp z

lemma exp_rem (N : ℕ) (z : ℝ) :
    |Real.exp z - expT N z| ≤ |z| ^ (N + 1) * Real.exp |z| := by
  have h2 := (hasSum_nat_add_iff' (N + 1)).mpr (exp_hasSum' z)
  have hb := (exp_hasSum' |z|).mul_left (|z| ^ (N + 1))
  rw [expT_real, ← Real.norm_eq_abs]
  refine h2.norm_le_of_bounded hb fun k => ?_
  rw [Real.norm_eq_abs, abs_div, abs_pow, Nat.abs_cast, pow_add]
  have hf : (k.factorial : ℝ) ≤ ((k + (N + 1)).factorial : ℝ) := by
    exact_mod_cast Nat.factorial_le (by omega)
  have hk : 0 < (k.factorial : ℝ) := by positivity
  calc |z| ^ k * |z| ^ (N + 1) / ((k + (N + 1)).factorial : ℝ)
      ≤ |z| ^ k * |z| ^ (N + 1) / (k.factorial : ℝ) :=
        div_le_div_of_nonneg_left (by positivity) hk hf
    _ = _ := by ring

/-- `λ log y + log(1-y)` -/
noncomputable def Lf (lam y : ℝ) : ℝ := lam * Real.log y + Real.log (1 - y)

/-- the Beta weight `y^{p-1}(1-y)^{q-1}` -/
noncomputable def wB (p q : ℕ) (y : ℝ) : ℝ := y ^ (p - 1) * (1 - y) ^ (q - 1)

lemma mom_fun (p q : ℕ) (lam : ℝ) (k : ℕ) :
    (fun y => wB p q y * Lf lam y ^ k) = fun y => ∑ i ∈ Finset.range (k + 1),
      (lam ^ i * (k.choose i : ℝ)) * V05.ibf p q i (k - i) y := by
  funext y
  rw [Lf, add_pow, Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  unfold wB V05.ibf; ring

lemma ii_sum {ι : Type*} (s : Finset ι) {f : ι → ℝ → ℝ}
    (h : ∀ i ∈ s, IntervalIntegrable (f i) volume 0 1) :
    IntervalIntegrable (fun y => ∑ i ∈ s, f i y) volume 0 1 := by
  have := IntervalIntegrable.sum s h
  convert this using 1; funext y; simp [Finset.sum_apply]

lemma mom_int (p q : ℕ) (lam : ℝ) (k : ℕ) :
    IntervalIntegrable (fun y => wB p q y * Lf lam y ^ k) volume 0 1 := by
  rw [mom_fun]
  exact ii_sum _ fun i _ => (V05.ibf_int _ _ _ _).const_mul _

lemma moment_eq (p q : ℕ) (lam : ℝ) (k : ℕ) :
    ∫ y in (0:ℝ)..1, wB p q y * Lf lam y ^ k =
      ∑ i ∈ Finset.range (k + 1), lam ^ i * (k.choose i : ℝ) * V05.Ib p q i (k - i) := by
  rw [mom_fun, intervalIntegral.integral_finsetSum fun i _ => (V05.ibf_int _ _ _ _).const_mul _]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [intervalIntegral.integral_const_mul]; rfl

lemma ls_eq {y : ℝ} (hy0 : 0 < y) (hy1 : y < 1) : V05.ls y = -Real.log (y * (1 - y)) := by
  unfold V05.ls
  rw [abs_of_neg (Real.log_neg hy0 hy1), abs_of_neg (Real.log_neg (by linarith) (by linarith)),
    Real.log_mul hy0.ne' (by linarith)]
  ring

lemma exp_ls {y c : ℝ} (hy0 : 0 < y) (hy1 : y < 1) :
    Real.exp (c * V05.ls y) = (y * (1 - y)) ^ (-c) := by
  rw [ls_eq hy0 hy1, Real.rpow_def_of_pos (mul_pos hy0 (by linarith))]
  ring_nf

lemma wB_bd (p q : ℕ) {y : ℝ} (hy0 : 0 < y) (hy1 : y < 1) : 0 ≤ wB p q y ∧ wB p q y ≤ 1 := by
  unfold wB
  have h1 : 0 < 1 - y := by linarith
  refine ⟨by positivity, ?_⟩
  exact mul_le_one₀ (pow_le_one₀ hy0.le hy1.le) (by positivity) (pow_le_one₀ h1.le (by linarith))

lemma Lf_bd (lam : ℝ) (y : ℝ) : |Lf lam y| ≤ (|lam| + 1) * V05.ls y := by
  unfold Lf V05.ls
  have := abs_add_le (lam * Real.log y) (Real.log (1 - y))
  rw [abs_mul] at this
  have h1 := abs_nonneg (Real.log y)
  have h2 := abs_nonneg (Real.log (1 - y))
  have h3 := abs_nonneg lam
  nlinarith

lemma int_taylor (p q : ℕ) (lam : ℝ) (N : ℕ) :
    (fun u : ℝ => (∫ y in (0:ℝ)..1, wB p q y * Real.exp (u * Lf lam y)) -
      ∑ k ∈ Finset.range (N + 1),
        (∫ y in (0:ℝ)..1, wB p q y * Lf lam y ^ k) / (k.factorial : ℝ) * u ^ k)
      =O[𝓝 0] (fun u => u ^ (N + 1)) := by
  set δ : ℝ := 1 / (4 * (|lam| + 1)) with hδ
  have hl : 0 < |lam| + 1 := by positivity
  have hδ0 : 0 < δ := by positivity
  set C : ℝ := (|lam| + 1) ^ (N + 1) * (4 * (((N + 1 : ℕ) : ℝ) + 1)) ^ (N + 1) with hC
  have hC0 : 0 ≤ C := by positivity
  set g : ℝ → ℝ := fun y => (y * (1 - y)) ^ (-(3 / 4 : ℝ)) with hg
  have hgi : IntervalIntegrable g volume 0 1 := V05.ii_tt _ (by norm_num)
  refine IsBigO.of_bound (C * ∫ y in (0:ℝ)..1, g y) ?_
  filter_upwards [Metric.ball_mem_nhds (0 : ℝ) hδ0] with u hu
  rw [Metric.mem_ball, dist_zero_right, Real.norm_eq_abs] at hu
  -- pointwise facts
  have hpt : ∀ y ∈ Ioo (0:ℝ) 1, |u * Lf lam y| ≤ V05.ls y / 4 := by
    intro y hy
    rw [abs_mul]
    have h1 := Lf_bd lam y
    have h2 := V05.ls_nonneg y
    have h3 : |u| * (|lam| + 1) ≤ 1 / 4 := by
      rw [hδ] at hu
      rw [lt_div_iff₀ (by positivity)] at hu
      nlinarith
    calc |u| * |Lf lam y| ≤ |u| * ((|lam| + 1) * V05.ls y) :=
          mul_le_mul_of_nonneg_left h1 (abs_nonneg u)
      _ = (|u| * (|lam| + 1)) * V05.ls y := by ring
      _ ≤ 1 / 4 * V05.ls y := mul_le_mul_of_nonneg_right h3 h2
      _ = V05.ls y / 4 := by ring
  have hexp : ∀ y ∈ Ioo (0:ℝ) 1, Real.exp |u * Lf lam y| ≤ (y * (1 - y)) ^ (-(1 / 4 : ℝ)) := by
    intro y hy
    rw [← exp_ls hy.1 hy.2]
    exact Real.exp_le_exp.mpr (by linarith [hpt y hy])
  have hint1 : IntervalIntegrable (fun y => wB p q y * Real.exp (u * Lf lam y)) volume 0 1 := by
    refine IntervalIntegrable.mono_fun' (V05.ii_tt (-(1 / 4 : ℝ)) (by norm_num))
      (by unfold wB Lf; fun_prop : Measurable _).aestronglyMeasurable ?_
    filter_upwards [V05.ae_Ioo01] with y hy
    obtain ⟨h0, h1⟩ := wB_bd p q hy.1 hy.2
    rw [Real.norm_eq_abs, abs_mul, abs_of_nonneg h0, abs_of_pos (Real.exp_pos _)]
    calc wB p q y * Real.exp (u * Lf lam y) ≤ 1 * Real.exp |u * Lf lam y| := by
          gcongr
          exact le_abs_self _
      _ ≤ _ := by rw [one_mul]; exact hexp y hy
  have hrw : (∫ y in (0:ℝ)..1, wB p q y * Real.exp (u * Lf lam y)) -
      ∑ k ∈ Finset.range (N + 1),
        (∫ y in (0:ℝ)..1, wB p q y * Lf lam y ^ k) / (k.factorial : ℝ) * u ^ k =
      ∫ y in (0:ℝ)..1, wB p q y * (Real.exp (u * Lf lam y) - expT N (u * Lf lam y)) := by
    have e : (fun y => wB p q y * (Real.exp (u * Lf lam y) - expT N (u * Lf lam y))) =
        fun y => wB p q y * Real.exp (u * Lf lam y) - ∑ k ∈ Finset.range (N + 1),
          (u ^ k / (k.factorial : ℝ)) * (wB p q y * Lf lam y ^ k) := by
      funext y
      rw [expT_real, mul_sub, Finset.mul_sum]
      congr 1
      refine Finset.sum_congr rfl fun k _ => ?_
      rw [mul_pow]; ring
    rw [e, intervalIntegral.integral_sub hint1
      (ii_sum _ fun k _ => (mom_int p q lam k).const_mul _),
      intervalIntegral.integral_finsetSum fun k _ => (mom_int p q lam k).const_mul _]
    congr 1
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [intervalIntegral.integral_const_mul]; ring
  rw [hrw, Real.norm_eq_abs (u ^ (N + 1)), abs_pow]
  have hfin : (∫ y in (0:ℝ)..1, |u| ^ (N + 1) * C * g y) =
      C * (∫ y in (0:ℝ)..1, g y) * |u| ^ (N + 1) := by
    rw [intervalIntegral.integral_const_mul]; ring
  rw [← hfin]
  refine intervalIntegral.norm_integral_le_of_norm_le zero_le_one ?_
    (hgi.const_mul _)
  filter_upwards [V05.ae_ne_one] with y hy1 hy
  have hy' : y ∈ Ioo (0:ℝ) 1 := ⟨hy.1, lt_of_le_of_ne hy.2 hy1⟩
  have hp : 0 < y * (1 - y) := mul_pos hy'.1 (by linarith [hy'.2])
  obtain ⟨h0, h1⟩ := wB_bd p q hy'.1 hy'.2
  have hr := exp_rem N (u * Lf lam y)
  have hls := V05.ls_pow_le hy' (N + 1)
  have hLp : |Lf lam y| ^ (N + 1) ≤ (|lam| + 1) ^ (N + 1) * V05.ls y ^ (N + 1) := by
    rw [← mul_pow]; exact pow_le_pow_left₀ (abs_nonneg _) (Lf_bd lam y) _
  have hA : |u * Lf lam y| ^ (N + 1) ≤ |u| ^ (N + 1) * (C * (y * (1 - y)) ^ (-(1 / 2 : ℝ))) := by
    rw [abs_mul, mul_pow]
    refine mul_le_mul_of_nonneg_left ?_ (by positivity)
    refine hLp.trans ?_
    rw [hC, mul_assoc]
    exact mul_le_mul_of_nonneg_left hls (by positivity)
  have hE := hexp y hy'
  have hq : (y * (1 - y)) ^ (-(1 / 2 : ℝ)) * (y * (1 - y)) ^ (-(1 / 4 : ℝ)) = g y := by
    rw [hg, ← Real.rpow_add hp]; norm_num
  rw [Real.norm_eq_abs, abs_mul, abs_of_nonneg h0]
  calc wB p q y * |Real.exp (u * Lf lam y) - expT N (u * Lf lam y)|
      ≤ 1 * (|u * Lf lam y| ^ (N + 1) * Real.exp |u * Lf lam y|) := by
        gcongr
    _ ≤ 1 * ((|u| ^ (N + 1) * (C * (y * (1 - y)) ^ (-(1 / 2 : ℝ)))) *
          (y * (1 - y)) ^ (-(1 / 4 : ℝ))) := by
        gcongr
    _ = |u| ^ (N + 1) * C * g y := by rw [← hq]; ring

end Q04

namespace Q04
open Asymptotics Filter Topology MeasureTheory Set

lemma I_beta (p q : ℕ) (hp : 1 ≤ p) (hq : 1 ≤ q) (lam u : ℝ) (h1 : 0 < (p : ℝ) + lam * u)
    (h2 : 0 < (q : ℝ) + u) :
    (∫ y in (0:ℝ)..1, wB p q y * Real.exp (u * Lf lam y)) *
        Real.Gamma ((p : ℝ) + q + (lam + 1) * u) =
      Real.Gamma ((p : ℝ) + lam * u) * Real.Gamma ((q : ℝ) + u) := by
  have key := X05.beta_real h1 h2
  rw [show (p : ℝ) + lam * u + ((q : ℝ) + u) = (p : ℝ) + q + (lam + 1) * u by ring] at key
  rw [← key]
  congr 1
  refine intervalIntegral.integral_congr_uIoo fun t ht => ?_
  rw [uIoo_of_le zero_le_one] at ht
  have ht0 := ht.1
  have ht1 : 0 < 1 - t := by linarith [ht.2]
  rw [show (p : ℝ) + lam * u - 1 = ((p - 1 : ℕ) : ℝ) + lam * u by
      push_cast [Nat.cast_sub hp]; ring,
    show (q : ℝ) + u - 1 = ((q - 1 : ℕ) : ℝ) + u by push_cast [Nat.cast_sub hq]; ring,
    Real.rpow_add ht0, Real.rpow_add ht1, Real.rpow_natCast, Real.rpow_natCast,
    Real.rpow_def_of_pos ht0, Real.rpow_def_of_pos ht1]
  unfold wB Lf
  rw [show u * (lam * Real.log t + Real.log (1 - t)) = Real.log t * (lam * u) +
    Real.log (1 - t) * u by ring, Real.exp_add]
  ring

lemma Gamma_shift (P : ℕ) {s : ℝ} (hs : -1 < s) :
    Real.Gamma ((P : ℝ) + 1 + s) =
      (P.factorial : ℝ) * Real.Gamma (1 + s) * ∏ j ∈ Finset.range P, (1 + s / ((j : ℝ) + 1)) := by
  induction P with
  | zero => simp
  | succ P ih =>
    have hP : (0 : ℝ) < (P : ℝ) + 1 + s := by have := P.cast_nonneg (α := ℝ); linarith
    rw [show ((P + 1 : ℕ) : ℝ) + 1 + s = ((P : ℝ) + 1 + s) + 1 by push_cast; ring,
      Real.Gamma_add_one hP.ne', ih, Finset.prod_range_succ, Nat.factorial_succ]
    have : (P : ℝ) + 1 ≠ 0 := by positivity
    push_cast
    field_simp

/-- `log Γ(1+s) + γ s` -/
noncomputable def Zf (s : ℝ) : ℝ := Real.log (Real.Gamma (1 + s)) + Real.eulerMascheroniConstant * s

/-- `Σ_{j<P} log(1 + s/(j+1))` -/
noncomputable def Hs (P : ℕ) (s : ℝ) : ℝ := ∑ j ∈ Finset.range P, Real.log (1 + s / ((j : ℝ) + 1))

lemma Gamma_exp (P : ℕ) {s : ℝ} (hs : -1 < s) :
    Real.Gamma ((P : ℝ) + 1 + s) = (P.factorial : ℝ) *
      Real.exp (Zf s - Real.eulerMascheroniConstant * s + Hs P s) := by
  rw [Gamma_shift P hs, Real.exp_add, Hs, Real.exp_sum, mul_assoc]
  congr 1
  have hG : 0 < Real.Gamma (1 + s) := Real.Gamma_pos_of_pos (by linarith)
  congr 1
  · rw [Zf, add_sub_cancel_right, Real.exp_log hG]
  · refine Finset.prod_congr rfl fun j _ => ?_
    have hj : (1 : ℝ) ≤ (j : ℝ) + 1 := by have := j.cast_nonneg (α := ℝ); linarith
    have : -1 < s / ((j : ℝ) + 1) := by
      rw [lt_div_iff₀ (by linarith)]; nlinarith
    rw [Real.exp_log (by linarith)]

/-- the logarithm of `B(a+1+λu, b+1+u)/B(a+1,b+1)` -/
noncomputable def DeltaF (a b : ℕ) (lam u : ℝ) : ℝ :=
  Zf (u * lam) + Hs a (u * lam) + Zf u + Hs b u - Zf (u * (lam + 1)) - Hs (a + b + 1) (u * (lam + 1))

lemma I_eq (a b : ℕ) (lam u : ℝ) (h1 : -1 < u * lam) (h2 : -1 < u) (h3 : -1 < u * (lam + 1)) :
    (∫ y in (0:ℝ)..1, wB (a + 1) (b + 1) y * Real.exp (u * Lf lam y)) =
      V05.Bet (a + 1) (b + 1) * Real.exp (DeltaF a b lam u) := by
  have key := I_beta (a + 1) (b + 1) (by omega) (by omega) lam u (by push_cast; linarith)
    (by push_cast; linarith)
  rw [show ((a + 1 : ℕ) : ℝ) + ((b + 1 : ℕ) : ℝ) + (lam + 1) * u =
      ((a + b + 1 : ℕ) : ℝ) + 1 + u * (lam + 1) by push_cast; ring,
    show ((a + 1 : ℕ) : ℝ) + lam * u = (a : ℝ) + 1 + u * lam by push_cast; ring,
    show ((b + 1 : ℕ) : ℝ) + u = (b : ℝ) + 1 + u by push_cast; ring,
    Gamma_exp _ h3, Gamma_exp _ h1, Gamma_exp _ h2] at key
  have hB : V05.Bet (a + 1) (b + 1) =
      (a.factorial * b.factorial : ℝ) / ((a + b + 1).factorial : ℝ) := by
    unfold V05.Bet
    rw [show a + 1 + (b + 1) - 1 = a + b + 1 by omega]
    simp
  rw [hB, DeltaF]
  have hE : 0 < Real.exp (Zf (u * (lam + 1)) - Real.eulerMascheroniConstant * (u * (lam + 1)) +
      Hs (a + b + 1) (u * (lam + 1))) := Real.exp_pos _
  have hF : (0 : ℝ) < ((a + b + 1).factorial : ℝ) := by positivity
  set I := ∫ y in (0:ℝ)..1, wB (a + 1) (b + 1) y * Real.exp (u * Lf lam y)
  have e : Zf (u * lam) + Hs a (u * lam) + Zf u + Hs b u - Zf (u * (lam + 1)) -
      Hs (a + b + 1) (u * (lam + 1)) =
      (Zf (u * lam) - Real.eulerMascheroniConstant * (u * lam) + Hs a (u * lam)) +
      (Zf u - Real.eulerMascheroniConstant * u + Hs b u) -
      (Zf (u * (lam + 1)) - Real.eulerMascheroniConstant * (u * (lam + 1)) +
        Hs (a + b + 1) (u * (lam + 1))) := by ring
  rw [e, Real.exp_sub, Real.exp_add]
  have hne : ((a + b + 1).factorial : ℝ) * Real.exp (Zf (u * (lam + 1)) -
      Real.eulerMascheroniConstant * (u * (lam + 1)) + Hs (a + b + 1) (u * (lam + 1))) ≠ 0 := by
    positivity
  rw [eq_div_of_mul_eq hne key]
  field_simp

lemma sum_Icc_two {M : Type*} [AddCommMonoid M] (g : ℕ → M) (N : ℕ) :
    ∑ x ∈ Finset.Icc 2 N, g x = ∑ k ∈ Finset.range (N - 1), g (k + 2) := by
  induction N with
  | zero => simp
  | succ N ih =>
    rcases Nat.eq_zero_or_pos N with h | h
    · subst h; simp
    · rw [Finset.sum_Icc_succ_top (by omega), ih, show N + 1 - 1 = (N - 1) + 1 by omega,
        Finset.sum_range_succ, show N - 1 + 2 = N + 1 by omega]

/-- the truncated zeta part `Σ_{w=2}^N (-1)^w ζ(w) s^w / w` -/
noncomputable def ZT (N : ℕ) (s : ℝ) : ℝ :=
  ∑ w ∈ Finset.Icc 2 N, (-1 : ℝ) ^ w * X05.z w * s ^ w / (w : ℝ)

lemma Z_trunc (N : ℕ) : (fun s => Zf s - ZT N s) =O[𝓝 0] (fun s => s ^ (N + 1)) := by
  refine IsBigO.of_bound 4 ?_
  filter_upwards [Metric.ball_mem_nhds (0 : ℝ) (by norm_num : (0 : ℝ) < 1 / 2)] with s hs
  rw [Metric.mem_ball, dist_zero_right, Real.norm_eq_abs] at hs
  have hS := Q04.proof_1 s (by linarith)
  set f : ℕ → ℝ := fun w => (-1 : ℝ) ^ (w + 2) * X05.z (w + 2) * s ^ (w + 2) / ((w : ℝ) + 2)
  have hZ : ZT N s = ∑ k ∈ Finset.range (N - 1), f k := by
    rw [ZT, sum_Icc_two]
    refine Finset.sum_congr rfl fun k _ => ?_
    simp [f]
  have h2 := (hasSum_nat_add_iff' (N - 1)).mpr hS
  have hb := (hasSum_geometric_two.mul_left (2 * |s| ^ (N + 1)))
  rw [hZ, Real.norm_eq_abs (s ^ (N + 1)), abs_pow, Zf]
  have := h2.norm_le_of_bounded hb fun k => ?_
  · simp only [Real.norm_eq_abs] at this ⊢; linarith
  · simp only [f, Real.norm_eq_abs]
    have hz := X05.z_bd (n := k + (N - 1) + 2) (by omega)
    have hk : (2 : ℝ) ≤ ((k + (N - 1) : ℕ) : ℝ) + 2 := by
      have := (k + (N - 1)).cast_nonneg (α := ℝ); linarith
    rw [abs_div, abs_mul, abs_mul, abs_pow, abs_pow, abs_neg, abs_one, one_pow, one_mul,
      abs_of_pos (by linarith : (0:ℝ) < ((k + (N - 1) : ℕ) : ℝ) + 2)]
    have hs1 : |s| ≤ 1 := by linarith
    have e1 : |s| ^ (k + (N - 1) + 2) ≤ |s| ^ (N + 1 + k) :=
      pow_le_pow_of_le_one (abs_nonneg s) hs1 (by omega)
    have e2 : |s| ^ k ≤ (1 / 2) ^ k := pow_le_pow_left₀ (abs_nonneg s) hs.le k
    have e3 : |s| ^ (N + 1 + k) ≤ |s| ^ (N + 1) * (1 / 2) ^ k := by
      rw [pow_add]; exact mul_le_mul_of_nonneg_left e2 (by positivity)
    have hp : 0 ≤ |s| ^ (k + (N - 1) + 2) := by positivity
    rw [div_le_iff₀ (by linarith)]
    calc |X05.z (k + (N - 1) + 2)| * |s| ^ (k + (N - 1) + 2)
        ≤ 2 * (|s| ^ (N + 1) * (1 / 2) ^ k) := by
          have := mul_le_mul hz (e1.trans e3) hp (by norm_num)
          linarith
      _ ≤ 2 * |s| ^ (N + 1) * (1 / 2) ^ k * (((k + (N - 1) : ℕ) : ℝ) + 2) := by
          have : 0 ≤ 2 * |s| ^ (N + 1) * (1 / 2) ^ k := by positivity
          nlinarith

/-- the truncated harmonic part `Σ_{r=1}^N (-1)^{r-1}/r H^{(r)}_P s^r` -/
noncomputable def HT (N P : ℕ) (s : ℝ) : ℝ :=
  ∑ r ∈ Finset.Icc 1 N, (-1 : ℝ) ^ (r - 1) / (r : ℝ) * V05.Hr r P * s ^ r

lemma HT_eq (N P : ℕ) (s : ℝ) :
    HT N P s = ∑ j ∈ Finset.range P, logT N (s * ((j : ℝ) + 1)⁻¹) := by
  unfold HT logT V05.Hr
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun r _ => ?_
  rw [Finset.mul_sum, Finset.sum_mul]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [Rat.smul_def]
  push_cast
  rw [mul_pow, inv_pow]
  ring

lemma Hs_trunc (N P : ℕ) : (fun s => Hs P s - HT N P s) =O[𝓝 0] (fun s => s ^ (N + 1)) := by
  have e : (fun s => Hs P s - HT N P s) = fun s => ∑ j ∈ Finset.range P,
      (Real.log (1 + s * ((j : ℝ) + 1)⁻¹) - logT N (s * ((j : ℝ) + 1)⁻¹)) := by
    funext s
    rw [HT_eq, Hs, ← Finset.sum_sub_distrib]
    simp only [div_eq_mul_inv]
  rw [e]
  have := IsBigO.sum (s := Finset.range P)
    (A := fun (j : ℕ) (s : ℝ) => Real.log (1 + s * ((j : ℝ) + 1)⁻¹) - logT N (s * ((j : ℝ) + 1)⁻¹))
    (g := fun s => s ^ (N + 1)) fun j _ => by
      have ht : Tendsto (fun s : ℝ => s * ((j : ℝ) + 1)⁻¹) (𝓝 0) (𝓝 0) := by
        simpa using (tendsto_id (x := 𝓝 (0 : ℝ))).mul_const ((j : ℝ) + 1)⁻¹
      exact ((log_sub_logT N).comp_tendsto ht).trans (mulpow _ (N + 1))
  refine this.congr_left fun s => ?_
  simp [Finset.sum_apply]

end Q04

namespace Q04
open Asymptotics Filter Topology MeasureTheory Set

/-- `x ↦ λu`, `x' ↦ u`, `𝔷_w ↦ ζ(w)`, as a polynomial in `u` -/
noncomputable def phiL (lam : ℝ) : S2 →+* Polynomial ℝ :=
  MvPolynomial.eval₂Hom (Polynomial.C.comp (MvPolynomial.aeval X05.z).toRingHom)
    ![Polynomial.C lam * Polynomial.X, Polynomial.X]

noncomputable def psiL (lam u : ℝ) : S2 →+* ℝ :=
  MvPolynomial.eval₂Hom (MvPolynomial.aeval X05.z).toRingHom ![u * lam, u]

lemma eval_phiL (lam u : ℝ) (F : S2) : (phiL lam F).eval u = psiL lam u F := by
  unfold phiL psiL
  rw [← Polynomial.coe_evalRingHom, MvPolynomial.coe_eval₂Hom, MvPolynomial.eval₂_comp_left,
    MvPolynomial.coe_eval₂Hom]
  congr 1
  · exact RingHom.ext fun c => by
      simp only [RingHom.comp_apply, Polynomial.coe_evalRingHom, Polynomial.eval_C]
  · funext i
    fin_cases i
    · show Polynomial.eval u (Polynomial.C lam * Polynomial.X) = u * lam
      rw [Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_X, mul_comm]
    · show Polynomial.eval u Polynomial.X = u
      rw [Polynomial.eval_X]

/-- the truncation of `DeltaF` -/
noncomputable def TF (a b N : ℕ) (lam u : ℝ) : ℝ :=
  ZT N (u * lam) + HT N a (u * lam) + ZT N u + HT N b u - ZT N (u * (lam + 1)) -
    HT N (a + b + 1) (u * (lam + 1))

lemma psiL_DeltaT (a b N : ℕ) (lam u : ℝ) :
    psiL lam u (DeltaT N (a + 1) (b + 1)) = TF a b N lam u := by
  unfold DeltaT
  rw [show a + 1 + (b + 1) - 1 = a + b + 1 by omega]
  simp only [map_add, map_sub, map_sum, map_mul, map_pow, map_rat_smul, psiL,
    MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C, Nat.add_sub_cancel]
  unfold TF ZT HT
  rw [show ∀ A B C D E F : ℝ, A + B + C + D - E - F = (B + D - F) + (A + C - E) by intros; ring]
  rw [← Finset.sum_add_distrib, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib,
    ← Finset.sum_sub_distrib]
  congr 1
  · refine Finset.sum_congr rfl fun r _ => ?_
    simp only [Rat.smul_def, V05.HrQ_cast, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.head_cons]
    push_cast
    ring
  · refine Finset.sum_congr rfl fun w _ => ?_
    simp only [Rat.smul_def, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons,
      AlgHom.toRingHom_eq_coe, RingHom.coe_coe, MvPolynomial.aeval_X]
    push_cast
    ring

lemma tailR (P : Polynomial ℝ) (K : ℕ) :
    (fun u : ℝ => P.eval u - ∑ k ∈ Finset.range (K + 1), P.coeff k * u ^ k) =O[𝓝 0]
      (fun u => u ^ (K + 1)) := by
  set N := P.natDegree + K + 2
  have hN : P.natDegree < N := by omega
  have hKN : K + 1 ≤ N := by omega
  have hg : Continuous (fun u : ℝ => ∑ k ∈ Finset.Ico (K + 1) N, P.coeff k * u ^ (k - (K + 1))) := by
    fun_prop
  have hb := (hg.tendsto 0).isBigO_one ℝ
  have hm := (isBigO_refl (fun u : ℝ => u ^ (K + 1)) (𝓝 0)).mul hb
  simp only [mul_one] at hm
  refine hm.congr_left fun u => ?_
  rw [Polynomial.eval_eq_sum_range' hN, ← Finset.sum_range_add_sum_Ico _ hKN, add_sub_cancel_left,
    Finset.mul_sum]
  refine Finset.sum_congr rfl fun k hk => ?_
  have := (Finset.mem_Ico.mp hk).1
  rw [mul_left_comm, ← pow_add, Nat.add_sub_cancel' this]

/-- the index `x^i x'^{k-i}` -/
noncomputable def ix (k i : ℕ) : Fin 2 →₀ ℕ := Finsupp.single 0 i + Finsupp.single 1 (k - i)

lemma phiL_coeff (lam : ℝ) (F : S2) (k : ℕ) :
    (phiL lam F).coeff k = ∑ i ∈ Finset.range (k + 1),
      MvPolynomial.aeval X05.z (F.coeff (ix k i)) * lam ^ i := by
  induction F using MvPolynomial.induction_on' with
  | monomial t c =>
    have hm : phiL lam (MvPolynomial.monomial t c) =
        Polynomial.C (MvPolynomial.aeval X05.z c * lam ^ t 0) * Polynomial.X ^ (t 0 + t 1) := by
      rw [MvPolynomial.monomial_eq, Finsupp.prod_fintype _ _ (fun i => by simp)]
      simp only [map_mul, map_prod, map_pow, Fin.prod_univ_two, phiL, MvPolynomial.eval₂Hom_C,
        MvPolynomial.eval₂Hom_X', Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons,
        RingHom.coe_comp, Function.comp_apply, Polynomial.C_mul, Polynomial.C_pow, mul_pow,
        pow_add]
      simp only [AlgHom.toRingHom_eq_coe, RingHom.coe_coe]
      ring
    rw [hm, Polynomial.coeff_C_mul_X_pow]
    simp only [MvPolynomial.coeff_monomial]
    split_ifs with hk
    · rw [Finset.sum_eq_single (t 0)]
      · have : t = ix k (t 0) := by
          ext j; fin_cases j <;> simp [ix] <;> omega
        rw [if_pos this]
      · intro i _ hi
        rw [if_neg]
        · simp
        · intro h; apply hi; rw [h]; simp [ix]
      · intro h; exfalso; apply h; rw [Finset.mem_range]; omega
    · symm
      refine Finset.sum_eq_zero fun i hi => ?_
      rw [if_neg]
      · simp
      · intro h; apply hk; rw [h]; simp [ix]; rw [Finset.mem_range] at hi; omega
  | add F G hF hG =>
    rw [map_add, Polynomial.coeff_add, hF, hG, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [MvPolynomial.coeff_add, map_add, add_mul]

end Q04

namespace Q04
open Asymptotics Filter Topology MeasureTheory Set

lemma comp_mul {F G : ℝ → ℝ} {N : ℕ} (h : (fun s => F s - G s) =O[𝓝 0] (fun s => s ^ (N + 1)))
    (c : ℝ) : (fun u => F (u * c) - G (u * c)) =O[𝓝 0] (fun u => u ^ (N + 1)) := by
  have ht : Tendsto (fun u : ℝ => u * c) (𝓝 0) (𝓝 0) := by
    simpa using (tendsto_id (x := 𝓝 (0 : ℝ))).mul_const c
  exact (h.comp_tendsto ht).trans (mulpow c (N + 1))

lemma Delta_trunc (a b N : ℕ) (lam : ℝ) :
    (fun u => DeltaF a b lam u - TF a b N lam u) =O[𝓝 0] (fun u => u ^ (N + 1)) := by
  have h1 := comp_mul (Z_trunc N) lam
  have h2 := comp_mul (Hs_trunc N a) lam
  have h3 := Z_trunc N
  have h4 := Hs_trunc N b
  have h5 := comp_mul (Z_trunc N) (lam + 1)
  have h6 := comp_mul (Hs_trunc N (a + b + 1)) (lam + 1)
  refine (((((h1.add h2).add h3).add h4).sub h5).sub h6).congr_left fun u => ?_
  simp only [DeltaF, TF]
  ring

lemma TF_zero (a b : ℕ) (lam u : ℝ) : TF a b 0 lam u = 0 := by
  simp [TF, ZT, HT]

/-- the formal truncated `exp(Δ_⋆)` at `(p,q) = (a+1,b+1)` -/
noncomputable def GT (a b N : ℕ) : S2 := expT N (DeltaT N (a + 1) (b + 1))

lemma main_bigO (a b N : ℕ) (lam : ℝ) :
    (fun u : ℝ => ∑ k ∈ Finset.range (N + 1),
      ((∫ y in (0:ℝ)..1, wB (a + 1) (b + 1) y * Lf lam y ^ k) / (k.factorial : ℝ) -
        V05.Bet (a + 1) (b + 1) * (phiL lam (GT a b N)).coeff k) * u ^ k)
      =O[𝓝 0] (fun u => u ^ (N + 1)) := by
  have hI := int_taylor (a + 1) (b + 1) lam N
  have hEq : ∀ᶠ u in 𝓝 (0 : ℝ),
      (∫ y in (0:ℝ)..1, wB (a + 1) (b + 1) y * Real.exp (u * Lf lam y)) =
        V05.Bet (a + 1) (b + 1) * Real.exp (DeltaF a b lam u) := by
    have hl : 0 < |lam| + 2 := by positivity
    filter_upwards [Metric.ball_mem_nhds (0 : ℝ) (by positivity : (0 : ℝ) < 1 / (|lam| + 2))]
      with u hu
    rw [Metric.mem_ball, dist_zero_right, Real.norm_eq_abs, lt_div_iff₀ hl] at hu
    have hu1 : |u| < 1 := by nlinarith [abs_nonneg u, abs_nonneg lam]
    have hb1 : |u * lam| < 1 := by rw [abs_mul]; nlinarith [abs_nonneg u, abs_nonneg lam]
    have hb2 : |u * (lam + 1)| < 1 := by
      rw [abs_mul]
      have := abs_add_le lam 1
      rw [abs_one] at this
      nlinarith [abs_nonneg u, abs_nonneg lam, abs_nonneg (lam + 1)]
    exact I_eq a b lam u (by linarith [neg_abs_le (u * lam)]) (by linarith [neg_abs_le u])
      (by linarith [neg_abs_le (u * (lam + 1))])
  have hD := Delta_trunc a b N lam
  have hD0 : (fun u => DeltaF a b lam u) =O[𝓝 0] (fun u => u) := by
    have := Delta_trunc a b 0 lam
    simpa [TF_zero] using this
  have huO1 : (fun u : ℝ => u) =O[𝓝 0] (fun _ => (1 : ℝ)) :=
    (tendsto_id (x := 𝓝 (0 : ℝ))).isBigO_one ℝ
  have hD1 := hD0.trans huO1
  have hT1 : (fun u => TF a b N lam u) =O[𝓝 0] (fun _ => (1 : ℝ)) := by
    have hu1 : (fun u : ℝ => u ^ (N + 1)) =O[𝓝 0] (fun _ => (1 : ℝ)) := by
      have : Tendsto (fun u : ℝ => u ^ (N + 1)) (𝓝 0) (𝓝 0) := by
        have := ((continuous_pow (N + 1)).tendsto (0 : ℝ)); simpa using this
      exact this.isBigO_one ℝ
    exact (hD1.sub (hD.trans hu1)).congr_left fun u => by ring
  have hE := ((exp_sub_expT N).comp_tendsto (hD0.trans_tendsto tendsto_id)).trans (hD0.pow (N + 1))
  have hE2 := expT_diff hD hD1 hT1 N
  have hP : ∀ u : ℝ, expT N (TF a b N lam u) = (phiL lam (GT a b N)).eval u := by
    intro u
    rw [eval_phiL, GT, expT_map (A := S2) (B := ℝ) (psiL lam u), psiL_DeltaT]
  have htail := tailR (phiL lam (GT a b N)) N
  have hsum := ((hI.neg_left).add (((hE.add hE2).add htail).const_mul_left
    (V05.Bet (a + 1) (b + 1))))
  refine hsum.congr' ?_ EventuallyEq.rfl
  filter_upwards [hEq] with u hu
  simp only [Function.comp, Pi.add_apply, Pi.neg_apply]
  rw [← hP u, hu]
  simp only [sub_mul, Finset.sum_sub_distrib, mul_assoc, ← Finset.mul_sum]
  ring

lemma coeff_ident (a b N : ℕ) (lam : ℝ) :
    (∫ y in (0:ℝ)..1, wB (a + 1) (b + 1) y * Lf lam y ^ N) / (N.factorial : ℝ) =
      V05.Bet (a + 1) (b + 1) * (phiL lam (GT a b N)).coeff N := by
  have h := (main_bigO a b N lam).trans_isLittleO (isLittleO_pow_pow (Nat.lt_succ_self N))
  have := poly_littleO_zero N _ (h.mono nhdsWithin_le_nhds) N le_rfl
  linarith

lemma Ib_formula (a b i m : ℕ) :
    V05.Ib (a + 1) (b + 1) i m =
      V05.Bet (a + 1) (b + 1) * MvPolynomial.aeval X05.z (DsharpF i m (a + 1) (b + 1)) := by
  set N := i + m with hN
  set e : ℕ → ℝ := fun j => (N.choose j : ℝ) * V05.Ib (a + 1) (b + 1) j (N - j) /
    (N.factorial : ℝ) - V05.Bet (a + 1) (b + 1) *
      MvPolynomial.aeval X05.z ((GT a b N).coeff (ix N j)) with he
  have key : ∀ lam : ℝ, ∑ j ∈ Finset.range (N + 1), e j * lam ^ j = 0 := by
    intro lam
    have h := coeff_ident a b N lam
    rw [moment_eq, phiL_coeff] at h
    simp only [he, sub_mul, Finset.sum_sub_distrib]
    rw [sub_eq_zero]
    calc _ = (∑ j ∈ Finset.range (N + 1), lam ^ j * (N.choose j : ℝ) *
          V05.Ib (a + 1) (b + 1) j (N - j)) / (N.factorial : ℝ) := by
          rw [Finset.sum_div]; exact Finset.sum_congr rfl fun j _ => by ring
      _ = _ := h
      _ = _ := by rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun j _ => by ring
  have hz := poly_littleO_zero N e ((isLittleO_zero (fun u : ℝ => u ^ N) (𝓝[≠] 0)).congr_left
    fun u => (key u).symm) i (by omega)
  simp only [he, sub_eq_zero] at hz
  rw [show N - i = m by omega] at hz
  have hix : ix N i = Finsupp.single 0 i + Finsupp.single 1 m := by
    simp [ix, show N - i = m by omega]
  rw [hix] at hz
  have hD : MvPolynomial.aeval X05.z (DsharpF i m (a + 1) (b + 1)) =
      (i.factorial : ℝ) * (m.factorial : ℝ) *
        MvPolynomial.aeval X05.z ((GT a b N).coeff (Finsupp.single 0 i + Finsupp.single 1 m)) := by
    rw [DsharpF, map_rat_smul, Rat.smul_def, GT]
    push_cast
    ring
  have hc : ((N.choose i : ℕ) : ℝ) * (i.factorial : ℝ) * (m.factorial : ℝ) = (N.factorial : ℝ) := by
    have := Nat.choose_mul_factorial_mul_factorial (show i ≤ N by omega)
    rw [show N - i = m by omega] at this
    exact_mod_cast this
  have hF : (0 : ℝ) < (N.factorial : ℝ) := by positivity
  rw [hD]
  rw [div_eq_iff hF.ne'] at hz
  rw [← hc] at hz
  have hi : (0 : ℝ) < (i.factorial : ℝ) * (m.factorial : ℝ) := by positivity
  have hC : (0 : ℝ) < ((N.choose i : ℕ) : ℝ) := by
    have := Nat.choose_pos (show i ≤ N by omega); exact_mod_cast this
  apply mul_left_cancel₀ hC.ne'
  linear_combination hz

end Q04

theorem Q04.proof_4 : Q04.stmt_4 := by
  intro p q hp hq i m
  obtain ⟨a, rfl⟩ : ∃ a, p = a + 1 := ⟨p - 1, by omega⟩
  obtain ⟨b, rfl⟩ : ∃ b, q = b + 1 := ⟨q - 1, by omega⟩
  exact ⟨V05.ibf_int _ _ _ _, Q04.Ib_formula a b i m, Q04.DsharpF_weight i m _ _⟩
namespace Q04
open MeasureTheory Set Filter Topology

lemma rhoQ_cast (j n : ℕ) : ((V05.rhoQ j n : ℚ) : ℝ) = V05.rho j n := by
  unfold V05.rhoQ V05.rho; push_cast; rfl

lemma Ib_eq {p q : ℕ} (hp : 1 ≤ p) (hq : 1 ≤ q) (i m : ℕ) :
    V05.Ib p q i m = V05.Bet p q * MvPolynomial.aeval X05.z (DsharpF i m p q) :=
  (Q04.proof_4 p q hp hq i m).2.1

lemma Jk_int {N : ℕ} (hN : 1 ≤ N) (k : ℕ) :
    (∫ y in (0:ℝ)..1, (y * (1 - y)) ^ (N - 1) *
        (iteratedDeriv k (fun u => V05.E N u y) 0 / (k.factorial : ℝ))) =
      ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
        (eF k N j i (k - i) : ℝ) * V05.Ib (N + j) N i (k - i) := by
  have e : (∫ y in (0:ℝ)..1, (y * (1 - y)) ^ (N - 1) *
        (iteratedDeriv k (fun u => V05.E N u y) 0 / (k.factorial : ℝ))) =
      ∫ y in (0:ℝ)..1, ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
        (eF k N j i (k - i) : ℝ) * V05.ibf (N + j) N i (k - i) y := by
    refine intervalIntegral.integral_congr_uIoo fun y hy => ?_
    rw [uIoo_of_le zero_le_one] at hy
    rw [stmt3_c k N hy.1 hy.2, Finset.mul_sum]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    unfold V05.ibf
    rw [show N + j - 1 = (N - 1) + j by omega, pow_add, mul_pow]
    ring
  rw [e, intervalIntegral.integral_finsetSum fun j _ =>
    ii_sum _ fun i _ => (V05.ibf_int _ _ _ _).const_mul _]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [intervalIntegral.integral_finsetSum fun i _ => (V05.ibf_int _ _ _ _).const_mul _]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [intervalIntegral.integral_const_mul]; rfl

lemma aeval_Pf (k N : ℕ) : MvPolynomial.aeval X05.z (Pf k N) =
    -(1 / ((N : ℝ) * (2 * (N : ℝ) - 1))) * ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
      ((eF k N j i (k - i) : ℝ) * V05.rho j N) *
        MvPolynomial.aeval X05.z (DsharpF i (k - i) (N + j) N) := by
  unfold Pf
  simp only [map_rat_smul, map_sum, Rat.smul_def]
  push_cast
  simp only [rhoQ_cast]

lemma term_eq {N : ℕ} (hN : 1 ≤ N) (k : ℕ) :
    V05.an N * (iteratedDeriv k (V05.J N) 0 / (k.factorial : ℝ)) =
      -MvPolynomial.aeval X05.z (Pf k N) := by
  rw [(V05.proof_5 N hN k).2, Jk_int hN, aeval_Pf]
  have hS : ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
      (eF k N j i (k - i) : ℝ) * V05.Ib (N + j) N i (k - i) =
      V05.Bet N N * ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
        ((eF k N j i (k - i) : ℝ) * V05.rho j N) *
          MvPolynomial.aeval X05.z (DsharpF i (k - i) (N + j) N) := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Ib_eq (by omega) hN, V05.Bet_rho hN]
    ring
  rw [hS, ← V05.an_Bet hN]
  ring

end Q04

theorem Q04.proof_5 : Q04.stmt_5 := by
  refine ⟨Q04.Pf_weight, fun α hα0 hα1 k => ?_⟩
  have hf : ∀ n : ℕ, MvPolynomial.aeval X05.z (Q04.Pf k (n + 1)) * D04.c α ^ (2 * (n + 1)) =
      -(V05.an (n + 1) * D04.c α ^ (2 * (n + 1)) *
        (iteratedDeriv k (V05.J (n + 1)) 0 / (k.factorial : ℝ))) := by
    intro n
    have := Q04.term_eq (N := n + 1) (by omega) k
    linear_combination (D04.c α ^ (2 * (n + 1))) * this
  obtain ⟨δ, hδ0, hδ1, -, hS4⟩ := V05.proof_4 α hα0 hα1
  refine ⟨?_, ?_⟩
  · obtain ⟨hc0, hc1, -, -, -⟩ := D04.cs_facts hα0 hα1
    obtain ⟨δ', hδ0', hδ1', hr⟩ := V05.exists_delta hc0 hc1
    obtain ⟨C, hC0, hb⟩ := V05.trm_bd k
    have hJ : ∀ (n : ℕ) (u : ℝ), |u| < 1 →
        V05.J (n + 1) u = V05.comb (V05.PhiF (n + 1)) (V05.LL (n + 1) 0) u := by
      intro n u hu
      rw [show V05.LL (n + 1) 0 = [((1 : ℝ) / 2, (0, 0, 0))] from rfl, V05.comb_L0,
        V05.proof_3 (n + 1) (by omega) u hu]
      rfl
    have hJk : ∀ n : ℕ, iteratedDeriv k (V05.J (n + 1)) 0
        = V05.comb (V05.PhiF (n + 1)) (V05.LL (n + 1) k) 0 := by
      intro n
      have hevJ : V05.J (n + 1) =ᶠ[nhds 0] V05.comb (V05.PhiF (n + 1)) (V05.LL (n + 1) 0) :=
        Filter.eventually_of_mem (Ioo_mem_nhds (show (-1 : ℝ) < 0 by norm_num)
          (show (0 : ℝ) < 1 by norm_num)) fun u hu => hJ n u (abs_lt.mpr hu)
      rw [hevJ.iteratedDeriv_eq]
      have h0' : (0 : ℝ) ∈ Set.Ioo (-V05.dl (n + 1)) (V05.dl (n + 1)) := by
        unfold V05.dl; constructor <;> norm_num
      exact V05.comb_iter_local (n + 1) isOpen_Ioo
        (fun z u hu => V05.PhiF_hasDerivAt (by omega) z hu) _ k 0 h0'
    have hr0 : 0 ≤ D04.c α ^ 2 * V05.rho4 δ' := by unfold V05.rho4; positivity
    have hk : (0 : ℝ) < (k.factorial : ℝ) := by positivity
    refine (V05.maj_summable k hr0 hr (C / (k.factorial : ℝ))).of_nonneg_of_le
      (fun n => abs_nonneg _) fun n => ?_
    rw [hf n, abs_neg, hJk n, mul_div_assoc', abs_div, abs_of_pos hk]
    have := hb (D04.c α) δ' hc0 hδ0'.le hδ1' n 0 (by rw [abs_zero]; exact hδ0'.le)
    unfold V05.trm at this
    rw [div_le_iff₀ hk]
    calc _ ≤ _ := this
      _ = _ := by field_simp
  · have h := (hS4 k).neg
    rw [neg_neg] at h
    exact h.congr_fun fun n => by rw [hf n]

namespace Q04
open Finsupp

lemma mono_bound (s : ℕ →₀ ℕ) (w : ℕ) : s w * w ≤ s.sum (fun w e => e * w) := by
  by_cases hw : w ∈ s.support
  · exact Finset.single_le_sum (f := fun w => s w * w) (fun _ _ => Nat.zero_le _) hw
  · rw [Finsupp.notMem_support_iff.mp hw]; simp

lemma cases5 (s : ℕ →₀ ℕ) (h1 : ∀ w ∈ s.support, 2 ≤ w) (h2 : (s.sum fun w e => e * w) ≤ 5) :
    s ∈ ({0, single 2 1, single 3 1, single 4 1, single 5 1, single 2 2,
      single 2 1 + single 3 1} : Finset (ℕ →₀ ℕ)) := by
  have hb : ∀ w, s w * w ≤ 5 := fun w => (mono_bound s w).trans h2
  have hz : ∀ w, (w < 2 ∨ 5 < w) → s w = 0 := by
    intro w hw
    by_contra hne
    have h3 := h1 w (Finsupp.mem_support_iff.mpr hne)
    have h4 := hb w
    have h5 : 1 ≤ s w := Nat.one_le_iff_ne_zero.mpr hne
    rcases hw with hw | hw
    · omega
    · nlinarith
  have hsub : s.support ⊆ ({2, 3, 4, 5} : Finset ℕ) := fun w hw => by
    have h3 := h1 w hw
    have hne := Finsupp.mem_support_iff.mp hw
    have : w ≤ 5 := by by_contra hc; exact hne (hz w (Or.inr (by omega)))
    simp only [Finset.mem_insert, Finset.mem_singleton]; omega
  have hsum : (s.sum fun w e => e * w) = s 2 * 2 + s 3 * 3 + s 4 * 4 + s 5 * 5 := by
    rw [Finsupp.sum_of_support_subset s hsub (fun w e => e * w) (fun _ _ => by simp)]
    simp [Finset.sum_insert]; ring
  have hs : s = single 2 (s 2) + single 3 (s 3) + single 4 (s 4) + single 5 (s 5) := by
    ext w
    by_cases hw : w < 2 ∨ 5 < w
    · rw [hz w hw]; simp only [Finsupp.coe_add, Pi.add_apply, Finsupp.single_apply]
      split_ifs <;> omega
    · have : w = 2 ∨ w = 3 ∨ w = 4 ∨ w = 5 := by omega
      rcases this with rfl | rfl | rfl | rfl <;> simp [Finsupp.single_apply]
  rw [hsum] at h2
  rw [hs]
  generalize s 2 = a at h2
  generalize s 3 = b at h2
  generalize s 4 = c at h2
  generalize s 5 = d at h2
  have ha : a ≤ 2 := by omega
  have hb' : b ≤ 1 := by omega
  have hc : c ≤ 1 := by omega
  have hd : d ≤ 1 := by omega
  interval_cases a <;> interval_cases b <;> interval_cases c <;> interval_cases d <;>
    first | omega | simp

lemma fs_ne {f g : ℕ →₀ ℕ} (h : f 2 ≠ g 2 ∨ f 3 ≠ g 3 ∨ f 4 ≠ g 4 ∨ f 5 ≠ g 5) : f ≠ g := by
  intro e; subst e; simp at h

lemma aeval_S5 (P : MvPolynomial ℕ ℚ)
    (h : ∀ s ∈ P.support, (∀ w ∈ s.support, 2 ≤ w) ∧ (s.sum fun w e => e * w) ≤ 5) :
    MvPolynomial.aeval X05.z P = (P.coeff 0 : ℝ) + (P.coeff (single 2 1) : ℝ) * X05.z 2 +
      (P.coeff (single 3 1) : ℝ) * X05.z 3 + (P.coeff (single 4 1) : ℝ) * X05.z 4 +
      (P.coeff (single 5 1) : ℝ) * X05.z 5 + (P.coeff (single 2 2) : ℝ) * X05.z 2 ^ 2 +
      (P.coeff (single 2 1 + single 3 1) : ℝ) * (X05.z 2 * X05.z 3) := by
  rw [MvPolynomial.aeval_eq_eval₂Hom, MvPolynomial.coe_eval₂Hom, MvPolynomial.eval₂_eq]
  rw [Finset.sum_subset (fun s hs => cases5 s (h s hs).1 (h s hs).2)
    (fun s _ hs => by rw [MvPolynomial.notMem_support_iff.mp hs]; simp)]
  repeat rw [Finset.sum_insert (by
    simp only [Finset.mem_insert, Finset.mem_singleton, not_or]
    repeat' apply And.intro
    all_goals exact fs_ne (by simp [Finsupp.single_apply]))]
  rw [Finset.sum_singleton]
  simp [Finsupp.support_single, Finsupp.support_add_eq, Finsupp.single_apply]
  ring

lemma wsum_single (w e : ℕ) : ((Finsupp.single w e).sum fun w e => e * w) = e * w := by
  rw [Finsupp.sum_single_index (by simp)]

end Q04

namespace Q04
/-- a submultiplicative seminorm compatible with the `ℚ`-action -/
structure IsSN {R : Type*} [CommRing R] [Algebra ℚ R] (ν : R → ℝ) : Prop where
  nonneg : ∀ a, 0 ≤ ν a
  zero : ν 0 = 0
  add : ∀ a b, ν (a + b) ≤ ν a + ν b
  mul : ∀ a b, ν (a * b) ≤ ν a * ν b
  smul : ∀ (c : ℚ) a, ν (c • a) ≤ |(c : ℝ)| * ν a
  one : ν 1 ≤ 1

section
variable {R : Type*} [CommRing R] [Algebra ℚ R] {ν : R → ℝ}

lemma IsSN.sum (h : IsSN ν) {ι : Type*} (s : Finset ι) (f : ι → R) :
    ν (∑ i ∈ s, f i) ≤ ∑ i ∈ s, ν (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [h.zero]
  | insert a s ha ih =>
    rw [Finset.sum_insert ha, Finset.sum_insert ha]
    exact (h.add _ _).trans (by linarith)

lemma IsSN.pow (h : IsSN ν) (a : R) (r : ℕ) : ν (a ^ r) ≤ ν a ^ r := by
  induction r with
  | zero => simpa using h.one
  | succ r ih =>
    rw [pow_succ, pow_succ]
    exact (h.mul _ _).trans (mul_le_mul_of_nonneg_right ih (h.nonneg _))

lemma IsSN.neg (h : IsSN ν) (a : R) : ν (-a) ≤ ν a := by
  have := h.smul (-1) a
  simpa using this

lemma IsSN.sub (h : IsSN ν) (a b : R) : ν (a - b) ≤ ν a + ν b := by
  rw [sub_eq_add_neg]; exact (h.add _ _).trans (by linarith [h.neg b])

lemma IsSN.expT_le (h : IsSN ν) (K : ℕ) (a : R) :
    ν (expT K a) ≤ ∑ r ∈ Finset.range (K + 1), ν a ^ r := by
  unfold expT
  refine (h.sum _ _).trans (Finset.sum_le_sum fun r _ => ?_)
  refine (h.smul _ _).trans ?_
  have : |((((r.factorial : ℚ)⁻¹ : ℚ)) : ℝ)| ≤ 1 := by
    push_cast
    rw [abs_of_nonneg (by positivity)]
    exact inv_le_one_of_one_le₀ (by exact_mod_cast r.factorial_pos)
  calc _ ≤ 1 * ν (a ^ r) := mul_le_mul_of_nonneg_right this (h.nonneg _)
    _ ≤ _ := by rw [one_mul]; exact h.pow a r

lemma IsSN.logT_le (h : IsSN ν) (K : ℕ) (a : R) :
    ν (logT K a) ≤ ∑ m ∈ Finset.Icc 1 K, ν a ^ m := by
  unfold logT
  refine (h.sum _ _).trans (Finset.sum_le_sum fun m hm => ?_)
  refine (h.smul _ _).trans ?_
  have hm1 : (1 : ℝ) ≤ m := by exact_mod_cast (Finset.mem_Icc.mp hm).1
  have : |((((-1 : ℚ) ^ (m - 1) / (m : ℚ) : ℚ)) : ℝ)| ≤ 1 := by
    push_cast
    rw [abs_div, abs_pow, abs_neg, abs_one, one_pow, abs_of_pos (by linarith)]
    exact (div_le_one (by linarith)).mpr hm1
  calc _ ≤ 1 * ν (a ^ m) := mul_le_mul_of_nonneg_right this (h.nonneg _)
    _ ≤ _ := by rw [one_mul]; exact h.pow a m

end

lemma geom_le {x : ℝ} (hx : 1 ≤ x) (K : ℕ) :
    ∑ r ∈ Finset.range (K + 1), x ^ r ≤ ((K : ℝ) + 1) * x ^ K := by
  calc ∑ r ∈ Finset.range (K + 1), x ^ r ≤ ∑ r ∈ Finset.range (K + 1), x ^ K :=
        Finset.sum_le_sum fun r hr => pow_le_pow_right₀ hx (Nat.lt_succ_iff.mp (Finset.mem_range.mp hr))
    _ = _ := by simp

/-- the `ℓ¹` norm of an `MvPolynomial` over a seminormed coefficient ring -/
noncomputable def Nmv {R σ : Type*} [CommRing R] (ν : R → ℝ) (p : MvPolynomial σ R) : ℝ :=
  ∑ s ∈ p.support, ν (p.coeff s)

section
variable {R σ : Type*} [CommRing R] [Algebra ℚ R] {ν : R → ℝ}

lemma Nmv_coeff (h : IsSN ν) (p : MvPolynomial σ R) (s : σ →₀ ℕ) : ν (p.coeff s) ≤ Nmv ν p := by
  by_cases hs : s ∈ p.support
  · exact Finset.single_le_sum (f := fun t => ν (p.coeff t)) (fun t _ => h.nonneg _) hs
  · rw [MvPolynomial.notMem_support_iff.mp hs, h.zero]
    exact Finset.sum_nonneg fun _ _ => h.nonneg _

lemma Nmv_nonneg (h : IsSN ν) (p : MvPolynomial σ R) : 0 ≤ Nmv ν p :=
  Finset.sum_nonneg fun _ _ => h.nonneg _

lemma Nmv_add (h : IsSN ν) (p q : MvPolynomial σ R) : Nmv ν (p + q) ≤ Nmv ν p + Nmv ν q := by
  classical
  unfold Nmv
  have e1 : ∑ s ∈ p.support, ν (p.coeff s) = ∑ s ∈ p.support ∪ q.support, ν (p.coeff s) :=
    Finset.sum_subset Finset.subset_union_left fun s _ hs => by
      rw [MvPolynomial.notMem_support_iff.mp hs, h.zero]
  have e2 : ∑ s ∈ q.support, ν (q.coeff s) = ∑ s ∈ p.support ∪ q.support, ν (q.coeff s) :=
    Finset.sum_subset Finset.subset_union_right fun s _ hs => by
      rw [MvPolynomial.notMem_support_iff.mp hs, h.zero]
  rw [e1, e2, ← Finset.sum_add_distrib]
  calc ∑ s ∈ (p + q).support, ν ((p + q).coeff s)
      ≤ ∑ s ∈ p.support ∪ q.support, ν ((p + q).coeff s) :=
        Finset.sum_le_sum_of_subset_of_nonneg MvPolynomial.support_add fun _ _ _ => h.nonneg _
    _ ≤ _ := Finset.sum_le_sum fun s _ => by rw [MvPolynomial.coeff_add]; exact h.add _ _

lemma Nmv_zero : Nmv ν (0 : MvPolynomial σ R) = 0 := by simp [Nmv]

lemma Nmv_sum (h : IsSN ν) {ι : Type*} (s : Finset ι) (f : ι → MvPolynomial σ R) :
    Nmv ν (∑ i ∈ s, f i) ≤ ∑ i ∈ s, Nmv ν (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [Nmv_zero]
  | insert a s ha ih =>
    rw [Finset.sum_insert ha, Finset.sum_insert ha]
    exact (Nmv_add h _ _).trans (by linarith)

lemma Nmv_monomial (h : IsSN ν) (s : σ →₀ ℕ) (a : R) : Nmv ν (MvPolynomial.monomial s a) ≤ ν a := by
  classical
  unfold Nmv
  rw [MvPolynomial.support_monomial]
  split_ifs with ha
  · simp [ha, h.zero]
  · simp

lemma Nmv_mul (h : IsSN ν) (p q : MvPolynomial σ R) : Nmv ν (p * q) ≤ Nmv ν p * Nmv ν q := by
  classical
  conv_lhs => rw [p.as_sum, q.as_sum, Finset.sum_mul_sum]
  simp only [MvPolynomial.monomial_mul]
  refine (Nmv_sum h _ _).trans ?_
  unfold Nmv
  rw [Finset.sum_mul_sum]
  refine Finset.sum_le_sum fun a _ => (Nmv_sum h _ _).trans (Finset.sum_le_sum fun b _ => ?_)
  exact (Nmv_monomial h _ _).trans (h.mul _ _)

lemma Nmv_smul (h : IsSN ν) (c : ℚ) (p : MvPolynomial σ R) : Nmv ν (c • p) ≤ |(c : ℝ)| * Nmv ν p := by
  unfold Nmv
  rw [Finset.mul_sum]
  calc ∑ s ∈ (c • p).support, ν ((c • p).coeff s) ≤ ∑ s ∈ p.support, ν ((c • p).coeff s) :=
        Finset.sum_le_sum_of_subset_of_nonneg MvPolynomial.support_smul fun _ _ _ => h.nonneg _
    _ ≤ _ := Finset.sum_le_sum fun s _ => by rw [MvPolynomial.coeff_smul]; exact h.smul _ _

lemma Nmv_isSN (h : IsSN ν) : IsSN (Nmv ν : MvPolynomial σ R → ℝ) where
  nonneg := Nmv_nonneg h
  zero := Nmv_zero
  add := Nmv_add h
  mul := Nmv_mul h
  smul := Nmv_smul h
  one := by
    rw [← MvPolynomial.C_1, MvPolynomial.C_apply]
    exact (Nmv_monomial h _ _).trans h.one

lemma Nmv_C (h : IsSN ν) (a : R) : Nmv ν (MvPolynomial.C a : MvPolynomial σ R) ≤ ν a := by
  rw [MvPolynomial.C_apply]; exact Nmv_monomial h _ _

lemma Nmv_X (h : IsSN ν) (i : σ) : Nmv ν (MvPolynomial.X i : MvPolynomial σ R) ≤ 1 := by
  rw [MvPolynomial.X]; exact (Nmv_monomial h _ _).trans h.one

end

/-- the `ℓ¹` norm of a `Polynomial` -/
noncomputable def Npo {R : Type*} [CommRing R] (ν : R → ℝ) (p : Polynomial R) : ℝ :=
  ∑ s ∈ p.support, ν (p.coeff s)

section
variable {R : Type*} [CommRing R] [Algebra ℚ R] {ν : R → ℝ}

lemma Npo_coeff (h : IsSN ν) (p : Polynomial R) (s : ℕ) : ν (p.coeff s) ≤ Npo ν p := by
  by_cases hs : s ∈ p.support
  · exact Finset.single_le_sum (f := fun t => ν (p.coeff t)) (fun t _ => h.nonneg _) hs
  · rw [Polynomial.notMem_support_iff.mp hs, h.zero]
    exact Finset.sum_nonneg fun _ _ => h.nonneg _

lemma Npo_add (h : IsSN ν) (p q : Polynomial R) : Npo ν (p + q) ≤ Npo ν p + Npo ν q := by
  classical
  unfold Npo
  have e1 : ∑ s ∈ p.support, ν (p.coeff s) = ∑ s ∈ p.support ∪ q.support, ν (p.coeff s) :=
    Finset.sum_subset Finset.subset_union_left fun s _ hs => by
      rw [Polynomial.notMem_support_iff.mp hs, h.zero]
  have e2 : ∑ s ∈ q.support, ν (q.coeff s) = ∑ s ∈ p.support ∪ q.support, ν (q.coeff s) :=
    Finset.sum_subset Finset.subset_union_right fun s _ hs => by
      rw [Polynomial.notMem_support_iff.mp hs, h.zero]
  rw [e1, e2, ← Finset.sum_add_distrib]
  calc ∑ s ∈ (p + q).support, ν ((p + q).coeff s)
      ≤ ∑ s ∈ p.support ∪ q.support, ν ((p + q).coeff s) :=
        Finset.sum_le_sum_of_subset_of_nonneg Polynomial.support_add fun _ _ _ => h.nonneg _
    _ ≤ _ := Finset.sum_le_sum fun s _ => by rw [Polynomial.coeff_add]; exact h.add _ _

lemma Npo_sum (h : IsSN ν) {ι : Type*} (s : Finset ι) (f : ι → Polynomial R) :
    Npo ν (∑ i ∈ s, f i) ≤ ∑ i ∈ s, Npo ν (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [Npo]
  | insert a s ha ih =>
    rw [Finset.sum_insert ha, Finset.sum_insert ha]
    exact (Npo_add h _ _).trans (by linarith)

lemma Npo_monomial (h : IsSN ν) (s : ℕ) (a : R) : Npo ν (Polynomial.monomial s a) ≤ ν a := by
  classical
  unfold Npo
  by_cases ha : a = 0
  · simp [ha, h.zero]
  · rw [Polynomial.support_monomial s ha]; simp

lemma Npo_mul (h : IsSN ν) (p q : Polynomial R) : Npo ν (p * q) ≤ Npo ν p * Npo ν q := by
  classical
  conv_lhs => rw [p.as_sum_support, q.as_sum_support, Finset.sum_mul_sum]
  simp only [Polynomial.monomial_mul_monomial]
  refine (Npo_sum h _ _).trans ?_
  unfold Npo
  rw [Finset.sum_mul_sum]
  refine Finset.sum_le_sum fun a _ => (Npo_sum h _ _).trans (Finset.sum_le_sum fun b _ => ?_)
  exact (Npo_monomial h _ _).trans (h.mul _ _)

lemma Npo_isSN (h : IsSN ν) : IsSN (Npo ν : Polynomial R → ℝ) where
  nonneg := fun p => Finset.sum_nonneg fun _ _ => h.nonneg _
  zero := by simp [Npo]
  add := Npo_add h
  mul := Npo_mul h
  smul := fun c p => by
    unfold Npo
    rw [Finset.mul_sum]
    calc ∑ s ∈ (c • p).support, ν ((c • p).coeff s) ≤ ∑ s ∈ p.support, ν ((c • p).coeff s) :=
          Finset.sum_le_sum_of_subset_of_nonneg (Polynomial.support_smul c p)
            fun _ _ _ => h.nonneg _
      _ ≤ _ := Finset.sum_le_sum fun s _ => by rw [Polynomial.coeff_smul]; exact h.smul _ _
  one := by
    rw [← Polynomial.monomial_zero_one]
    exact (Npo_monomial h _ _).trans h.one

lemma Npo_C (h : IsSN ν) (a : R) : Npo ν (Polynomial.C a) ≤ ν a := by
  rw [← Polynomial.monomial_zero_left]; exact Npo_monomial h _ _

lemma Npo_X (h : IsSN ν) : Npo ν (Polynomial.X : Polynomial R) ≤ 1 := by
  rw [← Polynomial.monomial_one_one_eq_X]; exact (Npo_monomial h _ _).trans h.one

end

/-- the absolute value on `ℚ` -/
noncomputable def nQ (q : ℚ) : ℝ := |(q : ℝ)|

lemma nQ_isSN : IsSN nQ where
  nonneg := fun _ => abs_nonneg _
  zero := by simp [nQ]
  add := fun a b => by unfold nQ; push_cast; exact abs_add_le _ _
  mul := fun a b => by unfold nQ; push_cast; rw [abs_mul]
  smul := fun c a => by unfold nQ; rw [smul_eq_mul]; push_cast; rw [abs_mul]
  one := by simp [nQ]


lemma IsSN.mul_le {R : Type*} [CommRing R] [Algebra ℚ R] {ν : R → ℝ} (h : IsSN ν) {a b : R}
    {A B : ℝ} (ha : ν a ≤ A) (hb : ν b ≤ B) : ν (a * b) ≤ A * B :=
  (h.mul a b).trans (mul_le_mul ha hb (h.nonneg _) ((h.nonneg _).trans ha))

lemma IsSN.smul_le {R : Type*} [CommRing R] [Algebra ℚ R] {ν : R → ℝ} (h : IsSN ν) (c : ℚ) {a : R}
    {A C : ℝ} (hc : |(c : ℝ)| ≤ C) (ha : ν a ≤ A) : ν (c • a) ≤ C * A :=
  (h.smul c a).trans (mul_le_mul hc ha (h.nonneg _) ((abs_nonneg _).trans hc))

lemma sum_pow_mono {x y : ℝ} (hx : 0 ≤ x) (hxy : x ≤ y) (S : Finset ℕ) :
    ∑ m ∈ S, x ^ m ≤ ∑ m ∈ S, y ^ m :=
  Finset.sum_le_sum fun m _ => pow_le_pow_left₀ hx hxy m

noncomputable abbrev nu3 : MvPolynomial (Fin 3) ℚ → ℝ := Nmv nQ
noncomputable abbrev nuE : Polynomial (MvPolynomial (Fin 3) ℚ) → ℝ := Npo nu3
lemma nu3_SN : IsSN nu3 := Nmv_isSN nQ_isSN
lemma nuE_SN : IsSN nuE := Npo_isSN nu3_SN

noncomputable def Wc (K : ℕ) : ℝ := ∑ r ∈ Finset.Icc 1 K, (2 : ℝ) ^ r
noncomputable def Sc (K : ℕ) : ℝ := ∑ m ∈ Finset.Icc 1 K, Wc K ^ m
noncomputable def Bc (K : ℕ) : ℝ := 2 + 2 * Sc K

lemma inv_fact_le (r : ℕ) : |((((r.factorial : ℚ)⁻¹ : ℚ)) : ℝ)| ≤ 1 := by
  push_cast
  rw [abs_of_nonneg (by positivity)]
  exact inv_le_one_of_one_le₀ (by exact_mod_cast r.factorial_pos)

lemma alt_le (m : ℕ) : |((((-1 : ℚ) ^ (m + 1) / (m : ℚ) : ℚ)) : ℝ)| ≤ 1 := by
  push_cast
  rw [abs_div, abs_pow, abs_neg, abs_one, one_pow]
  rcases Nat.eq_zero_or_pos m with h | h
  · subst h; simp
  · have : (1 : ℝ) ≤ m := by exact_mod_cast h
    rw [abs_of_pos (by linarith)]; exact (div_le_one (by linarith)).mpr this

lemma alt_le' (m : ℕ) : |((((-1 : ℚ) ^ (m - 1) / (m : ℚ) : ℚ)) : ℝ)| ≤ 1 := by
  push_cast
  rw [abs_div, abs_pow, abs_neg, abs_one, one_pow]
  rcases Nat.eq_zero_or_pos m with h | h
  · subst h; simp
  · have : (1 : ℝ) ≤ m := by exact_mod_cast h
    rw [abs_of_pos (by linarith)]; exact (div_le_one (by linarith)).mpr this

lemma Eser_bd (n K : ℕ) (hn : 1 ≤ n) : nuE (Eser n K) ≤ ((K : ℝ) + 1) * (Bc K * n) ^ K := by
  have h := nuE_SN
  have h3 := nu3_SN
  have hn' : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hWc0 : 0 ≤ Wc K := Finset.sum_nonneg fun _ _ => by positivity
  have hSc0 : 0 ≤ Sc K := Finset.sum_nonneg fun _ _ => by positivity
  have hX : ∀ i : Fin 3, nu3 (MvPolynomial.X i) ≤ 1 := fun i => Nmv_X nQ_isSN i
  have hl : nu3 (MvPolynomial.X 2 - MvPolynomial.X 1) ≤ 2 :=
    (h3.sub _ _).trans (by linarith [hX 2, hX 1])
  have hPX : nuE (Polynomial.X) ≤ 1 := Npo_X h3
  have hCX : nuE (Polynomial.C (MvPolynomial.X 2 - MvPolynomial.X 1) * Polynomial.X) ≤ 2 := by
    have := h.mul_le ((Npo_C h3 _).trans hl) hPX
    linarith
  have hw : nuE (Polynomial.C (MvPolynomial.X 0 : MvPolynomial (Fin 3) ℚ) *
      ∑ r ∈ Finset.Icc 1 K, ((r.factorial : ℚ)⁻¹) •
        (Polynomial.C (MvPolynomial.X 2 - MvPolynomial.X 1) * Polynomial.X) ^ r) ≤ Wc K := by
    have hs : nuE (∑ r ∈ Finset.Icc 1 K, ((r.factorial : ℚ)⁻¹) •
        (Polynomial.C (MvPolynomial.X 2 - MvPolynomial.X 1) * Polynomial.X) ^ r) ≤ Wc K := by
      refine (h.sum _ _).trans (Finset.sum_le_sum fun r _ => ?_)
      have := h.smul_le _ (inv_fact_le r) ((h.pow _ r).trans (pow_le_pow_left₀ (h.nonneg _) hCX r))
      linarith
    have := h.mul_le ((Npo_C h3 _).trans (hX 0)) hs
    linarith
  have hlog := (h.logT_le K _).trans (sum_pow_mono (h.nonneg _) hw _)
  have hA : nu3 ((-((n : ℚ) - 1)) • (MvPolynomial.X 1 : MvPolynomial (Fin 3) ℚ) +
      (n : ℚ) • MvPolynomial.X 2) ≤ 2 * n := by
    have e1 : |(((-((n : ℚ) - 1)) : ℚ) : ℝ)| ≤ n := by
      push_cast; rw [abs_of_nonpos (by linarith)]; linarith
    have e2 : |(((n : ℚ)) : ℝ)| ≤ n := by push_cast; rw [abs_of_nonneg (by linarith)]
    have := h3.add ((-((n : ℚ) - 1)) • (MvPolynomial.X 1 : MvPolynomial (Fin 3) ℚ))
      ((n : ℚ) • MvPolynomial.X 2)
    have a1 := h3.smul_le _ e1 (hX 1)
    have a2 := h3.smul_le _ e2 (hX 2)
    linarith
  have hCA := h.mul_le ((Npo_C h3 _).trans hA) hPX
  have e3 : |(((1 - 2 * (n : ℚ)) : ℚ) : ℝ)| ≤ 2 * n := by
    push_cast; rw [abs_of_nonpos (by linarith)]; linarith
  have hS := h.smul_le _ e3 hlog
  unfold Eser
  dsimp only
  have hL := (h.add _ _).trans (add_le_add hCA hS)
  have hL' : nuE (Polynomial.C ((-((n : ℚ) - 1)) • (MvPolynomial.X 1 : MvPolynomial (Fin 3) ℚ) +
      (n : ℚ) • MvPolynomial.X 2) * Polynomial.X + (1 - 2 * (n : ℚ)) • logT K
      (Polynomial.C (MvPolynomial.X 0 : MvPolynomial (Fin 3) ℚ) *
      ∑ r ∈ Finset.Icc 1 K, ((r.factorial : ℚ)⁻¹) •
        (Polynomial.C (MvPolynomial.X 2 - MvPolynomial.X 1) * Polynomial.X) ^ r)) ≤ Bc K * n := by
    refine hL.trans ?_
    unfold Bc Sc; nlinarith
  have hB1 : 1 ≤ Bc K * n := by unfold Bc; nlinarith
  exact (h.expT_le K _).trans ((sum_pow_mono (h.nonneg _) hL' _).trans (geom_le hB1 K))

lemma eF_bd (k n j i m : ℕ) (hn : 1 ≤ n) :
    |(eF k n j i m : ℝ)| ≤ ((k : ℝ) + 1) * (Bc k * n) ^ k := by
  have h1 : |(eF k n j i m : ℝ)| ≤ nu3 (Ecoef k n) :=
    Nmv_coeff nQ_isSN (Ecoef k n) (Finsupp.single 0 j + Finsupp.single 1 i + Finsupp.single 2 m)
  have h2 : nu3 (Ecoef k n) ≤ nuE (Eser n k) := Npo_coeff nu3_SN (Eser n k) k
  exact h1.trans (h2.trans (Eser_bd n k hn))

noncomputable abbrev nu1 : MvPolynomial ℕ ℚ → ℝ := Nmv nQ
noncomputable abbrev nu2 : S2 → ℝ := Nmv nu1
lemma nu1_SN : IsSN nu1 := Nmv_isSN nQ_isSN
lemma nu2_SN : IsSN nu2 := Nmv_isSN nu1_SN

lemma Hr_le (r P : ℕ) (hr : 1 ≤ r) {M : ℝ} (hP : (P : ℝ) ≤ M) (hM : 1 ≤ M) :
    0 ≤ V05.Hr r P ∧ V05.Hr r P ≤ 1 + Real.log M := by
  have h0 : 0 ≤ V05.Hr r P := Finset.sum_nonneg fun _ _ => by positivity
  refine ⟨h0, ?_⟩
  have h1 : V05.Hr r P ≤ V05.Hr 1 P := Finset.sum_le_sum fun l _ => by
    have hl : (1 : ℝ) ≤ (l : ℝ) + 1 := by have := l.cast_nonneg (α := ℝ); linarith
    rw [pow_one]
    exact one_div_le_one_div_of_le (by positivity) (le_self_pow₀ hl (by omega))
  have h2 : V05.Hr 1 P = (harmonic P : ℝ) := by
    unfold V05.Hr harmonic
    push_cast
    refine Finset.sum_congr rfl fun l _ => by ring
  have hlog : 0 ≤ Real.log M := Real.log_nonneg hM
  rcases Nat.eq_zero_or_pos P with hP0 | hP0
  · subst hP0; simp [V05.Hr] at h1 ⊢; linarith
  · have h3 := harmonic_le_one_add_log P
    have h4 : Real.log P ≤ Real.log M :=
      Real.log_le_log (by exact_mod_cast hP0) hP
    linarith

noncomputable def Dc (K : ℕ) : ℝ :=
  1 + ∑ r ∈ Finset.Icc 1 K, 3 * (2 : ℝ) ^ r + ∑ w ∈ Finset.Icc 2 K, 3 * (2 : ℝ) ^ w

lemma Dc_one (K : ℕ) : 1 ≤ Dc K := by
  unfold Dc
  have := Finset.sum_nonneg (s := Finset.Icc 1 K) fun r _ => (by positivity : (0:ℝ) ≤ 3 * 2 ^ r)
  have := Finset.sum_nonneg (s := Finset.Icc 2 K) fun r _ => (by positivity : (0:ℝ) ≤ 3 * 2 ^ r)
  linarith

lemma DeltaT_bd (K p q : ℕ) {M : ℝ} (hM : 1 ≤ M) (hpq : ((p + q - 1 : ℕ) : ℝ) ≤ M) :
    nu2 (DeltaT K p q) ≤ Dc K * (1 + Real.log M) := by
  have h := nu2_SN
  have hlog : 0 ≤ Real.log M := Real.log_nonneg hM
  set Hm := 1 + Real.log M with hHm
  have hHm1 : 1 ≤ Hm := by linarith
  have hx : ∀ i : Fin 2, nu2 (MvPolynomial.X i) ≤ 1 := fun i => Nmv_X nu1_SN i
  have hxx : nu2 (MvPolynomial.X 0 + MvPolynomial.X 1 : S2) ≤ 2 :=
    (h.add _ _).trans (by linarith [hx 0, hx 1])
  have hpow1 : ∀ (i : Fin 2) (r : ℕ), nu2 ((MvPolynomial.X i : S2) ^ r) ≤ 1 := fun i r =>
    (h.pow _ r).trans (pow_le_one₀ (h.nonneg _) (hx i))
  have hpow2 : ∀ r : ℕ, nu2 ((MvPolynomial.X 0 + MvPolynomial.X 1 : S2) ^ r) ≤ 2 ^ r := fun r =>
    (h.pow _ r).trans (pow_le_pow_left₀ (h.nonneg _) hxx r)
  have hH : ∀ r P : ℕ, 1 ≤ r → P ≤ p + q - 1 → |((V05.HrQ r P : ℚ) : ℝ)| ≤ Hm := by
    intro r P hr hP
    rw [V05.HrQ_cast]
    have hP' : (P : ℝ) ≤ M := le_trans (by exact_mod_cast hP) hpq
    obtain ⟨a, b⟩ := Hr_le r P hr hP' hM
    rw [abs_of_nonneg a]; exact b
  unfold DeltaT
  dsimp only
  refine (h.add _ _).trans ?_
  have hA : nu2 (∑ r ∈ Finset.Icc 1 K, ((-1 : ℚ) ^ (r - 1) / (r : ℚ)) •
      (V05.HrQ r (p - 1) • (MvPolynomial.X 0 : S2) ^ r + V05.HrQ r (q - 1) • MvPolynomial.X 1 ^ r
        - V05.HrQ r (p + q - 1) • (MvPolynomial.X 0 + MvPolynomial.X 1) ^ r)) ≤
      ∑ r ∈ Finset.Icc 1 K, 3 * (2 : ℝ) ^ r * Hm := by
    refine (h.sum _ _).trans (Finset.sum_le_sum fun r hr => ?_)
    have hr1 := (Finset.mem_Icc.mp hr).1
    have b1 := h.smul_le _ (hH r (p - 1) hr1 (by omega)) (hpow1 0 r)
    have b2 := h.smul_le _ (hH r (q - 1) hr1 (by omega)) (hpow1 1 r)
    have b3 := h.smul_le _ (hH r (p + q - 1) hr1 le_rfl) (hpow2 r)
    have b4 := (h.sub _ _).trans (add_le_add ((h.add _ _).trans (add_le_add b1 b2)) b3)
    have b5 := h.smul_le _ (alt_le' r) b4
    have : (1 : ℝ) ≤ 2 ^ r := one_le_pow₀ (by norm_num)
    nlinarith
  have hB : nu2 (∑ w ∈ Finset.Icc 2 K, ((-1 : ℚ) ^ (w + 1) / (w : ℚ)) •
      (MvPolynomial.C (MvPolynomial.X w) * ((MvPolynomial.X 0 + MvPolynomial.X 1 : S2) ^ w -
        MvPolynomial.X 0 ^ w - MvPolynomial.X 1 ^ w))) ≤ ∑ w ∈ Finset.Icc 2 K, 3 * (2 : ℝ) ^ w * Hm := by
    refine (h.sum _ _).trans (Finset.sum_le_sum fun w _ => ?_)
    have hC : nu2 (MvPolynomial.C (MvPolynomial.X w) : S2) ≤ 1 :=
      (Nmv_C nu1_SN _).trans (Nmv_X nQ_isSN w)
    have b4 := (h.sub _ _).trans (add_le_add ((h.sub _ _).trans (add_le_add (hpow2 w) (hpow1 0 w)))
      (hpow1 1 w))
    have b5 := h.smul_le _ (alt_le w) (h.mul_le hC b4)
    have : (1 : ℝ) ≤ 2 ^ w := one_le_pow₀ (by norm_num)
    nlinarith
  have e : Dc K * Hm = Hm + (∑ r ∈ Finset.Icc 1 K, 3 * (2 : ℝ) ^ r * Hm) +
      ∑ w ∈ Finset.Icc 2 K, 3 * (2 : ℝ) ^ w * Hm := by
    unfold Dc; rw [add_mul, add_mul, Finset.sum_mul, Finset.sum_mul]; ring
  rw [e]; linarith

lemma DsharpF_bd (i m p q : ℕ) {M : ℝ} (hM : 1 ≤ M) (hpq : ((p + q - 1 : ℕ) : ℝ) ≤ M) :
    nu1 (DsharpF i m p q) ≤ ((i.factorial * m.factorial : ℕ) : ℝ) *
      (((i + m : ℕ) : ℝ) + 1) * (Dc (i + m) * (1 + Real.log M)) ^ (i + m) := by
  have hD := DeltaT_bd (i + m) p q hM hpq
  have hlog : 0 ≤ Real.log M := Real.log_nonneg hM
  have hDc : 1 ≤ Dc (i + m) := by
    unfold Dc
    have := Finset.sum_nonneg (s := Finset.Icc 1 (i + m)) fun r _ => (by positivity : (0:ℝ) ≤ 3 * 2 ^ r)
    have := Finset.sum_nonneg (s := Finset.Icc 2 (i + m)) fun r _ => (by positivity : (0:ℝ) ≤ 3 * 2 ^ r)
    linarith
  have hB1 : 1 ≤ Dc (i + m) * (1 + Real.log M) := by nlinarith
  have hE := (nu2_SN.expT_le (i + m) _).trans
    ((sum_pow_mono (nu2_SN.nonneg _) hD _).trans (geom_le hB1 (i + m)))
  have hc := Nmv_coeff nu1_SN (expT (i + m) (DeltaT (i + m) p q))
    (Finsupp.single 0 i + Finsupp.single 1 m)
  unfold DsharpF
  have habs : |((((i.factorial * m.factorial : ℕ) : ℚ)) : ℝ)| = ((i.factorial * m.factorial : ℕ) : ℝ) := by
    push_cast; rw [abs_of_nonneg (by positivity)]
  refine (nu1_SN.smul_le _ habs.le (hc.trans hE)).trans ?_
  push_cast
  ring_nf
  exact le_rfl

end Q04

theorem Q04.proof_6 : Q04.stmt_6 := by
  intro k
  refine ⟨((k : ℝ) + 1) ^ 2 * (((k : ℝ) + 1) * Q04.Bc k ^ k) *
    (((k.factorial * k.factorial : ℕ) : ℝ) * (((k : ℝ) + 1) * Q04.Dc k ^ k)), fun n hn s => ?_⟩
  have hn' : (1 : ℝ) ≤ n := by exact_mod_cast hn
  set M : ℝ := 2 * (n : ℝ) + k with hMdef
  have hM : 1 ≤ M := by have := k.cast_nonneg (α := ℝ); linarith
  have hlog : 0 ≤ Real.log M := Real.log_nonneg hM
  set Hm := 1 + Real.log M with hHm
  have h1 := Q04.nu1_SN
  have hc := Q04.Nmv_coeff Q04.nQ_isSN (Q04.Pf k n) s
  unfold Q04.nQ at hc
  refine hc.trans ?_
  have e0 : |((-(1 / ((n : ℚ) * (2 * (n : ℚ) - 1))) : ℚ) : ℝ)| ≤ 1 := by
    have hp : (0 : ℝ) < (n : ℝ) * (2 * (n : ℝ) - 1) := by nlinarith
    push_cast
    rw [abs_neg, abs_of_pos (div_pos one_pos hp)]
    rw [div_le_one hp]; nlinarith
  have hterm : ∀ j ∈ Finset.range (k + 1), ∀ i ∈ Finset.range (k + 1),
      Q04.nu1 ((Q04.eF k n j i (k - i) * V05.rhoQ j n) • Q04.DsharpF i (k - i) (n + j) n) ≤
        (((k : ℝ) + 1) * (Q04.Bc k * n) ^ k) *
          (((k.factorial * k.factorial : ℕ) : ℝ) * (((k : ℝ) + 1) * (Q04.Dc k * Hm) ^ k)) := by
    intro j hj i hi
    have hj' := Nat.lt_succ_iff.mp (Finset.mem_range.mp hj)
    have hi' := Nat.lt_succ_iff.mp (Finset.mem_range.mp hi)
    have hr : |(((Q04.eF k n j i (k - i) * V05.rhoQ j n) : ℚ) : ℝ)| ≤
        ((k : ℝ) + 1) * (Q04.Bc k * n) ^ k := by
      push_cast
      rw [abs_mul, Q04.rhoQ_cast]
      obtain ⟨r0, r1⟩ := V05.rho_bd j n
      rw [abs_of_nonneg r0]
      have := Q04.eF_bd k n j i (k - i) hn
      nlinarith [abs_nonneg ((Q04.eF k n j i (k - i) : ℚ) : ℝ)]
    have hpq : ((n + j + n - 1 : ℕ) : ℝ) ≤ M := by
      rw [hMdef]
      have : n + j + n - 1 ≤ 2 * n + k := by omega
      have := (Nat.cast_le (α := ℝ)).mpr this
      push_cast at this ⊢; linarith
    have hD := Q04.DsharpF_bd i (k - i) (n + j) n hM hpq
    rw [show i + (k - i) = k by omega] at hD
    have hfac : ((i.factorial * (k - i).factorial : ℕ) : ℝ) ≤ ((k.factorial * k.factorial : ℕ) : ℝ) := by
      exact_mod_cast Nat.mul_le_mul (Nat.factorial_le hi') (Nat.factorial_le (Nat.sub_le k i))
    have hD' : Q04.nu1 (Q04.DsharpF i (k - i) (n + j) n) ≤
        ((k.factorial * k.factorial : ℕ) : ℝ) * (((k : ℝ) + 1) * (Q04.Dc k * Hm) ^ k) := by
      refine hD.trans ?_
      rw [mul_assoc]
      have := Q04.Dc_one k
      exact mul_le_mul_of_nonneg_right hfac (by positivity)
    exact h1.smul_le _ hr hD'
  unfold Q04.Pf
  refine (h1.smul_le _ e0 ((h1.sum _ _).trans (Finset.sum_le_sum fun j hj =>
    (h1.sum _ _).trans (Finset.sum_le_sum fun i hi => hterm j hj i hi)))).trans ?_
  simp only [Finset.sum_const, Finset.card_range, nsmul_eq_mul]
  rw [mul_pow, mul_pow]
  push_cast
  ring_nf
  exact le_rfl

theorem Q04.proof_7 : Q04.stmt_7 := by
  have z1 : ∀ n, (Q04.Pf 4 n).coeff (Finsupp.single 5 1) = 0 := by
    intro n; by_contra hne
    have h := (Q04.Pf_weight 4 n _ (MvPolynomial.mem_support_iff.mpr hne)).2
    rw [Q04.wsum_single] at h; omega
  have z2 : ∀ n, (Q04.Pf 4 n).coeff (Finsupp.single 2 1 + Finsupp.single 3 1) = 0 := by
    intro n; by_contra hne
    have h := (Q04.Pf_weight 4 n _ (MvPolynomial.mem_support_iff.mpr hne)).2
    rw [Finsupp.sum_add_index' (fun _ => by simp) (fun _ _ _ => by ring), Q04.wsum_single,
      Q04.wsum_single] at h
    omega
  refine ⟨fun n => (Q04.Pf 4 n).coeff 0, fun n => (Q04.Pf 4 n).coeff (Finsupp.single 2 1),
    fun n => (Q04.Pf 4 n).coeff (Finsupp.single 3 1), fun n => (Q04.Pf 4 n).coeff (Finsupp.single 4 1),
    fun n => (Q04.Pf 4 n).coeff (Finsupp.single 2 2), fun α hα0 hα1 => ?_⟩
  have h := (Q04.proof_5.2 α hα0 hα1 4).2
  refine h.congr_fun fun n => ?_
  rw [Q04.aeval_S5 _ (fun s hs => ⟨(Q04.Pf_weight 4 (n + 1) s hs).1,
    (Q04.Pf_weight 4 (n + 1) s hs).2.trans (by norm_num)⟩), z1, z2]
  push_cast; ring
theorem Q04.proof_8 : Q04.stmt_8 := by
  refine ⟨fun n => (Q04.Pf 5 n).coeff 0, fun n => (Q04.Pf 5 n).coeff (Finsupp.single 2 1),
    fun n => (Q04.Pf 5 n).coeff (Finsupp.single 3 1), fun n => (Q04.Pf 5 n).coeff (Finsupp.single 4 1),
    fun n => (Q04.Pf 5 n).coeff (Finsupp.single 5 1), fun n => (Q04.Pf 5 n).coeff (Finsupp.single 2 2),
    fun n => (Q04.Pf 5 n).coeff (Finsupp.single 2 1 + Finsupp.single 3 1), fun α hα0 hα1 => ?_⟩
  have h := (Q04.proof_5.2 α hα0 hα1 5).2
  refine h.congr_fun fun n => ?_
  rw [Q04.aeval_S5 _ (Q04.Pf_weight 5 (n + 1))]
namespace Q04
open Polynomial

/-- `∂_a + ∂_b` on `ℚ[y,a,b]` -/
noncomputable def Dd : Derivation ℚ (MvPolynomial (Fin 3) ℚ) (MvPolynomial (Fin 3) ℚ) :=
  MvPolynomial.pderiv 1 + MvPolynomial.pderiv 2

/-- `∂_a + ∂_b` applied to each coefficient of a polynomial in `u` -/
noncomputable def Dp (F : Polynomial (MvPolynomial (Fin 3) ℚ)) : Polynomial (MvPolynomial (Fin 3) ℚ) :=
  ∑ k ∈ Finset.range (F.natDegree + 1), Polynomial.C (Dd (F.coeff k)) * Polynomial.X ^ k

lemma Dp_coeff (F : Polynomial (MvPolynomial (Fin 3) ℚ)) (k : ℕ) : (Dp F).coeff k = Dd (F.coeff k) := by
  unfold Dp
  rw [Polynomial.finsetSum_coeff]
  simp only [Polynomial.coeff_C_mul_X_pow]
  rw [Finset.sum_ite_eq]
  split_ifs with hk
  · rfl
  · rw [Polynomial.coeff_eq_zero_of_natDegree_lt (by simp at hk; omega), map_zero]

lemma Dp_add (F G : Polynomial (MvPolynomial (Fin 3) ℚ)) : Dp (F + G) = Dp F + Dp G := by
  ext k; simp [Dp_coeff, map_add]

lemma Dp_smul (c : ℚ) (F : Polynomial (MvPolynomial (Fin 3) ℚ)) : Dp (c • F) = c • Dp F := by
  ext k; simp [Dp_coeff, Polynomial.coeff_smul]

lemma Dp_mul (F G : Polynomial (MvPolynomial (Fin 3) ℚ)) : Dp (F * G) = Dp F * G + F * Dp G := by
  refine Polynomial.ext fun k => ?_
  rw [Dp_coeff, Polynomial.coeff_add, Polynomial.coeff_mul, Polynomial.coeff_mul,
    Polynomial.coeff_mul, map_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl (fun x _ => ?_)
  rw [Derivation.leibniz, Dp_coeff, Dp_coeff, smul_eq_mul, smul_eq_mul]
  ring

lemma Dp_sum {ι : Type*} (S : Finset ι) (f : ι → Polynomial (MvPolynomial (Fin 3) ℚ)) :
    Dp (∑ i ∈ S, f i) = ∑ i ∈ S, Dp (f i) := by
  ext k; simp [Dp_coeff, Polynomial.finsetSum_coeff, map_sum]

lemma Dp_C (a : MvPolynomial (Fin 3) ℚ) : Dp (Polynomial.C a) = Polynomial.C (Dd a) := by
  ext k; rw [Dp_coeff, Polynomial.coeff_C, Polynomial.coeff_C]; split_ifs <;> simp

lemma Dp_X : Dp (Polynomial.X : Polynomial (MvPolynomial (Fin 3) ℚ)) = 0 := by
  ext k; rw [Dp_coeff, Polynomial.coeff_X]; split_ifs <;> simp

lemma Dp_one : Dp (1 : Polynomial (MvPolynomial (Fin 3) ℚ)) = 0 := by
  rw [← Polynomial.C_1, Dp_C, Derivation.map_one_eq_zero, map_zero]

lemma Dp_pow (F : Polynomial (MvPolynomial (Fin 3) ℚ)) (r : ℕ) :
    Dp (F ^ (r + 1)) = ((r : ℚ) + 1) • (F ^ r * Dp F) := by
  induction r with
  | zero => simp
  | succ r ih =>
    rw [pow_succ, Dp_mul, ih]
    rw [smul_mul_assoc, Algebra.smul_def, Algebra.smul_def]
    push_cast
    try simp only [map_add, map_one, map_natCast]
    ring

lemma Dp_expT (K : ℕ) (F : Polynomial (MvPolynomial (Fin 3) ℚ)) :
    Dp (expT (K + 1) F) = expT K F * Dp F := by
  unfold expT
  rw [Finset.sum_range_succ', Dp_add, Dp_sum, pow_zero, Dp_smul, Dp_one, smul_zero, add_zero,
    Finset.sum_mul]
  refine Finset.sum_congr rfl fun r _ => ?_
  rw [Dp_smul, Dp_pow, smul_smul, smul_mul_assoc]
  congr 1
  rw [Nat.factorial_succ]
  push_cast
  have : (r.factorial : ℚ) ≠ 0 := by positivity
  field_simp

lemma Dp_logT (K : ℕ) {w : Polynomial (MvPolynomial (Fin 3) ℚ)} (hw : Dp w = 0) :
    Dp (logT K w) = 0 := by
  unfold logT
  rw [Dp_sum]
  refine Finset.sum_eq_zero fun m hm => ?_
  obtain ⟨r, rfl⟩ : ∃ r, m = r + 1 := ⟨m - 1, by have := (Finset.mem_Icc.mp hm).1; omega⟩
  rw [Dp_smul, Dp_pow, hw, mul_zero, smul_zero, smul_zero]

lemma Dd_X (i : Fin 3) : Dd (MvPolynomial.X i) = if i = 0 then 0 else 1 := by
  unfold Dd
  rw [Derivation.add_apply, MvPolynomial.pderiv_X, MvPolynomial.pderiv_X]
  fin_cases i <;> simp

lemma Dd_C (c : ℚ) : Dd (MvPolynomial.C c) = 0 := by
  unfold Dd; rw [Derivation.add_apply]; simp

lemma Dd_ell : Dd (MvPolynomial.X 2 - MvPolynomial.X 1 : MvPolynomial (Fin 3) ℚ) = 0 := by
  rw [map_sub, Dd_X, Dd_X]; simp

/-- the exponent `L_K` of `Eser n K` -/
noncomputable def LK (n K : ℕ) : Polynomial (MvPolynomial (Fin 3) ℚ) :=
  Polynomial.C ((-((n : ℚ) - 1)) • (MvPolynomial.X 1 : MvPolynomial (Fin 3) ℚ) +
      (n : ℚ) • MvPolynomial.X 2) * Polynomial.X + (1 - 2 * (n : ℚ)) • logT K
    (Polynomial.C (MvPolynomial.X 0 : MvPolynomial (Fin 3) ℚ) *
      ∑ r ∈ Finset.Icc 1 K, ((r.factorial : ℚ)⁻¹) •
        (Polynomial.C (MvPolynomial.X 2 - MvPolynomial.X 1) * Polynomial.X) ^ r)

/-- the inner series `w_K` -/
noncomputable def wK (K : ℕ) : Polynomial (MvPolynomial (Fin 3) ℚ) :=
  Polynomial.C (MvPolynomial.X 0 : MvPolynomial (Fin 3) ℚ) *
      ∑ r ∈ Finset.Icc 1 K, ((r.factorial : ℚ)⁻¹) •
        (Polynomial.C (MvPolynomial.X 2 - MvPolynomial.X 1) * Polynomial.X) ^ r

lemma Eser_eq (n K : ℕ) : Eser n K = expT K (LK n K) := rfl

lemma LK_eq (n K : ℕ) : LK n K = Polynomial.C ((-((n : ℚ) - 1)) • (MvPolynomial.X 1 : MvPolynomial (Fin 3) ℚ) +
      (n : ℚ) • MvPolynomial.X 2) * Polynomial.X + (1 - 2 * (n : ℚ)) • logT K (wK K) := rfl

lemma Dp_wK (K : ℕ) : Dp (wK K) = 0 := by
  unfold wK
  rw [Dp_mul, Dp_C, Dd_X, if_pos rfl, map_zero, zero_mul, zero_add, Dp_sum]
  rw [Finset.sum_eq_zero fun r hr => ?_, mul_zero]
  obtain ⟨m, rfl⟩ : ∃ m, r = m + 1 := ⟨r - 1, by have := (Finset.mem_Icc.mp hr).1; omega⟩
  rw [Dp_smul, Dp_pow, Dp_mul, Dp_C, Dp_X, Dd_ell]
  simp

lemma Dp_LK (n K : ℕ) : Dp (LK n K) = Polynomial.X := by
  rw [LK_eq, Dp_add, Dp_smul, Dp_logT K (Dp_wK K), smul_zero, add_zero, Dp_mul, Dp_C, Dp_X,
    mul_zero, add_zero, map_add, Dd.map_smul, Dd.map_smul, Dd_X, Dd_X]
  simp only [Fin.isValue, one_ne_zero, if_false, show (2 : Fin 3) ≠ 0 by decide]
  rw [← add_smul, show (-((n : ℚ) - 1) + n) = 1 by ring, one_smul, Polynomial.C_1, one_mul]

lemma X_dvd_wK (K : ℕ) : Polynomial.X ∣ wK K := by
  unfold wK
  refine dvd_mul_of_dvd_right (Finset.dvd_sum fun r hr => ?_) _
  rw [Algebra.smul_def, mul_pow]
  exact dvd_mul_of_dvd_right (dvd_mul_of_dvd_right (dvd_pow_self _
    (by have := (Finset.mem_Icc.mp hr).1; omega)) _) _

lemma qsmul_dvd {a F : Polynomial (MvPolynomial (Fin 3) ℚ)} (c : ℚ) (h : a ∣ F) : a ∣ c • F := by
  rw [Algebra.smul_def]; exact dvd_mul_of_dvd_right h _

lemma wK_succ (k : ℕ) : Polynomial.X ^ (k + 1) ∣ wK (k + 1) - wK k := by
  unfold wK
  rw [Finset.sum_Icc_succ_top (by omega), mul_add, add_sub_cancel_left]
  refine dvd_mul_of_dvd_right (qsmul_dvd _ ?_) _
  rw [mul_pow]; exact dvd_mul_left _ _

lemma LK_succ (n k : ℕ) : Polynomial.X ^ (k + 1) ∣ LK n (k + 1) - LK n k := by
  rw [LK_eq, LK_eq, add_sub_add_left_eq_sub, ← smul_sub]
  refine qsmul_dvd _ ?_
  unfold logT
  rw [Finset.sum_Icc_succ_top (by omega), add_sub_right_comm, ← Finset.sum_sub_distrib]
  refine dvd_add (Finset.dvd_sum fun m _ => ?_) (qsmul_dvd _ (pow_dvd_pow_of_dvd (X_dvd_wK _) _))
  rw [← smul_sub]
  exact qsmul_dvd _ ((wK_succ k).trans (sub_dvd_pow_sub_pow _ _ _))

lemma Ecoef_trunc (n k : ℕ) : Ecoef k n = (expT k (LK n (k + 1))).coeff k := by
  have h : Polynomial.X ^ (k + 1) ∣ expT k (LK n (k + 1)) - expT k (LK n k) := by
    unfold expT
    rw [← Finset.sum_sub_distrib]
    refine Finset.dvd_sum fun r _ => ?_
    rw [← smul_sub]
    exact qsmul_dvd _ ((LK_succ n k).trans (sub_dvd_pow_sub_pow _ _ _))
  have := (Polynomial.X_pow_dvd_iff.mp h) k (by omega)
  rw [Polynomial.coeff_sub, sub_eq_zero] at this
  rw [Ecoef, Eser_eq, this]

lemma Dd_Ecoef (n k : ℕ) : Dd (Ecoef (k + 1) n) = Ecoef k n := by
  have h := congrArg (fun F => F.coeff (k + 1)) (Dp_expT k (LK n (k + 1)))
  simp only [Dp_coeff, Dp_LK, Polynomial.coeff_mul_X] at h
  rw [Ecoef, Eser_eq, h, Ecoef_trunc]

lemma eF_DE (k n j i m : ℕ) :
    ((i : ℚ) + 1) * eF (k + 1) n j (i + 1) m + ((m : ℚ) + 1) * eF (k + 1) n j i (m + 1) =
      eF k n j i m := by
  have h := congrArg (fun P => P.coeff (Finsupp.single 0 j + Finsupp.single 1 i + Finsupp.single 2 m))
    (Dd_Ecoef n k)
  rw [show eF k n j i m = (Ecoef k n).coeff
    (Finsupp.single 0 j + Finsupp.single 1 i + Finsupp.single 2 m) from rfl, ← h]
  unfold Dd
  rw [Derivation.add_apply, MvPolynomial.coeff_add, MvPolynomial.coeff_pderiv,
    MvPolynomial.coeff_pderiv]
  have e1 : Finsupp.single (0 : Fin 3) j + Finsupp.single 1 i + Finsupp.single 2 m +
      Finsupp.single 1 1 = Finsupp.single 0 j + Finsupp.single 1 (i + 1) + Finsupp.single 2 m := by
    ext v; fin_cases v <;> simp
  have e2 : Finsupp.single (0 : Fin 3) j + Finsupp.single 1 i + Finsupp.single 2 m +
      Finsupp.single 2 1 = Finsupp.single 0 j + Finsupp.single 1 i + Finsupp.single 2 (m + 1) := by
    ext v; fin_cases v <;> simp
  rw [e1, e2]
  simp [eF]
  ring

end Q04

namespace Q04

abbrev Q2 := MvPolynomial (Fin 2) ℚ

/-- `𝔷_w ↦ 0` on the coefficients -/
noncomputable def K0 : S2 →+* Q2 := MvPolynomial.map MvPolynomial.constantCoeff

/-- `[𝔷_w] P` as `∂_{𝔷_w} P` at `𝔷 = 0` -/
noncomputable def dw (w : ℕ) (P : MvPolynomial ℕ ℚ) : ℚ :=
  MvPolynomial.constantCoeff (MvPolynomial.pderiv w P)

lemma dw_eq (w : ℕ) (P : MvPolynomial ℕ ℚ) : dw w P = P.coeff (Finsupp.single w 1) := by
  unfold dw; rw [MvPolynomial.constantCoeff_eq, MvPolynomial.coeff_pderiv]; simp

lemma dw_mul (w : ℕ) (P Q : MvPolynomial ℕ ℚ) :
    dw w (P * Q) = dw w P * MvPolynomial.constantCoeff Q + MvPolynomial.constantCoeff P * dw w Q := by
  unfold dw; rw [Derivation.leibniz, map_add, smul_eq_mul, smul_eq_mul, map_mul, map_mul]; ring

lemma dw_add (w : ℕ) (P Q : MvPolynomial ℕ ℚ) : dw w (P + Q) = dw w P + dw w Q := by
  unfold dw; rw [map_add, map_add]

lemma dw_sum (w : ℕ) {ι : Type*} (S : Finset ι) (f : ι → MvPolynomial ℕ ℚ) :
    dw w (∑ i ∈ S, f i) = ∑ i ∈ S, dw w (f i) := by
  unfold dw; rw [map_sum, map_sum]

lemma dw_smul (w : ℕ) (c : ℚ) (P : MvPolynomial ℕ ℚ) : dw w (c • P) = c * dw w P := by
  rw [dw_eq, dw_eq, MvPolynomial.coeff_smul, smul_eq_mul]

lemma dw_zero (w : ℕ) : dw w 0 = 0 := by unfold dw; simp

/-- `[𝔷_w]` applied to each coefficient -/
noncomputable def psiW (w : ℕ) (F : S2) : Q2 :=
  ∑ t ∈ F.support, MvPolynomial.monomial t (dw w (F.coeff t))

lemma psiW_coeff (w : ℕ) (F : S2) (t : Fin 2 →₀ ℕ) : (psiW w F).coeff t = dw w (F.coeff t) := by
  classical
  unfold psiW
  rw [MvPolynomial.coeff_sum]
  simp only [MvPolynomial.coeff_monomial]
  rw [Finset.sum_ite_eq']
  split_ifs with h
  · rfl
  · rw [MvPolynomial.notMem_support_iff.mp h, dw_zero]

lemma K0_coeff (F : S2) (t : Fin 2 →₀ ℕ) : (K0 F).coeff t = MvPolynomial.constantCoeff (F.coeff t) := by
  unfold K0; rw [MvPolynomial.coeff_map]

lemma psiW_add (w : ℕ) (F G : S2) : psiW w (F + G) = psiW w F + psiW w G := by
  ext t; simp [psiW_coeff, dw_add]

lemma psiW_smul (w : ℕ) (c : ℚ) (F : S2) : psiW w (c • F) = c • psiW w F := by
  ext t; simp [psiW_coeff, dw_smul]

lemma psiW_sum (w : ℕ) {ι : Type*} (S : Finset ι) (f : ι → S2) :
    psiW w (∑ i ∈ S, f i) = ∑ i ∈ S, psiW w (f i) := by
  ext t; simp [psiW_coeff, MvPolynomial.coeff_sum, dw_sum]

lemma psiW_mul (w : ℕ) (F G : S2) : psiW w (F * G) = psiW w F * K0 G + K0 F * psiW w G := by
  ext t
  rw [psiW_coeff, MvPolynomial.coeff_mul, dw_sum, MvPolynomial.coeff_add, MvPolynomial.coeff_mul,
    MvPolynomial.coeff_mul, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun x _ => ?_
  rw [dw_mul, psiW_coeff, psiW_coeff, K0_coeff, K0_coeff]

lemma psiW_C (w : ℕ) (c : MvPolynomial ℕ ℚ) : psiW w (MvPolynomial.C c) = MvPolynomial.C (dw w c) := by
  ext t; rw [psiW_coeff, MvPolynomial.coeff_C, MvPolynomial.coeff_C]; split_ifs <;> simp [dw_zero]

lemma psiW_one (w : ℕ) : psiW w 1 = 0 := by
  rw [← MvPolynomial.C_1, psiW_C]; unfold dw; simp

lemma psiW_X (w : ℕ) (i : Fin 2) : psiW w (MvPolynomial.X i) = 0 := by
  ext t; rw [psiW_coeff, MvPolynomial.coeff_X']; split_ifs <;> simp [dw, dw_zero]

lemma psiW_pow (w : ℕ) (F : S2) (r : ℕ) :
    psiW w (F ^ (r + 1)) = ((r : ℚ) + 1) • ((K0 F) ^ r * psiW w F) := by
  induction r with
  | zero => simp
  | succ r ih =>
    rw [pow_succ, psiW_mul, ih, map_pow]
    rw [smul_mul_assoc, Algebra.smul_def, Algebra.smul_def]
    push_cast
    try simp only [map_add, map_one, map_natCast]
    ring

lemma psiW_Xpow (w : ℕ) (i : Fin 2) (r : ℕ) : psiW w ((MvPolynomial.X i : S2) ^ r) = 0 := by
  rcases r with _ | r
  · rw [pow_zero, psiW_one]
  · rw [psiW_pow, psiW_X, mul_zero, smul_zero]

lemma psiW_XXpow (w : ℕ) (r : ℕ) : psiW w ((MvPolynomial.X 0 + MvPolynomial.X 1 : S2) ^ r) = 0 := by
  rcases r with _ | r
  · rw [pow_zero, psiW_one]
  · rw [psiW_pow, psiW_add, psiW_X, psiW_X, add_zero, mul_zero, smul_zero]

lemma psiW_expT (w K : ℕ) (F : S2) : psiW w (expT (K + 1) F) = expT K (K0 F) * psiW w F := by
  unfold expT
  rw [Finset.sum_range_succ', psiW_add, psiW_sum, pow_zero, psiW_smul, psiW_one, smul_zero, add_zero,
    Finset.sum_mul]
  refine Finset.sum_congr rfl fun r _ => ?_
  rw [psiW_smul, psiW_pow, smul_smul, smul_mul_assoc]
  congr 1
  rw [Nat.factorial_succ]
  push_cast
  have : (r.factorial : ℚ) ≠ 0 := by positivity
  field_simp

/-- the harmonic part of `Δ_⋆` over `ℚ` -/
noncomputable def Hb (N p q : ℕ) : Q2 :=
  ∑ r ∈ Finset.Icc 1 N, ((-1 : ℚ) ^ (r - 1) / (r : ℚ)) •
      (V05.HrQ r (p - 1) • (MvPolynomial.X 0 : Q2) ^ r + V05.HrQ r (q - 1) • MvPolynomial.X 1 ^ r
        - V05.HrQ r (p + q - 1) • (MvPolynomial.X 0 + MvPolynomial.X 1) ^ r)

/-- `𝖯_w(x,x')` -/
noncomputable def Qb (w : ℕ) : Q2 :=
  ((-1 : ℚ) ^ (w + 1) / (w : ℚ)) •
    ((MvPolynomial.X 0 + MvPolynomial.X 1 : Q2) ^ w - MvPolynomial.X 0 ^ w - MvPolynomial.X 1 ^ w)

lemma K0_X (i : Fin 2) : K0 (MvPolynomial.X i) = MvPolynomial.X i := by
  unfold K0; rw [MvPolynomial.map_X]

lemma K0_CX (w : ℕ) : K0 (MvPolynomial.C (MvPolynomial.X w)) = 0 := by
  unfold K0; rw [MvPolynomial.map_C, MvPolynomial.constantCoeff_X, map_zero]

lemma K0_Delta (N p q : ℕ) : K0 (DeltaT N p q) = Hb N p q := by
  unfold DeltaT Hb
  simp only [map_add, map_sub, map_sum, map_rat_smul, map_pow, map_mul, K0_X, K0_CX, zero_mul,
    smul_zero, Finset.sum_const_zero, add_zero]

lemma psiW_sub (w : ℕ) (F G : S2) : psiW w (F - G) = psiW w F - psiW w G := by
  ext t; simp [psiW_coeff, dw_eq]

lemma psiW_Delta (N p q w : ℕ) (hw : 2 ≤ w) (hwN : w ≤ N) : psiW w (DeltaT N p q) = Qb w := by
  unfold DeltaT
  simp only [psiW_add, psiW_sub, psiW_sum, psiW_smul, psiW_Xpow, psiW_XXpow, psiW_mul, psiW_C,
    K0_CX, zero_mul, add_zero, sub_zero, smul_zero, Finset.sum_const_zero, zero_add, dw_eq,
    map_sub, map_pow, map_add, K0_X, mul_zero]
  rw [Finset.sum_eq_single w]
  · rw [MvPolynomial.coeff_X, if_pos rfl, MvPolynomial.C_1, one_mul]; rfl
  · intro v _ hv
    rw [MvPolynomial.coeff_X', if_neg, MvPolynomial.C_0, zero_mul, smul_zero]
    intro h; apply hv
    have := DFunLike.congr_fun h v
    simp at this; by_contra hc; simp [Finsupp.single_apply, hc] at this
  · intro h; exfalso; exact h (Finset.mem_Icc.mpr ⟨hw, hwN⟩)

end Q04

namespace Q04

lemma Xw_D (w i m p q : ℕ) (hw : 2 ≤ w) (hN : w ≤ i + m) :
    (DsharpF i m p q).coeff (Finsupp.single w 1) =
      ((i.factorial * m.factorial : ℕ) : ℚ) *
        (expT (i + m - 1) (Hb (i + m) p q) * Qb w).coeff (Finsupp.single 0 i + Finsupp.single 1 m) := by
  unfold DsharpF
  rw [MvPolynomial.coeff_smul, smul_eq_mul, ← dw_eq, ← psiW_coeff]
  obtain ⟨N', hN'⟩ : ∃ N', i + m = N' + 1 := ⟨i + m - 1, by omega⟩
  rw [hN', psiW_expT, K0_Delta, psiW_Delta _ _ _ _ hw (by omega), Nat.add_sub_cancel]

/-- every monomial has total degree `≥ r` -/
def Ord (r : ℕ) (F : Q2) : Prop := ∀ t : Fin 2 →₀ ℕ, F.coeff t ≠ 0 → r ≤ t 0 + t 1

lemma Ord_zero' (F : Q2) : Ord 0 F := fun _ _ => Nat.zero_le _

lemma Ord_add {r : ℕ} {F G : Q2} (hF : Ord r F) (hG : Ord r G) : Ord r (F + G) := by
  intro t h
  rw [MvPolynomial.coeff_add] at h
  by_cases h1 : F.coeff t = 0
  · rw [h1, zero_add] at h; exact hG t h
  · exact hF t h1

lemma Ord_smul {r : ℕ} (c : ℚ) {F : Q2} (hF : Ord r F) : Ord r (c • F) := by
  intro t h
  rw [MvPolynomial.coeff_smul] at h
  exact hF t (right_ne_zero_of_mul h)

lemma Ord_sub {r : ℕ} {F G : Q2} (hF : Ord r F) (hG : Ord r G) : Ord r (F - G) := by
  rw [sub_eq_add_neg, ← neg_one_smul ℚ G]; exact Ord_add hF (Ord_smul _ hG)

lemma Ord_sum {ι : Type*} (S : Finset ι) {r : ℕ} (f : ι → Q2) (h : ∀ i ∈ S, Ord r (f i)) :
    Ord r (∑ i ∈ S, f i) := by
  classical
  induction S using Finset.induction_on with
  | empty => intro t ht; simp at ht
  | insert a S ha ih =>
    rw [Finset.sum_insert ha]
    exact Ord_add (h a (Finset.mem_insert_self _ _)) (ih fun i hi => h i (Finset.mem_insert_of_mem hi))

lemma Ord_mul {a b : ℕ} {F G : Q2} (hF : Ord a F) (hG : Ord b G) : Ord (a + b) (F * G) := by
  intro t h
  rw [MvPolynomial.coeff_mul] at h
  obtain ⟨⟨u, v⟩, huv, hne⟩ := Finset.exists_ne_zero_of_sum_ne_zero h
  rw [Finset.HasAntidiagonal.mem_antidiagonal] at huv
  have h1 := hF u (left_ne_zero_of_mul hne)
  have h2 := hG v (right_ne_zero_of_mul hne)
  subst huv
  simp only [Finsupp.add_apply]
  omega

lemma Ord_mono {a b : ℕ} {F : Q2} (hF : Ord a F) (h : b ≤ a) : Ord b F :=
  fun t ht => h.trans (hF t ht)

lemma Ord_pow {a : ℕ} {F : Q2} (hF : Ord a F) (r : ℕ) : Ord (r * a) (F ^ r) := by
  induction r with
  | zero => simpa using Ord_zero' 1
  | succ r ih =>
    rw [pow_succ, show (r + 1) * a = r * a + a by ring]; exact Ord_mul ih hF

lemma Ord_X (i : Fin 2) : Ord 1 (MvPolynomial.X i : Q2) := by
  intro t h
  rw [MvPolynomial.coeff_X'] at h
  split_ifs at h with h0
  · subst h0; fin_cases i <;> simp
  · exact absurd rfl h

lemma Ord_pow_sub {a : ℕ} {F G : Q2} (h : Ord a (F - G)) (r : ℕ) : Ord a (F ^ r - G ^ r) := by
  obtain ⟨c, hc⟩ := sub_dvd_pow_sub_pow F G r
  rw [hc]; simpa using Ord_mul h (Ord_zero' c)

lemma Ord_expT_sub {a : ℕ} {F G : Q2} (h : Ord a (F - G)) (d : ℕ) :
    Ord a (expT d F - expT d G) := by
  unfold expT
  rw [← Finset.sum_sub_distrib]
  exact Ord_sum _ _ fun r _ => by rw [← smul_sub]; exact Ord_smul _ (Ord_pow_sub h r)

lemma Ord_expT_tail {H : Q2} (hH : Ord 1 H) {d K : ℕ} (hdK : d ≤ K) :
    Ord (d + 1) (expT K H - expT d H) := by
  unfold expT
  rw [Finset.range_eq_Ico, Finset.range_eq_Ico,
    ← Finset.sum_Ico_consecutive _ (Nat.zero_le (d + 1)) (by omega : d + 1 ≤ K + 1),
    add_sub_cancel_left]
  refine Ord_sum _ _ fun r hr => Ord_smul _ (Ord_mono (Ord_pow hH r) ?_)
  have := (Finset.mem_Ico.mp hr).1; omega

lemma Hb_term_Ord (r p q : ℕ) :
    Ord r (((-1 : ℚ) ^ (r - 1) / (r : ℚ)) •
      (V05.HrQ r (p - 1) • (MvPolynomial.X 0 : Q2) ^ r + V05.HrQ r (q - 1) • MvPolynomial.X 1 ^ r
        - V05.HrQ r (p + q - 1) • (MvPolynomial.X 0 + MvPolynomial.X 1) ^ r)) := by
  have h0 := Ord_pow (Ord_X 0) r
  have h1 := Ord_pow (Ord_X 1) r
  have h2 := Ord_pow (Ord_add (Ord_X 0) (Ord_X 1)) r
  simp only [mul_one] at h0 h1 h2
  exact Ord_smul _ (Ord_sub (Ord_add (Ord_smul _ h0) (Ord_smul _ h1)) (Ord_smul _ h2))

lemma Hb_Ord1 (N p q : ℕ) : Ord 1 (Hb N p q) :=
  Ord_sum _ _ fun r hr => Ord_mono (Hb_term_Ord r p q) (Finset.mem_Icc.mp hr).1

lemma Hb_tail {d N : ℕ} (p q : ℕ) (hdN : d ≤ N) : Ord (d + 1) (Hb N p q - Hb d p q) := by
  unfold Hb
  rw [← Finset.Ico_add_one_right_eq_Icc, ← Finset.Ico_add_one_right_eq_Icc,
    ← Finset.sum_Ico_consecutive _ (by omega : 1 ≤ d + 1) (by omega : d + 1 ≤ N + 1),
    add_sub_cancel_left]
  refine Ord_sum _ _ fun r hr => Ord_mono (Hb_term_Ord r p q) ?_
  have := (Finset.mem_Ico.mp hr).1; omega

lemma coeff_eq_of_Ord {d : ℕ} {F G : Q2} (h : Ord (d + 1) (F - G)) (t : Fin 2 →₀ ℕ)
    (ht : t 0 + t 1 ≤ d) : F.coeff t = G.coeff t := by
  by_contra hne
  have := h t (by rw [MvPolynomial.coeff_sub]; exact sub_ne_zero.mpr hne)
  omega

/-- `[x^I x'^M] exp(Δ^H_⋆)` truncated at its own degree -/
noncomputable def gH (p q I M : ℕ) : ℚ :=
  (expT (I + M) (Hb (I + M) p q)).coeff (Finsupp.single 0 I + Finsupp.single 1 M)

lemma Hb_trunc (K N p q I M : ℕ) (hK : I + M ≤ K) (hN : I + M ≤ N) :
    (expT K (Hb N p q)).coeff (Finsupp.single 0 I + Finsupp.single 1 M) = gH p q I M := by
  have h1 := Ord_expT_tail (Hb_Ord1 N p q) hK
  have h2 := Ord_expT_sub (Hb_tail p q hN) (I + M)
  have h3 := Ord_add h1 h2
  rw [sub_add_sub_cancel] at h3
  exact coeff_eq_of_Ord h3 _ (by simp)

lemma coeff_mul_mon (G : Q2) (a b i m : ℕ) :
    (G * MvPolynomial.monomial (Finsupp.single 0 a + Finsupp.single 1 b) 1).coeff
        (Finsupp.single 0 i + Finsupp.single 1 m) =
      if a ≤ i ∧ b ≤ m then G.coeff (Finsupp.single 0 (i - a) + Finsupp.single 1 (m - b)) else 0 := by
  rw [MvPolynomial.coeff_mul_monomial', mul_one]
  have hle : (Finsupp.single (0 : Fin 2) a + Finsupp.single 1 b ≤
      Finsupp.single 0 i + Finsupp.single 1 m) ↔ (a ≤ i ∧ b ≤ m) := by
    rw [Finsupp.le_def]
    constructor
    · intro h; exact ⟨by simpa using h 0, by simpa using h 1⟩
    · intro h j; fin_cases j <;> simp [h.1, h.2]
  by_cases h : a ≤ i ∧ b ≤ m
  · rw [if_pos (hle.mpr h), if_pos h]
    congr 1
    ext j; fin_cases j <;> simp
  · rw [if_neg (fun h' => h (hle.mp h')), if_neg h]

lemma mon_eq (a b : ℕ) : (MvPolynomial.X 0 : Q2) ^ a * MvPolynomial.X 1 ^ b =
    MvPolynomial.monomial (Finsupp.single 0 a + Finsupp.single 1 b) 1 := by
  rw [MvPolynomial.X_pow_eq_monomial, MvPolynomial.X_pow_eq_monomial, MvPolynomial.monomial_mul,
    one_mul]

lemma Qb2_eq : Qb 2 = -(MvPolynomial.monomial (Finsupp.single 0 1 + Finsupp.single 1 1) (1 : ℚ)) := by
  rw [← mon_eq, Qb, MvPolynomial.smul_eq_C_mul]
  have e : (MvPolynomial.X 0 + MvPolynomial.X 1 : Q2) ^ 2 - MvPolynomial.X 0 ^ 2 -
      MvPolynomial.X 1 ^ 2 = MvPolynomial.C 2 * (MvPolynomial.X 0 ^ 1 * MvPolynomial.X 1 ^ 1) := by
    rw [show (MvPolynomial.C 2 : Q2) = 2 from map_ofNat MvPolynomial.C 2]; ring
  rw [e, ← mul_assoc, ← MvPolynomial.C_mul]
  norm_num

lemma Qb3_eq : Qb 3 = MvPolynomial.monomial (Finsupp.single 0 2 + Finsupp.single 1 1) (1 : ℚ) +
    MvPolynomial.monomial (Finsupp.single 0 1 + Finsupp.single 1 2) 1 := by
  rw [← mon_eq, ← mon_eq, Qb, MvPolynomial.smul_eq_C_mul]
  have e : (MvPolynomial.X 0 + MvPolynomial.X 1 : Q2) ^ 3 - MvPolynomial.X 0 ^ 3 -
      MvPolynomial.X 1 ^ 3 = MvPolynomial.C 3 * (MvPolynomial.X 0 ^ 2 * MvPolynomial.X 1 ^ 1 +
        MvPolynomial.X 0 ^ 1 * MvPolynomial.X 1 ^ 2) := by
    rw [show (MvPolynomial.C 3 : Q2) = 3 from map_ofNat MvPolynomial.C 3]; ring
  rw [e, ← mul_assoc, ← MvPolynomial.C_mul]
  norm_num

end Q04

namespace Q04

lemma X2D_mid (K I p q : ℕ) (hI : I ≤ K) :
    (DsharpF (I + 1) (K + 1 - I) p q).coeff (Finsupp.single 2 1) =
      -((((I + 1).factorial * (K + 1 - I).factorial : ℕ) : ℚ) * gH p q I (K - I)) := by
  rw [Xw_D 2 (I + 1) (K + 1 - I) p q le_rfl (by omega), show I + 1 + (K + 1 - I) = K + 2 by omega,
    Qb2_eq, mul_neg, MvPolynomial.coeff_neg, coeff_mul_mon, if_pos (by omega),
    show I + 1 - 1 = I by omega, show K + 1 - I - 1 = K - I by omega,
    Hb_trunc (K + 2 - 1) (K + 2) p q I (K - I) (by omega) (by omega)]
  ring

lemma X2D_edge0 (K p q : ℕ) : (DsharpF 0 (K + 2) p q).coeff (Finsupp.single 2 1) = 0 := by
  rw [Xw_D 2 0 (K + 2) p q le_rfl (by omega), Qb2_eq, mul_neg, MvPolynomial.coeff_neg, coeff_mul_mon,
    if_neg (by omega)]
  simp

lemma X2D_edge1 (K p q : ℕ) : (DsharpF (K + 2) 0 p q).coeff (Finsupp.single 2 1) = 0 := by
  rw [Xw_D 2 (K + 2) 0 p q le_rfl (by omega), Qb2_eq, mul_neg, MvPolynomial.coeff_neg, coeff_mul_mon,
    if_neg (by omega)]
  simp

lemma X3D (K i p q : ℕ) (hi : i ≤ K + 3) :
    (DsharpF i (K + 3 - i) p q).coeff (Finsupp.single 3 1) =
      (((i.factorial * (K + 3 - i).factorial : ℕ) : ℚ)) *
        ((if 2 ≤ i ∧ 1 ≤ K + 3 - i then gH p q (i - 2) (K + 3 - i - 1) else 0) +
          (if 1 ≤ i ∧ 2 ≤ K + 3 - i then gH p q (i - 1) (K + 3 - i - 2) else 0)) := by
  rw [Xw_D 3 i (K + 3 - i) p q (by norm_num) (by omega), show i + (K + 3 - i) = K + 3 by omega,
    Qb3_eq, mul_add, MvPolynomial.coeff_add, coeff_mul_mon, coeff_mul_mon]
  congr 1
  congr 1
  · split_ifs with h
    · exact Hb_trunc _ _ _ _ _ _ (by omega) (by omega)
    · rfl
  · split_ifs with h
    · exact Hb_trunc _ _ _ _ _ _ (by omega) (by omega)
    · rfl

lemma Pf_coeff (k n : ℕ) (s : ℕ →₀ ℕ) : (Pf k n).coeff s =
    -(1 / ((n : ℚ) * (2 * (n : ℚ) - 1))) * ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
      (eF k n j i (k - i) * V05.rhoQ j n) * (DsharpF i (k - i) (n + j) n).coeff s := by
  unfold Pf
  simp only [MvPolynomial.coeff_smul, MvPolynomial.coeff_sum, smul_eq_mul]

lemma X2_inner (K n j : ℕ) :
    ∑ i ∈ Finset.range (K + 2 + 1), (eF (K + 2) n j i (K + 2 - i) * V05.rhoQ j n) *
        (DsharpF i (K + 2 - i) (n + j) n).coeff (Finsupp.single 2 1) =
      -(V05.rhoQ j n * ∑ I ∈ Finset.range (K + 1), eF (K + 2) n j (I + 1) (K + 1 - I) *
        ((((I + 1).factorial * (K + 1 - I).factorial : ℕ) : ℚ) * gH (n + j) n I (K - I))) := by
  rw [Finset.sum_range_succ', Finset.sum_range_succ]
  rw [show K + 2 - 0 = K + 2 by omega, X2D_edge0, show K + 1 + 1 = K + 2 by omega,
    Nat.sub_self, X2D_edge1, mul_zero, mul_zero, add_zero, add_zero, Finset.mul_sum,
    ← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun I hI => ?_
  have hI' := Nat.lt_succ_iff.mp (Finset.mem_range.mp hI)
  rw [show K + 2 - (I + 1) = K + 1 - I by omega, X2D_mid K I _ _ hI']
  ring

lemma X3_inner (K n j : ℕ) :
    ∑ i ∈ Finset.range (K + 3 + 1), (eF (K + 3) n j i (K + 3 - i) * V05.rhoQ j n) *
        (DsharpF i (K + 3 - i) (n + j) n).coeff (Finsupp.single 3 1) =
      V05.rhoQ j n * ∑ I ∈ Finset.range (K + 1),
        (eF (K + 3) n j (I + 2) (K + 1 - I) * (((I + 2).factorial * (K + 1 - I).factorial : ℕ) : ℚ) +
          eF (K + 3) n j (I + 1) (K + 2 - I) * (((I + 1).factorial * (K + 2 - I).factorial : ℕ) : ℚ)) *
            gH (n + j) n I (K - I) := by
  rw [Finset.sum_congr rfl fun i hi => by
    rw [X3D K i (n + j) n (Nat.lt_succ_iff.mp (Finset.mem_range.mp hi)), mul_add, mul_add]]
  rw [Finset.sum_add_distrib]
  -- first sum: `2 ≤ i ≤ K + 2`
  have hA : ∑ i ∈ Finset.range (K + 3 + 1), (eF (K + 3) n j i (K + 3 - i) * V05.rhoQ j n) *
      ((((i.factorial * (K + 3 - i).factorial : ℕ) : ℚ)) *
        (if 2 ≤ i ∧ 1 ≤ K + 3 - i then gH (n + j) n (i - 2) (K + 3 - i - 1) else 0)) =
      V05.rhoQ j n * ∑ I ∈ Finset.range (K + 1),
        eF (K + 3) n j (I + 2) (K + 1 - I) * (((I + 2).factorial * (K + 1 - I).factorial : ℕ) : ℚ) *
          gH (n + j) n I (K - I) := by
    rw [Finset.sum_range_succ', Finset.sum_range_succ', show K + 3 = K + 2 + 1 by ring,
      Finset.sum_range_succ]
    rw [if_neg (by omega), if_neg (by omega), if_neg (by omega)]
    simp only [mul_zero, add_zero, Finset.mul_sum]
    refine Finset.sum_congr rfl fun I hI => ?_
    have hI' := Nat.lt_succ_iff.mp (Finset.mem_range.mp hI)
    rw [if_pos (by omega), show I + 1 + 1 - 2 = I by omega, show K + 2 + 1 - (I + 1 + 1) = K + 1 - I by omega,
      show K + 1 - I - 1 = K - I by omega, show I + 1 + 1 = I + 2 by ring]
    ring
  -- second sum: `1 ≤ i ≤ K + 1`
  have hB : ∑ i ∈ Finset.range (K + 3 + 1), (eF (K + 3) n j i (K + 3 - i) * V05.rhoQ j n) *
      ((((i.factorial * (K + 3 - i).factorial : ℕ) : ℚ)) *
        (if 1 ≤ i ∧ 2 ≤ K + 3 - i then gH (n + j) n (i - 1) (K + 3 - i - 2) else 0)) =
      V05.rhoQ j n * ∑ I ∈ Finset.range (K + 1),
        eF (K + 3) n j (I + 1) (K + 2 - I) * (((I + 1).factorial * (K + 2 - I).factorial : ℕ) : ℚ) *
          gH (n + j) n I (K - I) := by
    rw [Finset.sum_range_succ', show K + 3 = K + 1 + 1 + 1 by ring, Finset.sum_range_succ,
      Finset.sum_range_succ]
    rw [if_neg (by omega), if_neg (by omega), if_neg (by omega)]
    simp only [mul_zero, add_zero, Finset.mul_sum]
    refine Finset.sum_congr rfl fun I hI => ?_
    have hI' := Nat.lt_succ_iff.mp (Finset.mem_range.mp hI)
    rw [if_pos (by omega), show I + 1 - 1 = I by omega, show K + 1 + 1 + 1 - (I + 1) = K + 2 - I by omega,
      show K + 2 - I - 2 = K - I by omega]
    ring
  rw [hA, hB, ← mul_add, ← Finset.sum_add_distrib]
  congr 1
  refine Finset.sum_congr rfl fun I _ => ?_
  ring

lemma inner_DE (K n j : ℕ) :
    ∑ I ∈ Finset.range (K + 1),
        (eF (K + 3) n j (I + 2) (K + 1 - I) * (((I + 2).factorial * (K + 1 - I).factorial : ℕ) : ℚ) +
          eF (K + 3) n j (I + 1) (K + 2 - I) * (((I + 1).factorial * (K + 2 - I).factorial : ℕ) : ℚ)) *
            gH (n + j) n I (K - I) =
      ∑ I ∈ Finset.range (K + 1), eF (K + 2) n j (I + 1) (K + 1 - I) *
        ((((I + 1).factorial * (K + 1 - I).factorial : ℕ) : ℚ) * gH (n + j) n I (K - I)) := by
  refine Finset.sum_congr rfl fun I hI => ?_
  have hI' := Nat.lt_succ_iff.mp (Finset.mem_range.mp hI)
  obtain ⟨M, rfl⟩ : ∃ M, K = I + M := ⟨K - I, by omega⟩
  have hde := eF_DE (I + M + 2) n j (I + 1) (M + 1)
  rw [show I + M + 2 + 1 = I + M + 3 by ring, show M + 1 + 1 = M + 2 by ring,
    show I + 1 + 1 = I + 2 by ring] at hde
  rw [show I + M + 1 - I = M + 1 by omega, show I + M + 2 - I = M + 2 by omega,
    ← hde, Nat.factorial_succ (I + 1), show M + 2 = (M + 1) + 1 by ring, Nat.factorial_succ (M + 1)]
  push_cast
  ring

theorem stmt9_formal (k n : ℕ) :
    (Pf (k + 1) n).coeff (Finsupp.single 3 1) = -(Pf k n).coeff (Finsupp.single 2 1) := by
  rcases Nat.lt_or_ge k 2 with hk | hk
  · have h1 : (Pf (k + 1) n).coeff (Finsupp.single 3 1) = 0 := by
      by_contra hne
      have h := (Pf_weight (k + 1) n _ (MvPolynomial.mem_support_iff.mpr hne)).2
      rw [wsum_single] at h; omega
    have h2 : (Pf k n).coeff (Finsupp.single 2 1) = 0 := by
      by_contra hne
      have h := (Pf_weight k n _ (MvPolynomial.mem_support_iff.mpr hne)).2
      rw [wsum_single] at h; omega
    rw [h1, h2, neg_zero]
  obtain ⟨K, rfl⟩ : ∃ K, k = K + 2 := ⟨k - 2, by omega⟩
  rw [Pf_coeff, Pf_coeff, show K + 2 + 1 = K + 3 by ring]
  simp only [X2_inner, X3_inner, inner_DE]
  rw [Finset.sum_range_succ (f := fun j => V05.rhoQ j n * _)]
  rw [Finset.sum_eq_zero (s := Finset.range (K + 1)) (f := fun I => eF (K + 2) n (K + 3) (I + 1) (K + 1 - I) * _)
    fun I _ => by
      have : eF (K + 2) n (K + 3) (I + 1) (K + 1 - I) = 0 := by
        by_contra hne; have := (eF_supp _ _ _ _ _ hne).1; omega
      rw [this, zero_mul]]
  rw [mul_zero, add_zero, show K + 3 = K + 2 + 1 by ring, Finset.sum_neg_distrib]
  ring

end Q04

theorem Q04.proof_9 : Q04.stmt_9 := by
  refine ⟨fun k n _ => Q04.stmt9_formal k n, fun k α hα0 hα1 => ?_⟩
  obtain ⟨hc0, hc1, -, -, -⟩ := D04.cs_facts hα0 hα1
  obtain ⟨C, hC⟩ := Q04.proof_6 k
  have hr0 : 0 ≤ D04.c α ^ 2 := by positivity
  have hr1 : D04.c α ^ 2 < 1 := by nlinarith
  have hsum : Summable (fun n : ℕ => (((Q04.Pf k (n + 1)).coeff (Finsupp.single 2 1) : ℚ) : ℝ)
      * D04.c α ^ (2 * (n + 1))) := by
    refine (V05.maj_summable k hr0 hr1 (|C| * ((k : ℝ) + 2) ^ k)).of_norm_bounded fun n => ?_
    have hm : (1 : ℝ) ≤ (n : ℝ) + 1 := by have := n.cast_nonneg (α := ℝ); linarith
    have hb := hC (n + 1) (by omega) (Finsupp.single 2 1)
    push_cast at hb
    have hL0 : 0 ≤ Real.log (2 * ((n : ℝ) + 1) + k) := Real.log_nonneg (by
      have := k.cast_nonneg (α := ℝ); linarith)
    have hL1 : 1 + Real.log (2 * ((n : ℝ) + 1) + k) ≤ ((k : ℝ) + 2) * ((n : ℝ) + 1) := by
      have := Real.log_le_sub_one_of_pos (show (0 : ℝ) < 2 * ((n : ℝ) + 1) + k by positivity)
      have := k.cast_nonneg (α := ℝ)
      nlinarith
    have hP : ((n : ℝ) + 1) ^ k * (1 + Real.log (2 * ((n : ℝ) + 1) + k)) ^ k ≤
        ((k : ℝ) + 2) ^ k * ((n : ℝ) + 1) ^ (2 * k) := by
      calc ((n : ℝ) + 1) ^ k * (1 + Real.log (2 * ((n : ℝ) + 1) + k)) ^ k
          ≤ ((n : ℝ) + 1) ^ k * (((k : ℝ) + 2) * ((n : ℝ) + 1)) ^ k :=
            mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (by linarith) hL1 k) (by positivity)
        _ = _ := by rw [mul_pow, two_mul, pow_add]; ring
    have hb2 : |(((Q04.Pf k (n + 1)).coeff (Finsupp.single 2 1) : ℚ) : ℝ)| ≤
        |C| * (((k : ℝ) + 2) ^ k * ((n : ℝ) + 1) ^ (2 * k)) := by
      refine hb.trans ?_
      rw [mul_assoc]
      exact mul_le_mul (le_abs_self C) hP (by positivity) (abs_nonneg C)
    rw [Real.norm_eq_abs, abs_mul, abs_pow, pow_mul, abs_of_pos hc0]
    calc _ ≤ |C| * (((k : ℝ) + 2) ^ k * ((n : ℝ) + 1) ^ (2 * k)) * (D04.c α ^ 2) ^ (n + 1) :=
          mul_le_mul_of_nonneg_right hb2 (by positivity)
      _ = _ := by ring
  refine ⟨_, hsum.hasSum, ?_⟩
  have h := hsum.hasSum.neg
  refine h.congr_fun fun n => ?_
  rw [Q04.stmt9_formal]
  push_cast
  ring

namespace Q04

lemma z2_eq : X05.z 2 = V05.z2 := by rw [X05.z_two, V05.z2]

lemma z3_eq : X05.z 3 = V05.z3 := by
  unfold X05.z V05.z3; norm_num

/-- the witness of `thm:flat`: `P_{k,n}` for `k ≠ 2, 3`, and the explicit forms of
`V05.proof_12` for `k = 2, 3` -/
noncomputable def muW (r r' q' : ℕ → ℚ) (k n : ℕ) : MvPolynomial ℕ ℚ :=
  if k = 2 then MvPolynomial.C (r n) + MvPolynomial.C (1 / (4 * (n : ℚ) ^ 2 - 1)) * MvPolynomial.X 2
  else if k = 3 then MvPolynomial.C (r' n) + MvPolynomial.C (q' n) * MvPolynomial.X 2
    - MvPolynomial.C (1 / (4 * (n : ℚ) ^ 2 - 1)) * MvPolynomial.X 3
  else Pf k n

lemma muW_weight (r r' q' : ℕ → ℚ) (k n : ℕ) : WeightLE k (muW r r' q' k n) := by
  unfold muW
  split_ifs with h2 h3
  · subst h2
    exact WL_add (WL_mono (WL_C _) (by norm_num)) (by simpa using WL_mul (WL_C _) (WL_X le_rfl))
  · subst h3
    rw [sub_eq_add_neg]
    refine WL_add (WL_add (WL_mono (WL_C _) (by norm_num)) ?_) (WL_neg ?_)
    · exact WL_mono (by simpa using WL_mul (WL_C _) (WL_X le_rfl)) (by norm_num)
    · simpa using WL_mul (WL_C _) (WL_X (le_refl 2 |>.trans (by norm_num : 2 ≤ 3)))
  · exact Pf_weight k n

end Q04

theorem Q04.proof_10 : Q04.stmt_10 := by
  obtain ⟨r, r', q', -, hser⟩ := V05.proof_12
  refine ⟨Q04.muW r r' q', Q04.muW_weight r r' q', fun n _ => ⟨r n, r' n, q' n, ?_, ?_⟩,
    fun α hα0 hα1 => ?_⟩
  · simp [Q04.muW]
  · simp [Q04.muW]
  obtain ⟨-, -, -, h2, h3, hr, Ψ, hq, hr'⟩ := hser α hα0 hα1
  have a2 : ∀ n : ℕ, MvPolynomial.aeval X05.z (Q04.muW r r' q' 2 (n + 1)) =
      (r (n + 1) : ℝ) + V05.z2 / (4 * ((n : ℝ) + 1) ^ 2 - 1) := by
    intro n
    simp [Q04.muW, Q04.z2_eq]
    ring
  have a3 : ∀ n : ℕ, MvPolynomial.aeval X05.z (Q04.muW r r' q' 3 (n + 1)) =
      (r' (n + 1) : ℝ) + (q' (n + 1) : ℝ) * V05.z2 - V05.z3 / (4 * ((n : ℝ) + 1) ^ 2 - 1) := by
    intro n
    simp [Q04.muW, Q04.z2_eq, Q04.z3_eq]
    ring
  refine ⟨fun k => ?_, ?_, Ψ, ?_, ?_⟩
  · by_cases hk2 : k = 2
    · subst hk2
      exact h2.congr_fun fun n => by rw [a2]
    · by_cases hk3 : k = 3
      · subst hk3
        exact h3.congr_fun fun n => by rw [a3]
      · have := (Q04.proof_5.2 α hα0 hα1 k).2
        refine this.congr_fun fun n => ?_
        simp [Q04.muW, hk2, hk3]
  · rw [Q04.z2_eq]
    exact hr.congr_fun fun n => by simp [Q04.muW]
  · have n1 : (0 : ℕ →₀ ℕ) ≠ Finsupp.single 2 1 := (Finsupp.single_ne_zero.mpr one_ne_zero).symm
    have n2 : Finsupp.single 3 1 ≠ Finsupp.single 2 1 := by
      intro h; have := DFunLike.congr_fun h 3; simp at this
    exact hq.congr_fun fun n => by simp [Q04.muW, MvPolynomial.coeff_X', n1, n2]
  · rw [Q04.z2_eq, Q04.z3_eq]
    exact hr'.congr_fun fun n => by simp [Q04.muW, MvPolynomial.coeff_X']

