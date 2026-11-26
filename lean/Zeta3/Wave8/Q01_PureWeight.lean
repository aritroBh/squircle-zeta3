import Mathlib
import Zeta3.Wave4.X01_RhombusDerivs
import Zeta3.Wave4.X05_DiamondA4A5
/-!
# Q01 — pure weight (`thm:pw`) to every order in `τ`

Target id: Q01.
Paper labels: `thm:pw` (paper/ieee/main.tex lines 1183–1191, proof sketch lines 1193–1195,
`τ = 1 - 1/p` from line 1181); detailed proof paper/ieee/proofs/sec_pureweight.tex
(Lemmas `lem:pw:pure` (lines 43–49), `lem:pw:logA` (70–80), the expansion of `g` (82–83),
`lem:pw:remainder` (85–94), `lem:pw:asymp` (98–111), `prop:pw:rhombus` (113–131), the table of
Section `sec:pw:series` (133–153), the rhombus check (155)).

Quoted verbatim.

main.tex 1184–1190 (`thm:pw`):

  To every finite order in $\tau$ (the full series diverges, by Theorem~\ref{thm:zero}),
  $\log\frac{P(p)}{4\sqrt2}=\sum_{k\ge2}w_k\,\tau^k$, where each $w_k$ is a polynomial with
  rational coefficients in $\zeta(2),\dots,\zeta(k)$, homogeneous of weight $k$. Hence the
  coefficient of $\tau^n$ in $P$ has pure weight $n$:
  \frac{P}{\sqrt2}=4+\zeta(2)\tau^2+\zeta(3)\tau^3-\frac{33}{16}\zeta(4)\tau^4\\
  -\Bigl(\frac34\zeta(5)+\frac74\zeta(2)\zeta(3)\Bigr)\tau^5+O(\tau^6).
  The same holds for every rhombus (with polynomials in $q^2$ as coefficients) and, after
  division by the rational factor $Y_N(\tau)=...$, ... for the surface area of
  $|x_1|^p+\dots+|x_N|^p=1$ in every dimension $N$.

main.tex 1194 (proof sketch):

  Put $s=\eta\sigma$ in Step~2; by Step~3, $\log A(\eta,\eta\sigma)=\sum_{n\ge2}
  \frac{(-1)^n\zeta(n)}{n}\bigl[(\frac12+\sigma)^n+(\frac12-\sigma)^n-1\bigr]\eta^n$:
  $\zeta(n)$ always arrives together with $\eta^n$, ...

sec_pureweight.tex 39 (meaning of weight):

  A real number is ``a polynomial with rational coefficients in $\zeta(2),\dots,\zeta(k)$,
  homogeneous of weight $k$'' if it is the value of such a $c$ at $Z_j=\zeta(j)$.

sec_pureweight.tex 44 (`lem:pw:pure`):

  Pure series form a subring of $\mathcal R$. If $f$ is pure with $c_0=1$, then $1/f$ and
  $\log f$ are pure. If $f$ is pure with $c_0=0$, then $\exp f$ is pure.

sec_pureweight.tex 75 (`lem:pw:logA`, consequence):

  Consequently $\Psi(\eta,s)=\sum_{a,b\ge0}\Psi_{ab}(\zeta)\,\eta^as^b$, where
  $\Psi_{ab}\in\QQ[Z_2,\dots,Z_{a+b}]$ is homogeneous of weight $a+b$, $\Psi_{00}=1$, and
  $\Psi_{ab}=0$ if $a+b=1$.

sec_pureweight.tex 83 (expansion of `g`):

  $g(z)=\sum_mg_m(q)z^m$ with $g_m\in\QQ[q]$ ... $g_m(-q)=(-1)^mg_m(q)$ ... The first ones
  are $g_0=1$, $g_1=\frac q2$, $g_2=\frac14-\frac{q^2}8$ and
  $g_4=-\frac1{96}+\frac{5q^2}{96}-\frac{5q^4}{128}$.

sec_pureweight.tex 87–89 (`lem:pw:remainder`):

  For $z\in\RR$ and every $M\ge1$ there is $C_M$ (depending on $q$ only through $|q|\le1$) with
  \Bigl|g(z)-\sum_{m<M}g_mz^m\Bigr|\le C_M|z|^M\eu^{|z|/2}.

sec_pureweight.tex 99–102 (`lem:pw:asymp`):

  For every $M\ge1$ and $|\eta|\le\frac14$,
  \frac{P_{A,B}(p)}{4r}=\sum_{m<M,\ m\ \mathrm{even}}g_m(q)\,\eta^m\,\partial_s^m\Psi(\eta,0)+O(|\eta|^M).

sec_pureweight.tex 114–120 (`prop:pw:rhombus`):

  There are $\hat c_N\in\QQ[q^2][Z_2,\dots,Z_N]$, homogeneous of weight $N$, with $\hat c_0=1$,
  $\hat c_1=0$, such that for every $K\ge0$, as $\tau\to0^+$,
  \frac{P_{A,B}(p)}{4r}=\sum_{N\le K}(-1)^N\hat c_N(\zeta)\,\tau^N+O(\tau^{K+1}),
  \log\frac{P_{A,B}(p)}{4r}=\sum_{k=2}^Kw_k(q)\tau^k+O(\tau^{K+1}),

sec_pureweight.tex 153 (the logarithmic coefficients):

  we get $w_2=\tilde c_2=\frac14Z_2$, $w_3=\tilde c_3=\frac14Z_3$,
  $w_4=\tilde c_4-\frac12\tilde c_2^{\,2}=-\frac14Z_2^2+\frac1{32}Z_4$ and
  $w_5=\tilde c_5-\tilde c_2\tilde c_3=...=-\frac12Z_2Z_3-\frac3{16}Z_5$.

sec_pureweight.tex 155 (rhombus check), 249 (checks):

  For general $q$ the same table gives $c_2=...=\frac{(1-q^2)}4Z_2$ and
  $c_3=...=-\frac{1-q^2}4Z_3$ ... For the rhombus: $w_2=(1-q^2)\zeta(2)/4$,
  $w_3=(1-q^2)\zeta(3)/4$.

Formalization choices.
* `ζ(n)` is `X05.z n = (riemannZeta n).re` (imported, not redefined). For `n ≥ 2` it is real
  and equals the series `∑_{k≥1} 1/k^n` (Mathlib `zeta_nat_eq_tsum_of_gt_one`).
* Weight (supplement line 39): `IsZetaPoly k x` says that `x` is the value at `Z_j = ζ(j)` of a
  rational polynomial in `Z_2, Z_3, …` all of whose monomials have weight `k` (`Z_j` has weight
  `j`). It is written as a finite list of monomials `(c, [n_1,…,n_r])` with `c ∈ ℚ`, every
  `n_i ≥ 2` and `n_1+…+n_r = k`, standing for `c·ζ(n_1)⋯ζ(n_r)` (so automatically `n_i ≤ k`;
  the empty list of indices is a rational constant, of weight `0`). `IsZetaPolyQ2 k f` is the
  same with coefficients in `ℚ[q²]`: `f(q) = Σ c · q^{2j} · ζ(n_1)⋯ζ(n_r)`, for every real `q`.
  Since `ζ(2j) ∈ ℚ π^{2j}`, weight is a property of the representation, not of the number; the
  statements only assert that some representation of the stated weight exists, as the paper does.
* `P` is the library perimeter `X01.P A B` (the paper's `eq:Pint` with semi-axes `A, B`), and for
  the squircle `X01.P 1 1`. It is the true Euclidean length of the curve by `X01.proof_3`
  (`X01.perimeter A B p = ENNReal.ofReal (X01.P A B p)`), so all statements are about the true
  perimeter. `r = X01.r A B = √(A²+B²)`, `q = X01.q A B = (A²-B²)/(A²+B²)`; for `A = B = 1`,
  `4r = 4√2` and `q = 0`. `X01.I A B η = P(1/(1+η))/(4r)`.
* `A(η,s)` (supplement: `Ψ(η,s)`) is `X05.A η s`, defined by the Beta integral.
* `τ = tau p = 1 - 1/p`. The expansions of `thm:pw` are stated for `p > 1` and `τ ≤ τ₀`
  (the paper's `τ → 0⁺`), with an explicit `C`: `|error| ≤ C τ^{K+1}`. The coefficients are one
  sequence `w : ℕ → …`, chosen before `K` ("the `w_k` are the same for every `K`").
* `P > 0` for `p > 0`, so `Real.log (P/(4r))` is a genuine logarithm.
* `g_m(q)` is defined as the Taylor coefficient
  `iteratedDeriv m (z ↦ √(cosh z + q sinh z)) 0 / m!`; stmt_2 records that this function is
  `C^n` at `0` for every `n` (so the coefficients are genuine). Likewise stmt_4 records the
  analyticity that makes the `s`- and `η`-derivatives of `A` genuine.
* `lem:pw:pure` is about formal series; stmt_1 states its content at the level of the numbers
  `IsZetaPoly`: closure under `+`, rational multiples, products (weights add), and the
  coefficient form of `log` and `exp` of a series `Σ c_n τ^n`, `c_0 = 1`, `Σ w_n τ^n`, `w_0 = 0`:
  `f = exp(w)` is equivalent (for formal series) to `f' = w' f`, i.e.
  `n c_n = Σ_{k=1}^{n} k w_k c_{n-k}` for `n ≥ 1`.

Already proved elsewhere and imported, not restated:
* `eq:pw:P` (Step 1, the exact formula) is `X05.proof_1`; its rhombus form `eq:pw:rhombus` is
  `X01.proof_4`.
* Step 2 (the Beta moments) is `X05.proof_2` / `X01.proof_5`.
* Step 3, the generating identity `log A(η,s) = Σ_{n≥2} (-1)^n ζ(n)/n [(η/2+s)^n + (η/2-s)^n
  - η^n]` for `|η|, |s| < 1/2`, is `X05.proof_3`; the paper's form with `s = ησ` is this identity
  at `s = ησ` (for `|η| < 1/2`, `|ησ| < 1/2`), since `(η/2 ± ησ)^n = (1/2 ± σ)^n η^n`.
* The displayed series of `P/√2` through `τ⁵` with `O(τ⁶)` is `X05.proof_7` (eq:Pint, two-sided)
  and `V01.proof_9` (true perimeter).
* The case `M = 6`, `q = 0` of `lem:pw:remainder`/`lem:pw:asymp` is `X05.proof_4`/`X05.proof_5`.

Statements.
* stmt_1: `lem:pw:pure` at the level of numbers (ring closure; `log` and `exp` in coefficient form).
* stmt_2: the coefficients `g_m(q)`: genuine Taylor coefficients, rational polynomials in `q` of
  parity `m`, and `g_0, g_1, g_2, g_4`.
* stmt_3: `lem:pw:remainder` for every `M ≥ 1`, uniformly in `|q| ≤ 1`.
* stmt_4: `lem:pw:logA`, consequence: the pairing of `ζ(n)` with total degree: the `η^a` Taylor
  coefficient of `∂_s^m A(η,0)` is of weight `a+m`; odd `s`-moments vanish.
* stmt_5: `lem:pw:asymp` for every `M ≥ 1` and every rhombus.
* stmt_6: `prop:pw:rhombus`, first expansion: `P_{A,B}/(4r)` to every order in `τ`, coefficients
  in `ℚ[q²][ζ]` of pure weight.
* stmt_7: `prop:pw:rhombus`, second expansion = third sentence of `thm:pw` for the rhombus:
  `log(P_{A,B}/(4r))` to every order, `w_k ∈ ℚ[q²][ζ]` of weight `k`, `w_2, w_3` explicit.
* stmt_8: the "Hence" of `thm:pw`: the coefficient of `τ^n` in `P/√2` is of pure weight `n`,
  to every order, with constant term `4`.
* stmt_9: `w_2, …, w_5` (supplement line 153) for the squircle, as an explicit `O(τ⁶)` bound.
* stmt_10: `thm:pw`, first sentence, to every order, for the squircle (the final claim).

NOT COVERED: the parenthetical "the full series diverges" (that is `thm:zero`, proved in Y06);
the `N`-dimensional surface-area part with the factor `Y_N` (excluded from this target); the
formal-series (indeterminate `Z_j`) level of `lem:pw:pure` and `lem:pw:logA` (only their
consequences for the numbers are stated); the statement `w_k ∈ ℚ[q²][ζ]` is stated, but not
the stronger formal identity `Σ w_k τ^k = log(Σ (-1)^N ĉ_N τ^N)` in `ℚ[q²][Z][[τ]]`; the
numerical checks of the supplement (quadrature, `τ^{12}`); the explicit values of `w_k` for
`k ≥ 6`.
-/
-- BEGIN STATEMENTS
namespace Q01

open Real

/-- `τ = 1 - 1/p` (main.tex 1181). -/
noncomputable def tau (p : ℝ) : ℝ := 1 - 1 / p

/-- Supplement line 39: `x` is the value at `Z_j = ζ(j)` of a rational polynomial in
`Z_2, Z_3, …` homogeneous of weight `k`; a monomial `(c, [n_1,…,n_r])` stands for
`c · ζ(n_1) ⋯ ζ(n_r)` with all `n_i ≥ 2` and `n_1 + … + n_r = k`. -/
def IsZetaPoly (k : ℕ) (x : ℝ) : Prop :=
  ∃ L : List (ℚ × List ℕ),
    (∀ e ∈ L, (∀ n ∈ e.2, 2 ≤ n) ∧ e.2.sum = k) ∧
    x = (L.map fun e => (e.1 : ℝ) * (e.2.map X05.z).prod).sum

/-- Supplement line 25 (`w_k ∈ ℚ[q²][ζ(2),…,ζ(k)]` of weight `k`): for every real `q`,
`f(q) = Σ c · q^{2j} · ζ(n_1) ⋯ ζ(n_r)` (finite sum, `c ∈ ℚ`, all `n_i ≥ 2`, `Σ n_i = k`). -/
def IsZetaPolyQ2 (k : ℕ) (f : ℝ → ℝ) : Prop :=
  ∃ L : List (ℚ × ℕ × List ℕ),
    (∀ e ∈ L, (∀ n ∈ e.2.2, 2 ≤ n) ∧ e.2.2.sum = k) ∧
    ∀ x : ℝ, f x = (L.map fun e => (e.1 : ℝ) * x ^ (2 * e.2.1) * (e.2.2.map X05.z).prod).sum

/-- `g_m(q)`, the Taylor coefficients at `z = 0` of `g(z) = √(cosh z + q sinh z)`
(supplement line 83). -/
noncomputable def gcoef (m : ℕ) (q : ℝ) : ℝ :=
  iteratedDeriv m (fun z : ℝ => Real.sqrt (Real.cosh z + q * Real.sinh z)) 0 / (m.factorial : ℝ)

/-- `lem:pw:pure` (supplement line 44) for the numbers: values of weight-homogeneous rational
zeta polynomials are closed under sums, rational multiples and products (weights add); and if
`Σ c_n τ^n = exp(Σ w_n τ^n)` coefficientwise (`c_0 = 1`, `w_0 = 0`,
`n c_n = Σ_{k=1}^n k w_k c_{n-k}` for `n ≥ 1`), then the `c_n` are pure iff the `w_n` are
(`log` and `exp` of pure series are pure). -/
def stmt_1 : Prop :=
  (∀ k : ℕ, ∀ x y : ℝ, IsZetaPoly k x → IsZetaPoly k y → IsZetaPoly k (x + y)) ∧
  (∀ k : ℕ, ∀ x : ℝ, ∀ c : ℚ, IsZetaPoly k x → IsZetaPoly k ((c : ℝ) * x)) ∧
  (∀ a b : ℕ, ∀ x y : ℝ, IsZetaPoly a x → IsZetaPoly b y → IsZetaPoly (a + b) (x * y)) ∧
  (∀ c w : ℕ → ℝ, c 0 = 1 → w 0 = 0 →
    (∀ n : ℕ, 1 ≤ n →
      (n : ℝ) * c n = ∑ k ∈ Finset.Icc 1 n, (k : ℝ) * w k * c (n - k)) →
    ((∀ n, IsZetaPoly n (c n)) → ∀ n, IsZetaPoly n (w n)) ∧
    ((∀ n, IsZetaPoly n (w n)) → ∀ n, IsZetaPoly n (c n)))

/-- Supplement line 83: `g(z) = √(cosh z + q sinh z)` is smooth at `z = 0`; its Taylor
coefficients `g_m` are polynomials in `q` with rational coefficients, with
`g_m(-q) = (-1)^m g_m(q)` (so `g_m ∈ ℚ[q²]` for even `m`); `g_0 = 1`, `g_1 = q/2`,
`g_2 = 1/4 - q²/8`, `g_4 = -1/96 + 5q²/96 - 5q⁴/128`. -/
def stmt_2 : Prop :=
  (∀ q : ℝ, ∀ n : ℕ,
    ContDiffAt ℝ n (fun z : ℝ => Real.sqrt (Real.cosh z + q * Real.sinh z)) 0) ∧
  (∀ m : ℕ, ∃ Q : Polynomial ℚ, ∀ q : ℝ, gcoef m q = Polynomial.aeval q Q) ∧
  (∀ m : ℕ, ∀ q : ℝ, gcoef m (-q) = (-1) ^ m * gcoef m q) ∧
  (∀ q : ℝ, gcoef 0 q = 1 ∧ gcoef 1 q = q / 2 ∧ gcoef 2 q = 1 / 4 - q ^ 2 / 8 ∧
    gcoef 4 q = -1 / 96 + 5 * q ^ 2 / 96 - 5 * q ^ 4 / 128)

/-- `lem:pw:remainder` (supplement lines 87–89): for every `M ≥ 1` there is `C_M`, depending on
`q` only through `|q| ≤ 1`, with `|g(z) - Σ_{m<M} g_m z^m| ≤ C_M |z|^M e^{|z|/2}` for real `z`. -/
def stmt_3 : Prop :=
  ∀ M : ℕ, 1 ≤ M → ∃ C : ℝ, ∀ q : ℝ, |q| ≤ 1 → ∀ z : ℝ,
    |Real.sqrt (Real.cosh z + q * Real.sinh z) - ∑ m ∈ Finset.range M, gcoef m q * z ^ m|
      ≤ C * |z| ^ M * Real.exp (|z| / 2)

/-- `lem:pw:logA`, "Consequently" (supplement line 75) and the paper's "$\zeta(n)$ always arrives
together with $\eta^n$" (main.tex 1194): for `|η| < 1/2`, `s ↦ A(η,s)` is analytic at `0`; for
every `m`, `η ↦ ∂_s^m A(η,0)` is analytic at `0` and its `η^a` Taylor coefficient
`m! Ψ_{a,m}` is a rational zeta polynomial homogeneous of weight `a + m`; and
`∂_s^m A(η,0) = 0` for odd `m` (supplement line 68). -/
def stmt_4 : Prop :=
  (∀ η : ℝ, |η| < 1 / 2 → AnalyticAt ℝ (fun s => X05.A η s) 0) ∧
  ∀ m : ℕ,
    AnalyticAt ℝ (fun η => iteratedDeriv m (fun s => X05.A η s) 0) 0 ∧
    (∀ a : ℕ, IsZetaPoly (a + m)
      (iteratedDeriv a (fun η => iteratedDeriv m (fun s => X05.A η s) 0) 0 /
        (a.factorial : ℝ))) ∧
    (Odd m → ∀ η : ℝ, |η| < 1 / 2 → iteratedDeriv m (fun s => X05.A η s) 0 = 0)

/-- `lem:pw:asymp` (supplement lines 99–102): for every rhombus and every `M ≥ 1`, for
`|η| ≤ 1/4` (`p = 1/(1+η)`),
`P_{A,B}(p)/(4r) = Σ_{m<M, m even} g_m(q) η^m ∂_s^m A(η,0) + O(|η|^M)`. -/
def stmt_5 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B → ∀ M : ℕ, 1 ≤ M → ∃ C : ℝ, ∀ η : ℝ, |η| ≤ 1 / 4 →
    |X01.I A B η - ∑ m ∈ (Finset.range M).filter Even,
        gcoef m (X01.q A B) * η ^ m * iteratedDeriv m (fun s => X05.A η s) 0|
      ≤ C * |η| ^ M

/-- `prop:pw:rhombus`, first expansion (supplement lines 114–117): there are coefficients
`ĉ_N ∈ ℚ[q²][ζ(2),…,ζ(N)]`, homogeneous of weight `N` (`ĉ_0 = 1`, `ĉ_1 = 0`), the same for every
rhombus and every order, with `P_{A,B}(p)/(4r) = Σ_{N≤K} ĉ_N(q) τ^N + O(τ^{K+1})` as `τ → 0⁺`
(the sign `(-1)^N` of the supplement is absorbed into `ĉ_N`). -/
def stmt_6 : Prop :=
  ∃ c : ℕ → ℝ → ℝ, (∀ x, c 0 x = 1) ∧ (∀ x, c 1 x = 0) ∧ (∀ N, IsZetaPolyQ2 N (c N)) ∧
    ∀ A B : ℝ, 0 < A → 0 < B → ∀ K : ℕ, ∃ C τ₀ : ℝ, 0 < τ₀ ∧
      ∀ p : ℝ, 1 < p → tau p ≤ τ₀ →
        |X01.P A B p / (4 * X01.r A B) -
            ∑ N ∈ Finset.range (K + 1), c N (X01.q A B) * tau p ^ N|
          ≤ C * tau p ^ (K + 1)

/-- `thm:pw`, "The same holds for every rhombus (with polynomials in $q^2$ as coefficients)"
(supplement `prop:pw:rhombus`, line 118): there are `w_k ∈ ℚ[q²][ζ(2),…,ζ(k)]` homogeneous of
weight `k`, one sequence for all rhombi and all orders, with
`log(P_{A,B}(p)/(4r)) = Σ_{k=2}^K w_k(q) τ^k + O(τ^{K+1})` as `τ → 0⁺`; and
`w_2 = (1-q²)ζ(2)/4`, `w_3 = (1-q²)ζ(3)/4` (supplement lines 155, 249). -/
def stmt_7 : Prop :=
  ∃ w : ℕ → ℝ → ℝ, (∀ k, IsZetaPolyQ2 k (w k)) ∧
    (∀ x, w 2 x = (1 - x ^ 2) * X05.z 2 / 4 ∧ w 3 x = (1 - x ^ 2) * X05.z 3 / 4) ∧
    ∀ A B : ℝ, 0 < A → 0 < B → ∀ K : ℕ, 2 ≤ K → ∃ C τ₀ : ℝ, 0 < τ₀ ∧
      ∀ p : ℝ, 1 < p → tau p ≤ τ₀ →
        |Real.log (X01.P A B p / (4 * X01.r A B)) -
            ∑ k ∈ Finset.Icc 2 K, w k (X01.q A B) * tau p ^ k|
          ≤ C * tau p ^ (K + 1)

/-- `thm:pw`, "Hence the coefficient of $\tau^n$ in $P$ has pure weight $n$", to every order:
there is one sequence `c_n`, `c_0 = 4`, each `c_n` a rational zeta polynomial homogeneous of
weight `n`, with `P(p)/√2 = Σ_{n≤K} c_n τ^n + O(τ^{K+1})` as `τ → 0⁺`, for every `K`. -/
def stmt_8 : Prop :=
  ∃ c : ℕ → ℝ, c 0 = 4 ∧ (∀ n, IsZetaPoly n (c n)) ∧
    ∀ K : ℕ, ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → tau p ≤ τ₀ →
      |X01.P 1 1 p / Real.sqrt 2 - ∑ n ∈ Finset.range (K + 1), c n * tau p ^ n|
        ≤ C * tau p ^ (K + 1)

/-- Supplement line 153, the logarithmic coefficients of the squircle through `τ⁵`:
`w_2 = ζ(2)/4`, `w_3 = ζ(3)/4`, `w_4 = -ζ(2)²/4 + ζ(4)/32`, `w_5 = -ζ(2)ζ(3)/2 - 3ζ(5)/16`,
i.e. `|log(P/(4√2)) - (w_2τ² + w_3τ³ + w_4τ⁴ + w_5τ⁵)| ≤ C τ⁶` for `p > 1`, `τ ≤ τ₀`. -/
def stmt_9 : Prop :=
  ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → tau p ≤ τ₀ →
    |Real.log (X01.P 1 1 p / (4 * Real.sqrt 2)) -
        (X05.z 2 / 4 * tau p ^ 2 + X05.z 3 / 4 * tau p ^ 3
          + (-(X05.z 2) ^ 2 / 4 + X05.z 4 / 32) * tau p ^ 4
          + (-(X05.z 2 * X05.z 3) / 2 - 3 * X05.z 5 / 16) * tau p ^ 5)|
      ≤ C * tau p ^ 6

/-- `thm:pw`, first sentence (main.tex 1184–1185), to every finite order, for the squircle:
there is one sequence `w_k`, each a rational polynomial in `ζ(2),…,ζ(k)` homogeneous of weight
`k`, such that for every `K ≥ 2` there are `C`, `τ₀ > 0` with
`|log(P(p)/(4√2)) - Σ_{k=2}^K w_k τ^k| ≤ C τ^{K+1}` for all `p > 1` with `τ = 1 - 1/p ≤ τ₀`. -/
def stmt_10 : Prop :=
  ∃ w : ℕ → ℝ, (∀ k, IsZetaPoly k (w k)) ∧
    ∀ K : ℕ, 2 ≤ K → ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → tau p ≤ τ₀ →
      |Real.log (X01.P 1 1 p / (4 * Real.sqrt 2)) -
          ∑ k ∈ Finset.Icc 2 K, w k * tau p ^ k|
        ≤ C * tau p ^ (K + 1)

end Q01
-- END STATEMENTS

namespace Q01

lemma zp_zero (k : ℕ) : IsZetaPoly k 0 := ⟨[], by simp, by simp⟩

lemma zp_const (c : ℚ) : IsZetaPoly 0 (c : ℝ) := ⟨[(c, [])], by simp, by simp⟩

lemma zp_add {k : ℕ} {x y : ℝ} (hx : IsZetaPoly k x) (hy : IsZetaPoly k y) :
    IsZetaPoly k (x + y) := by
  obtain ⟨L1, h1, rfl⟩ := hx
  obtain ⟨L2, h2, rfl⟩ := hy
  refine ⟨L1 ++ L2, fun e he => ?_, by simp⟩
  rcases List.mem_append.1 he with h | h
  · exact h1 e h
  · exact h2 e h

lemma zp_smul {k : ℕ} {x : ℝ} (c : ℚ) (hx : IsZetaPoly k x) : IsZetaPoly k ((c : ℝ) * x) := by
  obtain ⟨L, h1, rfl⟩ := hx
  refine ⟨L.map fun e => (c * e.1, e.2), fun e he => ?_, ?_⟩
  · obtain ⟨a, ha, rfl⟩ := List.mem_map.1 he
    exact h1 a ha
  · rw [← List.sum_map_mul_left]
    simp only [List.map_map]
    congr 1
    apply List.map_congr_left
    intro e _
    simp only [Function.comp, Rat.cast_mul]; ring

lemma zp_mono_mul {b : ℕ} {y : ℝ} (c : ℚ) (ns : List ℕ) (hns : ∀ n ∈ ns, 2 ≤ n)
    (hy : IsZetaPoly b y) :
    IsZetaPoly (ns.sum + b) ((c : ℝ) * (ns.map X05.z).prod * y) := by
  obtain ⟨L, h1, rfl⟩ := hy
  refine ⟨L.map fun e => (c * e.1, ns ++ e.2), fun e he => ?_, ?_⟩
  · obtain ⟨a, ha, rfl⟩ := List.mem_map.1 he
    refine ⟨fun n hn => ?_, ?_⟩
    · rcases List.mem_append.1 hn with h | h
      · exact hns n h
      · exact (h1 a ha).1 n h
    · simp [(h1 a ha).2]
  · rw [← List.sum_map_mul_left]
    simp only [List.map_map]
    congr 1
    apply List.map_congr_left
    intro e _
    simp only [Function.comp, Rat.cast_mul, List.map_append, List.prod_append]; ring

lemma zp_mul {a b : ℕ} {x y : ℝ} (hx : IsZetaPoly a x) (hy : IsZetaPoly b y) :
    IsZetaPoly (a + b) (x * y) := by
  obtain ⟨L, h1, rfl⟩ := hx
  induction L with
  | nil => simpa using zp_zero (a + b)
  | cons e L ih =>
    simp only [List.map_cons, List.sum_cons, add_mul]
    refine zp_add ?_ (ih fun e' he' => h1 e' (List.mem_cons_of_mem _ he'))
    have := zp_mono_mul e.1 e.2 (h1 e (List.mem_cons_self)).1 hy
    rwa [(h1 e (List.mem_cons_self)).2] at this

lemma zp_sum {k : ℕ} {ι : Type*} (s : Finset ι) (f : ι → ℝ) (hf : ∀ i ∈ s, IsZetaPoly k (f i)) :
    IsZetaPoly k (∑ i ∈ s, f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using zp_zero k
  | insert i s hi ih =>
    rw [Finset.sum_insert hi]
    exact zp_add (hf i (Finset.mem_insert_self _ _))
      (ih fun j hj => hf j (Finset.mem_insert_of_mem hj))

lemma zp_cast {k k' : ℕ} {x : ℝ} (h : k = k') (hx : IsZetaPoly k x) : IsZetaPoly k' x := h ▸ hx

/-- the term `k w_k c_{n-k}` is of weight `n` -/
lemma zp_term {n k : ℕ} (hk : k ≤ n) {w c : ℝ} (hw : IsZetaPoly k w) (hc : IsZetaPoly (n - k) c) :
    IsZetaPoly n ((k : ℝ) * w * c) := by
  have := zp_mul (zp_smul (k : ℚ) hw) hc
  rw [Nat.add_sub_cancel' hk] at this
  simpa [mul_assoc] using this

end Q01

theorem Q01.proof_1 : Q01.stmt_1 := by
  refine ⟨fun k x y => Q01.zp_add, fun k x c => Q01.zp_smul c,
    fun a b x y => Q01.zp_mul, fun c w hc0 hw0 hrec => ⟨fun hc => ?_, fun hw => ?_⟩⟩
  · intro n
    induction n using Nat.strong_induction_on with
    | _ n ih =>
      rcases Nat.eq_zero_or_pos n with rfl | hn
      · rw [hw0]; exact Q01.zp_zero 0
      · have h := hrec n hn
        have hsplit : ∑ k ∈ Finset.Icc 1 n, (k : ℝ) * w k * c (n - k) =
            ∑ k ∈ Finset.Ico 1 n, (k : ℝ) * w k * c (n - k) + (n : ℝ) * w n := by
          rw [← Finset.Ico_add_one_right_eq_Icc, Finset.sum_Ico_succ_top (by omega),
            Nat.sub_self, hc0, mul_one]
        have hn' : (n : ℝ) ≠ 0 := by positivity
        have hwn : w n = ((1 / n : ℚ) : ℝ) * ((n : ℚ) * c n
            + ((-1 : ℚ) : ℝ) * ∑ k ∈ Finset.Ico 1 n, (k : ℝ) * w k * c (n - k)) := by
          push_cast
          field_simp
          linarith
        rw [hwn]
        refine Q01.zp_smul _ (Q01.zp_add (Q01.zp_smul _ (hc n)) (Q01.zp_smul _ ?_))
        refine Q01.zp_sum _ _ fun k hk => ?_
        obtain ⟨hk1, hk2⟩ := Finset.mem_Ico.1 hk
        exact Q01.zp_term hk2.le (ih k hk2) (hc _)
  · intro n
    induction n using Nat.strong_induction_on with
    | _ n ih =>
      rcases Nat.eq_zero_or_pos n with rfl | hn
      · rw [hc0]; simpa using Q01.zp_const 1
      · have h := hrec n hn
        have hn' : (n : ℝ) ≠ 0 := by positivity
        have hcn : c n = ((1 / n : ℚ) : ℝ) * ∑ k ∈ Finset.Icc 1 n, (k : ℝ) * w k * c (n - k) := by
          push_cast
          field_simp
          linarith
        rw [hcn]
        refine Q01.zp_smul _ (Q01.zp_sum _ _ fun k hk => ?_)
        obtain ⟨hk1, hk2⟩ := Finset.mem_Icc.1 hk
        exact Q01.zp_term hk2 (hw k) (ih _ (by omega))

namespace Q01

noncomputable def gf (q z : ℝ) : ℝ := Real.sqrt (Real.cosh z + q * Real.sinh z)
noncomputable def hf (q z : ℝ) : ℝ := Real.cosh z + q * Real.sinh z
noncomputable def hf' (q z : ℝ) : ℝ := Real.sinh z + q * Real.cosh z
noncomputable def uf (q z : ℝ) : ℝ := hf' q z / hf q z

open Polynomial in
noncomputable def Qp : ℕ → Polynomial ℚ
  | 0 => 1
  | m + 1 => derivative (Qp m) * (1 - X ^ 2) + Qp m * (C (1 / 2) * X)

lemma hf_deriv (q z : ℝ) : HasDerivAt (hf q) (hf' q z) z := by
  unfold hf hf'
  exact (Real.hasDerivAt_cosh z).add ((Real.hasDerivAt_sinh z).const_mul q)

lemma hf'_deriv (q z : ℝ) : HasDerivAt (hf' q) (hf q z) z := by
  unfold hf hf'
  convert (Real.hasDerivAt_sinh z).add ((Real.hasDerivAt_cosh z).const_mul q) using 1

lemma hf_sq (q z : ℝ) : hf q z ^ 2 - hf' q z ^ 2 = 1 - q ^ 2 := by
  unfold hf hf'
  have := Real.cosh_sq_sub_sinh_sq z
  linear_combination (1 - q ^ 2) * this

lemma uf_deriv {q z : ℝ} (hz : 0 < hf q z) : HasDerivAt (uf q) (1 - uf q z ^ 2) z := by
  have h := (hf'_deriv q z).div (hf_deriv q z) hz.ne'
  unfold uf
  convert h using 1
  field_simp

lemma gf_eq (q : ℝ) : gf q = fun z => Real.sqrt (hf q z) := rfl

lemma gf_deriv {q z : ℝ} (hz : 0 < hf q z) : HasDerivAt (gf q) (uf q z / 2 * gf q z) z := by
  have h := (hf_deriv q z).sqrt hz.ne'
  rw [gf_eq]
  convert h using 1
  unfold uf
  have hs : Real.sqrt (hf q z) ^ 2 = hf q z := Real.sq_sqrt hz.le
  have hs0 : Real.sqrt (hf q z) ≠ 0 := (Real.sqrt_pos.2 hz).ne'
  field_simp
  rw [hs]

lemma hf_cont (q : ℝ) : Continuous (hf q) := by unfold hf; fun_prop

lemma gf_struct (q : ℝ) (m : ℕ) : ∀ z, 0 < hf q z →
    iteratedDeriv m (gf q) z = Polynomial.aeval (uf q z) (Qp m) * gf q z := by
  induction m with
  | zero => intro z _; simp [Qp]
  | succ m ih =>
    intro z hz
    have hU : IsOpen {y | 0 < hf q y} := isOpen_lt continuous_const (hf_cont q)
    have hev : iteratedDeriv m (gf q) =ᶠ[nhds z]
        fun y => Polynomial.aeval (uf q y) (Qp m) * gf q y :=
      Filter.eventually_of_mem (hU.mem_nhds hz) fun y hy => ih y hy
    rw [iteratedDeriv_succ, hev.deriv_eq]
    have h1 := (((Qp m).hasDerivAt_aeval (uf q z)).comp z (uf_deriv hz)).mul (gf_deriv hz)
    have h2 : HasDerivAt (fun y => Polynomial.aeval (uf q y) (Qp m) * gf q y) _ z := h1
    rw [h2.deriv]
    simp only [Qp, map_add, map_mul, map_sub, map_one, map_pow, Polynomial.aeval_X,
      Polynomial.aeval_C, Function.comp]
    simp only [eq_ratCast, Rat.cast_div, Rat.cast_one, Rat.cast_ofNat]
    ring

lemma hf_zero (q : ℝ) : hf q 0 = 1 := by simp [hf]
lemma uf_zero (q : ℝ) : uf q 0 = q := by simp [uf, hf, hf']
lemma gf_zero (q : ℝ) : gf q 0 = 1 := by simp [gf]

lemma gcoef_eq (m : ℕ) (q : ℝ) : gcoef m q = Polynomial.aeval q (Qp m) / (m.factorial : ℝ) := by
  unfold gcoef
  have := gf_struct q m 0 (by rw [hf_zero]; norm_num)
  rw [uf_zero, gf_zero, mul_one] at this
  rw [← this]; rfl

end Q01

theorem Q01.proof_2 : Q01.stmt_2 := by
  refine ⟨fun q n => ?_, fun m => ?_, fun m q => ?_, fun q => ?_⟩
  · have h : ContDiffAt ℝ n (fun z : ℝ => Real.cosh z + q * Real.sinh z) 0 :=
      (Real.contDiff_cosh.add (contDiff_const.mul Real.contDiff_sinh)).contDiffAt
    exact h.sqrt (by simp)
  · refine ⟨Q01.Qp m * Polynomial.C (1 / (m.factorial : ℚ)), fun q => ?_⟩
    rw [Q01.gcoef_eq, map_mul, Polynomial.aeval_C]
    simp only [eq_ratCast, Rat.cast_div, Rat.cast_one, Rat.cast_natCast]
    ring
  · have he : Q01.gf (-q) = fun z => Q01.gf q (-z) := by
      funext z; simp [Q01.gf, Real.cosh_neg, Real.sinh_neg]
    unfold Q01.gcoef
    change iteratedDeriv m (Q01.gf (-q)) 0 / _ = (-1) ^ m * (iteratedDeriv m (Q01.gf q) 0 / _)
    rw [he, iteratedDeriv_comp_neg, neg_zero, smul_eq_mul]
    ring
  · simp only [Q01.gcoef_eq, Q01.Qp, map_add, map_mul, map_sub, map_one, map_pow,
      Polynomial.aeval_X, Polynomial.aeval_C, Polynomial.derivative_mul, Polynomial.derivative_one, 
      Polynomial.derivative_C, Polynomial.derivative_X, Polynomial.derivative_pow,
      eq_ratCast, Rat.cast_div, Rat.cast_one, Rat.cast_ofNat, map_natCast, map_zero]
    norm_num [Nat.factorial]
    refine ⟨by ring, by ring, by ring⟩

namespace Q01

lemma hf_exp (q z : ℝ) : hf q z = ((1 + q) * Real.exp z + (1 - q) * Real.exp (-z)) / 2 := by
  unfold hf; rw [Real.cosh_eq, Real.sinh_eq]; ring

lemma hf'_exp (q z : ℝ) : hf' q z = ((1 + q) * Real.exp z - (1 - q) * Real.exp (-z)) / 2 := by
  unfold hf'; rw [Real.cosh_eq, Real.sinh_eq]; ring

lemma hf_pos {q : ℝ} (hq : |q| ≤ 1) (z : ℝ) : 0 < hf q z := by
  rw [hf_exp]
  obtain ⟨h1, h2⟩ := abs_le.1 hq
  have e1 := Real.exp_pos z
  have e2 := Real.exp_pos (-z)
  rcases le_total q 0 with h | h
  · nlinarith [mul_nonneg (by linarith : (0:ℝ) ≤ 1 + q) e1.le]
  · nlinarith [mul_nonneg (by linarith : (0:ℝ) ≤ 1 - q) e2.le]

lemma hf_le {q : ℝ} (hq : |q| ≤ 1) (z : ℝ) : hf q z ≤ Real.exp |z| := by
  have h1 : q * Real.sinh z ≤ |Real.sinh z| := by
    calc q * Real.sinh z ≤ |q * Real.sinh z| := le_abs_self _
      _ = |q| * |Real.sinh z| := abs_mul _ _
      _ ≤ 1 * |Real.sinh z| := by gcongr
      _ = _ := one_mul _
  have h2 : Real.cosh z + |Real.sinh z| = Real.exp |z| := by
    rw [← Real.cosh_abs, Real.abs_sinh, Real.cosh_add_sinh]
  unfold hf; linarith

lemma uf_bd {q : ℝ} (hq : |q| ≤ 1) (z : ℝ) : |uf q z| ≤ 1 := by
  have hp := hf_pos hq z
  unfold uf
  rw [abs_div, abs_of_pos hp, div_le_one hp, abs_le]
  obtain ⟨h1, h2⟩ := abs_le.1 hq
  have e1 := Real.exp_pos z
  have e2 := Real.exp_pos (-z)
  rw [hf_exp, hf'_exp]
  constructor
  · nlinarith [mul_nonneg (by linarith : (0:ℝ) ≤ 1 + q) e1.le]
  · nlinarith [mul_nonneg (by linarith : (0:ℝ) ≤ 1 - q) e2.le]

lemma gf_le {q : ℝ} (hq : |q| ≤ 1) (z : ℝ) : gf q z ≤ Real.exp (|z| / 2) := by
  change Real.sqrt (hf q z) ≤ _
  rw [Real.sqrt_le_iff]
  refine ⟨(Real.exp_pos _).le, ?_⟩
  rw [sq, ← Real.exp_add, add_halves]
  exact hf_le hq z

lemma gf_nonneg (q z : ℝ) : 0 ≤ gf q z := Real.sqrt_nonneg _

lemma gf_contDiff {q : ℝ} (hq : |q| ≤ 1) (n : ℕ) : ContDiff ℝ n (gf q) :=
  (Real.contDiff_cosh.add (contDiff_const.mul Real.contDiff_sinh)).sqrt
    fun z => (hf_pos hq z).ne'

lemma abs_of_mem_uIoo {ξ z : ℝ} (h : ξ ∈ Set.uIoo 0 z) : |ξ| ≤ |z| := by
  simp only [Set.uIoo, Set.mem_Ioo] at h
  rcases le_total 0 z with hz | hz
  · rw [min_eq_left hz, max_eq_right hz] at h
    rw [abs_of_pos h.1, abs_of_nonneg hz]; exact h.2.le
  · rw [min_eq_right hz, max_eq_left hz] at h
    rw [abs_of_neg h.2, abs_of_nonpos hz]; linarith [h.1]

end Q01

theorem Q01.proof_3 : Q01.stmt_3 := by
  intro M hM
  obtain ⟨B, hB⟩ := (isCompact_Icc (a := (-1 : ℝ)) (b := 1)).exists_bound_of_continuousOn
    (Q01.Qp M).continuous_aeval.continuousOn
  refine ⟨max B 0 / (M.factorial : ℝ), fun q hq z => ?_⟩
  change |Q01.gf q z - ∑ m ∈ Finset.range M, Q01.gcoef m q * z ^ m| ≤ _
  have hfac : (0 : ℝ) < M.factorial := by positivity
  rcases eq_or_ne z 0 with rfl | hz
  · obtain ⟨M', rfl⟩ : ∃ M', M = M' + 1 := ⟨M - 1, by omega⟩
    rw [Finset.sum_range_succ']
    simp only [ne_eq, Nat.add_one_ne_zero, not_false_eq_true, zero_pow, mul_zero,
      Finset.sum_const_zero, pow_zero, mul_one, zero_add, Q01.gf_zero]
    have : Q01.gcoef 0 q = 1 := by simp [Q01.gcoef_eq, Q01.Qp]
    rw [this]; simp
  · obtain ⟨M', rfl⟩ : ∃ M', M = M' + 1 := ⟨M - 1, by omega⟩
    have hcd := Q01.gf_contDiff hq (M' + 1)
    obtain ⟨ξ, hξ, heq⟩ := taylor_mean_remainder_lagrange_iteratedDeriv (x₀ := 0) (x := z)
      (n := M') hz.symm (by exact_mod_cast hcd.contDiffOn)
    have htay : taylorWithinEval (Q01.gf q) M' (Set.uIcc 0 z) 0 z =
        ∑ m ∈ Finset.range (M' + 1), Q01.gcoef m q * z ^ m := by
      rw [taylor_within_apply]
      refine Finset.sum_congr rfl fun k _ => ?_
      rw [iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_uIcc hz.symm)
        ((Q01.gf_contDiff hq k).contDiffAt) (Set.left_mem_uIcc (a := (0:ℝ)) (b := z))]
      unfold Q01.gcoef
      simp only [smul_eq_mul, sub_zero]
      change _ = iteratedDeriv k (Q01.gf q) 0 / _ * _
      field_simp
    rw [← htay, heq, Q01.gf_struct q (M' + 1) ξ (Q01.hf_pos hq ξ), sub_zero]
    have h1 : |Polynomial.aeval (Q01.uf q ξ) (Q01.Qp (M' + 1))| ≤ max B 0 := by
      have := hB (Q01.uf q ξ) (abs_le.1 (Q01.uf_bd hq ξ))
      rw [Real.norm_eq_abs] at this
      exact this.trans (le_max_left _ _)
    have h2 : Q01.gf q ξ ≤ Real.exp (|z| / 2) :=
      (Q01.gf_le hq ξ).trans (Real.exp_le_exp.2 (by
        have := Q01.abs_of_mem_uIoo hξ; linarith))
    rw [abs_div, abs_mul, abs_mul, abs_pow, abs_of_nonneg (Q01.gf_nonneg q ξ),
      abs_of_pos hfac]
    have h3 : 0 ≤ |z| ^ (M' + 1) := by positivity
    calc |Polynomial.aeval (Q01.uf q ξ) (Q01.Qp (M' + 1))| * Q01.gf q ξ * |z| ^ (M' + 1)
          / ((M' + 1).factorial : ℝ)
        ≤ max B 0 * Real.exp (|z| / 2) * |z| ^ (M' + 1) / ((M' + 1).factorial : ℝ) := by
          have : (0:ℝ) ≤ max B 0 := le_max_right B 0
          have := Q01.gf_nonneg q ξ
          gcongr
      _ = _ := by ring

/-- helper for `stmt_4` (odd part): `logit (1 - t) = - logit t`. -/
lemma Q01.logit_symm (t : ℝ) : X05.logit (1 - t) = -X05.logit t := by
  unfold X05.logit; rw [sub_sub_cancel, ← Real.log_inv, inv_div]

/-- helper for `stmt_4`: `A(η,-s) = A(η,s)` for all `η, s` (substitution `t ↦ 1-t`). -/
lemma Q01.A_neg_all (η s : ℝ) : X05.A η (-s) = X05.A η s := by
  unfold X05.A
  refine congrArg (fun I => (1 + η) * I) ?_
  have h := intervalIntegral.integral_comp_sub_left
    (fun t : ℝ => (t * (1 - t)) ^ (η / 2) * Real.exp (s * X05.logit t)) (a := 0) (b := 1) 1
  simp only [sub_zero, sub_self] at h
  rw [← h]
  refine intervalIntegral.integral_congr fun t _ => ?_
  simp only [Q01.logit_symm, sub_sub_cancel]
  rw [mul_comm (1 - t) t]; ring_nf

/-- helper for `stmt_4`: odd `s`-derivatives of `A` at `s = 0` vanish. -/
lemma Q01.iter_A_odd (η : ℝ) {m : ℕ} (hm : Odd m) :
    iteratedDeriv m (fun s => X05.A η s) 0 = 0 := by
  have h := iteratedDeriv_comp_neg m (fun s => X05.A η s) 0
  simp only [Q01.A_neg_all, neg_zero, smul_eq_mul, hm.neg_one_pow] at h
  linarith

namespace Q01

open Filter Topology

/-- coefficients `c_n = (-1)^n ζ(n)/n` (`n ≥ 2`) of `ℓ(x) = log Γ(1+x) + γ x`. -/
noncomputable def cz (n : ℕ) : ℝ := if 2 ≤ n then (-1 : ℝ) ^ n * X05.z n / n else 0

lemma cz_bd (n : ℕ) : |cz n| ≤ 1 := by
  unfold cz
  split_ifs with h
  · have hz := X05.z_bd h
    have hn : (2 : ℝ) ≤ n := by exact_mod_cast h
    rw [abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
      abs_of_pos (by linarith : (0:ℝ) < n), div_le_one (by linarith)]
    linarith
  · simp

lemma cz_pure (n : ℕ) : IsZetaPoly n (cz n) := by
  unfold cz
  split_ifs with h
  · refine ⟨[((-1 : ℚ) ^ n / n, [n])], ?_, ?_⟩
    · simp [h]
    · simp; ring
  · exact zp_zero n

noncomputable def lser : FormalMultilinearSeries ℝ ℝ ℝ := FormalMultilinearSeries.ofScalars ℝ cz

lemma lser_radius : (1 : ENNReal) ≤ lser.radius := by
  refine ENNReal.le_of_forall_nnreal_lt fun r hr => ?_
  refine lser.le_radius_of_bound 1 fun n => ?_
  have hr1 : (r : ℝ) < 1 := by exact_mod_cast hr
  unfold lser
  rw [FormalMultilinearSeries.ofScalars_norm, Real.norm_eq_abs]
  calc |cz n| * (r : ℝ) ^ n ≤ 1 * 1 := by
        gcongr
        · exact cz_bd n
        · exact pow_le_one₀ r.2 hr1.le
    _ = 1 := one_mul 1

/-- `ℓ(x) = Σ_{n≥2} (-1)^n ζ(n)/n x^n`. -/
noncomputable def ell : ℝ → ℝ := lser.sum

lemma ell_hasFPS : HasFPowerSeriesOnBall ell lser 0 lser.radius :=
  lser.hasFPowerSeriesOnBall (lt_of_lt_of_le one_pos lser_radius)

lemma mem_ball1 {x : ℝ} (hx : |x| < 1) : x ∈ Metric.eball (0 : ℝ) lser.radius := by
  refine lt_of_lt_of_le ?_ lser_radius
  rw [edist_zero_right, Real.eNorm_eq_ofReal_abs]
  exact ENNReal.ofReal_lt_one.2 hx

lemma ell_hasSum {x : ℝ} (hx : |x| < 1) : HasSum (fun n => cz n * x ^ n) (ell x) := by
  have h := ell_hasFPS.hasSum (mem_ball1 hx)
  rw [zero_add] at h
  unfold lser at h
  simpa [FormalMultilinearSeries.ofScalars_apply_eq, mul_comm] using h

lemma ell_analytic {x : ℝ} (hx : |x| < 1) : AnalyticAt ℝ ell x :=
  ell_hasFPS.analyticAt_of_mem (by simpa using mem_ball1 hx)

lemma ell_coef (n : ℕ) : iteratedDeriv n ell 0 / (n.factorial : ℝ) = cz n := by
  have h1 := (ell_analytic (x := 0) (by simp)).hasFPowerSeriesAt
  have h2 := h1.eq_formalMultilinearSeries ell_hasFPS.hasFPowerSeriesAt
  unfold lser at h2
  have := FormalMultilinearSeries.ofScalars_series_injective ℝ ℝ h2
  exact congrFun this n

/-- pure weight-`d` germs: analytic at `0`, `a`-th Taylor coefficient of weight `a + d`. -/
def PF (d : ℕ) (f : ℝ → ℝ) : Prop :=
  AnalyticAt ℝ f 0 ∧ ∀ a : ℕ, IsZetaPoly (a + d) (iteratedDeriv a f 0 / (a.factorial : ℝ))

lemma pf_ell : PF 0 ell :=
  ⟨ell_analytic (by simp), fun a => by rw [ell_coef]; exact cz_pure a⟩

lemma pf_congr {d : ℕ} {f g : ℝ → ℝ} (h : f =ᶠ[𝓝 0] g) (hf : PF d f) : PF d g :=
  ⟨hf.1.congr h, fun a => by rw [← h.iteratedDeriv_eq a]; exact hf.2 a⟩

lemma pf_smul {d : ℕ} {f : ℝ → ℝ} (c : ℚ) (hf : PF d f) : PF d (fun x => (c : ℝ) * f x) := by
  refine ⟨analyticAt_const.mul hf.1, fun a => ?_⟩
  rw [iteratedDeriv_const_mul_field, mul_div_assoc]
  exact zp_smul c (hf.2 a)

lemma pf_add {d : ℕ} {f g : ℝ → ℝ} (hf : PF d f) (hg : PF d g) : PF d (fun x => f x + g x) := by
  refine ⟨hf.1.add hg.1, fun a => ?_⟩
  have := iteratedDeriv_add (n := a) hf.1.contDiffAt hg.1.contDiffAt
  rw [show (fun x => f x + g x) = f + g from rfl, this, add_div]
  exact zp_add (hf.2 a) (hg.2 a)

lemma pf_zero (d : ℕ) : PF d (fun _ => 0) :=
  ⟨analyticAt_const, fun a => by simp; exact zp_zero _⟩

lemma pf_sum {d : ℕ} {ι : Type*} (s : Finset ι) (f : ι → ℝ → ℝ) (hf : ∀ i ∈ s, PF d (f i)) :
    PF d (fun x => ∑ i ∈ s, f i x) := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using pf_zero d
  | insert i s hi ih =>
    simp only [Finset.sum_insert hi]
    exact pf_add (hf i (Finset.mem_insert_self _ _))
      (ih fun j hj => hf j (Finset.mem_insert_of_mem hj))

lemma pf_mul {d e : ℕ} {f g : ℝ → ℝ} (hf : PF d f) (hg : PF e g) :
    PF (d + e) (fun x => f x * g x) := by
  refine ⟨hf.1.mul hg.1, fun a => ?_⟩
  have := iteratedDeriv_mul (n := a) hf.1.contDiffAt hg.1.contDiffAt
  rw [show (fun x => f x * g x) = f * g from rfl, this, Finset.sum_div]
  refine zp_sum _ _ fun i hi => ?_
  have hi' : i ≤ a := Nat.lt_succ_iff.1 (Finset.mem_range.1 hi)
  have h := zp_mul (hf.2 i) (hg.2 (a - i))
  have hw : i + d + (a - i + e) = a + (d + e) := by omega
  rw [hw] at h
  convert h using 1
  rw [Nat.cast_choose ℝ hi']
  have h1 : (i.factorial : ℝ) ≠ 0 := by positivity
  have h2 : ((a - i).factorial : ℝ) ≠ 0 := by positivity
  have h3 : (a.factorial : ℝ) ≠ 0 := by positivity
  field_simp

lemma pf_deriv {d : ℕ} {f : ℝ → ℝ} (hf : PF d f) : PF (d + 1) (deriv f) := by
  refine ⟨hf.1.deriv, fun a => ?_⟩
  rw [← iteratedDeriv_succ']
  have h := zp_smul ((a + 1 : ℕ) : ℚ) (hf.2 (a + 1))
  rw [show a + 1 + d = a + (d + 1) by omega] at h
  convert h using 1
  rw [Nat.factorial_succ]; push_cast
  have : (a.factorial : ℝ) ≠ 0 := by positivity
  field_simp

lemma pf_iter {d : ℕ} {f : ℝ → ℝ} (hf : PF d f) (j : ℕ) : PF (d + j) (iteratedDeriv j f) := by
  induction j with
  | zero => simpa using hf
  | succ j ih => rw [iteratedDeriv_succ]; exact pf_deriv ih

lemma iter_comp_half {f : ℝ → ℝ} (hf : AnalyticAt ℝ f 0) (n : ℕ) :
    iteratedDeriv n (fun x => f ((1 / 2 : ℝ) * x)) 0 = (1 / 2 : ℝ) ^ n * iteratedDeriv n f 0 := by
  obtain ⟨R, hR, hfR⟩ := hf.exists_ball_analyticOnNhd
  have hs : IsOpen (Metric.ball (0 : ℝ) R) := Metric.isOpen_ball
  have h0 : (0 : ℝ) ∈ Metric.ball (0 : ℝ) R := Metric.mem_ball_self hR
  have hmaps : Set.MapsTo ((1 / 2 : ℝ) * ·) (Metric.ball (0 : ℝ) R) (Metric.ball (0 : ℝ) R) := by
    intro x hx
    simp only [Metric.mem_ball, dist_zero_right, Real.norm_eq_abs] at hx ⊢
    rw [abs_mul]; norm_num; linarith [abs_nonneg x]
  have hcd : ContDiffOn ℝ n f (Metric.ball (0 : ℝ) R) := (hfR.contDiffOn hs.uniqueDiffOn).of_le
    le_top
  have key := iteratedDerivWithin_comp_const_smul h0 hs.uniqueDiffOn hcd (1 / 2 : ℝ) hmaps
  rw [iteratedDerivWithin_of_isOpen hs h0, mul_zero,
    iteratedDerivWithin_of_isOpen hs h0, smul_eq_mul] at key
  exact key

lemma pf_half {d : ℕ} {f : ℝ → ℝ} (hf : PF d f) : PF d (fun x => f ((1 / 2 : ℝ) * x)) := by
  refine ⟨?_, fun a => ?_⟩
  · have : AnalyticAt ℝ (fun x : ℝ => (1 / 2 : ℝ) * x) 0 := analyticAt_const.mul analyticAt_id
    exact (by simpa using hf.1 : AnalyticAt ℝ f ((1 / 2 : ℝ) * 0)).comp this
  · rw [iter_comp_half hf.1, mul_div_assoc]
    have := zp_smul ((1 / 2 : ℚ) ^ a) (hf.2 a)
    simpa using this

/-- `exp` of a pure weight-0 germ vanishing at `0` is pure (by `proof_1`'s recursion). -/
lemma pf_exp {f : ℝ → ℝ} (hf : PF 0 f) (h0 : f 0 = 0) : PF 0 (fun x => Real.exp (f x)) := by
  have hE : AnalyticAt ℝ (fun x => Real.exp (f x)) 0 := hf.1.rexp'
  refine ⟨hE, fun a => ?_⟩
  set c : ℕ → ℝ := fun n => iteratedDeriv n (fun x => Real.exp (f x)) 0 / (n.factorial : ℝ)
  set w : ℕ → ℝ := fun n => iteratedDeriv n f 0 / (n.factorial : ℝ)
  have hc0 : c 0 = 1 := by simp [c, h0]
  have hw0 : w 0 = 0 := by simp [w, h0]
  -- derivative of `exp ∘ f` near 0
  have hder : deriv (fun x => Real.exp (f x)) =ᶠ[𝓝 0] fun x => deriv f x * Real.exp (f x) := by
    filter_upwards [hf.1.eventually_analyticAt] with x hx
    rw [(hx.differentiableAt.hasDerivAt.exp).deriv]; ring
  have hrec : ∀ n : ℕ, 1 ≤ n →
      (n : ℝ) * c n = ∑ k ∈ Finset.Icc 1 n, (k : ℝ) * w k * c (n - k) := by
    intro n hn
    obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
    have h1 : iteratedDeriv (m + 1) (fun x => Real.exp (f x)) 0 =
        ∑ i ∈ Finset.range (m + 1), (m.choose i : ℝ) * iteratedDeriv (i + 1) f 0 *
          iteratedDeriv (m - i) (fun x => Real.exp (f x)) 0 := by
      rw [iteratedDeriv_succ', hder.iteratedDeriv_eq]
      have := iteratedDeriv_mul (n := m) hf.1.deriv.contDiffAt hE.contDiffAt
      rw [show (fun x => deriv f x * Real.exp (f x)) = deriv f * (fun x => Real.exp (f x))
        from rfl, this]
      refine Finset.sum_congr rfl fun i _ => ?_
      rw [← iteratedDeriv_succ']
    have hIcc : Finset.Icc 1 (m + 1) = Finset.Ico 1 (m + 2) := rfl
    rw [hIcc, Finset.sum_Ico_eq_sum_range, show m + 2 - 1 = m + 1 by omega]
    simp only [c, w, h1, Finset.mul_sum, Finset.sum_div]
    refine Finset.sum_congr rfl fun i hi => ?_
    have hi' : i ≤ m := Nat.lt_succ_iff.1 (Finset.mem_range.1 hi)
    rw [show m + 1 - (1 + i) = m - i by omega, add_comm 1 i, Nat.cast_choose ℝ hi']
    rw [Nat.factorial_succ, Nat.factorial_succ]
    push_cast
    have h1 : (i.factorial : ℝ) ≠ 0 := by positivity
    have h2 : ((m - i).factorial : ℝ) ≠ 0 := by positivity
    have h3 : (m.factorial : ℝ) ≠ 0 := by positivity
    field_simp
  have := ((Q01.proof_1).2.2.2 c w hc0 hw0 hrec).2 (fun n => by simpa using hf.2 n) a
  simpa using this

end Q01

namespace Q01

open Filter Topology

lemma pf_cast {d d' : ℕ} {f : ℝ → ℝ} (h : d = d') (hf : PF d f) : PF d' f := h ▸ hf

lemma ell_zero : ell 0 = 0 := by
  have := ell_coef 0
  simpa [cz] using this

/-- `G(x) = exp ℓ(x) = Γ(1+x) e^{γx}`. -/
noncomputable def Gf (x : ℝ) : ℝ := Real.exp (ell x)
/-- `1/G`. -/
noncomputable def Hf (x : ℝ) : ℝ := Real.exp (-ell x)

lemma pf_G : PF 0 Gf := pf_exp pf_ell ell_zero

lemma pf_H : PF 0 Hf := by
  have h1 : PF 0 (fun x => ((-1 : ℚ) : ℝ) * ell x) := pf_smul (-1) pf_ell
  have h2 := pf_exp h1 (by simp [ell_zero])
  refine pf_congr (Eventually.of_forall fun x => ?_) h2
  simp [Hf]

lemma G_analytic {x : ℝ} (hx : |x| < 1) : AnalyticAt ℝ Gf x := (ell_analytic hx).rexp'

/-- Step 3 in product form: `A(η,s) = G(η/2+s) G(η/2-s) / G(η)`. -/
lemma A_eq_G {η s : ℝ} (hη : |η| < 1 / 2) (hs : |s| < 1 / 2) :
    X05.A η s = Gf ((1 / 2) * η + s) * Gf ((1 / 2) * η - s) * Hf η := by
  obtain ⟨hpos, hsum⟩ := X05.proof_3 η s hη hs
  have hu : |(1 / 2) * η + s| < 1 := by
    calc _ ≤ |(1 / 2) * η| + |s| := abs_add_le _ _
      _ < 1 := by rw [abs_mul]; norm_num; linarith
  have hv : |(1 / 2) * η - s| < 1 := by
    calc _ ≤ |(1 / 2) * η| + |s| := abs_sub _ _
      _ < 1 := by rw [abs_mul]; norm_num; linarith
  have he : |η| < 1 := by linarith
  have h3 := ((ell_hasSum hu).add (ell_hasSum hv)).sub (ell_hasSum he)
  rw [← hasSum_nat_add_iff' 2] at h3
  have h4 : HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * X05.z (k + 2) / ((k + 2 : ℕ) : ℝ) *
        ((η / 2 + s) ^ (k + 2) + (η / 2 - s) ^ (k + 2) - η ^ (k + 2)))
      (ell ((1 / 2) * η + s) + ell ((1 / 2) * η - s) - ell η) := by
    convert h3 using 1
    · funext k
      simp only [cz, show 2 ≤ k + 2 by omega, ite_true]
      ring
    · simp [Finset.sum_range_succ, cz]
  have hl := hsum.unique h4
  rw [← Real.exp_log hpos, hl, Gf, Gf, Hf, ← Real.exp_add, ← Real.exp_add]
  ring_nf

/-- the closed form of `∂_s^m A(η,0)`. -/
noncomputable def Bf (m : ℕ) (η : ℝ) : ℝ :=
  Hf η * ∑ j ∈ Finset.range (m + 1), (((m.choose j : ℚ) * (-1) ^ (m - j) : ℚ) : ℝ) *
    (iteratedDeriv j Gf ((1 / 2) * η) * iteratedDeriv (m - j) Gf ((1 / 2) * η))

lemma pf_B (m : ℕ) : PF m (Bf m) := by
  refine pf_cast (zero_add m) (pf_mul pf_H (pf_sum _ _ fun j hj => ?_))
  have hj' : j ≤ m := Nat.lt_succ_iff.1 (Finset.mem_range.1 hj)
  refine pf_smul _ (pf_cast (by omega) (pf_mul (pf_half (pf_iter pf_G j))
    (pf_half (pf_iter pf_G (m - j)))))

lemma iter_A_eq {η : ℝ} (hη : |η| < 1 / 2) (m : ℕ) :
    iteratedDeriv m (fun s => X05.A η s) 0 = Bf m η := by
  have hev : (fun s => X05.A η s) =ᶠ[𝓝 0]
      fun s => (fun s => Gf ((1 / 2) * η + s) * Gf ((1 / 2) * η - s)) s * Hf η := by
    have : {s : ℝ | |s| < 1 / 2} ∈ 𝓝 (0 : ℝ) := by
      have := Metric.ball_mem_nhds (0 : ℝ) (by norm_num : (0:ℝ) < 1 / 2)
      simpa [Metric.ball, Real.dist_eq] using this
    filter_upwards [this] with s hs
    exact A_eq_G hη hs
  rw [hev.iteratedDeriv_eq, iteratedDeriv_mul_const_field]
  have hx : |(1 / 2) * η| < 1 := by rw [abs_mul]; norm_num; linarith
  have c1 : ContDiffAt ℝ m (fun s => Gf ((1 / 2) * η + s)) 0 := by
    have : AnalyticAt ℝ Gf ((1 / 2) * η + 0) := by simpa using G_analytic hx
    exact (this.comp (analyticAt_const.add analyticAt_id)).contDiffAt
  have c2 : ContDiffAt ℝ m (fun s => Gf ((1 / 2) * η - s)) 0 := by
    have : AnalyticAt ℝ Gf ((1 / 2) * η - 0) := by simpa using G_analytic hx
    exact (this.comp (analyticAt_const.sub analyticAt_id)).contDiffAt
  have := iteratedDeriv_mul (n := m) c1 c2
  rw [show (fun s => Gf ((1 / 2) * η + s) * Gf ((1 / 2) * η - s)) =
    (fun s => Gf ((1 / 2) * η + s)) * (fun s => Gf ((1 / 2) * η - s)) from rfl, this]
  unfold Bf
  rw [mul_comm, Finset.mul_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [iteratedDeriv_comp_const_add j Gf, iteratedDeriv_comp_const_sub (m - j) Gf]
  simp only [add_zero, sub_zero, smul_eq_mul]
  push_cast
  ring

end Q01

theorem Q01.proof_4 : Q01.stmt_4 := by
  have hU : {η : ℝ | |η| < 1 / 2} ∈ nhds (0 : ℝ) := by
    have := Metric.ball_mem_nhds (0 : ℝ) (by norm_num : (0:ℝ) < 1 / 2)
    simpa [Metric.ball, Real.dist_eq] using this
  refine ⟨fun η hη => ?_, fun m => ?_⟩
  · have hx : |(1 / 2) * η| < 1 := by rw [abs_mul]; norm_num; linarith
    have hev : (fun s => X05.A η s) =ᶠ[nhds 0]
        fun s => Q01.Gf ((1 / 2) * η + s) * Q01.Gf ((1 / 2) * η - s) * Q01.Hf η := by
      filter_upwards [hU] with s hs
      exact Q01.A_eq_G hη hs
    refine AnalyticAt.congr ?_ hev.symm
    have a1 : AnalyticAt ℝ Q01.Gf ((1 / 2) * η + 0) := by simpa using Q01.G_analytic hx
    have a2 : AnalyticAt ℝ Q01.Gf ((1 / 2) * η - 0) := by simpa using Q01.G_analytic hx
    exact ((a1.comp (analyticAt_const.add analyticAt_id)).mul
      (a2.comp (analyticAt_const.sub analyticAt_id))).mul analyticAt_const
  · have hev : Q01.Bf m =ᶠ[nhds 0] fun η => iteratedDeriv m (fun s => X05.A η s) 0 := by
      filter_upwards [hU] with η hη
      exact (Q01.iter_A_eq hη m).symm
    have h := Q01.pf_congr hev (Q01.pf_B m)
    exact ⟨h.1, h.2, fun hm η _ => Q01.iter_A_odd η hm⟩

namespace Q01

open MeasureTheory Set Filter Topology

lemma log_bd' {x ε : ℝ} (hx0 : 0 < x) (hx1 : x ≤ 1) (hε : 0 < ε) :
    |Real.log x| ≤ x ^ (-ε) / ε := by
  have h := Real.neg_rpow_div_le_log hx0.le hε
  have hl : Real.log x ≤ 0 := Real.log_nonpos hx0.le hx1
  rw [abs_of_nonpos hl]
  have : -x ^ (-ε) / ε = -(x ^ (-ε) / ε) := by ring
  linarith

lemma logit_bd' {t ε : ℝ} (ht0 : 0 < t) (ht1 : t < 1) (hε : 0 < ε) :
    |X05.logit t| ≤ 2 / ε * (t * (1 - t)) ^ (-ε) := by
  have h1 : 0 < 1 - t := by linarith
  have hy : 0 < t * (1 - t) := mul_pos ht0 h1
  rw [X05.logit, Real.log_div ht0.ne' h1.ne']
  have a := log_bd' ht0 ht1.le hε
  have b := log_bd' h1 (by linarith) hε
  have a' : t ^ (-ε) ≤ (t * (1 - t)) ^ (-ε) :=
    Real.rpow_le_rpow_of_nonpos hy (by nlinarith) (by linarith)
  have b' : (1 - t) ^ (-ε) ≤ (t * (1 - t)) ^ (-ε) :=
    Real.rpow_le_rpow_of_nonpos hy (by nlinarith) (by linarith)
  have a'' : t ^ (-ε) / ε ≤ (t * (1 - t)) ^ (-ε) / ε := div_le_div_of_nonneg_right a' hε.le
  have b'' : (1 - t) ^ (-ε) / ε ≤ (t * (1 - t)) ^ (-ε) / ε := div_le_div_of_nonneg_right b' hε.le
  calc |Real.log t - Real.log (1 - t)| ≤ |Real.log t| + |Real.log (1 - t)| := abs_sub _ _
    _ ≤ (t * (1 - t)) ^ (-ε) / ε + (t * (1 - t)) ^ (-ε) / ε := by linarith
    _ = _ := by ring

/-- the constant `(16(k+1))^k` -/
noncomputable def Kc (k : ℕ) : ℝ := (16 * ((k : ℝ) + 1)) ^ k

lemma Kc_pos (k : ℕ) : 0 < Kc k := by unfold Kc; positivity

lemma logit_pow_bd' {t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) (k : ℕ) :
    |X05.logit t| ^ k ≤ Kc k * (t * (1 - t)) ^ (-(1 / 8 : ℝ)) := by
  have h1 : 0 < 1 - t := by linarith
  have hy : 0 < t * (1 - t) := mul_pos ht0 h1
  have hy1 : t * (1 - t) ≤ 1 := by nlinarith
  set ε : ℝ := 1 / (8 * ((k : ℝ) + 1)) with hε
  have hε0 : 0 < ε := by positivity
  have hb := logit_bd' ht0 ht1 hε0
  have h2 : 2 / ε = 16 * ((k : ℝ) + 1) := by rw [hε]; field_simp; ring
  rw [h2] at hb
  calc |X05.logit t| ^ k ≤ (16 * ((k : ℝ) + 1) * (t * (1 - t)) ^ (-ε)) ^ k :=
        pow_le_pow_left₀ (abs_nonneg _) hb k
    _ = Kc k * (t * (1 - t)) ^ (-ε * k) := by
        rw [mul_pow, Real.rpow_mul_natCast hy.le]; rfl
    _ ≤ Kc k * (t * (1 - t)) ^ (-(1 / 8 : ℝ)) := by
        gcongr ?_ * ?_
        · exact (Kc_pos k).le
        · refine Real.rpow_le_rpow_of_exponent_ge hy hy1 ?_
          rw [hε]
          have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
          rw [neg_mul, neg_le_neg_iff, div_mul_eq_mul_div, one_mul,
            div_le_div_iff₀ (by positivity) (by norm_num)]
          nlinarith

lemma main_bd' {t a b : ℝ} (ht0 : 0 < t) (ht1 : t < 1) (ha : |a| ≤ 1 / 4) (hb : |b| ≤ 1 / 4)
    (m : ℕ) : |X05.logit t| ^ m * (t ^ a * (1 - t) ^ b) ≤ Kc m * X05.q t := by
  have h1 : 0 < 1 - t := by linarith
  have hy : 0 < t * (1 - t) := mul_pos ht0 h1
  have e1 := logit_pow_bd' ht0 ht1 m
  have e2 := X05.pow_pair_bd ht0 ht1 ha hb
  calc |X05.logit t| ^ m * (t ^ a * (1 - t) ^ b)
      ≤ (Kc m * (t * (1 - t)) ^ (-(1 / 8 : ℝ))) * (t * (1 - t)) ^ (-(1 / 2 : ℝ)) :=
        mul_le_mul e1 e2 (by positivity) (mul_nonneg (Kc_pos m).le (by positivity))
    _ = Kc m * X05.q t := by
      rw [X05.q, mul_assoc, ← Real.rpow_add hy]; norm_num

lemma G_integrand_bd' {η s t : ℝ} (hη : |η| ≤ 1 / 4) (hs : |s| ≤ 1 / 8) (ht0 : 0 < t)
    (ht1 : t < 1) (k : ℕ) :
    ‖(t * (1 - t)) ^ (η / 2) * X05.logit t ^ k * Real.exp (s * X05.logit t)‖ ≤
      Kc k * X05.q t := by
  have e : (t * (1 - t)) ^ (η / 2) * X05.logit t ^ k * Real.exp (s * X05.logit t) =
      X05.logit t ^ k * ((t * (1 - t)) ^ (η / 2) * Real.exp (s * X05.logit t)) := by ring
  rw [e, X05.pointwise2 ht0 ht1, Real.norm_eq_abs, abs_mul, abs_pow,
    abs_of_pos (by have : 0 < 1 - t := by linarith
                   positivity : (0:ℝ) < t ^ (η / 2 + s) * (1 - t) ^ (η / 2 - s))]
  have h1 : |η / 2 + s| ≤ 1 / 4 := by
    rw [abs_le] at hη hs ⊢; constructor <;> linarith
  have h2 : |η / 2 - s| ≤ 1 / 4 := by
    rw [abs_le] at hη hs ⊢; constructor <;> linarith
  exact main_bd' ht0 ht1 h1 h2 k

lemma G_int' {η s : ℝ} (hη : |η| ≤ 1 / 4) (hs : |s| ≤ 1 / 8) (k : ℕ) :
    IntervalIntegrable
      (fun t => (t * (1 - t)) ^ (η / 2) * X05.logit t ^ k * Real.exp (s * X05.logit t))
      volume 0 1 := by
  refine IntervalIntegrable.mono_fun' (X05.q_int.const_mul (Kc k)) (by unfold X05.logit; fun_prop) ?_
  rw [uIoc_of_le zero_le_one, ← restrict_Ioo_eq_restrict_Ioc]
  exact ae_restrict_of_forall_mem measurableSet_Ioo fun t ht =>
    G_integrand_bd' hη hs ht.1 ht.2 k

lemma G_deriv' {η s : ℝ} (hη : |η| ≤ 1 / 4) (hs : |s| < 1 / 8) (k : ℕ) :
    HasDerivAt (X05.G η k) (X05.G η (k + 1) s) s := by
  have hS : Set.Ioo (-(1 / 8) : ℝ) (1 / 8) ∈ 𝓝 s :=
    Ioo_mem_nhds (by linarith [neg_abs_le s]) (by linarith [le_abs_self s])
  have := intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (μ := volume) (a := 0) (b := 1)
    (F := fun x t => (t * (1 - t)) ^ (η / 2) * X05.logit t ^ k * Real.exp (x * X05.logit t))
    (F' := fun x t => (t * (1 - t)) ^ (η / 2) * X05.logit t ^ (k + 1) *
      Real.exp (x * X05.logit t))
    (bound := fun t => Kc (k + 1) * X05.q t) hS
    (Filter.Eventually.of_forall fun x => by unfold X05.logit; fun_prop)
    (G_int' hη hs.le k) (by unfold X05.logit; fun_prop) ?_ (X05.q_int.const_mul _) ?_
  · exact this.2
  · filter_upwards [X05.ae_ne_one] with t ht1 htm x hx
    rw [uIoc_of_le zero_le_one] at htm
    have hx' : |x| ≤ 1 / 8 := by rw [abs_le]; exact ⟨hx.1.le, hx.2.le⟩
    exact G_integrand_bd' hη hx' htm.1 (lt_of_le_of_ne htm.2 ht1) (k + 1)
  · refine Filter.Eventually.of_forall fun t _ x _ => ?_
    have h := ((hasDerivAt_id' x).mul_const (X05.logit t)).exp.const_mul
      ((t * (1 - t)) ^ (η / 2) * X05.logit t ^ k)
    convert h using 1
    ring

lemma iter_A' {η : ℝ} (hη : |η| ≤ 1 / 4) :
    ∀ k, ∀ s, |s| < 1 / 8 → iteratedDeriv k (fun s => X05.A η s) s = (1 + η) * X05.G η k s := by
  intro k
  induction k with
  | zero => intro s _; simp [X05.A_eq_G]
  | succ k ih =>
    intro s hs
    rw [iteratedDeriv_succ]
    have hev : iteratedDeriv k (fun s => X05.A η s) =ᶠ[𝓝 s] fun s => (1 + η) * X05.G η k s := by
      have hS : Set.Ioo (-(1 / 8) : ℝ) (1 / 8) ∈ 𝓝 s :=
        Ioo_mem_nhds (by linarith [neg_abs_le s]) (by linarith [le_abs_self s])
      filter_upwards [hS] with x hx
      exact ih x (by rw [abs_lt]; exact hx)
    rw [hev.deriv_eq]
    exact ((G_deriv' hη hs k).const_mul (1 + η)).deriv

lemma iter_A0 {η : ℝ} (hη : |η| ≤ 1 / 4) (k : ℕ) :
    iteratedDeriv k (fun s => X05.A η s) 0 =
      (1 + η) * ∫ t in (0 : ℝ)..1, (t * (1 - t)) ^ (η / 2) * X05.logit t ^ k := by
  rw [iter_A' hη k 0 (by norm_num), X05.G]
  simp

lemma sqrtcosh_bd' {η t : ℝ} (hη : |η| ≤ 1 / 4) (ht0 : 0 < t) (ht1 : t < 1) (k : ℕ) :
    |X05.logit t| ^ k * ((t * (1 - t)) ^ (η / 2) * Real.exp (|η * X05.logit t| / 2)) ≤
      2 * (Kc k * X05.q t) := by
  have hη2 : |η / 2| ≤ 1 / 8 := by rw [abs_div]; norm_num; linarith
  have hη2' : |-(η / 2)| ≤ 1 / 8 := by rwa [abs_neg]
  have a := G_integrand_bd' hη hη2 ht0 ht1 k
  have b := G_integrand_bd' hη hη2' ht0 ht1 k
  have hw : 0 ≤ (t * (1 - t)) ^ (η / 2) := by
    have : 0 < 1 - t := by linarith
    positivity
  rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_pow, abs_of_nonneg hw,
    abs_of_pos (Real.exp_pos _)] at a b
  have he : Real.exp (|η * X05.logit t| / 2) ≤
      Real.exp (η / 2 * X05.logit t) + Real.exp (-(η / 2) * X05.logit t) := by
    rcases le_total 0 (η * X05.logit t) with h | h
    · rw [abs_of_nonneg h]
      have : Real.exp (η * X05.logit t / 2) = Real.exp (η / 2 * X05.logit t) := by ring_nf
      rw [this]; linarith [Real.exp_pos (-(η / 2) * X05.logit t)]
    · rw [abs_of_nonpos h]
      have : Real.exp (-(η * X05.logit t) / 2) = Real.exp (-(η / 2) * X05.logit t) := by ring_nf
      rw [this]; linarith [Real.exp_pos (η / 2 * X05.logit t)]
  have hl : 0 ≤ |X05.logit t| ^ k := by positivity
  calc |X05.logit t| ^ k * ((t * (1 - t)) ^ (η / 2) * Real.exp (|η * X05.logit t| / 2))
      ≤ |X05.logit t| ^ k * ((t * (1 - t)) ^ (η / 2) *
          (Real.exp (η / 2 * X05.logit t) + Real.exp (-(η / 2) * X05.logit t))) :=
        mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left he hw) hl
    _ = (t * (1 - t)) ^ (η / 2) * |X05.logit t| ^ k * Real.exp (η / 2 * X05.logit t) +
        (t * (1 - t)) ^ (η / 2) * |X05.logit t| ^ k * Real.exp (-(η / 2) * X05.logit t) := by
        ring
    _ ≤ _ := by linarith

end Q01

namespace Q01

open MeasureTheory Set Filter Topology

lemma q_abs_le {A B : ℝ} (hA : 0 < A) (hB : 0 < B) : |X01.q A B| ≤ 1 := by
  unfold X01.q
  have h : 0 < A ^ 2 + B ^ 2 := by positivity
  rw [abs_div, abs_of_pos h, div_le_one h, abs_le]
  constructor <;> nlinarith [sq_nonneg A, sq_nonneg B]

lemma I_repr {A B : ℝ} (hA : 0 < A) (hB : 0 < B) {η : ℝ} (hη : |η| ≤ 1 / 4) :
    X01.I A B η = (1 + η) *
      ∫ t in (0 : ℝ)..1, (t * (1 - t)) ^ (η / 2) * gf (X01.q A B) (η * X05.logit t) := by
  have h1 : 0 < 1 + η := by have := neg_abs_le η; linarith
  have hp : 0 < 1 / (1 + η) := by positivity
  have h := X01.proof_4 A B hA hB (1 / (1 + η)) hp
  have he : 1 / (1 / (1 + η)) - 1 = η := by field_simp; ring
  rw [he] at h
  have hr := X01.r_pos (B := B) hA
  unfold X01.I
  rw [h]
  field_simp
  rfl

lemma w_gf_int {q η : ℝ} (hq : |q| ≤ 1) (hη : |η| ≤ 1 / 4) :
    IntervalIntegrable (fun t => (t * (1 - t)) ^ (η / 2) * gf q (η * X05.logit t)) volume 0 1 := by
  refine IntervalIntegrable.mono_fun' (X05.q_int.const_mul (2 * Kc 0))
    (by unfold gf X05.logit; fun_prop) ?_
  rw [uIoc_of_le zero_le_one, ← restrict_Ioo_eq_restrict_Ioc]
  refine ae_restrict_of_forall_mem measurableSet_Ioo fun t ht => ?_
  have h1 : 0 < 1 - t := by linarith [ht.2]
  have hw : 0 ≤ (t * (1 - t)) ^ (η / 2) := by have := ht.1; positivity
  have hb := sqrtcosh_bd' hη ht.1 ht.2 0
  rw [pow_zero, one_mul] at hb
  have hg := gf_le hq (η * X05.logit t)
  dsimp only
  rw [Real.norm_eq_abs, abs_mul, abs_of_nonneg hw, abs_of_nonneg (gf_nonneg _ _)]
  calc (t * (1 - t)) ^ (η / 2) * gf q (η * X05.logit t)
      ≤ (t * (1 - t)) ^ (η / 2) * Real.exp (|η * X05.logit t| / 2) :=
        mul_le_mul_of_nonneg_left hg hw
    _ ≤ 2 * (Kc 0 * X05.q t) := hb
    _ = 2 * Kc 0 * X05.q t := by ring

lemma w_pow_int {η : ℝ} (hη : |η| ≤ 1 / 4) (k : ℕ) :
    IntervalIntegrable (fun t => (t * (1 - t)) ^ (η / 2) * X05.logit t ^ k) volume 0 1 := by
  have := G_int' hη (by norm_num : |(0 : ℝ)| ≤ 1 / 8) k
  simpa using this

end Q01

theorem Q01.proof_5 : Q01.stmt_5 := by
  intro A B hA hB M hM
  obtain ⟨C3, hC3⟩ := Q01.proof_3 M hM
  set q := X01.q A B with hqdef
  have hq : |q| ≤ 1 := Q01.q_abs_le hA hB
  set J := ∫ t in (0 : ℝ)..1, X05.q t with hJ
  refine ⟨5 / 4 * (|C3| * (2 * Q01.Kc M) * J), fun η hη => ?_⟩
  have h1η : 0 < 1 + η := by have := neg_abs_le η; linarith
  have h1η' : |1 + η| ≤ 5 / 4 := by
    rw [abs_of_pos h1η]; linarith [le_abs_self η]
  -- the odd terms vanish
  have hsum : ∑ m ∈ (Finset.range M).filter Even,
        Q01.gcoef m q * η ^ m * iteratedDeriv m (fun s => X05.A η s) 0 =
      ∑ m ∈ Finset.range M, Q01.gcoef m q * η ^ m * iteratedDeriv m (fun s => X05.A η s) 0 := by
    refine Finset.sum_filter_of_ne fun m _ hne => ?_
    by_contra hodd
    rw [Nat.not_even_iff_odd] at hodd
    exact hne (by rw [Q01.iter_A_odd η hodd, mul_zero])
  rw [hsum]
  simp only [Q01.iter_A0 hη]
  have hfs : ∑ m ∈ Finset.range M, Q01.gcoef m q * η ^ m *
        ((1 + η) * ∫ t in (0 : ℝ)..1, (t * (1 - t)) ^ (η / 2) * X05.logit t ^ m) =
      (1 + η) * ∫ t in (0 : ℝ)..1, ∑ m ∈ Finset.range M,
        Q01.gcoef m q * η ^ m * ((t * (1 - t)) ^ (η / 2) * X05.logit t ^ m) := by
    rw [intervalIntegral.integral_finsetSum (fun m _ => (Q01.w_pow_int hη m).const_mul _),
      Finset.mul_sum]
    refine Finset.sum_congr rfl fun m _ => ?_
    rw [intervalIntegral.integral_const_mul]; ring
  have hsi : IntervalIntegrable (fun t => ∑ m ∈ Finset.range M,
      Q01.gcoef m q * η ^ m * ((t * (1 - t)) ^ (η / 2) * X05.logit t ^ m)) MeasureTheory.volume 0 1 := by
    have := IntervalIntegrable.sum (Finset.range M)
      (fun m _ => (Q01.w_pow_int hη m).const_mul (Q01.gcoef m q * η ^ m))
    convert this using 1
    funext t; simp [Finset.sum_apply]
  rw [hfs, Q01.I_repr hA hB hη, ← hqdef, ← mul_sub, ← intervalIntegral.integral_sub
    (Q01.w_gf_int hq hη) hsi]
  rw [abs_mul]
  have hbd : ‖∫ t in (0 : ℝ)..1, ((t * (1 - t)) ^ (η / 2) * Q01.gf q (η * X05.logit t) -
        ∑ m ∈ Finset.range M, Q01.gcoef m q * η ^ m * ((t * (1 - t)) ^ (η / 2) * X05.logit t ^ m))‖
      ≤ ∫ t in (0 : ℝ)..1, |C3| * |η| ^ M * (2 * Q01.Kc M) * X05.q t := by
    refine intervalIntegral.norm_integral_le_of_norm_le zero_le_one ?_
      (X05.q_int.const_mul _)
    filter_upwards [X05.ae_ne_one] with t ht1 ht
    have ht0 : 0 < t := ht.1
    have ht1' : t < 1 := lt_of_le_of_ne ht.2 ht1
    have h1 : 0 < 1 - t := by linarith
    have hw : 0 ≤ (t * (1 - t)) ^ (η / 2) := by positivity
    set w := (t * (1 - t)) ^ (η / 2) with hwdef
    have e : w * Q01.gf q (η * X05.logit t) -
        ∑ m ∈ Finset.range M, Q01.gcoef m q * η ^ m * (w * X05.logit t ^ m)
        = w * (Q01.gf q (η * X05.logit t) -
          ∑ m ∈ Finset.range M, Q01.gcoef m q * (η * X05.logit t) ^ m) := by
      rw [mul_sub, Finset.mul_sum]
      congr 1
      refine Finset.sum_congr rfl fun m _ => ?_
      ring
    rw [e, Real.norm_eq_abs, abs_mul, abs_of_nonneg hw]
    have hr := hC3 q hq (η * X05.logit t)
    have hb := Q01.sqrtcosh_bd' hη ht0 ht1' M
    have hC : C3 ≤ |C3| := le_abs_self C3
    have hz : 0 ≤ |η * X05.logit t| ^ M * Real.exp (|η * X05.logit t| / 2) := by positivity
    calc w * |Q01.gf q (η * X05.logit t) -
          ∑ m ∈ Finset.range M, Q01.gcoef m q * (η * X05.logit t) ^ m|
        ≤ w * (|C3| * |η * X05.logit t| ^ M *
            Real.exp (|η * X05.logit t| / 2)) := by
          refine mul_le_mul_of_nonneg_left (hr.trans ?_) hw
          rw [mul_assoc, mul_assoc]
          exact mul_le_mul_of_nonneg_right hC hz
      _ = |C3| * |η| ^ M * (|X05.logit t| ^ M * (w *
            Real.exp (|η * X05.logit t| / 2))) := by
          rw [abs_mul, mul_pow]; ring
      _ ≤ |C3| * |η| ^ M * (2 * (Q01.Kc M * X05.q t)) := by
          gcongr
      _ = _ := by ring
  have hJ' : ∫ t in (0 : ℝ)..1, |C3| * |η| ^ M * (2 * Q01.Kc M) * X05.q t =
      |C3| * |η| ^ M * (2 * Q01.Kc M) * J := by
    rw [intervalIntegral.integral_const_mul]
  rw [hJ'] at hbd
  have hJ0 : 0 ≤ J := by
    rw [hJ]
    refine intervalIntegral.integral_nonneg zero_le_one fun t ht => ?_
    unfold X05.q; exact Real.rpow_nonneg (mul_nonneg ht.1 (by linarith [ht.2])) _
  rw [Real.norm_eq_abs] at hbd
  have hK := Q01.Kc_pos M
  calc |1 + η| * _ ≤ 5 / 4 * (|C3| * |η| ^ M * (2 * Q01.Kc M) * J) :=
        mul_le_mul h1η' hbd (abs_nonneg _) (by norm_num)
    _ = _ := by ring

namespace Q01

open Filter Topology

/-! ### `IsZetaPolyQ2` closure -/

lemma zq_zero (k : ℕ) : IsZetaPolyQ2 k (fun _ => 0) := ⟨[], by simp, by simp⟩

lemma zq_of_eq {k : ℕ} {f g : ℝ → ℝ} (h : ∀ x, f x = g x) (hf : IsZetaPolyQ2 k f) :
    IsZetaPolyQ2 k g := by
  have : f = g := funext h
  exact this ▸ hf

lemma zq_cast {k k' : ℕ} {f : ℝ → ℝ} (h : k = k') (hf : IsZetaPolyQ2 k f) :
    IsZetaPolyQ2 k' f := h ▸ hf

lemma zq_add {k : ℕ} {f g : ℝ → ℝ} (hf : IsZetaPolyQ2 k f) (hg : IsZetaPolyQ2 k g) :
    IsZetaPolyQ2 k (fun x => f x + g x) := by
  obtain ⟨L1, h1, e1⟩ := hf
  obtain ⟨L2, h2, e2⟩ := hg
  refine ⟨L1 ++ L2, fun e he => ?_, fun x => by simp [e1, e2]⟩
  rcases List.mem_append.1 he with h | h
  · exact h1 e h
  · exact h2 e h

lemma zq_smul {k : ℕ} {f : ℝ → ℝ} (c : ℚ) (hf : IsZetaPolyQ2 k f) :
    IsZetaPolyQ2 k (fun x => (c : ℝ) * f x) := by
  obtain ⟨L, h1, e1⟩ := hf
  refine ⟨L.map fun e => (c * e.1, e.2), fun e he => ?_, fun x => ?_⟩
  · obtain ⟨a, ha, rfl⟩ := List.mem_map.1 he
    exact h1 a ha
  · show (c : ℝ) * f x = _
    rw [e1, ← List.sum_map_mul_left]
    simp only [List.map_map]
    congr 1
    apply List.map_congr_left
    intro e _
    simp only [Function.comp, Rat.cast_mul]; ring

lemma zq_sum {k : ℕ} {ι : Type*} (s : Finset ι) (f : ι → ℝ → ℝ)
    (hf : ∀ i ∈ s, IsZetaPolyQ2 k (f i)) : IsZetaPolyQ2 k (fun x => ∑ i ∈ s, f i x) := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using zq_zero k
  | insert i s hi ih =>
    simp only [Finset.sum_insert hi]
    exact zq_add (hf i (Finset.mem_insert_self _ _))
      (ih fun j hj => hf j (Finset.mem_insert_of_mem hj))

lemma zq_mono (c : ℚ) (j : ℕ) : IsZetaPolyQ2 0 (fun x => (c : ℝ) * x ^ (2 * j)) :=
  ⟨[(c, j, [])], by simp, fun x => by simp⟩

lemma zq_const {k : ℕ} {y : ℝ} (hy : IsZetaPoly k y) : IsZetaPolyQ2 k (fun _ => y) := by
  obtain ⟨L, h1, rfl⟩ := hy
  refine ⟨L.map fun e => (e.1, 0, e.2), fun e he => ?_, fun x => ?_⟩
  · obtain ⟨a, ha, rfl⟩ := List.mem_map.1 he
    exact h1 a ha
  · simp only [List.map_map]
    congr 1
    apply List.map_congr_left
    intro e _
    simp

lemma zq_mono_mul {b : ℕ} {g : ℝ → ℝ} (c : ℚ) (j : ℕ) (ns : List ℕ) (hns : ∀ n ∈ ns, 2 ≤ n)
    (hg : IsZetaPolyQ2 b g) :
    IsZetaPolyQ2 (ns.sum + b) (fun x => (c : ℝ) * x ^ (2 * j) * (ns.map X05.z).prod * g x) := by
  obtain ⟨L, h1, e1⟩ := hg
  refine ⟨L.map fun e => (c * e.1, j + e.2.1, ns ++ e.2.2), fun e he => ?_, fun x => ?_⟩
  · obtain ⟨a, ha, rfl⟩ := List.mem_map.1 he
    refine ⟨fun n hn => ?_, ?_⟩
    · rcases List.mem_append.1 hn with h | h
      · exact hns n h
      · exact (h1 a ha).1 n h
    · simp [(h1 a ha).2]
  · show (c : ℝ) * x ^ (2 * j) * (ns.map X05.z).prod * g x = _
    rw [e1, ← List.sum_map_mul_left]
    simp only [List.map_map]
    congr 1
    apply List.map_congr_left
    intro e _
    simp only [Function.comp, Rat.cast_mul, List.map_append, List.prod_append, mul_add, pow_add]
    ring

lemma zq_mul {a b : ℕ} {f g : ℝ → ℝ} (hf : IsZetaPolyQ2 a f) (hg : IsZetaPolyQ2 b g) :
    IsZetaPolyQ2 (a + b) (fun x => f x * g x) := by
  obtain ⟨L, h1, e1⟩ := hf
  refine zq_of_eq (fun x => by rw [e1 x]) ?_
  clear e1
  induction L with
  | nil => simpa using zq_zero (a + b)
  | cons e L ih =>
    simp only [List.map_cons, List.sum_cons, add_mul]
    refine zq_add ?_ (ih fun e' he' => h1 e' (List.mem_cons_of_mem _ he'))
    have := zq_mono_mul e.1 e.2.1 e.2.2 (h1 e (List.mem_cons_self)).1 hg
    rwa [(h1 e (List.mem_cons_self)).2] at this

/-- even `g_m` lie in `ℚ[q²]`. -/
lemma gcoef_q2 {m : ℕ} (hm : Even m) : IsZetaPolyQ2 0 (gcoef m) := by
  obtain ⟨Q, hQ⟩ := (Q01.proof_2).2.1 m
  have hpar := (Q01.proof_2).2.2.1 m
  have key : ∀ x, gcoef m x = ∑ i ∈ Finset.range (Q.natDegree + 1),
      (if Even i then ((Q.coeff i : ℚ) : ℝ) * x ^ (2 * (i / 2)) else 0) := by
    intro x
    have h1 : gcoef m x = (gcoef m x + gcoef m (-x)) / 2 := by
      rw [hpar, hm.neg_one_pow]; ring
    rw [h1, hQ x, hQ (-x), Polynomial.aeval_eq_sum_range, Polynomial.aeval_eq_sum_range,
      ← Finset.sum_add_distrib, Finset.sum_div]
    refine Finset.sum_congr rfl fun i _ => ?_
    split_ifs with hi
    · obtain ⟨r, rfl⟩ := hi
      rw [show (r + r) / 2 = r by omega, show r + r = 2 * r by ring, Even.neg_pow ⟨r, by ring⟩]
      simp [Algebra.smul_def]
    · rw [Nat.not_even_iff_odd] at hi
      rw [hi.neg_pow]
      simp [Algebra.smul_def]
  refine zq_of_eq (fun x => (key x).symm) (zq_sum _ _ fun i _ => ?_)
  by_cases hi : Even i
  · simpa [hi] using zq_mono (Q.coeff i) (i / 2)
  · simpa [hi] using zq_zero 0

/-- weight one is empty: `IsZetaPoly 1 x → x = 0`. -/
lemma zp_one_eq_zero {x : ℝ} (hx : IsZetaPoly 1 x) : x = 0 := by
  obtain ⟨L, h1, rfl⟩ := hx
  have : ∀ e ∈ L, (e.1 : ℝ) * (e.2.map X05.z).prod = 0 := by
    intro e he
    exfalso
    obtain ⟨h2, hs⟩ := h1 e he
    rcases e with ⟨c, ns⟩
    cases ns with
    | nil => simp at hs
    | cons n ns =>
      have := h2 n (by simp)
      simp at hs; omega
  rw [List.map_congr_left this]; simp

/-- Taylor's theorem with `O`-remainder for an analytic germ. -/
lemma taylor_bd {f : ℝ → ℝ} (hf : AnalyticAt ℝ f 0) (n : ℕ) :
    ∃ δ > 0, ∃ D : ℝ, ∀ x : ℝ, |x| < δ →
      |f x - ∑ k ∈ Finset.range n, iteratedDeriv k f 0 / (k.factorial : ℝ) * x ^ k| ≤ D * |x| ^ n := by
  have h := (hf.hasFPowerSeriesAt.isBigO_sub_partialSum_pow n).bound
  obtain ⟨D, hD⟩ := h
  obtain ⟨δ, hδ, hball⟩ := Metric.eventually_nhds_iff.1 hD
  refine ⟨δ, hδ, D, fun x hx => ?_⟩
  have := hball (by simpa [Real.dist_eq] using hx)
  simp only [zero_add, FormalMultilinearSeries.partialSum,
    FormalMultilinearSeries.ofScalars_apply_eq, smul_eq_mul, Real.norm_eq_abs, abs_pow,
    abs_abs] at this
  exact this

/-- `β(m,a)`, the `η^a` coefficient of `∂_s^m A(η,0)`. -/
noncomputable def beta (m a : ℕ) : ℝ :=
  iteratedDeriv a (fun η => iteratedDeriv m (fun s => X05.A η s) 0) 0 / (a.factorial : ℝ)

lemma beta_pure (m a : ℕ) : IsZetaPoly (a + m) (beta m a) := ((Q01.proof_4).2 m).2.1 a

lemma beta_odd {m : ℕ} (hm : Odd m) (a : ℕ) : beta m a = 0 := by
  unfold beta
  have : (fun η => iteratedDeriv m (fun s => X05.A η s) 0) = fun _ => 0 :=
    funext fun η => iter_A_odd η hm
  rw [this]; simp

/-- `ĉ_N(q)` -/
noncomputable def chat (N : ℕ) (x : ℝ) : ℝ :=
  (-1 : ℝ) ^ N * ∑ m ∈ Finset.range (N + 1), gcoef m x * beta m (N - m)

lemma chat_pure (N : ℕ) : IsZetaPolyQ2 N (chat N) := by
  have h : IsZetaPolyQ2 N (fun x => (((-1) ^ N : ℚ) : ℝ) *
      ∑ m ∈ Finset.range (N + 1), gcoef m x * beta m (N - m)) := by
    refine zq_smul _ (zq_sum _ _ fun m hm => ?_)
    have hm' : m ≤ N := Nat.lt_succ_iff.1 (Finset.mem_range.1 hm)
    rcases Nat.even_or_odd m with he | ho
    · exact zq_cast (by omega) (zq_mul (gcoef_q2 he) (zq_const (beta_pure m (N - m))))
    · simpa [beta_odd ho] using zq_zero N
  exact zq_of_eq (fun x => by simp [chat]) h

lemma A00 : X05.A 0 0 = 1 := by simp [X05.A]

lemma chat_zero (x : ℝ) : chat 0 x = 1 := by
  simp [chat, beta, gcoef_eq, Qp, A00]

lemma chat_one (x : ℝ) : chat 1 x = 0 := by
  have h0 : beta 0 1 = 0 := zp_one_eq_zero (by simpa using beta_pure 0 1)
  have h1 : beta 1 0 = 0 := beta_odd odd_one 0
  simp [chat, Finset.sum_range_succ, h0, h1]

end Q01

lemma Q01.chat_asymp : ∀ A B : ℝ, 0 < A → 0 < B → ∀ K : ℕ, ∃ C τ₀ : ℝ, 0 < τ₀ ∧
    ∀ p : ℝ, 1 < p → Q01.tau p ≤ τ₀ →
      |X01.P A B p / (4 * X01.r A B) -
          ∑ N ∈ Finset.range (K + 1), Q01.chat N (X01.q A B) * Q01.tau p ^ N|
        ≤ C * Q01.tau p ^ (K + 1) := by
  intro A B hA hB K
  set q := X01.q A B with hqdef
  obtain ⟨C5, hC5⟩ := Q01.proof_5 A B hA hB (K + 1) (by omega)
  have hT := fun m => Q01.taylor_bd ((Q01.proof_4).2 m).1 (K + 1 - m)
  choose δf hδf Df hDf using hT
  have hne : (Finset.range (K + 1)).Nonempty := ⟨0, by simp⟩
  set δ := (Finset.range (K + 1)).inf' hne δf with hδ
  have hδpos : 0 < δ := by
    rw [hδ, Finset.lt_inf'_iff]; exact fun m _ => hδf m
  set C := C5 + ∑ m ∈ Finset.range (K + 1), |Q01.gcoef m q| * |Df m| with hC
  refine ⟨C, min (1 / 4) (δ / 2), lt_min (by norm_num) (by positivity), fun p hp hτ => ?_⟩
  set τ := Q01.tau p with hτdef
  have hp0 : 0 < p := by linarith
  have hτpos : 0 < τ := by
    rw [hτdef, Q01.tau]
    have : 1 / p < 1 := by rw [div_lt_one hp0]; exact hp
    linarith
  set η := -τ with hη
  have hηabs : |η| = τ := by rw [hη, abs_neg, abs_of_pos hτpos]
  have hη4 : |η| ≤ 1 / 4 := by rw [hηabs]; exact hτ.trans (min_le_left _ _)
  have hηδ : ∀ m ∈ Finset.range (K + 1), |η| < δf m := by
    intro m hm
    have : δ ≤ δf m := Finset.inf'_le _ hm
    rw [hηabs]
    have := hτ.trans (min_le_right _ _)
    linarith
  -- `P/(4r) = I(η)`
  have hI : X01.P A B p / (4 * X01.r A B) = X01.I A B η := by
    unfold X01.I
    congr 2
    rw [hη, hτdef, Q01.tau]; field_simp; ring
  -- the coefficient sum
  have hsum : ∑ N ∈ Finset.range (K + 1), Q01.chat N q * τ ^ N =
      ∑ m ∈ Finset.range (K + 1), Q01.gcoef m q * η ^ m *
        ∑ a ∈ Finset.range (K + 1 - m), Q01.beta m a * η ^ a := by
    have h1 : ∀ N, Q01.chat N q * τ ^ N =
        ∑ k ∈ Finset.range (N + 1), Q01.gcoef k q * Q01.beta k (N - k) * η ^ (k + (N - k)) := by
      intro N
      unfold Q01.chat
      rw [mul_comm, ← mul_assoc, Finset.mul_sum]
      refine Finset.sum_congr rfl fun k hk => ?_
      have hk' : k ≤ N := Nat.lt_succ_iff.1 (Finset.mem_range.1 hk)
      rw [Nat.add_sub_cancel' hk', hη, neg_pow]
      ring
    simp only [h1]
    rw [Finset.sum_range_diag_flip (K + 1)
      (fun k j => Q01.gcoef k q * Q01.beta k j * η ^ (k + j))]
    refine Finset.sum_congr rfl fun m _ => ?_
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun a _ => ?_
    ring
  -- odd terms of the `lem:pw:asymp` sum vanish
  have hfilt : ∑ m ∈ (Finset.range (K + 1)).filter Even,
        Q01.gcoef m q * η ^ m * iteratedDeriv m (fun s => X05.A η s) 0 =
      ∑ m ∈ Finset.range (K + 1),
        Q01.gcoef m q * η ^ m * iteratedDeriv m (fun s => X05.A η s) 0 := by
    refine Finset.sum_filter_of_ne fun m _ hne => ?_
    by_contra hodd
    rw [Nat.not_even_iff_odd] at hodd
    exact hne (by rw [Q01.iter_A_odd η hodd, mul_zero])
  have h5 := hC5 η hη4
  rw [hfilt] at h5
  rw [hI, hsum]
  have hsplit : X01.I A B η - ∑ m ∈ Finset.range (K + 1), Q01.gcoef m q * η ^ m *
        ∑ a ∈ Finset.range (K + 1 - m), Q01.beta m a * η ^ a =
      (X01.I A B η - ∑ m ∈ Finset.range (K + 1),
        Q01.gcoef m q * η ^ m * iteratedDeriv m (fun s => X05.A η s) 0) +
      ∑ m ∈ Finset.range (K + 1), Q01.gcoef m q * η ^ m *
        (iteratedDeriv m (fun s => X05.A η s) 0 -
          ∑ a ∈ Finset.range (K + 1 - m), Q01.beta m a * η ^ a) := by
    have : ∑ m ∈ Finset.range (K + 1), Q01.gcoef m q * η ^ m *
        (iteratedDeriv m (fun s => X05.A η s) 0 -
          ∑ a ∈ Finset.range (K + 1 - m), Q01.beta m a * η ^ a) =
        ∑ m ∈ Finset.range (K + 1),
          Q01.gcoef m q * η ^ m * iteratedDeriv m (fun s => X05.A η s) 0 -
        ∑ m ∈ Finset.range (K + 1), Q01.gcoef m q * η ^ m *
          ∑ a ∈ Finset.range (K + 1 - m), Q01.beta m a * η ^ a := by
      rw [← Finset.sum_sub_distrib]
      refine Finset.sum_congr rfl fun m _ => ?_
      ring
    rw [this]; ring
  rw [hsplit]
  have h2 : |∑ m ∈ Finset.range (K + 1), Q01.gcoef m q * η ^ m *
        (iteratedDeriv m (fun s => X05.A η s) 0 -
          ∑ a ∈ Finset.range (K + 1 - m), Q01.beta m a * η ^ a)| ≤
      (∑ m ∈ Finset.range (K + 1), |Q01.gcoef m q| * |Df m|) * τ ^ (K + 1) := by
    refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
    rw [Finset.sum_mul]
    refine Finset.sum_le_sum fun m hm => ?_
    have hm' : m ≤ K := Nat.lt_succ_iff.1 (Finset.mem_range.1 hm)
    have hD := hDf m η (hηδ m hm)
    have hb : iteratedDeriv m (fun s => X05.A η s) 0 -
        ∑ a ∈ Finset.range (K + 1 - m), Q01.beta m a * η ^ a =
        iteratedDeriv m (fun s => X05.A η s) 0 - ∑ k ∈ Finset.range (K + 1 - m),
          iteratedDeriv k (fun η => iteratedDeriv m (fun s => X05.A η s) 0) 0 /
            (k.factorial : ℝ) * η ^ k := rfl
    rw [hb, abs_mul, abs_mul, abs_pow, hηabs]
    have hDD : Df m * |η| ^ (K + 1 - m) ≤ |Df m| * τ ^ (K + 1 - m) := by
      rw [hηabs]; exact mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity)
    calc |Q01.gcoef m q| * τ ^ m * |_| ≤ |Q01.gcoef m q| * τ ^ m * (|Df m| * τ ^ (K + 1 - m)) :=
          mul_le_mul_of_nonneg_left (hD.trans hDD) (by positivity)
      _ = |Q01.gcoef m q| * |Df m| * (τ ^ m * τ ^ (K + 1 - m)) := by ring
      _ = |Q01.gcoef m q| * |Df m| * τ ^ (K + 1) := by
          rw [← pow_add, Nat.add_sub_cancel' (by omega)]
  rw [hηabs] at h5
  calc _ ≤ |X01.I A B η - ∑ m ∈ Finset.range (K + 1),
        Q01.gcoef m q * η ^ m * iteratedDeriv m (fun s => X05.A η s) 0| + |_| := abs_add_le _ _
    _ ≤ C5 * τ ^ (K + 1) + (∑ m ∈ Finset.range (K + 1), |Q01.gcoef m q| * |Df m|) *
          τ ^ (K + 1) := add_le_add h5 h2
    _ = C * τ ^ (K + 1) := by rw [hC]; ring

theorem Q01.proof_6 : Q01.stmt_6 :=
  ⟨Q01.chat, Q01.chat_zero, Q01.chat_one, Q01.chat_pure, Q01.chat_asymp⟩

namespace Q01

open Filter Topology

/-- logarithmic coefficients: `w_0 = 0`, `n w_n = n c_n - Σ_{k=1}^{n-1} k w_k c_{n-k}`. -/
noncomputable def wlog (c : ℕ → ℝ) (n : ℕ) : ℝ :=
  if n = 0 then 0 else
    ((n : ℝ) * c n - ∑ k : Fin n, if (k : ℕ) = 0 then 0 else
      ((k : ℕ) : ℝ) * wlog c k * c (n - k)) / n
decreasing_by exact k.isLt

lemma wlog_zero (c : ℕ → ℝ) : wlog c 0 = 0 := by rw [wlog]; simp

lemma wlog_eq (c : ℕ → ℝ) {n : ℕ} (hn : 1 ≤ n) :
    wlog c n = ((n : ℝ) * c n - ∑ k ∈ Finset.Ico 1 n, (k : ℝ) * wlog c k * c (n - k)) / n := by
  have hfin : ∑ k : Fin n, (if (k : ℕ) = 0 then 0 else ((k : ℕ) : ℝ) * wlog c k * c (n - k)) =
      ∑ k ∈ Finset.Ico 1 n, (k : ℝ) * wlog c k * c (n - k) := by
    rw [Fin.sum_univ_eq_sum_range (fun k => if k = 0 then 0 else (k : ℝ) * wlog c k * c (n - k)),
      Finset.range_eq_Ico, Finset.sum_eq_sum_Ico_succ_bot (by omega)]
    simp only [↓reduceIte, zero_add]
    exact Finset.sum_congr rfl fun k hk => by
      have := (Finset.mem_Ico.1 hk).1
      simp [show k ≠ 0 by omega]
  rw [wlog]; simp only [show n ≠ 0 from by omega, ↓reduceIte]; rw [hfin]

lemma wlog_rec (c : ℕ → ℝ) (hc0 : c 0 = 1) (n : ℕ) (hn : 1 ≤ n) :
    (n : ℝ) * c n = ∑ k ∈ Finset.Icc 1 n, (k : ℝ) * wlog c k * c (n - k) := by
  have hn' : (n : ℝ) ≠ 0 := by positivity
  rw [← Finset.Ico_add_one_right_eq_Icc, Finset.sum_Ico_succ_top hn, Nat.sub_self, hc0, mul_one,
    wlog_eq c hn]
  field_simp
  ring

lemma wlog_q2 (c : ℕ → ℝ → ℝ) (hc : ∀ N, IsZetaPolyQ2 N (c N)) :
    ∀ n, IsZetaPolyQ2 n (fun x => wlog (fun N => c N x) n) := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    rcases Nat.eq_zero_or_pos n with rfl | hn
    · simpa [wlog_zero] using zq_zero 0
    · have hn' : (n : ℝ) ≠ 0 := by positivity
      have key : ∀ x, ((1 / n : ℚ) : ℝ) * ((n : ℚ) * c n x + ((-1 : ℚ) : ℝ) *
          ∑ k ∈ Finset.Ico 1 n, ((k : ℚ) : ℝ) * (wlog (fun N => c N x) k * c (n - k) x)) =
          wlog (fun N => c N x) n := by
        intro x
        rw [wlog_eq _ hn]
        push_cast
        field_simp
        ring
      refine zq_of_eq key (zq_smul _ (zq_add (zq_smul _ (hc n)) (zq_smul _ ?_)))
      refine zq_sum _ _ fun k hk => ?_
      obtain ⟨hk1, hk2⟩ := Finset.mem_Ico.1 hk
      exact zq_smul _ (zq_cast (by omega) (zq_mul (ih k hk2) (hc (n - k))))

/-- the exp recursion for Taylor coefficients. -/
lemma exp_rec {f : ℝ → ℝ} (hf : AnalyticAt ℝ f 0) (n : ℕ) (hn : 1 ≤ n) :
    (n : ℝ) * (iteratedDeriv n (fun x => Real.exp (f x)) 0 / (n.factorial : ℝ)) =
      ∑ k ∈ Finset.Icc 1 n, (k : ℝ) * (iteratedDeriv k f 0 / (k.factorial : ℝ)) *
        (iteratedDeriv (n - k) (fun x => Real.exp (f x)) 0 / ((n - k).factorial : ℝ)) := by
  have hE : AnalyticAt ℝ (fun x => Real.exp (f x)) 0 := hf.rexp'
  have hder : deriv (fun x => Real.exp (f x)) =ᶠ[𝓝 0] fun x => deriv f x * Real.exp (f x) := by
    filter_upwards [hf.eventually_analyticAt] with x hx
    rw [(hx.differentiableAt.hasDerivAt.exp).deriv]; ring
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  have h1 : iteratedDeriv (m + 1) (fun x => Real.exp (f x)) 0 =
      ∑ i ∈ Finset.range (m + 1), (m.choose i : ℝ) * iteratedDeriv (i + 1) f 0 *
        iteratedDeriv (m - i) (fun x => Real.exp (f x)) 0 := by
    rw [iteratedDeriv_succ', hder.iteratedDeriv_eq]
    have := iteratedDeriv_mul (n := m) hf.deriv.contDiffAt hE.contDiffAt
    rw [show (fun x => deriv f x * Real.exp (f x)) = deriv f * (fun x => Real.exp (f x))
      from rfl, this]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [← iteratedDeriv_succ']
  have hIcc : Finset.Icc 1 (m + 1) = Finset.Ico 1 (m + 2) := rfl
  rw [hIcc, Finset.sum_Ico_eq_sum_range, show m + 2 - 1 = m + 1 by omega]
  simp only [h1, Finset.mul_sum, Finset.sum_div]
  refine Finset.sum_congr rfl fun i hi => ?_
  have hi' : i ≤ m := Nat.lt_succ_iff.1 (Finset.mem_range.1 hi)
  rw [show m + 1 - (1 + i) = m - i by omega, add_comm 1 i, Nat.cast_choose ℝ hi']
  rw [Nat.factorial_succ, Nat.factorial_succ]
  push_cast
  have h1 : (i.factorial : ℝ) ≠ 0 := by positivity
  have h2 : ((m - i).factorial : ℝ) ≠ 0 := by positivity
  have h3 : (m.factorial : ℝ) ≠ 0 := by positivity
  field_simp

lemma poly_iter (a : ℕ → ℝ) (n k : ℕ) :
    iteratedDeriv k (fun x => ∑ N ∈ Finset.range n, a N * x ^ N) 0 =
      if k < n then (k.factorial : ℝ) * a k else 0 := by
  induction n with
  | zero => simp
  | succ n ih =>
    simp only [Finset.sum_range_succ]
    have h1 : ContDiffAt ℝ k (fun x : ℝ => ∑ N ∈ Finset.range n, a N * x ^ N) 0 := by fun_prop
    have h2 : ContDiffAt ℝ k (fun x : ℝ => a n * x ^ n) 0 := by fun_prop
    have := iteratedDeriv_add (n := k) h1 h2
    rw [show (fun x => ∑ N ∈ Finset.range n, a N * x ^ N + a n * x ^ n) =
      (fun x : ℝ => ∑ N ∈ Finset.range n, a N * x ^ N) + (fun x : ℝ => a n * x ^ n) from rfl, this,
      ih, iteratedDeriv_const_mul_field, iteratedDeriv_fun_pow_zero]
    by_cases hk : k < n
    · simp [hk, show k < n + 1 by omega, show k ≠ n by omega]
    · by_cases hkn : k = n
      · subst hkn; simp; ring
      · simp [hk, hkn, show ¬ k < n + 1 by omega]

lemma log_lip {a b : ℝ} (ha : 1 / 2 ≤ a) (hb : 1 / 2 ≤ b) :
    |Real.log a - Real.log b| ≤ 2 * |a - b| := by
  have ha0 : 0 < a := by linarith
  have hb0 : 0 < b := by linarith
  rw [← Real.log_div ha0.ne' hb0.ne', abs_le]
  constructor
  · have h := Real.log_le_sub_one_of_pos (div_pos hb0 ha0)
    rw [Real.log_div hb0.ne' ha0.ne'] at h
    rw [Real.log_div ha0.ne' hb0.ne']
    have : b / a - 1 = (b - a) / a := by field_simp
    rw [this] at h
    have h2 : (b - a) / a ≤ 2 * |a - b| := by
      rw [div_le_iff₀ ha0]
      have := neg_abs_le (a - b)
      have := abs_nonneg (a - b)
      nlinarith [abs_sub_comm a b, le_abs_self (b - a)]
    linarith
  · have h := Real.log_le_sub_one_of_pos (div_pos ha0 hb0)
    have : a / b - 1 = (a - b) / b := by field_simp
    rw [this] at h
    have h2 : (a - b) / b ≤ 2 * |a - b| := by
      rw [div_le_iff₀ hb0]
      have := abs_nonneg (a - b)
      nlinarith [le_abs_self (a - b)]
    linarith

/-- `log` of an asymptotic series with constant term `1`. -/
lemma log_asymp (c : ℕ → ℝ) (hc0 : c 0 = 1) (K : ℕ) (C₀ : ℝ) :
    ∃ C τ₁ : ℝ, 0 < τ₁ ∧ ∀ τ F : ℝ, 0 < τ → τ ≤ τ₁ →
      |F - ∑ N ∈ Finset.range (K + 1), c N * τ ^ N| ≤ C₀ * τ ^ (K + 1) →
      |Real.log F - ∑ k ∈ Finset.range (K + 1), wlog c k * τ ^ k| ≤ C * τ ^ (K + 1) := by
  set S : ℝ → ℝ := fun τ => ∑ N ∈ Finset.range (K + 1), c N * τ ^ N with hSdef
  have hS0 : S 0 = 1 := by simp [hSdef, Finset.sum_range_succ', hc0]
  have hScont : Continuous S := by rw [hSdef]; fun_prop
  have hSan : AnalyticAt ℝ S 0 := by rw [hSdef]; fun_prop
  set L : ℝ → ℝ := fun τ => Real.log (S τ) with hLdef
  have hLan : AnalyticAt ℝ L 0 := hSan.log (by rw [hS0]; norm_num)
  have hpos : ∀ᶠ x in 𝓝 (0 : ℝ), 0 < S x :=
    continuousAt_const.eventually_lt hScont.continuousAt (by rw [hS0]; norm_num)
  have hexp : (fun x => Real.exp (L x)) =ᶠ[𝓝 0] S := by
    filter_upwards [hpos] with x hx
    simp [hLdef, Real.exp_log hx]
  have he : ∀ n, n ≤ K →
      iteratedDeriv n (fun x => Real.exp (L x)) 0 / (n.factorial : ℝ) = c n := by
    intro n hn
    rw [hexp.iteratedDeriv_eq, hSdef, poly_iter, ite_eq_left_of_eq_true _ _ (eq_true (by omega))]
    field_simp
  have hl : ∀ n, n ≤ K → iteratedDeriv n L 0 / (n.factorial : ℝ) = wlog c n := by
    intro n
    induction n using Nat.strong_induction_on with
    | _ n ih =>
      intro hnK
      rcases Nat.eq_zero_or_pos n with rfl | hn
      · simp [wlog_zero, hLdef, hS0]
      · have h1 := exp_rec hLan n hn
        have h2 := wlog_rec c hc0 n hn
        rw [← Finset.Ico_add_one_right_eq_Icc, Finset.sum_Ico_succ_top hn, Nat.sub_self] at h1 h2
        simp only [Nat.succ_eq_add_one, zero_add] at h1 h2
        rw [he n hnK, he 0 (by omega), hc0, mul_one] at h1
        rw [hc0, mul_one] at h2
        have hsum : ∑ k ∈ Finset.Ico 1 n, (k : ℝ) * (iteratedDeriv k L 0 / (k.factorial : ℝ)) *
            (iteratedDeriv (n - k) (fun x => Real.exp (L x)) 0 / ((n - k).factorial : ℝ)) =
            ∑ k ∈ Finset.Ico 1 n, (k : ℝ) * wlog c k * c (n - k) := by
          refine Finset.sum_congr rfl fun k hk => ?_
          obtain ⟨hk1, hk2⟩ := Finset.mem_Ico.1 hk
          rw [ih k hk2 (by omega), he (n - k) (by omega)]
        rw [hsum] at h1
        have hn' : (n : ℝ) ≠ 0 := by positivity
        have := h1.symm.trans h2
        have h3 : (n : ℝ) * (iteratedDeriv n L 0 / (n.factorial : ℝ)) = n * wlog c n := by linarith
        exact mul_left_cancel₀ hn' h3
  obtain ⟨δ, hδ, D, hD⟩ := taylor_bd hLan (K + 1)
  obtain ⟨δ2, hδ2, hS⟩ := Metric.continuousAt_iff.1 hScont.continuousAt (1 / 4) (by norm_num)
  refine ⟨2 * |C₀| + |D|, min (min (δ / 2) (δ2 / 2)) (min 1 (1 / (4 * (|C₀| + 1)))),
    by positivity, fun τ F hτ hτ1 hF => ?_⟩
  have hτδ : τ < δ := by
    have := hτ1.trans ((min_le_left _ _).trans (min_le_left _ _)); linarith
  have hτδ2 : τ < δ2 := by
    have := hτ1.trans ((min_le_left _ _).trans (min_le_right _ _)); linarith
  have hτ1' : τ ≤ 1 := hτ1.trans ((min_le_right _ _).trans (min_le_left _ _))
  have hτc : τ ≤ 1 / (4 * (|C₀| + 1)) := hτ1.trans ((min_le_right _ _).trans (min_le_right _ _))
  have hSτ : |S τ - 1| < 1 / 4 := by
    have := hS (x := τ) (by rw [Real.dist_eq, sub_zero, abs_of_pos hτ]; exact hτδ2)
    rwa [Real.dist_eq, hS0] at this
  have hpow : τ ^ (K + 1) ≤ τ := by
    calc τ ^ (K + 1) ≤ τ ^ 1 := pow_le_pow_of_le_one hτ.le hτ1' (by omega)
      _ = τ := pow_one τ
  have hF' : |F - S τ| ≤ |C₀| * τ ^ (K + 1) :=
    hF.trans (mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity))
  have hsmall : |C₀| * τ ^ (K + 1) ≤ 1 / 4 := by
    calc |C₀| * τ ^ (K + 1) ≤ |C₀| * τ := mul_le_mul_of_nonneg_left hpow (abs_nonneg _)
      _ ≤ |C₀| * (1 / (4 * (|C₀| + 1))) := mul_le_mul_of_nonneg_left hτc (abs_nonneg _)
      _ ≤ 1 / 4 := by
          rw [mul_one_div, div_le_div_iff₀ (by positivity) (by norm_num)]
          nlinarith [abs_nonneg C₀]
  have hSlo' : 3 / 4 < S τ := by rw [abs_lt] at hSτ; linarith
  have hSlo : 1 / 2 ≤ S τ := by linarith
  have hFlo : 1 / 2 ≤ F := by
    rw [abs_le] at hF'; linarith
  have hlip := log_lip hFlo hSlo
  have hT := hD τ (by rw [abs_of_pos hτ]; exact hτδ)
  have hsw : ∑ k ∈ Finset.range (K + 1), wlog c k * τ ^ k =
      ∑ k ∈ Finset.range (K + 1), iteratedDeriv k L 0 / (k.factorial : ℝ) * τ ^ k := by
    refine Finset.sum_congr rfl fun k hk => ?_
    rw [hl k (Nat.lt_succ_iff.1 (Finset.mem_range.1 hk))]
  rw [hsw]
  rw [abs_of_pos hτ] at hT
  have hT' : D * τ ^ (K + 1) ≤ |D| * τ ^ (K + 1) :=
    mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity)
  calc |Real.log F - ∑ k ∈ Finset.range (K + 1), iteratedDeriv k L 0 / (k.factorial : ℝ) * τ ^ k|
      ≤ |Real.log F - Real.log (S τ)| +
        |L τ - ∑ k ∈ Finset.range (K + 1), iteratedDeriv k L 0 / (k.factorial : ℝ) * τ ^ k| := by
        have := abs_add_le (Real.log F - Real.log (S τ))
          (L τ - ∑ k ∈ Finset.range (K + 1), iteratedDeriv k L 0 / (k.factorial : ℝ) * τ ^ k)
        simpa [hLdef] using this
    _ ≤ 2 * (|C₀| * τ ^ (K + 1)) + |D| * τ ^ (K + 1) := by
        have := mul_le_mul_of_nonneg_left hF' (by norm_num : (0:ℝ) ≤ 2)
        linarith
    _ = _ := by ring

end Q01

namespace Q01

open Filter Topology

/-- a polynomial that is `O(x^n)` as `x → 0⁺` vanishes. -/
lemma poly_zero (n : ℕ) : ∀ (d : ℕ → ℝ) (C δ : ℝ), 0 < δ →
    (∀ x : ℝ, 0 < x → x < δ → |∑ k ∈ Finset.range n, d k * x ^ k| ≤ C * x ^ n) →
    ∀ k < n, d k = 0 := by
  induction n with
  | zero => intro d C δ _ _ k hk; omega
  | succ n ih =>
    intro d C δ hδ h
    have hsplit : ∀ x : ℝ, ∑ k ∈ Finset.range (n + 1), d k * x ^ k =
        x * ∑ k ∈ Finset.range n, d (k + 1) * x ^ k + d 0 := by
      intro x
      rw [Finset.sum_range_succ', Finset.mul_sum]
      simp only [pow_zero, mul_one, pow_succ]
      congr 1
      refine Finset.sum_congr rfl fun k _ => by ring
    have hd0 : d 0 = 0 := by
      by_contra hne
      have hpos : 0 < |d 0| := abs_pos.2 hne
      set B := ∑ k ∈ Finset.range n, |d (k + 1)| with hB
      have hB0 : 0 ≤ B := Finset.sum_nonneg fun k _ => abs_nonneg _
      set x := min (δ / 2) (min (1 / 2) (|d 0| / (2 * (|C| + B + 1)))) with hx
      have hx0 : 0 < x := by positivity
      have hxδ : x < δ := by have := min_le_left (δ / 2) (min (1 / 2) (|d 0| / (2 * (|C| + B + 1)))); linarith
      have hx1 : x ≤ 1 / 2 := (min_le_right _ _).trans (min_le_left _ _)
      have hx2 : x ≤ |d 0| / (2 * (|C| + B + 1)) := (min_le_right _ _).trans (min_le_right _ _)
      have h1 := h x hx0 hxδ
      rw [hsplit] at h1
      have hrest : |∑ k ∈ Finset.range n, d (k + 1) * x ^ k| ≤ B := by
        refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun k _ => ?_)
        rw [abs_mul, abs_pow, abs_of_pos hx0]
        have : x ^ k ≤ 1 := pow_le_one₀ hx0.le (by linarith)
        nlinarith [abs_nonneg (d (k + 1))]
      have hCx : C * x ^ (n + 1) ≤ |C| * x := by
        have : x ^ (n + 1) ≤ x := by
          calc x ^ (n + 1) ≤ x ^ 1 := pow_le_pow_of_le_one hx0.le (by linarith) (by omega)
            _ = x := pow_one x
        calc C * x ^ (n + 1) ≤ |C| * x ^ (n + 1) :=
              mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity)
          _ ≤ |C| * x := mul_le_mul_of_nonneg_left this (abs_nonneg _)
      have hd : |d 0| ≤ (|C| + B) * x := by
        have e : d 0 = (x * ∑ k ∈ Finset.range n, d (k + 1) * x ^ k + d 0) -
            x * ∑ k ∈ Finset.range n, d (k + 1) * x ^ k := by ring
        rw [e]
        refine (abs_sub _ _).trans ?_
        rw [abs_mul, abs_of_pos hx0]
        nlinarith
      have : (|C| + B) * x < |d 0| := by
        calc (|C| + B) * x ≤ (|C| + B) * (|d 0| / (2 * (|C| + B + 1))) :=
              mul_le_mul_of_nonneg_left hx2 (by positivity)
          _ < |d 0| := by
              rw [mul_div_assoc', div_lt_iff₀ (by positivity)]
              nlinarith [abs_nonneg C]
      linarith
    have h' : ∀ x : ℝ, 0 < x → x < δ →
        |∑ k ∈ Finset.range n, d (k + 1) * x ^ k| ≤ C * x ^ n := by
      intro x hx hxδ
      have h1 := h x hx hxδ
      rw [hsplit, hd0, add_zero, abs_mul, abs_of_pos hx, pow_succ] at h1
      have : x * |∑ k ∈ Finset.range n, d (k + 1) * x ^ k| ≤ x * (C * x ^ n) := by linarith
      exact le_of_mul_le_mul_left this hx
    have ih' := ih (fun k => d (k + 1)) C δ hδ h'
    intro k hk
    rcases k with _ | k
    · exact hd0
    · exact ih' k (by omega)

lemma beta0_vals : beta 0 2 = -(X05.z 2) / 4 ∧ beta 0 3 = X05.z 3 / 4 := by
  obtain ⟨C, δ, hδ, h6⟩ := X05.proof_6
  obtain ⟨δ', hδ', D, hD⟩ := taylor_bd ((Q01.proof_4).2 0).1 4
  set a4 : ℝ := (Real.pi ^ 2 / 6) ^ 2 / 32 - 7 * (Real.pi ^ 4 / 90) / 32
  set a5 : ℝ := -(Real.pi ^ 2 / 6) * X05.z 3 / 16 + 3 * X05.z 5 / 16
  set d : ℕ → ℝ := fun k => beta 0 k -
    (if k = 0 then 1 else if k = 2 then -(Real.pi ^ 2 / 6) / 4 else if k = 3 then X05.z 3 / 4
      else 0) with hd
  have hz := poly_zero 4 d (|D| + |C| + |a4| + |a5|) (min (min δ δ') 1) (by positivity)
    (fun x hx hxδ => ?_)
  · have h2 := hz 2 (by norm_num)
    have h3 := hz 3 (by norm_num)
    simp only [hd] at h2 h3
    norm_num at h2 h3
    rw [X05.z_two]
    constructor <;> linarith
  · have hxδ1 : x ≤ δ := by have := min_le_left (min δ δ') 1; have := min_le_left δ δ'; linarith
    have hxδ2 : x < δ' := by
      have := min_le_left (min δ δ') 1; have := min_le_right δ δ'; linarith
    have hx1 : x ≤ 1 := by have := min_le_right (min δ δ') 1; linarith
    have hA := (h6 x (by rw [abs_of_pos hx]; exact hxδ1)).1
    have hT := hD x (by rw [abs_of_pos hx]; exact hxδ2)
    simp only [iteratedDeriv_zero] at hT
    rw [abs_of_pos hx] at hA hT
    have e : ∑ k ∈ Finset.range 4, d k * x ^ k =
        -(X05.A x 0 - ∑ k ∈ Finset.range 4,
          iteratedDeriv k (fun η => X05.A η 0) 0 / (k.factorial : ℝ) * x ^ k) +
        (X05.A x 0 - (1 - (Real.pi ^ 2 / 6) / 4 * x ^ 2 + X05.z 3 / 4 * x ^ 3 + a4 * x ^ 4 +
          a5 * x ^ 5)) + a4 * x ^ 4 + a5 * x ^ 5 := by
      simp only [hd, beta, Finset.sum_range_succ, Finset.sum_range_zero]
      norm_num
      ring
    rw [e]
    have h4 : x ^ 6 ≤ x ^ 4 := pow_le_pow_of_le_one hx.le hx1 (by norm_num)
    have h5 : x ^ 5 ≤ x ^ 4 := pow_le_pow_of_le_one hx.le hx1 (by norm_num)
    have hx4 : 0 ≤ x ^ 4 := by positivity
    have q1 := abs_add_le (-(X05.A x 0 - ∑ k ∈ Finset.range 4,
          iteratedDeriv k (fun η => X05.A η 0) 0 / (k.factorial : ℝ) * x ^ k) +
        (X05.A x 0 - (1 - (Real.pi ^ 2 / 6) / 4 * x ^ 2 + X05.z 3 / 4 * x ^ 3 + a4 * x ^ 4 +
          a5 * x ^ 5)) + a4 * x ^ 4) (a5 * x ^ 5)
    have q2 := abs_add_le (-(X05.A x 0 - ∑ k ∈ Finset.range 4,
          iteratedDeriv k (fun η => X05.A η 0) 0 / (k.factorial : ℝ) * x ^ k) +
        (X05.A x 0 - (1 - (Real.pi ^ 2 / 6) / 4 * x ^ 2 + X05.z 3 / 4 * x ^ 3 + a4 * x ^ 4 +
          a5 * x ^ 5))) (a4 * x ^ 4)
    have q3 := abs_add_le (-(X05.A x 0 - ∑ k ∈ Finset.range 4,
          iteratedDeriv k (fun η => X05.A η 0) 0 / (k.factorial : ℝ) * x ^ k))
        (X05.A x 0 - (1 - (Real.pi ^ 2 / 6) / 4 * x ^ 2 + X05.z 3 / 4 * x ^ 3 + a4 * x ^ 4 +
          a5 * x ^ 5))
    rw [abs_neg] at q3
    have r1 : |a4 * x ^ 4| = |a4| * x ^ 4 := by rw [abs_mul, abs_of_nonneg hx4]
    have r2 : |a5 * x ^ 5| = |a5| * x ^ 5 := by rw [abs_mul, abs_of_nonneg (pow_nonneg hx.le 5)]
    have s1 : D * x ^ 4 ≤ |D| * x ^ 4 := mul_le_mul_of_nonneg_right (le_abs_self _) hx4
    have s2 : C * x ^ 6 ≤ |C| * x ^ 4 :=
      (mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity)).trans
        (mul_le_mul_of_nonneg_left h4 (abs_nonneg _))
    have s3 : |a5| * x ^ 5 ≤ |a5| * x ^ 4 := mul_le_mul_of_nonneg_left h5 (abs_nonneg _)
    nlinarith

lemma beta2_vals : beta 2 0 = 2 * X05.z 2 ∧ beta 2 1 = -2 * X05.z 3 := by
  obtain ⟨C, δ, hδ, h6⟩ := X05.proof_6
  obtain ⟨δ', hδ', D, hD⟩ := taylor_bd ((Q01.proof_4).2 2).1 2
  set b2 : ℝ := -(Real.pi ^ 2 / 6) ^ 2 / 2 + 3 * (Real.pi ^ 4 / 90) / 2
  set b3 : ℝ := (Real.pi ^ 2 / 6) * X05.z 3 - X05.z 5
  set d : ℕ → ℝ := fun k => beta 2 k -
    (if k = 0 then 2 * (Real.pi ^ 2 / 6) else if k = 1 then -2 * X05.z 3 else 0) with hd
  have hz := poly_zero 2 d (|D| + |C| + |b2| + |b3|) (min (min δ δ') 1) (by positivity)
    (fun x hx hxδ => ?_)
  · have h0 := hz 0 (by norm_num)
    have h1 := hz 1 (by norm_num)
    simp only [hd] at h0 h1
    norm_num at h0 h1
    rw [X05.z_two]
    constructor <;> linarith
  · have hxδ1 : x ≤ δ := by have := min_le_left (min δ δ') 1; have := min_le_left δ δ'; linarith
    have hxδ2 : x < δ' := by
      have := min_le_left (min δ δ') 1; have := min_le_right δ δ'; linarith
    have hx1 : x ≤ 1 := by have := min_le_right (min δ δ') 1; linarith
    have hA := (h6 x (by rw [abs_of_pos hx]; exact hxδ1)).2.1
    have hT := hD x (by rw [abs_of_pos hx]; exact hxδ2)
    rw [abs_of_pos hx] at hA hT
    set B2 := iteratedDeriv 2 (fun s => X05.A x s) 0
    have e : ∑ k ∈ Finset.range 2, d k * x ^ k =
        -(B2 - ∑ k ∈ Finset.range 2,
          iteratedDeriv k (fun η => iteratedDeriv 2 (fun s => X05.A η s) 0) 0 /
            (k.factorial : ℝ) * x ^ k) +
        (B2 - (2 * (Real.pi ^ 2 / 6) - 2 * X05.z 3 * x + b2 * x ^ 2 + b3 * x ^ 3)) +
        b2 * x ^ 2 + b3 * x ^ 3 := by
      simp only [hd, beta, Finset.sum_range_succ, Finset.sum_range_zero]
      norm_num
      ring
    rw [e]
    have h4 : x ^ 4 ≤ x ^ 2 := pow_le_pow_of_le_one hx.le hx1 (by norm_num)
    have h3 : x ^ 3 ≤ x ^ 2 := pow_le_pow_of_le_one hx.le hx1 (by norm_num)
    have hx2 : 0 ≤ x ^ 2 := by positivity
    have q1 := abs_add_le (-(B2 - ∑ k ∈ Finset.range 2,
          iteratedDeriv k (fun η => iteratedDeriv 2 (fun s => X05.A η s) 0) 0 /
            (k.factorial : ℝ) * x ^ k) +
        (B2 - (2 * (Real.pi ^ 2 / 6) - 2 * X05.z 3 * x + b2 * x ^ 2 + b3 * x ^ 3)) +
        b2 * x ^ 2) (b3 * x ^ 3)
    have q2 := abs_add_le (-(B2 - ∑ k ∈ Finset.range 2,
          iteratedDeriv k (fun η => iteratedDeriv 2 (fun s => X05.A η s) 0) 0 /
            (k.factorial : ℝ) * x ^ k) +
        (B2 - (2 * (Real.pi ^ 2 / 6) - 2 * X05.z 3 * x + b2 * x ^ 2 + b3 * x ^ 3)))
        (b2 * x ^ 2)
    have q3 := abs_add_le (-(B2 - ∑ k ∈ Finset.range 2,
          iteratedDeriv k (fun η => iteratedDeriv 2 (fun s => X05.A η s) 0) 0 /
            (k.factorial : ℝ) * x ^ k))
        (B2 - (2 * (Real.pi ^ 2 / 6) - 2 * X05.z 3 * x + b2 * x ^ 2 + b3 * x ^ 3))
    rw [abs_neg] at q3
    have r1 : |b2 * x ^ 2| = |b2| * x ^ 2 := by rw [abs_mul, abs_of_nonneg hx2]
    have r2 : |b3 * x ^ 3| = |b3| * x ^ 3 := by rw [abs_mul, abs_of_nonneg (pow_nonneg hx.le 3)]
    have s1 : D * x ^ 2 ≤ |D| * x ^ 2 := mul_le_mul_of_nonneg_right (le_abs_self _) hx2
    have s2 : C * x ^ 4 ≤ |C| * x ^ 2 :=
      (mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity)).trans
        (mul_le_mul_of_nonneg_left h4 (abs_nonneg _))
    have s3 : |b3| * x ^ 3 ≤ |b3| * x ^ 2 := mul_le_mul_of_nonneg_left h3 (abs_nonneg _)
    nlinarith

lemma chat_two (x : ℝ) : chat 2 x = (1 - x ^ 2) * X05.z 2 / 4 := by
  obtain ⟨g0, g1, g2, g4⟩ := (Q01.proof_2).2.2.2 x
  simp only [chat, Finset.sum_range_succ, Finset.sum_range_zero, g0, g2,
    beta_odd odd_one, beta0_vals.1, beta2_vals.1]
  norm_num
  ring

lemma chat_three (x : ℝ) : chat 3 x = (1 - x ^ 2) * X05.z 3 / 4 := by
  obtain ⟨g0, g1, g2, g4⟩ := (Q01.proof_2).2.2.2 x
  simp only [chat, Finset.sum_range_succ, Finset.sum_range_zero, g0, g2,
    beta_odd odd_one, beta_odd (show Odd 3 by decide), beta0_vals.2, beta2_vals.2]
  norm_num
  ring

end Q01

namespace Q01

lemma wlog_one (x : ℝ) : wlog (fun N => chat N x) 1 = 0 := by
  rw [wlog_eq _ le_rfl]; simp [chat_one]

lemma wlog_two (x : ℝ) : wlog (fun N => chat N x) 2 = chat 2 x := by
  rw [wlog_eq _ (by norm_num), Finset.sum_Ico_succ_top (by norm_num), Finset.Ico_self]
  simp [wlog_one, chat_one]

lemma wlog_three (x : ℝ) : wlog (fun N => chat N x) 3 = chat 3 x := by
  rw [wlog_eq _ (by norm_num), Finset.sum_Ico_succ_top (by norm_num),
    Finset.sum_Ico_succ_top (by norm_num), Finset.Ico_self]
  simp [wlog_one, chat_one]

end Q01

theorem Q01.proof_7 : Q01.stmt_7 := by
  refine ⟨fun k x => Q01.wlog (fun N => Q01.chat N x) k, Q01.wlog_q2 Q01.chat Q01.chat_pure,
    fun x => ⟨by dsimp only; rw [Q01.wlog_two, Q01.chat_two],
      by dsimp only; rw [Q01.wlog_three, Q01.chat_three]⟩,
    fun A B hA hB K hK => ?_⟩
  obtain ⟨C6, τ₀, hτ₀, h6⟩ := Q01.chat_asymp A B hA hB K
  obtain ⟨C, τ₁, hτ₁, hlog⟩ :=
    Q01.log_asymp (fun N => Q01.chat N (X01.q A B)) (Q01.chat_zero _) K C6
  refine ⟨C, min τ₀ τ₁, lt_min hτ₀ hτ₁, fun p hp hτ => ?_⟩
  have hp0 : 0 < p := by linarith
  have hτpos : 0 < Q01.tau p := by
    unfold Q01.tau
    have : 1 / p < 1 := by rw [div_lt_one hp0]; exact hp
    linarith
  have h := hlog (Q01.tau p) _ hτpos (hτ.trans (min_le_right _ _))
    (h6 p hp (hτ.trans (min_le_left _ _)))
  have hs : ∑ k ∈ Finset.range (K + 1), Q01.wlog (fun N => Q01.chat N (X01.q A B)) k *
        Q01.tau p ^ k =
      ∑ k ∈ Finset.Icc 2 K, Q01.wlog (fun N => Q01.chat N (X01.q A B)) k * Q01.tau p ^ k := by
    rw [Finset.range_eq_Ico, Finset.sum_eq_sum_Ico_succ_bot (by omega),
      Finset.sum_eq_sum_Ico_succ_bot (by omega), ← Finset.Ico_add_one_right_eq_Icc,
      Q01.wlog_zero, Q01.wlog_one]
    simp
  rw [hs] at h
  exact h

namespace Q01

lemma zq_eval0 {k : ℕ} {f : ℝ → ℝ} (hf : IsZetaPolyQ2 k f) : IsZetaPoly k (f 0) := by
  obtain ⟨L, h1, e1⟩ := hf
  refine ⟨L.map fun e => (if e.2.1 = 0 then e.1 else 0, e.2.2), fun e he => ?_, ?_⟩
  · obtain ⟨a, ha, rfl⟩ := List.mem_map.1 he
    exact h1 a ha
  · rw [e1 0, List.map_map]
    congr 1
    apply List.map_congr_left
    intro e _
    by_cases h : e.2.1 = 0 <;> simp [h]

lemma r11 : X01.r 1 1 = Real.sqrt 2 := by unfold X01.r; norm_num
lemma q11 : X01.q 1 1 = 0 := by unfold X01.q; norm_num

end Q01

theorem Q01.proof_8 : Q01.stmt_8 := by
  obtain ⟨c, hc0, _, hcp, hc⟩ := Q01.proof_6
  refine ⟨fun n => 4 * c n 0, by simp [hc0], fun n => ?_, fun K => ?_⟩
  · have := Q01.zp_smul 4 (Q01.zq_eval0 (hcp n))
    simpa using this
  · obtain ⟨C, τ₀, hτ₀, h⟩ := hc 1 1 one_pos one_pos K
    refine ⟨4 * C, τ₀, hτ₀, fun p hp hτ => ?_⟩
    have h1 := h p hp hτ
    rw [Q01.r11, Q01.q11] at h1
    have hs : Real.sqrt 2 ≠ 0 := by positivity
    have e : X01.P 1 1 p / Real.sqrt 2 - ∑ n ∈ Finset.range (K + 1), 4 * c n 0 * Q01.tau p ^ n =
        4 * (X01.P 1 1 p / (4 * Real.sqrt 2) -
          ∑ N ∈ Finset.range (K + 1), c N 0 * Q01.tau p ^ N) := by
      rw [mul_sub, Finset.mul_sum]
      congr 1
      · field_simp
      · refine Finset.sum_congr rfl fun n _ => by ring
    rw [e, abs_mul, abs_of_pos (by norm_num : (0:ℝ) < 4), mul_assoc]
    exact mul_le_mul_of_nonneg_left h1 (by norm_num)


lemma Q01.P11 (p : ℝ) : X01.P 1 1 p = X05.P p := by
  unfold X01.P X05.P X01.integrand; simp

set_option maxHeartbeats 1000000 in
lemma Q01.alg9 {z2 z3 z4 z5 C T Pv : ℝ} (hz : z2 ^ 2 = 5 / 2 * z4)
    (h2 : |z2| ≤ 2) (h3 : |z3| ≤ 2) (h4 : |z4| ≤ 2) (h5 : |z5| ≤ 2)
    (hT0 : 0 < T) (hT1 : T ≤ 1) (hTs : T ≤ 1 / (2 * (6 + |C|)))
    (h : |4 * Pv - (4 + z2 * T ^ 2 + z3 * T ^ 3 - 33 / 16 * z4 * T ^ 4
        - (3 / 4 * z5 + 7 / 4 * z2 * z3) * T ^ 5)| ≤ C * T ^ 6) :
    |Real.log Pv - (z2 / 4 * T ^ 2 + z3 / 4 * T ^ 3 + (-z2 ^ 2 / 4 + z4 / 32) * T ^ 4
        + (-(z2 * z3) / 2 - 3 * z5 / 16) * T ^ 5)|
      ≤ (2 * (6 + |C|) ^ 3 + 7 * |C| + |C| ^ 2 + 18) * T ^ 6 := by
  set K := |C| with hK
  have hK0 : 0 ≤ K := abs_nonneg C
  set u := z2 / 4 + z3 / 4 * T with hu
  set v := -33 / 64 * z4 - (3 / 16 * z5 + 7 / 16 * z2 * z3) * T with hv
  set y := T ^ 2 * u + T ^ 4 * v with hy
  set x := Pv - 1 with hx
  set e := x - y with he
  clear_value e x y v u K
  have hT2 : T ^ 2 ≤ 1 := pow_le_one₀ hT0.le hT1
  have hT6 : 0 < T ^ 6 := by positivity
  have hT42 : T ^ 4 ≤ T ^ 2 := pow_le_pow_of_le_one hT0.le hT1 (by norm_num)
  have hT62 : T ^ 6 ≤ T ^ 2 := pow_le_pow_of_le_one hT0.le hT1 (by norm_num)
  have hT21 : T ^ 2 ≤ T := by nlinarith
  have he' : |e| ≤ K / 4 * T ^ 6 := by
    have : 4 * Pv - (4 + z2 * T ^ 2 + z3 * T ^ 3 - 33 / 16 * z4 * T ^ 4
        - (3 / 4 * z5 + 7 / 4 * z2 * z3) * T ^ 5) = 4 * e := by
      simp only [he, hx, hy, hu, hv]; ring
    rw [this, abs_mul] at h
    have : C * T ^ 6 ≤ K * T ^ 6 := by gcongr; rw [hK]; exact le_abs_self C
    norm_num at h; linarith
  have hz23 : |z2 * z3| ≤ 4 := by rw [abs_mul]; nlinarith [abs_nonneg z2, abs_nonneg z3]
  have hu' : |u| ≤ 1 := by
    have : |z3 / 4 * T| ≤ 1 / 2 := by
      rw [abs_mul, abs_of_pos hT0, abs_div]; norm_num
      nlinarith [abs_nonneg z3]
    have h2' : |z2 / 4| ≤ 1 / 2 := by rw [abs_div]; norm_num; linarith
    calc |u| ≤ |z2 / 4| + |z3 / 4 * T| := by rw [hu]; exact abs_add_le _ _
      _ ≤ 1 := by linarith
  have hv' : |v| ≤ 5 := by
    have h1 : |(3 / 16 * z5 + 7 / 16 * z2 * z3) * T| ≤ 3 := by
      rw [abs_mul, abs_of_pos hT0]
      have : |3 / 16 * z5 + 7 / 16 * z2 * z3| ≤ 3 := by
        calc _ ≤ |3 / 16 * z5| + |7 / 16 * (z2 * z3)| := by
                rw [mul_assoc]; exact abs_add_le _ _
          _ ≤ 3 := by
              have a1 : |3 / 16 * z5| ≤ 3 / 8 := by
                rw [abs_mul, abs_of_pos (by norm_num : (0:ℝ) < 3 / 16)]; linarith
              have a2 : |7 / 16 * (z2 * z3)| ≤ 7 / 4 := by
                rw [abs_mul, abs_of_pos (by norm_num : (0:ℝ) < 7 / 16)]; linarith
              linarith
      nlinarith [abs_nonneg (3 / 16 * z5 + 7 / 16 * z2 * z3)]
    have a1 : |-33 / 64 * z4| ≤ 33 / 32 := by rw [abs_mul]; norm_num; linarith
    calc |v| ≤ |-33 / 64 * z4| + |(3 / 16 * z5 + 7 / 16 * z2 * z3) * T| := by
          rw [hv]; exact abs_sub _ _
      _ ≤ 5 := by linarith
  have hy' : |y| ≤ 6 * T ^ 2 := by
    calc |y| ≤ |T ^ 2 * u| + |T ^ 4 * v| := by rw [hy]; exact abs_add_le _ _
      _ = T ^ 2 * |u| + T ^ 4 * |v| := by
          rw [abs_mul, abs_mul, abs_of_pos (by positivity : (0:ℝ) < T ^ 2),
            abs_of_pos (by positivity : (0:ℝ) < T ^ 4)]
      _ ≤ T ^ 2 * 1 + T ^ 2 * 5 := by
          gcongr
      _ = 6 * T ^ 2 := by ring
  have hx' : |x| ≤ (6 + K) * T ^ 2 := by
    have : x = y + e := by simp only [he]; ring
    rw [this]
    calc |y + e| ≤ |y| + |e| := abs_add_le _ _
      _ ≤ 6 * T ^ 2 + K * T ^ 2 := by
          have : K / 4 * T ^ 6 ≤ K * T ^ 2 := by
            have : K / 4 ≤ K := by linarith
            calc K / 4 * T ^ 6 ≤ K * T ^ 2 := by gcongr
              _ = _ := rfl
          linarith
      _ = _ := by ring
  have h5K : 0 < 6 + K := by linarith
  have hxh : |x| ≤ 1 / 2 := by
    calc |x| ≤ (6 + K) * T ^ 2 := hx'
      _ ≤ (6 + K) * T := by gcongr
      _ ≤ (6 + K) * (1 / (2 * (6 + K))) := by gcongr
      _ = 1 / 2 := by field_simp
  have hlog := Real.abs_log_sub_add_sum_range_le (x := -x) (by rw [abs_neg]; linarith) 2
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, abs_neg] at hlog
  have hPv : 1 - -x = Pv := by simp only [hx]; ring
  rw [hPv] at hlog
  have hx3 : |x| ^ (2 + 1) / (1 - |x|) ≤ 2 * ((6 + K) ^ 3 * T ^ 6) := by
    rw [div_le_iff₀ (by linarith)]
    have : |x| ^ 3 ≤ ((6 + K) * T ^ 2) ^ 3 := by gcongr
    have : ((6 + K) * T ^ 2) ^ 3 = (6 + K) ^ 3 * T ^ 6 := by ring
    have : 0 ≤ |x| ^ 3 := by positivity
    norm_num
    nlinarith
  -- the algebraic identity
  have hid : Real.log Pv - (z2 / 4 * T ^ 2 + z3 / 4 * T ^ 3 + (-z2 ^ 2 / 4 + z4 / 32) * T ^ 4
        + (-(z2 * z3) / 2 - 3 * z5 / 16) * T ^ 5) =
      (-x + x ^ 2 / 2 + Real.log Pv)
        + (e - y * e - e ^ 2 / 2 - T ^ 6 * ((z3 / 4) ^ 2 + 2 * u * v + T ^ 2 * v ^ 2) / 2) := by
    have hxe : x = y + e := by simp only [he]; ring
    rw [hxe]; simp only [hy, hu, hv]
    linear_combination ((7 / 32) * T ^ 4) * hz
  rw [hid]
  have hye : |y * e| ≤ 3 / 2 * K * T ^ 6 := by
    rw [abs_mul]
    calc |y| * |e| ≤ (6 * T ^ 2) * (K / 4 * T ^ 6) := by gcongr
      _ ≤ (6 * 1) * (K / 4 * T ^ 6) := by gcongr
      _ = _ := by ring
  have he2 : |e ^ 2 / 2| ≤ K ^ 2 * T ^ 6 := by
    rw [abs_div, abs_pow, abs_of_pos (two_pos : (0:ℝ) < 2)]
    have : |e| ^ 2 ≤ (K / 4 * T ^ 6) ^ 2 := by gcongr
    have : (K / 4 * T ^ 6) ^ 2 = K ^ 2 / 16 * T ^ 6 * T ^ 6 := by ring
    have : K ^ 2 / 16 * T ^ 6 * T ^ 6 ≤ K ^ 2 / 16 * T ^ 6 * 1 := by
      gcongr; exact pow_le_one₀ hT0.le hT1
    have : 0 ≤ K ^ 2 * T ^ 6 := by positivity
    linarith
  have hS : |T ^ 6 * ((z3 / 4) ^ 2 + 2 * u * v + T ^ 2 * v ^ 2) / 2| ≤ 18 * T ^ 6 := by
    have h1 : |(z3 / 4) ^ 2 + 2 * u * v + T ^ 2 * v ^ 2| ≤ 36 := by
      have a1 : |(z3 / 4) ^ 2| ≤ 1 / 4 := by
        rw [abs_pow]
        have : |z3 / 4| ≤ 1 / 2 := by
          rw [abs_div, abs_of_pos (by norm_num : (0:ℝ) < 4)]; linarith
        calc |z3 / 4| ^ 2 ≤ (1 / 2) ^ 2 := by gcongr
          _ = 1 / 4 := by norm_num
      have a2 : |2 * u * v| ≤ 10 := by
        rw [abs_mul, abs_mul, abs_two]
        calc 2 * |u| * |v| ≤ 2 * 1 * 5 := by gcongr
          _ = 10 := by norm_num
      have a3 : |T ^ 2 * v ^ 2| ≤ 25 := by
        rw [abs_mul, abs_pow, abs_pow, abs_of_pos hT0]
        calc T ^ 2 * |v| ^ 2 ≤ 1 * 5 ^ 2 := by gcongr
          _ = 25 := by norm_num
      calc _ ≤ |(z3 / 4) ^ 2 + 2 * u * v| + |T ^ 2 * v ^ 2| := abs_add_le _ _
        _ ≤ |(z3 / 4) ^ 2| + |2 * u * v| + |T ^ 2 * v ^ 2| := by
            gcongr; exact abs_add_le _ _
        _ ≤ 36 := by linarith
    rw [abs_div, abs_mul, abs_of_pos hT6, abs_two]
    calc T ^ 6 * |(z3 / 4) ^ 2 + 2 * u * v + T ^ 2 * v ^ 2| / 2 ≤ T ^ 6 * 36 / 2 := by gcongr
      _ = 18 * T ^ 6 := by ring
  have hA : |-x + x ^ 2 / 2 + Real.log Pv| ≤ 2 * ((6 + K) ^ 3 * T ^ 6) := by
    have := hlog.trans hx3
    convert this using 3; ring
  have hB1 := abs_sub (e - y * e - e ^ 2 / 2)
    (T ^ 6 * ((z3 / 4) ^ 2 + 2 * u * v + T ^ 2 * v ^ 2) / 2)
  have hB2 := abs_sub (e - y * e) (e ^ 2 / 2)
  have hB3 := abs_sub e (y * e)
  have hC := abs_add_le (-x + x ^ 2 / 2 + Real.log Pv)
    (e - y * e - e ^ 2 / 2 - T ^ 6 * ((z3 / 4) ^ 2 + 2 * u * v + T ^ 2 * v ^ 2) / 2)
  have hK3 : 0 ≤ K * T ^ 6 := by positivity
  nlinarith

theorem Q01.proof_9 : Q01.stmt_9 := by
  obtain ⟨C, δ, hδ, h7⟩ := X05.proof_7
  refine ⟨2 * (6 + |C|) ^ 3 + 7 * |C| + |C| ^ 2 + 18, min (min δ 1) (1 / (2 * (6 + |C|))),
    lt_min (lt_min hδ one_pos) (by positivity), fun p hp hτ => ?_⟩
  have hp0 : 0 < p := by linarith
  have hT : 1 - 1 / p = Q01.tau p := rfl
  have hTpos : 0 < Q01.tau p := by
    unfold Q01.tau
    have : 1 / p < 1 := by rw [div_lt_one hp0]; exact hp
    linarith
  have hTδ : Q01.tau p ≤ δ := hτ.trans ((min_le_left _ _).trans (min_le_left _ _))
  have hT1 : Q01.tau p ≤ 1 := hτ.trans ((min_le_left _ _).trans (min_le_right _ _))
  have hTs : Q01.tau p ≤ 1 / (2 * (6 + |C|)) := hτ.trans (min_le_right _ _)
  have h := h7 p hp0 (by rw [hT, abs_of_pos hTpos]; exact hTδ)
  rw [hT, abs_of_pos hTpos, ← X05.z_two, ← X05.z_four, ← Q01.P11] at h
  have hs2 : Real.sqrt 2 ≠ 0 := by positivity
  have hP : X01.P 1 1 p / Real.sqrt 2 = 4 * (X01.P 1 1 p / (4 * Real.sqrt 2)) := by
    field_simp
  rw [hP] at h
  have hz : X05.z 2 ^ 2 = 5 / 2 * X05.z 4 := by rw [X05.z_two, X05.z_four]; ring
  exact Q01.alg9 hz (X05.z_bd le_rfl) (X05.z_bd (by norm_num)) (X05.z_bd (by norm_num))
    (X05.z_bd (by norm_num)) hTpos hT1 hTs h

theorem Q01.proof_10 : Q01.stmt_10 := by
  obtain ⟨w, hwp, _, hw⟩ := Q01.proof_7
  refine ⟨fun k => w k 0, fun k => Q01.zq_eval0 (hwp k), fun K hK => ?_⟩
  obtain ⟨C, τ₀, hτ₀, h⟩ := hw 1 1 one_pos one_pos K hK
  refine ⟨C, τ₀, hτ₀, fun p hp hτ => ?_⟩
  have := h p hp hτ
  rwa [Q01.r11, Q01.q11] at this
