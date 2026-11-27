import Mathlib
import Zeta3.Wave10.S04_M4alphaZetaParts
import Zeta3.Wave10.S08_FieldRuleMain
/-!
# P04 -- `M_4(α)` at every angle, with the closed forms of `R_4` and `Ψ_2`

Target id: P04.
Paper label: `thm:M4alpha` (paper/ieee/main.tex, lines 1572-1583).
Supplement: paper/ieee/proofs/sec_g8_M4alpha.tex (setting 7-8; `g8:thm:M4alpha` 10-20; the data
files, line 22; the precise form `g8:eq:M4formal` 24-33; "Why a finite computation proves the
identity for every n" 37-51; M8, line 57; proof 74-76) and paper/ieee/proofs/sec_g8_M2alpha.tex
(the criterion section 37-99: `g8:lem:formal` 41-43, `g8:lem:shift` 49-60, `g8:lem:hcoef` 62-68,
`g8:lem:mu` 70-72, `g8:prop:crit` 76-91).
Data: paper/ieee/proofs/scripts/g8_terms_k4_partz2.json (`Ψ_2`) and g8_terms_k4_part1.json
(`R_4`); checked in the supplement by g8_allN.py (library g8_zsum.py), output g8_allN_k4.out.

This is the stretch continuation of `Zeta3.Wave10.S04_M4alphaZetaParts`, which proved the theorem
with `R_4 = [1]M_4` and `Ψ_2 = [ζ(2)]M_4` as SERIES parts. Here the two series are identified with
the finite combinations of harmonic polylogarithms of the data files.

Imported, not restated (all from frozen statement blocks): `M_k(α) = V05.Mk k α`,
`c = D04.c α = cos(α/2)`, the real zeta values `ζ(w) = X05.z w = (riemannZeta w).re`
(`zeta_nat_eq_tsum_of_gt_one` identifies `riemannZeta w` with `Σ_{n≥1} n^{-w}` for `w ≥ 2`; it is
real there), the formal coefficients `P_{k,n} = Q04.Pf k n`, the series parts
`[z]M_k(α) = R01.zM s k α` with `R01.AllSummable α` (absolute convergence), and from S03 the
harmonic polylogarithms `S03.H w x` (iterated integrals `∫_0^x f_a H_{w'}`), the nested sums
`S03.Z v M`, the dictionary data `S03.dd w = (ε, κ, 𝔠_w, v_w)` of `g6:lem:dict` (proved in
`S03.stmt_5`: `H_w(x) = Σ_{N≥1} 𝔠_w ε^N N^{-κ} Z_{v_w}(N-1) x^N` for `|x| < 1`), the elements of
`𝔇` as lists `S03.DTerm` of triples `(Q, e, w)` meaning `Q(N)/Π_p (N-p)^{e_p} [w]`, with
`S03.ratEval` and the evaluation `S03.ST T N = Σ R_w(N) Z_w(N)` (so `⟨T⟩_n = S03.ST T (2n)`).
`S04.stmt_7` (proved) is the theorem with the series parts; it is used, not restated.
`Zeta3.Wave10.S08_FieldRuleMain` is imported as the job asks (it uses the objects of S03); none of its
statements is used in a statement here.

Quoted verbatim (main.tex 1572-1583):

  \begin{theorem}[$M_4$ at every angle$^\dagger$]\label{thm:M4alpha}
  For $0<\alpha<\pi$, with $c=\cos\frac\alpha2$ and $\Psi$ the $\zeta(2)$ part of $M_3$
  (Theorem~\ref{thm:flat}),
  \[
  \begin{split}
  M_4(\alpha)={}&R_4(c)+\zeta(2)\Psi_2(c)-\zeta(3)\Psi(c)\\
  &+\frac{\pi^4}{2880}\Bigl[\Bigl(9c-\frac{30}c+\frac{21}{c^3}\Bigr)\operatorname{artanh}c+23-\frac{21}{c^2}\Bigr],
  \end{split}
  \]
  where $R_4$ is a rational combination of $109$ terms $c^eH_w(\pm c)$ with harmonic
  polylogarithms of weights two to five and $-3\le e\le1$, and $\Psi_2$ is a rational combination
  of $21$ such terms of weights one to three and the constant $\frac12$ (in our notes).
  \end{theorem}

Quoted verbatim (supplement sec_g8_M4alpha.tex, line 8): "A term $c^eH_w(\pm c)$ means the pair
$c^e\bigl(H_w(c)+s\,H_{\bar w}(c)\bigr)$, $s=(-1)^{e+\nu(w)}$, which is $c^eH_w(c)+(-c)^eH_w(-c)$".
(line 22) "The $109$ terms of $R_4$ are the file \script{g8_terms_k4_part1.json}, the terms of
$\Psi_2$ are \script{g8_terms_k4_partz2.json} (each line is $[\text{``}c^e\,\mathrm{He}[w]\text{''},q]$:
the term $q\,c^e(H_w+sH_{\bar w})$). ... In what follows $\Psi_2$ is the constant $\frac12$ plus
these $21$ terms".
(line 33) "where, more precisely, $R_4$, $\Psi_2$, $\Psi$, $F_1$, $F_2$ are the sums
$\sum_n[\gvs z]P_{4,n}c^{2n}$ for $\gvs z=1,\gvs2,-\gvs3,\gvs4,\gvs2^2$".

Quoted verbatim (sec_g8_M2alpha.tex, `g8:lem:shift`, 49-55): "For every index word $v$ and
integer $j\ge-2$ there is an element $S_j(v)\in\mathfrak D$ with $\gsub{S_j(v)}{n}=Z_v(2n+j)$ for
every $n\ge1$. It is defined by $S_j(\emptyset)=[\emptyset]$, $S_0(v)=[v]$ and, for $v=(a,v')$
with $\sigma=\gsgn(a)$,
  S_j(v)=[v]+\sum_{t=1}^{j}\frac{\sigma^t}{(N+t)^{|a|}}\,S_{t-1}(v')\quad(j>0),\qquad
  S_j(v)=[v]-\sum_{t=j+1}^{0}\frac{\sigma^t}{(N+t)^{|a|}}\,S_{t-1}(v')\quad(j<0),
where $N=2n$. All its denominators are powers of $N+t$ with $t\ge-1$, so the poles lie at
integers $\le1$."

Quoted verbatim (`g8:lem:hcoef`, 62-63): "Let $w$ be a nonempty regular word, with data
$(\varepsilon,\kappa,v_w,\gcf_w)$ as in Lemma~\ref{g6:lem:dict}. For every integer $\delta\ge-1$
put $h_{w,\delta}=\gcf_w\,\varepsilon^\delta\,(N+\delta)^{-\kappa}\,S_{\delta-1}(v_w)\in\mathfrak D$.
Then $[x^{2n+\delta}]H_w(x)=\gsub{h_{w,\delta}}{n}$ for every $n\ge1$."

Quoted verbatim (`g8:prop:crit`, 76-87): "Let $k\ge0$ and let $\gvs z$ be as above; write
$p_n=[\gvs z]P_{k,n}$. Let
  L(c)=\sum_i q_ic^{e_i}H_{w_i}(c)\;+\;\sum_jq'_jc^{e'_j}\bigl(H_{w'_j}(c)+s_jH_{\bar w'_j}(c)\bigr)
  \;+\;\sum_{n\ge1}\rho_nc^{2n}\;+\;\sum_lq''_lc^{e''_l}
be a finite sum of this form, plus finitely many terms $q''c^{e''}$ ($q''\in\QQ$, $e''\le1$
integers; ...), with rational $q,q'$, all exponents $e_i,e'_j\le1$, regular words, $\bar w$ the
word with every nonzero letter negated, $s_j=(-1)^{e'_j+\nu(w'_j)}$ ... For $\iota\in\{0,1\}$ let
$L_\iota\in\mathfrak D$ be the sum of the elements $q\,h_{w_i,\iota-e_i}$ and
$q'(h_{w'_j,\iota-e'_j}+s_jh_{\bar w'_j,\iota-e'_j})$ of Lemma~\ref{g8:lem:hcoef} ... Suppose
(i) $L_0=T_z$ and $L_1=0$ in $\mathfrak D$; (ii) the coefficients of $c^N$, for all $N\le1$, in
the Laurent series of $L-\sum_{n\ge1}\rho_nc^{2n}$ ... vanish ... Then
$L(c)=\sum_{n\ge1}p_nc^{2n}$ for $0<c<1$." (`T_z` is the element of `g8:lem:mu`:
"$[\gvs z]P_{k,n}=\gsub{T_z}{n}$ for every $n\ge1$".)

Quoted verbatim (sec_g8_M4alpha.tex, proof, 75): "Lemma~\ref{g8:prop:crit} (for $\Psi_2$ with its
constant term $\frac12$) with M8 gives the parts $1$ and $\gvs2$".

## Formalization choices

* The term lists are NOT printed in the paper or the supplement: they are copied exactly from the
  data files g8_terms_k4_partz2.json (21 terms and the constant `1/2`) and g8_terms_k4_part1.json
  (109 terms, no constant). A file line `["c^e He[w]", "q"]` is the triple `(q, e, w) : PTerm`
  (all lines of both files are of the paired kind `He`), standing for
  `q c^e (H_w(c) + s H_{w̄}(c))`, `s = (-1)^{e+ν(w)}` (`sgnP e w`), `w̄ = bar w` (every letter
  negated), `ν(w)` the number of nonzero letters (`nu w`). The line `["1", "1/2"]` is the constant
  term `(1/2, 0)` (a `q'' c^{e''}` term). `Psi2T`, `Psi2C`, `R4T`, `R4C = []` are these lists in
  file order; `Lfun T C x` is the function `Σ q'' x^{e''} + Σ q x^e (H_w(x) + s H_{w̄}(x))` (real
  integer powers `x ^ (e : ℤ)`; only used at `0 < x < 1`, so `x^{-3}` is not a junk value).
  `Psi2 = Lfun Psi2T Psi2C` and `R4 = Lfun R4T []`.
* `hc w N` is the Taylor coefficient `[x^N]H_w = 𝔠_w ε^N N^{-κ} Z_{v_w}(N-1)` of `g6:lem:dict`
  (`N ≥ 1`; `0` for `N = 0`, as `H_w(0) = 0`), written with S03's dictionary data `S03.dd`; by
  `S03.stmt_5` it is the Taylor coefficient of `S03.H w`. `coefL T C N` is then the coefficient
  of `x^N` (`N ∈ ℤ`) in the Laurent series of `Lfun T C`.
* The shift elements `S_j(v)` are `Sh v j : List S03.DTerm`, by the recursion of `g8:lem:shift`
  (`scl r p k τ` multiplies a term by `r/(N-p)^k`; `(N+t)^{|a|}` is `p = -t`). `h_{w,δ}` is
  `hW w δ`, and `L_ι` is `Lι T ι`. Equality in `𝔇` ("the same `R_w` for every `w`") is `Deq`: for
  every word `u` the sums of the coefficients `R` of the terms with word `u` agree at every
  rational `N > 1`. For elements whose poles lie at integers `≤ 1` (all `h_{w,δ}`, `δ ≥ -1`, by
  `stmt_2`) this is equality of the rational functions (they agree on an infinite set).
* The criterion (`stmt_3`) is stated in the form used for `k = 4`: paired terms and constant terms
  only (no unpaired terms `q c^e H_w` and no extra series `Σ ρ_n c^{2n}`, which occur only for
  `k = 2`), for every `k` and every monomial `𝔷^s`, with `T_z` any element with
  `[𝔷^s]P_{k,n} = ⟨T_z⟩_n` for all `n ≥ 1`. The words are nonempty regular words over `{0,1,-1}`
  (the data contain no empty word except through the constant terms).
* `stmt_4`/`stmt_6` (the content of M8) are stated directly on the coefficients: the coefficient
  of `c^{2n}` of the term list equals `[𝔷^s]P_{4,n}` for every `n ≥ 1`, the coefficient of
  `c^{2n+1}` vanishes, and all coefficients of `c^N`, `N ≤ 1`, vanish. They are statements about
  rational numbers for every `n`; by `stmt_2`, `stmt_3` they follow from the finitely many
  rational-function identities `Deq (Lι T 0) T_z`, `Deq (Lι T 1) []` that the supplement's script
  checks.
* `stmt_5`/`stmt_7`: as functions on `0 < c < 1` (`HasSum`, so the series converges). The parts
  are `[𝔷2]` (`Finsupp.single 2 1`) and `[1]` (`0`), as in S04.
* `Ψ = [ζ(2)]M_3(α) = R01.zM (single 2 1) 3 α` stays the series (as in S04; no closed form is
  asked for here).

## Statements

* stmt_1: `g8:lem:shift` (poles at integers `≤ 1`, and `⟨S_j(v)⟩_n = Z_v(2n+j)`, `j ≥ -2`).
* stmt_2: `g8:lem:hcoef` (`[x^{2n+δ}]H_w = ⟨h_{w,δ}⟩_n`, `δ ≥ -1`, poles at integers `≤ 1`).
* stmt_3: `g8:prop:crit` in the form used.
* stmt_4: M8 for `Ψ_2`: every Taylor coefficient of the 21 terms plus `1/2` is right.
* stmt_5: `Σ_{n≥1} [𝔷2]P_{4,n} c^{2n} = Ψ_2(c)` for `0 < c < 1`.
* stmt_6: M8 for `R_4`: every Taylor coefficient of the 109 terms is right.
* stmt_7: `Σ_{n≥1} [1]P_{4,n} c^{2n} = R_4(c)` for `0 < c < 1`.
* stmt_8: `thm:M4alpha` with the closed forms of `R_4` and `Ψ_2`.

## Numerical check before freezing

A Python script (scratch, not kept) mirroring the definitions above
checked: the Lean data lists, parsed back from this file, equal the two JSON files line by line
(21 + constant `1/2`, and 109); `stmt_1` for all index words over `{-3,-2,-1,1,2}` of length
`≤ 3`, `-2 ≤ j ≤ 4`, `n ≤ 5` (exact), and the poles; `stmt_2` and `hc` against the coefficient
recursion of the iterated integrals for all regular words of length `≤ 4` (exact). `stmt_4`,
`stmt_6`: `coefL` at `c^{2n}` equals `[𝔷2]P_{4,n}`, resp. `[1]P_{4,n}`, from the supplement's
independent routine g6_common.py (exact rationals), for `n = 1, …, 80`; the coefficients of
`c^{2n+1}` (`n ≤ 80`) and of `c^N`, `-10 ≤ N ≤ 1`, are `0`; `P_{4,n}` has only the monomials
`1, 𝔷2, 𝔷3, 𝔷4, 𝔷2²`. With mpmath (45 digits; `H_w` summed from the dictionary series, 220
terms), at `α = 2π/3` (`c = 1/2`) and `α = 2.5`: `Ψ_2(c) - Σ_{n≤80}[𝔷2]P_{4,n}c^{2n}` is
`2.5·10⁻⁴⁶`, `4.0·10⁻⁴⁶`; `R_4(c) - Σ_{n≤80}[1]P_{4,n}c^{2n}` is `-3.1·10⁻⁴⁶`, `-1.4·10⁻⁴⁶`
(`Ψ_2(1/2) = 0.067783306883119476707872587590`, `R_4(1/2) = -0.140447335287165935793247545488`);
the right side of `stmt_8` (with `Ψ` the series `Σ_{n≤80}[𝔷2]P_{3,n}c^{2n}`) differs from
`M_4(α)` computed from the definition of `D04.Mint` (Cauchy integral in `u`, radius `0.2`, 72
nodes) by `4.7·10⁻⁴²` and `7.1·10⁻⁴²` (`M_4(2π/3) = 0.0021017551526063123975140574239779`).

## NOT COVERED

* The general form of the criterion with unpaired terms and an extra series `Σ ρ_n c^{2n}` (only
  used for `k = 2`, Theorem `g8:thm:M2alpha`, done by hand in Z07).
* `g8:lem:formal` and `g8:lem:mu` as separate statements (`g8:lem:mu` is S08's open `stmt_4`; in
  `stmt_3` the element `T_z` is a hypothesis). The explicit `T_z` of the supplement's script
  (29 words for `𝔷 = 1`, 5 for `𝔷2`) is not written down.
* The counts in the theorem's wording (weights two to five, `1, 8, 27, 73` terms per weight;
  weights one to three for `Ψ_2`; `-3 ≤ e ≤ 1`) are properties of the data lists, visible in them,
  not separate statements; uniqueness or minimality of the term lists (gap box (ii)) is not
  claimed by the paper.
* A closed form of `Ψ` (the theorem uses `Ψ` from Theorem `thm:flat`; it stays the series).
* The depth-two remarks after the theorem and the numerical checks of the supplement.
-/
-- BEGIN STATEMENTS
namespace P04

open Real

/-- A paired term `(q, e, w)` of the data files: `q c^e (H_w(c) + s H_{w̄}(c))`. -/
abbrev PTerm := ℚ × ℤ × List ℤ

/-- `ν(w)`: the number of nonzero letters of `w`. -/
def nu (w : List ℤ) : ℕ := (w.filter fun a => a ≠ 0).length

/-- `w̄`: every letter negated. -/
def bar (w : List ℤ) : List ℤ := w.map fun a => -a

/-- The sign `s = (-1)^{e + ν(w)}` of the pair `c^e H_w(±c)`. -/
def sgnP (e : ℤ) (w : List ℤ) : ℚ := (-1 : ℚ) ^ (e + (nu w : ℤ))

/-- The function `Σ_l q''_l x^{e''_l} + Σ_j q_j x^{e_j} (H_{w_j}(x) + s_j H_{w̄_j}(x))` of a list
`T` of paired terms and a list `C` of constant terms `(q'', e'')`. -/
noncomputable def Lfun (T : List PTerm) (C : List (ℚ × ℤ)) (x : ℝ) : ℝ :=
  (C.map fun t => (t.1 : ℝ) * x ^ t.2).sum +
  (T.map fun t => (t.1 : ℝ) * x ^ t.2.1 *
    (S03.H t.2.2 x + (sgnP t.2.1 t.2.2 : ℝ) * S03.H (bar t.2.2) x)).sum

/-- `[x^N]H_w = 𝔠_w ε^N N^{-κ} Z_{v_w}(N-1)` for `N ≥ 1` (Lemma `g6:lem:dict`, data `S03.dd w`),
and `0` for `N = 0`. -/
def hc (w : List ℤ) (N : ℕ) : ℚ :=
  if N = 0 then 0 else
    ((S03.dd w).2.2.1 : ℚ) * ((S03.dd w).1 : ℚ) ^ N / (N : ℚ) ^ (S03.dd w).2.1 *
      S03.Z (S03.dd w).2.2.2 (N - 1)

/-- The coefficient of `x^N` (`N ∈ ℤ`) in the Laurent series of `Lfun T C`. -/
def coefL (T : List PTerm) (C : List (ℚ × ℤ)) (N : ℤ) : ℚ :=
  (C.map fun t => if t.2 = N then t.1 else 0).sum +
  (T.map fun t => if 1 ≤ N - t.2.1 then
      t.1 * (hc t.2.2 (N - t.2.1).toNat + sgnP t.2.1 t.2.2 * hc (bar t.2.2) (N - t.2.1).toNat)
    else 0).sum

/-- Multiply a term `R [w]` of `𝔇` by `r / (N - p)^k`. -/
noncomputable def scl (r : ℚ) (p : ℤ) (k : ℕ) (τ : S03.DTerm) : S03.DTerm :=
  (Polynomial.C r * τ.1, τ.2.1 + Finsupp.single p k, τ.2.2)

/-- The shift elements `S_j(v) ∈ 𝔇` of `g8:lem:shift`: `S_j(∅) = [∅]`, and for `v = (a, v')`,
`S_j(v) = [v] + Σ_{t=1}^{j} σ^t (N+t)^{-|a|} S_{t-1}(v')` (`j > 0`),
`S_j(v) = [v] - Σ_{t=j+1}^{0} σ^t (N+t)^{-|a|} S_{t-1}(v')` (`j ≤ 0`), `σ = sgn a`. -/
noncomputable def Sh : List ℤ → ℤ → List S03.DTerm
  | [], _ => [(1, 0, [])]
  | a :: v, j =>
      (1, 0, a :: v) ::
      (if 0 < j then
        (List.range j.toNat).flatMap fun i =>
          (Sh v (i : ℤ)).map (scl ((Int.sign a : ℚ) ^ ((i : ℤ) + 1)) (-((i : ℤ) + 1)) a.natAbs)
      else
        (List.range (-j).toNat).flatMap fun i =>
          (Sh v (-(i : ℤ) - 1)).map (scl (-(Int.sign a : ℚ) ^ (-(i : ℤ))) (i : ℤ) a.natAbs))

/-- `h_{w,δ} = 𝔠_w ε^δ (N+δ)^{-κ} S_{δ-1}(v_w)` (`g8:lem:hcoef`). -/
noncomputable def hW (w : List ℤ) (δ : ℤ) : List S03.DTerm :=
  (Sh (S03.dd w).2.2.2 (δ - 1)).map
    (scl (((S03.dd w).2.2.1 : ℚ) * ((S03.dd w).1 : ℚ) ^ δ) (-δ) (S03.dd w).2.1)

/-- `L_ι = Σ_j q_j (h_{w_j, ι-e_j} + s_j h_{w̄_j, ι-e_j})` (`g8:prop:crit`). -/
noncomputable def Lι (T : List PTerm) (ι : ℤ) : List S03.DTerm :=
  T.flatMap fun t =>
    (hW t.2.2 (ι - t.2.1)).map (scl t.1 0 0) ++
    (hW (bar t.2.2) (ι - t.2.1)).map (scl (t.1 * sgnP t.2.1 t.2.2) 0 0)

/-- The coefficient `R_u(N)` of the word `u` in `T ∈ 𝔇` (the sum over the terms with word `u`). -/
noncomputable def Dval (T : List S03.DTerm) (u : List ℤ) (N : ℚ) : ℚ :=
  ((T.filter fun τ => τ.2.2 = u).map fun τ => S03.ratEval τ.1 τ.2.1 N).sum

/-- Equality in `𝔇`: for every word the coefficients agree as functions of `N`, at every rational
`N > 1` (an infinite set, so as rational functions). -/
def Deq (T T' : List S03.DTerm) : Prop :=
  ∀ u : List ℤ, ∀ N : ℚ, 1 < N → Dval T u N = Dval T' u N

/-- The 21 paired terms of `Ψ_2`, copied from g8_terms_k4_partz2.json (file order). -/
def Psi2T : List PTerm := [
    ((-1 : ℚ), (-2 : ℤ), ([1] : List ℤ)), ((-1/4 : ℚ), (-1 : ℤ), ([1] : List ℤ)),
    ((5/4 : ℚ), (0 : ℤ), ([1] : List ℤ)), ((5/4 : ℚ), (-2 : ℤ), ([0, 1] : List ℤ)),
    ((1/4 : ℚ), (-1 : ℤ), ([0, 1] : List ℤ)), ((-1 : ℚ), (0 : ℤ), ([0, 1] : List ℤ)),
    ((1/2 : ℚ), (-2 : ℤ), ([1, -1] : List ℤ)), ((1/4 : ℚ), (-1 : ℤ), ([1, -1] : List ℤ)),
    ((-3/4 : ℚ), (0 : ℤ), ([1, -1] : List ℤ)), ((-1/4 : ℚ), (-1 : ℤ), ([0, 0, 1] : List ℤ)),
    ((1/4 : ℚ), (1 : ℤ), ([0, 0, 1] : List ℤ)), ((-1/4 : ℚ), (-1 : ℤ), ([0, 1, -1] : List ℤ)),
    ((1/4 : ℚ), (1 : ℤ), ([0, 1, -1] : List ℤ)), ((-3/8 : ℚ), (-3 : ℤ), ([1, 0, 1] : List ℤ)),
    ((1/4 : ℚ), (-1 : ℤ), ([1, 0, 1] : List ℤ)), ((1/8 : ℚ), (1 : ℤ), ([1, 0, 1] : List ℤ)),
    ((3/8 : ℚ), (-3 : ℤ), ([1, 0, -1] : List ℤ)), ((-1/2 : ℚ), (-1 : ℤ), ([1, 0, -1] : List ℤ)),
    ((1/8 : ℚ), (1 : ℤ), ([1, 0, -1] : List ℤ)), ((1/4 : ℚ), (-1 : ℤ), ([1, -1, 1] : List ℤ)),
    ((-1/4 : ℚ), (1 : ℤ), ([1, -1, 1] : List ℤ))]

/-- The constant term `1/2` of `Ψ_2` (the line `["1", "1/2"]` of g8_terms_k4_partz2.json). -/
def Psi2C : List (ℚ × ℤ) := [(1/2, 0)]

/-- The 109 paired terms of `R_4`, copied from g8_terms_k4_part1.json (file order). -/
def R4T : List PTerm := [
    ((-1 : ℚ), (0 : ℤ), ([0, 1] : List ℤ)), ((1 : ℚ), (0 : ℤ), ([0, 0, 1] : List ℤ)),
    ((1 : ℚ), (0 : ℤ), ([0, 1, -1] : List ℤ)), ((2 : ℚ), (-2 : ℤ), ([1, 0, 1] : List ℤ)),
    ((1/2 : ℚ), (-1 : ℤ), ([1, 0, 1] : List ℤ)), ((-5/2 : ℚ), (0 : ℤ), ([1, 0, 1] : List ℤ)),
    ((-2 : ℚ), (-2 : ℤ), ([1, 0, -1] : List ℤ)), ((-1/2 : ℚ), (-1 : ℤ), ([1, 0, -1] : List ℤ)),
    ((5/2 : ℚ), (0 : ℤ), ([1, 0, -1] : List ℤ)), ((-1 : ℚ), (0 : ℤ), ([0, 0, 0, 1] : List ℤ)),
    ((-1 : ℚ), (0 : ℤ), ([0, 0, 1, -1] : List ℤ)), ((-5/2 : ℚ), (-2 : ℤ), ([0, 1, 0, 1] : List ℤ)),
    ((-1/2 : ℚ), (-1 : ℤ), ([0, 1, 0, 1] : List ℤ)), ((2 : ℚ), (0 : ℤ), ([0, 1, 0, 1] : List ℤ)),
    ((5/2 : ℚ), (-2 : ℤ), ([0, 1, 0, -1] : List ℤ)), ((1/2 : ℚ), (-1 : ℤ), ([0, 1, 0, -1] : List ℤ)),
    ((-3 : ℚ), (0 : ℤ), ([0, 1, 0, -1] : List ℤ)), ((1 : ℚ), (0 : ℤ), ([0, 1, -1, 1] : List ℤ)),
    ((-1 : ℚ), (-2 : ℤ), ([1, 0, 0, 1] : List ℤ)), ((-1/2 : ℚ), (-1 : ℤ), ([1, 0, 0, 1] : List ℤ)),
    ((3/2 : ℚ), (0 : ℤ), ([1, 0, 0, 1] : List ℤ)), ((1 : ℚ), (-2 : ℤ), ([1, 0, 0, -1] : List ℤ)),
    ((1/2 : ℚ), (-1 : ℤ), ([1, 0, 0, -1] : List ℤ)), ((-3/2 : ℚ), (0 : ℤ), ([1, 0, 0, -1] : List ℤ)),
    ((-1 : ℚ), (-2 : ℤ), ([1, 0, 1, -1] : List ℤ)), ((-1/2 : ℚ), (-1 : ℤ), ([1, 0, 1, -1] : List ℤ)),
    ((3/2 : ℚ), (0 : ℤ), ([1, 0, 1, -1] : List ℤ)), ((-1 : ℚ), (-2 : ℤ), ([1, 0, -1, 1] : List ℤ)),
    ((-1/2 : ℚ), (-1 : ℤ), ([1, 0, -1, 1] : List ℤ)), ((3/2 : ℚ), (0 : ℤ), ([1, 0, -1, 1] : List ℤ)),
    ((-1 : ℚ), (-2 : ℤ), ([1, -1, 0, 1] : List ℤ)), ((-1/2 : ℚ), (-1 : ℤ), ([1, -1, 0, 1] : List ℤ)),
    ((3/2 : ℚ), (0 : ℤ), ([1, -1, 0, 1] : List ℤ)), ((1 : ℚ), (-2 : ℤ), ([1, -1, 0, -1] : List ℤ)),
    ((1/2 : ℚ), (-1 : ℤ), ([1, -1, 0, -1] : List ℤ)), ((-3/2 : ℚ), (0 : ℤ), ([1, -1, 0, -1] : List ℤ)),
    ((1 : ℚ), (0 : ℤ), ([0, 0, 0, 0, 1] : List ℤ)), ((-1 : ℚ), (1 : ℤ), ([0, 0, 0, 0, 1] : List ℤ)),
    ((1 : ℚ), (0 : ℤ), ([0, 0, 0, 1, -1] : List ℤ)), ((-1 : ℚ), (1 : ℤ), ([0, 0, 0, 1, -1] : List ℤ)),
    ((1/2 : ℚ), (-1 : ℤ), ([0, 0, 1, 0, 1] : List ℤ)), ((-1/2 : ℚ), (1 : ℤ), ([0, 0, 1, 0, 1] : List ℤ)),
    ((-1/2 : ℚ), (-1 : ℤ), ([0, 0, 1, 0, -1] : List ℤ)), ((1 : ℚ), (0 : ℤ), ([0, 0, 1, 0, -1] : List ℤ)),
    ((-1/2 : ℚ), (1 : ℤ), ([0, 0, 1, 0, -1] : List ℤ)), ((-1 : ℚ), (0 : ℤ), ([0, 0, 1, -1, 1] : List ℤ)),
    ((1 : ℚ), (1 : ℤ), ([0, 0, 1, -1, 1] : List ℤ)), ((1/2 : ℚ), (-1 : ℤ), ([0, 1, 0, 0, 1] : List ℤ)),
    ((-1/2 : ℚ), (1 : ℤ), ([0, 1, 0, 0, 1] : List ℤ)), ((-1/2 : ℚ), (-1 : ℤ), ([0, 1, 0, 0, -1] : List ℤ)),
    ((1 : ℚ), (0 : ℤ), ([0, 1, 0, 0, -1] : List ℤ)), ((-1/2 : ℚ), (1 : ℤ), ([0, 1, 0, 0, -1] : List ℤ)),
    ((1/2 : ℚ), (-1 : ℤ), ([0, 1, 0, 1, -1] : List ℤ)), ((-1/2 : ℚ), (1 : ℤ), ([0, 1, 0, 1, -1] : List ℤ)),
    ((1/2 : ℚ), (-1 : ℤ), ([0, 1, 0, -1, 1] : List ℤ)), ((-1 : ℚ), (0 : ℤ), ([0, 1, 0, -1, 1] : List ℤ)),
    ((1/2 : ℚ), (1 : ℤ), ([0, 1, 0, -1, 1] : List ℤ)), ((1/2 : ℚ), (-1 : ℤ), ([0, 1, -1, 0, 1] : List ℤ)),
    ((-1 : ℚ), (0 : ℤ), ([0, 1, -1, 0, 1] : List ℤ)), ((1/2 : ℚ), (1 : ℤ), ([0, 1, -1, 0, 1] : List ℤ)),
    ((-1/2 : ℚ), (-1 : ℤ), ([0, 1, -1, 0, -1] : List ℤ)), ((1/2 : ℚ), (1 : ℤ), ([0, 1, -1, 0, -1] : List ℤ)),
    ((-1 : ℚ), (0 : ℤ), ([0, 1, -1, 1, -1] : List ℤ)), ((1 : ℚ), (1 : ℤ), ([0, 1, -1, 1, -1] : List ℤ)),
    ((1/2 : ℚ), (-1 : ℤ), ([1, 0, 0, 0, 1] : List ℤ)), ((-1/2 : ℚ), (1 : ℤ), ([1, 0, 0, 0, 1] : List ℤ)),
    ((-1/2 : ℚ), (-1 : ℤ), ([1, 0, 0, 0, -1] : List ℤ)), ((1 : ℚ), (0 : ℤ), ([1, 0, 0, 0, -1] : List ℤ)),
    ((-1/2 : ℚ), (1 : ℤ), ([1, 0, 0, 0, -1] : List ℤ)), ((1/2 : ℚ), (-1 : ℤ), ([1, 0, 0, 1, -1] : List ℤ)),
    ((-1/2 : ℚ), (1 : ℤ), ([1, 0, 0, 1, -1] : List ℤ)), ((1/2 : ℚ), (-1 : ℤ), ([1, 0, 0, -1, 1] : List ℤ)),
    ((-1 : ℚ), (0 : ℤ), ([1, 0, 0, -1, 1] : List ℤ)), ((1/2 : ℚ), (1 : ℤ), ([1, 0, 0, -1, 1] : List ℤ)),
    ((3/4 : ℚ), (-3 : ℤ), ([1, 0, 1, 0, 1] : List ℤ)), ((-1/2 : ℚ), (-1 : ℤ), ([1, 0, 1, 0, 1] : List ℤ)),
    ((-1/4 : ℚ), (1 : ℤ), ([1, 0, 1, 0, 1] : List ℤ)), ((-3/4 : ℚ), (-3 : ℤ), ([1, 0, 1, 0, -1] : List ℤ)),
    ((1 : ℚ), (-1 : ℤ), ([1, 0, 1, 0, -1] : List ℤ)), ((-1/4 : ℚ), (1 : ℤ), ([1, 0, 1, 0, -1] : List ℤ)),
    ((-1/2 : ℚ), (-1 : ℤ), ([1, 0, 1, -1, 1] : List ℤ)), ((1/2 : ℚ), (1 : ℤ), ([1, 0, 1, -1, 1] : List ℤ)),
    ((-3/4 : ℚ), (-3 : ℤ), ([1, 0, -1, 0, 1] : List ℤ)), ((3/2 : ℚ), (-1 : ℤ), ([1, 0, -1, 0, 1] : List ℤ)),
    ((-1 : ℚ), (0 : ℤ), ([1, 0, -1, 0, 1] : List ℤ)), ((1/4 : ℚ), (1 : ℤ), ([1, 0, -1, 0, 1] : List ℤ)),
    ((3/4 : ℚ), (-3 : ℤ), ([1, 0, -1, 0, -1] : List ℤ)), ((-1 : ℚ), (-1 : ℤ), ([1, 0, -1, 0, -1] : List ℤ)),
    ((1/4 : ℚ), (1 : ℤ), ([1, 0, -1, 0, -1] : List ℤ)), ((1/2 : ℚ), (-1 : ℤ), ([1, 0, -1, 1, -1] : List ℤ)),
    ((-1 : ℚ), (0 : ℤ), ([1, 0, -1, 1, -1] : List ℤ)), ((1/2 : ℚ), (1 : ℤ), ([1, 0, -1, 1, -1] : List ℤ)),
    ((1/2 : ℚ), (-1 : ℤ), ([1, -1, 0, 0, 1] : List ℤ)), ((-1 : ℚ), (0 : ℤ), ([1, -1, 0, 0, 1] : List ℤ)),
    ((1/2 : ℚ), (1 : ℤ), ([1, -1, 0, 0, 1] : List ℤ)), ((-1/2 : ℚ), (-1 : ℤ), ([1, -1, 0, 0, -1] : List ℤ)),
    ((1/2 : ℚ), (1 : ℤ), ([1, -1, 0, 0, -1] : List ℤ)), ((1/2 : ℚ), (-1 : ℤ), ([1, -1, 0, 1, -1] : List ℤ)),
    ((-1 : ℚ), (0 : ℤ), ([1, -1, 0, 1, -1] : List ℤ)), ((1/2 : ℚ), (1 : ℤ), ([1, -1, 0, 1, -1] : List ℤ)),
    ((1/2 : ℚ), (-1 : ℤ), ([1, -1, 0, -1, 1] : List ℤ)), ((-1/2 : ℚ), (1 : ℤ), ([1, -1, 0, -1, 1] : List ℤ)),
    ((-1/2 : ℚ), (-1 : ℤ), ([1, -1, 1, 0, 1] : List ℤ)), ((1/2 : ℚ), (1 : ℤ), ([1, -1, 1, 0, 1] : List ℤ)),
    ((1/2 : ℚ), (-1 : ℤ), ([1, -1, 1, 0, -1] : List ℤ)), ((-1 : ℚ), (0 : ℤ), ([1, -1, 1, 0, -1] : List ℤ)),
    ((1/2 : ℚ), (1 : ℤ), ([1, -1, 1, 0, -1] : List ℤ)), ((1 : ℚ), (0 : ℤ), ([1, -1, 1, -1, 1] : List ℤ)),
    ((-1 : ℚ), (1 : ℤ), ([1, -1, 1, -1, 1] : List ℤ))]

/-- `Ψ_2(c)`: the constant `1/2` plus the 21 terms. -/
noncomputable def Psi2 (x : ℝ) : ℝ := Lfun Psi2T Psi2C x

/-- `R_4(c)`: the 109 terms. -/
noncomputable def R4 (x : ℝ) : ℝ := Lfun R4T [] x

/-- `g8:lem:shift`: for an index word `v` and `j ≥ -2`, the poles of `S_j(v)` are at integers
`p ≤ 1`, and `⟨S_j(v)⟩_n = Z_v(2n+j)` for every `n ≥ 1`. -/
def stmt_1 : Prop :=
  ∀ v : List ℤ, S03.IsIndex v → ∀ j : ℤ, -2 ≤ j →
    (∀ τ ∈ Sh v j, ∀ p ∈ τ.2.1.support, p ≤ 1) ∧
    ∀ n : ℕ, 1 ≤ n → S03.ST (Sh v j) (2 * n) = S03.Z v (2 * (n : ℤ) + j).toNat

/-- `g8:lem:hcoef`: for a nonempty regular word `w` over `{0,1,-1}` and `δ ≥ -1`, the poles of
`h_{w,δ}` are at integers `p ≤ 1`, and `[x^{2n+δ}]H_w = ⟨h_{w,δ}⟩_n` for every `n ≥ 1`. -/
def stmt_2 : Prop :=
  ∀ w : List ℤ, w ≠ [] → S03.IsLetters w → S03.Regular w → ∀ δ : ℤ, -1 ≤ δ →
    (∀ τ ∈ hW w δ, ∀ p ∈ τ.2.1.support, p ≤ 1) ∧
    ∀ n : ℕ, 1 ≤ n → hc w (2 * (n : ℤ) + δ).toNat = S03.ST (hW w δ) (2 * n)

/-- `g8:prop:crit` (paired and constant terms): let `[𝔷^s]P_{k,n} = ⟨T_z⟩_n` for all `n ≥ 1`, let
`T` be a finite list of paired terms `q c^e (H_w + s H_{w̄})`, `e ≤ 1`, `w` nonempty regular over
`{0,1,-1}`, and `C` a list of terms `q'' c^{e''}`, `e'' ≤ 1`. If (i) `L_0 = T_z` and `L_1 = 0` in
`𝔇` and (ii) the coefficients of `c^N`, `N ≤ 1`, vanish, then
`Lfun T C c = Σ_{n≥1} [𝔷^s]P_{k,n} c^{2n}` for `0 < c < 1`. -/
def stmt_3 : Prop :=
  ∀ (k : ℕ) (s : ℕ →₀ ℕ) (Tz : List S03.DTerm) (T : List PTerm) (C : List (ℚ × ℤ)),
    (∀ n : ℕ, 1 ≤ n → (Q04.Pf k n).coeff s = S03.ST Tz (2 * n)) →
    (∀ t ∈ T, t.2.2 ≠ [] ∧ S03.IsLetters t.2.2 ∧ S03.Regular t.2.2 ∧ t.2.1 ≤ 1) →
    (∀ t ∈ C, t.2 ≤ 1) →
    Deq (Lι T 0) Tz → Deq (Lι T 1) [] →
    (∀ N : ℤ, N ≤ 1 → coefL T C N = 0) →
    ∀ c : ℝ, 0 < c → c < 1 →
      HasSum (fun n : ℕ => (((Q04.Pf k (n + 1)).coeff s : ℚ) : ℝ) * c ^ (2 * (n + 1)))
        (Lfun T C c)

/-- M8 for `Ψ_2` (the part `𝔷2`): the coefficient of `c^{2n}` in `Ψ_2` (the 21 terms and `1/2`)
is `[𝔷2]P_{4,n}` for every `n ≥ 1`, the coefficient of `c^{2n+1}` is `0` for every `n ≥ 1`, and
the coefficients of `c^N`, `N ≤ 1`, are `0`. -/
def stmt_4 : Prop :=
  (∀ n : ℕ, 1 ≤ n → coefL Psi2T Psi2C (2 * (n : ℤ)) = (Q04.Pf 4 n).coeff (Finsupp.single 2 1)) ∧
  (∀ n : ℕ, 1 ≤ n → coefL Psi2T Psi2C (2 * (n : ℤ) + 1) = 0) ∧
  (∀ N : ℤ, N ≤ 1 → coefL Psi2T Psi2C N = 0)

/-- The `ζ(2)` part: `Σ_{n≥1} [𝔷2]P_{4,n} c^{2n} = Ψ_2(c)` for `0 < c < 1`. -/
def stmt_5 : Prop :=
  ∀ c : ℝ, 0 < c → c < 1 →
    HasSum (fun n : ℕ => (((Q04.Pf 4 (n + 1)).coeff (Finsupp.single 2 1) : ℚ) : ℝ)
      * c ^ (2 * (n + 1))) (Psi2 c)

/-- M8 for `R_4` (the part `1`): the coefficient of `c^{2n}` in the 109 terms is `[1]P_{4,n}` for
every `n ≥ 1`, the coefficient of `c^{2n+1}` is `0` for every `n ≥ 1`, and the coefficients of
`c^N`, `N ≤ 1`, are `0`. -/
def stmt_6 : Prop :=
  (∀ n : ℕ, 1 ≤ n → coefL R4T [] (2 * (n : ℤ)) = (Q04.Pf 4 n).coeff 0) ∧
  (∀ n : ℕ, 1 ≤ n → coefL R4T [] (2 * (n : ℤ) + 1) = 0) ∧
  (∀ N : ℤ, N ≤ 1 → coefL R4T [] N = 0)

/-- The rational part: `Σ_{n≥1} [1]P_{4,n} c^{2n} = R_4(c)` for `0 < c < 1`. -/
def stmt_7 : Prop :=
  ∀ c : ℝ, 0 < c → c < 1 →
    HasSum (fun n : ℕ => (((Q04.Pf 4 (n + 1)).coeff 0 : ℚ) : ℝ) * c ^ (2 * (n + 1))) (R4 c)

/-- `thm:M4alpha`: for `0 < α < π`, `c = cos(α/2)`, with `R_4` the 109 terms, `Ψ_2` the 21 terms
and `1/2`, and `Ψ = [ζ(2)]M_3(α)` (an absolutely convergent series),
`M_4(α) = R_4(c) + ζ(2)Ψ_2(c) - ζ(3)Ψ(c) + (π⁴/2880)[(9c - 30/c + 21/c³) artanh c + 23 - 21/c²]`. -/
def stmt_8 : Prop :=
  ∀ α : ℝ, 0 < α → α < π →
    R01.AllSummable α ∧
    V05.Mk 4 α =
      R4 (D04.c α) + X05.z 2 * Psi2 (D04.c α)
        - X05.z 3 * R01.zM (Finsupp.single 2 1) 3 α
        + π ^ 4 / 2880 *
          ((9 * D04.c α - 30 / D04.c α + 21 / D04.c α ^ 3) * Real.artanh (D04.c α)
            + 23 - 21 / D04.c α ^ 2)

end P04
-- END STATEMENTS

namespace P04

lemma ratEval_scl (r : ℚ) (p : ℤ) (k : ℕ) (τ : S03.DTerm) (N : ℚ) :
    S03.ratEval (scl r p k τ).1 (scl r p k τ).2.1 N = r / (N - p) ^ k * S03.ratEval τ.1 τ.2.1 N := by
  simp only [scl, S03.ratEval, Polynomial.eval_mul, Polynomial.eval_C]
  rw [Finsupp.prod_add_index' (by intro _; simp) (by intro _ _ _; exact pow_add _ _ _),
    Finsupp.prod_single_index (by simp)]
  rw [div_mul_div_comm, mul_comm (r : ℚ), mul_comm ((N - p) ^ k)]

lemma ST_nil (N : ℕ) : S03.ST [] N = 0 := by simp [S03.ST]

lemma ST_cons (τ : S03.DTerm) (L : List S03.DTerm) (N : ℕ) :
    S03.ST (τ :: L) N = S03.ratEval τ.1 τ.2.1 (N : ℚ) * S03.Z τ.2.2 N + S03.ST L N := by
  simp [S03.ST]

lemma ST_append (L L' : List S03.DTerm) (N : ℕ) :
    S03.ST (L ++ L') N = S03.ST L N + S03.ST L' N := by
  simp [S03.ST]

lemma ST_map_scl (r : ℚ) (p : ℤ) (k : ℕ) (L : List S03.DTerm) (N : ℕ) :
    S03.ST (L.map (scl r p k)) N = r / ((N : ℚ) - p) ^ k * S03.ST L N := by
  induction L with
  | nil => simp [S03.ST]
  | cons τ L ih =>
    rw [List.map_cons, ST_cons, ST_cons, ih, ratEval_scl]
    simp only [scl]; ring

lemma ST_flatMap {α : Type*} (l : List α) (f : α → List S03.DTerm) (N : ℕ) :
    S03.ST (l.flatMap f) N = (l.map fun x => S03.ST (f x) N).sum := by
  induction l with
  | nil => simp [S03.ST]
  | cons x l ih => rw [List.flatMap_cons, ST_append, ih, List.map_cons, List.sum_cons]

lemma list_range_sum {M : Type*} [AddCommMonoid M] (f : ℕ → M) (m : ℕ) :
    ((List.range m).map f).sum = ∑ i ∈ Finset.range m, f i := by
  induction m with
  | zero => simp
  | succ m ih => rw [List.range_succ, List.map_append, List.sum_append, ih, Finset.sum_range_succ]; simp

lemma Z_succ (a : ℤ) (v : List ℤ) (M : ℕ) :
    S03.Z (a :: v) (M + 1) = S03.Z (a :: v) M +
      (Int.sign a : ℚ) ^ (M + 1) / ((M + 1 : ℕ) : ℚ) ^ a.natAbs * S03.Z v M := by
  simp only [S03.Z]
  rw [Finset.sum_Icc_succ_top (by omega)]
  simp

lemma Z_add (a : ℤ) (v : List ℤ) (M m : ℕ) :
    S03.Z (a :: v) (M + m) = S03.Z (a :: v) M + ∑ i ∈ Finset.range m,
      (Int.sign a : ℚ) ^ (M + i + 1) / ((M + i + 1 : ℕ) : ℚ) ^ a.natAbs * S03.Z v (M + i) := by
  induction m with
  | zero => simp
  | succ m ih => rw [← add_assoc, Z_succ, ih, Finset.sum_range_succ]; ring

lemma sgn_sq (a : ℤ) (ha : a ≠ 0) : (Int.sign a : ℚ) ^ 2 = 1 := by
  rcases S03.sign_ne a ha with h | h <;> rw [h] <;> norm_num

lemma sgn_even (a : ℤ) (ha : a ≠ 0) (n i : ℕ) :
    (Int.sign a : ℚ) ^ (2 * n + i) = (Int.sign a : ℚ) ^ i := by
  rw [pow_add, pow_mul, sgn_sq a ha, one_pow, one_mul]

lemma ratEval_one (N : ℚ) : S03.ratEval 1 0 N = 1 := by simp [S03.ratEval]

lemma Sh_poles : ∀ (v : List ℤ) (j : ℤ), -2 ≤ j → ∀ τ ∈ Sh v j, ∀ p ∈ τ.2.1.support, p ≤ 1
  | [], j, _ => by intro τ hτ p hp; simp [Sh] at hτ; subst hτ; simp at hp
  | a :: v, j, hj => by
    intro τ hτ p hp
    rw [Sh] at hτ
    simp only [List.pure_def, List.bind_eq_flatMap] at hτ
    rcases List.mem_cons.mp hτ with h | h
    · subst h; simp at hp
    · split_ifs at h with hj0
      · obtain ⟨i, hi, h⟩ := List.mem_flatMap.mp h
        obtain ⟨τ', hτ', rfl⟩ := List.mem_map.mp h
        simp only [scl] at hp
        simp only [List.mem_flatMap, List.mem_range, List.mem_singleton] at hi
        obtain ⟨i', -, rfl⟩ := hi
        rcases Finset.mem_union.mp (Finsupp.support_add hp) with hp | hp
        · exact Sh_poles v i' (by omega) τ' hτ' p hp
        · have := Finsupp.support_single_subset hp
          simp at this; omega
      · obtain ⟨i, hi, h⟩ := List.mem_flatMap.mp h
        obtain ⟨τ', hτ', rfl⟩ := List.mem_map.mp h
        simp only [scl] at hp
        simp only [List.mem_flatMap, List.mem_range, List.mem_singleton] at hi
        obtain ⟨i', hi', rfl⟩ := hi
        rcases Finset.mem_union.mp (Finsupp.support_add hp) with hp | hp
        · exact Sh_poles v (-((i' : ℕ) : ℤ) - 1) (by omega) τ' hτ' p hp
        · have := Finsupp.support_single_subset hp
          simp at this; omega

lemma Sh_pos (a : ℤ) (v : List ℤ) (m : ℕ) (hm : 0 < m) :
    Sh (a :: v) (m : ℤ) = (1, 0, a :: v) :: (List.range m).flatMap (fun (i : ℕ) =>
      (Sh v (i : ℤ)).map (scl ((Int.sign a : ℚ) ^ ((i : ℤ) + 1)) (-((i : ℤ) + 1)) a.natAbs)) := by
  have h : (0 : ℤ) < m := by exact_mod_cast hm
  rw [Sh]; split_ifs
  simp only [List.pure_def, List.bind_eq_flatMap, Int.toNat_natCast, List.flatMap_assoc, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]

lemma Sh_nonpos (a : ℤ) (v : List ℤ) (m : ℕ) :
    Sh (a :: v) (-(m : ℤ)) = (1, 0, a :: v) :: (List.range m).flatMap (fun (i : ℕ) =>
      (Sh v (-(i : ℤ) - 1)).map (scl (-(Int.sign a : ℚ) ^ (-(i : ℤ))) (i : ℤ) a.natAbs)) := by
  have h : ¬ (0 : ℤ) < -(m : ℤ) := by omega
  rw [Sh]; split_ifs
  simp only [List.pure_def, List.bind_eq_flatMap, neg_neg, Int.toNat_natCast, List.flatMap_assoc,
    List.flatMap_cons, List.flatMap_nil, List.append_nil]

lemma Sh_val : ∀ v : List ℤ, S03.IsIndex v → ∀ j : ℤ, -2 ≤ j → ∀ n : ℕ, 1 ≤ n →
    S03.ST (Sh v j) (2 * n) = S03.Z v (2 * (n : ℤ) + j).toNat
  | [], _, j, _, n, _ => by simp [Sh, ST_cons, ST_nil, ratEval_one, S03.Z]
  | a :: v, hv, j, hj, n, hn => by
    have ha : a ≠ 0 := hv a (by simp)
    have hv' : S03.IsIndex v := fun b hb => hv b (by simp [hb])
    have hσ := S03.sign_ne a ha
    rcases lt_or_ge 0 j with hj0 | hj0
    · obtain ⟨m, rfl⟩ : ∃ m : ℕ, j = m := ⟨j.toNat, by omega⟩
      have hm : 0 < m := by exact_mod_cast hj0
      rw [Sh_pos a v m hm, ST_cons, ST_flatMap, list_range_sum]
      simp only [ratEval_one, one_mul]
      rw [show (2 * (n : ℤ) + m).toNat = 2 * n + m by omega, Z_add]
      congr 1
      apply Finset.sum_congr rfl
      intro i _
      rw [ST_map_scl, Sh_val v hv' i (by omega) n hn,
        show (2 * (n : ℤ) + i).toNat = 2 * n + i by omega]
      rw [show (i : ℤ) + 1 = ((i + 1 : ℕ) : ℤ) by push_cast; ring, zpow_natCast,
        show 2 * n + i + 1 = 2 * n + (i + 1) by ring, sgn_even a ha]
      push_cast; ring_nf
    · obtain ⟨n', rfl⟩ : ∃ n' : ℕ, n = n' + 1 := ⟨n - 1, by omega⟩
      have h3 : j = -((0 : ℕ) : ℤ) ∨ j = -((1 : ℕ) : ℤ) ∨ j = -((2 : ℕ) : ℤ) := by omega
      rcases h3 with rfl | rfl | rfl
      · rw [Sh_nonpos, show (2 * ((n' + 1 : ℕ) : ℤ) + -((0 : ℕ) : ℤ)).toNat = 2 * (n' + 1) by omega]
        simp [ST_cons, ST_nil, ratEval_one]
      · rw [Sh_nonpos]
        simp only [List.range_one, List.flatMap_singleton, ST_cons, ST_map_scl, ratEval_one, one_mul]
        rw [Sh_val v hv' _ (by norm_num) _ hn]
        rw [show (2 * ((n' + 1 : ℕ) : ℤ) + -((1 : ℕ) : ℤ)).toNat = 2 * n' + 1 by omega,
          show (2 * ((n' + 1 : ℕ) : ℤ) + (-((0 : ℕ) : ℤ) - 1)).toNat = 2 * n' + 1 by omega,
          show 2 * (n' + 1) = 2 * n' + 1 + 1 by ring, Z_succ]
        rw [show 2 * n' + 1 + 1 = 2 * (n' + 1) + 0 by ring, sgn_even a ha]
        simp; ring
      · rw [Sh_nonpos]
        simp only [show List.range 2 = [0, 1] from rfl, List.flatMap_cons, List.flatMap_nil,
          List.append_nil, ST_cons, ST_append, ST_map_scl, ratEval_one, one_mul]
        rw [Sh_val v hv' _ (by norm_num) _ hn, Sh_val v hv' _ (by norm_num) _ hn]
        rw [show (2 * ((n' + 1 : ℕ) : ℤ) + -((2 : ℕ) : ℤ)).toNat = 2 * n' by omega,
          show (2 * ((n' + 1 : ℕ) : ℤ) + (-((0 : ℕ) : ℤ) - 1)).toNat = 2 * n' + 1 by omega,
          show (2 * ((n' + 1 : ℕ) : ℤ) + (-((1 : ℕ) : ℤ) - 1)).toNat = 2 * n' by omega,
          show 2 * (n' + 1) = 2 * n' + 1 + 1 by ring, Z_succ, Z_succ]
        rw [show 2 * n' + 1 + 1 = 2 * (n' + 1) + 0 by ring, sgn_even a ha,
          show 2 * n' + 1 = 2 * n' + 1 by rfl, sgn_even a ha]
        rcases hσ with h | h <;> simp [h] <;> ring

lemma hW_poles (w : List ℤ) (δ : ℤ) (hδ : -1 ≤ δ) : ∀ τ ∈ hW w δ, ∀ p ∈ τ.2.1.support, p ≤ 1 := by
  intro τ hτ p hp
  obtain ⟨τ', hτ', rfl⟩ := List.mem_map.mp hτ
  simp only [scl] at hp
  rcases Finset.mem_union.mp (Finsupp.support_add hp) with hp | hp
  · exact Sh_poles _ (δ - 1) (by omega) τ' hτ' p hp
  · have := Finsupp.support_single_subset hp
    simp at this; omega

lemma pm_zpow_two (ε : ℤ) (hε : ε = 1 ∨ ε = -1) (n : ℕ) (δ : ℤ) :
    (ε : ℚ) ^ (2 * (n : ℤ) + δ) = (ε : ℚ) ^ δ := by
  have h0 : (ε : ℚ) ≠ 0 := by rcases hε with rfl | rfl <;> norm_num
  rw [zpow_add₀ h0, zpow_mul, show ((ε : ℚ) ^ (2 : ℤ)) = 1 by
    rcases hε with rfl | rfl <;> norm_num, one_zpow, one_mul]

lemma hc_val (w : List ℤ) (hne : w ≠ []) (hl : S03.IsLetters w) (hr : S03.Regular w)
    (δ : ℤ) (hδ : -1 ≤ δ) (n : ℕ) (hn : 1 ≤ n) :
    hc w (2 * (n : ℤ) + δ).toNat = S03.ST (hW w δ) (2 * n) := by
  obtain ⟨hε, -, -, hv, -, -⟩ := S03.proof_5.2.1 w hne hl hr
  unfold hW
  rw [ST_map_scl, Sh_val _ hv (δ - 1) (by omega) n hn]
  unfold hc
  split_ifs with h0
  · omega
  have hM : (((2 * (n : ℤ) + δ).toNat : ℕ) : ℤ) = 2 * (n : ℤ) + δ := by omega
  rw [← zpow_natCast, hM, pm_zpow_two _ hε,
    show (2 * (n : ℤ) + δ).toNat - 1 = (2 * (n : ℤ) + (δ - 1)).toNat by omega]
  have hM' : (((2 * (n : ℤ) + δ).toNat : ℕ) : ℚ) = ((2 * n : ℕ) : ℚ) - ((-δ : ℤ) : ℚ) := by
    have := congrArg (fun z : ℤ => (z : ℚ)) hM
    push_cast at this ⊢; linarith
  rw [hM']

end P04

theorem P04.proof_1 : P04.stmt_1 := fun v hv j hj =>
  ⟨P04.Sh_poles v j hj, fun n hn => P04.Sh_val v hv j hj n hn⟩

theorem P04.proof_2 : P04.stmt_2 := fun w hne hl hr δ hδ =>
  ⟨P04.hW_poles w δ hδ, fun n hn => P04.hc_val w hne hl hr δ hδ n hn⟩

namespace P04

open Real

lemma bar_ok (w : List ℤ) (hne : w ≠ []) (hl : S03.IsLetters w) (hr : S03.Regular w) :
    bar w ≠ [] ∧ S03.IsLetters (bar w) ∧ S03.Regular (bar w) := by
  refine ⟨by simpa [bar] using hne, ?_, ?_⟩
  · intro a ha
    obtain ⟨b, hb, rfl⟩ := List.mem_map.mp ha
    rcases hl b hb with h | h | h <;> simp [h]
  · intro a ha
    simp only [bar, List.getLast?_map, Option.mem_def, Option.map_eq_some_iff] at ha
    obtain ⟨b, hb, rfl⟩ := ha
    exact neg_ne_zero.mpr (hr b hb)

/-! ### Grouping by words -/

lemma Dval_cons (τ : S03.DTerm) (L : List S03.DTerm) (u : List ℤ) (N : ℚ) :
    Dval (τ :: L) u N = (if τ.2.2 = u then S03.ratEval τ.1 τ.2.1 N else 0) + Dval L u N := by
  unfold Dval
  rw [List.filter_cons]
  split_ifs with h1 h2 h2 <;> simp_all

lemma ST_group (L : List S03.DTerm) (U : Finset (List ℤ)) (hU : ∀ τ ∈ L, τ.2.2 ∈ U) (N : ℕ) :
    S03.ST L N = ∑ u ∈ U, Dval L u N * S03.Z u N := by
  induction L with
  | nil => simp [S03.ST, Dval]
  | cons τ L ih =>
    rw [ST_cons, ih (fun τ' h => hU τ' (by simp [h]))]
    simp only [Dval_cons, add_mul, Finset.sum_add_distrib]
    congr 1
    simp only [ite_mul, zero_mul]
    rw [Finset.sum_ite_eq U τ.2.2]
    rw [if_pos (hU τ (by simp))]

lemma ST_of_Deq (L L' : List S03.DTerm) (h : Deq L L') (N : ℕ) (hN : 1 < N) :
    S03.ST L N = S03.ST L' N := by
  classical
  let U : Finset (List ℤ) := (L.map fun τ => τ.2.2).toFinset ∪ (L'.map fun τ => τ.2.2).toFinset
  rw [ST_group L U (fun τ hτ => by simp [U]; left; exact ⟨τ.1, τ.2.1, hτ⟩) N,
    ST_group L' U (fun τ hτ => by simp [U]; right; exact ⟨τ.1, τ.2.1, hτ⟩) N]
  apply Finset.sum_congr rfl
  intro u _
  rw [h u N (by exact_mod_cast hN)]

/-! ### The coefficients of `c^{2n+ι}` -/

lemma ST_Lι (T : List PTerm) (ι : ℤ) (N : ℕ) :
    S03.ST (Lι T ι) N = (T.map fun t => t.1 * S03.ST (hW t.2.2 (ι - t.2.1)) N +
      t.1 * sgnP t.2.1 t.2.2 * S03.ST (hW (bar t.2.2) (ι - t.2.1)) N).sum := by
  unfold Lι
  rw [ST_flatMap]
  congr 1
  apply List.map_congr_left
  intro t _
  rw [ST_append, ST_map_scl, ST_map_scl]
  simp

lemma list_sum_congr {β : Type*} (L : List β) (f g : β → ℚ) (h : ∀ x ∈ L, f x = g x) :
    (L.map f).sum = (L.map g).sum := by
  rw [List.map_congr_left h]

lemma coefL_eq (T : List PTerm) (C : List (ℚ × ℤ))
    (hT : ∀ t ∈ T, t.2.2 ≠ [] ∧ S03.IsLetters t.2.2 ∧ S03.Regular t.2.2 ∧ t.2.1 ≤ 1)
    (hC : ∀ t ∈ C, t.2 ≤ 1) (ι : ℤ) (hι : 0 ≤ ι) (n : ℕ) (hn : 1 ≤ n) :
    coefL T C (2 * (n : ℤ) + ι) = S03.ST (Lι T ι) (2 * n) := by
  unfold coefL
  rw [ST_Lι]
  have h0 : (C.map fun t => if t.2 = 2 * (n : ℤ) + ι then t.1 else 0).sum = 0 := by
    apply List.sum_eq_zero
    intro x hx
    obtain ⟨t, ht, rfl⟩ := List.mem_map.mp hx
    have := hC t ht
    rw [if_neg (by omega)]
  rw [h0, zero_add]
  apply list_sum_congr
  intro t ht
  obtain ⟨hne, hl, hr, he⟩ := hT t ht
  obtain ⟨hne', hl', hr'⟩ := bar_ok _ hne hl hr
  rw [if_pos (by omega)]
  have e1 : 2 * (n : ℤ) + ι - t.2.1 = 2 * (n : ℤ) + (ι - t.2.1) := by ring
  rw [e1, hc_val _ hne hl hr _ (by omega) n hn, hc_val _ hne' hl' hr' _ (by omega) n hn]
  ring

/-! ### The Laurent expansion of `Lfun` -/

lemma hc_hasSum (w : List ℤ) (hne : w ≠ []) (hl : S03.IsLetters w) (hr : S03.Regular w)
    (x : ℝ) (hx : |x| < 1) :
    HasSum (fun m : ℕ => ((hc w (m + 1) : ℚ) : ℝ) * x ^ (m + 1)) (S03.H w x) := by
  have h := (S03.proof_5.2.1 w hne hl hr).2.2.2.2.2 x hx
  convert h using 2 with m
  simp [hc]

lemma pair_hasSum (t : PTerm) (hne : t.2.2 ≠ []) (hl : S03.IsLetters t.2.2)
    (hr : S03.Regular t.2.2) (c : ℝ) (hc0 : 0 < c) (hc1 : c < 1) :
    HasSum (fun N : ℤ => (((if 1 ≤ N - t.2.1 then t.1 * (hc t.2.2 (N - t.2.1).toNat +
        sgnP t.2.1 t.2.2 * hc (bar t.2.2) (N - t.2.1).toNat) else 0 : ℚ) : ℝ) * c ^ N))
      ((t.1 : ℝ) * c ^ t.2.1 *
        (S03.H t.2.2 c + (sgnP t.2.1 t.2.2 : ℝ) * S03.H (bar t.2.2) c)) := by
  obtain ⟨q, e, w⟩ := t
  simp only at hne hl hr ⊢
  obtain ⟨hne', hl', hr'⟩ := bar_ok _ hne hl hr
  have hx : |c| < 1 := by rw [abs_of_pos hc0]; exact hc1
  have h1 := hc_hasSum w hne hl hr c hx
  have h2 := hc_hasSum (bar w) hne' hl' hr' c hx
  have hG := ((h1.add (h2.mul_left (sgnP e w : ℝ))).mul_left ((q : ℝ) * c ^ e))
  let φ : ℕ → ℤ := fun m => (m : ℤ) + 1 + e
  have hφ : Function.Injective φ := by intro a b h; simp [φ] at h; exact h
  rw [← hφ.hasSum_iff]
  · convert hG using 1
    funext m
    simp only [Function.comp, φ]
    rw [if_pos (by omega), show (m : ℤ) + 1 + e - e = ((m + 1 : ℕ) : ℤ) by push_cast; ring,
      Int.toNat_natCast]
    rw [show (m : ℤ) + 1 + e = e + ((m + 1 : ℕ) : ℤ) by push_cast; ring,
      zpow_add₀ hc0.ne', zpow_natCast]
    push_cast; ring
  · intro N hN
    rw [if_neg]
    · simp
    · intro h
      exact hN ⟨(N - e - 1).toNat, by simp only [φ]; omega⟩

lemma const_hasSum (t : ℚ × ℤ) (c : ℝ) :
    HasSum (fun N : ℤ => (((if t.2 = N then t.1 else 0 : ℚ) : ℝ) * c ^ N)) ((t.1 : ℝ) * c ^ t.2) := by
  have := hasSum_single (f := fun N : ℤ => (((if t.2 = N then t.1 else 0 : ℚ) : ℝ) * c ^ N)) t.2
    (fun b hb => by rw [if_neg (Ne.symm hb)]; simp)
  simpa using this

lemma Lfun_hasSum (T : List PTerm) (C : List (ℚ × ℤ))
    (hT : ∀ t ∈ T, t.2.2 ≠ [] ∧ S03.IsLetters t.2.2 ∧ S03.Regular t.2.2 ∧ t.2.1 ≤ 1)
    (c : ℝ) (hc0 : 0 < c) (hc1 : c < 1) :
    HasSum (fun N : ℤ => ((coefL T C N : ℚ) : ℝ) * c ^ N) (Lfun T C c) := by
  have hCs : ∀ C : List (ℚ × ℤ), HasSum (fun N : ℤ =>
      (((C.map fun t => if t.2 = N then t.1 else 0).sum : ℚ) : ℝ) * c ^ N)
      (C.map fun t => (t.1 : ℝ) * c ^ t.2).sum := by
    intro C
    induction C with
    | nil => simp
    | cons t C ih =>
      simp only [List.map_cons, List.sum_cons, Rat.cast_add, add_mul]
      exact (const_hasSum t c).add ih
  have hTs : ∀ T : List PTerm,
      (∀ t ∈ T, t.2.2 ≠ [] ∧ S03.IsLetters t.2.2 ∧ S03.Regular t.2.2 ∧ t.2.1 ≤ 1) →
      HasSum (fun N : ℤ => (((T.map fun t => if 1 ≤ N - t.2.1 then
        t.1 * (hc t.2.2 (N - t.2.1).toNat + sgnP t.2.1 t.2.2 * hc (bar t.2.2) (N - t.2.1).toNat)
        else 0).sum : ℚ) : ℝ) * c ^ N)
      (T.map fun t => (t.1 : ℝ) * c ^ t.2.1 *
        (S03.H t.2.2 c + (sgnP t.2.1 t.2.2 : ℝ) * S03.H (bar t.2.2) c)).sum := by
    intro T
    induction T with
    | nil => intro _; simp
    | cons t T ih =>
      intro h
      simp only [List.map_cons, List.sum_cons, Rat.cast_add, add_mul]
      obtain ⟨h1, h2, h3, -⟩ := h t (by simp)
      exact (pair_hasSum t h1 h2 h3 c hc0 hc1).add (ih (fun t' ht' => h t' (by simp [ht'])))
  have := (hCs C).add (hTs T hT)
  unfold coefL Lfun
  simpa only [Rat.cast_add, add_mul] using this

end P04

theorem P04.proof_3 : P04.stmt_3 := by
  intro k s Tz T C hTz hT hC h0 h1 hneg c hc0 hc1
  have hL := P04.Lfun_hasSum T C hT c hc0 hc1
  let ψ : ℕ → ℤ := fun n => 2 * ((n : ℤ) + 1)
  have hψ : Function.Injective ψ := by intro a b h; simp [ψ] at h; exact h
  rw [← hψ.hasSum_iff] at hL
  · convert hL using 1
    funext n
    simp only [Function.comp, ψ]
    have e1 : 2 * ((n : ℤ) + 1) = 2 * (((n + 1 : ℕ) : ℤ)) + 0 := by push_cast; ring
    rw [e1, P04.coefL_eq T C hT hC 0 le_rfl (n + 1) (by omega),
      P04.ST_of_Deq _ _ h0 _ (by omega), ← hTz (n + 1) (by omega)]
    rw [show 2 * (((n + 1 : ℕ) : ℤ)) + 0 = ((2 * (n + 1) : ℕ) : ℤ) by push_cast; ring,
      zpow_natCast]
  · intro N hN
    have hcoef : P04.coefL T C N = 0 := by
      rcases le_or_gt N 1 with hN1 | hN1
      · exact hneg N hN1
      · obtain ⟨m, hm⟩ := Int.even_or_odd' N
        rcases hm with hm | hm
        · exact absurd ⟨(m - 1).toNat, by simp only [ψ]; omega⟩ hN
        · obtain ⟨n, rfl⟩ : ∃ n : ℕ, m = n := ⟨m.toNat, by omega⟩
          rw [hm, P04.coefL_eq T C hT hC 1 (by norm_num) n (by omega),
            P04.ST_of_Deq _ _ h1 _ (by omega)]
          simp [S03.ST]
    rw [hcoef]; simp

namespace P04

lemma series_of_coef (T : List PTerm) (C : List (ℚ × ℤ))
    (hT : ∀ t ∈ T, t.2.2 ≠ [] ∧ S03.IsLetters t.2.2 ∧ S03.Regular t.2.2 ∧ t.2.1 ≤ 1)
    (p : ℕ → ℚ) (heven : ∀ n : ℕ, 1 ≤ n → coefL T C (2 * (n : ℤ)) = p n)
    (hodd : ∀ n : ℕ, 1 ≤ n → coefL T C (2 * (n : ℤ) + 1) = 0)
    (hneg : ∀ N : ℤ, N ≤ 1 → coefL T C N = 0) (c : ℝ) (hc0 : 0 < c) (hc1 : c < 1) :
    HasSum (fun n : ℕ => ((p (n + 1) : ℚ) : ℝ) * c ^ (2 * (n + 1))) (Lfun T C c) := by
  have hL := Lfun_hasSum T C hT c hc0 hc1
  let ψ : ℕ → ℤ := fun n => 2 * ((n : ℤ) + 1)
  have hψ : Function.Injective ψ := by intro a b h; simp [ψ] at h; exact h
  rw [← hψ.hasSum_iff] at hL
  · convert hL using 1
    funext n
    simp only [Function.comp, ψ]
    rw [show 2 * ((n : ℤ) + 1) = 2 * (((n + 1 : ℕ) : ℤ)) by push_cast; ring,
      heven (n + 1) (by omega),
      show 2 * (((n + 1 : ℕ) : ℤ)) = ((2 * (n + 1) : ℕ) : ℤ) by push_cast; ring, zpow_natCast]
  · intro N hN
    have hcoef : coefL T C N = 0 := by
      rcases le_or_gt N 1 with hN1 | hN1
      · exact hneg N hN1
      · obtain ⟨m, hm⟩ := Int.even_or_odd' N
        rcases hm with hm | hm
        · exact absurd ⟨(m - 1).toNat, by simp only [ψ]; omega⟩ hN
        · obtain ⟨n, rfl⟩ : ∃ n : ℕ, m = n := ⟨m.toNat, by omega⟩
          rw [hm]; exact hodd n (by omega)
    rw [hcoef]; simp

lemma Psi2T_ok : ∀ t ∈ Psi2T, t.2.2 ≠ [] ∧ S03.IsLetters t.2.2 ∧ S03.Regular t.2.2 ∧ t.2.1 ≤ 1 := by
  simp only [Psi2T, S03.IsLetters, S03.Regular]
  decide

lemma R4T_ok : ∀ t ∈ R4T, t.2.2 ≠ [] ∧ S03.IsLetters t.2.2 ∧ S03.Regular t.2.2 ∧ t.2.1 ≤ 1 := by
  simp only [R4T, S03.IsLetters, S03.Regular]
  decide

end P04

namespace P04

lemma Z_zero (v : List ℤ) : S03.Z v 0 = if v = [] then 1 else 0 := by
  cases v <;> simp [S03.Z]

lemma sgn_lit : Int.sign 2 = 1 ∧ Int.sign 3 = 1 ∧ Int.sign (-2) = -1 ∧ Int.sign (-3) = -1 := by
  decide

lemma Icc_lit : Finset.Icc 1 2 = {1, 2} ∧ Finset.Icc 1 3 = {1, 2, 3} ∧
    Finset.Icc 1 4 = {1, 2, 3, 4} ∧ Finset.Icc 1 5 = {1, 2, 3, 4, 5} := by
  decide

lemma coefL_low (T : List PTerm) (C : List (ℚ × ℤ)) (N : ℤ)
    (hT : ∀ t ∈ T, N < t.2.1 + 1) (hC : ∀ t ∈ C, t.2 ≠ N) : coefL T C N = 0 := by
  unfold coefL
  rw [List.sum_eq_zero, List.sum_eq_zero, add_zero]
  · intro x hx
    obtain ⟨t, ht, rfl⟩ := List.mem_map.mp hx
    have := hT t ht
    rw [if_neg (by omega)]
  · intro x hx
    obtain ⟨t, ht, rfl⟩ := List.mem_map.mp hx
    rw [if_neg (hC t ht)]

lemma Psi2T_e : ∀ t ∈ Psi2T, -3 ≤ t.2.1 := by simp only [Psi2T]; decide
lemma R4T_e : ∀ t ∈ R4T, -3 ≤ t.2.1 := by simp only [R4T]; decide

lemma low_cases (P : ℤ → Prop) (h3 : ∀ N : ℤ, N ≤ -3 → P N) (hm2 : P (-2)) (hm1 : P (-1))
    (h0 : P 0) (h1 : P 1) : ∀ N : ℤ, N ≤ 1 → P N := by
  intro N hN
  rcases le_or_gt N (-3) with h | h
  · exact h3 N h
  · have : N = -2 ∨ N = -1 ∨ N = 0 ∨ N = 1 := by omega
    rcases this with rfl | rfl | rfl | rfl <;> assumption

lemma Psi2_low : ∀ N : ℤ, N ≤ 1 → coefL Psi2T Psi2C N = 0 := by
  apply low_cases
  · intro N hN
    apply coefL_low
    · intro t ht; have := Psi2T_e t ht; omega
    · intro t ht; simp [Psi2C] at ht; subst ht; simp; omega
  all_goals
    simp [coefL, Psi2T, Psi2C, hc, sgnP, nu, bar, S03.dd, Z_succ, Z_zero, S03.Z, sgn_lit, Icc_lit] <;> norm_num

lemma R4_low : ∀ N : ℤ, N ≤ 1 → coefL R4T [] N = 0 := by
  apply low_cases
  · intro N hN
    apply coefL_low
    · intro t ht; have := R4T_e t ht; omega
    · simp
  all_goals
    simp [coefL, R4T, hc, sgnP, nu, bar, S03.dd, Z_succ, Z_zero, S03.Z, sgn_lit, Icc_lit] <;> norm_num

end P04


namespace P04

/-! ### Explicit stuffle and products in `𝔇` -/

/-- the letter `a ⊕ b` of the quasi-shuffle -/
def opl (a b : ℤ) : ℤ := Int.sign a * Int.sign b * ((a.natAbs + b.natAbs : ℕ) : ℤ)

/-- the quasi-shuffle (stuffle) of two words, as a list of words -/
def stuf : List ℤ → List ℤ → List (List ℤ)
  | [], v => [v]
  | a :: u, [] => [a :: u]
  | a :: u, b :: v => (stuf u (b :: v)).map (a :: ·) ++ (stuf (a :: u) v).map (b :: ·) ++
      (stuf u v).map (opl a b :: ·)
termination_by u v => u.length + v.length

lemma stuf_nil_left (v : List ℤ) : stuf [] v = [v] := by rw [stuf]
lemma stuf_nil_right (a : ℤ) (u : List ℤ) : stuf (a :: u) [] = [a :: u] := by rw [stuf]
lemma stuf_cons (a b : ℤ) (u v : List ℤ) : stuf (a :: u) (b :: v) =
    (stuf u (b :: v)).map (a :: ·) ++ (stuf (a :: u) v).map (b :: ·) ++
      (stuf u v).map (opl a b :: ·) := by rw [stuf]

lemma Z_stuf : ∀ u v : List ℤ, S03.IsIndex u → S03.IsIndex v → ∀ M : ℕ,
    S03.Z u M * S03.Z v M = ((stuf u v).map fun w => S03.Z w M).sum := by
  intro u
  induction u with
  | nil => intro v _ _ M; simp [stuf_nil_left, S03.Z]
  | cons a u ihu =>
    intro v hu hv
    have ha : a ≠ 0 := hu a (by simp)
    have hu' : S03.IsIndex u := fun x hx => hu x (by simp [hx])
    induction v with
    | nil => intro M; simp [stuf_nil_right, S03.Z]
    | cons b v ihv =>
      intro M
      have hb : b ≠ 0 := hv b (by simp)
      have hv' : S03.IsIndex v := fun x hx => hv x (by simp [hx])
      have hg := fun i => (S03.gq_mul a b ha hb i).1
      rw [stuf_cons, List.map_append, List.map_append, List.sum_append, List.sum_append,
        List.map_map, List.map_map, List.map_map]
      rw [show ((fun w => S03.Z w M) ∘ (a :: ·)) = fun w => S03.Z (a :: w) M from rfl,
        show ((fun w => S03.Z w M) ∘ (b :: ·)) = fun w => S03.Z (b :: w) M from rfl,
        show ((fun w => S03.Z w M) ∘ (opl a b :: ·)) = fun w => S03.Z (opl a b :: w) M from rfl,
        S03.sum_map_cons, S03.sum_map_cons, S03.sum_map_cons, S03.stuffle_core]
      have e1 : ∀ i, S03.Z u (i - 1) * S03.Z (b :: v) (i - 1) =
          ((stuf u (b :: v)).map fun w => S03.Z w (i - 1)).sum := fun i => ihu (b :: v) hu' hv _
      have e2 : ∀ i, S03.Z (a :: u) (i - 1) * S03.Z v (i - 1) =
          ((stuf (a :: u) v).map fun w => S03.Z w (i - 1)).sum := fun i => ihv hv' _
      have e3 : ∀ i, S03.Z u (i - 1) * S03.Z v (i - 1) =
          ((stuf u v).map fun w => S03.Z w (i - 1)).sum := fun i => ihu v hu' hv' _
      simp_rw [e1, e2, e3, hg]
      rfl

lemma stuf_index : ∀ u v : List ℤ, S03.IsIndex u → S03.IsIndex v → ∀ w ∈ stuf u v, S03.IsIndex w := by
  intro u
  induction u with
  | nil => intro v _ hv w hw; rw [stuf_nil_left] at hw; simp at hw; subst hw; exact hv
  | cons a u ihu =>
    intro v hu hv
    have ha : a ≠ 0 := hu a (by simp)
    have hu' : S03.IsIndex u := fun x hx => hu x (by simp [hx])
    induction v with
    | nil => intro w hw; rw [stuf_nil_right] at hw; simp at hw; subst hw; exact hu
    | cons b v ihv =>
      intro w hw
      have hb : b ≠ 0 := hv b (by simp)
      have hv' : S03.IsIndex v := fun x hx => hv x (by simp [hx])
      have hc := (S03.gq_mul a b ha hb 0).2.2
      rw [stuf_cons] at hw
      simp only [List.mem_append, List.mem_map] at hw
      rcases hw with (⟨w', hw', rfl⟩ | ⟨w', hw', rfl⟩) | ⟨w', hw', rfl⟩
      · intro x hx; simp at hx; rcases hx with rfl | hx; exact ha; exact ihu _ hu' hv w' hw' x hx
      · intro x hx; simp at hx; rcases hx with rfl | hx; exact hb; exact ihv hv' w' hw' x hx
      · intro x hx; simp at hx; rcases hx with rfl | hx; exact hc; exact ihu _ hu' hv' w' hw' x hx

/-- the product in `𝔇` -/
noncomputable def mulL (L1 L2 : List S03.DTerm) : List S03.DTerm :=
  L1.flatMap fun τ1 => L2.flatMap fun τ2 =>
    (stuf τ1.2.2 τ2.2.2).map fun w => (τ1.1 * τ2.1, τ1.2.1 + τ2.2.1, w)

lemma ST_stufmap (Q : Polynomial ℚ) (e : ℤ →₀ ℕ) (u v : List ℤ) (hu : S03.IsIndex u)
    (hv : S03.IsIndex v) (N : ℕ) :
    S03.ST ((stuf u v).map fun w => (Q, e, w)) N = S03.ratEval Q e N * (S03.Z u N * S03.Z v N) := by
  rw [Z_stuf u v hu hv N]
  unfold S03.ST
  rw [List.map_map, ← List.sum_map_mul_left]
  rfl

lemma ST_mulL (L1 L2 : List S03.DTerm) (h1 : ∀ τ ∈ L1, S03.IsIndex τ.2.2)
    (h2 : ∀ τ ∈ L2, S03.IsIndex τ.2.2) (N : ℕ) :
    S03.ST (mulL L1 L2) N = S03.ST L1 N * S03.ST L2 N := by
  induction L1 with
  | nil => simp [mulL, S03.ST]
  | cons τ1 L1 ih =>
    have ih' := ih (fun τ h => h1 τ (by simp [h]))
    unfold mulL at ih' ⊢
    rw [List.flatMap_cons, ST_append, ih', ST_cons, add_mul]
    congr 1
    have hτ1 := h1 τ1 (by simp)
    clear ih ih' h1
    induction L2 with
    | nil => simp [S03.ST]
    | cons τ2 L2 ih2 =>
      rw [List.flatMap_cons, ST_append, ih2 (fun τ h => h2 τ (by simp [h])), ST_cons,
        ST_stufmap _ _ _ _ hτ1 (h2 τ2 (by simp)), S08.PA.ratEval_mul]
      ring

end P04

namespace P04

lemma Dval_nil (u : List ℤ) (N : ℚ) : Dval [] u N = 0 := by simp [Dval]

lemma Dval_append (L L' : List S03.DTerm) (u : List ℤ) (N : ℚ) :
    Dval (L ++ L') u N = Dval L u N + Dval L' u N := by
  simp [Dval, List.filter_append]

lemma Dval_flatMap {α : Type*} (l : List α) (f : α → List S03.DTerm) (u : List ℤ) (N : ℚ) :
    Dval (l.flatMap f) u N = (l.map fun x => Dval (f x) u N).sum := by
  induction l with
  | nil => simp [Dval]
  | cons x l ih => rw [List.flatMap_cons, Dval_append, ih, List.map_cons, List.sum_cons]

lemma scl_w (r : ℚ) (p : ℤ) (k : ℕ) (τ : S03.DTerm) : (scl r p k τ).2.2 = τ.2.2 := rfl

lemma Dval_map_scl (r : ℚ) (p : ℤ) (k : ℕ) (L : List S03.DTerm) (u : List ℤ) (N : ℚ) :
    Dval (L.map (scl r p k)) u N = r / (N - p) ^ k * Dval L u N := by
  induction L with
  | nil => simp [Dval]
  | cons τ L ih =>
    rw [List.map_cons, Dval_cons, Dval_cons, ih, ratEval_scl, scl_w]
    by_cases h : τ.2.2 = u <;> simp only [h, if_true, if_false] <;> ring

lemma Dval_stufmap (Q : Polynomial ℚ) (e : ℤ →₀ ℕ) (u v x : List ℤ) (N : ℚ) :
    Dval ((stuf u v).map fun w => (Q, e, w)) x N =
      S03.ratEval Q e N * ((stuf u v).map fun w => if w = x then (1 : ℚ) else 0).sum := by
  induction (stuf u v) with
  | nil => simp [Dval]
  | cons w L ih =>
    rw [List.map_cons, Dval_cons, ih, List.map_cons, List.sum_cons]
    simp only
    split_ifs <;> ring

lemma Dval_mulL (L1 L2 : List S03.DTerm) (x : List ℤ) (N : ℚ) :
    Dval (mulL L1 L2) x N = (L1.map fun τ1 => (L2.map fun τ2 =>
      S03.ratEval τ1.1 τ1.2.1 N * S03.ratEval τ2.1 τ2.2.1 N *
        ((stuf τ1.2.2 τ2.2.2).map fun w => if w = x then (1 : ℚ) else 0).sum).sum).sum := by
  unfold mulL
  rw [Dval_flatMap]
  congr 1
  apply List.map_congr_left
  intro τ1 _
  rw [Dval_flatMap]
  congr 1
  apply List.map_congr_left
  intro τ2 _
  rw [Dval_stufmap, S08.PA.ratEval_mul]


lemma Sh_nil' (j : ℤ) : Sh [] j = [(1, 0, [])] := by rw [Sh]

lemma Sh_0 (a : ℤ) (v : List ℤ) : Sh (a :: v) 0 = [(1, 0, a :: v)] := by
  have := Sh_nonpos a v 0; simpa using this

lemma Sh_1 (a : ℤ) (v : List ℤ) : Sh (a :: v) 1 = (1, 0, a :: v) ::
    (Sh v 0).map (scl (Int.sign a : ℚ) (-1) a.natAbs) := by
  have := Sh_pos a v 1 (by norm_num); simpa using this

lemma Sh_2 (a : ℤ) (v : List ℤ) : Sh (a :: v) 2 = (1, 0, a :: v) ::
    ((Sh v 0).map (scl (Int.sign a : ℚ) (-1) a.natAbs) ++
     (Sh v 1).map (scl ((Int.sign a : ℚ) ^ 2) (-2) a.natAbs)) := by
  have := Sh_pos a v 2 (by norm_num)
  simp only [show List.range 2 = [0, 1] from rfl, List.flatMap_cons, List.flatMap_nil,
    List.append_nil] at this
  simpa using this

lemma Sh_3 (a : ℤ) (v : List ℤ) : Sh (a :: v) 3 = (1, 0, a :: v) ::
    ((Sh v 0).map (scl (Int.sign a : ℚ) (-1) a.natAbs) ++
     ((Sh v 1).map (scl ((Int.sign a : ℚ) ^ 2) (-2) a.natAbs) ++
     (Sh v 2).map (scl ((Int.sign a : ℚ) ^ 3) (-3) a.natAbs))) := by
  have := Sh_pos a v 3 (by norm_num)
  simp only [show List.range 3 = [0, 1, 2] from rfl, List.flatMap_cons, List.flatMap_nil,
    List.append_nil] at this
  simpa using this

lemma Sh_m1 (a : ℤ) (v : List ℤ) : Sh (a :: v) (-1) = (1, 0, a :: v) ::
    (Sh v (-1)).map (scl (-1) 0 a.natAbs) := by
  have := Sh_nonpos a v 1; simpa using this

lemma Sh_m2 (a : ℤ) (v : List ℤ) : Sh (a :: v) (-2) = (1, 0, a :: v) ::
    ((Sh v (-1)).map (scl (-1) 0 a.natAbs) ++
     (Sh v (-2)).map (scl (-(Int.sign a : ℚ)⁻¹) 1 a.natAbs)) := by
  have := Sh_nonpos a v 2
  simp only [show List.range 2 = [0, 1] from rfl, List.flatMap_cons, List.flatMap_nil,
    List.append_nil] at this
  simpa using this

lemma Dval_base (Q : Polynomial ℚ) (e : ℤ →₀ ℕ) (w u : List ℤ) (L : List S03.DTerm) (N : ℚ) :
    Dval ((Q, e, w) :: L) u N = (if w = u then S03.ratEval Q e N else 0) + Dval L u N := Dval_cons _ _ _ _

lemma ratEval_10 (N : ℚ) : S03.ratEval 1 0 N = 1 := ratEval_one N

end P04

namespace P04

lemma Z_nil (M : ℕ) : S03.Z [] M = 1 := rfl

/-- the polynomial `c_0 + c_1 X + c_2 X^2 + ⋯` -/
noncomputable def polyOf : List ℚ → Polynomial ℚ
  | [] => 0
  | c :: cs => Polynomial.C c + Polynomial.X * polyOf cs

/-- Horner evaluation -/
def horner : List ℚ → ℚ → ℚ
  | [], _ => 0
  | c :: cs, N => c + N * horner cs N

lemma eval_polyOf (cs : List ℚ) (N : ℚ) : (polyOf cs).eval N = horner cs N := by
  induction cs with
  | nil => simp [polyOf, horner]
  | cons c cs ih => simp [polyOf, horner, ih]

/-- a rational term `P(N) / Π (N - p)^k` with the empty word, `P` given by its coefficients -/
noncomputable def rt (cs : List ℚ) (ps : List (ℤ × ℕ)) : S03.DTerm :=
  ps.foldr (fun pk τ => scl 1 pk.1 pk.2 τ) (polyOf cs, 0, [])

lemma rt_w (cs : List ℚ) (ps : List (ℤ × ℕ)) : (rt cs ps).2.2 = [] := by
  induction ps with
  | nil => rfl
  | cons pk ps ih => exact ih

lemma rt_eval (cs : List ℚ) (ps : List (ℤ × ℕ)) (N : ℚ) :
    S03.ratEval (rt cs ps).1 (rt cs ps).2.1 N =
      horner cs N / (ps.map fun pk => (N - pk.1) ^ pk.2).prod := by
  induction ps with
  | nil => simp [rt, S03.ratEval, eval_polyOf]
  | cons pk ps ih =>
    rw [show rt cs (pk :: ps) = scl 1 pk.1 pk.2 (rt cs ps) from rfl, ratEval_scl, ih, List.map_cons,
      List.prod_cons]
    ring

lemma Hn1 (r n : ℕ) (hr : 1 ≤ r) (hn : 1 ≤ n) (j : ℕ) : V05.HrQ r (n + j - 1) =
    (2 : ℚ) ^ (r - 1) * (S03.Z [(r : ℤ)] (2 * n) + S03.Z [-(r : ℤ)] (2 * n)) - 1 / (n : ℚ) ^ r +
      ∑ t ∈ Finset.range j, 1 / ((n : ℚ) + t) ^ r := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  rw [show m + 1 + j - 1 = m + j by omega, S08.PA.HrQ_shift r m j, ← S08.PA.hn_eq_HrQ]
  have h3 := (S03.proof3 r (m + 1) hr (by omega)).1
  rw [show m + 1 - 1 = m by omega] at h3
  rw [← h3]
  push_cast
  have : ∀ l ∈ Finset.range j, 1 / ((m : ℚ) + l + 1) ^ r = 1 / ((m : ℚ) + 1 + l) ^ r :=
    fun l _ => by ring_nf
  rw [Finset.sum_congr rfl this]

lemma Hn3 (r n : ℕ) (hr : 1 ≤ r) (hn : 1 ≤ n) (j : ℕ) : V05.HrQ r (n + j + n - 1) =
    S03.Z [(r : ℤ)] (2 * n) - 1 / (2 * (n : ℚ)) ^ r + ∑ l ∈ Finset.range j, 1 / (2 * (n : ℚ) + l) ^ r := by
  rw [show n + j + n - 1 = (2 * n - 1) + j by omega, S08.PA.HrQ_shift r (2 * n - 1) j, ← S08.PA.hn_eq_HrQ]
  have h3 := (S03.proof3 r n hr hn).2
  rw [← h3]
  have h2 : ((2 * n - 1 : ℕ) : ℚ) + 1 = 2 * (n : ℚ) := by
    rw [Nat.cast_sub (by omega)]; push_cast; ring
  have : ∀ l ∈ Finset.range j, 1 / (((2 * n - 1 : ℕ) : ℚ) + l + 1) ^ r = 1 / (2 * (n : ℚ) + l) ^ r :=
    fun l _ => by rw [add_right_comm, h2]
  rw [Finset.sum_congr rfl this]

lemma Hn2 (r n : ℕ) (hr : 1 ≤ r) (hn : 1 ≤ n) : V05.HrQ r (n - 1) =
    (2 : ℚ) ^ (r - 1) * (S03.Z [(r : ℤ)] (2 * n) + S03.Z [-(r : ℤ)] (2 * n)) - 1 / (n : ℚ) ^ r := by
  have := Hn1 r n hr hn 0
  simpa using this

/-- the atoms `α_r(j)` (k = 0), `β_r(j)` (k = 1), `γ_r(j)` (k = 2): word parts and one rational term -/
noncomputable def atm : ℕ → ℕ → ℕ → List S03.DTerm
  | 0, 1, 0 => [(1, 0, [-1]), rt [-1] [(0, 1)]]
  | 0, 1, 1 => [(1, 0, [-1]), rt [0] [(0, 1)]]
  | 0, 1, 2 => [(1, 0, [-1]), rt [0, 0, 1] [(-2, 1), (-1, 1), (0, 1)]]
  | 0, 1, 3 => [(1, 0, [-1]), rt [0, 0, 5, 2] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]]
  | 0, 1, 4 => [(1, 0, [-1]), rt [0, 0, 98, 95, 30, 3] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]]
  | 0, 2, 0 => [(1, 0, [2]), scl 2 0 0 (1, 0, [-2]), rt [-3] [(0, 2)]]
  | 0, 2, 1 => [(1, 0, [2]), scl 2 0 0 (1, 0, [-2]), rt [0] [(0, 2)]]
  | 0, 2, 2 => [(1, 0, [2]), scl 2 0 0 (1, 0, [-2]), rt [0, 0, 0, 4, 3] [(-2, 2), (-1, 2), (0, 2)]]
  | 0, 2, 3 => [(1, 0, [2]), scl 2 0 0 (1, 0, [-2]), rt [0, 0, 0, 72, 99, 42, 6] [(-4, 2), (-2, 2), (-1, 2), (0, 2)]]
  | 0, 2, 4 => [(1, 0, [2]), scl 2 0 0 (1, 0, [-2]), rt [0, 0, 0, 24096, 58284, 58476, 31899, 10272, 1956, 204, 9] [(-6, 2), (-4, 2), (-3, 2), (-2, 2), (-1, 2), (0, 2)]]
  | 0, 3, 0 => [scl 3 0 0 (1, 0, [3]), scl 4 0 0 (1, 0, [-3]), rt [-7] [(0, 3)]]
  | 0, 3, 1 => [scl 3 0 0 (1, 0, [3]), scl 4 0 0 (1, 0, [-3]), rt [0] [(0, 3)]]
  | 0, 3, 2 => [scl 3 0 0 (1, 0, [3]), scl 4 0 0 (1, 0, [-3]), rt [0, 0, 0, 0, 12, 18, 7] [(-2, 3), (-1, 3), (0, 3)]]
  | 0, 3, 3 => [scl 3 0 0 (1, 0, [3]), scl 4 0 0 (1, 0, [-3]), rt [0, 0, 0, 0, 816, 1908, 1715, 741, 159, 14] [(-4, 3), (-2, 3), (-1, 3), (0, 3)]]
  | 0, 3, 4 => [scl 3 0 0 (1, 0, [3]), scl 4 0 0 (1, 0, [-3]), rt [0, 0, 0, 0, 4814208, 18583776, 31857320, 32083044, 21164814, 9637481, 3100530, 706449, 111918, 11757, 738, 21] [(-6, 3), (-4, 3), (-3, 3), (-2, 3), (-1, 3), (0, 3)]]
  | 0, 4, 0 => [scl 7 0 0 (1, 0, [4]), scl 8 0 0 (1, 0, [-4]), rt [-15] [(0, 4)]]
  | 0, 4, 1 => [scl 7 0 0 (1, 0, [4]), scl 8 0 0 (1, 0, [-4]), rt [0] [(0, 4)]]
  | 0, 4, 2 => [scl 7 0 0 (1, 0, [4]), scl 8 0 0 (1, 0, [-4]), rt [0, 0, 0, 0, 0, 32, 72, 56, 15] [(-2, 4), (-1, 4), (0, 4)]]
  | 0, 4, 3 => [scl 7 0 0 (1, 0, [4]), scl 8 0 0 (1, 0, [-4]), rt [0, 0, 0, 0, 0, 8448, 27936, 38640, 28815, 12540, 3234, 468, 30] [(-4, 4), (-2, 4), (-1, 4), (0, 4)]]
  | 0, 4, 4 => [scl 7 0 0 (1, 0, [4]), scl 8 0 0 (1, 0, [-4]), rt [0, 0, 0, 0, 0, 890376192, 4733710848, 11545504512, 17132667888, 17310484512, 12632131848, 6890649432, 2867073183, 919403952, 227701188, 43280256, 6210984, 652416, 47424, 2136, 45] [(-6, 4), (-4, 4), (-3, 4), (-2, 4), (-1, 4), (0, 4)]]
  | 1, 1, 0 => [(1, 0, [-1]), rt [-1] [(0, 1)]]
  | 1, 1, 1 => [(1, 0, [-1]), rt [-2] [(0, 1)]]
  | 1, 1, 2 => [(1, 0, [-1]), rt [-2, -3] [(-1, 1), (0, 1)]]
  | 1, 1, 3 => [(1, 0, [-1]), rt [-4, -9, -4] [(-2, 1), (-1, 1), (0, 1)]]
  | 1, 1, 4 => [(1, 0, [-1]), rt [-12, -33, -24, -5] [(-3, 1), (-2, 1), (-1, 1), (0, 1)]]
  | 1, 2, 0 => [(1, 0, [2]), scl 2 0 0 (1, 0, [-2]), rt [-3] [(0, 2)]]
  | 1, 2, 1 => [(1, 0, [2]), scl 2 0 0 (1, 0, [-2]), rt [-4] [(0, 2)]]
  | 1, 2, 2 => [(1, 0, [2]), scl 2 0 0 (1, 0, [-2]), rt [-4, -8, -5] [(-1, 2), (0, 2)]]
  | 1, 2, 3 => [(1, 0, [2]), scl 2 0 0 (1, 0, [-2]), rt [-16, -48, -57, -30, -6] [(-2, 2), (-1, 2), (0, 2)]]
  | 1, 2, 4 => [(1, 0, [2]), scl 2 0 0 (1, 0, [-2]), rt [-144, -528, -821, -672, -304, -72, -7] [(-3, 2), (-2, 2), (-1, 2), (0, 2)]]
  | 1, 3, 0 => [scl 3 0 0 (1, 0, [3]), scl 4 0 0 (1, 0, [-3]), rt [-7] [(0, 3)]]
  | 1, 3, 1 => [scl 3 0 0 (1, 0, [3]), scl 4 0 0 (1, 0, [-3]), rt [-8] [(0, 3)]]
  | 1, 3, 2 => [scl 3 0 0 (1, 0, [3]), scl 4 0 0 (1, 0, [-3]), rt [-8, -24, -24, -9] [(-1, 3), (0, 3)]]
  | 1, 3, 3 => [scl 3 0 0 (1, 0, [3]), scl 4 0 0 (1, 0, [-3]), rt [-64, -288, -528, -513, -279, -81, -10] [(-2, 3), (-1, 3), (0, 3)]]
  | 1, 3, 4 => [scl 3 0 0 (1, 0, [3]), scl 4 0 0 (1, 0, [-3]), rt [-1728, -9504, -22608, -30771, -26460, -14931, -5544, -1311, -180, -11] [(-3, 3), (-2, 3), (-1, 3), (0, 3)]]
  | 1, 4, 0 => [scl 7 0 0 (1, 0, [4]), scl 8 0 0 (1, 0, [-4]), rt [-15] [(0, 4)]]
  | 1, 4, 1 => [scl 7 0 0 (1, 0, [4]), scl 8 0 0 (1, 0, [-4]), rt [-16] [(0, 4)]]
  | 1, 4, 2 => [scl 7 0 0 (1, 0, [4]), scl 8 0 0 (1, 0, [-4]), rt [-16, -64, -96, -64, -17] [(-1, 4), (0, 4)]]
  | 1, 4, 3 => [scl 7 0 0 (1, 0, [4]), scl 8 0 0 (1, 0, [-4]), rt [-256, -1536, -3968, -5760, -5153, -2916, -1022, -204, -18] [(-2, 4), (-1, 4), (0, 4)]]
  | 1, 4, 4 => [scl 7 0 0 (1, 0, [4]), scl 8 0 0 (1, 0, [-4]), rt [-20736, -152064, -501120, -981120, -1272449, -1153008, -749308, -352320, -119144, -28320, -4504, -432, -19] [(-3, 4), (-2, 4), (-1, 4), (0, 4)]]
  | 2, 2, 0 => [(1, 0, [2]), rt [-1] [(0, 2)]]
  | 2, 2, 1 => [(1, 0, [2]), rt [0] [(0, 2)]]
  | 2, 2, 2 => [(1, 0, [2]), rt [0, 0, 1] [(-1, 2), (0, 2)]]
  | 2, 2, 3 => [(1, 0, [2]), rt [0, 0, 5, 6, 2] [(-2, 2), (-1, 2), (0, 2)]]
  | 2, 2, 4 => [(1, 0, [2]), rt [0, 0, 49, 96, 72, 24, 3] [(-3, 2), (-2, 2), (-1, 2), (0, 2)]]
  | 2, 3, 0 => [(1, 0, [3]), rt [-1] [(0, 3)]]
  | 2, 3, 1 => [(1, 0, [3]), rt [0] [(0, 3)]]
  | 2, 3, 2 => [(1, 0, [3]), rt [0, 0, 0, 1] [(-1, 3), (0, 3)]]
  | 2, 3, 3 => [(1, 0, [3]), rt [0, 0, 0, 9, 15, 9, 2] [(-2, 3), (-1, 3), (0, 3)]]
  | 2, 3, 4 => [(1, 0, [3]), rt [0, 0, 0, 251, 684, 795, 504, 183, 36, 3] [(-3, 3), (-2, 3), (-1, 3), (0, 3)]]
  | 2, 4, 0 => [(1, 0, [4]), rt [-1] [(0, 4)]]
  | 2, 4, 1 => [(1, 0, [4]), rt [0] [(0, 4)]]
  | 2, 4, 2 => [(1, 0, [4]), rt [0, 0, 0, 0, 1] [(-1, 4), (0, 4)]]
  | 2, 4, 3 => [(1, 0, [4]), rt [0, 0, 0, 0, 17, 36, 30, 12, 2] [(-2, 4), (-1, 4), (0, 4)]]
  | 2, 4, 4 => [(1, 0, [4]), rt [0, 0, 0, 0, 1393, 4848, 7484, 6720, 3848, 1440, 344, 48, 3] [(-3, 4), (-2, 4), (-1, 4), (0, 4)]]
  | _, _, _ => []

/-- the prefactors `-(1/(n(2n-1))) e^{(4)}_{j,i,4-i}(n) ρ_j(n)` at `N = 2n` -/
noncomputable def fz : ℕ → ℕ → S03.DTerm
  | 0, 0 => rt [0, 0, 0, 0, -1 / 192] [(0, 1), (1, 1)]
  | 0, 1 => rt [0, 0, 0, -1 / 24, 1 / 48] [(0, 1), (1, 1)]
  | 0, 2 => rt [0, 0, -1 / 8, 1 / 8, -1 / 32] [(0, 1), (1, 1)]
  | 0, 3 => rt [0, -1 / 6, 1 / 4, -1 / 8, 1 / 48] [(0, 1), (1, 1)]
  | 0, 4 => rt [-1 / 12, 1 / 6, -1 / 8, 1 / 24, -1 / 192] [(0, 1), (1, 1)]
  | 1, 0 => rt [0, -1 / 24, -1 / 24, 1 / 48, 1 / 24, 1 / 48] [(0, 2), (1, 1)]
  | 1, 1 => rt [0, 0, 1 / 12, 1 / 24, -1 / 24, -1 / 12] [(0, 2), (1, 1)]
  | 1, 2 => rt [0, 0, 0, 0, -1 / 8, 1 / 8] [(0, 2), (1, 1)]
  | 1, 3 => rt [0, 0, 1 / 12, -5 / 24, 5 / 24, -1 / 12] [(0, 2), (1, 1)]
  | 1, 4 => rt [0, 1 / 24, -1 / 8, 7 / 48, -1 / 12, 1 / 48] [(0, 2), (1, 1)]
  | 2, 0 => rt [0, 0, 7 / 24, 5 / 48, -5 / 24, -5 / 32, -1 / 32] [(-1, 1), (0, 2), (1, 1)]
  | 2, 1 => rt [0, 0, -2 / 3, -5 / 12, 11 / 24, 1 / 2, 1 / 8] [(-1, 1), (0, 2), (1, 1)]
  | 2, 2 => rt [0, 0, 1 / 2, 1 / 2, -1 / 4, -9 / 16, -3 / 16] [(-1, 1), (0, 2), (1, 1)]
  | 2, 3 => rt [0, 0, -1 / 6, -1 / 6, -1 / 24, 1 / 4, 1 / 8] [(-1, 1), (0, 2), (1, 1)]
  | 2, 4 => rt [0, 0, 1 / 24, -1 / 48, 1 / 24, -1 / 32, -1 / 32] [(-1, 1), (0, 2), (1, 1)]
  | 3, 0 => rt [0, 0, -1 / 2, -13 / 24, 5 / 16, 25 / 48, 3 / 16, 1 / 48] [(-2, 1), (-1, 1), (0, 2), (1, 1)]
  | 3, 1 => rt [0, 0, 5 / 3, 23 / 12, -23 / 24, -11 / 6, -17 / 24, -1 / 12] [(-2, 1), (-1, 1), (0, 2), (1, 1)]
  | 3, 2 => rt [0, 0, -2, -5 / 2, 1, 19 / 8, 1, 1 / 8] [(-2, 1), (-1, 1), (0, 2), (1, 1)]
  | 3, 3 => rt [0, 0, 1, 17 / 12, -3 / 8, -4 / 3, -5 / 8, -1 / 12] [(-2, 1), (-1, 1), (0, 2), (1, 1)]
  | 3, 4 => rt [0, 0, -1 / 6, -7 / 24, 1 / 48, 13 / 48, 7 / 48, 1 / 48] [(-2, 1), (-1, 1), (0, 2), (1, 1)]
  | 4, 0 => rt [0, 0, 1 / 2, 17 / 24, -7 / 48, -61 / 96, -67 / 192, -7 / 96, -1 / 192] [(-3, 1), (-2, 1), (-1, 1), (0, 2), (1, 1)]
  | 4, 1 => rt [0, 0, -2, -17 / 6, 7 / 12, 61 / 24, 67 / 48, 7 / 24, 1 / 48] [(-3, 1), (-2, 1), (-1, 1), (0, 2), (1, 1)]
  | 4, 2 => rt [0, 0, 3, 17 / 4, -7 / 8, -61 / 16, -67 / 32, -7 / 16, -1 / 32] [(-3, 1), (-2, 1), (-1, 1), (0, 2), (1, 1)]
  | 4, 3 => rt [0, 0, -2, -17 / 6, 7 / 12, 61 / 24, 67 / 48, 7 / 24, 1 / 48] [(-3, 1), (-2, 1), (-1, 1), (0, 2), (1, 1)]
  | 4, 4 => rt [0, 0, 1 / 2, 17 / 24, -7 / 48, -61 / 96, -67 / 192, -7 / 96, -1 / 192] [(-3, 1), (-2, 1), (-1, 1), (0, 2), (1, 1)]
  | _, _ => (0, 0, [])


lemma atmA_ST (r j n : ℕ) (hr : 1 ≤ r) (hr4 : r ≤ 4) (hj : j ≤ 4) (hn : 1 ≤ n) :
    S03.ST (atm 0 r j) (2 * n) = V05.HrQ r (n + j - 1) - V05.HrQ r (n + j + n - 1) := by
  rw [Hn1 r n hr hn, Hn3 r n hr hn]
  have hn' : (1 : ℚ) ≤ n := by exact_mod_cast hn
  have h0 : (n : ℚ) ≠ 0 := by positivity
  interval_cases r <;> interval_cases j <;>
    (simp [atm, ST_cons, ST_nil, rt_eval, horner, ratEval_scl, ratEval_one, Finset.sum_range_succ,
      scl_w, rt_w, Z_nil] <;> (try field_simp) <;> ring)

end P04

namespace P04

lemma atmB_ST (r j n : ℕ) (hr : 1 ≤ r) (hr4 : r ≤ 4) (hj : j ≤ 4) (hn : 1 ≤ n) :
    S03.ST (atm 1 r j) (2 * n) = V05.HrQ r (n - 1) - V05.HrQ r (n + j + n - 1) := by
  rw [Hn2 r n hr hn, Hn3 r n hr hn]
  have hn' : (1 : ℚ) ≤ n := by exact_mod_cast hn
  have h0 : (n : ℚ) ≠ 0 := by positivity
  interval_cases r <;> interval_cases j <;>
    (simp [atm, ST_cons, ST_nil, rt_eval, horner, ratEval_scl, ratEval_one, Finset.sum_range_succ,
      scl_w, rt_w, Z_nil] <;> (try field_simp) <;> ring)

lemma atmG_ST (r j n : ℕ) (hr : 2 ≤ r) (hr4 : r ≤ 4) (hj : j ≤ 4) (hn : 1 ≤ n) :
    S03.ST (atm 2 r j) (2 * n) = V05.HrQ r (n + j + n - 1) := by
  rw [Hn3 r n (by omega) hn]
  have hn' : (1 : ℚ) ≤ n := by exact_mod_cast hn
  have h0 : (n : ℚ) ≠ 0 := by positivity
  interval_cases r <;> interval_cases j <;>
    (simp [atm, ST_cons, ST_nil, rt_eval, horner, ratEval_scl, ratEval_one, Finset.sum_range_succ,
      scl_w, rt_w, Z_nil] <;> (try field_simp) <;> ring)

lemma atm_idx (k r j : ℕ) : ∀ τ ∈ atm k r j, S03.IsIndex τ.2.2 := by
  unfold atm
  split <;> intro τ hτ <;> simp only [List.mem_cons, List.mem_nil_iff, or_false] at hτ <;>
    rcases hτ with rfl | rfl | rfl <;> simp [S03.IsIndex, scl_w, rt_w]

lemma fz_w (j i : ℕ) : (fz j i).2.2 = [] := by
  unfold fz; split <;> simp [rt_w]

lemma fz_eval (j i n : ℕ) (hj : j ≤ 4) (hi : i ≤ 4) (hn : 1 ≤ n) :
    S03.ratEval (fz j i).1 (fz j i).2.1 ((2 * n : ℕ) : ℚ) =
      -(1 / ((n : ℚ) * (2 * (n : ℚ) - 1))) * (S04.cf4 n j i * V05.rhoQ j n) := by
  have hn' : (1 : ℚ) ≤ n := by exact_mod_cast hn
  have h0 : (n : ℚ) ≠ 0 := by positivity
  have h1 : 2 * (n : ℚ) - 1 ≠ 0 := by linarith
  have h1' : (n : ℚ) * 2 - 1 ≠ 0 := by linarith
  interval_cases j <;> interval_cases i <;>
    (simp [fz, rt_eval, horner, S04.cf4, V05.rhoQ, Finset.prod_range_succ] <;> (try field_simp) <;> ring)

end P04

namespace P04
section gHsec
open MvPolynomial

lemma coeff_XX (a b i m : ℕ) : ((X 0 : Q04.Q2) ^ a * X 1 ^ b).coeff
    (Finsupp.single 0 i + Finsupp.single 1 m) = if a = i ∧ b = m then 1 else 0 := by
  rw [Q04.mon_eq, coeff_monomial]
  by_cases h : a = i ∧ b = m
  · rw [if_pos h, if_pos (by rw [h.1, h.2])]
  · rw [if_neg h, if_neg]; intro h'; apply h
    have h0 := congrArg (fun f => f 0) h'; have h1 := congrArg (fun f => f 1) h'
    simp at h0 h1; exact ⟨h0, h1⟩

lemma coeff_sXX (c : ℚ) (a b i m : ℕ) : (c • ((X 0 : Q04.Q2) ^ a * X 1 ^ b)).coeff
    (Finsupp.single 0 i + Finsupp.single 1 m) = if a = i ∧ b = m then c else 0 := by
  rw [coeff_smul, coeff_XX, smul_eq_mul]; split_ifs <;> simp

lemma gH_two (p q : ℕ) :
    Q04.gH p q 2 0 = (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) ^ 2 / 2
      - (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) / 2 ∧
    Q04.gH p q 1 1 = (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) *
      (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + V05.HrQ 2 (p + q - 1) ∧
    Q04.gH p q 0 2 = (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) ^ 2 / 2
      - (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) / 2 := by
  set T1 : Q04.Q2 := ((-1 : ℚ) ^ (1 - 1) / (1 : ℕ)) •
      (V05.HrQ 1 (p - 1) • (X 0 : Q04.Q2) ^ 1 + V05.HrQ 1 (q - 1) • X 1 ^ 1
        - V05.HrQ 1 (p + q - 1) • (X 0 + X 1) ^ 1) with hT1
  set T2 : Q04.Q2 := ((-1 : ℚ) ^ (2 - 1) / (2 : ℕ)) •
      (V05.HrQ 2 (p - 1) • (X 0 : Q04.Q2) ^ 2 + V05.HrQ 2 (q - 1) • X 1 ^ 2
        - V05.HrQ 2 (p + q - 1) • (X 0 + X 1) ^ 2) with hT2
  have o1 : Q04.Ord 1 T1 := Q04.Hb_term_Ord 1 p q
  have o2 : Q04.Ord 2 T2 := Q04.Hb_term_Ord 2 p q
  have hHb : Q04.Hb 2 p q = T1 + T2 := by
    unfold Q04.Hb
    rw [show Finset.Icc 1 2 = {1, 2} by decide, Finset.sum_insert (by decide), Finset.sum_singleton]
  have hE : Q04.expT 2 (T1 + T2) = (1 + T1 + T2 + (1/2 : ℚ) • T1 ^ 2) +
      (T1 * T2 + (1/2 : ℚ) • T2 ^ 2) := by
    unfold Q04.expT
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial]
    norm_num
    algebra
  have hR : Q04.Ord 3 (T1 * T2 + (1/2 : ℚ) • T2 ^ 2) :=
    Q04.Ord_add (Q04.Ord_mul o1 o2) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_pow o2 2) (by norm_num)))
  have hG : 1 + T1 + T2 + (1/2 : ℚ) • T1 ^ 2 =
      (1 : ℚ) • ((X 0 : Q04.Q2) ^ 0 * X 1 ^ 0)
      + (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) • ((X 0 : Q04.Q2) ^ 1 * X 1 ^ 0)
      + (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) • ((X 0 : Q04.Q2) ^ 0 * X 1 ^ 1)
      + ((V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) ^ 2 / 2
          - (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) / 2) • ((X 0 : Q04.Q2) ^ 2 * X 1 ^ 0)
      + ((V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) *
          (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + V05.HrQ 2 (p + q - 1)) •
          ((X 0 : Q04.Q2) ^ 1 * X 1 ^ 1)
      + ((V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) ^ 2 / 2
          - (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) / 2) • ((X 0 : Q04.Q2) ^ 0 * X 1 ^ 2) := by
    rw [hT1, hT2]
    norm_num
    algebra
  have key : ∀ i m : ℕ, i + m = 2 → Q04.gH p q i m =
      (1 + T1 + T2 + (1/2 : ℚ) • T1 ^ 2).coeff (Finsupp.single 0 i + Finsupp.single 1 m) := by
    intro i m him
    unfold Q04.gH
    rw [him, hHb, hE, _root_.MvPolynomial.coeff_add, Q04.coeff_eq_of_Ord (d := 2) (F := T1 * T2 + (1/2 : ℚ) • T2 ^ 2)
      (G := 0) (by simpa using hR) _ (by simp; omega)]
    simp
  refine ⟨?_, ?_, ?_⟩
  · rw [key 2 0 rfl, hG]; simp only [_root_.MvPolynomial.coeff_add, coeff_sXX]; norm_num
  · rw [key 1 1 rfl, hG]; simp only [_root_.MvPolynomial.coeff_add, coeff_sXX]; norm_num
  · rw [key 0 2 rfl, hG]; simp only [_root_.MvPolynomial.coeff_add, coeff_sXX]; norm_num

end gHsec

/-- the lists for `-(i!(4-i)!)·[𝔷2]`-free parts: `β₁² - β₂`, `α₁β₁ + γ₂`, `α₁² - α₂` -/
noncomputable def Dz2 (j : ℕ) : ℕ → List S03.DTerm
  | 1 => mulL (atm 1 1 j) (atm 1 1 j) ++ (atm 1 2 j).map (scl (-1) 0 0)
  | 2 => mulL (atm 0 1 j) (atm 1 1 j) ++ atm 2 2 j
  | 3 => mulL (atm 0 1 j) (atm 0 1 j) ++ (atm 0 2 j).map (scl (-1) 0 0)
  | _ => []

/-- the factors `-3, -4, -3` -/
def cz2 : ℕ → ℚ
  | 1 => -3
  | 2 => -4
  | 3 => -3
  | _ => 0

/-- the element `T_z` of `g8:lem:mu` for `𝔷 = 𝔷2`, `k = 4` -/
noncomputable def Tz2 : List S03.DTerm :=
  (List.range 5).flatMap fun j => (List.range 5).flatMap fun i =>
    mulL [fz j i] ((Dz2 j i).map (scl (cz2 i) 0 0))

lemma ST_single_nil (τ : S03.DTerm) (hτ : τ.2.2 = []) (L : List S03.DTerm)
    (hL : ∀ τ' ∈ L, S03.IsIndex τ'.2.2) (N : ℕ) :
    S03.ST (mulL [τ] L) N = S03.ratEval τ.1 τ.2.1 (N : ℚ) * S03.ST L N := by
  rw [ST_mulL _ _ (by intro τ' h; simp at h; subst h; rw [hτ]; intro a ha; simp at ha) hL,
    ST_cons, ST_nil, hτ, Z_nil]
  ring

lemma idx_append {L L' : List S03.DTerm} (h : ∀ τ ∈ L, S03.IsIndex τ.2.2)
    (h' : ∀ τ ∈ L', S03.IsIndex τ.2.2) : ∀ τ ∈ L ++ L', S03.IsIndex τ.2.2 := by
  intro τ hτ; rcases List.mem_append.mp hτ with h1 | h1; exact h τ h1; exact h' τ h1

lemma idx_map_scl {L : List S03.DTerm} (h : ∀ τ ∈ L, S03.IsIndex τ.2.2) (r : ℚ) (p : ℤ) (k : ℕ) :
    ∀ τ ∈ L.map (scl r p k), S03.IsIndex τ.2.2 := by
  intro τ hτ; obtain ⟨τ', h1, rfl⟩ := List.mem_map.mp hτ; rw [scl_w]; exact h τ' h1

lemma idx_mulL {L L' : List S03.DTerm} (h : ∀ τ ∈ L, S03.IsIndex τ.2.2)
    (h' : ∀ τ ∈ L', S03.IsIndex τ.2.2) : ∀ τ ∈ mulL L L', S03.IsIndex τ.2.2 := by
  intro τ hτ
  unfold mulL at hτ
  obtain ⟨τ1, h1, hτ⟩ := List.mem_flatMap.mp hτ
  obtain ⟨τ2, h2, hτ⟩ := List.mem_flatMap.mp hτ
  obtain ⟨w, hw, rfl⟩ := List.mem_map.mp hτ
  exact stuf_index _ _ (h τ1 h1) (h' τ2 h2) w hw

lemma Dz2_idx (j i : ℕ) : ∀ τ ∈ Dz2 j i, S03.IsIndex τ.2.2 := by
  unfold Dz2
  split
  · exact idx_append (idx_mulL (atm_idx _ _ _) (atm_idx _ _ _)) (idx_map_scl (atm_idx _ _ _) _ _ _)
  · exact idx_append (idx_mulL (atm_idx _ _ _) (atm_idx _ _ _)) (atm_idx _ _ _)
  · exact idx_append (idx_mulL (atm_idx _ _ _) (atm_idx _ _ _)) (idx_map_scl (atm_idx _ _ _) _ _ _)
  · intro τ h; simp at h

lemma Dz2_val (j i n : ℕ) (hj : j ≤ 4) (hi : i ≤ 4) (hn : 1 ≤ n) :
    cz2 i * S03.ST (Dz2 j i) (2 * n) =
      (Q04.DsharpF i (4 - i) (n + j) n).coeff (Finsupp.single 2 1) := by
  interval_cases i
  · rw [show (4 : ℕ) - 0 = 2 + 2 from rfl, Q04.X2D_edge0]; simp [Dz2, S03.ST]
  · have h := Q04.X2D_mid 2 0 (n + j) n (by norm_num)
    rw [show (4 : ℕ) - 1 = 2 + 1 - 0 from rfl, h]
    rw [(gH_two (n + j) n).2.2]
    simp only [Dz2, cz2, ST_append, ST_map_scl, ST_mulL _ _ (atm_idx _ _ _) (atm_idx _ _ _),
      atmB_ST 1 j n le_rfl (by norm_num) hj hn, atmB_ST 2 j n (by norm_num) (by norm_num) hj hn]
    rw [show n + j + n - 1 = n + j + n - 1 from rfl]
    norm_num <;> ring
  · have h := Q04.X2D_mid 2 1 (n + j) n (by norm_num)
    rw [show (4 : ℕ) - 2 = 2 + 1 - 1 from rfl, show (2 : ℕ) = 1 + 1 from rfl, h]
    rw [(gH_two (n + j) n).2.1]
    simp only [Dz2, cz2, ST_append, ST_map_scl, ST_mulL _ _ (atm_idx _ _ _) (atm_idx _ _ _),
      atmA_ST 1 j n le_rfl (by norm_num) hj hn, atmB_ST 1 j n le_rfl (by norm_num) hj hn,
      atmG_ST 2 j n le_rfl (by norm_num) hj hn]
    norm_num <;> ring
  · have h := Q04.X2D_mid 2 2 (n + j) n (by norm_num)
    rw [show (4 : ℕ) - 3 = 2 + 1 - 2 from rfl, show (3 : ℕ) = 2 + 1 from rfl, h]
    rw [(gH_two (n + j) n).1]
    simp only [Dz2, cz2, ST_append, ST_map_scl, ST_mulL _ _ (atm_idx _ _ _) (atm_idx _ _ _),
      atmA_ST 1 j n le_rfl (by norm_num) hj hn, atmA_ST 2 j n (by norm_num) (by norm_num) hj hn]
    norm_num <;> ring
  · rw [show (4 : ℕ) - 4 = 0 from rfl, show (4 : ℕ) = 2 + 2 from rfl, Q04.X2D_edge1]; simp [Dz2, S03.ST]

lemma Tz2_ok (n : ℕ) (hn : 1 ≤ n) :
    (Q04.Pf 4 n).coeff (Finsupp.single 2 1) = S03.ST Tz2 (2 * n) := by
  unfold Q04.Pf Tz2
  rw [MvPolynomial.coeff_smul, smul_eq_mul, ST_flatMap, list_range_sum, MvPolynomial.coeff_sum,
    Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro j hj
  rw [ST_flatMap, list_range_sum, MvPolynomial.coeff_sum, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i hi
  have hj' : j ≤ 4 := by simp at hj; omega
  have hi' : i ≤ 4 := by simp at hi; omega
  rw [ST_single_nil _ (fz_w j i) _ (idx_map_scl (Dz2_idx j i) _ _ _), fz_eval j i n hj' hi' hn,
    ST_map_scl, MvPolynomial.coeff_smul, smul_eq_mul, ← Dz2_val j i n hj' hi' hn,
    S04.eF4 n j i hj' hi']
  simp
  ring

end P04

namespace P04

lemma Dval_zero_of (L : List S03.DTerm) (u : List ℤ) (N : ℚ) (h : ∀ τ ∈ L, τ.2.2 ≠ u) :
    Dval L u N = 0 := by
  unfold Dval
  rw [List.filter_eq_nil_iff.mpr (fun τ hτ => by simpa using h τ hτ)]
  simp

lemma Dval_mulL_fz (j i : ℕ) (L : List S03.DTerm) (u : List ℤ) (N : ℚ) :
    Dval (mulL [fz j i] L) u N = S03.ratEval (fz j i).1 (fz j i).2.1 N * Dval L u N := by
  unfold mulL
  rw [List.flatMap_singleton, Dval_flatMap, fz_w]
  induction L with
  | nil => simp [Dval]
  | cons τ L ih =>
    rw [List.map_cons, List.sum_cons, ih, Dval_cons, stuf_nil_left, List.map_singleton, Dval_cons,
      Dval_nil]
    simp only [S08.PA.ratEval_mul]
    split_ifs <;> ring

lemma sgn_pos_lit : Int.sign 1 = 1 ∧ Int.sign 2 = 1 ∧ Int.sign 3 = 1 ∧ Int.sign 4 = 1 ∧
    Int.sign (-1) = -1 ∧ Int.sign (-2) = -1 ∧ Int.sign (-3) = -1 ∧ Int.sign (-4) = -1 := by decide

lemma range5 : List.range 5 = [0, 1, 2, 3, 4] := rfl

lemma opl_eval (a b : ℤ) : opl a b = Int.sign a * Int.sign b * ((a.natAbs + b.natAbs : ℕ) : ℤ) := rfl

end P04

namespace P04

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 20000 in
lemma P2w0_0 (N : ℚ) (hN : 1 < N) : Dval (Lι Psi2T 0) [] N = Dval Tz2 [] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, Psi2T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz2, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [Dz2, cz2, Dval_map_scl, Dval_append, Dval_mulL, Dval_nil, List.map_nil]
  simp only [atm, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, rt_w]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, rt_w]
  field_simp
  ring

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 20000 in
lemma P2w1_0 (N : ℚ) (hN : 1 < N) : Dval (Lι Psi2T 1) [] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, Psi2T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, rt_w]
  field_simp
  ring

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 20000 in
lemma P2w0_1 (N : ℚ) (hN : 1 < N) : Dval (Lι Psi2T 0) [-2] N = Dval Tz2 [-2] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, Psi2T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz2, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [Dz2, cz2, Dval_map_scl, Dval_append, Dval_mulL, Dval_nil, List.map_nil]
  simp only [atm, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, rt_w]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, rt_w]
  field_simp
  ring

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 20000 in
lemma P2w1_1 (N : ℚ) (hN : 1 < N) : Dval (Lι Psi2T 1) [-2] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, Psi2T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, rt_w]
  field_simp
  ring

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 20000 in
lemma P2w0_2 (N : ℚ) (hN : 1 < N) : Dval (Lι Psi2T 0) [-1] N = Dval Tz2 [-1] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, Psi2T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz2, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [Dz2, cz2, Dval_map_scl, Dval_append, Dval_mulL, Dval_nil, List.map_nil]
  simp only [atm, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, rt_w]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, rt_w]
  field_simp
  ring

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 20000 in
lemma P2w1_2 (N : ℚ) (hN : 1 < N) : Dval (Lι Psi2T 1) [-1] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, Psi2T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, rt_w]
  field_simp
  ring

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 20000 in
lemma P2w0_3 (N : ℚ) (hN : 1 < N) : Dval (Lι Psi2T 0) [-1, -1] N = Dval Tz2 [-1, -1] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, Psi2T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz2, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [Dz2, cz2, Dval_map_scl, Dval_append, Dval_mulL, Dval_nil, List.map_nil]
  simp only [atm, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, rt_w]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, rt_w]
  field_simp
  ring

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 20000 in
lemma P2w1_3 (N : ℚ) (hN : 1 < N) : Dval (Lι Psi2T 1) [-1, -1] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, Psi2T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, rt_w]
  field_simp
  ring

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 20000 in
lemma P2w0_4 (N : ℚ) (hN : 1 < N) : Dval (Lι Psi2T 0) [2] N = Dval Tz2 [2] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, Psi2T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz2, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [Dz2, cz2, Dval_map_scl, Dval_append, Dval_mulL, Dval_nil, List.map_nil]
  simp only [atm, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, rt_w]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, rt_w]
  field_simp
  ring

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 20000 in
lemma P2w1_4 (N : ℚ) (hN : 1 < N) : Dval (Lι Psi2T 1) [2] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, Psi2T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, rt_w]
  field_simp
  ring

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 20000 in
lemma P2o_L0 (u : List ℤ) (N : ℚ) (h0 : ¬ ([] = u)) (h1 : ¬ ([-2] = u)) (h2 : ¬ ([-1] = u)) (h3 : ¬ ([-1, -1] = u)) (h4 : ¬ ([2] = u)) : Dval (Lι Psi2T 0) u N = 0 := by
  simp only [Lι, Psi2T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, rt_w, h0, h1, h2, h3, h4]
  try simp

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 20000 in
lemma P2o_L1 (u : List ℤ) (N : ℚ) (h0 : ¬ ([] = u)) (h1 : ¬ ([-2] = u)) (h2 : ¬ ([-1] = u)) (h3 : ¬ ([-1, -1] = u)) (h4 : ¬ ([2] = u)) : Dval (Lι Psi2T 1) u N = 0 := by
  simp only [Lι, Psi2T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, rt_w, h0, h1, h2, h3, h4]
  try simp

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 20000 in
lemma P2o_T (u : List ℤ) (N : ℚ) (h0 : ¬ ([] = u)) (h1 : ¬ ([-2] = u)) (h2 : ¬ ([-1] = u)) (h3 : ¬ ([-1, -1] = u)) (h4 : ¬ ([2] = u)) : Dval Tz2 u N = 0 := by
  simp only [Tz2, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [Dz2, cz2, Dval_map_scl, Dval_append, Dval_mulL, Dval_nil, List.map_nil]
  simp only [atm, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, rt_w]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, rt_w, h0, h1, h2, h3, h4]
  try simp

lemma Deq_P2_0 : Deq (Lι Psi2T 0) Tz2 := by
  intro u N hN
  by_cases e0 : [] = u
  · subst e0; exact P2w0_0 N hN
  by_cases e1 : [-2] = u
  · subst e1; exact P2w0_1 N hN
  by_cases e2 : [-1] = u
  · subst e2; exact P2w0_2 N hN
  by_cases e3 : [-1, -1] = u
  · subst e3; exact P2w0_3 N hN
  by_cases e4 : [2] = u
  · subst e4; exact P2w0_4 N hN
  rw [P2o_L0 u N e0 e1 e2 e3 e4, P2o_T u N e0 e1 e2 e3 e4]

lemma Deq_P2_1 : Deq (Lι Psi2T 1) [] := by
  intro u N hN
  by_cases e0 : [] = u
  · subst e0; rw [P2w1_0 N hN, Dval_nil]
  by_cases e1 : [-2] = u
  · subst e1; rw [P2w1_1 N hN, Dval_nil]
  by_cases e2 : [-1] = u
  · subst e2; rw [P2w1_2 N hN, Dval_nil]
  by_cases e3 : [-1, -1] = u
  · subst e3; rw [P2w1_3 N hN, Dval_nil]
  by_cases e4 : [2] = u
  · subst e4; rw [P2w1_4 N hN, Dval_nil]
  rw [P2o_L1 u N e0 e1 e2 e3 e4, Dval_nil]

end P04

namespace P04

lemma Psi2C_ok : ∀ t ∈ Psi2C, t.2 ≤ 1 := by simp [Psi2C]

lemma Psi2_even (n : ℕ) (hn : 1 ≤ n) :
    coefL Psi2T Psi2C (2 * (n : ℤ)) = (Q04.Pf 4 n).coeff (Finsupp.single 2 1) := by
  rw [show 2 * (n : ℤ) = 2 * (n : ℤ) + 0 by ring, coefL_eq Psi2T Psi2C Psi2T_ok Psi2C_ok 0 le_rfl n hn,
    ST_of_Deq _ _ Deq_P2_0 _ (by omega), ← Tz2_ok n hn]

lemma Psi2_odd (n : ℕ) (hn : 1 ≤ n) : coefL Psi2T Psi2C (2 * (n : ℤ) + 1) = 0 := by
  rw [coefL_eq Psi2T Psi2C Psi2T_ok Psi2C_ok 1 (by norm_num) n hn,
    ST_of_Deq _ _ Deq_P2_1 _ (by omega), ST_nil]

end P04

/-- M8 for `Ψ_2`: (i), (ii) from `coefL_eq`, the element `Tz2` (`Tz2_ok`) and the identities
`Deq_P2_0`, `Deq_P2_1`; (iii) from `Psi2_low`. -/
theorem P04.proof_4 : P04.stmt_4 :=
  ⟨P04.Psi2_even, P04.Psi2_odd, P04.Psi2_low⟩

theorem P04.proof_5 : P04.stmt_5 := fun c hc0 hc1 =>
  P04.series_of_coef P04.Psi2T P04.Psi2C P04.Psi2T_ok _ P04.proof_4.1 P04.proof_4.2.1
    P04.proof_4.2.2 c hc0 hc1


namespace P04
section gH4sec
open MvPolynomial

set_option maxHeartbeats 10000000 in
lemma gH_four (p q : ℕ) (i : ℕ) (hi : i ≤ 4) :
    Q04.gH p q i (4 - i) = (((1 : ℚ)) • ((X 0 : Q04.Q2) ^ 0 * X 1 ^ 0) + ((1 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1))) • ((X 0 : Q04.Q2) ^ 0 * X 1 ^ 1) + ((1 / 2 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 2 : ℚ) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1))) • ((X 0 : Q04.Q2) ^ 0 * X 1 ^ 2) + ((1 / 6 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 2 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) + (1 / 3 : ℚ) * (V05.HrQ 3 (q - 1) - V05.HrQ 3 (p + q - 1))) • ((X 0 : Q04.Q2) ^ 0 * X 1 ^ 3) + ((1 / 24 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 4 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) + (1 / 3 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 3 (q - 1) - V05.HrQ 3 (p + q - 1)) + (1 / 8 : ℚ) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) + (-1 / 4 : ℚ) * (V05.HrQ 4 (q - 1) - V05.HrQ 4 (p + q - 1))) • ((X 0 : Q04.Q2) ^ 0 * X 1 ^ 4) + ((1 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1))) • ((X 0 : Q04.Q2) ^ 1 * X 1 ^ 0) + ((1 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (1 : ℚ) * V05.HrQ 2 (p + q - 1)) • ((X 0 : Q04.Q2) ^ 1 * X 1 ^ 1) + ((1 / 2 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 2 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) + (1 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 2 (p + q - 1) + (-1 : ℚ) * V05.HrQ 3 (p + q - 1)) • ((X 0 : Q04.Q2) ^ 1 * X 1 ^ 2) + ((1 / 6 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 2 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) + (1 / 3 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 3 (q - 1) - V05.HrQ 3 (p + q - 1)) + (1 / 2 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 2 (p + q - 1) + (-1 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 3 (p + q - 1) + (-1 / 2 : ℚ) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) * V05.HrQ 2 (p + q - 1) + (1 : ℚ) * V05.HrQ 4 (p + q - 1)) • ((X 0 : Q04.Q2) ^ 1 * X 1 ^ 3) + ((1 / 2 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 2 : ℚ) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1))) • ((X 0 : Q04.Q2) ^ 2 * X 1 ^ 0) + ((1 / 2 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (1 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 2 (p + q - 1) + (-1 / 2 : ℚ) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (-1 : ℚ) * V05.HrQ 3 (p + q - 1)) • ((X 0 : Q04.Q2) ^ 2 * X 1 ^ 1) + ((1 / 4 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 4 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) + (1 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 2 (p + q - 1) + (-1 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 3 (p + q - 1) + (-1 / 4 : ℚ) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (1 / 4 : ℚ) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) + (-1 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 3 (p + q - 1) + (1 / 2 : ℚ) * V05.HrQ 2 (p + q - 1) * V05.HrQ 2 (p + q - 1) + (3 / 2 : ℚ) * V05.HrQ 4 (p + q - 1)) • ((X 0 : Q04.Q2) ^ 2 * X 1 ^ 2) + ((1 / 6 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 2 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) + (1 / 3 : ℚ) * (V05.HrQ 3 (p - 1) - V05.HrQ 3 (p + q - 1))) • ((X 0 : Q04.Q2) ^ 3 * X 1 ^ 0) + ((1 / 6 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (1 / 2 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 2 (p + q - 1) + (-1 / 2 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (-1 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 3 (p + q - 1) + (-1 / 2 : ℚ) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) * V05.HrQ 2 (p + q - 1) + (1 / 3 : ℚ) * (V05.HrQ 3 (p - 1) - V05.HrQ 3 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (1 : ℚ) * V05.HrQ 4 (p + q - 1)) • ((X 0 : Q04.Q2) ^ 3 * X 1 ^ 1) + ((1 / 24 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 4 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) + (1 / 3 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 3 (p - 1) - V05.HrQ 3 (p + q - 1)) + (1 / 8 : ℚ) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) + (-1 / 4 : ℚ) * (V05.HrQ 4 (p - 1) - V05.HrQ 4 (p + q - 1))) • ((X 0 : Q04.Q2) ^ 4 * X 1 ^ 0)).coeff (Finsupp.single 0 i + Finsupp.single 1 (4 - i)) := by
  set T1 : Q04.Q2 := ((-1 : ℚ) ^ (1 - 1) / (1 : ℕ)) •
      (V05.HrQ 1 (p - 1) • (X 0 : Q04.Q2) ^ 1 + V05.HrQ 1 (q - 1) • X 1 ^ 1
        - V05.HrQ 1 (p + q - 1) • (X 0 + X 1) ^ 1) with hT1
  have o1 : Q04.Ord 1 T1 := Q04.Hb_term_Ord 1 p q
  set T2 : Q04.Q2 := ((-1 : ℚ) ^ (2 - 1) / (2 : ℕ)) •
      (V05.HrQ 2 (p - 1) • (X 0 : Q04.Q2) ^ 2 + V05.HrQ 2 (q - 1) • X 1 ^ 2
        - V05.HrQ 2 (p + q - 1) • (X 0 + X 1) ^ 2) with hT2
  have o2 : Q04.Ord 2 T2 := Q04.Hb_term_Ord 2 p q
  set T3 : Q04.Q2 := ((-1 : ℚ) ^ (3 - 1) / (3 : ℕ)) •
      (V05.HrQ 3 (p - 1) • (X 0 : Q04.Q2) ^ 3 + V05.HrQ 3 (q - 1) • X 1 ^ 3
        - V05.HrQ 3 (p + q - 1) • (X 0 + X 1) ^ 3) with hT3
  have o3 : Q04.Ord 3 T3 := Q04.Hb_term_Ord 3 p q
  set T4 : Q04.Q2 := ((-1 : ℚ) ^ (4 - 1) / (4 : ℕ)) •
      (V05.HrQ 4 (p - 1) • (X 0 : Q04.Q2) ^ 4 + V05.HrQ 4 (q - 1) • X 1 ^ 4
        - V05.HrQ 4 (p + q - 1) • (X 0 + X 1) ^ 4) with hT4
  have o4 : Q04.Ord 4 T4 := Q04.Hb_term_Ord 4 p q
  have hHb : Q04.Hb 4 p q = T1 + T2 + T3 + T4 := by
    unfold Q04.Hb
    rw [show Finset.Icc 1 4 = {1, 2, 3, 4} by decide, Finset.sum_insert (by decide),
      Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_singleton]
    abel
  obtain ⟨R, hR, hE⟩ : ∃ R : Q04.Q2, Q04.Ord 5 R ∧
      Q04.expT 4 (T1 + T2 + T3 + T4) = ((1 : ℚ) • (1) + (1 : ℚ) • (T4) + (1 : ℚ) • (T3) + (1 : ℚ) • (T2) + (1 / 2 : ℚ) • (T2 ^ 2) + (1 : ℚ) • (T1) + (1 : ℚ) • (T1 * T3) + (1 : ℚ) • (T1 * T2) + (1 / 2 : ℚ) • (T1 ^ 2) + (1 / 2 : ℚ) • (T1 ^ 2 * T2) + (1 / 6 : ℚ) • (T1 ^ 3) + (1 / 24 : ℚ) • (T1 ^ 4)) + R := by
    refine ⟨(1 / 2 : ℚ) • (T4 ^ 2) + (1 / 6 : ℚ) • (T4 ^ 3) + (1 / 24 : ℚ) • (T4 ^ 4) + (1 : ℚ) • (T3 * T4) + (1 / 2 : ℚ) • (T3 * T4 ^ 2) + (1 / 6 : ℚ) • (T3 * T4 ^ 3) + (1 / 2 : ℚ) • (T3 ^ 2) + (1 / 2 : ℚ) • (T3 ^ 2 * T4) + (1 / 4 : ℚ) • (T3 ^ 2 * T4 ^ 2) + (1 / 6 : ℚ) • (T3 ^ 3) + (1 / 6 : ℚ) • (T3 ^ 3 * T4) + (1 / 24 : ℚ) • (T3 ^ 4) + (1 : ℚ) • (T2 * T4) + (1 / 2 : ℚ) • (T2 * T4 ^ 2) + (1 / 6 : ℚ) • (T2 * T4 ^ 3) + (1 : ℚ) • (T2 * T3) + (1 : ℚ) • (T2 * T3 * T4) + (1 / 2 : ℚ) • (T2 * T3 * T4 ^ 2) + (1 / 2 : ℚ) • (T2 * T3 ^ 2) + (1 / 2 : ℚ) • (T2 * T3 ^ 2 * T4) + (1 / 6 : ℚ) • (T2 * T3 ^ 3) + (1 / 2 : ℚ) • (T2 ^ 2 * T4) + (1 / 4 : ℚ) • (T2 ^ 2 * T4 ^ 2) + (1 / 2 : ℚ) • (T2 ^ 2 * T3) + (1 / 2 : ℚ) • (T2 ^ 2 * T3 * T4) + (1 / 4 : ℚ) • (T2 ^ 2 * T3 ^ 2) + (1 / 6 : ℚ) • (T2 ^ 3) + (1 / 6 : ℚ) • (T2 ^ 3 * T4) + (1 / 6 : ℚ) • (T2 ^ 3 * T3) + (1 / 24 : ℚ) • (T2 ^ 4) + (1 : ℚ) • (T1 * T4) + (1 / 2 : ℚ) • (T1 * T4 ^ 2) + (1 / 6 : ℚ) • (T1 * T4 ^ 3) + (1 : ℚ) • (T1 * T3 * T4) + (1 / 2 : ℚ) • (T1 * T3 * T4 ^ 2) + (1 / 2 : ℚ) • (T1 * T3 ^ 2) + (1 / 2 : ℚ) • (T1 * T3 ^ 2 * T4) + (1 / 6 : ℚ) • (T1 * T3 ^ 3) + (1 : ℚ) • (T1 * T2 * T4) + (1 / 2 : ℚ) • (T1 * T2 * T4 ^ 2) + (1 : ℚ) • (T1 * T2 * T3) + (1 : ℚ) • (T1 * T2 * T3 * T4) + (1 / 2 : ℚ) • (T1 * T2 * T3 ^ 2) + (1 / 2 : ℚ) • (T1 * T2 ^ 2) + (1 / 2 : ℚ) • (T1 * T2 ^ 2 * T4) + (1 / 2 : ℚ) • (T1 * T2 ^ 2 * T3) + (1 / 6 : ℚ) • (T1 * T2 ^ 3) + (1 / 2 : ℚ) • (T1 ^ 2 * T4) + (1 / 4 : ℚ) • (T1 ^ 2 * T4 ^ 2) + (1 / 2 : ℚ) • (T1 ^ 2 * T3) + (1 / 2 : ℚ) • (T1 ^ 2 * T3 * T4) + (1 / 4 : ℚ) • (T1 ^ 2 * T3 ^ 2) + (1 / 2 : ℚ) • (T1 ^ 2 * T2 * T4) + (1 / 2 : ℚ) • (T1 ^ 2 * T2 * T3) + (1 / 4 : ℚ) • (T1 ^ 2 * T2 ^ 2) + (1 / 6 : ℚ) • (T1 ^ 3 * T4) + (1 / 6 : ℚ) • (T1 ^ 3 * T3) + (1 / 6 : ℚ) • (T1 ^ 3 * T2), (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_add (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_pow o4 2) (by norm_num : (5 : ℕ) ≤ 2 * 4))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_pow o4 3) (by norm_num : (5 : ℕ) ≤ 3 * 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_pow o4 4) (by norm_num : (5 : ℕ) ≤ 4 * 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul o3 o4) (by norm_num : (5 : ℕ) ≤ 3 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul o3 (Q04.Ord_pow o4 2)) (by norm_num : (5 : ℕ) ≤ 3 + 2 * 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul o3 (Q04.Ord_pow o4 3)) (by norm_num : (5 : ℕ) ≤ 3 + 3 * 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_pow o3 2) (by norm_num : (5 : ℕ) ≤ 2 * 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o3 2) o4) (by norm_num : (5 : ℕ) ≤ 2 * 3 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o3 2) (Q04.Ord_pow o4 2)) (by norm_num : (5 : ℕ) ≤ 2 * 3 + 2 * 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_pow o3 3) (by norm_num : (5 : ℕ) ≤ 3 * 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o3 3) o4) (by norm_num : (5 : ℕ) ≤ 3 * 3 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_pow o3 4) (by norm_num : (5 : ℕ) ≤ 4 * 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul o2 o4) (by norm_num : (5 : ℕ) ≤ 2 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul o2 (Q04.Ord_pow o4 2)) (by norm_num : (5 : ℕ) ≤ 2 + 2 * 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul o2 (Q04.Ord_pow o4 3)) (by norm_num : (5 : ℕ) ≤ 2 + 3 * 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul o2 o3) (by norm_num : (5 : ℕ) ≤ 2 + 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul o2 o3) o4) (by norm_num : (5 : ℕ) ≤ 2 + 3 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul o2 o3) (Q04.Ord_pow o4 2)) (by norm_num : (5 : ℕ) ≤ 2 + 3 + 2 * 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul o2 (Q04.Ord_pow o3 2)) (by norm_num : (5 : ℕ) ≤ 2 + 2 * 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul o2 (Q04.Ord_pow o3 2)) o4) (by norm_num : (5 : ℕ) ≤ 2 + 2 * 3 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul o2 (Q04.Ord_pow o3 3)) (by norm_num : (5 : ℕ) ≤ 2 + 3 * 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o2 2) o4) (by norm_num : (5 : ℕ) ≤ 2 * 2 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o2 2) (Q04.Ord_pow o4 2)) (by norm_num : (5 : ℕ) ≤ 2 * 2 + 2 * 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o2 2) o3) (by norm_num : (5 : ℕ) ≤ 2 * 2 + 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul (Q04.Ord_pow o2 2) o3) o4) (by norm_num : (5 : ℕ) ≤ 2 * 2 + 3 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o2 2) (Q04.Ord_pow o3 2)) (by norm_num : (5 : ℕ) ≤ 2 * 2 + 2 * 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_pow o2 3) (by norm_num : (5 : ℕ) ≤ 3 * 2)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o2 3) o4) (by norm_num : (5 : ℕ) ≤ 3 * 2 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o2 3) o3) (by norm_num : (5 : ℕ) ≤ 3 * 2 + 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_pow o2 4) (by norm_num : (5 : ℕ) ≤ 4 * 2)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul o1 o4) (by norm_num : (5 : ℕ) ≤ 1 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul o1 (Q04.Ord_pow o4 2)) (by norm_num : (5 : ℕ) ≤ 1 + 2 * 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul o1 (Q04.Ord_pow o4 3)) (by norm_num : (5 : ℕ) ≤ 1 + 3 * 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul o1 o3) o4) (by norm_num : (5 : ℕ) ≤ 1 + 3 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul o1 o3) (Q04.Ord_pow o4 2)) (by norm_num : (5 : ℕ) ≤ 1 + 3 + 2 * 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul o1 (Q04.Ord_pow o3 2)) (by norm_num : (5 : ℕ) ≤ 1 + 2 * 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul o1 (Q04.Ord_pow o3 2)) o4) (by norm_num : (5 : ℕ) ≤ 1 + 2 * 3 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul o1 (Q04.Ord_pow o3 3)) (by norm_num : (5 : ℕ) ≤ 1 + 3 * 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul o1 o2) o4) (by norm_num : (5 : ℕ) ≤ 1 + 2 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul o1 o2) (Q04.Ord_pow o4 2)) (by norm_num : (5 : ℕ) ≤ 1 + 2 + 2 * 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul o1 o2) o3) (by norm_num : (5 : ℕ) ≤ 1 + 2 + 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul (Q04.Ord_mul o1 o2) o3) o4) (by norm_num : (5 : ℕ) ≤ 1 + 2 + 3 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul o1 o2) (Q04.Ord_pow o3 2)) (by norm_num : (5 : ℕ) ≤ 1 + 2 + 2 * 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul o1 (Q04.Ord_pow o2 2)) (by norm_num : (5 : ℕ) ≤ 1 + 2 * 2)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul o1 (Q04.Ord_pow o2 2)) o4) (by norm_num : (5 : ℕ) ≤ 1 + 2 * 2 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul o1 (Q04.Ord_pow o2 2)) o3) (by norm_num : (5 : ℕ) ≤ 1 + 2 * 2 + 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul o1 (Q04.Ord_pow o2 3)) (by norm_num : (5 : ℕ) ≤ 1 + 3 * 2)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o1 2) o4) (by norm_num : (5 : ℕ) ≤ 2 * 1 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o1 2) (Q04.Ord_pow o4 2)) (by norm_num : (5 : ℕ) ≤ 2 * 1 + 2 * 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o1 2) o3) (by norm_num : (5 : ℕ) ≤ 2 * 1 + 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul (Q04.Ord_pow o1 2) o3) o4) (by norm_num : (5 : ℕ) ≤ 2 * 1 + 3 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o1 2) (Q04.Ord_pow o3 2)) (by norm_num : (5 : ℕ) ≤ 2 * 1 + 2 * 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul (Q04.Ord_pow o1 2) o2) o4) (by norm_num : (5 : ℕ) ≤ 2 * 1 + 2 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_mul (Q04.Ord_pow o1 2) o2) o3) (by norm_num : (5 : ℕ) ≤ 2 * 1 + 2 + 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o1 2) (Q04.Ord_pow o2 2)) (by norm_num : (5 : ℕ) ≤ 2 * 1 + 2 * 2)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o1 3) o4) (by norm_num : (5 : ℕ) ≤ 3 * 1 + 4)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o1 3) o3) (by norm_num : (5 : ℕ) ≤ 3 * 1 + 3)))) (Q04.Ord_smul _ (Q04.Ord_mono (Q04.Ord_mul (Q04.Ord_pow o1 3) o2) (by norm_num : (5 : ℕ) ≤ 3 * 1 + 2)))), ?_⟩
    unfold Q04.expT
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial]
    norm_num
    algebra
  have hG : (1 : ℚ) • (1) + (1 : ℚ) • (T4) + (1 : ℚ) • (T3) + (1 : ℚ) • (T2) + (1 / 2 : ℚ) • (T2 ^ 2) + (1 : ℚ) • (T1) + (1 : ℚ) • (T1 * T3) + (1 : ℚ) • (T1 * T2) + (1 / 2 : ℚ) • (T1 ^ 2) + (1 / 2 : ℚ) • (T1 ^ 2 * T2) + (1 / 6 : ℚ) • (T1 ^ 3) + (1 / 24 : ℚ) • (T1 ^ 4) = ((1 : ℚ)) • ((X 0 : Q04.Q2) ^ 0 * X 1 ^ 0) + ((1 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1))) • ((X 0 : Q04.Q2) ^ 0 * X 1 ^ 1) + ((1 / 2 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 2 : ℚ) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1))) • ((X 0 : Q04.Q2) ^ 0 * X 1 ^ 2) + ((1 / 6 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 2 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) + (1 / 3 : ℚ) * (V05.HrQ 3 (q - 1) - V05.HrQ 3 (p + q - 1))) • ((X 0 : Q04.Q2) ^ 0 * X 1 ^ 3) + ((1 / 24 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 4 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) + (1 / 3 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 3 (q - 1) - V05.HrQ 3 (p + q - 1)) + (1 / 8 : ℚ) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) + (-1 / 4 : ℚ) * (V05.HrQ 4 (q - 1) - V05.HrQ 4 (p + q - 1))) • ((X 0 : Q04.Q2) ^ 0 * X 1 ^ 4) + ((1 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1))) • ((X 0 : Q04.Q2) ^ 1 * X 1 ^ 0) + ((1 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (1 : ℚ) * V05.HrQ 2 (p + q - 1)) • ((X 0 : Q04.Q2) ^ 1 * X 1 ^ 1) + ((1 / 2 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 2 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) + (1 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 2 (p + q - 1) + (-1 : ℚ) * V05.HrQ 3 (p + q - 1)) • ((X 0 : Q04.Q2) ^ 1 * X 1 ^ 2) + ((1 / 6 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 2 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) + (1 / 3 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 3 (q - 1) - V05.HrQ 3 (p + q - 1)) + (1 / 2 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 2 (p + q - 1) + (-1 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 3 (p + q - 1) + (-1 / 2 : ℚ) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) * V05.HrQ 2 (p + q - 1) + (1 : ℚ) * V05.HrQ 4 (p + q - 1)) • ((X 0 : Q04.Q2) ^ 1 * X 1 ^ 3) + ((1 / 2 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 2 : ℚ) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1))) • ((X 0 : Q04.Q2) ^ 2 * X 1 ^ 0) + ((1 / 2 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (1 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 2 (p + q - 1) + (-1 / 2 : ℚ) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (-1 : ℚ) * V05.HrQ 3 (p + q - 1)) • ((X 0 : Q04.Q2) ^ 2 * X 1 ^ 1) + ((1 / 4 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 4 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) + (1 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 2 (p + q - 1) + (-1 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 3 (p + q - 1) + (-1 / 4 : ℚ) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (1 / 4 : ℚ) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) * (V05.HrQ 2 (q - 1) - V05.HrQ 2 (p + q - 1)) + (-1 : ℚ) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 3 (p + q - 1) + (1 / 2 : ℚ) * V05.HrQ 2 (p + q - 1) * V05.HrQ 2 (p + q - 1) + (3 / 2 : ℚ) * V05.HrQ 4 (p + q - 1)) • ((X 0 : Q04.Q2) ^ 2 * X 1 ^ 2) + ((1 / 6 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 2 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) + (1 / 3 : ℚ) * (V05.HrQ 3 (p - 1) - V05.HrQ 3 (p + q - 1))) • ((X 0 : Q04.Q2) ^ 3 * X 1 ^ 0) + ((1 / 6 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (1 / 2 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 2 (p + q - 1) + (-1 / 2 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (-1 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * V05.HrQ 3 (p + q - 1) + (-1 / 2 : ℚ) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) * V05.HrQ 2 (p + q - 1) + (1 / 3 : ℚ) * (V05.HrQ 3 (p - 1) - V05.HrQ 3 (p + q - 1)) * (V05.HrQ 1 (q - 1) - V05.HrQ 1 (p + q - 1)) + (1 : ℚ) * V05.HrQ 4 (p + q - 1)) • ((X 0 : Q04.Q2) ^ 3 * X 1 ^ 1) + ((1 / 24 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) + (-1 / 4 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) + (1 / 3 : ℚ) * (V05.HrQ 1 (p - 1) - V05.HrQ 1 (p + q - 1)) * (V05.HrQ 3 (p - 1) - V05.HrQ 3 (p + q - 1)) + (1 / 8 : ℚ) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) * (V05.HrQ 2 (p - 1) - V05.HrQ 2 (p + q - 1)) + (-1 / 4 : ℚ) * (V05.HrQ 4 (p - 1) - V05.HrQ 4 (p + q - 1))) • ((X 0 : Q04.Q2) ^ 4 * X 1 ^ 0) := by
    rw [hT1, hT2, hT3, hT4]
    norm_num
    algebra
  unfold Q04.gH
  rw [show i + (4 - i) = 4 by omega, hHb, hE, hG, _root_.MvPolynomial.coeff_add,
    Q04.coeff_eq_of_Ord (d := 4) (F := R) (G := 0) (by rw [sub_zero]; exact hR) _ (by simp; omega)]
  simp
end gH4sec

/-- the lists for `[1]D^♯_{i,4-i}(n+j,n) = i!(4-i)! gH` in the atoms -/
noncomputable def Dz1 (j : ℕ) : ℕ → List S03.DTerm
  | 0 => (atm 1 4 j).map (scl (-6 : ℚ) 0 0) ++
      (mulL (atm 1 1 j) (atm 1 3 j)).map (scl (8 : ℚ) 0 0) ++
      (mulL (atm 1 2 j) (atm 1 2 j)).map (scl (3 : ℚ) 0 0) ++
      (mulL (atm 1 1 j) (mulL (atm 1 1 j) (atm 1 2 j))).map (scl (-6 : ℚ) 0 0) ++
      (mulL (atm 1 1 j) (mulL (atm 1 1 j) (mulL (atm 1 1 j) (atm 1 1 j)))).map (scl (1 : ℚ) 0 0)
  | 1 => (atm 2 4 j).map (scl (6 : ℚ) 0 0) ++
      (mulL (atm 0 1 j) (atm 1 3 j)).map (scl (2 : ℚ) 0 0) ++
      (mulL (atm 1 1 j) (atm 2 3 j)).map (scl (-6 : ℚ) 0 0) ++
      (mulL (atm 1 2 j) (atm 2 2 j)).map (scl (-3 : ℚ) 0 0) ++
      (mulL (atm 0 1 j) (mulL (atm 1 1 j) (atm 1 2 j))).map (scl (-3 : ℚ) 0 0) ++
      (mulL (atm 1 1 j) (mulL (atm 1 1 j) (atm 2 2 j))).map (scl (3 : ℚ) 0 0) ++
      (mulL (atm 0 1 j) (mulL (atm 1 1 j) (mulL (atm 1 1 j) (atm 1 1 j)))).map (scl (1 : ℚ) 0 0)
  | 2 => (atm 2 4 j).map (scl (6 : ℚ) 0 0) ++
      (mulL (atm 0 1 j) (atm 2 3 j)).map (scl (-4 : ℚ) 0 0) ++
      (mulL (atm 1 1 j) (atm 2 3 j)).map (scl (-4 : ℚ) 0 0) ++
      (mulL (atm 0 2 j) (atm 1 2 j)).map (scl (1 : ℚ) 0 0) ++
      (mulL (atm 2 2 j) (atm 2 2 j)).map (scl (2 : ℚ) 0 0) ++
      (mulL (atm 0 1 j) (mulL (atm 0 1 j) (atm 1 2 j))).map (scl (-1 : ℚ) 0 0) ++
      (mulL (atm 0 1 j) (mulL (atm 1 1 j) (atm 2 2 j))).map (scl (4 : ℚ) 0 0) ++
      (mulL (atm 0 2 j) (mulL (atm 1 1 j) (atm 1 1 j))).map (scl (-1 : ℚ) 0 0) ++
      (mulL (atm 0 1 j) (mulL (atm 0 1 j) (mulL (atm 1 1 j) (atm 1 1 j)))).map (scl (1 : ℚ) 0 0)
  | 3 => (atm 2 4 j).map (scl (6 : ℚ) 0 0) ++
      (mulL (atm 0 1 j) (atm 2 3 j)).map (scl (-6 : ℚ) 0 0) ++
      (mulL (atm 0 3 j) (atm 1 1 j)).map (scl (2 : ℚ) 0 0) ++
      (mulL (atm 0 2 j) (atm 2 2 j)).map (scl (-3 : ℚ) 0 0) ++
      (mulL (atm 0 1 j) (mulL (atm 0 1 j) (atm 2 2 j))).map (scl (3 : ℚ) 0 0) ++
      (mulL (atm 0 1 j) (mulL (atm 0 2 j) (atm 1 1 j))).map (scl (-3 : ℚ) 0 0) ++
      (mulL (atm 0 1 j) (mulL (atm 0 1 j) (mulL (atm 0 1 j) (atm 1 1 j)))).map (scl (1 : ℚ) 0 0)
  | 4 => (atm 0 4 j).map (scl (-6 : ℚ) 0 0) ++
      (mulL (atm 0 1 j) (atm 0 3 j)).map (scl (8 : ℚ) 0 0) ++
      (mulL (atm 0 2 j) (atm 0 2 j)).map (scl (3 : ℚ) 0 0) ++
      (mulL (atm 0 1 j) (mulL (atm 0 1 j) (atm 0 2 j))).map (scl (-6 : ℚ) 0 0) ++
      (mulL (atm 0 1 j) (mulL (atm 0 1 j) (mulL (atm 0 1 j) (atm 0 1 j)))).map (scl (1 : ℚ) 0 0)
  | _ => []

lemma idx_mulL3 {A B C : List S03.DTerm} (ha : ∀ τ ∈ A, S03.IsIndex τ.2.2)
    (hb : ∀ τ ∈ B, S03.IsIndex τ.2.2) (hc : ∀ τ ∈ C, S03.IsIndex τ.2.2) :
    ∀ τ ∈ mulL A (mulL B C), S03.IsIndex τ.2.2 := idx_mulL ha (idx_mulL hb hc)

lemma ST_m2 (k r j k' r' j' : ℕ) (N : ℕ) :
    S03.ST (mulL (atm k r j) (atm k' r' j')) N = S03.ST (atm k r j) N * S03.ST (atm k' r' j') N :=
  ST_mulL _ _ (atm_idx _ _ _) (atm_idx _ _ _) N

lemma ST_m3 (k r j k' r' j' k'' r'' j'' : ℕ) (N : ℕ) :
    S03.ST (mulL (atm k r j) (mulL (atm k' r' j') (atm k'' r'' j''))) N =
      S03.ST (atm k r j) N * (S03.ST (atm k' r' j') N * S03.ST (atm k'' r'' j'') N) := by
  rw [ST_mulL _ _ (atm_idx _ _ _) (idx_mulL (atm_idx _ _ _) (atm_idx _ _ _)), ST_m2]

lemma ST_m4 (k r j k' r' j' k'' r'' j'' k3 r3 j3 : ℕ) (N : ℕ) :
    S03.ST (mulL (atm k r j) (mulL (atm k' r' j') (mulL (atm k'' r'' j'') (atm k3 r3 j3)))) N =
      S03.ST (atm k r j) N * (S03.ST (atm k' r' j') N *
        (S03.ST (atm k'' r'' j'') N * S03.ST (atm k3 r3 j3) N)) := by
  rw [ST_mulL _ _ (atm_idx _ _ _) (idx_mulL (atm_idx _ _ _) (idx_mulL (atm_idx _ _ _) (atm_idx _ _ _))),
    ST_m3]

lemma Dz1_idx (j i : ℕ) : ∀ τ ∈ Dz1 j i, S03.IsIndex τ.2.2 := by
  unfold Dz1
  split <;>
  repeat (first | exact atm_idx _ _ _ | apply idx_append | apply idx_map_scl | apply idx_mulL |
    (intro τ h; simp at h))

set_option maxHeartbeats 4000000 in
lemma Dz1_val (j i n : ℕ) (hj : j ≤ 4) (hi : i ≤ 4) (hn : 1 ≤ n) :
    S03.ST (Dz1 j i) (2 * n) = (Q04.DsharpF i (4 - i) (n + j) n).coeff 0 := by
  rw [R01.Dsharp_coeff0, gH_four (n + j) n i hi]
  have a1 := atmA_ST 1 j n le_rfl (by norm_num) hj hn
  have a2 := atmA_ST 2 j n (by norm_num) (by norm_num) hj hn
  have a3 := atmA_ST 3 j n (by norm_num) (by norm_num) hj hn
  have a4 := atmA_ST 4 j n (by norm_num) (by norm_num) hj hn
  have b1 := atmB_ST 1 j n le_rfl (by norm_num) hj hn
  have b2 := atmB_ST 2 j n (by norm_num) (by norm_num) hj hn
  have b3 := atmB_ST 3 j n (by norm_num) (by norm_num) hj hn
  have b4 := atmB_ST 4 j n (by norm_num) (by norm_num) hj hn
  have c2 := atmG_ST 2 j n le_rfl (by norm_num) hj hn
  have c3 := atmG_ST 3 j n (by norm_num) (by norm_num) hj hn
  have c4 := atmG_ST 4 j n (by norm_num) (by norm_num) hj hn
  interval_cases i <;>
  · simp only [_root_.MvPolynomial.coeff_add, coeff_sXX]
    simp only [Dz1, ST_append, ST_map_scl, ST_m2, ST_m3, ST_m4, a1, a2, a3, a4, b1, b2, b3, b4, c2, c3, c4]
    norm_num
    ring

/-- the element `T_z` of `g8:lem:mu` for `𝔷 = 1`, `k = 4` -/
noncomputable def Tz1 : List S03.DTerm :=
  (List.range 5).flatMap fun j => (List.range 5).flatMap fun i => mulL [fz j i] (Dz1 j i)

lemma Tz1_ok (n : ℕ) (hn : 1 ≤ n) : (Q04.Pf 4 n).coeff 0 = S03.ST Tz1 (2 * n) := by
  unfold Q04.Pf Tz1
  rw [MvPolynomial.coeff_smul, smul_eq_mul, ST_flatMap, list_range_sum, MvPolynomial.coeff_sum,
    Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro j hj
  rw [ST_flatMap, list_range_sum, MvPolynomial.coeff_sum, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i hi
  have hj' : j ≤ 4 := by simp at hj; omega
  have hi' : i ≤ 4 := by simp at hi; omega
  rw [ST_single_nil _ (fz_w j i) _ (Dz1_idx j i), fz_eval j i n hj' hi' hn,
    MvPolynomial.coeff_smul, smul_eq_mul, ← Dz1_val j i n hj' hi' hn, S04.eF4 n j i hj' hi']
  push_cast
  ring

end P04

namespace P04

/-- the product `[a_1] ⋯ [a_k]` in `𝔇` -/
noncomputable def sp : List ℤ → List S03.DTerm
  | [] => [(1, 0, [])]
  | a :: m => mulL [(1, 0, [a])] (sp m)

lemma sp_idx : ∀ m : List ℤ, (∀ a ∈ m, a ≠ 0) → ∀ τ ∈ sp m, S03.IsIndex τ.2.2
  | [], _ => by intro τ h; simp [sp] at h; subst h; intro a ha; simp at ha
  | a :: m, hm => by
    refine idx_mulL ?_ (sp_idx m (fun b hb => hm b (by simp [hb])))
    intro τ h; simp at h; subst h; intro b hb; simp at hb; rw [hb]; exact hm a (by simp)

lemma ST_sp : ∀ m : List ℤ, (∀ a ∈ m, a ≠ 0) → ∀ N : ℕ,
    S03.ST (sp m) N = (m.map fun a => S03.Z [a] N).prod
  | [], _, N => by simp [sp, S03.ST, ratEval_one, Z_nil]
  | a :: m, hm, N => by
    rw [sp, ST_mulL _ _ ?_ (sp_idx m (fun b hb => hm b (by simp [hb]))),
      ST_sp m (fun b hb => hm b (by simp [hb])) N]
    · simp [S03.ST, ratEval_one]
    · intro τ h; simp at h; subst h; intro b hb; simp at hb; rw [hb]; exact hm a (by simp)

lemma ST_rt_sp (cs : List ℚ) (ps : List (ℤ × ℕ)) (m : List ℤ) (hm : ∀ a ∈ m, a ≠ 0) (N : ℕ) :
    S03.ST (mulL [rt cs ps] (sp m)) N =
      S03.ratEval (rt cs ps).1 (rt cs ps).2.1 (N : ℚ) * (m.map fun a => S03.Z [a] N).prod := by
  rw [ST_single_nil _ (rt_w cs ps) _ (sp_idx m hm), ST_sp m hm]

lemma Dval_mulL_w0 (τ : S03.DTerm) (hτ : τ.2.2 = []) (L : List S03.DTerm) (u : List ℤ) (N : ℚ) :
    Dval (mulL [τ] L) u N = S03.ratEval τ.1 τ.2.1 N * Dval L u N := by
  unfold mulL
  rw [List.flatMap_singleton, Dval_flatMap, hτ]
  induction L with
  | nil => simp [Dval]
  | cons τ' L ih =>
    rw [List.map_cons, List.sum_cons, ih, Dval_cons, stuf_nil_left, List.map_singleton, Dval_cons,
      Dval_nil]
    simp only [S08.PA.ratEval_mul]
    split_ifs <;> ring

lemma Dval_mulL_rt (cs : List ℚ) (ps : List (ℤ × ℕ)) (L : List S03.DTerm) (u : List ℤ) (N : ℚ) :
    Dval (mulL [rt cs ps] L) u N = S03.ratEval (rt cs ps).1 (rt cs ps).2.1 N * Dval L u N :=
  Dval_mulL_w0 _ (rt_w cs ps) L u N

/-- the blocks `i!(4-i)! [x^i x'^{4-i}] exp(Δ)` expanded in single-letter nested sums -/
noncomputable def E1 : ℕ → ℕ → List S03.DTerm
  | 0, 0 => mulL [rt [192] [(0, 4)]] (sp []) ++
      mulL [rt [-48] []] (sp [-4]) ++
      mulL [rt [-32] [(0, 1)]] (sp [-3]) ++
      mulL [rt [32] []] (sp [-3, -1]) ++
      mulL [rt [-48] [(0, 2)]] (sp [-2]) ++
      mulL [rt [12] []] (sp [-2, -2]) ++
      mulL [rt [24] [(0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-12] []] (sp [-2, -1, -1]) ++
      mulL [rt [12] []] (sp [-2, 2]) ++
      mulL [rt [-96] [(0, 3)]] (sp [-1]) ++
      mulL [rt [24] [(0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-4] [(0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [-6] []] (sp [-1, -1, 2]) ++
      mulL [rt [12] [(0, 1)]] (sp [-1, 2]) ++
      mulL [rt [24] []] (sp [-1, 3]) ++
      mulL [rt [-24] [(0, 2)]] (sp [2]) ++
      mulL [rt [3] []] (sp [2, 2]) ++
      mulL [rt [-24] [(0, 1)]] (sp [3]) ++
      mulL [rt [-42] []] (sp [4])
  | 0, 1 => mulL [rt [-8] [(0, 1)]] (sp [-3]) ++
      mulL [rt [8] []] (sp [-3, -1]) ++
      mulL [rt [12] [(0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-6] []] (sp [-2, -1, -1]) ++
      mulL [rt [-6] []] (sp [-2, 2]) ++
      mulL [rt [-24] [(0, 3)]] (sp [-1]) ++
      mulL [rt [12] [(0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-4] [(0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [12] [(0, 2)]] (sp [2]) ++
      mulL [rt [-3] []] (sp [2, 2]) ++
      mulL [rt [6] []] (sp [4])
  | 0, 2 => mulL [rt [-16] [(0, 2)]] (sp [-2]) ++
      mulL [rt [4] []] (sp [-2, -2]) ++
      mulL [rt [8] [(0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-4] []] (sp [-2, -1, -1]) ++
      mulL [rt [4] []] (sp [-2, 2]) ++
      mulL [rt [8] [(0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-4] [(0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [2] []] (sp [-1, -1, 2]) ++
      mulL [rt [-4] [(0, 1)]] (sp [-1, 2]) ++
      mulL [rt [-8] []] (sp [-1, 3]) ++
      mulL [rt [-8] [(0, 2)]] (sp [2]) ++
      mulL [rt [3] []] (sp [2, 2]) ++
      mulL [rt [8] [(0, 1)]] (sp [3]) ++
      mulL [rt [6] []] (sp [4])
  | 0, 3 => mulL [rt [-8] [(0, 1)]] (sp [-3]) ++
      mulL [rt [8] []] (sp [-3, -1]) ++
      mulL [rt [12] [(0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-6] []] (sp [-2, -1, -1]) ++
      mulL [rt [-6] []] (sp [-2, 2]) ++
      mulL [rt [-24] [(0, 3)]] (sp [-1]) ++
      mulL [rt [12] [(0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-4] [(0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [12] [(0, 2)]] (sp [2]) ++
      mulL [rt [-3] []] (sp [2, 2]) ++
      mulL [rt [6] []] (sp [4])
  | 0, 4 => mulL [rt [192] [(0, 4)]] (sp []) ++
      mulL [rt [-48] []] (sp [-4]) ++
      mulL [rt [-32] [(0, 1)]] (sp [-3]) ++
      mulL [rt [32] []] (sp [-3, -1]) ++
      mulL [rt [-48] [(0, 2)]] (sp [-2]) ++
      mulL [rt [12] []] (sp [-2, -2]) ++
      mulL [rt [24] [(0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-12] []] (sp [-2, -1, -1]) ++
      mulL [rt [12] []] (sp [-2, 2]) ++
      mulL [rt [-96] [(0, 3)]] (sp [-1]) ++
      mulL [rt [24] [(0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-4] [(0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [-6] []] (sp [-1, -1, 2]) ++
      mulL [rt [12] [(0, 1)]] (sp [-1, 2]) ++
      mulL [rt [24] []] (sp [-1, 3]) ++
      mulL [rt [-24] [(0, 2)]] (sp [2]) ++
      mulL [rt [3] []] (sp [2, 2]) ++
      mulL [rt [-24] [(0, 1)]] (sp [3]) ++
      mulL [rt [-42] []] (sp [4])
  | 1, 0 => mulL [rt [384] [(0, 4)]] (sp []) ++
      mulL [rt [-48] []] (sp [-4]) ++
      mulL [rt [-64] [(0, 1)]] (sp [-3]) ++
      mulL [rt [32] []] (sp [-3, -1]) ++
      mulL [rt [-96] [(0, 2)]] (sp [-2]) ++
      mulL [rt [12] []] (sp [-2, -2]) ++
      mulL [rt [48] [(0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-12] []] (sp [-2, -1, -1]) ++
      mulL [rt [12] []] (sp [-2, 2]) ++
      mulL [rt [-192] [(0, 3)]] (sp [-1]) ++
      mulL [rt [48] [(0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-8] [(0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [-6] []] (sp [-1, -1, 2]) ++
      mulL [rt [24] [(0, 1)]] (sp [-1, 2]) ++
      mulL [rt [24] []] (sp [-1, 3]) ++
      mulL [rt [-48] [(0, 2)]] (sp [2]) ++
      mulL [rt [3] []] (sp [2, 2]) ++
      mulL [rt [-48] [(0, 1)]] (sp [3]) ++
      mulL [rt [-42] []] (sp [4])
  | 1, 1 => mulL [rt [8] []] (sp [-3, -1]) ++
      mulL [rt [12] [(0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-6] []] (sp [-2, -1, -1]) ++
      mulL [rt [-6] []] (sp [-2, 2]) ++
      mulL [rt [-48] [(0, 3)]] (sp [-1]) ++
      mulL [rt [24] [(0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-6] [(0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [-6] [(0, 1)]] (sp [-1, 2]) ++
      mulL [rt [24] [(0, 2)]] (sp [2]) ++
      mulL [rt [-3] []] (sp [2, 2]) ++
      mulL [rt [12] [(0, 1)]] (sp [3]) ++
      mulL [rt [6] []] (sp [4])
  | 1, 2 => mulL [rt [-16] [(0, 2)]] (sp [-2]) ++
      mulL [rt [4] []] (sp [-2, -2]) ++
      mulL [rt [8] [(0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-4] []] (sp [-2, -1, -1]) ++
      mulL [rt [4] []] (sp [-2, 2]) ++
      mulL [rt [8] [(0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-4] [(0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [2] []] (sp [-1, -1, 2]) ++
      mulL [rt [-4] [(0, 1)]] (sp [-1, 2]) ++
      mulL [rt [-8] []] (sp [-1, 3]) ++
      mulL [rt [-8] [(0, 2)]] (sp [2]) ++
      mulL [rt [3] []] (sp [2, 2]) ++
      mulL [rt [8] [(0, 1)]] (sp [3]) ++
      mulL [rt [6] []] (sp [4])
  | 1, 3 => mulL [rt [-16] [(0, 1)]] (sp [-3]) ++
      mulL [rt [8] []] (sp [-3, -1]) ++
      mulL [rt [12] [(0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-6] []] (sp [-2, -1, -1]) ++
      mulL [rt [-6] []] (sp [-2, 2]) ++
      mulL [rt [-2] [(0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [6] [(0, 1)]] (sp [-1, 2]) ++
      mulL [rt [-3] []] (sp [2, 2]) ++
      mulL [rt [-12] [(0, 1)]] (sp [3]) ++
      mulL [rt [6] []] (sp [4])
  | 1, 4 => mulL [rt [-48] []] (sp [-4]) ++
      mulL [rt [32] []] (sp [-3, -1]) ++
      mulL [rt [12] []] (sp [-2, -2]) ++
      mulL [rt [-12] []] (sp [-2, -1, -1]) ++
      mulL [rt [12] []] (sp [-2, 2]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [-6] []] (sp [-1, -1, 2]) ++
      mulL [rt [24] []] (sp [-1, 3]) ++
      mulL [rt [3] []] (sp [2, 2]) ++
      mulL [rt [-42] []] (sp [4])
  | 2, 0 => mulL [rt [384, 1728, 2976, 2352, 744] [(-1, 4), (0, 4)]] (sp []) ++
      mulL [rt [-48] []] (sp [-4]) ++
      mulL [rt [-64, -96] [(-1, 1), (0, 1)]] (sp [-3]) ++
      mulL [rt [32] []] (sp [-3, -1]) ++
      mulL [rt [-96, -240, -168] [(-1, 2), (0, 2)]] (sp [-2]) ++
      mulL [rt [12] []] (sp [-2, -2]) ++
      mulL [rt [48, 72] [(-1, 1), (0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-12] []] (sp [-2, -1, -1]) ++
      mulL [rt [12] []] (sp [-2, 2]) ++
      mulL [rt [-192, -672, -816, -360] [(-1, 3), (0, 3)]] (sp [-1]) ++
      mulL [rt [48, 120, 84] [(-1, 2), (0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-8, -12] [(-1, 1), (0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [-6] []] (sp [-1, -1, 2]) ++
      mulL [rt [24, 36] [(-1, 1), (0, 1)]] (sp [-1, 2]) ++
      mulL [rt [24] []] (sp [-1, 3]) ++
      mulL [rt [-48, -120, -84] [(-1, 2), (0, 2)]] (sp [2]) ++
      mulL [rt [3] []] (sp [2, 2]) ++
      mulL [rt [-48, -72] [(-1, 1), (0, 1)]] (sp [3]) ++
      mulL [rt [-42] []] (sp [4])
  | 2, 1 => mulL [rt [0, 0, 0, 0, 0, -24] [(-2, 1), (-1, 4), (0, 4)]] (sp []) ++
      mulL [rt [0, 0, 8] [(-2, 1), (-1, 1), (0, 1)]] (sp [-3]) ++
      mulL [rt [8] []] (sp [-3, -1]) ++
      mulL [rt [0, 0, 0, 12] [(-2, 1), (-1, 2), (0, 2)]] (sp [-2]) ++
      mulL [rt [24, 48, 12] [(-2, 1), (-1, 1), (0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-6] []] (sp [-2, -1, -1]) ++
      mulL [rt [-6] []] (sp [-2, 2]) ++
      mulL [rt [-96, -384, -576, -384, -72] [(-2, 1), (-1, 3), (0, 3)]] (sp [-1]) ++
      mulL [rt [48, 144, 144, 36] [(-2, 1), (-1, 2), (0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-12, -24, -8] [(-2, 1), (-1, 1), (0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [-12, -24, -12] [(-2, 1), (-1, 1), (0, 1)]] (sp [-1, 2]) ++
      mulL [rt [48, 144, 144, 48] [(-2, 1), (-1, 2), (0, 2)]] (sp [2]) ++
      mulL [rt [-3] []] (sp [2, 2]) ++
      mulL [rt [24, 48, 24] [(-2, 1), (-1, 1), (0, 1)]] (sp [3]) ++
      mulL [rt [6] []] (sp [4])
  | 2, 2 => mulL [rt [0, 0, 0, 0, 0, -48, -24] [(-2, 2), (-1, 4), (0, 4)]] (sp []) ++
      mulL [rt [-64, -224, -288, -144, -24] [(-2, 2), (-1, 2), (0, 2)]] (sp [-2]) ++
      mulL [rt [4] []] (sp [-2, -2]) ++
      mulL [rt [16, 32, 8] [(-2, 1), (-1, 1), (0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-4] []] (sp [-2, -1, -1]) ++
      mulL [rt [4] []] (sp [-2, 2]) ++
      mulL [rt [0, 0, 0, 0, 48, 24] [(-2, 2), (-1, 3), (0, 3)]] (sp [-1]) ++
      mulL [rt [32, 112, 144, 56, 4] [(-2, 2), (-1, 2), (0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-8, -16, -4] [(-2, 1), (-1, 1), (0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [2] []] (sp [-1, -1, 2]) ++
      mulL [rt [-8, -16, -4] [(-2, 1), (-1, 1), (0, 1)]] (sp [-1, 2]) ++
      mulL [rt [-8] []] (sp [-1, 3]) ++
      mulL [rt [-32, -112, -144, -88, -20] [(-2, 2), (-1, 2), (0, 2)]] (sp [2]) ++
      mulL [rt [3] []] (sp [2, 2]) ++
      mulL [rt [16, 32, 8] [(-2, 1), (-1, 1), (0, 1)]] (sp [3]) ++
      mulL [rt [6] []] (sp [4])
  | 2, 3 => mulL [rt [0, 0, 0, 0, 0, -96, -96, -24] [(-2, 3), (-1, 4), (0, 4)]] (sp []) ++
      mulL [rt [-16, -24] [(-1, 1), (0, 1)]] (sp [-3]) ++
      mulL [rt [8] []] (sp [-3, -1]) ++
      mulL [rt [0, 0, 0, 12] [(-2, 1), (-1, 2), (0, 2)]] (sp [-2]) ++
      mulL [rt [24, 48, 12] [(-2, 1), (-1, 1), (0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-6] []] (sp [-2, -1, -1]) ++
      mulL [rt [-6] []] (sp [-2, 2]) ++
      mulL [rt [0, 0, 0, 0, 96, 96, 24] [(-2, 3), (-1, 3), (0, 3)]] (sp [-1]) ++
      mulL [rt [0, 0, 0, -24, -12] [(-2, 2), (-1, 2), (0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-4, -8] [(-2, 1), (-1, 1), (0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [12, 24, 12] [(-2, 1), (-1, 1), (0, 1)]] (sp [-1, 2]) ++
      mulL [rt [-3] []] (sp [2, 2]) ++
      mulL [rt [-24, -48, -24] [(-2, 1), (-1, 1), (0, 1)]] (sp [3]) ++
      mulL [rt [6] []] (sp [4])
  | 2, 4 => mulL [rt [0, 0, 0, 0, 0, -192, -288, -144, -24] [(-2, 4), (-1, 4), (0, 4)]] (sp []) ++
      mulL [rt [-48] []] (sp [-4]) ++
      mulL [rt [0, 0, 32] [(-2, 1), (-1, 1), (0, 1)]] (sp [-3]) ++
      mulL [rt [32] []] (sp [-3, -1]) ++
      mulL [rt [0, 0, 0, 48, 24] [(-2, 2), (-1, 2), (0, 2)]] (sp [-2]) ++
      mulL [rt [12] []] (sp [-2, -2]) ++
      mulL [rt [0, 0, -24] [(-2, 1), (-1, 1), (0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-12] []] (sp [-2, -1, -1]) ++
      mulL [rt [12] []] (sp [-2, 2]) ++
      mulL [rt [0, 0, 0, 0, 96, 96, 24] [(-2, 3), (-1, 3), (0, 3)]] (sp [-1]) ++
      mulL [rt [0, 0, 0, -24, -12] [(-2, 2), (-1, 2), (0, 2)]] (sp [-1, -1]) ++
      mulL [rt [0, 0, 4] [(-2, 1), (-1, 1), (0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [-6] []] (sp [-1, -1, 2]) ++
      mulL [rt [0, 0, -12] [(-2, 1), (-1, 1), (0, 1)]] (sp [-1, 2]) ++
      mulL [rt [24] []] (sp [-1, 3]) ++
      mulL [rt [0, 0, 0, 24, 12] [(-2, 2), (-1, 2), (0, 2)]] (sp [2]) ++
      mulL [rt [3] []] (sp [2, 2]) ++
      mulL [rt [0, 0, 24] [(-2, 1), (-1, 1), (0, 1)]] (sp [3]) ++
      mulL [rt [-42] []] (sp [4])
  | 3, 0 => mulL [rt [6144, 41472, 121728, 203040, 210744, 139320, 57264, 13392, 1368] [(-2, 4), (-1, 4), (0, 4)]] (sp []) ++
      mulL [rt [-48] []] (sp [-4]) ++
      mulL [rt [-128, -288, -128] [(-2, 1), (-1, 1), (0, 1)]] (sp [-3]) ++
      mulL [rt [32] []] (sp [-3, -1]) ++
      mulL [rt [-384, -1440, -2040, -1224, -264] [(-2, 2), (-1, 2), (0, 2)]] (sp [-2]) ++
      mulL [rt [12] []] (sp [-2, -2]) ++
      mulL [rt [96, 216, 96] [(-2, 1), (-1, 1), (0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-12] []] (sp [-2, -1, -1]) ++
      mulL [rt [12] []] (sp [-2, 2]) ++
      mulL [rt [-1536, -8064, -17568, -20376, -13152, -4464, -624] [(-2, 3), (-1, 3), (0, 3)]] (sp [-1]) ++
      mulL [rt [192, 720, 1020, 612, 132] [(-2, 2), (-1, 2), (0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-16, -36, -16] [(-2, 1), (-1, 1), (0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [-6] []] (sp [-1, -1, 2]) ++
      mulL [rt [48, 108, 48] [(-2, 1), (-1, 1), (0, 1)]] (sp [-1, 2]) ++
      mulL [rt [24] []] (sp [-1, 3]) ++
      mulL [rt [-192, -720, -1020, -612, -132] [(-2, 2), (-1, 2), (0, 2)]] (sp [2]) ++
      mulL [rt [3] []] (sp [2, 2]) ++
      mulL [rt [-96, -216, -96] [(-2, 1), (-1, 1), (0, 1)]] (sp [3]) ++
      mulL [rt [-42] []] (sp [4])
  | 3, 1 => mulL [rt [0, 0, 0, 0, 0, -1332, -2820, -2280, -840, -120] [(-4, 1), (-2, 4), (-1, 4), (0, 4)]] (sp []) ++
      mulL [rt [0, 0, 40, 16] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-3]) ++
      mulL [rt [8] []] (sp [-3, -1]) ++
      mulL [rt [0, 0, 0, 144, 144, 36] [(-4, 1), (-2, 2), (-1, 2), (0, 2)]] (sp [-2]) ++
      mulL [rt [96, 240, 120, 12] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-6] []] (sp [-2, -1, -1]) ++
      mulL [rt [-6] []] (sp [-2, 2]) ++
      mulL [rt [-1536, -8448, -19584, -24768, -17772, -7008, -1344, -84] [(-4, 1), (-2, 3), (-1, 3), (0, 3)]] (sp [-1]) ++
      mulL [rt [384, 1536, 2400, 1662, 498, 48] [(-4, 1), (-2, 2), (-1, 2), (0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-48, -120, -70, -10] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [-48, -120, -90, -18] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, 2]) ++
      mulL [rt [384, 1536, 2400, 1806, 642, 84] [(-4, 1), (-2, 2), (-1, 2), (0, 2)]] (sp [2]) ++
      mulL [rt [-3] []] (sp [2, 2]) ++
      mulL [rt [96, 240, 180, 36] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [3]) ++
      mulL [rt [6] []] (sp [4])
  | 3, 2 => mulL [rt [0, 0, 0, 0, 0, -4608, -10368, -9216, -4032, -864, -72] [(-4, 2), (-2, 4), (-1, 4), (0, 4)]] (sp []) ++
      mulL [rt [-1024, -4352, -7424, -6080, -2528, -512, -40] [(-4, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [-2]) ++
      mulL [rt [4] []] (sp [-2, -2]) ++
      mulL [rt [64, 160, 80, 8] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-4] []] (sp [-2, -1, -1]) ++
      mulL [rt [4] []] (sp [-2, 2]) ++
      mulL [rt [0, 0, 0, 0, 2304, 4032, 2592, 720, 72] [(-4, 2), (-2, 3), (-1, 3), (0, 3)]] (sp [-1]) ++
      mulL [rt [512, 2176, 3712, 2656, 784, 64, -4] [(-4, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-32, -80, -40, -4] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [2] []] (sp [-1, -1, 2]) ++
      mulL [rt [-32, -80, -40, -4] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, 2]) ++
      mulL [rt [-8] []] (sp [-1, 3]) ++
      mulL [rt [-512, -2176, -3712, -3424, -1744, -448, -44] [(-4, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [2]) ++
      mulL [rt [3] []] (sp [2, 2]) ++
      mulL [rt [64, 160, 80, 8] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [3]) ++
      mulL [rt [6] []] (sp [4])
  | 3, 3 => mulL [rt [0, 0, 0, 0, 0, -15552, -36384, -34092, -16188, -4056, -504, -24] [(-4, 3), (-2, 4), (-1, 4), (0, 4)]] (sp []) ++
      mulL [rt [-32, -72, -32] [(-2, 1), (-1, 1), (0, 1)]] (sp [-3]) ++
      mulL [rt [8] []] (sp [-3, -1]) ++
      mulL [rt [0, 0, 0, 144, 144, 36] [(-4, 1), (-2, 2), (-1, 2), (0, 2)]] (sp [-2]) ++
      mulL [rt [96, 240, 120, 12] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-6] []] (sp [-2, -1, -1]) ++
      mulL [rt [-6] []] (sp [-2, 2]) ++
      mulL [rt [0, 0, 0, 0, 7872, 14688, 10476, 3552, 576, 36] [(-4, 3), (-2, 3), (-1, 3), (0, 3)]] (sp [-1]) ++
      mulL [rt [0, 0, 0, -504, -582, -210, -24] [(-4, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-16, -40, -10, 2] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [48, 120, 90, 18] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, 2]) ++
      mulL [rt [0, 0, 0, 72, 138, 78, 12] [(-4, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [2]) ++
      mulL [rt [-3] []] (sp [2, 2]) ++
      mulL [rt [-96, -240, -180, -36] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [3]) ++
      mulL [rt [6] []] (sp [4])
  | 3, 4 => mulL [rt [0, 0, 0, 0, 0, -50688, -119424, -110496, -49080, -9528, -48, 240, 24] [(-4, 4), (-2, 4), (-1, 4), (0, 4)]] (sp []) ++
      mulL [rt [-48] []] (sp [-4]) ++
      mulL [rt [0, 0, 160, 64] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-3]) ++
      mulL [rt [32] []] (sp [-3, -1]) ++
      mulL [rt [0, 0, 0, 864, 888, 264, 24] [(-4, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [-2]) ++
      mulL [rt [12] []] (sp [-2, -2]) ++
      mulL [rt [0, 0, -120, -48] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-12] []] (sp [-2, -1, -1]) ++
      mulL [rt [12] []] (sp [-2, 2]) ++
      mulL [rt [0, 0, 0, 0, 6528, 10944, 6552, 1632, 144] [(-4, 3), (-2, 3), (-1, 3), (0, 3)]] (sp [-1]) ++
      mulL [rt [0, 0, 0, -432, -444, -132, -12] [(-4, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [-1, -1]) ++
      mulL [rt [0, 0, 20, 8] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [-6] []] (sp [-1, -1, 2]) ++
      mulL [rt [0, 0, -60, -24] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, 2]) ++
      mulL [rt [24] []] (sp [-1, 3]) ++
      mulL [rt [0, 0, 0, 432, 444, 132, 12] [(-4, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [2]) ++
      mulL [rt [3] []] (sp [2, 2]) ++
      mulL [rt [0, 0, 120, 48] [(-4, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [3]) ++
      mulL [rt [-42] []] (sp [4])
  | 4, 0 => mulL [rt [497664, 4105728, 15327360, 34241760, 51004104, 53354880, 40175928, 21932640, 8613408, 2373120, 435480, 47808, 2376] [(-3, 4), (-2, 4), (-1, 4), (0, 4)]] (sp []) ++
      mulL [rt [-48] []] (sp [-4]) ++
      mulL [rt [-384, -1056, -768, -160] [(-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-3]) ++
      mulL [rt [32] []] (sp [-3, -1]) ++
      mulL [rt [-3456, -15840, -29832, -28512, -14520, -3744, -384] [(-3, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [-2]) ++
      mulL [rt [12] []] (sp [-2, -2]) ++
      mulL [rt [288, 792, 576, 120] [(-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-12] []] (sp [-2, -1, -1]) ++
      mulL [rt [12] []] (sp [-2, 2]) ++
      mulL [rt [-41472, -266112, -747936, -1209240, -1233792, -820992, -355680, -96696, -14976, -1008] [(-3, 3), (-2, 3), (-1, 3), (0, 3)]] (sp [-1]) ++
      mulL [rt [1728, 7920, 14916, 14256, 7260, 1872, 192] [(-3, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-48, -132, -96, -20] [(-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [-6] []] (sp [-1, -1, 2]) ++
      mulL [rt [144, 396, 288, 60] [(-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, 2]) ++
      mulL [rt [24] []] (sp [-1, 3]) ++
      mulL [rt [-1728, -7920, -14916, -14256, -7260, -1872, -192] [(-3, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [2]) ++
      mulL [rt [3] []] (sp [2, 2]) ++
      mulL [rt [-288, -792, -576, -120] [(-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [3]) ++
      mulL [rt [-42] []] (sp [4])
  | 4, 1 => mulL [rt [0, 0, 0, 0, 0, -874608, -3357888, -5722896, -5679552, -3613824, -1527048, -427920, -76560, -7920, -360] [(-6, 1), (-4, 1), (-3, 4), (-2, 4), (-1, 4), (0, 4)]] (sp []) ++
      mulL [rt [0, 0, 784, 760, 240, 24] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-3]) ++
      mulL [rt [8] []] (sp [-3, -1]) ++
      mulL [rt [0, 0, 0, 9480, 18660, 14424, 5460, 1008, 72] [(-6, 1), (-4, 1), (-3, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [-2]) ++
      mulL [rt [1728, 5472, 4920, 1788, 264, 12] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-6] []] (sp [-2, -1, -1]) ++
      mulL [rt [-6] []] (sp [-2, 2]) ++
      mulL [rt [-248832, -1700352, -5163264, -9191808, -10513008, -8021088, -4132080, -1425672, -319056, -43152, -3024, -72] [(-6, 1), (-4, 1), (-3, 3), (-2, 3), (-1, 3), (0, 3)]] (sp [-1]) ++
      mulL [rt [20736, 103680, 219456, 244872, 156528, 58680, 12564, 1392, 60] [(-6, 1), (-4, 1), (-3, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-864, -2736, -2656, -1084, -192, -12] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [-864, -2736, -3048, -1464, -312, -24] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, 2]) ++
      mulL [rt [20736, 103680, 219456, 254352, 175188, 73104, 18024, 2400, 132] [(-6, 1), (-4, 1), (-3, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [2]) ++
      mulL [rt [-3] []] (sp [2, 2]) ++
      mulL [rt [1728, 5472, 6096, 2928, 624, 48] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [3]) ++
      mulL [rt [6] []] (sp [4])
  | 4, 2 => mulL [rt [0, 0, 0, 0, 0, -16497792, -66875424, -121865760, -131653880, -93582976, -45897000, -15822272, -3826336, -634560, -68520, -4320, -120] [(-6, 2), (-4, 2), (-3, 4), (-2, 4), (-1, 4), (0, 4)]] (sp []) ++
      mulL [rt [-331776, -1797120, -4216320, -5477760, -4366752, -2246688, -760568, -168256, -23400, -1856, -64] [(-6, 2), (-4, 2), (-3, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [-2]) ++
      mulL [rt [4] []] (sp [-2, -2]) ++
      mulL [rt [1152, 3648, 3280, 1192, 176, 8] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-4] []] (sp [-2, -1, -1]) ++
      mulL [rt [4] []] (sp [-2, 2]) ++
      mulL [rt [0, 0, 0, 0, 2759040, 8906592, 12596832, 10233736, 5256320, 1767456, 388384, 53672, 4224, 144] [(-6, 2), (-4, 2), (-3, 3), (-2, 3), (-1, 3), (0, 3)]] (sp [-1]) ++
      mulL [rt [165888, 898560, 2108160, 2587200, 1821616, 761840, 184324, 21984, 188, -224, -16] [(-6, 2), (-4, 2), (-3, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-576, -1824, -1640, -596, -88, -4] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [2] []] (sp [-1, -1, 2]) ++
      mulL [rt [-576, -1824, -1640, -596, -88, -4] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, 2]) ++
      mulL [rt [-8] []] (sp [-1, 3]) ++
      mulL [rt [-165888, -898560, -2108160, -2890560, -2545136, -1484848, -576244, -146272, -23212, -2080, -80] [(-6, 2), (-4, 2), (-3, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [2]) ++
      mulL [rt [3] []] (sp [2, 2]) ++
      mulL [rt [1152, 3648, 3280, 1192, 176, 8] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [3]) ++
      mulL [rt [6] []] (sp [4])
  | 4, 3 => mulL [rt [0, 0, 0, 0, 0, -302220288, -1272282624, -2419374912, -2744497728, -2065241424, -1083724800, -405238896, -108141024, -20192256, -2499336, -177456, -3312, 432, 24] [(-6, 3), (-4, 3), (-3, 4), (-2, 4), (-1, 4), (0, 4)]] (sp []) ++
      mulL [rt [-96, -264, -192, -40] [(-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-3]) ++
      mulL [rt [8] []] (sp [-3, -1]) ++
      mulL [rt [0, 0, 0, 9480, 18660, 14424, 5460, 1008, 72] [(-6, 1), (-4, 1), (-3, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [-2]) ++
      mulL [rt [1728, 5472, 4920, 1788, 264, 12] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-6] []] (sp [-2, -1, -1]) ++
      mulL [rt [-6] []] (sp [-2, 2]) ++
      mulL [rt [0, 0, 0, 0, 51452928, 176196096, 268061376, 238895424, 138385584, 54606144, 14934768, 2819160, 357648, 28752, 1296, 24] [(-6, 3), (-4, 3), (-3, 3), (-2, 3), (-1, 3), (0, 3)]] (sp [-1]) ++
      mulL [rt [0, 0, 0, -186048, -417360, -390696, -197952, -58560, -10092, -936, -36] [(-6, 2), (-4, 2), (-3, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [-1, -1]) ++
      mulL [rt [-288, -912, -624, -108, 16, 4] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [864, 2736, 3048, 1464, 312, 24] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, 2]) ++
      mulL [rt [0, 0, 0, 41472, 125280, 151560, 95988, 34656, 7176, 792, 36] [(-6, 2), (-4, 2), (-3, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [2]) ++
      mulL [rt [-3] []] (sp [2, 2]) ++
      mulL [rt [-1728, -5472, -6096, -2928, -624, -48] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [3]) ++
      mulL [rt [6] []] (sp [4])
  | 4, 4 => mulL [rt [0, 0, 0, 0, 0, -5342257152, -22886074368, -44006588928, -49853814144, -36539831040, -17685594816, -5292242976, -602182776, 252949632, 156303672, 44243232, 7959648, 955776, 74712, 3456, 72] [(-6, 4), (-4, 4), (-3, 4), (-2, 4), (-1, 4), (0, 4)]] (sp []) ++
      mulL [rt [-48] []] (sp [-4]) ++
      mulL [rt [0, 0, 3136, 3040, 960, 96] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-3]) ++
      mulL [rt [32] []] (sp [-3, -1]) ++
      mulL [rt [0, 0, 0, 289152, 584160, 478272, 203928, 47808, 5832, 288] [(-6, 2), (-4, 2), (-3, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [-2]) ++
      mulL [rt [12] []] (sp [-2, -2]) ++
      mulL [rt [0, 0, -2352, -2280, -720, -72] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-2, -1]) ++
      mulL [rt [-12] []] (sp [-2, -1, -1]) ++
      mulL [rt [12] []] (sp [-2, 2]) ++
      mulL [rt [0, 0, 0, 0, 38513664, 120333312, 162611904, 123726816, 57363792, 15983976, 2182848, -101856, -97968, -17400, -1440, -48] [(-6, 3), (-4, 3), (-3, 3), (-2, 3), (-1, 3), (0, 3)]] (sp [-1]) ++
      mulL [rt [0, 0, 0, -144576, -292080, -239136, -101964, -23904, -2916, -144] [(-6, 2), (-4, 2), (-3, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [-1, -1]) ++
      mulL [rt [0, 0, 392, 380, 120, 12] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, -1, -1]) ++
      mulL [rt [1] []] (sp [-1, -1, -1, -1]) ++
      mulL [rt [-6] []] (sp [-1, -1, 2]) ++
      mulL [rt [0, 0, -1176, -1140, -360, -36] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [-1, 2]) ++
      mulL [rt [24] []] (sp [-1, 3]) ++
      mulL [rt [0, 0, 0, 144576, 292080, 239136, 101964, 23904, 2916, 144] [(-6, 2), (-4, 2), (-3, 2), (-2, 2), (-1, 2), (0, 2)]] (sp [2]) ++
      mulL [rt [3] []] (sp [2, 2]) ++
      mulL [rt [0, 0, 2352, 2280, 720, 72] [(-6, 1), (-4, 1), (-3, 1), (-2, 1), (-1, 1), (0, 1)]] (sp [3]) ++
      mulL [rt [-42] []] (sp [4])
  | _, _ => []

lemma idx_rt_sp (cs : List ℚ) (ps : List (ℤ × ℕ)) (m : List ℤ) (hm : ∀ a ∈ m, a ≠ 0) :
    ∀ τ ∈ mulL [rt cs ps] (sp m), S03.IsIndex τ.2.2 :=
  idx_mulL (by intro τ h; simp at h; subst h; rw [rt_w]; intro a ha; simp at ha) (sp_idx m hm)

lemma ST_mulL_w0 (τ : S03.DTerm) (hτ : τ.2.2 = []) (L : List S03.DTerm) (N : ℕ) :
    S03.ST (mulL [τ] L) N = S03.ratEval τ.1 τ.2.1 (N : ℚ) * S03.ST L N := by
  unfold mulL
  rw [List.flatMap_singleton, hτ]
  induction L with
  | nil => simp [S03.ST]
  | cons τ' L ih =>
    rw [List.flatMap_cons, ST_append, ih, stuf_nil_left, List.map_singleton, ST_cons, ST_cons, ST_nil,
      S08.PA.ratEval_mul]
    ring

set_option maxHeartbeats 10000000 in
lemma E1_val (j i n : ℕ) (hj : j ≤ 4) (hi : i ≤ 4) (hn : 1 ≤ n) :
    S03.ST (E1 j i) (2 * n) = S03.ST (Dz1 j i) (2 * n) := by
  have hn' : (1 : ℚ) ≤ n := by exact_mod_cast hn
  have h0 : (n : ℚ) ≠ 0 := by positivity
  have h2 : (2 * n : ℕ) ≠ 0 := by omega
  interval_cases j <;> interval_cases i <;>
  · simp only [E1, Dz1, ST_append, ST_map_scl, ST_m2, ST_m3, ST_m4]
    simp only [ST_rt_sp, List.mem_cons, List.not_mem_nil, ne_eq, forall_eq_or_imp, IsEmpty.forall_iff,
      implies_true, and_true, List.forall_mem_cons, not_false_eq_true, reduceCtorEq, Int.reduceNeg,
      Int.reduceEq, OfNat.ofNat_ne_zero, one_ne_zero]
    simp [atm, ST_cons, ST_nil, rt_eval, horner, ratEval_scl, ratEval_one, scl_w, rt_w, Z_nil]
    try field_simp
    ring

lemma E1_ST (j i n : ℕ) (hj : j ≤ 4) (hi : i ≤ 4) (hn : 1 ≤ n) :
    S03.ST (mulL [fz j i] (E1 j i)) (2 * n) = S03.ST (mulL [fz j i] (Dz1 j i)) (2 * n) := by
  rw [ST_mulL_w0 _ (fz_w j i), ST_mulL_w0 _ (fz_w j i),
    E1_val j i n hj hi hn]

/-- `T_z` for `𝔷 = 1`, `k = 4`, in single-letter form -/
noncomputable def Tz1E : List S03.DTerm :=
  (List.range 5).flatMap fun j => (List.range 5).flatMap fun i => mulL [fz j i] (E1 j i)

lemma Tz1E_ok (n : ℕ) (hn : 1 ≤ n) : (Q04.Pf 4 n).coeff 0 = S03.ST Tz1E (2 * n) := by
  rw [Tz1_ok n hn]
  unfold Tz1 Tz1E
  rw [ST_flatMap, ST_flatMap, list_range_sum, list_range_sum]
  apply Finset.sum_congr rfl
  intro j hj
  rw [ST_flatMap, ST_flatMap, list_range_sum, list_range_sum]
  apply Finset.sum_congr rfl
  intro i hi
  exact (E1_ST j i n (by simp at hj; omega) (by simp at hi; omega) hn).symm

end P04

namespace P04

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_0 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [] N = Dval Tz1E [] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_0 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_1 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-4] N = Dval Tz1E [-4] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_1 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-4] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_2 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-3] N = Dval Tz1E [-3] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_2 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-3] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_3 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-3, -1] N = Dval Tz1E [-3, -1] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_3 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-3, -1] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_4 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-2] N = Dval Tz1E [-2] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_4 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-2] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_5 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-2, -2] N = Dval Tz1E [-2, -2] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_5 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-2, -2] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_6 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-2, -1] N = Dval Tz1E [-2, -1] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_6 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-2, -1] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_7 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-2, -1, -1] N = Dval Tz1E [-2, -1, -1] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_7 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-2, -1, -1] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_8 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-2, 2] N = Dval Tz1E [-2, 2] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_8 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-2, 2] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_9 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-1] N = Dval Tz1E [-1] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_9 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-1] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_10 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-1, -3] N = Dval Tz1E [-1, -3] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_10 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-1, -3] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_11 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-1, -2] N = Dval Tz1E [-1, -2] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_11 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-1, -2] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_12 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-1, -2, -1] N = Dval Tz1E [-1, -2, -1] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_12 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-1, -2, -1] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_13 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-1, -1] N = Dval Tz1E [-1, -1] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_13 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-1, -1] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_14 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-1, -1, -2] N = Dval Tz1E [-1, -1, -2] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_14 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-1, -1, -2] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_15 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-1, -1, -1] N = Dval Tz1E [-1, -1, -1] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_15 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-1, -1, -1] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_16 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-1, -1, -1, -1] N = Dval Tz1E [-1, -1, -1, -1] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_16 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-1, -1, -1, -1] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_17 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-1, -1, 2] N = Dval Tz1E [-1, -1, 2] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_17 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-1, -1, 2] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_18 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-1, 2] N = Dval Tz1E [-1, 2] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_18 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-1, 2] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_19 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-1, 2, -1] N = Dval Tz1E [-1, 2, -1] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_19 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-1, 2, -1] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_20 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [-1, 3] N = Dval Tz1E [-1, 3] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_20 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [-1, 3] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_21 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [2] N = Dval Tz1E [2] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_21 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [2] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_22 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [2, -2] N = Dval Tz1E [2, -2] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_22 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [2, -2] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_23 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [2, -1] N = Dval Tz1E [2, -1] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_23 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [2, -1] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_24 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [2, -1, -1] N = Dval Tz1E [2, -1, -1] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_24 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [2, -1, -1] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_25 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [2, 2] N = Dval Tz1E [2, 2] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_25 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [2, 2] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_26 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [3] N = Dval Tz1E [3] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_26 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [3] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_27 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [3, -1] N = Dval Tz1E [3, -1] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_27 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [3, -1] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w0_28 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 0) [4] N = Dval Tz1E [4] N := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4w1_28 (N : ℚ) (hN : 1 < N) : Dval (Lι R4T 1) [4] N = 0 := by
  have h0 : N ≠ 0 := by positivity
  have h1 : N - 1 ≠ 0 := by linarith
  have g1 : N + 1 ≠ 0 := by linarith
  have g1' : 1 + N ≠ 0 := by linarith
  have g1'' : N - -1 ≠ 0 := by linarith
  have g2 : N + 2 ≠ 0 := by linarith
  have g2' : 2 + N ≠ 0 := by linarith
  have g2'' : N - -2 ≠ 0 := by linarith
  have g3 : N + 3 ≠ 0 := by linarith
  have g3' : 3 + N ≠ 0 := by linarith
  have g3'' : N - -3 ≠ 0 := by linarith
  have g4 : N + 4 ≠ 0 := by linarith
  have g4' : 4 + N ≠ 0 := by linarith
  have g4'' : N - -4 ≠ 0 := by linarith
  have g5 : N + 5 ≠ 0 := by linarith
  have g5' : 5 + N ≠ 0 := by linarith
  have g5'' : N - -5 ≠ 0 := by linarith
  have g6 : N + 6 ≠ 0 := by linarith
  have g6' : 6 + N ≠ 0 := by linarith
  have g6'' : N - -6 ≠ 0 := by linarith
  have g7 : N + 7 ≠ 0 := by linarith
  have g7' : 7 + N ≠ 0 := by linarith
  have g7'' : N - -7 ≠ 0 := by linarith
  have g8 : N + 8 ≠ 0 := by linarith
  have g8' : 8 + N ≠ 0 := by linarith
  have g8'' : N - -8 ≠ 0 := by linarith
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  field_simp
  ring

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4o_L0 (u : List ℤ) (N : ℚ) (h0 : ¬ ([] = u)) (h1 : ¬ ([-4] = u)) (h2 : ¬ ([-3] = u)) (h3 : ¬ ([-3, -1] = u)) (h4 : ¬ ([-2] = u)) (h5 : ¬ ([-2, -2] = u)) (h6 : ¬ ([-2, -1] = u)) (h7 : ¬ ([-2, -1, -1] = u)) (h8 : ¬ ([-2, 2] = u)) (h9 : ¬ ([-1] = u)) (h10 : ¬ ([-1, -3] = u)) (h11 : ¬ ([-1, -2] = u)) (h12 : ¬ ([-1, -2, -1] = u)) (h13 : ¬ ([-1, -1] = u)) (h14 : ¬ ([-1, -1, -2] = u)) (h15 : ¬ ([-1, -1, -1] = u)) (h16 : ¬ ([-1, -1, -1, -1] = u)) (h17 : ¬ ([-1, -1, 2] = u)) (h18 : ¬ ([-1, 2] = u)) (h19 : ¬ ([-1, 2, -1] = u)) (h20 : ¬ ([-1, 3] = u)) (h21 : ¬ ([2] = u)) (h22 : ¬ ([2, -2] = u)) (h23 : ¬ ([2, -1] = u)) (h24 : ¬ ([2, -1, -1] = u)) (h25 : ¬ ([2, 2] = u)) (h26 : ¬ ([3] = u)) (h27 : ¬ ([3, -1] = u)) (h28 : ¬ ([4] = u)) : Dval (Lι R4T 0) u N = 0 := by
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul, h0, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14, h15, h16, h17, h18, h19, h20, h21, h22, h23, h24, h25, h26, h27, h28]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  try simp

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4o_L1 (u : List ℤ) (N : ℚ) (h0 : ¬ ([] = u)) (h1 : ¬ ([-4] = u)) (h2 : ¬ ([-3] = u)) (h3 : ¬ ([-3, -1] = u)) (h4 : ¬ ([-2] = u)) (h5 : ¬ ([-2, -2] = u)) (h6 : ¬ ([-2, -1] = u)) (h7 : ¬ ([-2, -1, -1] = u)) (h8 : ¬ ([-2, 2] = u)) (h9 : ¬ ([-1] = u)) (h10 : ¬ ([-1, -3] = u)) (h11 : ¬ ([-1, -2] = u)) (h12 : ¬ ([-1, -2, -1] = u)) (h13 : ¬ ([-1, -1] = u)) (h14 : ¬ ([-1, -1, -2] = u)) (h15 : ¬ ([-1, -1, -1] = u)) (h16 : ¬ ([-1, -1, -1, -1] = u)) (h17 : ¬ ([-1, -1, 2] = u)) (h18 : ¬ ([-1, 2] = u)) (h19 : ¬ ([-1, 2, -1] = u)) (h20 : ¬ ([-1, 3] = u)) (h21 : ¬ ([2] = u)) (h22 : ¬ ([2, -2] = u)) (h23 : ¬ ([2, -1] = u)) (h24 : ¬ ([2, -1, -1] = u)) (h25 : ¬ ([2, 2] = u)) (h26 : ¬ ([3] = u)) (h27 : ¬ ([3, -1] = u)) (h28 : ¬ ([4] = u)) : Dval (Lι R4T 1) u N = 0 := by
  simp only [Lι, R4T, hW, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_map_scl,
    S03.dd, bar, List.map_cons, List.map_nil, Int.reduceNeg, Int.reduceSub, ite_true,
    ite_false, Int.reduceEq, Int.reduceMul, Nat.cast_one, Int.reduceAdd, Nat.cast_ofNat]
  simp only [Sh_nil', Sh_0, Sh_1, Sh_2, Sh_3, Sh_m1, Sh_m2, Dval_append, Dval_map_scl, List.map_cons,
    List.map_nil, List.cons_append, List.nil_append]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul, h0, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14, h15, h16, h17, h18, h19, h20, h21, h22, h23, h24, h25, h26, h27, h28]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  try simp

set_option maxHeartbeats 40000000 in
set_option maxRecDepth 100000 in
lemma R4o_T (u : List ℤ) (N : ℚ) (h0 : ¬ ([] = u)) (h1 : ¬ ([-4] = u)) (h2 : ¬ ([-3] = u)) (h3 : ¬ ([-3, -1] = u)) (h4 : ¬ ([-2] = u)) (h5 : ¬ ([-2, -2] = u)) (h6 : ¬ ([-2, -1] = u)) (h7 : ¬ ([-2, -1, -1] = u)) (h8 : ¬ ([-2, 2] = u)) (h9 : ¬ ([-1] = u)) (h10 : ¬ ([-1, -3] = u)) (h11 : ¬ ([-1, -2] = u)) (h12 : ¬ ([-1, -2, -1] = u)) (h13 : ¬ ([-1, -1] = u)) (h14 : ¬ ([-1, -1, -2] = u)) (h15 : ¬ ([-1, -1, -1] = u)) (h16 : ¬ ([-1, -1, -1, -1] = u)) (h17 : ¬ ([-1, -1, 2] = u)) (h18 : ¬ ([-1, 2] = u)) (h19 : ¬ ([-1, 2, -1] = u)) (h20 : ¬ ([-1, 3] = u)) (h21 : ¬ ([2] = u)) (h22 : ¬ ([2, -2] = u)) (h23 : ¬ ([2, -1] = u)) (h24 : ¬ ([2, -1, -1] = u)) (h25 : ¬ ([2, 2] = u)) (h26 : ¬ ([3] = u)) (h27 : ¬ ([3, -1] = u)) (h28 : ¬ ([4] = u)) : Dval Tz1E u N = 0 := by
  simp only [Tz1E, range5, List.flatMap_cons, List.flatMap_nil, Dval_append, Dval_nil, Dval_mulL_fz]
  simp only [E1, Dval_append, Dval_mulL_rt]
  simp only [sp, mulL, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.cons_append, List.nil_append, List.flatMap_singleton]
  simp only [stuf_cons, stuf_nil_left, stuf_nil_right, List.map_cons, List.map_nil, List.sum_cons,
    List.sum_nil, opl_eval, List.cons_append, List.nil_append, List.map_append, List.sum_append, List.flatMap_cons,
    List.flatMap_nil, List.append_nil]
  norm_num only [Int.sign_neg, Int.sign_one, Int.natAbs_neg, Int.natAbs_one, Int.reduceNeg,
    List.cons.injEq, ratEval_one, rt_eval, horner, scl_w, ratEval_scl, reduceCtorEq, ite_true,
    ite_false, true_and, and_true, false_and, and_false, List.prod_cons, List.prod_nil, List.map_cons,
    List.map_nil, fz, Dval_cons, Dval_nil, sgnP, nu, List.filter_cons, List.filter_nil, List.length_cons,
    List.length_nil, ratEval_10, Dval_append, Dval_map_scl, sgn_pos_lit, decide_true, decide_false,
    Int.sign_zero, Int.reduceAbs, rt_w, S08.PA.ratEval_mul, h0, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14, h15, h16, h17, h18, h19, h20, h21, h22, h23, h24, h25, h26, h27, h28]
  try simp only [S03.ratEval, Polynomial.eval_one, Polynomial.eval_mul, add_zero, zero_add, Finsupp.prod_zero_index, div_one, mul_one, one_mul]
  try simp

lemma Deq_R4_0 : Deq (Lι R4T 0) Tz1E := by
  intro u N hN
  by_cases e0 : [] = u
  · subst e0; exact R4w0_0 N hN
  by_cases e1 : [-4] = u
  · subst e1; exact R4w0_1 N hN
  by_cases e2 : [-3] = u
  · subst e2; exact R4w0_2 N hN
  by_cases e3 : [-3, -1] = u
  · subst e3; exact R4w0_3 N hN
  by_cases e4 : [-2] = u
  · subst e4; exact R4w0_4 N hN
  by_cases e5 : [-2, -2] = u
  · subst e5; exact R4w0_5 N hN
  by_cases e6 : [-2, -1] = u
  · subst e6; exact R4w0_6 N hN
  by_cases e7 : [-2, -1, -1] = u
  · subst e7; exact R4w0_7 N hN
  by_cases e8 : [-2, 2] = u
  · subst e8; exact R4w0_8 N hN
  by_cases e9 : [-1] = u
  · subst e9; exact R4w0_9 N hN
  by_cases e10 : [-1, -3] = u
  · subst e10; exact R4w0_10 N hN
  by_cases e11 : [-1, -2] = u
  · subst e11; exact R4w0_11 N hN
  by_cases e12 : [-1, -2, -1] = u
  · subst e12; exact R4w0_12 N hN
  by_cases e13 : [-1, -1] = u
  · subst e13; exact R4w0_13 N hN
  by_cases e14 : [-1, -1, -2] = u
  · subst e14; exact R4w0_14 N hN
  by_cases e15 : [-1, -1, -1] = u
  · subst e15; exact R4w0_15 N hN
  by_cases e16 : [-1, -1, -1, -1] = u
  · subst e16; exact R4w0_16 N hN
  by_cases e17 : [-1, -1, 2] = u
  · subst e17; exact R4w0_17 N hN
  by_cases e18 : [-1, 2] = u
  · subst e18; exact R4w0_18 N hN
  by_cases e19 : [-1, 2, -1] = u
  · subst e19; exact R4w0_19 N hN
  by_cases e20 : [-1, 3] = u
  · subst e20; exact R4w0_20 N hN
  by_cases e21 : [2] = u
  · subst e21; exact R4w0_21 N hN
  by_cases e22 : [2, -2] = u
  · subst e22; exact R4w0_22 N hN
  by_cases e23 : [2, -1] = u
  · subst e23; exact R4w0_23 N hN
  by_cases e24 : [2, -1, -1] = u
  · subst e24; exact R4w0_24 N hN
  by_cases e25 : [2, 2] = u
  · subst e25; exact R4w0_25 N hN
  by_cases e26 : [3] = u
  · subst e26; exact R4w0_26 N hN
  by_cases e27 : [3, -1] = u
  · subst e27; exact R4w0_27 N hN
  by_cases e28 : [4] = u
  · subst e28; exact R4w0_28 N hN
  rw [R4o_L0 u N e0 e1 e2 e3 e4 e5 e6 e7 e8 e9 e10 e11 e12 e13 e14 e15 e16 e17 e18 e19 e20 e21 e22 e23 e24 e25 e26 e27 e28, R4o_T u N e0 e1 e2 e3 e4 e5 e6 e7 e8 e9 e10 e11 e12 e13 e14 e15 e16 e17 e18 e19 e20 e21 e22 e23 e24 e25 e26 e27 e28]

lemma Deq_R4_1 : Deq (Lι R4T 1) [] := by
  intro u N hN
  by_cases e0 : [] = u
  · subst e0; rw [R4w1_0 N hN, Dval_nil]
  by_cases e1 : [-4] = u
  · subst e1; rw [R4w1_1 N hN, Dval_nil]
  by_cases e2 : [-3] = u
  · subst e2; rw [R4w1_2 N hN, Dval_nil]
  by_cases e3 : [-3, -1] = u
  · subst e3; rw [R4w1_3 N hN, Dval_nil]
  by_cases e4 : [-2] = u
  · subst e4; rw [R4w1_4 N hN, Dval_nil]
  by_cases e5 : [-2, -2] = u
  · subst e5; rw [R4w1_5 N hN, Dval_nil]
  by_cases e6 : [-2, -1] = u
  · subst e6; rw [R4w1_6 N hN, Dval_nil]
  by_cases e7 : [-2, -1, -1] = u
  · subst e7; rw [R4w1_7 N hN, Dval_nil]
  by_cases e8 : [-2, 2] = u
  · subst e8; rw [R4w1_8 N hN, Dval_nil]
  by_cases e9 : [-1] = u
  · subst e9; rw [R4w1_9 N hN, Dval_nil]
  by_cases e10 : [-1, -3] = u
  · subst e10; rw [R4w1_10 N hN, Dval_nil]
  by_cases e11 : [-1, -2] = u
  · subst e11; rw [R4w1_11 N hN, Dval_nil]
  by_cases e12 : [-1, -2, -1] = u
  · subst e12; rw [R4w1_12 N hN, Dval_nil]
  by_cases e13 : [-1, -1] = u
  · subst e13; rw [R4w1_13 N hN, Dval_nil]
  by_cases e14 : [-1, -1, -2] = u
  · subst e14; rw [R4w1_14 N hN, Dval_nil]
  by_cases e15 : [-1, -1, -1] = u
  · subst e15; rw [R4w1_15 N hN, Dval_nil]
  by_cases e16 : [-1, -1, -1, -1] = u
  · subst e16; rw [R4w1_16 N hN, Dval_nil]
  by_cases e17 : [-1, -1, 2] = u
  · subst e17; rw [R4w1_17 N hN, Dval_nil]
  by_cases e18 : [-1, 2] = u
  · subst e18; rw [R4w1_18 N hN, Dval_nil]
  by_cases e19 : [-1, 2, -1] = u
  · subst e19; rw [R4w1_19 N hN, Dval_nil]
  by_cases e20 : [-1, 3] = u
  · subst e20; rw [R4w1_20 N hN, Dval_nil]
  by_cases e21 : [2] = u
  · subst e21; rw [R4w1_21 N hN, Dval_nil]
  by_cases e22 : [2, -2] = u
  · subst e22; rw [R4w1_22 N hN, Dval_nil]
  by_cases e23 : [2, -1] = u
  · subst e23; rw [R4w1_23 N hN, Dval_nil]
  by_cases e24 : [2, -1, -1] = u
  · subst e24; rw [R4w1_24 N hN, Dval_nil]
  by_cases e25 : [2, 2] = u
  · subst e25; rw [R4w1_25 N hN, Dval_nil]
  by_cases e26 : [3] = u
  · subst e26; rw [R4w1_26 N hN, Dval_nil]
  by_cases e27 : [3, -1] = u
  · subst e27; rw [R4w1_27 N hN, Dval_nil]
  by_cases e28 : [4] = u
  · subst e28; rw [R4w1_28 N hN, Dval_nil]
  rw [R4o_L1 u N e0 e1 e2 e3 e4 e5 e6 e7 e8 e9 e10 e11 e12 e13 e14 e15 e16 e17 e18 e19 e20 e21 e22 e23 e24 e25 e26 e27 e28, Dval_nil]

end P04

namespace P04

lemma R4C_ok : ∀ t ∈ ([] : List (ℚ × ℤ)), t.2 ≤ 1 := by simp

lemma R4_even (n : ℕ) (hn : 1 ≤ n) :
    coefL R4T [] (2 * (n : ℤ)) = (Q04.Pf 4 n).coeff 0 := by
  rw [show 2 * (n : ℤ) = 2 * (n : ℤ) + 0 by ring, coefL_eq R4T [] R4T_ok R4C_ok 0 le_rfl n hn,
    ST_of_Deq _ _ Deq_R4_0 _ (by omega), ← Tz1E_ok n hn]

lemma R4_odd (n : ℕ) (hn : 1 ≤ n) : coefL R4T [] (2 * (n : ℤ) + 1) = 0 := by
  rw [coefL_eq R4T [] R4T_ok R4C_ok 1 (by norm_num) n hn,
    ST_of_Deq _ _ Deq_R4_1 _ (by omega), ST_nil]

end P04

/-- M8 for `R_4`: (i), (ii) from `coefL_eq`, the element `Tz1E` (`Tz1E_ok`, via `Tz1_ok` and
`gH_four`) and the identities `Deq_R4_0`, `Deq_R4_1`; (iii) from `R4_low`. -/
theorem P04.proof_6 : P04.stmt_6 :=
  ⟨P04.R4_even, P04.R4_odd, P04.R4_low⟩

theorem P04.proof_7 : P04.stmt_7 := fun c hc0 hc1 =>
  P04.series_of_coef P04.R4T [] P04.R4T_ok _ P04.proof_6.1 P04.proof_6.2.1
    P04.proof_6.2.2 c hc0 hc1

theorem P04.proof_8 : P04.stmt_8 := by
  intro α hα0 hα1
  obtain ⟨hA, hM⟩ := S04.proof_7 α hα0 hα1
  obtain ⟨hc0, hc1⟩ := S08.PD.c_pos α hα0 hα1
  have e1 : R01.zM 0 4 α = P04.R4 (D04.c α) := by
    have := R01.zM_nat 0 4 α
    simp only [Nat.cast_ofNat] at this
    rw [this]
    exact (P04.proof_7 _ hc0 hc1).tsum_eq
  have e2 : R01.zM (Finsupp.single 2 1) 4 α = P04.Psi2 (D04.c α) := by
    have := R01.zM_nat (Finsupp.single 2 1) 4 α
    simp only [Nat.cast_ofNat] at this
    rw [this]
    exact (P04.proof_5 _ hc0 hc1).tsum_eq
  refine ⟨hA, ?_⟩
  rw [hM, e1, e2]
