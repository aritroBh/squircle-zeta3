import Mathlib
import Zeta3.Wave1.T12_EdgeLaw
import Zeta3.Wave2.D00_Dilog
import Zeta3.Wave4.X04_PolyFirstDeriv
import Zeta3.Wave6.Z03_NonlocalValues
/-!
# S05 -- `P'''(1)` of every hexagon (thm:hexagon)

Target id: S05.
Paper label: `thm:hexagon` (main.tex 1709-1725, Theorem "P'''(1) of every hexagon").
Detailed proof: paper/ieee/proofs/sec_g4_hexagon.tex (Theorem 25 = `g4:thm:hexagon`, lines 8-19;
"What is proved" 21-22; "One edge" 26-34; Lemma `g4:lem:precise` 39-52; Lemma `g4:lem:euler`
70-81; Lemma `g4:lem:Jseries` 83-93; "A worked entry" 99-107; Lemma `g4:lem:NK` 116-126;
Lemma `g4:lem:assembly` 128-138; "Proof of the five sentences" 143-145).
Related: `prop:edgelaw` (main.tex 1671-1684; Lean T12), `prop:nonlocal` (main.tex 1697-1705;
Lean Z03), the setting of `thm:poly` (main.tex 1645, 1669; Lean X04 for the regular case).

## Quoted verbatim

main.tex 1709 (setting): "We use the weighted convention ... $K_p=\bigl\{\sum_{j=1}^3w_j|k_j
\langle v_j,x\rangle|^p\le1\bigr\}$, where the weights $w_j>0$ and the factors $k_j>0$ do not
change with $p$ and are chosen so that $K_1$ is the given hexagon ... On an edge $e$ keep the
notation of Proposition~\ref{prop:nonlocal}, with the three terms of
$\sum_jw_j|k_j\langle v_j,x\rangle|$ called $A_0t$, $A_1(1-t)$ and $B(1+\kappa t)$."
main.tex 1695: "$y(t)=V_0+t(V_1-V_0)$, $L$, $h$ as in Proposition~\ref{prop:edgelaw},
$q_0=V_0\cdot(V_1-V_0)=-LL_0$ ... so $A_1+B=1$, $A_0=A_1-B\kappa$, and $r=1+\kappa$."

main.tex 1711-1718 (Theorem thm:hexagon):
  For every non-degenerate centrally symmetric convex hexagon and all weights $w_j>0$,
  $P'''(1)=2\sum_eE_e$, summed over three edges, one from each pair of opposite edges. With
  $x=1-r$ and $\ell=-\log r=-\log(1-x)$, each $E_e$ is an explicit linear combination of
  \Li_3(x),\ \Li_3\bigl(\tfrac x{x-1}\bigr),\ \Li_2(x)\,\ell,\ \ell^3,\ \Li_2(x),\ \ell^2,
  \ell,\ \zeta(3),\ \zeta(2)\,\ell,\ \zeta(2),\ 1,
  whose coefficients are polynomials in $A_1$, $B$, $r$, $q_0$, $|V_0|^2$, $L^2$ and the
  logarithms of the three corner term values $A_0$, $A_1$, $B$ (each divided by its weight),
  divided by $L^5$. Written this way, the coefficient of $\zeta(3)$ is exactly $Z_e$, for every
  $r$ (Proposition~\ref{prop:edgelaw}). The non-local $\zeta(3)$ of
  Proposition~\ref{prop:nonlocal} appears at $r=2,\frac12$ through $\Li_3(1-r)$ and
  $\Li_3\bigl(\frac{r-1}r\bigr)$, with the coefficients $2\Phi_A+\Phi_L$ and $\Phi_A$. The
  theorem reduces to Proposition~\ref{prop:nonlocal} at $r=2,\frac12$ and, when $r=1$ on every
  edge and all $w_j=1$, to the hexagon law of Proposition~\ref{prop:edgelaw}.

Supplement 22 ("What is proved"): "``Non-degenerate'' means that $K_1$ has six vertices, i.e.
$v_1,v_2,v_3$ are pairwise non-parallel; then $0<B<1$, $A_0,A_1>0$, $r>0$.
$E_e=3!\int_0^1E_3\dd t$ is the contribution of the edge $e$ to $P'''(1)$. ``Written this way''
means: the eleven functions are the formal symbols, $\zeta(2),\zeta(3)$ are the numbers, and the
coefficient of $\zeta(3)$ is read off the formula as written. We do not claim that the eleven
functions are linearly independent ..."

Supplement 28 ("One edge"): "$a_{j_0}=A_0t$, $a_{j_1}=A_1(1-t)$, $a_{j_2}=B(1+\kappa t)$ with
$A_1=1-B$ ($\mathcal N(V_0)=1$) and $A_0+Br=1$ ($\mathcal N(V_1)=1$) ... Hence
$P'''(1)=6\sum_{\text{6 edges}}\int_0^1E_3\dd t=2\sum_{e}E_e$ over three edges, one from each
pair, with $E_e:=6\int_0^1E_3\dd t$."

Supplement 39-51 (Lemma g4:lem:precise):
  There are polynomials $c_1,\dots,c_{11}\in\ZZ[B,r,q_0,V_{02},L^2][\mathrm{l}_0,\mathrm{l}_1,
  \mathrm{l}_B]$ such that for every $r>0$, with $x=1-r$, $\ell=-\log r$,
  L^5E_e=c_1\Li_3(x)+c_2\Li_3\bigl(\tfrac x{x-1}\bigr)+c_3\Li_2(x)\ell+c_4\ell^3+c_5\Li_2(x)
  +c_6\ell^2+c_7\ell+c_8\zeta(3)+c_9\zeta(2)\ell+c_{10}\zeta(2)+c_{11}
  (for $r=1$ read the value $x=0$). The five blocks that matter for $\zeta(3)$, with
  $\Delta=q_0^2-L^2V_{02}=-L^2h^2$, are
  c_8&=-6(B-1)(Br-1)\,\Delta\,\bigl(L^2-3BL^2+3Bq_0(r-1)\bigr),
  c_1&=-6B\,\Delta\,\bigl(-3BL^2+3Bq_0(r-1)+(4-r)L^2-3q_0(r-1)\bigr),
  c_2&=-6B\,\Delta\,\bigl(-3BL^2r+3Bq_0r(r-1)+(r+2)L^2-3q_0(r-1)\bigr),
  c_3&=-6B(Br-1)\,\Delta\,\bigl(-2L^2+3q_0(r-1)\bigr),\qquad c_9=18B\,q_0\,r(B-1)(Br-1)\,\Delta ,
  and the other six blocks are listed in \path{scripts/g4_formula_table.txt}.
  (Supplement 54: "The printed $\ell$ is $-\log r$.")
  Supplement 30: "$\mathrm{l}_0=\log(A_0/w_{j_0})$, $\mathrm{l}_1=\log(A_1/w_{j_1})$,
  $\mathrm{l}_B=\log(B/w_{j_2})$".

Supplement 70-77 (Lemma g4:lem:euler): For $|x|<1$, with $H_n=\sum_{j\le n}\frac1j$,
  $H_n^{(2)}=\sum_{j\le n}\frac1{j^2}$:
  \sum_{n\ge1}H_nx^n=\frac{\mathcal L_1}{1-x},\quad \sum_{n\ge1}H^{(2)}_nx^n=
  \frac{\mathcal L_2}{1-x},\quad \sum_{n\ge1}H_n^2x^n=\frac{\mathcal L_1^2+\mathcal L_2}{1-x},
  \sum_{n\ge1}\frac{H_n}nx^n=\mathcal L_2+\tfrac12\mathcal L_1^2,\qquad
  \sum_{n\ge1}\frac{H_n^2}nx^n=\mathcal L_3+\mathcal L_1\mathcal L_2+\tfrac13\mathcal L_1^3,
  \sum_{n\ge1}\frac{H_n^{(2)}}nx^n=3\mathcal L_3-2\mathcal A+\mathcal L_1\mathcal L_2 .
  (Supplement 68: "$\mathcal L_1=-\log(1-x)$, $\mathcal L_2=\Li_2(x)$, $\mathcal L_3=\Li_3(x)$,
  and $\mathcal A=A(x)=\sum_{n\ge1}H_nx^n/n^2$".)

Supplement 64-66, 83-89 (Lemma g4:lem:Jseries): "for $|x|<1$ ...
  $J(n;a,b,c)=\int_0^1t^nX^aY^b\log^c(1-xt)\dd t$, $a+b+c\le3$" and
  J(n;a,b,c)=\sum_{q\ge1}x^q\,\lambda_c(q)\,\mathcal I_{n+q}(a,b),
  where $\log^c(1-z)=\sum_{q\ge1}\lambda_c(q)z^q$ and $\mathcal I_N(a,b)=\int_0^1s^N\log^as
  \log^b(1-s)\dd s$. Explicitly $\lambda_c(q)=c!\,[v^c]\Bigl(-\frac vq\exp\bigl(-\sum_{i\ge1}
  \frac{v^i}iH^{(i)}_{q-1}\bigr)\Bigr)$ ...

Supplement 100-107 (worked entries):
  J(0;1,0,1)=\mathcal L_1-\frac{\mathcal L_1-x}x-\frac{\mathcal L_2-x}x
   =\frac{-\mathcal L_1-\mathcal L_2+x(\mathcal L_1+2)}x .
  "Two other entries: $J(0;0,0,1)=-\mathcal L_1+\frac{\mathcal L_1}x-1$ ... and
  $J(1;0,0,1)=-\frac{\mathcal L_1}2+\frac{\mathcal L_1}{2x^2}-\frac14-\frac1{2x}$. The entries
  with $a+b+c=3$ contain $\mathcal A$, for instance"
  J(0;0,2,1)=\frac{1-x}{3x}\Bigl(12\mathcal L_3-6\mathcal A+6\mathcal L_1\mathcal L_2
   +6\mathcal L_2+\mathcal L_1^3+3\mathcal L_1^2+6\mathcal L_1\Bigr)-6 .

Supplement 141 (M6, the checks that sentences 3 and 4 rest on): "the coefficient of $\zeta(3)$
equals $Z_e=6A_0A_1\frac{h^2}L\bigl(3\frac{A_0L_0+A_1L_1}L-2\bigr)$ ... the coefficient of
$\Li_3(\frac x{x-1})$ equals $\Phi_A$ and that of $\Li_3(x)$ equals $2\Phi_A+\Phi_L$".
main.tex 1698-1703 (prop:nonlocal): "with $\gamma_e=Bh^2/L^3$,
  \Phi_A&=-6\gamma_e\bigl(3A_0\kappa q_0-L^2(2A_0+rA_1)\bigr),
  \Phi_L&=-18\gamma_e\Bigl(\kappa\bigl(-A_1^2+A_1B(2\kappa-1)+2B^2\kappa\bigr)q_0
   +L^2(A_0B+rA_1^2)\Bigr)."

Supplement 145 (fifth sentence, r = 1): "then $a_{j_2}$ is constant on every edge, so the
hexagon is fully local, hence affine-regular ... and $2\sum_eE_e=6\sum_{\text{edges}}\int E_3$ is
the quantity of Lemma~\ref{g4:lem:laws}(b)", i.e. the $d^3$ term of the hexagon law of
prop:edgelaw (main.tex 1680-1683):
  2^{\frac{d}{1+d}}\Bigl[P(1)+\frac{\pi^2\mathcal E}{24}d^2-\frac{\mathcal E}{8}
  \Bigl(\zeta(3)+\frac{2\pi^2}{3}\Bigr)d^3\Bigr]+O(d^4).

## Conventions (how the paper's objects are modelled)

* The plane is `X04.Plane` (Euclidean `ℝ²`) with `X04.dot`; `det2` is the 2×2 determinant.
* Weights and factors: only the products `u_j = k_j v_j` enter `K_p`, so a hexagon with weights
  is the data `u : Fin 3 → Plane`, `w : Fin 3 → ℝ`, and
  `F_p(x) = Σ_j w_j |⟨u_j,x⟩|^p`, `K_p = {F_p ≤ 1}` (this is the paper's
  `Σ_j w_j |k_j⟨v_j,x⟩|^p`). The term value is `a_j(x) = w_j|⟨u_j,x⟩|` (supplement:
  `a_j = w_jk_j|⟨v_j,x⟩|`), so "the corner term value divided by its weight" is `|⟨u_j,·⟩|`.
* `IsHex`: weights positive, `u_j` pairwise non-parallel (non-degenerate), and the corners
  `C_0, C_1, C_2` (the other three are `-C_j`) labelled so that `C_j` lies on the kink line
  `⟨u_j,x⟩=0` and on `F_1 = 1`, in counter-clockwise order `C_0, C_1, C_2, -C_0, ...`
  ("the weights are listed in the order of the corners, $w_j$ going with the kink line through
  the $j$-th corner"). stmt_5 says every such hexagon with every choice of weights arises, and
  stmt_6 that `K_1` is then the convex hull of the six corners.
* The three edges are `C_0 → C_1`, `C_1 → C_2`, `C_2 → -C_0` (one from each opposite pair); on
  edge `k` the kink lines are `j₀ = k` (through `V_0`), `j₁ = k+1` (through `V_1`), `j₂ = k+2`
  (in `Fin 3`). `A_0 = a_{j₀}(V_1)`, `A_1 = a_{j₁}(V_0)`, `B = a_{j₂}(V_0)`,
  `r = a_{j₂}(V_1)/B` (supplement line 9: "$B=a_{j_2}(V_0)$ ..., $r=a_{j_2}(V_1)/a_{j_2}(V_0)$").
* `P(p)` (`Perim`) is the Euclidean length (total variation, as in X04/W04) of the closed polar
  curve `θ ↦ F_p(e_θ)^{-1/p} e_θ`, `θ ∈ [0,2π]`, whose image is `{F_p = 1} = ∂K_p` (stmt_7). For
  `p < 1` this is the length of the same closed curve (the reading of checklist note 8).
* `E_e` (`Ee`) is the explicit edge integral `∫_0^1 ∂_p^3|x_p'(t)| |_{p=1} dt`, where
  `x_p(t) = F_p(y(t))^{-1/p} y(t)` is the arc of `∂K_p` over the edge. Since
  `|x_p'(t)| = G(t,d)` (Lemma g4:lem:arc) and `E_3` is the `d³` coefficient of `G`, this is the
  paper's `E_e = 3!∫_0^1E_3 dt`. Integrability of the integrand is stated explicitly (stmt_9),
  so Lean's convention "non-integrable ⇒ 0" plays no role.
* `ζ(3)` is `Z03.z3 = (riemannZeta 3).re` (Z03 stmt_1 relates it to `riemannZeta 3` and to
  `Σ 1/(n+1)³`); `ζ(2)` is written `π²/6` (Mathlib `riemannZeta_two`), as in Z03.
  `Li₂ = D00.Li2`, `Li₃ = Z03.Li3` (both genuine on `(-∞,1]`, which contains `x = 1-r` and
  `x/(x-1) = (r-1)/r` for every `r > 0`), `A = Z03.A` (the series, used for `|x| < 1` only).
* `H_n` is `harmonic n` (cast to `ℝ`), `H^{(2)}_n` is `H2 n`; series over `n ≥ 1` are indexed
  by `n+1`, `n : ℕ`.
* The six blocks `c_4, c_5, c_6, c_7, c_10, c_11` that the supplement does not print are
  existentially quantified integer polynomials in the eight edge variables
  `B, r, q_0, V_{02}, L², l_0, l_1, l_B` (in this order, `Fin 8`); the five printed blocks are
  written out.

## Already proved elsewhere (imported, not restated)

* Lemma g4:lem:NK (Nielsen–Kölbig with Landen) and the regrouping into `Φ_A A + Φ_L Li₃`: Z03
  stmt_10, stmt_11, stmt_12; `A' = (Li₂ + ½log²(1-x))/x`: Z03 stmt_8.
* The values of prop:nonlocal at `x = -1, ½`: Z03 stmt_13 to stmt_17.
* `ΣZ_e = -6𝓔/8` for affine-regular hexagons: T12 stmt_3, stmt_4.

## NOT COVERED

* The six coefficient blocks `c_4, c_5, c_6, c_7, c_10, c_11` as explicit polynomials: the
  supplement does not print them ("listed in scripts/g4_formula_table.txt", the output of
  machine step M6); they appear only existentially in stmt_10. (A numerical check of the full
  formula with the table's blocks was run on a random non-regular hexagon; see below.)
* Lemma g4:lem:assembly as a formula in the 80 monomial coefficients `c_{n,a,b,c}` of machine
  step M4 (not printed), Lemma g4:lem:struct, the 40 local Beta values (M2) and the 40 closed
  forms of `J` (M5, not printed; only the four worked entries, stmt_4, are stated). The
  continuation to every `r > 0` is part of stmt_10 ("for every `r > 0`").
* The remark that no closed form `J` contains `ζ(3)`: it is a statement about the shape of
  machine output and about which constants a reduction "introduces"; not stated.
* The `a + b ≤ 2` description of `𝓘_N(a,b)` in g4:lem:Jseries ("a polynomial in $H_N$,
  $H_N^{(2)}$, $\frac1{N+1}$ and $\fz_2$, divided by $N+1$"): not a precise formula; not stated.
  The general-`c` formula for `λ_c(q)` is stated only for `c = 1, 2, 3`, the values used
  (`c ≤ 3` since `a+b+c ≤ 3`).
* "The theorem reduces to Proposition 8 at r = 2, ½": prop:nonlocal speaks of "the ζ(3) of
  P'''(1) beyond Z_e", a basis-dependent coefficient; the needed facts are the Φ identities
  (stmt_12) together with Z03 stmt_12 to stmt_17. Not restated.
* "Written this way, the coefficient of ζ(3)" is formalised as "the printed block `c_8`" (the
  paper says it does not claim linear independence of the eleven functions); stmt_11 says
  `c_8 = L⁵ Z_e`.

## Numerical check before freezing

(mpmath, 30 digits, own code.) Random non-regular hexagon `V = (1.3,-0.2), (1.1,1.4),
(-0.5,1.2)`, weights `1.7, 0.6, 2.9` (edge ratios `r ≈ 1.384, 0.716, 1.010`): for every edge
`E_e` computed from the definition here (`∫∂_p³|x_p'|`, numerical differentiation) agrees with the
closed form of stmt_10 (printed blocks plus the table's six blocks) to `2·10⁻²⁶`; `2ΣE_e =
-21.05325812923…` agrees with 13-point finite differences of the polar perimeter to `10⁻¹¹`.
The paper's first example (corners `(1,0),(2,2),(0,1)`, weights `1,3,1`) gives `19.90026636…`
as printed. The printed blocks equal the table's (exactly, SymPy); `c_8 = L⁵Z_e`,
`c_2 = L⁵Φ_A`, `c_1 = L⁵(2Φ_A+Φ_L)` exactly; an affine-regular hexagon with `w = 1` gives
`P'''(1)` equal to the hexagon law to `10⁻²⁹`; the six Euler sums, `λ_1, λ_2, λ_3`, four
instances of the `J` series and the four worked entries (at `x = -0.63, 0.37`) to `10⁻²⁶`;
stmt_5 (positive `u`-scalings) on 2000 random convex hexagons with random weights.

Statements 10 and 14 moved to `Zeta3.Wave10.S07_HexagonFormula` (continuation of S05); this file keeps 1-9, 11, 12, 13.
-/

-- BEGIN STATEMENTS
namespace S05

open MeasureTheory Real Set

/-- The determinant of two plane vectors (nonzero iff they are not parallel). -/
def det2 (x y : X04.Plane) : ℝ := x 0 * y 1 - x 1 * y 0

/-- The weighted `F_p(x) = Σ_j w_j |⟨u_j, x⟩|^p` (the paper's `Σ_j w_j|k_j⟨v_j,x⟩|^p`,
`u_j = k_j v_j`). -/
noncomputable def Fw (u : Fin 3 → X04.Plane) (w : Fin 3 → ℝ) (p : ℝ) (x : X04.Plane) : ℝ :=
  ∑ j : Fin 3, w j * |X04.dot (u j) x| ^ p

/-- `K_p = {x : F_p(x) ≤ 1}`. -/
def Kw (u : Fin 3 → X04.Plane) (w : Fin 3 → ℝ) (p : ℝ) : Set X04.Plane := {x | Fw u w p x ≤ 1}

/-- The term value `a_j(x) = w_j |⟨u_j, x⟩|` (`= w_jk_j|⟨v_j,x⟩|`). -/
noncomputable def term (u : Fin 3 → X04.Plane) (w : Fin 3 → ℝ) (j : Fin 3) (x : X04.Plane) :
    ℝ :=
  w j * |X04.dot (u j) x|

/-- A weighted hexagon: kink normals `u`, weights `w`, corners `C` (the other three are `-C`). -/
structure Hex where
  u : Fin 3 → X04.Plane
  w : Fin 3 → ℝ
  C : Fin 3 → X04.Plane

/-- Non-degenerate weighted hexagon with labelled corners: `w_j > 0`; the `u_j` pairwise
non-parallel; `C_j` on the kink line `⟨u_j,x⟩ = 0` and on `F_1 = 1`; the corners
`C_0, C_1, C_2, -C_0, -C_1, -C_2` in counter-clockwise order. -/
def IsHex (H : Hex) : Prop :=
  (∀ j, 0 < H.w j) ∧ (∀ i j, i ≠ j → det2 (H.u i) (H.u j) ≠ 0) ∧
  (∀ j, X04.dot (H.u j) (H.C j) = 0) ∧ (∀ j, Fw H.u H.w 1 (H.C j) = 1) ∧
  0 < det2 (H.C 0) (H.C 1) ∧ 0 < det2 (H.C 1) (H.C 2) ∧ 0 < det2 (H.C 0) (H.C 2)

/-- Start corner `V_0` of edge `k`: `C_k`. -/
def eV0 (H : Hex) (k : Fin 3) : X04.Plane := H.C k

/-- End corner `V_1` of edge `k`: `C_{k+1}` for `k = 0, 1`, and `-C_0` for `k = 2`. -/
def eV1 (H : Hex) (k : Fin 3) : X04.Plane := if k = 2 then -H.C 0 else H.C (k + 1)

/-- The point `y(t) = V_0 + t(V_1 - V_0)` of edge `k`. -/
noncomputable def ept (H : Hex) (k : Fin 3) (t : ℝ) : X04.Plane :=
  eV0 H k + t • (eV1 H k - eV0 H k)

/-- The arc `x_p(t) = F_p(y(t))^{-1/p} y(t)` of `∂K_p` over edge `k`. -/
noncomputable def arc (H : Hex) (k : Fin 3) (p t : ℝ) : X04.Plane :=
  (Fw H.u H.w p (ept H k t)) ^ (-(1 / p)) • ept H k t

/-- The length element `|x_p'(t)|` (the supplement's `G(t,d)`, `d = p - 1`). -/
noncomputable def Gl (H : Hex) (k : Fin 3) (p t : ℝ) : ℝ := ‖deriv (arc H k p) t‖

/-- `E_e = 3!∫_0^1 E_3 dt = ∫_0^1 ∂_p^3 G(t, p-1)|_{p=1} dt` for edge `k`. -/
noncomputable def Ee (H : Hex) (k : Fin 3) : ℝ :=
  ∫ t in (0 : ℝ)..1, iteratedDeriv 3 (fun p => Gl H k p t) 1

/-- The unit vector `e_θ = (cos θ, sin θ)`. -/
noncomputable def ray (θ : ℝ) : X04.Plane := !₂[Real.cos θ, Real.sin θ]

/-- The polar curve `θ ↦ F_p(e_θ)^{-1/p} e_θ` (the boundary `{F_p = 1}`). -/
noncomputable def polar (H : Hex) (p θ : ℝ) : X04.Plane :=
  (Fw H.u H.w p (ray θ)) ^ (-(1 / p)) • ray θ

/-- `P(p)`: the Euclidean length of the closed curve `∂K_p` (total variation of the polar
parametrisation on `[0, 2π]`). -/
noncomputable def Perim (H : Hex) (p : ℝ) : ℝ :=
  (eVariationOn (polar H p) (Icc 0 (2 * π))).toReal

/-- Edge data `A_0 = a_{j₀}(V_1)`. -/
noncomputable def A0 (H : Hex) (k : Fin 3) : ℝ := term H.u H.w k (eV1 H k)
/-- Edge data `A_1 = a_{j₁}(V_0)`. -/
noncomputable def A1 (H : Hex) (k : Fin 3) : ℝ := term H.u H.w (k + 1) (eV0 H k)
/-- Edge data `B = a_{j₂}(V_0)`. -/
noncomputable def Bv (H : Hex) (k : Fin 3) : ℝ := term H.u H.w (k + 2) (eV0 H k)
/-- Edge data `r = a_{j₂}(V_1)/a_{j₂}(V_0)`. -/
noncomputable def rr (H : Hex) (k : Fin 3) : ℝ := term H.u H.w (k + 2) (eV1 H k) / Bv H k
/-- `L² = |V_1 - V_0|²`. -/
noncomputable def Lsq (H : Hex) (k : Fin 3) : ℝ := X04.dot (eV1 H k - eV0 H k) (eV1 H k - eV0 H k)
/-- `L = |V_1 - V_0|`. -/
noncomputable def Len (H : Hex) (k : Fin 3) : ℝ := Real.sqrt (Lsq H k)
/-- `q_0 = ⟨V_0, V_1 - V_0⟩`. -/
noncomputable def q0 (H : Hex) (k : Fin 3) : ℝ := X04.dot (eV0 H k) (eV1 H k - eV0 H k)
/-- `V_{02} = |V_0|²`. -/
noncomputable def V02 (H : Hex) (k : Fin 3) : ℝ := X04.dot (eV0 H k) (eV0 H k)
/-- `l_0 = log(A_0/w_{j₀})`. -/
noncomputable def l0 (H : Hex) (k : Fin 3) : ℝ := Real.log (A0 H k / H.w k)
/-- `l_1 = log(A_1/w_{j₁})`. -/
noncomputable def l1 (H : Hex) (k : Fin 3) : ℝ := Real.log (A1 H k / H.w (k + 1))
/-- `l_B = log(B/w_{j₂})`. -/
noncomputable def lB (H : Hex) (k : Fin 3) : ℝ := Real.log (Bv H k / H.w (k + 2))

/-- The eight edge variables `(B, r, q_0, V_{02}, L², l_0, l_1, l_B)`. -/
noncomputable def evars (H : Hex) (k : Fin 3) : Fin 8 → ℝ :=
  ![Bv H k, rr H k, q0 H k, V02 H k, Lsq H k, l0 H k, l1 H k, lB H k]

/-- Evaluation of an integer polynomial in the eight edge variables. -/
noncomputable def ev (H : Hex) (k : Fin 3) (c : MvPolynomial (Fin 8) ℤ) : ℝ :=
  MvPolynomial.eval₂ (Int.castRingHom ℝ) (evars H k) c

/-- `Δ = q_0² - L²V_{02}`. -/
def Dl (q V L2 : ℝ) : ℝ := q ^ 2 - L2 * V

/-- Printed block `c_8` (coefficient of `ζ(3)`). -/
def c8 (B r q V L2 : ℝ) : ℝ :=
  -6 * (B - 1) * (B * r - 1) * Dl q V L2 * (L2 - 3 * B * L2 + 3 * B * q * (r - 1))
/-- Printed block `c_1` (coefficient of `Li₃(x)`). -/
def c1 (B r q V L2 : ℝ) : ℝ :=
  -6 * B * Dl q V L2 * (-3 * B * L2 + 3 * B * q * (r - 1) + (4 - r) * L2 - 3 * q * (r - 1))
/-- Printed block `c_2` (coefficient of `Li₃(x/(x-1))`). -/
def c2 (B r q V L2 : ℝ) : ℝ :=
  -6 * B * Dl q V L2 *
    (-3 * B * L2 * r + 3 * B * q * r * (r - 1) + (r + 2) * L2 - 3 * q * (r - 1))
/-- Printed block `c_3` (coefficient of `Li₂(x)ℓ`). -/
def c3 (B r q V L2 : ℝ) : ℝ :=
  -6 * B * (B * r - 1) * Dl q V L2 * (-2 * L2 + 3 * q * (r - 1))
/-- Printed block `c_9` (coefficient of `ζ(2)ℓ`). -/
def c9 (B r q V L2 : ℝ) : ℝ :=
  18 * B * q * r * (B - 1) * (B * r - 1) * Dl q V L2

/-- The right side of Lemma g4:lem:precise for edge `k`, with the five printed blocks and the
six unprinted ones given as integer polynomials. `x = 1 - r`, `ℓ = -log r`, `ζ(2) = π²/6`. -/
noncomputable def closedForm (H : Hex) (k : Fin 3)
    (c4 c5 c6 c7 c10 c11 : MvPolynomial (Fin 8) ℤ) : ℝ :=
  let B := Bv H k
  let r := rr H k
  let q := q0 H k
  let V := V02 H k
  let L2 := Lsq H k
  let x := 1 - r
  let ℓ := -Real.log r
  c1 B r q V L2 * Z03.Li3 x + c2 B r q V L2 * Z03.Li3 (x / (x - 1))
    + c3 B r q V L2 * D00.Li2 x * ℓ + ev H k c4 * ℓ ^ 3 + ev H k c5 * D00.Li2 x
    + ev H k c6 * ℓ ^ 2 + ev H k c7 * ℓ + c8 B r q V L2 * Z03.z3
    + c9 B r q V L2 * (π ^ 2 / 6) * ℓ + ev H k c10 * (π ^ 2 / 6) + ev H k c11

/-- A plane vector as a pair (for T12's edge data). -/
def toPair (x : X04.Plane) : ℝ × ℝ := (x 0, x 1)

/-- Edge `k` in T12's format: the kink vectors through `V_0`, `V_1` scaled by their weights
(so T12's `ξ_0 = A_0`, `ξ_1 = A_1`), and the two corners. -/
def t12Edge (H : Hex) (k : Fin 3) : T12.Edge :=
  ⟨toPair (H.w k • H.u k), toPair (H.w (k + 1) • H.u (k + 1)), toPair (eV0 H k),
    toPair (eV1 H k)⟩

/-- `Φ_A` of prop:nonlocal for edge `k` (`γ_e = Bh²/L³`, `κ = r - 1`). -/
noncomputable def PhiA (H : Hex) (k : Fin 3) : ℝ :=
  let L := Len H k
  let h := T12.edgeH (toPair (eV0 H k)) (toPair (eV1 H k))
  let gam := Bv H k * h ^ 2 / L ^ 3
  let kap := rr H k - 1
  (-6) * gam * (3 * A0 H k * kap * q0 H k - L ^ 2 * (2 * A0 H k + rr H k * A1 H k))

/-- `Φ_L` of prop:nonlocal for edge `k`. -/
noncomputable def PhiL (H : Hex) (k : Fin 3) : ℝ :=
  let L := Len H k
  let h := T12.edgeH (toPair (eV0 H k)) (toPair (eV1 H k))
  let gam := Bv H k * h ^ 2 / L ^ 3
  let kap := rr H k - 1
  let a0 := A0 H k
  let a1 := A1 H k
  let B := Bv H k
  (-18) * gam * (kap * (-a1 ^ 2 + a1 * B * (2 * kap - 1) + 2 * B ^ 2 * kap) * q0 H k
    + L ^ 2 * (a0 * B + rr H k * a1 ^ 2))

/-- `𝓔 = Σ_e h_e²/L_e` over the six edges (opposite edges have the same `h`, `L`). -/
noncomputable def calE (H : Hex) : ℝ := 2 * ∑ k : Fin 3, T12.hsqOverL (t12Edge H k)

/-- `H^{(2)}_n = Σ_{j=1}^n 1/j²`. -/
noncomputable def H2 (n : ℕ) : ℝ := ∑ j ∈ Finset.Icc 1 n, 1 / (j : ℝ) ^ 2

/-- `𝓛_1 = -log(1-x)`. -/
noncomputable def cL1 (x : ℝ) : ℝ := -Real.log (1 - x)

/-- `λ_c(q)` for `c = 1, 2, 3` from `c![v^c](-(v/q)exp(-Σ_i v^iH^{(i)}_{q-1}/i))`:
`λ_1 = -1/q`, `λ_2 = 2H_{q-1}/q`, `λ_3 = -3(H_{q-1}² - H^{(2)}_{q-1})/q` (used for `q ≥ 1`). -/
noncomputable def lam : ℕ → ℕ → ℝ
  | 1, q => -1 / (q : ℝ)
  | 2, q => 2 * (harmonic (q - 1) : ℝ) / q
  | 3, q => -3 * ((harmonic (q - 1) : ℝ) ^ 2 - H2 (q - 1)) / q
  | _, _ => 0

/-- `J(n;a,b,c)(x) = ∫_0^1 t^n log^a t log^b(1-t) log^c(1-xt) dt`. -/
noncomputable def Jint (n a b c : ℕ) (x : ℝ) : ℝ :=
  ∫ t in (0 : ℝ)..1,
    t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b * Real.log (1 - x * t) ^ c

/-- `𝓘_N(a,b) = ∫_0^1 s^N log^a s log^b(1-s) ds`. -/
noncomputable def Iint (N a b : ℕ) : ℝ :=
  ∫ s in (0 : ℝ)..1, s ^ N * Real.log s ^ a * Real.log (1 - s) ^ b

/-- Lemma g4:lem:euler: the six generating functions, for `|x| < 1`. -/
def stmt_1 : Prop :=
  ∀ x : ℝ, |x| < 1 →
    HasSum (fun n : ℕ => (harmonic (n + 1) : ℝ) * x ^ (n + 1)) (cL1 x / (1 - x)) ∧
    HasSum (fun n : ℕ => H2 (n + 1) * x ^ (n + 1)) (D00.Li2 x / (1 - x)) ∧
    HasSum (fun n : ℕ => (harmonic (n + 1) : ℝ) ^ 2 * x ^ (n + 1))
      ((cL1 x ^ 2 + D00.Li2 x) / (1 - x)) ∧
    HasSum (fun n : ℕ => (harmonic (n + 1) : ℝ) / ((n : ℝ) + 1) * x ^ (n + 1))
      (D00.Li2 x + (1 / 2) * cL1 x ^ 2) ∧
    HasSum (fun n : ℕ => (harmonic (n + 1) : ℝ) ^ 2 / ((n : ℝ) + 1) * x ^ (n + 1))
      (Z03.Li3 x + cL1 x * D00.Li2 x + (1 / 3) * cL1 x ^ 3) ∧
    HasSum (fun n : ℕ => H2 (n + 1) / ((n : ℝ) + 1) * x ^ (n + 1))
      (3 * Z03.Li3 x - 2 * Z03.A x + cL1 x * D00.Li2 x)

/-- Lemma g4:lem:Jseries, first part: "$\log^c(1-z)=\sum_{q\ge1}\lambda_c(q)z^q$" with the
explicit `λ_c`, for `c = 1, 2, 3` and `|z| < 1`. -/
def stmt_2 : Prop :=
  ∀ c : ℕ, 1 ≤ c → c ≤ 3 → ∀ z : ℝ, |z| < 1 →
    HasSum (fun q : ℕ => lam c (q + 1) * z ^ (q + 1)) (Real.log (1 - z) ^ c)

/-- Lemma g4:lem:Jseries: for `|x| < 1`, `c ≥ 1`, `a + b + c ≤ 3`, the integrand of `J` is
integrable and `J(n;a,b,c) = Σ_{q≥1} x^q λ_c(q) 𝓘_{n+q}(a,b)`. -/
def stmt_3 : Prop :=
  ∀ x : ℝ, |x| < 1 → ∀ n a b c : ℕ, 1 ≤ c → a + b + c ≤ 3 →
    IntervalIntegrable
        (fun t : ℝ => t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b * Real.log (1 - x * t) ^ c)
        volume 0 1 ∧
      HasSum (fun q : ℕ => x ^ (q + 1) * lam c (q + 1) * Iint (n + q + 1) a b) (Jint n a b c x)

/-- The worked entries after Lemma g4:lem:Jseries, for `0 < |x| < 1` (they divide by `x`;
checklist note 24): `J(0;1,0,1)`, `J(0;0,0,1)`, `J(1;0,0,1)`, `J(0;0,2,1)`. -/
def stmt_4 : Prop :=
  ∀ x : ℝ, |x| < 1 → x ≠ 0 →
    Jint 0 1 0 1 x = (-cL1 x - D00.Li2 x + x * (cL1 x + 2)) / x ∧
    Jint 0 0 0 1 x = -cL1 x + cL1 x / x - 1 ∧
    Jint 1 0 0 1 x = -cL1 x / 2 + cL1 x / (2 * x ^ 2) - 1 / 4 - 1 / (2 * x) ∧
    Jint 0 0 2 1 x =
      (1 - x) / (3 * x) * (12 * Z03.Li3 x - 6 * Z03.A x + 6 * cL1 x * D00.Li2 x
        + 6 * D00.Li2 x + cL1 x ^ 3 + 3 * cL1 x ^ 2 + 6 * cL1 x) - 6

/-- Setting (main.tex 1709): every strictly convex centrally symmetric hexagon with corners
`V_0, V_1, V_2, -V_0, -V_1, -V_2` (counter-clockwise), with every choice of weights `w_j > 0`,
is `K_1` for suitable `u_j = k_jv_j` with `V_j` on the kink line `⟨u_j,x⟩ = 0`. -/
def stmt_5 : Prop :=
  ∀ (V : Fin 3 → X04.Plane) (w : Fin 3 → ℝ), (∀ j, 0 < w j) →
    0 < det2 (V 0) (V 1) → 0 < det2 (V 1) (V 2) → 0 < det2 (V 0) (V 2) →
    0 < det2 (V 1 - V 0) (V 2 - V 1) → 0 < det2 (V 2 - V 1) (-V 0 - V 2) →
    0 < det2 (-V 0 - V 2) (-V 1 + V 0) →
    ∃ u : Fin 3 → X04.Plane, IsHex ⟨u, w, V⟩

/-- Setting: `K_1` is the hexagon with the six corners `±C_j`. -/
def stmt_6 : Prop :=
  ∀ H : Hex, IsHex H →
    Kw H.u H.w 1 = convexHull ℝ (Set.range H.C ∪ Set.range (fun j => -H.C j))

/-- Lemma g4:lem:arc (setting): for `p > 0` the polar curve is a simple closed curve of finite
length whose image is `{F_p = 1} = ∂K_p`. -/
def stmt_7 : Prop :=
  ∀ H : Hex, IsHex H → ∀ p : ℝ, 0 < p →
    ContinuousOn (polar H p) (Icc 0 (2 * π)) ∧ polar H p (2 * π) = polar H p 0 ∧
    InjOn (polar H p) (Ico 0 (2 * π)) ∧
    polar H p '' Icc 0 (2 * π) = {x | Fw H.u H.w p x = 1} ∧
    frontier (Kw H.u H.w p) = {x | Fw H.u H.w p x = 1} ∧
    eVariationOn (polar H p) (Icc 0 (2 * π)) ≠ ⊤

/-- "One edge" and "What is proved": on each of the three edges the terms are
`A_0t`, `A_1(1-t)`, `B(1+κt)` with `κ = r - 1`; `A_1 + B = 1`, `A_0 + Br = 1`; `0 < B < 1`,
`A_0, A_1 > 0`, `r > 0`. -/
def stmt_8 : Prop :=
  ∀ H : Hex, IsHex H → ∀ k : Fin 3,
    (∀ t ∈ Icc (0 : ℝ) 1,
      term H.u H.w k (ept H k t) = A0 H k * t ∧
      term H.u H.w (k + 1) (ept H k t) = A1 H k * (1 - t) ∧
      term H.u H.w (k + 2) (ept H k t) = Bv H k * (1 + (rr H k - 1) * t)) ∧
    A1 H k + Bv H k = 1 ∧ A0 H k + Bv H k * rr H k = 1 ∧
    0 < Bv H k ∧ Bv H k < 1 ∧ 0 < A0 H k ∧ 0 < A1 H k ∧ 0 < rr H k

/-- First sentence of thm:hexagon: `P` is three times differentiable at `1`, the edge
integrands are integrable, and `P'''(1) = 2Σ_e E_e` over the three edges
`C_0C_1`, `C_1C_2`, `C_2(-C_0)`. -/
def stmt_9 : Prop :=
  ∀ H : Hex, IsHex H →
    ContDiffAt ℝ 3 (Perim H) 1 ∧
    (∀ k : Fin 3,
      IntervalIntegrable (fun t => iteratedDeriv 3 (fun p => Gl H k p t) 1) volume 0 1) ∧
    iteratedDeriv 3 (Perim H) 1 = 2 * ∑ k : Fin 3, Ee H k

/-- Third sentence: the coefficient of `ζ(3)`, `c_8/L⁵`, is exactly `Z_e` of prop:edgelaw
(T12's `Ze`, with `ξ_0 = A_0`, `ξ_1 = A_1`), for every edge and every `r`. -/
def stmt_11 : Prop :=
  ∀ H : Hex, IsHex H → ∀ k : Fin 3,
    c8 (Bv H k) (rr H k) (q0 H k) (V02 H k) (Lsq H k) = Len H k ^ 5 * T12.Ze (t12Edge H k)

/-- Fourth sentence: the coefficients of `Li₃(x)` and `Li₃(x/(x-1))` are `2Φ_A + Φ_L` and
`Φ_A`: `c_1 = L⁵(2Φ_A + Φ_L)`, `c_2 = L⁵Φ_A`. -/
def stmt_12 : Prop :=
  ∀ H : Hex, IsHex H → ∀ k : Fin 3,
    c1 (Bv H k) (rr H k) (q0 H k) (V02 H k) (Lsq H k)
        = Len H k ^ 5 * (2 * PhiA H k + PhiL H k) ∧
      c2 (Bv H k) (rr H k) (q0 H k) (V02 H k) (Lsq H k) = Len H k ^ 5 * PhiA H k

/-- Fifth sentence, `r = 1` on every edge and all `w_j = 1`: `P'''(1)` is `3!` times the `d³`
coefficient of the hexagon law of prop:edgelaw,
`2^{d/(1+d)}[P(1) + π²𝓔d²/24 - (𝓔/8)(ζ(3) + 2π²/3)d³]` (T12's `hexDisp` with `z3 = ζ(3)`). -/
def stmt_13 : Prop :=
  ∀ H : Hex, IsHex H → (∀ j, H.w j = 1) → (∀ k, rr H k = 1) →
    ContDiffAt ℝ 3 (Perim H) 1 ∧
    iteratedDeriv 3 (Perim H) 1
      = iteratedDeriv 3 (fun d => T12.hexDisp (Perim H 1) (calE H) Z03.z3 d) 0

end S05
-- END STATEMENTS

namespace S05
open Real Set

lemma dot_ept (u : X04.Plane) (H : Hex) (k : Fin 3) (t : ℝ) :
    X04.dot u (ept H k t) = (1 - t) * X04.dot u (eV0 H k) + t * X04.dot u (eV1 H k) := by
  simp only [ept, X04.dot, PiLp.add_apply, PiLp.smul_apply, PiLp.sub_apply, smul_eq_mul]; ring

lemma dot_det (u C x y : X04.Plane) (h : X04.dot u C = 0) :
    X04.dot u x * det2 C y = X04.dot u y * det2 C x := by
  unfold X04.dot det2 at *; linear_combination (x 0 * y 1 - x 1 * y 0) * h

lemma dot_ne (u u' C x : X04.Plane) (hu : det2 u u' ≠ 0) (h : X04.dot u C = 0)
    (hx : det2 C x ≠ 0) : X04.dot u x ≠ 0 := by
  intro h0
  have key : ∀ y, X04.dot u y = 0 := by
    intro y
    have e := dot_det u C x y h
    rw [h0, zero_mul] at e
    rcases mul_eq_zero.1 e.symm with h1 | h1
    · exact h1
    · exact absurd h1 hx
  have e0 := key (EuclideanSpace.single 0 1)
  have e1 := key (EuclideanSpace.single 1 1)
  simp [X04.dot] at e0 e1
  apply hu; simp [det2, e0, e1]

lemma sameSign (u C x y : X04.Plane) (h : X04.dot u C = 0) (hxy : 0 < det2 C x * det2 C y) :
    0 ≤ X04.dot u x * X04.dot u y := by
  have e := dot_det u C x y h
  have h2 : X04.dot u x * X04.dot u y * (det2 C x * det2 C y) = (X04.dot u x) ^ 2 * (det2 C y) ^ 2 := by
    linear_combination (-(X04.dot u x * det2 C y)) * e
  have : 0 ≤ X04.dot u x * X04.dot u y * (det2 C x * det2 C y) := by rw [h2]; positivity
  exact nonneg_of_mul_nonneg_left this hxy

lemma abs_comb (a b t : ℝ) (hab : 0 ≤ a * b) (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    |(1 - t) * a + t * b| = (1 - t) * |a| + t * |b| := by
  rcases le_total 0 a with ha | ha <;> rcases le_total 0 b with hb | hb
  · rw [abs_of_nonneg ha, abs_of_nonneg hb, abs_of_nonneg (by nlinarith)]
  · have : a * b ≤ 0 := mul_nonpos_of_nonneg_of_nonpos ha hb
    have hab0 : a * b = 0 := le_antisymm this hab
    rcases mul_eq_zero.1 hab0 with h | h <;> subst h <;> simp [abs_mul, abs_of_nonneg ht0, abs_of_nonneg (sub_nonneg.2 ht1)]
  · have : a * b ≤ 0 := mul_nonpos_of_nonpos_of_nonneg ha hb
    have hab0 : a * b = 0 := le_antisymm this hab
    rcases mul_eq_zero.1 hab0 with h | h <;> subst h <;> simp [abs_mul, abs_of_nonneg ht0, abs_of_nonneg (sub_nonneg.2 ht1)]
  · rw [abs_of_nonpos ha, abs_of_nonpos hb, abs_of_nonpos (by nlinarith)]; ring

lemma Fw_one (u : Fin 3 → X04.Plane) (w : Fin 3 → ℝ) (x : X04.Plane) :
    Fw u w 1 x = w 0 * |X04.dot (u 0) x| + w 1 * |X04.dot (u 1) x| + w 2 * |X04.dot (u 2) x| := by
  simp [Fw, Fin.sum_univ_three]

/-- Core of stmt_8 on one edge, with the roles `j₀, j₁, j₂` abstracted. -/
lemma edge_core (w0 w1 w2 a0 a1 b0 b1 : ℝ) (hw0 : 0 < w0) (hw1 : 0 < w1) (hw2 : 0 < w2)
    (ha0 : a0 ≠ 0) (ha1 : a1 ≠ 0) (hb0 : b0 ≠ 0) (hb1 : b1 ≠ 0)
    (hF0 : w1 * |a1| + w2 * |b0| = 1) (hF1 : w0 * |a0| + w2 * |b1| = 1) :
    w1 * |a1| + w2 * |b0| = 1 ∧ w0 * |a0| + w2 * |b0| * (w2 * |b1| / (w2 * |b0|)) = 1 ∧
    0 < w2 * |b0| ∧ w2 * |b0| < 1 ∧ 0 < w0 * |a0| ∧ 0 < w1 * |a1| ∧
    0 < w2 * |b1| / (w2 * |b0|) := by
  have p0 : 0 < |a0| := abs_pos.2 ha0
  have p1 : 0 < |a1| := abs_pos.2 ha1
  have pb0 : 0 < |b0| := abs_pos.2 hb0
  have pb1 : 0 < |b1| := abs_pos.2 hb1
  have hB : 0 < w2 * |b0| := by positivity
  refine ⟨hF0, ?_, hB, ?_, by positivity, by positivity, by positivity⟩
  · rw [mul_div_cancel₀ _ hB.ne']; exact hF1
  · nlinarith [mul_pos hw1 p1]


lemma edge_gen (H : Hex) (k j1 j2 : Fin 3) (V0 V1 C : X04.Plane) (hj1 : k + 1 = j1) (hj2 : k + 2 = j2)
    (hV0 : eV0 H k = V0) (hV1 : eV1 H k = V1) (hw : ∀ j, 0 < H.w j)
    (h0 : X04.dot (H.u k) V0 = 0) (h1 : X04.dot (H.u j1) V1 = 0) (hC : X04.dot (H.u j2) C = 0)
    (hCs : 0 < det2 C V0 * det2 C V1)
    (n0 : X04.dot (H.u k) V1 ≠ 0) (n1 : X04.dot (H.u j1) V0 ≠ 0) (n2 : X04.dot (H.u j2) V0 ≠ 0)
    (n3 : X04.dot (H.u j2) V1 ≠ 0)
    (F0 : H.w k * |X04.dot (H.u k) V0| + H.w j1 * |X04.dot (H.u j1) V0|
      + H.w j2 * |X04.dot (H.u j2) V0| = 1)
    (F1 : H.w k * |X04.dot (H.u k) V1| + H.w j1 * |X04.dot (H.u j1) V1|
      + H.w j2 * |X04.dot (H.u j2) V1| = 1) :
    (∀ t ∈ Icc (0 : ℝ) 1,
      term H.u H.w k (ept H k t) = A0 H k * t ∧
      term H.u H.w (k + 1) (ept H k t) = A1 H k * (1 - t) ∧
      term H.u H.w (k + 2) (ept H k t) = Bv H k * (1 + (rr H k - 1) * t)) ∧
    A1 H k + Bv H k = 1 ∧ A0 H k + Bv H k * rr H k = 1 ∧
    0 < Bv H k ∧ Bv H k < 1 ∧ 0 < A0 H k ∧ 0 < A1 H k ∧ 0 < rr H k := by
  have hss := sameSign (H.u j2) C V0 V1 hC hCs
  simp only [A0, A1, Bv, rr, term, hj1, hj2, hV0, hV1]
  rw [h0, abs_zero, mul_zero, zero_add] at F0
  rw [h1, abs_zero, mul_zero, add_zero] at F1
  obtain ⟨e1, e2, e3, e4, e5, e6, e7⟩ := edge_core (H.w k) (H.w j1) (H.w j2) _ _ _ _ (hw k) (hw j1)
    (hw j2) n0 n1 n2 n3 F0 F1
  refine ⟨?_, e1, e2, e3, e4, e5, e6, e7⟩
  intro t ht
  rw [dot_ept, dot_ept, dot_ept, hV0, hV1, h0, h1]
  have ht0 := ht.1; have ht1 := ht.2
  refine ⟨?_, ?_, ?_⟩
  · rw [mul_zero, zero_add, abs_mul, abs_of_nonneg ht0]; ring
  · rw [mul_zero, add_zero, abs_mul, abs_of_nonneg (sub_nonneg.2 ht1)]; ring
  · rw [abs_comb _ _ _ hss ht0 ht1]
    have hB : 0 < H.w j2 * |X04.dot (H.u j2) V0| := e3
    have := (hw j2).ne'
    have := (abs_pos.2 n2).ne'
    field_simp
    ring

end S05

theorem S05.proof_8 : S05.stmt_8 := by
  intro H hH k
  obtain ⟨hw, hdet, hker, hF, h01, h12, h02⟩ := hH
  have nz : ∀ i j : Fin 3, i ≠ j → ∀ x, S05.det2 (H.C i) x ≠ 0 → X04.dot (H.u i) x ≠ 0 :=
    fun i j hij x hx => S05.dot_ne _ _ _ _ (hdet i j hij) (hker i) hx
  have dneg : ∀ x y : X04.Plane, S05.det2 x (-y) = -S05.det2 x y := by
    intro x y; simp [S05.det2]; ring
  have dsw : ∀ x y : X04.Plane, S05.det2 x y = -S05.det2 y x := by
    intro x y; simp [S05.det2]; ring
  have Fs : ∀ j, H.w 0 * |X04.dot (H.u 0) (H.C j)| + H.w 1 * |X04.dot (H.u 1) (H.C j)|
      + H.w 2 * |X04.dot (H.u 2) (H.C j)| = 1 := fun j => by rw [← S05.Fw_one]; exact hF j
  have dotneg : ∀ u x : X04.Plane, X04.dot u (-x) = -X04.dot u x := by
    intro u x; simp [X04.dot]; ring
  obtain rfl | rfl | rfl : k = 0 ∨ k = 1 ∨ k = 2 := by fin_cases k <;> simp
  · refine S05.edge_gen H 0 1 2 (H.C 0) (H.C 1) (H.C 2) rfl rfl rfl (by simp [S05.eV1]) hw
      (hker 0) (hker 1) (hker 2) ?_ (nz 0 1 (by decide) _ h01.ne') (nz 1 0 (by decide) _ ?_)
      (nz 2 0 (by decide) _ ?_) (nz 2 0 (by decide) _ ?_) ?_ ?_
    · rw [dsw (H.C 2) (H.C 0), dsw (H.C 2) (H.C 1)]; nlinarith
    · rw [dsw]; linarith
    · rw [dsw]; linarith
    · rw [dsw]; linarith
    · linarith [Fs 0]
    · linarith [Fs 1]
  · refine S05.edge_gen H 1 2 0 (H.C 1) (H.C 2) (H.C 0) rfl rfl rfl (by simp [S05.eV1]) hw
      (hker 1) (hker 2) (hker 0) ?_ (nz 1 0 (by decide) _ h12.ne') (nz 2 0 (by decide) _ ?_)
      (nz 0 1 (by decide) _ h01.ne') (nz 0 1 (by decide) _ h02.ne') ?_ ?_
    · nlinarith
    · rw [dsw]; linarith
    · linarith [Fs 1]
    · linarith [Fs 2]
  · refine S05.edge_gen H 2 0 1 (H.C 2) (-H.C 0) (H.C 1) rfl rfl rfl (by simp [S05.eV1]) hw
      (hker 2) (by rw [dotneg, hker 0, neg_zero]) (hker 1) ?_ ?_ (nz 0 1 (by decide) _ ?_)
      (nz 1 0 (by decide) _ h12.ne') (nz 1 0 (by decide) _ ?_) ?_ ?_
    · rw [dneg, dsw (H.C 1) (H.C 0)]; nlinarith
    · rw [dotneg]; exact neg_ne_zero.2 (nz 2 0 (by decide) _ (by rw [dsw (H.C 2) (H.C 0)]; linarith))
    · exact h02.ne'
    · rw [dneg, dsw]; linarith
    · linarith [Fs 2]
    · have := Fs 0; simp only [dotneg, abs_neg]; linarith

namespace S05
open Real Set

lemma t12dot_smul (c : ℝ) (x y : X04.Plane) :
    T12.dot (toPair (c • x)) (toPair y) = c * X04.dot x y := by
  simp [T12.dot, toPair, X04.dot]; ring

lemma toPair_sub (x y : X04.Plane) : toPair (x - y) = toPair x - toPair y := by
  ext <;> simp [toPair]

lemma t12dot (x y : X04.Plane) : T12.dot (toPair x) (toPair y) = X04.dot x y := rfl

lemma edgeLen_eq (H : Hex) (k : Fin 3) :
    T12.edgeLen (toPair (eV0 H k)) (toPair (eV1 H k)) = Len H k := by
  unfold T12.edgeLen Len Lsq; rw [← toPair_sub, t12dot]

lemma Lsq_pos (H : Hex) (hH : IsHex H) (k : Fin 3) : 0 < Lsq H k := by
  have h8 := (S05.proof_8 H hH k)
  have hA0 := h8.2.2.2.2.2.1
  have hne : eV1 H k - eV0 H k ≠ 0 := by
    intro h
    have : eV1 H k = eV0 H k := sub_eq_zero.1 h
    have e := (h8.1 0 ⟨le_rfl, zero_le_one⟩).1
    simp only [mul_zero] at e
    have e2 : term H.u H.w k (ept H k 0) = A0 H k := by
      simp [A0, ept, this]
    linarith
  unfold Lsq X04.dot
  by_contra hc
  rw [not_lt] at hc
  have h0 : (eV1 H k - eV0 H k) 0 = 0 := by nlinarith [sq_nonneg ((eV1 H k - eV0 H k) 0), sq_nonneg ((eV1 H k - eV0 H k) 1)]
  have h1 : (eV1 H k - eV0 H k) 1 = 0 := by nlinarith [sq_nonneg ((eV1 H k - eV0 H k) 0), sq_nonneg ((eV1 H k - eV0 H k) 1)]
  apply hne
  ext i; fin_cases i <;> simp [h0, h1]

lemma foot_alg (a0 a1 b0 b1 L : ℝ) (hL : 0 < L)
    (hL2 : L ^ 2 = (b0 - a0) * (b0 - a0) + (b1 - a1) * (b1 - a1)) :
    (a0 + (-(a0 * (b0 - a0) + a1 * (b1 - a1)) / L / L) * (b0 - a0)) *
        (a0 + (-(a0 * (b0 - a0) + a1 * (b1 - a1)) / L / L) * (b0 - a0)) +
      (a1 + (-(a0 * (b0 - a0) + a1 * (b1 - a1)) / L / L) * (b1 - a1)) *
        (a1 + (-(a0 * (b0 - a0) + a1 * (b1 - a1)) / L / L) * (b1 - a1))
      = (a0 * a0 + a1 * a1) - (a0 * (b0 - a0) + a1 * (b1 - a1)) ^ 2 / L ^ 2 := by
  have hL0 : L ≠ 0 := hL.ne'
  field_simp
  linear_combination (-(a0 * (b0 - a0) + a1 * (b1 - a1)) ^ 2) * hL2

lemma edgeH_sq (H : Hex) (hH : IsHex H) (k : Fin 3) :
    T12.edgeH (toPair (eV0 H k)) (toPair (eV1 H k)) ^ 2 = V02 H k - q0 H k ^ 2 / Lsq H k := by
  have hL := Lsq_pos H hH k
  have hLn : Len H k ^ 2 = Lsq H k := Real.sq_sqrt hL.le
  have hLp : 0 < Len H k := Real.sqrt_pos.2 hL
  unfold T12.edgeH
  rw [Real.sq_sqrt (by unfold T12.dot; exact add_nonneg (mul_self_nonneg _) (mul_self_nonneg _))]
  unfold T12.foot T12.edgeL0
  rw [edgeLen_eq, ← toPair_sub, t12dot]
  rw [← hLn]
  unfold Lsq at hLn
  simp only [T12.dot, toPair, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
    Prod.fst_sub, Prod.snd_sub, smul_eq_mul, V02, q0, X04.dot, PiLp.sub_apply] at hLn ⊢
  exact foot_alg _ _ _ _ _ hLp hLn

lemma xi0_eq (H : Hex) (hH : IsHex H) (k : Fin 3) :
    |T12.dot (toPair (H.w k • H.u k)) (toPair (eV1 H k))| = A0 H k := by
  rw [t12dot_smul, abs_mul, abs_of_pos (hH.1 k)]; rfl

lemma xi1_eq (H : Hex) (hH : IsHex H) (k : Fin 3) :
    |T12.dot (toPair (H.w (k + 1) • H.u (k + 1))) (toPair (eV0 H k))| = A1 H k := by
  rw [t12dot_smul, abs_mul, abs_of_pos (hH.1 (k + 1))]; rfl

lemma q0_t12 (H : Hex) (k : Fin 3) :
    T12.dot (toPair (eV0 H k)) (toPair (eV1 H k) - toPair (eV0 H k)) = q0 H k := by
  rw [← toPair_sub, t12dot]; rfl

end S05

theorem S05.proof_11 : S05.stmt_11 := by
  intro H hH k
  obtain ⟨_, hAB, hAr, hB0, hB1, hA0, hA1, hr⟩ := S05.proof_8 H hH k
  have hL := S05.Lsq_pos H hH k
  have hLn : S05.Len H k ^ 2 = S05.Lsq H k := Real.sq_sqrt hL.le
  have hLp : 0 < S05.Len H k := Real.sqrt_pos.2 hL
  have hh := S05.edgeH_sq H hH k
  unfold T12.Ze S05.t12Edge
  simp only
  rw [S05.xi0_eq H hH, S05.xi1_eq H hH, S05.edgeLen_eq, T12.edgeL1, T12.edgeL0, S05.edgeLen_eq,
    S05.q0_t12, hh]
  have e1 : S05.A1 H k = 1 - S05.Bv H k := by linarith
  have e0 : S05.A0 H k = 1 - S05.Bv H k * S05.rr H k := by linarith
  rw [e0, e1, ← hLn]
  unfold S05.c8 S05.Dl
  have := hLp.ne'
  field_simp
  ring

theorem S05.proof_12 : S05.stmt_12 := by
  intro H hH k
  obtain ⟨_, hAB, hAr, hB0, hB1, hA0, hA1, hr⟩ := S05.proof_8 H hH k
  have hL := S05.Lsq_pos H hH k
  have hLn : S05.Len H k ^ 2 = S05.Lsq H k := Real.sq_sqrt hL.le
  have hLp : 0 < S05.Len H k := Real.sqrt_pos.2 hL
  have hh := S05.edgeH_sq H hH k
  have e1 : S05.A1 H k = 1 - S05.Bv H k := by linarith
  have e0 : S05.A0 H k = 1 - S05.Bv H k * S05.rr H k := by linarith
  unfold S05.PhiA S05.PhiL
  simp only
  rw [hh, e0, e1, ← hLn]
  unfold S05.c1 S05.c2 S05.Dl
  have := hLp.ne'
  constructor
  · field_simp
    ring
  · field_simp
    ring


namespace S05
open Real Set

/-- The rotated, scaled kink normal `s • (-x₁, x₀)`, so that `⟨perpS s x, y⟩ = s · det(x, y)`. -/
noncomputable def perpS (s : ℝ) (x : X04.Plane) : X04.Plane := s • !₂[-(x 1), x 0]

lemma dot_perpS (s : ℝ) (x y : X04.Plane) : X04.dot (perpS s x) y = s * det2 x y := by
  simp [perpS, X04.dot, det2]; ring

lemma det_perpS (s t : ℝ) (x y : X04.Plane) :
    det2 (perpS s x) (perpS t y) = s * t * det2 x y := by
  simp [perpS, det2]; ring

end S05

theorem S05.proof_5 : S05.stmt_5 := by
  intro V w hw h01 h12 h02 hc1 hc2 hc3
  set a := S05.det2 (V 1) (V 2) with ha
  set b := S05.det2 (V 0) (V 2) with hb
  set c := S05.det2 (V 0) (V 1) with hc
  have e1 : S05.det2 (V 1 - V 0) (V 2 - V 1) = a - b + c := by
    simp only [ha, hb, hc, S05.det2, PiLp.sub_apply]; ring
  have e2 : S05.det2 (V 2 - V 1) (-V 0 - V 2) = a + b - c := by
    simp only [ha, hb, hc, S05.det2, PiLp.sub_apply, PiLp.neg_apply]; ring
  have e3 : S05.det2 (-V 0 - V 2) (-V 1 + V 0) = b + c - a := by
    simp only [ha, hb, hc, S05.det2, PiLp.sub_apply, PiLp.neg_apply, PiLp.add_apply]; ring
  rw [e1] at hc1; rw [e2] at hc2; rw [e3] at hc3
  set y : Fin 3 → ℝ := ![(b + c - a) / (2 * b * c), (a + c - b) / (2 * a * c),
    (a + b - c) / (2 * a * b)] with hy
  have hy0 : 0 < y 0 := by simp only [hy]; simp; exact div_pos (by linarith) (by positivity)
  have hy1 : 0 < y 1 := by simp only [hy]; simp; exact div_pos (by linarith) (by positivity)
  have hy2 : 0 < y 2 := by simp only [hy]; simp; exact div_pos (by linarith) (by positivity)
  have hyp : ∀ j, 0 < y j := by intro j; fin_cases j <;> assumption
  refine ⟨fun j => S05.perpS (y j / w j) (V j), ?_⟩
  have dsw : ∀ x z : X04.Plane, S05.det2 x z = -S05.det2 z x := by
    intro x z; simp [S05.det2]; ring
  have dself : ∀ x : X04.Plane, S05.det2 x x = 0 := by intro x; simp [S05.det2]; ring
  have hwy : ∀ j, w j * (y j / w j) = y j := fun j => by field_simp [(hw j).ne']
  refine ⟨hw, ?_, ?_, ?_, h01, h12, h02⟩
  · intro i j hij
    simp only [S05.det_perpS]
    have hs : ∀ j, 0 < y j / w j := fun j => div_pos (hyp j) (hw j)
    refine mul_ne_zero (mul_pos (hs i) (hs j)).ne' ?_
    fin_cases i <;> fin_cases j <;> simp at hij ⊢ <;>
      first
      | exact h01.ne' | exact h12.ne' | exact h02.ne'
      | (rw [dsw]; simpa using h01.ne') | (rw [dsw]; simpa using h12.ne')
      | (rw [dsw]; simpa using h02.ne')
  · intro j; simp only [S05.dot_perpS, dself, mul_zero]
  · have q1 : |S05.det2 (V 0) (V 1)| = c := abs_of_pos h01
    have q2 : |S05.det2 (V 1) (V 0)| = c := by rw [dsw, abs_neg]; exact abs_of_pos h01
    have q3 : |S05.det2 (V 0) (V 2)| = b := abs_of_pos h02
    have q4 : |S05.det2 (V 2) (V 0)| = b := by rw [dsw, abs_neg]; exact abs_of_pos h02
    have q5 : |S05.det2 (V 1) (V 2)| = a := abs_of_pos h12
    have q6 : |S05.det2 (V 2) (V 1)| = a := by rw [dsw, abs_neg]; exact abs_of_pos h12
    intro i
    simp only [S05.Fw, Fin.sum_univ_three, Real.rpow_one, S05.dot_perpS, abs_mul,
      abs_of_pos (div_pos (hyp _) (hw _)), ← mul_assoc, hwy]
    obtain rfl | rfl | rfl : i = 0 ∨ i = 1 ∨ i = 2 := by fin_cases i <;> simp
    all_goals
      simp only [dself, abs_zero, mul_zero, zero_add, add_zero, q1, q2, q3, q4, q5, q6, hy,
        Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
        Matrix.tail_cons]
      field_simp
      ring

namespace S05
open Real Set Finset

lemma hR_succ (n : ℕ) : (harmonic (n + 1) : ℝ) = harmonic n + 1 / ((n : ℝ) + 1) := by
  rw [harmonic_succ]; push_cast; ring

lemma hR_sum (n : ℕ) : (harmonic n : ℝ) = ∑ i ∈ range n, 1 / ((i : ℝ) + 1) := by
  induction n with
  | zero => simp
  | succ n ih => rw [hR_succ, ih, sum_range_succ]

lemma H2_succ (n : ℕ) : H2 (n + 1) = H2 n + 1 / ((n : ℝ) + 1) ^ 2 := by
  unfold H2; rw [sum_Icc_succ_top (by omega)]; push_cast; ring

lemma H2_zero : H2 0 = 0 := by simp [H2]

/-- `Σ_{k=0}^n 1/((k+1)(n-k+1)) = 2H_{n+1}/(n+2)`. -/
lemma conv2 (n : ℕ) :
    ∑ k ∈ range (n + 1), 1 / (((k : ℝ) + 1) * (((n - k : ℕ) : ℝ) + 1))
      = 2 * (harmonic (n + 1) : ℝ) / ((n : ℝ) + 2) := by
  have h : ∀ k ∈ range (n + 1), 1 / (((k : ℝ) + 1) * (((n - k : ℕ) : ℝ) + 1))
      = (1 / ((n : ℝ) + 2)) * (1 / ((k : ℝ) + 1) + 1 / (((n - k : ℕ) : ℝ) + 1)) := by
    intro k hk
    have hkn : k ≤ n := Nat.lt_succ_iff.1 (mem_range.1 hk)
    rw [Nat.cast_sub hkn]
    have h1 : ((n : ℝ) - k + 1) ≠ 0 := by
      have : (k : ℝ) ≤ n := by exact_mod_cast hkn
      linarith
    field_simp
    ring
  rw [sum_congr rfl h, ← mul_sum, sum_add_distrib]
  have hr : ∑ k ∈ range (n + 1), 1 / (((n - k : ℕ) : ℝ) + 1)
      = ∑ k ∈ range (n + 1), 1 / ((k : ℝ) + 1) := by
    have := sum_range_reflect (fun k => 1 / ((k : ℝ) + 1)) (n + 1)
    simpa using this
  rw [hr, ← hR_sum]
  field_simp
  ring

/-- `Σ_{k=0}^n H_{k+1}/(k+2) = (H_{n+2}² - H^{(2)}_{n+2})/2`. -/
lemma sumP (n : ℕ) :
    ∑ k ∈ range (n + 1), (harmonic (k + 1) : ℝ) / ((k : ℝ) + 2)
      = ((harmonic (n + 2) : ℝ) ^ 2 - H2 (n + 2)) / 2 := by
  induction n with
  | zero =>
    simp [H2_succ, H2_zero, hR_succ]; norm_num
  | succ n ih =>
    rw [sum_range_succ, ih, show n + 1 + 2 = (n + 2) + 1 by ring, hR_succ (n + 2), H2_succ (n + 2)]
    push_cast
    field_simp
    ring

/-- `Σ_{k=0}^n H_{n+1-k}/(k+1) = H_{n+2}² - H^{(2)}_{n+2}`. -/
lemma sumR (n : ℕ) :
    ∑ k ∈ range (n + 1), (harmonic (n + 1 - k) : ℝ) / ((k : ℝ) + 1)
      = (harmonic (n + 2) : ℝ) ^ 2 - H2 (n + 2) := by
  induction n with
  | zero => simp [H2_succ, H2_zero, hR_succ]; norm_num
  | succ n ih =>
    have h : ∀ k ∈ range (n + 2), (harmonic (n + 1 + 1 - k) : ℝ) / ((k : ℝ) + 1)
        = (harmonic (n + 1 - k) : ℝ) / ((k : ℝ) + 1)
          + 1 / (((k : ℝ) + 1) * (((n + 1 - k : ℕ) : ℝ) + 1)) := by
      intro k hk
      have hkn : k ≤ n + 1 := Nat.lt_succ_iff.1 (mem_range.1 hk)
      rw [show n + 1 + 1 - k = (n + 1 - k) + 1 by omega, hR_succ]
      field_simp
    rw [sum_congr rfl h, sum_add_distrib, sum_range_succ, Nat.sub_self, ih, conv2]
    rw [show n + 1 + 2 = (n + 2) + 1 by ring, hR_succ (n + 2), H2_succ (n + 2)]
    simp only [harmonic_zero, Rat.cast_zero, zero_div, add_zero]
    push_cast
    field_simp
    ring

/-- `Σ_{k=0}^n 2H_{k+1}/((k+2)(n-k+1)) = 3(H_{n+2}² - H^{(2)}_{n+2})/(n+3)`. -/
lemma conv3 (n : ℕ) :
    ∑ k ∈ range (n + 1), 2 * (harmonic (k + 1) : ℝ) / (((k : ℝ) + 2) * (((n - k : ℕ) : ℝ) + 1))
      = 3 * ((harmonic (n + 2) : ℝ) ^ 2 - H2 (n + 2)) / ((n : ℝ) + 3) := by
  have h : ∀ k ∈ range (n + 1),
      2 * (harmonic (k + 1) : ℝ) / (((k : ℝ) + 2) * (((n - k : ℕ) : ℝ) + 1))
      = (2 / ((n : ℝ) + 3)) * ((harmonic (k + 1) : ℝ) / ((k : ℝ) + 2)
        + (harmonic (k + 1) : ℝ) / (((n - k : ℕ) : ℝ) + 1)) := by
    intro k hk
    have hkn : k ≤ n := Nat.lt_succ_iff.1 (mem_range.1 hk)
    rw [Nat.cast_sub hkn]
    have h1 : ((n : ℝ) - k + 1) ≠ 0 := by
      have : (k : ℝ) ≤ n := by exact_mod_cast hkn
      linarith
    field_simp
    ring
  rw [sum_congr rfl h, ← mul_sum, sum_add_distrib, sumP]
  have hr : ∑ k ∈ range (n + 1), (harmonic (k + 1) : ℝ) / (((n - k : ℕ) : ℝ) + 1)
      = ∑ k ∈ range (n + 1), (harmonic (n + 1 - k) : ℝ) / ((k : ℝ) + 1) := by
    rw [← sum_range_reflect]
    apply sum_congr rfl
    intro k hk
    have hkn : k < n + 1 := mem_range.1 hk
    rw [show n + 1 - 1 - k + 1 = n + 1 - k by omega, show n - (n + 1 - 1 - k) = k by omega]
  rw [hr, sumR]
  field_simp
  ring


/-- The series of `log(1-z)`, from index `0`. -/
noncomputable def ser1 (z : ℝ) (n : ℕ) : ℝ := -(z ^ (n + 1) / ((n : ℝ) + 1))
/-- The series of `log²(1-z)`, from index `0` (term of `z^{n+2}`). -/
noncomputable def ser2 (z : ℝ) (n : ℕ) : ℝ :=
  2 * (harmonic (n + 1) : ℝ) / ((n : ℝ) + 2) * z ^ (n + 2)
/-- The series of `log³(1-z)`, from index `0` (term of `z^{n+3}`). -/
noncomputable def ser3 (z : ℝ) (n : ℕ) : ℝ :=
  -3 * ((harmonic (n + 2) : ℝ) ^ 2 - H2 (n + 2)) / ((n : ℝ) + 3) * z ^ (n + 3)

lemma harm_le (n : ℕ) : (harmonic n : ℝ) ≤ n := by
  rw [hR_sum]
  calc ∑ i ∈ range n, 1 / ((i : ℝ) + 1) ≤ ∑ i ∈ range n, (1 : ℝ) := by
        apply sum_le_sum; intro i _
        rw [div_le_one (by positivity)]; linarith [(Nat.cast_nonneg i : (0:ℝ) ≤ i)]
    _ = n := by simp

lemma harm_nonneg (n : ℕ) : (0 : ℝ) ≤ harmonic n := by
  rw [hR_sum]; apply sum_nonneg; intro i _; positivity

lemma geom_summ {z : ℝ} (hz : |z| < 1) (C : ℝ) :
    Summable (fun n : ℕ => C * |z| ^ n) :=
  (summable_geometric_of_lt_one (abs_nonneg z) hz).mul_left C

lemma pow_le_pow_abs {z : ℝ} (hz : |z| < 1) (n m : ℕ) : |z| ^ (n + m) ≤ |z| ^ n := by
  rw [pow_add]; exact mul_le_of_le_one_right (by positivity) (pow_le_one₀ (abs_nonneg z) hz.le)

lemma l1_norm {z : ℝ} (hz : |z| < 1) : Summable (fun n => ‖ser1 z n‖) := by
  refine Summable.of_nonneg_of_le (fun _ => norm_nonneg _) (fun n => ?_) (geom_summ hz 1)
  rw [Real.norm_eq_abs, ser1, abs_neg, abs_div, abs_pow, one_mul, abs_of_pos (by positivity : (0:ℝ) < n + 1),
    div_le_iff₀ (by positivity)]
  have := pow_le_pow_abs hz n 1
  have : (0:ℝ) ≤ |z| ^ n := by positivity
  nlinarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]

lemma l2_norm {z : ℝ} (hz : |z| < 1) : Summable (fun n => ‖ser2 z n‖) := by
  refine Summable.of_nonneg_of_le (fun _ => norm_nonneg _) (fun n => ?_) (geom_summ hz 2)
  have hH := harm_le (n + 1)
  have hH0 := harm_nonneg (n + 1)
  push_cast at hH
  rw [Real.norm_eq_abs, ser2, abs_mul, abs_pow, abs_div, abs_of_pos (by positivity : (0:ℝ) < n + 2),
    abs_mul, abs_of_nonneg hH0, abs_two]
  have h1 : 2 * (harmonic (n + 1) : ℝ) / ((n : ℝ) + 2) ≤ 2 := by
    rw [div_le_iff₀ (by positivity)]; linarith
  have := pow_le_pow_abs hz n 2
  have : (0:ℝ) ≤ |z| ^ (n + 2) := by positivity
  calc 2 * (harmonic (n + 1) : ℝ) / ((n : ℝ) + 2) * |z| ^ (n + 2) ≤ 2 * |z| ^ (n + 2) :=
        mul_le_mul_of_nonneg_right h1 this
    _ ≤ 2 * |z| ^ n := by linarith

lemma hasSum_l1 {z : ℝ} (hz : |z| < 1) : HasSum (ser1 z) (Real.log (1 - z)) := by
  have := (Real.hasSum_pow_div_log_of_abs_lt_one hz).neg
  simp only [neg_neg] at this
  exact this

lemma cauchy (f g : ℕ → ℝ) (F G : ℝ) (hf : HasSum f F) (hg : HasSum g G)
    (hfn : Summable (fun n => ‖f n‖)) (hgn : Summable (fun n => ‖g n‖)) :
    HasSum (fun n => ∑ k ∈ range (n + 1), f k * g (n - k)) (F * G) := by
  have := hasSum_sum_range_mul_of_summable_norm hfn hgn
  rwa [hf.tsum_eq, hg.tsum_eq] at this

lemma hasSum_l2 {z : ℝ} (hz : |z| < 1) : HasSum (ser2 z) (Real.log (1 - z) ^ 2) := by
  have h := cauchy _ _ _ _ (hasSum_l1 hz) (hasSum_l1 hz) (l1_norm hz) (l1_norm hz)
  rw [← sq] at h
  convert h using 1
  funext n
  have e : ∀ k ∈ range (n + 1), ser1 z k * ser1 z (n - k)
      = z ^ (n + 2) * (1 / (((k : ℝ) + 1) * (((n - k : ℕ) : ℝ) + 1))) := by
    intro k hk
    have hkn : k ≤ n := Nat.lt_succ_iff.1 (mem_range.1 hk)
    have hp : z ^ (k + 1) * z ^ (n - k + 1) = z ^ (n + 2) := by
      rw [← pow_add]; congr 1; omega
    unfold ser1
    rw [← hp]
    field_simp
  rw [sum_congr rfl e, ← mul_sum, conv2, ser2]
  ring

lemma hasSum_l3 {z : ℝ} (hz : |z| < 1) : HasSum (ser3 z) (Real.log (1 - z) ^ 3) := by
  have h := cauchy _ _ _ _ (hasSum_l2 hz) (hasSum_l1 hz) (l2_norm hz) (l1_norm hz)
  rw [← pow_succ] at h
  convert h using 1
  funext n
  have e : ∀ k ∈ range (n + 1), ser2 z k * ser1 z (n - k)
      = -z ^ (n + 3) * (2 * (harmonic (k + 1) : ℝ) /
          (((k : ℝ) + 2) * (((n - k : ℕ) : ℝ) + 1))) := by
    intro k hk
    have hkn : k ≤ n := Nat.lt_succ_iff.1 (mem_range.1 hk)
    have hp : z ^ (k + 2) * z ^ (n - k + 1) = z ^ (n + 3) := by
      rw [← pow_add]; congr 1; omega
    unfold ser1 ser2
    rw [← hp]
    field_simp
  rw [sum_congr rfl e, ← mul_sum, conv3, ser3]
  ring

end S05

theorem S05.proof_2 : S05.stmt_2 := by
  intro c hc1 hc3 z hz
  obtain rfl | rfl | rfl : c = 1 ∨ c = 2 ∨ c = 3 := by omega
  · rw [pow_one]; convert S05.hasSum_l1 hz using 1
    funext q
    simp only [S05.lam, S05.ser1]; push_cast; ring
  · rw [← hasSum_nat_add_iff' 1]
    simp only [Finset.range_one, Finset.sum_singleton]
    simp only [S05.lam, Nat.cast_one, zero_add, Nat.sub_self, harmonic_zero, Rat.cast_zero,
      mul_zero, zero_div, zero_mul, sub_zero]
    convert S05.hasSum_l2 hz using 1
    funext q
    simp only [S05.ser2, Nat.add_sub_cancel]; push_cast; ring
  · rw [← hasSum_nat_add_iff' 2]
    simp only [Finset.sum_range_succ, Finset.sum_range_zero]
    have e1 : S05.lam 3 (0 + 1) = 0 := by simp [S05.lam, S05.H2_zero]
    have e2 : S05.lam 3 (1 + 1) = 0 := by
      simp [S05.lam, S05.H2]
    rw [e1, e2]
    simp only [zero_mul, add_zero, sub_zero]
    convert S05.hasSum_l3 hz using 1
    funext q
    simp only [S05.lam, S05.ser3, show q + 2 + 1 - 1 = q + 2 by omega]; push_cast; ring

namespace S05
open Real Set Finset


lemma H2_nonneg (n : ℕ) : 0 ≤ H2 n := by
  unfold H2; apply sum_nonneg; intro i _; positivity

lemma H2_le (n : ℕ) : H2 n ≤ n := by
  induction n with
  | zero => simp [H2_zero]
  | succ n ih =>
    rw [H2_succ]; push_cast
    have : 1 / ((n : ℝ) + 1) ^ 2 ≤ 1 := by
      rw [div_le_one (by positivity)]; nlinarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
    linarith

lemma harm_sq_le (n : ℕ) : (harmonic (n + 1) : ℝ) ^ 2 ≤ 9 * ((n : ℝ) + 1) := by
  have h := Z03.harmonic_le_sqrt n
  have h0 := harm_nonneg (n + 1)
  have hs : (((n : ℝ) + 1) ^ (1 / 2 : ℝ)) ^ 2 = (n : ℝ) + 1 := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul (by positivity)]; norm_num
  nlinarith

/-- Shift trick: if `c (n+1) - c n = d n` and `Σ d_n x^{n+1} = D`, `c 0 = 0`, then
`Σ c_{n+1} x^{n+1} = D/(1-x)`. -/
lemma shift_trick (c d : ℕ → ℝ) (x D : ℝ) (hx : |x| < 1) (hc0 : c 0 = 0)
    (hcd : ∀ n, c (n + 1) - c n = d n)
    (hs : Summable (fun n : ℕ => c (n + 1) * x ^ (n + 1)))
    (hD : HasSum (fun n : ℕ => d n * x ^ (n + 1)) D) :
    HasSum (fun n : ℕ => c (n + 1) * x ^ (n + 1)) (D / (1 - x)) := by
  have hs' := hs.hasSum
  set s := ∑' n : ℕ, c (n + 1) * x ^ (n + 1)
  have h1 : HasSum (fun n : ℕ => c n * x ^ (n + 1)) (x * s) := by
    rw [← hasSum_nat_add_iff' 1]
    simp only [Finset.range_one, Finset.sum_singleton, hc0, zero_mul, sub_zero]
    convert hs'.mul_left x using 1
    funext n; ring
  have h2 : HasSum (fun n : ℕ => d n * x ^ (n + 1)) (s - x * s) := by
    convert hs'.sub h1 using 1
    funext n; rw [← hcd]; ring
  have h4 := h2.unique hD
  have hx1 : 1 - x ≠ 0 := by linarith [le_abs_self x]
  rw [show D / (1 - x) = s by rw [← h4]; field_simp]
  exact hs'

lemma summ_bound (c : ℕ → ℝ) (x : ℝ) (hx : |x| < 1) (C : ℝ)
    (hc : ∀ n, |c n| ≤ C * ((n : ℝ) + 1)) :
    Summable (fun n : ℕ => c n * x ^ (n + 1)) := by
  apply Summable.of_norm_bounded ((Z03.summable_geom_mul (abs_nonneg x) hx).mul_left C)
  intro n
  rw [Real.norm_eq_abs, abs_mul, abs_pow]
  have h1 : |x| ^ (n + 1) ≤ |x| ^ n := by
    rw [pow_succ]; exact mul_le_of_le_one_right (by positivity) hx.le
  calc |c n| * |x| ^ (n + 1) ≤ (C * ((n : ℝ) + 1)) * |x| ^ n :=
        mul_le_mul (hc n) h1 (by positivity) (le_trans (abs_nonneg _) (hc n))
    _ = C * (((n : ℝ) + 1) * |x| ^ n) := by ring

lemma euler1 {x : ℝ} (hx : |x| < 1) :
    HasSum (fun n : ℕ => (harmonic (n + 1) : ℝ) * x ^ (n + 1)) (cL1 x / (1 - x)) :=
  Z03.hasSum_C hx

lemma euler2 {x : ℝ} (hx : |x| < 1) :
    HasSum (fun n : ℕ => H2 (n + 1) * x ^ (n + 1)) (D00.Li2 x / (1 - x)) := by
  refine shift_trick H2 (fun n => 1 / ((n : ℝ) + 1) ^ 2) x _ hx H2_zero
    (fun n => by rw [H2_succ]; ring) ?_ ?_
  · apply summ_bound (fun n => H2 (n + 1)) x hx 1
    intro n; rw [abs_of_nonneg (H2_nonneg _)]; have := H2_le (n + 1); push_cast at this; linarith
  · convert D00.hasSum_Li2 hx.le using 1; funext n; ring

lemma euler4 {x : ℝ} (hx : |x| < 1) :
    HasSum (fun n : ℕ => (harmonic (n + 1) : ℝ) / ((n : ℝ) + 1) * x ^ (n + 1))
      (D00.Li2 x + (1 / 2) * cL1 x ^ 2) := by
  have h := (Z03.summable_B hx).hasSum
  have e : (∑' n : ℕ, (harmonic (n + 1) : ℝ) * x ^ (n + 1) / ((n : ℝ) + 1)) = Z03.B x := rfl
  rw [e, Z03.B_eq hx] at h
  unfold cL1; rw [neg_sq]
  convert h using 1; funext n; ring

lemma euler3 {x : ℝ} (hx : |x| < 1) :
    HasSum (fun n : ℕ => (harmonic (n + 1) : ℝ) ^ 2 * x ^ (n + 1))
      ((cL1 x ^ 2 + D00.Li2 x) / (1 - x)) := by
  refine shift_trick (fun n => (harmonic n : ℝ) ^ 2)
    (fun n => 2 * ((harmonic (n + 1) : ℝ) / ((n : ℝ) + 1)) - 1 / ((n : ℝ) + 1) ^ 2) x _ hx
    (by simp) (fun n => by rw [hR_succ]; field_simp; ring) ?_ ?_
  · apply summ_bound (fun n => (harmonic (n + 1) : ℝ) ^ 2) x hx 9
    intro n; rw [abs_of_nonneg (sq_nonneg _)]; exact harm_sq_le n
  · have := ((euler4 hx).mul_left 2).sub (D00.hasSum_Li2 hx.le)
    convert this using 1
    · funext n; ring
    · ring


lemma hasDerivAt_cL1 {y : ℝ} (hy : y < 1) : HasDerivAt cL1 (1 / (1 - y)) y := by
  have hne : 1 - y ≠ 0 := by linarith
  have := (((hasDerivAt_id' y).const_sub 1).log hne).fun_neg
  rw [show (1 / (1 - y)) = -(-1 / (1 - y)) by ring]
  exact this

lemma continuousAt_cL1 {y : ℝ} (hy : y < 1) : ContinuousAt cL1 y :=
  (hasDerivAt_cL1 hy).continuousAt

lemma abs_lt_of_uIcc {x z : ℝ} (hx : |x| < 1) (hz : z ∈ uIcc 0 x) : |z| < 1 :=
  lt_of_le_of_lt (Z03.abs_le_of_mem_uIcc hz) hx

/-- Generic uniqueness: power series `F = Σ a_n z^{n+1}/(n+1)` whose derivative series equals
`g' ` and `g` with `g(0) = 0`. -/
lemma powser_eq (a : ℕ → ℝ) (ha : ∀ n, |a n| ≤ (n : ℝ) + 1) (g g' : ℝ → ℝ) {x : ℝ} (hx : |x| < 1)
    (hg0 : g 0 = 0) (hgc : ∀ z, |z| < 1 → ContinuousAt g z)
    (hg : ∀ z, |z| < 1 → z ≠ 0 → HasDerivAt g (g' z) z)
    (hser : ∀ z, |z| < 1 → z ≠ 0 → HasSum (fun n => a n * z ^ n) (g' z)) :
    HasSum (fun n => a n * x ^ (n + 1) / ((n : ℝ) + 1)) (g x) := by
  set F := fun z => ∑' n : ℕ, a n * z ^ (n + 1) / ((n : ℝ) + 1) with hF
  have hsum : Summable (fun n => a n * x ^ (n + 1) / ((n : ℝ) + 1)) := by
    have := summ_bound (fun n => a n / ((n : ℝ) + 1)) x hx 1 (fun n => by
      rw [abs_div, abs_of_pos (by positivity : (0:ℝ) < n + 1), div_le_iff₀ (by positivity)]
      nlinarith [ha n, (Nat.cast_nonneg n : (0:ℝ) ≤ n)])
    convert this using 1; funext n; ring
  have key := D00.eq_of_deriv_zero' (f := fun z => F z - g z) (a := 0) (b := x)
    (by
      intro z hz
      have hz1 := abs_lt_of_uIcc hx hz
      exact ((Z03.hasDerivAt_powser a ha hz1).continuousAt.sub (hgc z hz1)).continuousWithinAt)
    (by
      intro z hz
      obtain ⟨hz0, hzx⟩ := D00.mem_uIoo_zero hz
      have hz1 : |z| < 1 := lt_trans hzx hx
      have h1 := Z03.hasDerivAt_powser a ha hz1
      rw [(hser z hz1 hz0).tsum_eq] at h1
      have := h1.fun_sub (hg z hz1 hz0)
      rw [sub_self] at this
      exact this)
  have hF0 : F 0 = 0 := by simp [hF]
  simp only [hF0, hg0, sub_zero] at key
  have : F x = g x := by linarith
  rw [← this]
  exact hsum.hasSum

lemma euler5 {x : ℝ} (hx : |x| < 1) :
    HasSum (fun n : ℕ => (harmonic (n + 1) : ℝ) ^ 2 / ((n : ℝ) + 1) * x ^ (n + 1))
      (Z03.Li3 x + cL1 x * D00.Li2 x + (1 / 3) * cL1 x ^ 3) := by
  have h := powser_eq (fun n => (harmonic (n + 1) : ℝ) ^ 2 / 9)
    (fun n => by
      rw [abs_of_nonneg (by positivity)]; have := harm_sq_le n; linarith)
    (fun z => (Z03.Li3 z + cL1 z * D00.Li2 z + (1 / 3) * cL1 z ^ 3) / 9)
    (fun z => (cL1 z ^ 2 + D00.Li2 z) / (z * (1 - z)) / 9) hx
    (by simp [cL1])
    (fun z hz => by
      have hz1 : z < 1 := lt_of_le_of_lt (le_abs_self z) hz
      exact ((Z03.continuous_Li3.continuousAt.add ((continuousAt_cL1 hz1).mul
        D00.continuous_Li2.continuousAt)).add (continuousAt_const.mul
        ((continuousAt_cL1 hz1).pow 3))).div_const 9)
    (fun z hz hz0 => by
      have hz1 : z < 1 := lt_of_le_of_lt (le_abs_self z) hz
      have hne : 1 - z ≠ 0 := by linarith
      have := (((Z03.hasDerivAt_Li3 z).add ((hasDerivAt_cL1 hz1).mul
        (D00.hasDerivAt_Li2 hz1 hz0))).add (((hasDerivAt_cL1 hz1).pow 3).const_mul (1 / 3))).div_const 9
      convert this using 1
      rw [Z03.H_of_ne hz0]; unfold cL1; field_simp; ring)
    (fun z hz hz0 => by
      have hne : 1 - z ≠ 0 := by linarith [le_abs_self z]
      have := ((euler3 hz).div_const (9 * z))
      convert this using 1
      · funext n; rw [pow_succ]; field_simp; ring
      · field_simp)
  have := h.mul_left 9
  convert this using 1
  · funext n; field_simp
  · field_simp

lemma euler6 {x : ℝ} (hx : |x| < 1) :
    HasSum (fun n : ℕ => H2 (n + 1) / ((n : ℝ) + 1) * x ^ (n + 1))
      (3 * Z03.Li3 x - 2 * Z03.A x + cL1 x * D00.Li2 x) := by
  have hAc : ∀ z, |z| < 1 → ContinuousAt Z03.A z := by
    intro z hz
    apply Z03.continuousOn_A.continuousAt
    apply Icc_mem_nhds <;> linarith [neg_abs_le z, le_abs_self z]
  have h := powser_eq (fun n => H2 (n + 1))
    (fun n => by
      rw [abs_of_nonneg (H2_nonneg _)]; have := H2_le (n + 1); push_cast at this; linarith)
    (fun z => 3 * Z03.Li3 z - 2 * Z03.A z + cL1 z * D00.Li2 z)
    (fun z => D00.Li2 z / (z * (1 - z))) hx
    (by simp [cL1, Z03.A, Z03.Aterm])
    (fun z hz => by
      have hz1 : z < 1 := lt_of_le_of_lt (le_abs_self z) hz
      exact ((continuousAt_const.mul Z03.continuous_Li3.continuousAt).sub
        (continuousAt_const.mul (hAc z hz))).add ((continuousAt_cL1 hz1).mul
        D00.continuous_Li2.continuousAt))
    (fun z hz hz0 => by
      have hz1 : z < 1 := lt_of_le_of_lt (le_abs_self z) hz
      have hne : 1 - z ≠ 0 := by linarith
      have := ((((Z03.hasDerivAt_Li3 z).const_mul 3).sub
        ((Z03.hasDerivAt_A hz hz0).const_mul 2)).add ((hasDerivAt_cL1 hz1).mul
        (D00.hasDerivAt_Li2 hz1 hz0)))
      convert this using 1
      rw [Z03.H_of_ne hz0]; unfold cL1; field_simp; ring)
    (fun z hz hz0 => by
      have hne : 1 - z ≠ 0 := by linarith [le_abs_self z]
      have := ((euler2 hz).div_const z)
      convert this using 1
      · funext n; rw [pow_succ]; field_simp
      · field_simp)
  convert h using 1
  funext n; ring

end S05

theorem S05.proof_1 : S05.stmt_1 := fun _ hx =>
  ⟨S05.euler1 hx, S05.euler2 hx, S05.euler3 hx, S05.euler4 hx, S05.euler5 hx, S05.euler6 hx⟩

namespace S05
open Real Set Finset


lemma dot_add_smul (u P Q : X04.Plane) (a b : ℝ) :
    X04.dot u (a • P + b • Q) = a * X04.dot u P + b * X04.dot u Q := by
  simp only [X04.dot, PiLp.add_apply, PiLp.smul_apply, smul_eq_mul]; ring

lemma dot_neg' (u x : X04.Plane) : X04.dot u (-x) = -X04.dot u x := by
  simp only [X04.dot, PiLp.neg_apply]; ring

lemma abs_comb2 (a b α β : ℝ) (hab : 0 ≤ a * b) (hα : 0 ≤ α) (hβ : 0 ≤ β) :
    |α * a + β * b| = α * |a| + β * |b| := by
  rcases le_total 0 a with ha | ha <;> rcases le_total 0 b with hb | hb
  · rw [abs_of_nonneg ha, abs_of_nonneg hb, abs_of_nonneg (by positivity)]
  · have hab0 : a * b = 0 := le_antisymm (mul_nonpos_of_nonneg_of_nonpos ha hb) hab
    rcases mul_eq_zero.1 hab0 with h | h <;> subst h <;> simp [abs_mul, abs_of_nonneg hα, abs_of_nonneg hβ]
  · have hab0 : a * b = 0 := le_antisymm (mul_nonpos_of_nonpos_of_nonneg ha hb) hab
    rcases mul_eq_zero.1 hab0 with h | h <;> subst h <;> simp [abs_mul, abs_of_nonneg hα, abs_of_nonneg hβ]
  · rw [abs_of_nonpos ha, abs_of_nonpos hb, abs_of_nonpos (by nlinarith)]; ring

lemma Fw_neg (u : Fin 3 → X04.Plane) (w : Fin 3 → ℝ) (x : X04.Plane) :
    Fw u w 1 (-x) = Fw u w 1 x := by
  simp only [Fw_one, dot_neg', abs_neg]

lemma Fw_comb (u : Fin 3 → X04.Plane) (w : Fin 3 → ℝ) (P Q : X04.Plane) (α β : ℝ)
    (hα : 0 ≤ α) (hβ : 0 ≤ β) (hs : ∀ j, 0 ≤ X04.dot (u j) P * X04.dot (u j) Q) :
    Fw u w 1 (α • P + β • Q) = α * Fw u w 1 P + β * Fw u w 1 Q := by
  simp only [Fw_one, dot_add_smul, abs_comb2 _ _ _ _ (hs _) hα hβ]; ring

lemma cramer (P Q x : X04.Plane) (h : det2 P Q ≠ 0) :
    x = (det2 x Q / det2 P Q) • P + (det2 P x / det2 P Q) • Q := by
  unfold det2 at h ⊢
  have h' : Q 1 * P 0 - Q 0 * P 1 ≠ 0 := by intro h2; apply h; linarith
  ext i
  fin_cases i <;> simp <;> field_simp <;> ring

lemma sector (u : Fin 3 → X04.Plane) (w : Fin 3 → ℝ) (S : Set X04.Plane) (hS : Convex ℝ S)
    (P Q x : X04.Plane) (hPQ : 0 < det2 P Q) (hx1 : 0 ≤ det2 P x) (hx2 : 0 ≤ det2 x Q)
    (hs : ∀ j, 0 ≤ X04.dot (u j) P * X04.dot (u j) Q) (hFP : Fw u w 1 P = 1)
    (hFQ : Fw u w 1 Q = 1) (hFx : Fw u w 1 x ≤ 1) (hP : P ∈ S) (hQ : Q ∈ S) (h0 : 0 ∈ S) :
    x ∈ S := by
  set α := det2 x Q / det2 P Q
  set β := det2 P x / det2 P Q
  have hα : 0 ≤ α := div_nonneg hx2 hPQ.le
  have hβ : 0 ≤ β := div_nonneg hx1 hPQ.le
  have hx := cramer P Q x hPQ.ne'
  have hsum : α + β ≤ 1 := by
    have := Fw_comb u w P Q α β hα hβ hs
    rw [← hx, hFP, hFQ] at this; linarith
  -- x = α P + β Q + (1 - α - β) 0
  have hseg : α • P + β • Q ∈ S := by
    by_cases hab : α + β = 0
    · have ha0 : α = 0 := by linarith
      have hb0 : β = 0 := by linarith
      simp [ha0, hb0, h0]
    · have hpos : 0 < α + β := lt_of_le_of_ne (by linarith) (Ne.symm hab)
      have hm : (α / (α + β)) • P + (β / (α + β)) • Q ∈ S :=
        hS hP hQ (div_nonneg hα hpos.le) (div_nonneg hβ hpos.le) (by field_simp)
      have := hS hm h0 hpos.le (by linarith : 0 ≤ 1 - (α + β)) (by ring)
      convert this using 1
      simp only [smul_add, smul_smul, smul_zero, add_zero]
      congr 2 <;> field_simp
  rw [hx]; exact hseg

end S05

namespace S05

lemma convex_Kw (u : Fin 3 → X04.Plane) (w : Fin 3 → ℝ) (hw : ∀ j, 0 < w j) :
    Convex ℝ (Kw u w 1) := by
  intro x hx y hy a b ha hb hab
  simp only [Kw, Set.mem_ofPred_eq] at hx hy ⊢
  rw [Fw_one] at hx hy ⊢
  simp only [dot_add_smul]
  have e : ∀ j, w j * |a * X04.dot (u j) x + b * X04.dot (u j) y|
      ≤ a * (w j * |X04.dot (u j) x|) + b * (w j * |X04.dot (u j) y|) := by
    intro j
    have := abs_add_le (a * X04.dot (u j) x) (b * X04.dot (u j) y)
    rw [abs_mul, abs_mul, abs_of_nonneg ha, abs_of_nonneg hb] at this
    nlinarith [hw j]
  nlinarith [e 0, e 1, e 2]

end S05

theorem S05.proof_6 : S05.stmt_6 := by
  intro H hH
  obtain ⟨hw, hdet, hker, hF, h01, h12, h02⟩ := hH
  set S := convexHull ℝ (Set.range H.C ∪ Set.range (fun j => -H.C j)) with hSdef
  have hS : Convex ℝ S := convex_convexHull ℝ _
  have hC : ∀ j, H.C j ∈ S := fun j => subset_convexHull ℝ _ (Or.inl ⟨j, rfl⟩)
  have hnC : ∀ j, -H.C j ∈ S := fun j => subset_convexHull ℝ _ (Or.inr ⟨j, rfl⟩)
  have h0 : (0 : X04.Plane) ∈ S := by
    have := hS (hC 0) (hnC 0) (by norm_num : (0:ℝ) ≤ 1 / 2) (by norm_num : (0:ℝ) ≤ 1 / 2)
      (by norm_num)
    simpa using this
  have dsw : ∀ x z : X04.Plane, S05.det2 x z = -S05.det2 z x := by
    intro x z; simp [S05.det2]; ring
  have dn1 : ∀ x z : X04.Plane, S05.det2 (-x) z = -S05.det2 x z := by
    intro x z; simp [S05.det2]; ring
  have dn2 : ∀ x z : X04.Plane, S05.det2 x (-z) = -S05.det2 x z := by
    intro x z; simp [S05.det2]; ring
  have hFn : ∀ j, S05.Fw H.u H.w 1 (-H.C j) = 1 := fun j => by rw [S05.Fw_neg]; exact hF j
  apply Set.Subset.antisymm
  · intro x hx
    simp only [S05.Kw, Set.mem_ofPred_eq] at hx
    -- sign facts on the three sector pairs
    have s01 : ∀ j, 0 ≤ X04.dot (H.u j) (H.C 0) * X04.dot (H.u j) (H.C 1) := by
      intro j
      obtain rfl | rfl | rfl : j = 0 ∨ j = 1 ∨ j = 2 := by fin_cases j <;> simp
      · rw [hker 0, zero_mul]
      · rw [hker 1, mul_zero]
      · exact S05.sameSign _ _ _ _ (hker 2) (by rw [dsw (H.C 2), dsw (H.C 2)]; nlinarith)
    have s12 : ∀ j, 0 ≤ X04.dot (H.u j) (H.C 1) * X04.dot (H.u j) (H.C 2) := by
      intro j
      obtain rfl | rfl | rfl : j = 0 ∨ j = 1 ∨ j = 2 := by fin_cases j <;> simp
      · exact S05.sameSign _ _ _ _ (hker 0) (by nlinarith)
      · rw [hker 1, zero_mul]
      · rw [hker 2, mul_zero]
    have s20 : ∀ j, 0 ≤ X04.dot (H.u j) (H.C 2) * X04.dot (H.u j) (-H.C 0) := by
      intro j
      obtain rfl | rfl | rfl : j = 0 ∨ j = 1 ∨ j = 2 := by fin_cases j <;> simp
      · rw [S05.dot_neg', hker 0, neg_zero, mul_zero]
      · exact S05.sameSign _ _ _ _ (hker 1) (by rw [dn2, dsw (H.C 1) (H.C 0)]; nlinarith)
      · rw [hker 2, zero_mul]
    have neg2 : ∀ (P Q : X04.Plane), (∀ j, 0 ≤ X04.dot (H.u j) P * X04.dot (H.u j) Q) →
        ∀ j, 0 ≤ X04.dot (H.u j) (-P) * X04.dot (H.u j) (-Q) := by
      intro P Q h j; rw [S05.dot_neg', S05.dot_neg', neg_mul_neg]; exact h j
    have s20' : ∀ j, 0 ≤ X04.dot (H.u j) (-H.C 2) * X04.dot (H.u j) (H.C 0) := by
      intro j; have := neg2 _ _ s20 j; rwa [neg_neg] at this
    set d0 := S05.det2 (H.C 0) x
    set d1 := S05.det2 (H.C 1) x
    set d2 := S05.det2 (H.C 2) x
    rcases le_total 0 d0 with g0 | g0 <;> rcases le_total 0 d1 with g1 | g1 <;>
      rcases le_total 0 d2 with g2 | g2
    · -- S3
      exact S05.sector H.u H.w S hS _ _ x (by rw [dn2, dsw]; linarith) g2
        (by rw [dn2, dsw]; linarith) s20 (hF 2) (hFn 0) hx (hC 2) (hnC 0) h0
    · -- S2
      exact S05.sector H.u H.w S hS _ _ x h12 g1 (by rw [dsw]; linarith) s12 (hF 1) (hF 2) hx
        (hC 1) (hC 2) h0
    · exact S05.sector H.u H.w S hS _ _ x h01 g0 (by rw [dsw]; linarith) s01 (hF 0) (hF 1) hx
        (hC 0) (hC 1) h0
    · exact S05.sector H.u H.w S hS _ _ x h01 g0 (by rw [dsw]; linarith) s01 (hF 0) (hF 1) hx
        (hC 0) (hC 1) h0
    · -- S4
      exact S05.sector H.u H.w S hS _ _ x (by rw [dn1, dn2, neg_neg]; exact h01)
        (by rw [dn1]; linarith) (by rw [dn2, dsw, neg_neg]; exact g1) (neg2 _ _ s01)
        (hFn 0) (hFn 1) hx (hnC 0) (hnC 1) h0
    · exact S05.sector H.u H.w S hS _ _ x h12 g1 (by rw [dsw]; linarith) s12 (hF 1) (hF 2) hx
        (hC 1) (hC 2) h0
    · -- S5
      exact S05.sector H.u H.w S hS _ _ x (by rw [dn1, dn2, neg_neg]; exact h12)
        (by rw [dn1]; linarith) (by rw [dn2, dsw, neg_neg]; exact g2) (neg2 _ _ s12)
        (hFn 1) (hFn 2) hx (hnC 1) (hnC 2) h0
    · -- S6
      exact S05.sector H.u H.w S hS _ _ x (by rw [dn1, dsw]; linarith)
        (by rw [dn1]; linarith) (by rw [dsw]; linarith) s20'
        (hFn 2) (hF 0) hx (hnC 2) (hC 0) h0
  · apply convexHull_min _ (S05.convex_Kw H.u H.w hw)
    rintro y (⟨j, rfl⟩ | ⟨j, rfl⟩)
    · simp only [S05.Kw, Set.mem_ofPred_eq, hF j, le_refl]
    · simp only [S05.Kw, Set.mem_ofPred_eq, hFn j, le_refl]

namespace S05
open MeasureTheory Set Filter Real Interval Topology

lemma ae_Ioo01 : ∀ᵐ t ∂(volume.restrict (Ι (0:ℝ) 1)), t ∈ Ioo (0:ℝ) 1 := by
  rw [uIoc_of_le zero_le_one, ← Measure.restrict_congr_set Ioo_ae_eq_Ioc]
  exact ae_restrict_mem measurableSet_Ioo

lemma ii_tt (c : ℝ) (hc : -1 < c) :
    IntervalIntegrable (fun t : ℝ => (t * (1 - t)) ^ c) volume 0 1 := by
  have h1 : IntervalIntegrable (fun t : ℝ => (t * (1 - t)) ^ c) volume 0 (1/2) := by
    have : IntervalIntegrable (fun t : ℝ => t ^ c * (1 - t) ^ c) volume 0 (1/2) := by
      apply IntervalIntegrable.mul_continuousOn (intervalIntegral.intervalIntegrable_rpow' hc)
      apply ContinuousOn.rpow_const (by fun_prop)
      intro x hx; left; rw [uIcc_of_le (by norm_num)] at hx; linarith [hx.2]
    refine this.congr ?_
    intro t ht
    rw [uIoc_of_le (by norm_num)] at ht
    exact (mul_rpow ht.1.le (by linarith [ht.2])).symm
  have h2 := (h1.comp_sub_left 1).symm
  refine h1.trans ?_
  convert h2 using 1
  · funext t; ring_nf
  · norm_num
  · norm_num

lemma abs_log_le {x θ : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) (hθ : 0 < θ) :
    |log x| ≤ x ^ (-θ) / θ := by
  have := abs_log_mul_self_rpow_lt x θ h0 h1 hθ
  have hp := rpow_pos_of_pos h0 θ
  rw [abs_mul, abs_of_pos hp] at this
  rw [rpow_neg h0.le, show (x ^ θ)⁻¹ / θ = (1 / θ) / x ^ θ by field_simp, le_div_iff₀ hp]
  linarith

lemma logsum_le {t θ : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) (hθ : 0 < θ) :
    |log t| + |log (1 - t)| ≤ 2 / θ * (t * (1 - t)) ^ (-θ) := by
  obtain ⟨h0, h1⟩ := ht
  have hu : 0 < 1 - t := by linarith
  have hp : 0 < t * (1 - t) := mul_pos h0 hu
  have e1 := abs_log_le h0 h1.le hθ
  have e2 := abs_log_le hu (by linarith) hθ
  have f1 : t ^ (-θ) ≤ (t * (1 - t)) ^ (-θ) :=
    rpow_le_rpow_of_nonpos hp (by nlinarith) (by linarith)
  have f2 : (1 - t) ^ (-θ) ≤ (t * (1 - t)) ^ (-θ) :=
    rpow_le_rpow_of_nonpos hp (by nlinarith) (by linarith)
  have : t ^ (-θ) / θ + (1 - t) ^ (-θ) / θ ≤ 2 / θ * (t * (1 - t)) ^ (-θ) := by
    rw [← add_div, div_le_iff₀ hθ]
    have : 2 / θ * (t * (1 - t)) ^ (-θ) * θ = 2 * (t * (1 - t)) ^ (-θ) := by field_simp
    rw [this]; linarith
  linarith

/-- the log-sum `s(t) = |log t| + |log(1-t)|` -/
noncomputable def ls (t : ℝ) : ℝ := |log t| + |log (1 - t)|

lemma ls_nonneg (t : ℝ) : 0 ≤ ls t := by unfold ls; positivity

lemma ls_pow_le {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) (m : ℕ) :
    ls t ^ m ≤ (4 * (m + 1)) ^ m * (t * (1 - t)) ^ (-(1/2:ℝ)) := by
  have hp : 0 < t * (1 - t) := mul_pos ht.1 (by linarith [ht.2])
  have hθ : (0:ℝ) < 1 / (2 * (m + 1)) := by positivity
  have h := logsum_le ht hθ
  have h' : ls t ^ m ≤ (2 / (1 / (2 * (m + 1))) * (t * (1 - t)) ^ (-(1 / (2 * ((m:ℝ) + 1))))) ^ m :=
    pow_le_pow_left₀ (ls_nonneg t) h m
  rw [mul_pow, ← rpow_mul_natCast hp.le] at h'
  have hle1 : t * (1 - t) ≤ 1 := by nlinarith [ht.1, ht.2]
  have : (t * (1 - t)) ^ (-(1 / (2 * ((m:ℝ) + 1))) * (m:ℝ)) ≤ (t * (1 - t)) ^ (-(1/2:ℝ)) := by
    apply rpow_le_rpow_of_exponent_ge hp hle1
    rw [neg_mul, neg_le_neg_iff, div_mul_eq_mul_div, one_mul, div_le_iff₀ (by positivity)]
    linarith
  have e : (2 / (1 / (2 * ((m:ℝ) + 1)))) = 4 * (m + 1) := by field_simp; ring
  rw [e] at h'
  calc ls t ^ m ≤ _ := h'
    _ ≤ _ := by gcongr

/-- Integrability of anything bounded by `K * s(t)^m` -/
lemma ii_of_ls {f : ℝ → ℝ} (hf : AEStronglyMeasurable f (volume.restrict (Ι (0:ℝ) 1)))
    (K : ℝ) (m : ℕ) (h : ∀ t ∈ Ioo (0:ℝ) 1, |f t| ≤ K * ls t ^ m) :
    IntervalIntegrable f volume 0 1 := by
  refine IntervalIntegrable.mono_fun'
    (((ii_tt (-(1/2)) (by norm_num)).const_mul (|K| * (4 * (m + 1)) ^ m))) hf ?_
  filter_upwards [ae_Ioo01] with t ht
  have h1 := h t ht
  have h2 := ls_pow_le ht m
  have hK : K ≤ |K| := le_abs_self K
  have hl := pow_nonneg (ls_nonneg t) m
  simp only [Real.norm_eq_abs]
  calc |f t| ≤ K * ls t ^ m := h1
    _ ≤ |K| * ls t ^ m := by gcongr
    _ ≤ |K| * ((4 * (m + 1)) ^ m * (t * (1 - t)) ^ (-(1/2:ℝ))) := by gcongr
    _ = _ := by ring


end S05

namespace S05
open MeasureTheory Set Filter Real Interval Topology Finset

lemma H2_le' (n : ℕ) : H2 n ≤ n := by
  induction n with
  | zero => simp [H2_zero]
  | succ n ih =>
    rw [H2_succ]; push_cast
    have : 1 / ((n : ℝ) + 1) ^ 2 ≤ 1 := by
      rw [div_le_one (by positivity)]; nlinarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
    linarith

lemma H2_nonneg' (n : ℕ) : 0 ≤ H2 n := by
  unfold H2; apply Finset.sum_nonneg; intro i _; positivity

lemma lam_bound (c q : ℕ) : |lam c (q + 1)| ≤ 3 * ((q : ℝ) + 1) := by
  have hq : (0:ℝ) < (q : ℝ) + 1 := by positivity
  have hH := harm_le q
  have hH0 := harm_nonneg q
  have h2 := H2_le' q
  have h20 := H2_nonneg' q
  have hq0 : (0:ℝ) ≤ q := Nat.cast_nonneg q
  match c with
  | 0 => simp [lam]; positivity
  | 1 =>
    simp only [lam]; push_cast
    rw [abs_div, abs_neg, abs_one, abs_of_pos hq, div_le_iff₀ hq]; nlinarith
  | 2 =>
    simp only [lam, Nat.add_sub_cancel]; push_cast
    rw [abs_div, abs_of_pos hq, div_le_iff₀ hq, abs_mul, abs_two, abs_of_nonneg hH0]; nlinarith
  | 3 =>
    simp only [lam, Nat.add_sub_cancel]; push_cast
    rw [abs_div, abs_of_pos hq, div_le_iff₀ hq, abs_mul, abs_neg, abs_of_pos (by norm_num : (0:ℝ) < 3)]
    have : |(harmonic q : ℝ) ^ 2 - H2 q| ≤ (q : ℝ) ^ 2 + q := by
      rw [abs_le]; constructor <;> nlinarith
    nlinarith
  | c + 4 => simp [lam]; positivity

lemma abs_log_le_sum {y : ℝ} (hy : 0 < y) : |Real.log y| ≤ y + 1 / y := by
  rw [abs_le]; constructor
  · have := Real.log_le_sub_one_of_pos (inv_pos.2 hy)
    rw [Real.log_inv] at this
    have : 0 < y := hy
    rw [one_div]; linarith
  · have := Real.log_le_sub_one_of_pos hy
    have : 0 < 1 / y := by positivity
    linarith

lemma abs_log_lin_le {x t : ℝ} (hx : |x| < 1) (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    |Real.log (1 - x * t)| ≤ 2 + 1 / (1 - |x|) := by
  have hxt : |x * t| ≤ |x| := by
    rw [abs_mul, abs_of_nonneg ht0]; exact mul_le_of_le_one_right (abs_nonneg _) ht1
  have h1 : 1 - |x| ≤ 1 - x * t := by linarith [le_abs_self (x * t)]
  have h2 : 1 - x * t ≤ 2 := by linarith [neg_abs_le (x * t), abs_nonneg x]
  have hp : 0 < 1 - |x| := by linarith
  have hp2 : 0 < 1 - x * t := by linarith
  have := abs_log_le_sum hp2
  have : 1 / (1 - x * t) ≤ 1 / (1 - |x|) := one_div_le_one_div_of_le hp h1
  linarith

lemma abs_log_le_ls (t : ℝ) : |Real.log t| ≤ ls t := by unfold ls; linarith [abs_nonneg (Real.log (1 - t))]
lemma abs_log1_le_ls (t : ℝ) : |Real.log (1 - t)| ≤ ls t := by unfold ls; linarith [abs_nonneg (Real.log t)]

lemma g_bound {n a b : ℕ} {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    |t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b| ≤ 1 * ls t ^ (a + b) := by
  rw [abs_mul, abs_mul, abs_pow, abs_pow, abs_pow, one_mul, pow_add]
  have h1 : |t| ^ n ≤ 1 := pow_le_one₀ (abs_nonneg t) (by rw [abs_of_pos ht.1]; exact ht.2.le)
  have a1 := pow_le_pow_left₀ (abs_nonneg _) (abs_log_le_ls t) a
  have a2 := pow_le_pow_left₀ (abs_nonneg _) (abs_log1_le_ls t) b
  have := ls_nonneg t
  have p1 : 0 ≤ |Real.log t| ^ a := by positivity
  have p2 : 0 ≤ |Real.log (1 - t)| ^ b := by positivity
  calc |t| ^ n * |Real.log t| ^ a * |Real.log (1 - t)| ^ b
      ≤ 1 * |Real.log t| ^ a * |Real.log (1 - t)| ^ b := by gcongr
    _ ≤ 1 * ls t ^ a * ls t ^ b := by gcongr
    _ = ls t ^ a * ls t ^ b := by ring

lemma ii_g (n a b : ℕ) :
    IntervalIntegrable (fun t : ℝ => t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b) volume 0 1 :=
  ii_of_ls (by fun_prop : Measurable (fun t : ℝ => t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b)).aestronglyMeasurable
    1 (a + b) (fun t ht => g_bound ht)

lemma ii_J (x : ℝ) (hx : |x| < 1) (n a b c : ℕ) :
    IntervalIntegrable
        (fun t : ℝ => t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b * Real.log (1 - x * t) ^ c)
        volume 0 1 := by
  refine ii_of_ls (by fun_prop : Measurable (fun t : ℝ =>
      t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b * Real.log (1 - x * t) ^ c)).aestronglyMeasurable
    ((2 + 1 / (1 - |x|)) ^ c) (a + b) (fun t ht => ?_)
  rw [abs_mul, abs_pow]
  have h1 := g_bound (n := n) (a := a) (b := b) ht
  have h2 := pow_le_pow_left₀ (abs_nonneg _) (abs_log_lin_le hx ht.1.le ht.2.le) c
  have := pow_nonneg (ls_nonneg t) (a + b)
  calc |t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b| * |Real.log (1 - x * t)| ^ c
      ≤ (1 * ls t ^ (a + b)) * (2 + 1 / (1 - |x|)) ^ c := by gcongr
    _ = _ := by ring

end S05

theorem S05.proof_3 : S05.stmt_3 := by
  intro x hx n a b c hc1 hc3
  refine ⟨S05.ii_J x hx n a b c, ?_⟩
  set g := fun t : ℝ => t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b with hg
  have hgi : IntervalIntegrable g MeasureTheory.volume 0 1 := S05.ii_g n a b
  have hS : Summable (fun q : ℕ => |S05.lam c (q + 1)| * |x| ^ (q + 1)) := by
    apply Summable.of_nonneg_of_le (fun q => by positivity) (fun q => ?_)
      ((Z03.summable_geom_mul (abs_nonneg x) hx).mul_left 3)
    have h1 : |x| ^ (q + 1) ≤ |x| ^ q := by
      rw [pow_succ]; exact mul_le_of_le_one_right (by positivity) hx.le
    calc |S05.lam c (q + 1)| * |x| ^ (q + 1) ≤ (3 * ((q : ℝ) + 1)) * |x| ^ q :=
          mul_le_mul (S05.lam_bound c q) h1 (by positivity) (by positivity)
      _ = _ := by ring
  have key := intervalIntegral.hasSum_integral_of_dominated_convergence
    (μ := MeasureTheory.volume) (a := 0) (b := 1)
    (F := fun (q : ℕ) (t : ℝ) => g t * (S05.lam c (q + 1) * (x * t) ^ (q + 1)))
    (f := fun t => g t * Real.log (1 - x * t) ^ c)
    (fun q t => |S05.lam c (q + 1)| * |x| ^ (q + 1) * |g t|)
    (fun q => (by fun_prop : Measurable (fun t : ℝ =>
      g t * (S05.lam c (q + 1) * (x * t) ^ (q + 1)))).aestronglyMeasurable)
    (fun q => by
      refine Filter.Eventually.of_forall (fun t ht => ?_)
      rw [Set.uIoc_of_le zero_le_one] at ht
      have hgt : ∀ G L : ℝ, ‖G * (L * (x * t) ^ (q + 1))‖ = |L| * |x| ^ (q + 1) * |G| * t ^ (q + 1) := by
        intro G L
        rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_pow, abs_mul, abs_of_pos ht.1]; ring
      rw [hgt]
      have : t ^ (q + 1) ≤ 1 := pow_le_one₀ ht.1.le ht.2
      have h0 : 0 ≤ |S05.lam c (q + 1)| * |x| ^ (q + 1) * |g t| := by positivity
      calc |S05.lam c (q + 1)| * |x| ^ (q + 1) * |g t| * t ^ (q + 1)
          ≤ |S05.lam c (q + 1)| * |x| ^ (q + 1) * |g t| * 1 := by gcongr
        _ = _ := by ring)
    (Filter.Eventually.of_forall (fun t _ => hS.mul_right _))
    (by
      simp_rw [tsum_mul_right]
      exact hgi.norm.const_mul _)
    (Filter.Eventually.of_forall (fun t ht => by
      rw [Set.uIoc_of_le zero_le_one] at ht
      have hxt : |x * t| < 1 := by
        rw [abs_mul, abs_of_pos ht.1]
        calc |x| * t ≤ |x| * 1 := by gcongr; exact ht.2
          _ < 1 := by linarith
      exact (S05.proof_2 c hc1 (by omega) (x * t) hxt).mul_left (g t)))
  have e1 : (∫ t in (0:ℝ)..1, g t * Real.log (1 - x * t) ^ c) = S05.Jint n a b c x := by
    unfold S05.Jint; rfl
  rw [e1] at key
  convert key using 1
  funext q
  unfold S05.Iint
  rw [← intervalIntegral.integral_const_mul]
  congr 1; funext t
  simp only [hg]
  rw [mul_pow, show n + q + 1 = n + (q + 1) by ring, pow_add]
  ring

namespace S05
open MeasureTheory Set Filter Real Interval Topology Finset

lemma hs_L1 {x : ℝ} (hx : |x| < 1) :
    HasSum (fun q : ℕ => x ^ (q + 1) / ((q : ℝ) + 1)) (cL1 x) :=
  Real.hasSum_pow_div_log_of_abs_lt_one hx

lemma hs_L1_shift1 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) :
    HasSum (fun q : ℕ => x ^ (q + 1) / ((q : ℝ) + 2)) ((cL1 x - x) / x) := by
  have h := (hasSum_nat_add_iff' 1).2 (hs_L1 hx)
  simp only [Finset.range_one, Finset.sum_singleton] at h
  have := h.div_const x
  convert this using 1
  · funext q; push_cast; rw [pow_succ x (q + 1)]; field_simp; ring
  · simp

lemma hs_L1_shift2 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) :
    HasSum (fun q : ℕ => x ^ (q + 1) / ((q : ℝ) + 3)) ((cL1 x - x - x ^ 2 / 2) / x ^ 2) := by
  have h := (hasSum_nat_add_iff' 2).2 (hs_L1 hx)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero] at h
  have := h.div_const (x ^ 2)
  convert this using 1
  · funext q; push_cast; rw [show q + 2 + 1 = (q + 1) + 2 by ring, pow_add]; field_simp; ring
  · norm_num; ring

lemma hs_Li2_shift1 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) :
    HasSum (fun q : ℕ => x ^ (q + 1) / ((q : ℝ) + 2) ^ 2) ((D00.Li2 x - x) / x) := by
  have h := (hasSum_nat_add_iff' 1).2 (D00.hasSum_Li2 hx.le)
  simp only [Finset.range_one, Finset.sum_singleton] at h
  have := h.div_const x
  convert this using 1
  · funext q; push_cast; rw [pow_succ x (q + 1)]; field_simp; ring
  · simp

lemma Iint00 (N : ℕ) : Iint N 0 0 = 1 / ((N : ℝ) + 1) := by
  simp [Iint, integral_pow]

lemma t_log_pow_tendsto (k : ℕ) (hk : k ≠ 0) :
    Tendsto (fun t : ℝ => t * log t ^ k) (𝓝[>] 0) (𝓝 0) := by
  have hr : (0:ℝ) < (k:ℝ)⁻¹ := by positivity
  have h := (tendsto_log_mul_rpow_nhdsGT_zero hr).pow k
  rw [zero_pow hk] at h
  refine h.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with t ht
  rw [mul_pow, rpow_inv_natCast_pow (le_of_lt ht) hk, mul_comm]

lemma G_cont (n k : ℕ) (hk : k ≠ 0) :
    ContinuousOn (fun t : ℝ => t ^ (n+1) * log t ^ k) (Icc 0 1) := by
  intro x hx
  rcases eq_or_lt_of_le hx.1 with h | h
  · subst h
    apply ContinuousWithinAt.mono _ Set.Icc_subset_Ici_self
    rw [← continuousWithinAt_Ioi_iff_Ici, ContinuousWithinAt]
    have := (((continuous_pow n).tendsto 0).mono_left nhdsWithin_le_nhds).mul
      (t_log_pow_tendsto k hk)
    rw [mul_zero] at this
    simp only [ne_eq, Nat.add_eq_zero_iff, one_ne_zero, and_false, not_false_eq_true, zero_pow,
      zero_mul]
    refine this.congr fun t => ?_
    ring
  · exact ((continuous_pow _).continuousAt.mul ((continuousAt_log h.ne').pow k)).continuousWithinAt


lemma Iint10 (N : ℕ) : Iint N 1 0 = -1 / ((N : ℝ) + 1) ^ 2 := by
  have hN : (0:ℝ) < (N : ℝ) + 1 := by positivity
  have hc : ContinuousOn (fun t : ℝ => t ^ (N + 1) * Real.log t / ((N : ℝ) + 1)
      - t ^ (N + 1) / ((N : ℝ) + 1) ^ 2) (Icc 0 1) := by
    have := G_cont N 1 one_ne_zero
    simp only [pow_one] at this
    exact (this.div_const _).sub ((continuous_pow _).continuousOn.div_const _)
  have hd : ∀ t ∈ Ioo (0:ℝ) 1, HasDerivAt (fun t : ℝ => t ^ (N + 1) * Real.log t / ((N : ℝ) + 1)
      - t ^ (N + 1) / ((N : ℝ) + 1) ^ 2) (t ^ N * Real.log t) t := by
    intro t ht
    have h1 := (((hasDerivAt_pow (N + 1) t).mul (Real.hasDerivAt_log ht.1.ne')).div_const
      ((N : ℝ) + 1)).sub ((hasDerivAt_pow (N + 1) t).div_const (((N : ℝ) + 1) ^ 2))
    convert h1 using 1
    have := ht.1.ne'
    push_cast
    field_simp
    ring
  have hi : IntervalIntegrable (fun t : ℝ => t ^ N * Real.log t) volume 0 1 := by
    have := ii_g N 1 0
    simpa using this
  have := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one hc hd hi
  unfold Iint
  simp only [pow_one, pow_zero, mul_one]
  rw [this]
  simp
  field_simp

end S05

theorem S05.proof_4_parts (x : ℝ) (hx : |x| < 1) (h0 : x ≠ 0) :
    S05.Jint 0 1 0 1 x = (-S05.cL1 x - D00.Li2 x + x * (S05.cL1 x + 2)) / x ∧
    S05.Jint 0 0 0 1 x = -S05.cL1 x + S05.cL1 x / x - 1 ∧
    S05.Jint 1 0 0 1 x = -S05.cL1 x / 2 + S05.cL1 x / (2 * x ^ 2) - 1 / 4 - 1 / (2 * x) := by
  refine ⟨?_, ?_, ?_⟩
  · have h := (S05.proof_3 x hx 0 1 0 1 le_rfl (by norm_num)).2
    have h2 := ((S05.hs_L1 hx).sub (S05.hs_L1_shift1 hx h0)).sub (S05.hs_Li2_shift1 hx h0)
    refine h.unique ?_
    convert h2 using 1
    · funext q
      simp only [S05.lam, S05.Iint10]
      push_cast
      field_simp
      ring
    · field_simp
      ring
  · have h := (S05.proof_3 x hx 0 0 0 1 le_rfl (by norm_num)).2
    have h2 := (S05.hs_L1_shift1 hx h0).sub (S05.hs_L1 hx)
    refine h.unique ?_
    convert h2 using 1
    · funext q
      simp only [S05.lam, S05.Iint00]
      push_cast
      field_simp
      ring
    · field_simp
      ring
  · have h := (S05.proof_3 x hx 1 0 0 1 le_rfl (by norm_num)).2
    have h2 := ((S05.hs_L1_shift2 hx h0).sub (S05.hs_L1 hx)).div_const 2
    refine h.unique ?_
    convert h2 using 1
    · funext q
      simp only [S05.lam, S05.Iint00]
      push_cast
      field_simp
      ring
    · field_simp
      ring

namespace S05
open MeasureTheory Set Filter Real Interval Topology Finset

/-- `p_k(u) = Σ_{j≤k} (1-u)^j/(k+1)`, so that `u p_k(u) = (1-(1-u)^{k+1})/(k+1)`. -/
noncomputable def pp (k : ℕ) (u : ℝ) : ℝ := (∑ j ∈ range (k + 1), (1 - u) ^ j) / ((k : ℝ) + 1)

/-- `q_k(u) = (1-(1-u)^{k+1})/(k+1)`, a primitive of `(1-u)^k` vanishing at `0`. -/
noncomputable def qq (k : ℕ) (u : ℝ) : ℝ := (1 - (1 - u) ^ (k + 1)) / ((k : ℝ) + 1)

/-- A primitive of `p_k` vanishing at `0`. -/
noncomputable def rq (k : ℕ) (u : ℝ) : ℝ :=
  (∑ j ∈ range (k + 1), (1 - (1 - u) ^ (j + 1)) / ((j : ℝ) + 1)) / ((k : ℝ) + 1)

/-- `F_k(u) = q_k(u) log u - r_k(u)`, a primitive of `(1-u)^k log u`. -/
noncomputable def FF (k : ℕ) (u : ℝ) : ℝ := qq k u * Real.log u - rq k u

lemma qq_eq (k : ℕ) (u : ℝ) : qq k u = pp k u * u := by
  unfold qq pp
  have := geom_sum_mul (1 - u) (k + 1)
  rw [div_mul_eq_mul_div]
  congr 1
  linear_combination this

lemma hasDerivAt_qq (k : ℕ) (u : ℝ) : HasDerivAt (qq k) ((1 - u) ^ k) u := by
  have h := ((((hasDerivAt_id' u).const_sub 1).fun_pow (k + 1)).const_sub 1).div_const
    ((k : ℝ) + 1)
  unfold qq
  refine h.congr_deriv ?_
  rw [Nat.add_sub_cancel]
  push_cast
  field_simp

lemma hasDerivAt_rq (k : ℕ) (u : ℝ) : HasDerivAt (rq k) (pp k u) u := by
  unfold rq pp
  apply HasDerivAt.div_const
  have : ∀ j ∈ range (k + 1), HasDerivAt (fun u => (1 - (1 - u) ^ (j + 1)) / ((j : ℝ) + 1))
      ((1 - u) ^ j) u := fun j _ => hasDerivAt_qq j u
  have h := HasDerivAt.sum this
  convert h using 1
  funext v
  simp [Finset.sum_apply]

lemma hasDerivAt_FF (k : ℕ) {u : ℝ} (hu : 0 < u) :
    HasDerivAt (FF k) ((1 - u) ^ k * Real.log u) u := by
  have h := ((hasDerivAt_qq k u).fun_mul (Real.hasDerivAt_log hu.ne')).fun_sub (hasDerivAt_rq k u)
  unfold FF
  refine h.congr_deriv ?_
  rw [qq_eq]
  field_simp
  ring

lemma FF_eq (k : ℕ) (u : ℝ) : FF k u = pp k u * (u * Real.log u) - rq k u := by
  unfold FF; rw [qq_eq]; ring

lemma continuous_pp (k : ℕ) : Continuous (pp k) := by unfold pp; fun_prop
lemma continuous_rq (k : ℕ) : Continuous (rq k) := by unfold rq; fun_prop

lemma continuousOn_FF (k : ℕ) : ContinuousOn (FF k) (Icc 0 1) := by
  have h := G_cont 0 1 one_ne_zero
  simp only [zero_add, pow_one] at h
  have : FF k = fun u => pp k u * (u * Real.log u) - rq k u := funext (FF_eq k)
  rw [this]
  exact ((continuous_pp k).continuousOn.mul h).sub (continuous_rq k).continuousOn

lemma FF_zero (k : ℕ) : FF k 0 = 0 := by simp [FF, qq, rq]

lemma FF_one (k : ℕ) : FF k 1 = -(harmonic (k + 1) : ℝ) / ((k : ℝ) + 1) := by
  simp only [FF, Real.log_one, mul_zero, zero_sub, rq, sub_self]
  rw [hR_sum]
  simp only [ne_eq, Nat.add_eq_zero_iff, one_ne_zero, and_false, not_false_eq_true, zero_pow,
    sub_zero]
  ring

lemma ii_pow_log (k m : ℕ) :
    IntervalIntegrable (fun u : ℝ => (1 - u) ^ k * Real.log u ^ m) volume 0 1 := by
  refine ii_of_ls (by fun_prop : Measurable (fun u : ℝ => (1 - u) ^ k * Real.log u ^ m)).aestronglyMeasurable
    1 m (fun t ht => ?_)
  rw [abs_mul, abs_pow, abs_pow, one_mul]
  have h1 : |1 - t| ^ k ≤ 1 :=
    pow_le_one₀ (abs_nonneg _) (by rw [abs_of_pos (by linarith [ht.2])]; linarith [ht.1])
  have a1 := pow_le_pow_left₀ (abs_nonneg _) (abs_log_le_ls t) m
  have := pow_nonneg (abs_nonneg (Real.log t)) m
  nlinarith

/-- `∫_0^1 (1-u)^k log u du = -H_{k+1}/(k+1)`. -/
lemma intP (k : ℕ) :
    ∫ u in (0:ℝ)..1, (1 - u) ^ k * Real.log u = -(harmonic (k + 1) : ℝ) / ((k : ℝ) + 1) := by
  have hi := ii_pow_log k 1
  simp only [pow_one] at hi
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one (continuousOn_FF k)
    (fun u hu => hasDerivAt_FF k hu.1) hi, FF_one, FF_zero, sub_zero]

/-- `G_N(u) = q_N(u) log²u - (2/(N+1)) Σ_{k≤N} F_k(u)`, a primitive of `(1-u)^N log²u`. -/
noncomputable def GG (N : ℕ) (u : ℝ) : ℝ :=
  qq N u * Real.log u ^ 2 - 2 / ((N : ℝ) + 1) * ∑ k ∈ range (N + 1), FF k u

lemma hasDerivAt_GG (N : ℕ) {u : ℝ} (hu : 0 < u) :
    HasDerivAt (GG N) ((1 - u) ^ N * Real.log u ^ 2) u := by
  have hs : HasDerivAt (fun v => ∑ k ∈ range (N + 1), FF k v)
      (∑ k ∈ range (N + 1), (1 - u) ^ k * Real.log u) u := by
    have := HasDerivAt.sum (u := range (N + 1)) (fun k _ => hasDerivAt_FF k hu)
    convert this using 1
    funext v; simp [Finset.sum_apply]
  have h := ((hasDerivAt_qq N u).fun_mul ((Real.hasDerivAt_log hu.ne').fun_pow 2)).fun_sub
    (hs.const_mul (2 / ((N : ℝ) + 1)))
  unfold GG
  refine h.congr_deriv ?_
  rw [qq_eq, ← Finset.sum_mul]
  unfold pp
  field_simp
  ring

lemma continuousOn_GG (N : ℕ) : ContinuousOn (GG N) (Icc 0 1) := by
  have h := G_cont 0 2 two_ne_zero
  simp only [zero_add, pow_one] at h
  have : GG N = fun u => pp N u * (u * Real.log u ^ 2)
      - 2 / ((N : ℝ) + 1) * ∑ k ∈ range (N + 1), FF k u := by
    funext u; unfold GG; rw [qq_eq]; ring
  rw [this]
  exact ((continuous_pp N).continuousOn.mul h).sub (continuousOn_const.mul
    (continuousOn_finsetSum (range (N + 1)) (fun k _ => continuousOn_FF k)))

lemma sum_harm_div (N : ℕ) :
    ∑ k ∈ range (N + 1), (harmonic (k + 1) : ℝ) / ((k : ℝ) + 1)
      = ((harmonic (N + 1) : ℝ) ^ 2 + H2 (N + 1)) / 2 := by
  induction N with
  | zero => simp [H2_succ, H2_zero]
  | succ N ih =>
    rw [sum_range_succ, ih, hR_succ (N + 1), H2_succ (N + 1)]
    push_cast
    field_simp
    ring

/-- `𝓘_N(0,2) = ∫_0^1 s^N log²(1-s) ds = (H_{N+1}² + H^{(2)}_{N+1})/(N+1)`. -/
lemma Iint02 (N : ℕ) :
    Iint N 0 2 = ((harmonic (N + 1) : ℝ) ^ 2 + H2 (N + 1)) / ((N : ℝ) + 1) := by
  have hsub : Iint N 0 2 = ∫ u in (0:ℝ)..1, (1 - u) ^ N * Real.log u ^ 2 := by
    unfold Iint
    have := intervalIntegral.integral_comp_sub_left
      (fun u : ℝ => (1 - u) ^ N * Real.log u ^ 2) (a := 0) (b := 1) 1
    simp only [sub_sub_cancel, sub_zero, sub_self] at this
    rw [← this]
    congr 1; funext s; ring
  rw [hsub, intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one (continuousOn_GG N)
    (fun u hu => hasDerivAt_GG N hu.1) (ii_pow_log N 2)]
  simp only [GG, Real.log_one, Real.log_zero]
  simp only [ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, mul_zero, zero_sub,
    FF_zero, Finset.sum_const_zero, sub_zero, FF_one]
  rw [show ∑ k ∈ range (N + 1), -(harmonic (k + 1) : ℝ) / ((k : ℝ) + 1)
      = -∑ k ∈ range (N + 1), (harmonic (k + 1) : ℝ) / ((k : ℝ) + 1) by
    rw [← Finset.sum_neg_distrib]; congr 1; funext k; ring]
  rw [sum_harm_div]
  field_simp

end S05

namespace S05
open MeasureTheory Set Filter Real Interval Topology Finset

lemma shift_div (g : ℕ → ℝ) {x F : ℝ} (h0 : x ≠ 0) (h : HasSum (fun q => g q * x ^ (q + 1)) F) :
    HasSum (fun q => g (q + 1) * x ^ (q + 1)) ((F - g 0 * x) / x) := by
  have h1 := (hasSum_nat_add_iff' 1).2 h
  simp only [Finset.range_one, Finset.sum_singleton, zero_add, pow_one] at h1
  have := h1.div_const x
  convert this using 1
  funext q; rw [pow_succ x (q + 1)]; field_simp

/-- The worked entry `J(0;0,2,1)`. -/
lemma J0021 (x : ℝ) (hx : |x| < 1) (h0 : x ≠ 0) :
    Jint 0 0 2 1 x =
      (1 - x) / (3 * x) * (12 * Z03.Li3 x - 6 * Z03.A x + 6 * cL1 x * D00.Li2 x
        + 6 * D00.Li2 x + cL1 x ^ 3 + 3 * cL1 x ^ 2 + 6 * cL1 x) - 6 := by
  have h := (S05.proof_3 x hx 0 0 2 1 le_rfl (by norm_num)).2
  have hE1 := euler4 hx
  have hE2 := euler5 hx
  have hE3 := euler6 hx
  have sE1 := shift_div (fun q => (harmonic (q + 1) : ℝ) / ((q : ℝ) + 1)) h0 hE1
  have sE2 := shift_div (fun q => (harmonic (q + 1) : ℝ) ^ 2 / ((q : ℝ) + 1)) h0 hE2
  have sE3 := shift_div (fun q => H2 (q + 1) / ((q : ℝ) + 1)) h0 hE3
  have hL := hs_L1 hx
  have hLs := hs_L1_shift1 hx h0
  have hc := (((((((hE2.add hE3).neg.sub (hE1.mul_left 2)).sub (hL.mul_left 2)).add
    (hLs.mul_left 2)).add (sE1.mul_left 2)).add sE2).add sE3)
  refine h.unique ?_
  convert hc using 1
  · funext q
    simp only [lam, Iint02, zero_add]
    push_cast
    rw [show (q : ℝ) + 1 + 1 = (q : ℝ) + 2 by ring]
    rw [show q + 1 + 1 = (q + 1) + 1 by ring, hR_succ (q + 1), H2_succ (q + 1)]
    push_cast
    field_simp
    ring
  · have e1 : (harmonic (0 + 1) : ℝ) = 1 := by simp
    have e2 : H2 (0 + 1) = 1 := by simp [H2]
    rw [e1, e2]
    field_simp
    ring

end S05

theorem S05.proof_4 : S05.stmt_4 := by
  intro x hx h0
  obtain ⟨e1, e2, e3⟩ := S05.proof_4_parts x hx h0
  exact ⟨e1, e2, e3, S05.J0021 x hx h0⟩

namespace S05
open MeasureTheory Set Filter Real Topology

lemma eVar_add_le (f g : ℝ → ℝ) (s : Set ℝ) :
    eVariationOn (fun x => f x + g x) s ≤ eVariationOn f s + eVariationOn g s := by
  apply iSup_le
  rintro ⟨n, ⟨u, hu, us⟩⟩
  calc ∑ i ∈ Finset.range n, edist (f (u (i + 1)) + g (u (i + 1))) (f (u i) + g (u i))
      ≤ ∑ i ∈ Finset.range n, (edist (f (u (i + 1))) (f (u i)) + edist (g (u (i + 1))) (g (u i))) := by
        gcongr with i
        exact edist_add_add_le _ _ _ _
    _ = ∑ i ∈ Finset.range n, edist (f (u (i + 1))) (f (u i))
        + ∑ i ∈ Finset.range n, edist (g (u (i + 1))) (g (u i)) := Finset.sum_add_distrib
    _ ≤ eVariationOn f s + eVariationOn g s := by
        gcongr
        · exact eVariationOn.sum_le (f := f) (n := n) hu us
        · exact eVariationOn.sum_le (f := g) (n := n) hu us

lemma bv_add {f g : ℝ → ℝ} {s : Set ℝ} (hf : BoundedVariationOn f s)
    (hg : BoundedVariationOn g s) : BoundedVariationOn (fun x => f x + g x) s :=
  ne_top_of_le_ne_top (ENNReal.add_ne_top.2 ⟨hf, hg⟩) (eVar_add_le f g s)

lemma bv_const_mul {f : ℝ → ℝ} {s : Set ℝ} (c : ℝ) (hf : BoundedVariationOn f s) :
    BoundedVariationOn (fun x => c * f x) s := by
  have hc : BoundedVariationOn (fun _ : ℝ => c) s := by
    unfold BoundedVariationOn
    rw [eVariationOn.constant_on (by simp [Set.Subsingleton])]
    exact ENNReal.zero_ne_top
  exact hc.fun_mul hf

/-- `|cos|^p` has bounded variation on `[0, π]`. -/
lemma bv_abs_cos_pi {p : ℝ} (hp : 0 < p) :
    BoundedVariationOn (fun s => |Real.cos s| ^ p) (Icc 0 π) := by
  have hpi := Real.pi_pos
  have split := eVariationOn.Icc_add_Icc (fun s => |Real.cos s| ^ p) (s := univ)
    (a := 0) (b := π / 2) (c := π) (by linarith) (by linarith) (mem_univ _)
  simp only [univ_inter] at split
  unfold BoundedVariationOn
  rw [← split]
  refine ENNReal.add_ne_top.2 ⟨?_, ?_⟩
  · -- on `[0, π/2]`: reflect `s ↦ π/2 - s`, getting `sin^p`, monotone
    have hφ : AntitoneOn (fun s : ℝ => π / 2 - s) (Icc 0 (π / 2)) := fun a _ b _ h => by
      simp only; linarith
    have himg : (fun s : ℝ => π / 2 - s) '' Icc 0 (π / 2) = Icc 0 (π / 2) := by
      rw [Set.image_const_sub_Icc]; simp
    have := eVariationOn.comp_eq_of_antitoneOn (fun s => |Real.cos s| ^ p) _ hφ
      (t := Icc 0 (π / 2))
    rw [himg] at this
    rw [← this]
    have hmono : MonotoneOn ((fun s => |Real.cos s| ^ p) ∘ fun s : ℝ => π / 2 - s)
        (Icc 0 (π / 2)) := by
      intro a ha b hb hab
      simp only [Function.comp, Real.cos_pi_div_two_sub]
      have ha0 : 0 ≤ Real.sin a := Real.sin_nonneg_of_nonneg_of_le_pi ha.1 (by linarith [ha.2])
      have hb0 : 0 ≤ Real.sin b := Real.sin_nonneg_of_nonneg_of_le_pi hb.1 (by linarith [hb.2])
      rw [abs_of_nonneg ha0, abs_of_nonneg hb0]
      apply Real.rpow_le_rpow ha0 _ hp.le
      exact Real.sin_le_sin_of_le_of_le_pi_div_two (by linarith [ha.1]) hb.2 hab
    have h2 := hmono.eVariationOn_eq (a := 0) (b := π / 2) ⟨le_rfl, by linarith⟩
      ⟨by linarith, le_rfl⟩
    simp only [Set.inter_self] at h2
    rw [h2]; exact ENNReal.ofReal_ne_top
  · have hmono : MonotoneOn (fun s => |Real.cos s| ^ p) (Icc (π / 2) π) := by
      intro a ha b hb hab
      have ha0 : Real.cos a ≤ 0 := Real.cos_nonpos_of_pi_div_two_le_of_le ha.1 (by linarith [ha.2])
      have hb0 : Real.cos b ≤ 0 := Real.cos_nonpos_of_pi_div_two_le_of_le hb.1 (by linarith [hb.2])
      simp only
      rw [abs_of_nonpos ha0, abs_of_nonpos hb0]
      apply Real.rpow_le_rpow (by linarith) _ hp.le
      have := Real.cos_le_cos_of_nonneg_of_le_pi (by linarith [ha.1]) hb.2 hab
      linarith
    have := hmono.eVariationOn_eq (a := π / 2) (b := π) ⟨le_rfl, by linarith⟩ ⟨by linarith, le_rfl⟩
    simp only [Set.inter_self] at this
    rw [this]; exact ENNReal.ofReal_ne_top

end S05

namespace S05
open MeasureTheory Set Filter Real Topology

lemma bv_periodic (h : ℝ → ℝ) {T : ℝ} (hT : 0 < T) (hper : Function.Periodic h T)
    (h0 : BoundedVariationOn h (Icc 0 T)) (a b : ℝ) : BoundedVariationOn h (Icc a b) := by
  -- one period at `kT`
  have one : ∀ k : ℤ, BoundedVariationOn h (Icc (k * T) ((k + 1) * T)) := by
    intro k
    have hφ : MonotoneOn (fun s : ℝ => s + k * T) (Icc 0 T) := fun x _ y _ hxy => by
      simp only; linarith
    have himg : (fun s : ℝ => s + k * T) '' Icc 0 T = Icc (k * T) ((k + 1) * T) := by
      rw [Set.image_add_const_Icc]; congr 1 <;> ring
    have := eVariationOn.comp_eq_of_monotoneOn h _ hφ (t := Icc 0 T)
    rw [himg] at this
    unfold BoundedVariationOn
    rw [← this]
    have hc : h ∘ (fun s : ℝ => s + k * T) = h := by
      funext s; simp only [Function.comp]; exact (hper.int_mul k) s
    rw [hc]; exact h0
  -- `n` periods from `kT`
  have many : ∀ (k : ℤ) (n : ℕ), BoundedVariationOn h (Icc (k * T) ((k + n) * T)) := by
    intro k n
    induction n with
    | zero =>
      unfold BoundedVariationOn
      rw [eVariationOn.subsingleton]
      · exact ENNReal.zero_ne_top
      · simp only [Nat.cast_zero, add_zero, Icc_self]; exact subsingleton_singleton
    | succ n ih =>
      have := eVariationOn.Icc_add_Icc h (s := univ) (a := k * T) (b := (k + n) * T)
        (c := (k + (n + 1 : ℕ)) * T) (by nlinarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)])
        (by push_cast; nlinarith) (mem_univ _)
      simp only [univ_inter] at this
      unfold BoundedVariationOn
      rw [← this]
      refine ENNReal.add_ne_top.2 ⟨ih, ?_⟩
      have h1 := one (k + n)
      push_cast at h1 ⊢
      convert h1 using 2; ring
  -- cover `[a, b]`
  rcases le_or_gt a b with hab | hab
  · set k : ℤ := ⌊a / T⌋
    obtain ⟨n, hn⟩ := exists_nat_gt ((b - k * T) / T)
    have hk : (k : ℝ) * T ≤ a := by
      have := Int.floor_le (a / T)
      rwa [le_div_iff₀ hT] at this
    have hn' : b ≤ (k + n) * T := by
      rw [div_lt_iff₀ hT] at hn; nlinarith
    exact (many k n).mono (Icc_subset_Icc hk hn')
  · unfold BoundedVariationOn
    rw [eVariationOn.subsingleton]
    · exact ENNReal.zero_ne_top
    · rw [Icc_eq_empty (not_le.2 hab)]; exact subsingleton_empty

lemma bv_abs_cos (p : ℝ) (hp : 0 < p) (a b : ℝ) :
    BoundedVariationOn (fun s => |Real.cos s| ^ p) (Icc a b) :=
  bv_periodic _ Real.pi_pos (fun s => by simp [Real.cos_add_pi, abs_neg]) (bv_abs_cos_pi hp) a b

end S05

namespace S05
open MeasureTheory Set Filter Real Topology

lemma dot_ray (u : X04.Plane) (θ : ℝ) : X04.dot u (ray θ) = u 0 * Real.cos θ + u 1 * Real.sin θ := by
  simp [X04.dot, ray]

lemma bv_dot_ray (u : X04.Plane) {p : ℝ} (hp : 0 < p) (a b : ℝ) :
    BoundedVariationOn (fun θ => |X04.dot u (ray θ)| ^ p) (Icc a b) := by
  set z : ℂ := ⟨u 0, u 1⟩
  have h0 : ‖z‖ * Real.cos (Complex.arg z) = u 0 := Complex.norm_mul_cos_arg z
  have h1 : ‖z‖ * Real.sin (Complex.arg z) = u 1 := Complex.norm_mul_sin_arg z
  set α := Complex.arg z
  have e : (fun θ => |X04.dot u (ray θ)| ^ p)
      = fun θ => ‖z‖ ^ p * ((fun s => |Real.cos s| ^ p) ∘ fun s => s - α) θ := by
    funext θ
    simp only [Function.comp, dot_ray, ← h0, ← h1]
    rw [show ‖z‖ * Real.cos α * Real.cos θ + ‖z‖ * Real.sin α * Real.sin θ
      = ‖z‖ * Real.cos (θ - α) by rw [Real.cos_sub]; ring]
    rw [abs_mul, abs_of_nonneg (norm_nonneg _), Real.mul_rpow (norm_nonneg _) (abs_nonneg _)]
  rw [e]
  apply bv_const_mul
  have hφ : MonotoneOn (fun s : ℝ => s - α) (Icc a b) := fun x _ y _ hxy => by
    simp only; linarith
  unfold BoundedVariationOn
  rw [eVariationOn.comp_eq_of_monotoneOn _ _ hφ, Set.image_sub_const_Icc]
  exact bv_abs_cos p hp _ _

end S05

namespace S05
open MeasureTheory Set Filter Real Topology

lemma dot_cont' (u : X04.Plane) : Continuous (fun x : X04.Plane => X04.dot u x) := by
  unfold X04.dot; fun_prop

lemma Fw_cont (u : Fin 3 → X04.Plane) (w : Fin 3 → ℝ) {p : ℝ} (hp : 0 < p) :
    Continuous (Fw u w p) := by
  unfold Fw
  exact continuous_finsetSum _ (fun j _ => continuous_const.mul
    ((dot_cont' _).abs.rpow_const (fun _ => Or.inr hp.le)))

lemma Fw_smul (u : Fin 3 → X04.Plane) (w : Fin 3 → ℝ) (p c : ℝ) (x : X04.Plane) :
    Fw u w p (c • x) = |c| ^ p * Fw u w p x := by
  unfold Fw
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  have : X04.dot (u j) (c • x) = c * X04.dot (u j) x := by
    simp only [X04.dot, PiLp.smul_apply, smul_eq_mul]; ring
  rw [this, abs_mul, Real.mul_rpow (abs_nonneg _) (abs_nonneg _)]; ring

lemma Fw_pos (H : Hex) (hH : IsHex H) {p : ℝ} {x : X04.Plane} (hx : x ≠ 0) :
    0 < Fw H.u H.w p x := by
  obtain ⟨hw, hdet, -⟩ := hH
  have hne : ∃ j, X04.dot (H.u j) x ≠ 0 := by
    by_contra hc
    push Not at hc
    have h0 := hc 0
    have h1 := hc 1
    have hd := hdet 0 1 (by decide)
    unfold X04.dot at h0 h1
    unfold det2 at hd
    have e0 : (H.u 0 0 * H.u 1 1 - H.u 0 1 * H.u 1 0) * x 0 = 0 := by
      linear_combination (H.u 1 1) * h0 - (H.u 0 1) * h1
    have e1 : (H.u 0 0 * H.u 1 1 - H.u 0 1 * H.u 1 0) * x 1 = 0 := by
      linear_combination (H.u 0 0) * h1 - (H.u 1 0) * h0
    apply hx
    ext i; fin_cases i
    · simpa [hd] using e0
    · simpa [hd] using e1
  obtain ⟨j, hj⟩ := hne
  unfold Fw
  exact Finset.sum_pos' (fun i _ => mul_nonneg (hw i).le (Real.rpow_nonneg (abs_nonneg _) _))
    ⟨j, Finset.mem_univ _, mul_pos (hw j) (Real.rpow_pos_of_pos (abs_pos.2 hj) _)⟩

lemma Fw_zero (u : Fin 3 → X04.Plane) (w : Fin 3 → ℝ) {p : ℝ} (hp : 0 < p) :
    Fw u w p 0 = 0 := by
  simp [Fw, X04.dot, Real.zero_rpow hp.ne']

lemma ray_ne (θ : ℝ) : ray θ ≠ 0 := by
  intro h
  have h0 : Real.cos θ = 0 := by
    have := congrArg (fun v : X04.Plane => v 0) h; simpa [ray] using this
  have h1 : Real.sin θ = 0 := by
    have := congrArg (fun v : X04.Plane => v 1) h; simpa [ray] using this
  have := Real.sin_sq_add_cos_sq θ
  rw [h0, h1] at this; norm_num at this

lemma Phi_pos (H : Hex) (hH : IsHex H) {p : ℝ} (θ : ℝ) : 0 < Fw H.u H.w p (ray θ) :=
  Fw_pos H hH (ray_ne θ)

lemma continuous_ray : Continuous ray := by
  unfold ray
  fun_prop

lemma polar_cont (H : Hex) (hH : IsHex H) {p : ℝ} (hp : 0 < p) : Continuous (polar H p) := by
  have h1 : Continuous (fun θ => (Fw H.u H.w p (ray θ)) ^ (-(1 / p))) :=
    ((Fw_cont H.u H.w hp).comp continuous_ray).rpow_const
      (fun θ => Or.inl (Phi_pos H hH θ).ne')
  exact h1.smul continuous_ray

lemma F_polar (H : Hex) (hH : IsHex H) {p : ℝ} (hp : 0 < p) (θ : ℝ) :
    Fw H.u H.w p (polar H p θ) = 1 := by
  have hF := Phi_pos H hH (p := p) θ
  unfold polar
  rw [Fw_smul, abs_of_pos (Real.rpow_pos_of_pos hF _), ← Real.rpow_mul hF.le,
    show -(1 / p) * p = -1 by field_simp, Real.rpow_neg_one, inv_mul_cancel₀ hF.ne']

lemma polar_coords (x : X04.Plane) : ∃ r θ : ℝ, 0 ≤ r ∧ 0 ≤ θ ∧ θ < 2 * π ∧ x = r • ray θ := by
  set z : ℂ := ⟨x 0, x 1⟩
  have h0 : ‖z‖ * Real.cos (Complex.arg z) = x 0 := Complex.norm_mul_cos_arg z
  have h1 : ‖z‖ * Real.sin (Complex.arg z) = x 1 := Complex.norm_mul_sin_arg z
  have ha := Complex.neg_pi_lt_arg z
  have hb := Complex.arg_le_pi z
  have hpi := Real.pi_pos
  by_cases hc : Complex.arg z < 0
  · refine ⟨‖z‖, Complex.arg z + 2 * π, norm_nonneg _, by linarith, by linarith, ?_⟩
    ext i; fin_cases i
    · simp [ray, Real.cos_add_two_pi, h0]
    · simp [ray, Real.sin_add_two_pi, h1]
  · refine ⟨‖z‖, Complex.arg z, norm_nonneg _, by linarith, by linarith, ?_⟩
    ext i; fin_cases i
    · simp [ray, h0]
    · simp [ray, h1]

end S05

namespace S05
open MeasureTheory Set Filter Real Topology

lemma polar_inj (H : Hex) (hH : IsHex H) {p : ℝ} :
    InjOn (polar H p) (Ico 0 (2 * π)) := by
  intro a ha b hb hab
  have hpa := Real.rpow_pos_of_pos (Phi_pos H hH (p := p) a) (-(1 / p))
  have hpb := Real.rpow_pos_of_pos (Phi_pos H hH (p := p) b) (-(1 / p))
  set ra := (Fw H.u H.w p (ray a)) ^ (-(1 / p))
  set rb := (Fw H.u H.w p (ray b)) ^ (-(1 / p))
  have r0 : ∀ θ, ray θ 0 = Real.cos θ := fun θ => by simp [ray]
  have r1 : ∀ θ, ray θ 1 = Real.sin θ := fun θ => by simp [ray]
  have e0 : ra * Real.cos a = rb * Real.cos b := by
    have := congrArg (fun v : X04.Plane => v 0) hab
    simp only [polar, PiLp.smul_apply, smul_eq_mul, r0] at this; exact this
  have e1 : ra * Real.sin a = rb * Real.sin b := by
    have := congrArg (fun v : X04.Plane => v 1) hab
    simp only [polar, PiLp.smul_apply, smul_eq_mul, r1] at this; exact this
  have sa := Real.sin_sq_add_cos_sq a
  have sb := Real.sin_sq_add_cos_sq b
  have hr : ra = rb := by
    have : ra ^ 2 = rb ^ 2 := by
      have : ra ^ 2 * (Real.sin a ^ 2 + Real.cos a ^ 2) = rb ^ 2 * (Real.sin b ^ 2 + Real.cos b ^ 2) := by
        linear_combination (ra * Real.sin a + rb * Real.sin b) * e1 + (ra * Real.cos a + rb * Real.cos b) * e0
      rw [sa, sb] at this; linarith
    nlinarith
  rw [hr] at e0 e1
  have hc : Real.cos a = Real.cos b := mul_left_cancel₀ hpb.ne' e0
  have hs : Real.sin a = Real.sin b := mul_left_cancel₀ hpb.ne' e1
  obtain ⟨k, hk⟩ := Real.Angle.angle_eq_iff_two_pi_dvd_sub.1 (Real.Angle.cos_sin_inj hc hs)
  have hpi := Real.pi_pos
  have h1 : (k : ℝ) < 1 := by
    by_contra h; push Not at h; nlinarith [ha.1, ha.2, hb.1, hb.2]
  have h2 : (-1 : ℝ) < k := by
    by_contra h; push Not at h; nlinarith [ha.1, ha.2, hb.1, hb.2]
  have hk1 : k < 1 := by exact_mod_cast h1
  have hk2 : -1 < k := by exact_mod_cast h2
  have hk0 : k = 0 := by omega
  rw [hk0] at hk; simp at hk; linarith

lemma polar_image (H : Hex) (hH : IsHex H) {p : ℝ} (hp : 0 < p) :
    polar H p '' Icc 0 (2 * π) = {x | Fw H.u H.w p x = 1} := by
  apply Subset.antisymm
  · rintro _ ⟨θ, -, rfl⟩; exact F_polar H hH hp θ
  · intro x hx
    simp only [Set.mem_ofPred_eq] at hx
    obtain ⟨r, θ, hr, hθ0, hθ1, rfl⟩ := polar_coords x
    refine ⟨θ, ⟨hθ0, hθ1.le⟩, ?_⟩
    have hΦ := Phi_pos H hH (p := p) θ
    rw [Fw_smul, abs_of_nonneg hr] at hx
    have hrp : r ^ p = (Fw H.u H.w p (ray θ))⁻¹ := by
      field_simp; linarith
    have : r = (Fw H.u H.w p (ray θ)) ^ (-(1 / p)) := by
      rw [← Real.rpow_rpow_inv hr hp.ne', hrp, Real.inv_rpow hΦ.le, Real.rpow_neg hΦ.le,
        one_div]
    unfold polar; rw [this]

lemma frontier_Kw (H : Hex) (hH : IsHex H) {p : ℝ} (hp : 0 < p) :
    frontier (Kw H.u H.w p) = {x | Fw H.u H.w p x = 1} := by
  apply subset_antisymm
  · exact frontier_le_subset_eq (Fw_cont H.u H.w hp) continuous_const
  · intro x hx
    simp only [Set.mem_ofPred_eq] at hx
    unfold Kw
    rw [frontier, (isClosed_le (Fw_cont H.u H.w hp) continuous_const).closure_eq]
    refine ⟨le_of_eq hx, fun hint => ?_⟩
    have hK : {y | Fw H.u H.w p y ≤ 1} ∈ 𝓝 x := mem_interior_iff_mem_nhds.1 hint
    have ht : Tendsto (fun c : ℝ => c • x) (𝓝 1) (𝓝 x) := by
      have := (tendsto_id (x := 𝓝 (1:ℝ))).smul_const x
      simpa using this
    have hev := (ht.eventually hK).filter_mono (nhdsWithin_le_nhds (s := Ioi (1:ℝ)))
    obtain ⟨c, hc1, hc2⟩ := (hev.and self_mem_nhdsWithin).exists
    simp only [Fw_smul, hx, mul_one] at hc1
    have hc' : 1 < c := hc2
    have : 1 < |c| ^ p := Real.one_lt_rpow (by rw [abs_of_pos (by linarith)]; exact hc') hp
    linarith

lemma lipschitz_ray : LipschitzWith 2 ray := by
  apply LipschitzWith.of_dist_le_mul
  intro a b
  rw [EuclideanSpace.dist_eq, Fin.sum_univ_two]
  simp only [ray, Real.dist_eq]
  simp
  have hc := Real.abs_cos_sub_cos_le a b
  have hs := Real.abs_sin_sub_sin_le a b
  have hc2 : (Real.cos a - Real.cos b) ^ 2 ≤ (a - b) ^ 2 := by
    rw [← sq_abs, ← sq_abs (a - b)]; exact pow_le_pow_left₀ (abs_nonneg _) hc 2
  have hs2 : (Real.sin a - Real.sin b) ^ 2 ≤ (a - b) ^ 2 := by
    rw [← sq_abs, ← sq_abs (a - b)]; exact pow_le_pow_left₀ (abs_nonneg _) hs 2
  rw [Real.sqrt_le_iff]
  refine ⟨by positivity, ?_⟩
  rw [mul_pow, sq_abs]
  nlinarith [sq_nonneg (a - b)]

lemma bv_Phi (H : Hex) {p : ℝ} (hp : 0 < p) :
    BoundedVariationOn (fun θ => Fw H.u H.w p (ray θ)) (Icc 0 (2 * π)) := by
  have e : (fun θ => Fw H.u H.w p (ray θ)) = fun θ =>
      (H.w 0 * |X04.dot (H.u 0) (ray θ)| ^ p + H.w 1 * |X04.dot (H.u 1) (ray θ)| ^ p)
        + H.w 2 * |X04.dot (H.u 2) (ray θ)| ^ p := by
    funext θ; simp [Fw, Fin.sum_univ_three]
  rw [e]
  exact bv_add (bv_add (bv_const_mul _ (bv_dot_ray _ hp _ _)) (bv_const_mul _ (bv_dot_ray _ hp _ _)))
    (bv_const_mul _ (bv_dot_ray _ hp _ _))

lemma bv_polar (H : Hex) (hH : IsHex H) {p : ℝ} (hp : 0 < p) :
    BoundedVariationOn (polar H p) (Icc 0 (2 * π)) := by
  have hpi := Real.pi_pos
  set Φ := fun θ => Fw H.u H.w p (ray θ) with hΦdef
  have hΦc : Continuous Φ := (Fw_cont H.u H.w hp).comp continuous_ray
  obtain ⟨θm, hθm, hmin⟩ := isCompact_Icc.exists_isMinOn
    (nonempty_Icc.2 (by linarith : (0:ℝ) ≤ 2 * π)) hΦc.continuousOn
  obtain ⟨θM, hθM, hmax⟩ := isCompact_Icc.exists_isMaxOn
    (nonempty_Icc.2 (by linarith : (0:ℝ) ≤ 2 * π)) hΦc.continuousOn
  set m := Φ θm
  set M := Φ θM
  have hm : 0 < m := Phi_pos H hH θm
  have hmaps : MapsTo Φ (Icc 0 (2 * π)) (Icc m M) := fun θ hθ => ⟨hmin hθ, hmax hθ⟩
  set K := (1 / p) * m ^ (-(1 / p) - 1)
  have hK : 0 ≤ K := by positivity
  have hlip : LipschitzOnWith K.toNNReal (fun y : ℝ => y ^ (-(1 / p))) (Icc m M) := by
    refine Convex.lipschitzOnWith_of_nnNorm_hasDerivWithin_le (f' := fun y => -(1 / p) * y ^ (-(1 / p) - 1))
      (convex_Icc m M) (fun y hy => ?_) (fun y hy => ?_)
    · have hy0 : 0 < y := lt_of_lt_of_le hm hy.1
      exact (Real.hasDerivAt_rpow_const (Or.inl hy0.ne')).hasDerivWithinAt
    · have hy0 : 0 < y := lt_of_lt_of_le hm hy.1
      rw [← NNReal.coe_le_coe, coe_nnnorm, Real.coe_toNNReal _ hK, Real.norm_eq_abs, abs_mul, abs_neg,
        abs_of_pos (by positivity : (0:ℝ) < 1 / p), abs_of_pos (Real.rpow_pos_of_pos hy0 _)]
      have : y ^ (-(1 / p) - 1) ≤ m ^ (-(1 / p) - 1) :=
        Real.rpow_le_rpow_of_nonpos hm hy.1 (by have h1 : 0 < 1 / p := (by positivity); linarith)
      exact mul_le_mul_of_nonneg_left this (by positivity)
  have hr : BoundedVariationOn (fun θ => (Φ θ) ^ (-(1 / p))) (Icc 0 (2 * π)) :=
    hlip.comp_boundedVariationOn hmaps (bv_Phi H hp)
  have hray : BoundedVariationOn ray (Icc 0 (2 * π)) := by
    have := (lipschitz_ray.lipschitzOnWith (s := univ)).comp_boundedVariationOn (g := id)
      (mapsTo_univ _ _) (BoundedVariationOn.id_Icc 0 (2 * π))
    simpa using this
  exact hr.fun_smul hray

end S05

theorem S05.proof_7 : S05.stmt_7 := by
  intro H hH p hp
  refine ⟨(S05.polar_cont H hH hp).continuousOn, ?_, S05.polar_inj H hH, S05.polar_image H hH hp,
    S05.frontier_Kw H hH hp, S05.bv_polar H hH hp⟩
  simp [S05.polar, S05.ray]

namespace S05
open Set Real Filter Topology


/-- `|∂^n f(p)|` within `Ioi 0`. -/
noncomputable def nd (n : ℕ) (f : ℝ → ℝ) (p : ℝ) : ℝ := |iteratedDerivWithin n f (Ioi 0) p|

lemma uS : UniqueDiffOn ℝ (Ioi (0:ℝ)) := isOpen_Ioi.uniqueDiffOn

lemma nd_F (n : ℕ) (f : ℝ → ℝ) (p : ℝ) :
    ‖iteratedFDerivWithin ℝ n f (Ioi 0) p‖ = nd n f p := by
  rw [norm_iteratedFDerivWithin_eq_norm_iteratedDerivWithin, Real.norm_eq_abs]; rfl

lemma nd_nonneg (n : ℕ) (f : ℝ → ℝ) (p : ℝ) : 0 ≤ nd n f p := abs_nonneg _

lemma nd_zero (f : ℝ → ℝ) (p : ℝ) : nd 0 f p = |f p| := by simp [nd]

lemma nd_exp {h : ℝ → ℝ} (hh : ContDiffOn ℝ 3 h (Ioi 0)) {p : ℝ} (hp : 0 < p) {n : ℕ}
    (hn : n ≤ 3) {D : ℝ} (hD : ∀ i, 1 ≤ i → i ≤ n → nd i h p ≤ D ^ i) :
    nd n (fun q => Real.exp (h q)) p ≤ (n.factorial : ℝ) * Real.exp (h p) * D ^ n := by
  have := norm_iteratedFDerivWithin_comp_le (g := Real.exp) (f := h) (t := univ) (s := Ioi 0)
    (x := p) (n := n) (N := 3) (contDiff_exp.contDiffOn) hh (by exact_mod_cast hn)
    uniqueDiffOn_univ uS (mapsTo_univ _ _) hp (C := Real.exp (h p)) (D := D) ?_
    (fun i h1 h2 => by rw [nd_F]; exact hD i h1 h2)
  · rw [← nd_F]; exact this
  · intro i _
    rw [iteratedFDerivWithin_univ, norm_iteratedFDeriv_eq_norm_iteratedDeriv,
      iteratedDeriv_eq_iterate, Real.iter_deriv_exp, Real.norm_eq_abs, abs_of_pos (exp_pos _)]

lemma log_derivs (i : ℕ) (hi : i ≤ 3) : |iteratedDeriv i Real.log 1| ≤ 2 := by
  rcases i with _ | i
  · simp
  · rw [iteratedDeriv_succ', Real.deriv_log', iteratedDeriv_eq_iterate, iter_deriv_inv]
    simp only [one_zpow, mul_one, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul]
    have : i ≤ 2 := by omega
    interval_cases i <;> simp [Nat.factorial]

lemma nd_log {f : ℝ → ℝ} (hf : ContDiffOn ℝ 3 f (Ioi 0)) (hpos : ∀ q ∈ Ioi (0:ℝ), 0 < f q)
    {p : ℝ} (hp : 0 < p) {n : ℕ} (hn1 : 1 ≤ n) (hn : n ≤ 3) {D : ℝ}
    (hD : ∀ i, 1 ≤ i → i ≤ n → nd i f p ≤ f p * D ^ i) :
    nd n (fun q => Real.log (f q)) p ≤ (n.factorial : ℝ) * 2 * D ^ n := by
  have hc : 0 < f p := hpos p hp
  set φ : ℝ → ℝ := fun q => (f p)⁻¹ * f q with hφ
  have hEq : EqOn (fun q => Real.log (f q)) (fun q => Real.log (f p) + Real.log (φ q)) (Ioi 0) := by
    intro q hq
    simp only [hφ]
    rw [Real.log_mul (inv_ne_zero hc.ne') (hpos q hq).ne', Real.log_inv]; ring
  have h1 : nd n (fun q => Real.log (f q)) p = nd n (fun q => Real.log (φ q)) p := by
    unfold nd
    rw [iteratedDerivWithin_congr hEq hp, iteratedDerivWithin_const_add (by omega)]
  rw [h1]
  have hφc : ContDiffOn ℝ 3 φ (Ioi 0) := contDiffOn_const.mul hf
  have hφD : ∀ i, 1 ≤ i → i ≤ n → nd i φ p ≤ D ^ i := by
    intro i hi1 hi2
    have : nd i φ p = (f p)⁻¹ * nd i f p := by
      unfold nd
      rw [hφ, iteratedDerivWithin_const_mul_field, abs_mul, abs_of_pos (inv_pos.2 hc)]
    rw [this, inv_mul_le_iff₀ hc]; exact hD i hi1 hi2
  have := norm_iteratedFDerivWithin_comp_le (g := Real.log) (f := φ) (t := Ioi 0) (s := Ioi 0)
    (x := p) (n := n) (N := 3) (Real.contDiffOn_log.mono (fun x hx => ne_of_gt hx)) hφc
    (by exact_mod_cast hn) uS uS (fun q hq => by
      simp only [hφ, mem_Ioi]; exact mul_pos (inv_pos.2 hc) (hpos q hq)) hp (C := 2) (D := D) ?_
    (fun i h1 h2 => by rw [nd_F]; exact hφD i h1 h2)
  · rw [← nd_F]; simpa [Function.comp_def, mul_assoc] using this
  · intro i hi
    have h1 : φ p = 1 := by simp [hφ, inv_mul_cancel₀ hc.ne']
    rw [h1, iteratedFDerivWithin_of_isOpen i isOpen_Ioi (show (1:ℝ) ∈ Ioi 0 by norm_num),
      norm_iteratedFDeriv_eq_norm_iteratedDeriv, Real.norm_eq_abs]
    exact log_derivs i (hi.trans hn)

lemma choose_sum (n : ℕ) : ∑ i ∈ Finset.range (n + 1), (n.choose i : ℝ) = 2 ^ n := by
  exact_mod_cast Nat.sum_range_choose n

lemma nd_mul {f g : ℝ → ℝ} (hf : ContDiffOn ℝ 3 f (Ioi 0)) (hg : ContDiffOn ℝ 3 g (Ioi 0))
    {p : ℝ} (hp : 0 < p) {n : ℕ} (hn : n ≤ 3) {a b D : ℝ}
    (hfD : ∀ i ≤ n, nd i f p ≤ a * D ^ i) (hgD : ∀ i ≤ n, nd i g p ≤ b * D ^ i) :
    nd n (fun q => f q * g q) p ≤ 2 ^ n * (a * b) * D ^ n := by
  have := norm_iteratedFDerivWithin_mul_le (N := 3) hf hg uS hp (n := n) (by exact_mod_cast hn)
  rw [nd_F] at this
  simp only [nd_F] at this
  refine this.trans ?_
  calc ∑ i ∈ Finset.range (n + 1), (n.choose i : ℝ) * nd i f p * nd (n - i) g p
      ≤ ∑ i ∈ Finset.range (n + 1), (n.choose i : ℝ) * (a * b * D ^ n) := by
        refine Finset.sum_le_sum fun i hi => ?_
        have hi' : i ≤ n := Nat.lt_succ_iff.1 (Finset.mem_range.1 hi)
        have e : a * b * D ^ n = (a * D ^ i) * (b * D ^ (n - i)) := by
          rw [show D ^ n = D ^ i * D ^ (n - i) by rw [← pow_add, Nat.add_sub_cancel' hi']]; ring
        rw [mul_assoc, e]
        exact mul_le_mul_of_nonneg_left (mul_le_mul (hfD i hi') (hgD _ (Nat.sub_le _ _))
          (nd_nonneg _ _ _) ((nd_nonneg _ _ _).trans (hfD i hi'))) (Nat.cast_nonneg _)
    _ = 2 ^ n * (a * b) * D ^ n := by rw [← Finset.sum_mul, choose_sum]; ring

lemma nd_add {f g : ℝ → ℝ} (hf : ContDiffOn ℝ 3 f (Ioi 0)) (hg : ContDiffOn ℝ 3 g (Ioi 0))
    {p : ℝ} (hp : 0 < p) {n : ℕ} (hn : n ≤ 3) :
    nd n (fun q => f q + g q) p ≤ nd n f p + nd n g p := by
  unfold nd
  rw [iteratedDerivWithin_fun_add hp uS ((hf p hp).of_le (by exact_mod_cast hn))
    ((hg p hp).of_le (by exact_mod_cast hn))]
  exact abs_add_le _ _

lemma nd_cmul (c : ℝ) (f : ℝ → ℝ) (n : ℕ) (p : ℝ) :
    nd n (fun q => c * f q) p = |c| * nd n f p := by
  unfold nd; rw [iteratedDerivWithin_const_mul_field, abs_mul]

lemma nd_base (a c : ℝ) (n : ℕ) {p : ℝ} (hp : 0 < p) :
    nd n (fun q => a * Real.exp (c * q)) p = |a| * |c| ^ n * Real.exp (c * p) := by
  rw [nd_cmul]
  unfold nd
  rw [iteratedDerivWithin_of_isOpen isOpen_Ioi hp, iteratedDeriv_exp_const_mul, abs_mul, abs_pow,
    abs_of_pos (exp_pos _)]
  ring

lemma nd_inv_part (i : ℕ) (hi : i ≤ 3) {p : ℝ} (hp : 1 / 2 ≤ p) :
    nd i (fun q => -(1 + q⁻¹)) p ≤ 96 := by
  have hp0 : 0 < p := by linarith
  have e : (fun q : ℝ => -(1 + q⁻¹)) = fun q => (-1) * (1 + q⁻¹) := by ext q; ring
  rw [e, nd_cmul]
  simp only [abs_neg, abs_one, one_mul]
  unfold nd
  rw [iteratedDerivWithin_of_isOpen isOpen_Ioi hp0]
  have hinv : p⁻¹ ≤ 2 := by rw [inv_le_comm₀ hp0 (by norm_num)]; linarith
  rcases Nat.eq_zero_or_pos i with rfl | hi0
  · simp only [iteratedDeriv_zero]
    rw [abs_le]; constructor <;> nlinarith [inv_pos.2 hp0]
  · rw [iteratedDeriv_const_add hi0, iteratedDeriv_eq_iterate, iter_deriv_inv, abs_mul, abs_mul,
      abs_pow, abs_neg, abs_one, one_pow, one_mul, Nat.abs_cast]
    have hz : |p ^ (-1 - (i : ℤ))| = (p⁻¹) ^ (i + 1) := by
      rw [abs_of_pos (zpow_pos hp0 _), show (-1 - (i : ℤ)) = -((i + 1 : ℕ) : ℤ) by push_cast; ring,
        zpow_neg, zpow_natCast, inv_pow]
    rw [hz]
    have h1 : (p⁻¹) ^ (i + 1) ≤ 2 ^ (i + 1) := pow_le_pow_left₀ (inv_pos.2 hp0).le hinv _
    have h2 : (2:ℝ) ^ (i + 1) ≤ 16 := by
      calc (2:ℝ) ^ (i + 1) ≤ 2 ^ 4 := pow_le_pow_right₀ (by norm_num) (by omega)
        _ = 16 := by norm_num
    have h3 : (i.factorial : ℝ) ≤ 6 := by
      have : i.factorial ≤ Nat.factorial 3 := Nat.factorial_le hi
      exact_mod_cast this.trans (by decide)
    calc (i.factorial : ℝ) * (p⁻¹) ^ (i + 1) ≤ 6 * 16 :=
          mul_le_mul h3 (h1.trans h2) (by positivity) (by norm_num)
      _ = 96 := by norm_num

section Core
variable (w α β : Fin 3 → ℝ) (E0 E1 Y0 Y1 : ℝ)

/-- `F(p) = Σ_j w_j α_j^p`. -/
noncomputable def cF (q : ℝ) : ℝ := ∑ j, w j * Real.exp (Real.log (α j) * q)
/-- `S(p) = Σ_j w_j β_j α_j^{p-1}`. -/
noncomputable def cS (q : ℝ) : ℝ := ∑ j, (w j * β j / α j) * Real.exp (Real.log (α j) * q)
/-- `Σ_j |w_j β_j| α_j^{p-1}`. -/
noncomputable def cSt (q : ℝ) : ℝ := ∑ j, |w j * β j / α j| * Real.exp (Real.log (α j) * q)
/-- `|F E - S Y|²`. -/
noncomputable def cQ (q : ℝ) : ℝ :=
  (cF w α q * E0 - cS w α β q * Y0) ^ 2 + (cF w α q * E1 - cS w α β q * Y1) ^ 2
/-- `log G`. -/
noncomputable def ch (q : ℝ) : ℝ :=
  -(1 + q⁻¹) * Real.log (cF w α q) + (1 / 2) * Real.log (cQ w α β E0 E1 Y0 Y1 q)
/-- `G = F^{-1-1/p} √Q`. -/
noncomputable def cG (q : ℝ) : ℝ := Real.exp (ch w α β E0 E1 Y0 Y1 q)

lemma cF_eq : cF w α = fun q => w 0 * Real.exp (Real.log (α 0) * q) +
    w 1 * Real.exp (Real.log (α 1) * q) + w 2 * Real.exp (Real.log (α 2) * q) := by
  ext q; simp [cF, Fin.sum_univ_three]

lemma cS_eq : cS w α β = fun q => (w 0 * β 0 / α 0) * Real.exp (Real.log (α 0) * q) +
    (w 1 * β 1 / α 1) * Real.exp (Real.log (α 1) * q) +
    (w 2 * β 2 / α 2) * Real.exp (Real.log (α 2) * q) := by
  ext q; simp [cS, Fin.sum_univ_three]

lemma cF_cd : ContDiff ℝ 3 (cF w α) := by rw [cF_eq]; fun_prop
lemma cS_cd : ContDiff ℝ 3 (cS w α β) := by rw [cS_eq]; fun_prop

variable {w α β E0 E1 Y0 Y1}

lemma cF_pos (hw : ∀ j, 0 < w j) (q : ℝ) : 0 < cF w α q := by
  rw [cF_eq]; have := hw 0; have := hw 1; have := hw 2; positivity

lemma cQ_pos (hw : ∀ j, 0 < w j) (hdet : E0 * Y1 - E1 * Y0 ≠ 0) (q : ℝ) :
    0 < cQ w α β E0 E1 Y0 Y1 q := by
  unfold cQ
  set F := cF w α q
  set S := cS w α β q
  have hF := cF_pos (α := α) hw q
  by_contra h
  have h0 : F * E0 - S * Y0 = 0 := by nlinarith [sq_nonneg (F * E0 - S * Y0), sq_nonneg (F * E1 - S * Y1)]
  have h1 : F * E1 - S * Y1 = 0 := by nlinarith [sq_nonneg (F * E0 - S * Y0), sq_nonneg (F * E1 - S * Y1)]
  apply hdet
  have : F * (E0 * Y1 - E1 * Y0) = 0 := by linear_combination Y1 * h0 - Y0 * h1
  rcases mul_eq_zero.1 this with h | h
  · exact absurd h hF.ne'
  · exact h

lemma nd_three {f0 f1 f2 : ℝ → ℝ} (h0 : ContDiff ℝ 3 f0) (h1 : ContDiff ℝ 3 f1)
    (h2 : ContDiff ℝ 3 f2) {p : ℝ} (hp : 0 < p) {n : ℕ} (hn : n ≤ 3) :
    nd n (fun q => f0 q + f1 q + f2 q) p ≤ nd n f0 p + nd n f1 p + nd n f2 p := by
  refine (nd_add (h0.add h1).contDiffOn h2.contDiffOn hp hn).trans ?_
  have := nd_add h0.contDiffOn h1.contDiffOn hp hn
  linarith

lemma nd_cF (hw : ∀ j, 0 < w j) {D : ℝ} (hD : ∀ j, |Real.log (α j)| ≤ D) {p : ℝ} (hp : 0 < p)
    {n : ℕ} (hn : n ≤ 3) : nd n (cF w α) p ≤ cF w α p * D ^ n := by
  rw [cF_eq]
  refine (nd_three (by fun_prop) (by fun_prop) (by fun_prop) hp hn).trans ?_
  simp only [nd_base _ _ _ hp]
  have k : ∀ j, |w j| * |Real.log (α j)| ^ n * Real.exp (Real.log (α j) * p) ≤
      w j * Real.exp (Real.log (α j) * p) * D ^ n := by
    intro j
    rw [abs_of_pos (hw j)]
    have := pow_le_pow_left₀ (abs_nonneg _) (hD j) n
    have := hw j
    calc w j * |Real.log (α j)| ^ n * Real.exp (Real.log (α j) * p)
        ≤ w j * D ^ n * Real.exp (Real.log (α j) * p) := by gcongr
      _ = _ := by ring
  linarith [k 0, k 1, k 2]

lemma nd_cS {D : ℝ} (hD : ∀ j, |Real.log (α j)| ≤ D) {p : ℝ} (hp : 0 < p)
    {n : ℕ} (hn : n ≤ 3) : nd n (cS w α β) p ≤ cSt w α β p * D ^ n := by
  rw [cS_eq]
  refine (nd_three (by fun_prop) (by fun_prop) (by fun_prop) hp hn).trans ?_
  simp only [nd_base _ _ _ hp]
  have k : ∀ j, |w j * β j / α j| * |Real.log (α j)| ^ n * Real.exp (Real.log (α j) * p) ≤
      |w j * β j / α j| * Real.exp (Real.log (α j) * p) * D ^ n := by
    intro j
    have := pow_le_pow_left₀ (abs_nonneg _) (hD j) n
    calc |w j * β j / α j| * |Real.log (α j)| ^ n * Real.exp (Real.log (α j) * p)
        ≤ |w j * β j / α j| * D ^ n * Real.exp (Real.log (α j) * p) := by gcongr
      _ = _ := by ring
  have e : cSt w α β p = |w 0 * β 0 / α 0| * Real.exp (Real.log (α 0) * p) +
      |w 1 * β 1 / α 1| * Real.exp (Real.log (α 1) * p) +
      |w 2 * β 2 / α 2| * Real.exp (Real.log (α 2) * p) := by simp [cSt, Fin.sum_univ_three]
  rw [e]
  linarith [k 0, k 1, k 2]

lemma cSt_nonneg (q : ℝ) : 0 ≤ cSt w α β q := by
  unfold cSt; exact Finset.sum_nonneg fun j _ => by positivity

lemma nd_cQ (hw : ∀ j, 0 < w j) {D : ℝ} (hD : ∀ j, |Real.log (α j)| ≤ D) {p : ℝ} (hp : 0 < p)
    {n : ℕ} (hn : n ≤ 3) :
    nd n (cQ w α β E0 E1 Y0 Y1) p ≤ 2 ^ n * ((|E0 ^ 2 + E1 ^ 2| + 2 * |E0 * Y0 + E1 * Y1| +
      |Y0 ^ 2 + Y1 ^ 2|) * (cF w α p + cSt w α β p) ^ 2) * D ^ n := by
  set F := cF w α
  set S := cS w α β
  have e : cQ w α β E0 E1 Y0 Y1 = fun q => (E0 ^ 2 + E1 ^ 2) * (F q * F q) +
      (-2 * (E0 * Y0 + E1 * Y1)) * (F q * S q) + (Y0 ^ 2 + Y1 ^ 2) * (S q * S q) := by
    ext q; simp only [cQ, F, S]; ring
  have hFc : ContDiff ℝ 3 F := cF_cd w α
  have hSc : ContDiff ℝ 3 S := cS_cd w α β
  rw [e]
  refine (nd_three (by fun_prop) (by fun_prop) (by fun_prop) hp hn).trans ?_
  simp only [nd_cmul]
  have hF : ∀ i ≤ n, nd i F p ≤ F p * D ^ i := fun i hi => nd_cF hw hD hp (hi.trans hn)
  have hS : ∀ i ≤ n, nd i S p ≤ cSt w α β p * D ^ i := fun i hi => nd_cS hD hp (hi.trans hn)
  have a1 := nd_mul hFc.contDiffOn hFc.contDiffOn hp hn hF hF
  have a2 := nd_mul hFc.contDiffOn hSc.contDiffOn hp hn hF hS
  have a3 := nd_mul hSc.contDiffOn hSc.contDiffOn hp hn hS hS
  have hFp : 0 < F p := cF_pos hw p
  have hSt := cSt_nonneg (w := w) (α := α) (β := β) p
  have hDn : 0 ≤ (2:ℝ) ^ n * D ^ n := by
    have : 0 ≤ nd 0 F p := nd_nonneg _ _ _
    have h0 := hF 0 (Nat.zero_le _)
    simp only [pow_zero, mul_one, nd_zero] at h0
    have hD0 : 0 ≤ D := (abs_nonneg _).trans (hD 0)
    positivity
  have := abs_nonneg (E0 ^ 2 + E1 ^ 2)
  have := abs_nonneg (E0 * Y0 + E1 * Y1)
  have := abs_nonneg (Y0 ^ 2 + Y1 ^ 2)
  rw [abs_mul, abs_neg, abs_two] at *
  calc |E0 ^ 2 + E1 ^ 2| * nd n (fun q => F q * F q) p +
        2 * |E0 * Y0 + E1 * Y1| * nd n (fun q => F q * S q) p +
        |Y0 ^ 2 + Y1 ^ 2| * nd n (fun q => S q * S q) p
      ≤ |E0 ^ 2 + E1 ^ 2| * (2 ^ n * (F p * F p) * D ^ n) +
        2 * |E0 * Y0 + E1 * Y1| * (2 ^ n * (F p * cSt w α β p) * D ^ n) +
        |Y0 ^ 2 + Y1 ^ 2| * (2 ^ n * (cSt w α β p * cSt w α β p) * D ^ n) := by gcongr
    _ ≤ _ := by
        have e2 : ∀ X : ℝ, 2 ^ n * X * D ^ n = X * (2 ^ n * D ^ n) := fun X => by ring
        rw [e2, e2, e2, e2]
        have : |E0 ^ 2 + E1 ^ 2| * (F p * F p) + 2 * |E0 * Y0 + E1 * Y1| * (F p * cSt w α β p) +
            |Y0 ^ 2 + Y1 ^ 2| * (cSt w α β p * cSt w α β p) ≤
            (|E0 ^ 2 + E1 ^ 2| + 2 * |E0 * Y0 + E1 * Y1| + |Y0 ^ 2 + Y1 ^ 2|) *
              (F p + cSt w α β p) ^ 2 := by
          nlinarith [mul_nonneg hFp.le hSt, mul_nonneg hSt hSt, mul_nonneg hFp.le hFp.le]
        nlinarith

/-- the constant `M` of `nd_cQ`. -/
noncomputable def cM (E0 E1 Y0 Y1 : ℝ) : ℝ :=
  |E0 ^ 2 + E1 ^ 2| + 2 * |E0 * Y0 + E1 * Y1| + |Y0 ^ 2 + Y1 ^ 2|

lemma cM_nonneg : 0 ≤ cM E0 E1 Y0 Y1 := by unfold cM; positivity

lemma cQ_cd (hw : ∀ j, 0 < w j) : ContDiff ℝ 3 (cQ w α β E0 E1 Y0 Y1) := by
  have := cF_cd w α; have := cS_cd w α β
  unfold cQ; fun_prop

lemma nd_cQ_rel (hw : ∀ j, 0 < w j) {D : ℝ} (hD1 : 1 ≤ D) (hD : ∀ j, |Real.log (α j)| ≤ D)
    {p : ℝ} (hp : 0 < p) {K0 : ℝ} (hK0 : 0 ≤ K0)
    (hQ : (cF w α p + cSt w α β p) ^ 2 ≤ K0 * cQ w α β E0 E1 Y0 Y1 p)
    {i : ℕ} (hi1 : 1 ≤ i) (hi : i ≤ 3) :
    nd i (cQ w α β E0 E1 Y0 Y1) p ≤
      cQ w α β E0 E1 Y0 Y1 p * ((2 * (cM E0 E1 Y0 Y1 * K0 + 1)) * D) ^ i := by
  refine (nd_cQ hw hD hp hi).trans ?_
  have hM := cM_nonneg (E0 := E0) (E1 := E1) (Y0 := Y0) (Y1 := Y1)
  set Q := cQ w α β E0 E1 Y0 Y1 p
  set X := (cF w α p + cSt w α β p) ^ 2
  set M := cM E0 E1 Y0 Y1
  have h1 : M * K0 ≤ (M * K0 + 1) ^ i := by
    calc M * K0 ≤ M * K0 + 1 := by linarith
      _ = (M * K0 + 1) ^ 1 := (pow_one _).symm
      _ ≤ (M * K0 + 1) ^ i := pow_le_pow_right₀ (by nlinarith) hi1
  have hD0 : 0 ≤ D ^ i := by positivity
  rw [mul_pow, mul_pow]
  change 2 ^ i * (M * X) * D ^ i ≤ _
  have hQ0 : 0 ≤ Q := by simp only [Q, cQ]; positivity
  have : M * X ≤ Q * (M * K0 + 1) ^ i := by
    calc M * X ≤ M * (K0 * Q) := mul_le_mul_of_nonneg_left hQ hM
      _ = Q * (M * K0) := by ring
      _ ≤ Q * (M * K0 + 1) ^ i := mul_le_mul_of_nonneg_left h1 hQ0
  have h2 : (0:ℝ) ≤ 2 ^ i := by positivity
  calc 2 ^ i * (M * X) * D ^ i ≤ 2 ^ i * (Q * (M * K0 + 1) ^ i) * D ^ i := by gcongr
    _ = _ := by ring

lemma ch_cd (hw : ∀ j, 0 < w j) (hdet : E0 * Y1 - E1 * Y0 ≠ 0) :
    ContDiffOn ℝ 3 (ch w α β E0 E1 Y0 Y1) (Ioi 0) := by
  unfold ch
  have hF := cF_cd w α
  have hQ := cQ_cd (α := α) (β := β) (E0 := E0) (E1 := E1) (Y0 := Y0) (Y1 := Y1) hw
  refine ContDiffOn.add ?_ ?_
  · refine ContDiffOn.mul ?_ ?_
    · exact (contDiffOn_const.add (contDiffOn_inv ℝ |>.mono fun x hx => ne_of_gt hx)).neg
    · exact (hF.log fun q => (cF_pos hw q).ne').contDiffOn
  · exact contDiffOn_const.mul (hQ.log fun q => (cQ_pos hw hdet q).ne').contDiffOn

/-- the constant `K₃` in the bound for `log G`. -/
noncomputable def cK3 (b K2 : ℝ) : ℝ := 8 * 96 * (|b| + 12) + 6 * K2 ^ 3 + 1

lemma nd_ch (hw : ∀ j, 0 < w j) (hdet : E0 * Y1 - E1 * Y0 ≠ 0) {D : ℝ} (hD1 : 1 ≤ D)
    (hD : ∀ j, |Real.log (α j)| ≤ D) {p : ℝ} (hp : 1 / 2 ≤ p) {K0 : ℝ} (hK0 : 0 ≤ K0)
    (hQ : (cF w α p + cSt w α β p) ^ 2 ≤ K0 * cQ w α β E0 E1 Y0 Y1 p)
    {n : ℕ} (hn1 : 1 ≤ n) (hn : n ≤ 3) :
    nd n (ch w α β E0 E1 Y0 Y1) p ≤
      cK3 (Real.log (cF w α p)) (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) * D ^ n := by
  have hp0 : 0 < p := by linarith
  set K2 := 2 * (cM E0 E1 Y0 Y1 * K0 + 1)
  have hK2 : 1 ≤ K2 := by have := cM_nonneg (E0 := E0) (E1 := E1) (Y0 := Y0) (Y1 := Y1); nlinarith
  have hF := cF_cd w α
  have hQc := cQ_cd (α := α) (β := β) (E0 := E0) (E1 := E1) (Y0 := Y0) (Y1 := Y1) hw
  have hlF : ContDiffOn ℝ 3 (fun q => Real.log (cF w α q)) (Ioi 0) :=
    (hF.log fun q => (cF_pos hw q).ne').contDiffOn
  have hlQ : ContDiffOn ℝ 3 (fun q => Real.log (cQ w α β E0 E1 Y0 Y1 q)) (Ioi 0) :=
    (hQc.log fun q => (cQ_pos hw hdet q).ne').contDiffOn
  have hinv : ContDiffOn ℝ 3 (fun q : ℝ => -(1 + q⁻¹)) (Ioi 0) :=
    (contDiffOn_const.add (contDiffOn_inv ℝ |>.mono fun x hx => ne_of_gt hx)).neg
  have hDpow : ∀ i : ℕ, 1 ≤ D ^ i := fun i => one_le_pow₀ hD1
  unfold ch
  refine (nd_add (hinv.mul hlF) (contDiffOn_const.mul hlQ) hp0 hn).trans ?_
  rw [nd_cmul]
  -- first part
  have hA : ∀ i ≤ n, nd i (fun q : ℝ => -(1 + q⁻¹)) p ≤ 96 * D ^ i := fun i hi => by
    have := nd_inv_part i (hi.trans hn) hp
    nlinarith [hDpow i]
  have hB : ∀ i ≤ n, nd i (fun q => Real.log (cF w α q)) p ≤ (|Real.log (cF w α p)| + 12) * D ^ i := by
    intro i hi
    rcases Nat.eq_zero_or_pos i with rfl | hi0
    · rw [nd_zero]; simp
    · have := nd_log hF.contDiffOn (fun q _ => cF_pos hw q) hp0 hi0 (hi.trans hn)
        (fun j _ hj => nd_cF hw hD hp0 (hj.trans (hi.trans hn)))
      have hf : (i.factorial : ℝ) ≤ 6 := by
        have : i.factorial ≤ Nat.factorial 3 := Nat.factorial_le (hi.trans hn)
        exact_mod_cast this.trans (by decide)
      have : (i.factorial : ℝ) * 2 * D ^ i ≤ 12 * D ^ i := by nlinarith [hDpow i]
      nlinarith [abs_nonneg (Real.log (cF w α p)), hDpow i]
  have h1 := nd_mul hinv hlF hp0 hn hA hB
  have h2 := nd_log hQc.contDiffOn (fun q _ => cQ_pos hw hdet q) hp0 hn1 hn
    (fun i hi1 hi => nd_cQ_rel hw hD1 hD hp0 hK0 hQ hi1 (hi.trans hn))
  have hf : (n.factorial : ℝ) ≤ 6 := by
    have : n.factorial ≤ Nat.factorial 3 := Nat.factorial_le hn
    exact_mod_cast this.trans (by decide)
  have h2' : |(1:ℝ) / 2| * nd n (fun q => Real.log (cQ w α β E0 E1 Y0 Y1 q)) p ≤
      6 * K2 ^ 3 * D ^ n := by
    rw [mul_pow] at h2
    have hk : K2 ^ n ≤ K2 ^ 3 := pow_le_pow_right₀ hK2 hn
    have : (n.factorial : ℝ) * 2 * (K2 ^ n * D ^ n) ≤ 12 * K2 ^ 3 * D ^ n := by
      have := hDpow n
      have : 0 ≤ K2 ^ n := by positivity
      nlinarith [mul_le_mul hf hk (by positivity) (by norm_num)]
    rw [abs_of_pos (by norm_num : (0:ℝ) < 1 / 2)]
    linarith
  have h1' : 2 ^ n * (96 * (|Real.log (cF w α p)| + 12)) * D ^ n ≤
      8 * 96 * (|Real.log (cF w α p)| + 12) * D ^ n := by
    have : (2:ℝ) ^ n ≤ 8 := by
      calc (2:ℝ) ^ n ≤ 2 ^ 3 := pow_le_pow_right₀ (by norm_num) hn
        _ = 8 := by norm_num
    have hX : 0 ≤ 96 * (|Real.log (cF w α p)| + 12) * D ^ n := by
      have := hDpow n; positivity
    calc 2 ^ n * (96 * (|Real.log (cF w α p)| + 12)) * D ^ n
        = 2 ^ n * (96 * (|Real.log (cF w α p)| + 12) * D ^ n) := by ring
      _ ≤ 8 * (96 * (|Real.log (cF w α p)| + 12) * D ^ n) := mul_le_mul_of_nonneg_right this hX
      _ = _ := by ring
  unfold cK3
  nlinarith [hDpow n]

lemma cK3_ge (b K2 : ℝ) : 1 ≤ cK3 b K2 ∨ True := Or.inr trivial

lemma nd_cG (hw : ∀ j, 0 < w j) (hdet : E0 * Y1 - E1 * Y0 ≠ 0) {D : ℝ} (hD1 : 1 ≤ D)
    (hD : ∀ j, |Real.log (α j)| ≤ D) {p : ℝ} (hp : 1 / 2 ≤ p) {K0 : ℝ} (hK0 : 0 ≤ K0)
    (hQ : (cF w α p + cSt w α β p) ^ 2 ≤ K0 * cQ w α β E0 E1 Y0 Y1 p)
    {n : ℕ} (hn : n ≤ 3) :
    nd n (cG w α β E0 E1 Y0 Y1) p ≤ 6 * cG w α β E0 E1 Y0 Y1 p *
      (cK3 (Real.log (cF w α p)) (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) * D) ^ n := by
  have hp0 : 0 < p := by linarith
  set K3 := cK3 (Real.log (cF w α p)) (2 * (cM E0 E1 Y0 Y1 * K0 + 1))
  have hK3 : 1 ≤ K3 := by
    have : 0 ≤ (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) := by
      have := cM_nonneg (E0 := E0) (E1 := E1) (Y0 := Y0) (Y1 := Y1); nlinarith
    have h3 : 0 ≤ (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) ^ 3 := pow_nonneg this 3
    simp only [K3, cK3]
    nlinarith [abs_nonneg (Real.log (cF w α p))]
  have := nd_exp (ch_cd hw hdet) hp0 hn (D := K3 * D) (fun i hi1 hi => by
    refine (nd_ch hw hdet hD1 hD hp hK0 hQ hi1 (hi.trans hn)).trans ?_
    rw [mul_pow]
    have : K3 ≤ K3 ^ i := by
      calc K3 = K3 ^ 1 := (pow_one _).symm
        _ ≤ K3 ^ i := pow_le_pow_right₀ hK3 hi1
    have : 0 ≤ D ^ i := by positivity
    nlinarith)
  have hf : (n.factorial : ℝ) ≤ 6 := by
    have : n.factorial ≤ Nat.factorial 3 := Nat.factorial_le hn
    exact_mod_cast this.trans (by decide)
  unfold cG
  refine this.trans ?_
  have h0 : 0 ≤ (K3 * D) ^ n := pow_nonneg (mul_nonneg (by linarith) (by linarith)) n
  have := exp_pos (ch w α β E0 E1 Y0 Y1 p)
  exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right hf this.le) h0

end Core

end S05

namespace S05
open Set Real Filter Topology MeasureTheory Interval

/-! ### The length element on one edge -/

/-- weights in the edge order `j₀, j₁, j₂`. -/
noncomputable def ew (H : Hex) (k : Fin 3) : Fin 3 → ℝ := ![H.w k, H.w (k + 1), H.w (k + 2)]
/-- `|⟨u_j, y(t)⟩|` in the edge order. -/
noncomputable def eα (H : Hex) (k : Fin 3) (t : ℝ) : Fin 3 → ℝ :=
  ![A0 H k * t / H.w k, A1 H k * (1 - t) / H.w (k + 1),
    Bv H k * (1 + (rr H k - 1) * t) / H.w (k + 2)]
/-- their `t`-derivatives. -/
noncomputable def eβ (H : Hex) (k : Fin 3) : Fin 3 → ℝ :=
  ![A0 H k / H.w k, -(A1 H k / H.w (k + 1)), Bv H k * (rr H k - 1) / H.w (k + 2)]
/-- `V_1 - V_0`. -/
noncomputable def eE (H : Hex) (k : Fin 3) : X04.Plane := eV1 H k - eV0 H k

/-- the explicit length element `F^{-1-1/p}|F(V_1-V_0) - S y(t)|`. -/
noncomputable def Gf (H : Hex) (k : Fin 3) (p t : ℝ) : ℝ :=
  cG (ew H k) (eα H k t) (eβ H k) (eE H k 0) (eE H k 1) (ept H k t 0) (ept H k t 1) p

lemma sum_rot (f : Fin 3 → ℝ) (k : Fin 3) : ∑ j, f j = f k + f (k + 1) + f (k + 2) := by
  rw [Fin.sum_univ_three]
  fin_cases k <;> simp <;> ring_nf

lemma ew_pos (H : Hex) (hH : IsHex H) (k : Fin 3) (j : Fin 3) : 0 < ew H k j := by
  fin_cases j <;> simp [ew, hH.1 _]

lemma eα_pos (H : Hex) (hH : IsHex H) (k : Fin 3) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) (j : Fin 3) :
    0 < eα H k t j := by
  obtain ⟨-, -, -, hB, hB1, hA0, hA1, hr⟩ := S05.proof_8 H hH k
  have hw := hH.1
  have h2 : 0 < 1 + (rr H k - 1) * t := by nlinarith [ht.1, ht.2]
  have := ht.1; have : 0 < 1 - t := by linarith [ht.2]
  fin_cases j <;> simp [eα] <;> have := hw k <;> have := hw (k + 1) <;> have := hw (k + 2) <;>
    positivity

lemma abs_dot_edge (H : Hex) (hH : IsHex H) (k : Fin 3) {t : ℝ} (ht : t ∈ Icc (0:ℝ) 1) :
    |X04.dot (H.u k) (ept H k t)| = eα H k t 0 ∧
    |X04.dot (H.u (k + 1)) (ept H k t)| = eα H k t 1 ∧
    |X04.dot (H.u (k + 2)) (ept H k t)| = eα H k t 2 := by
  obtain ⟨h8, -⟩ := S05.proof_8 H hH k
  obtain ⟨e0, e1, e2⟩ := h8 t ht
  simp only [term] at e0 e1 e2
  have hw := hH.1
  simp only [eα, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.head_cons, Matrix.tail_cons]
  refine ⟨?_, ?_, ?_⟩ <;> rw [eq_div_iff (hw _).ne'] <;> linarith

lemma Fw_edge (H : Hex) (hH : IsHex H) (k : Fin 3) (p : ℝ) {t : ℝ} (ht : t ∈ Icc (0:ℝ) 1) :
    Fw H.u H.w p (ept H k t) = ∑ j, ew H k j * eα H k t j ^ p := by
  obtain ⟨a0, a1, a2⟩ := abs_dot_edge H hH k ht
  rw [Fw, sum_rot _ k, Fin.sum_univ_three, a0, a1, a2]
  simp [ew]

lemma Fw_edge_cF (H : Hex) (hH : IsHex H) (k : Fin 3) (p : ℝ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    Fw H.u H.w p (ept H k t) = cF (ew H k) (eα H k t) p := by
  rw [Fw_edge H hH k p (Ioo_subset_Icc_self ht), cF]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [Real.rpow_def_of_pos (eα_pos H hH k ht j)]

lemma hasDerivAt_eα (H : Hex) (k : Fin 3) (t : ℝ) (j : Fin 3) :
    HasDerivAt (fun s => eα H k s j) (eβ H k j) t := by
  fin_cases j <;> simp only [eα, eβ, Fin.zero_eta, Fin.mk_one, Fin.reduceFinMk,
    Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
    Matrix.tail_cons]
  · simpa using ((hasDerivAt_id t).const_mul (A0 H k)).div_const (H.w k)
  · exact ((((hasDerivAt_id t).const_sub 1).const_mul (A1 H k)).div_const
      (H.w (k + 1))).congr_deriv (by ring)
  · exact (((((hasDerivAt_id t).const_mul (rr H k - 1)).const_add 1).const_mul
      (Bv H k)).div_const (H.w (k + 2))).congr_deriv (by ring)

lemma hasDerivAt_ept (H : Hex) (k : Fin 3) (t : ℝ) : HasDerivAt (ept H k) (eE H k) t := by
  show HasDerivAt (fun s => eV0 H k + s • (eV1 H k - eV0 H k)) _ t
  exact (((hasDerivAt_id t).smul_const _).const_add _).congr_deriv (by simp [eE])

lemma edge_det_pos (H : Hex) (hH : IsHex H) (k : Fin 3) : 0 < det2 (eV0 H k) (eV1 H k) := by
  obtain ⟨-, -, -, -, h01, h12, h02⟩ := hH
  fin_cases k
  · simpa [eV0, eV1] using h01
  · simpa [eV0, eV1] using h12
  · simp only [eV0, eV1, det2]; simp only [det2] at h02
    simp; linarith

lemma edge_det_ne (H : Hex) (hH : IsHex H) (k : Fin 3) (t : ℝ) :
    eE H k 0 * ept H k t 1 - eE H k 1 * ept H k t 0 ≠ 0 := by
  have := edge_det_pos H hH k
  simp only [det2] at this
  simp only [eE, ept, PiLp.add_apply, PiLp.sub_apply, PiLp.smul_apply, smul_eq_mul]
  intro h; nlinarith

lemma sqrt_Q_eq (a E0 E1 y0 y1 : ℝ) :
    Real.sqrt ((a * E0) ^ 2 + (a * E1) ^ 2) = |a| * Real.sqrt (E0 ^ 2 + E1 ^ 2) := by
  rw [show (a * E0) ^ 2 + (a * E1) ^ 2 = a ^ 2 * (E0 ^ 2 + E1 ^ 2) by ring, Real.sqrt_mul (sq_nonneg _),
    Real.sqrt_sq_eq_abs]

lemma arc_deriv (H : Hex) (hH : IsHex H) (k : Fin 3) {p : ℝ} (hp : 0 < p) {t : ℝ}
    (ht : t ∈ Ioo (0:ℝ) 1) :
    HasDerivAt (arc H k p) (deriv (arc H k p) t) t ∧ Gl H k p t = Gf H k p t := by
  set Φ : ℝ → ℝ := fun s => ∑ j, ew H k j * eα H k s j ^ p with hΦ
  have hev : arc H k p =ᶠ[𝓝 t] fun s => Φ s ^ (-(1 / p)) • ept H k s := by
    filter_upwards [isOpen_Ioo.mem_nhds ht] with s hs
    simp only [arc, hΦ, Fw_edge H hH k p (Ioo_subset_Icc_self hs)]
  have hα := eα_pos H hH k ht
  set S := cS (ew H k) (eα H k t) (eβ H k) p
  have hΦd : HasDerivAt Φ (p * S) t := by
    have := HasDerivAt.fun_sum (u := Finset.univ) (fun j _ =>
      ((hasDerivAt_eα H k t j).rpow_const (p := p) (Or.inl (hα j).ne')).const_mul (ew H k j))
    convert this using 1
    simp only [S, cS, Finset.mul_sum]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [Real.rpow_sub_one (hα j).ne', Real.rpow_def_of_pos (hα j)]
    field_simp
  have hΦt : Φ t = cF (ew H k) (eα H k t) p := by
    simp only [hΦ, cF]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [Real.rpow_def_of_pos (hα j)]
  have hFpos : 0 < Φ t := by rw [hΦt]; exact cF_pos (ew_pos H hH k) p
  have hd := (hΦd.rpow_const (p := -(1 / p)) (Or.inl hFpos.ne')).smul (hasDerivAt_ept H k t)
  have hd' := hd.congr_of_eventuallyEq hev
  have hD : deriv (arc H k p) t = _ := hd'.deriv
  refine ⟨hD ▸ hd', ?_⟩
  rw [Gl, hD]
  set F := Φ t
  have e1 : F ^ (-(1 / p)) = F * F ^ (-(1 / p) - 1) := by
    rw [Real.rpow_sub_one hFpos.ne']; field_simp
  set c := F ^ (-(1 / p) - 1)
  have hc : 0 < c := Real.rpow_pos_of_pos hFpos _
  have hv : F ^ (-(1 / p)) • eE H k + (p * S * -(1 / p) * c) • ept H k t =
      c • (F • eE H k - S • ept H k t) := by
    rw [e1]
    ext i
    simp only [PiLp.add_apply, PiLp.smul_apply, PiLp.sub_apply, smul_eq_mul]
    field_simp
    ring
  rw [hv, norm_smul, Real.norm_eq_abs, abs_of_pos hc, X04.norm_plane]
  simp only [PiLp.sub_apply, PiLp.smul_apply, smul_eq_mul]
  rw [Gf, cG, ch, Real.exp_add, ← hΦt]
  have hQ : 0 < cQ (ew H k) (eα H k t) (eβ H k) (eE H k 0) (eE H k 1) (ept H k t 0)
      (ept H k t 1) p := by
    exact cQ_pos (ew_pos H hH k) (edge_det_ne H hH k t) p
  rw [mul_comm (-(1 + p⁻¹)) (Real.log F), ← Real.rpow_def_of_pos hFpos,
    mul_comm (1 / 2 : ℝ), ← Real.rpow_def_of_pos hQ, ← Real.sqrt_eq_rpow]
  have : c = F ^ (-(1 + p⁻¹)) := by simp only [c]; congr 1; ring
  rw [this, cQ, ← hΦt]

/-! ### uniform bounds -/

lemma St_le {x0 x1 x2 t c c2 : ℝ} (h0 : 0 ≤ x0) (h1 : 0 ≤ x1) (h2 : 0 ≤ x2)
    (ht : t ∈ Ioo (0:ℝ) 1) (hc : |c| ≤ c2) :
    x0 / t + x1 / (1 - t) + |c| * x2 ≤
      |x0 / t - x1 / (1 - t) + c * x2| + (4 + 2 * c2) * (x0 + x1 + x2) := by
  have ht0 := ht.1
  have ht1 : 0 < 1 - t := by linarith [ht.2]
  have hA : 0 ≤ x0 / t := div_nonneg h0 ht0.le
  have hB : 0 ≤ x1 / (1 - t) := div_nonneg h1 ht1.le
  have hC1 : |c * x2| ≤ c2 * x2 := by rw [abs_mul, abs_of_nonneg h2]; exact mul_le_mul_of_nonneg_right hc h2
  have hC2 : |c| * x2 = |c * x2| := by rw [abs_mul, abs_of_nonneg h2]
  have hc2 : 0 ≤ c2 := (abs_nonneg _).trans hc
  have l1 := le_abs_self (x0 / t - x1 / (1 - t) + c * x2)
  have l2 := neg_abs_le (x0 / t - x1 / (1 - t) + c * x2)
  have l3 := le_abs_self (c * x2)
  have l4 := neg_abs_le (c * x2)
  rw [hC2]
  rcases le_total t (1 / 2) with h | h
  · have : x1 / (1 - t) ≤ 2 * x1 := by rw [div_le_iff₀ ht1]; nlinarith
    nlinarith
  · have : x0 / t ≤ 2 * x0 := by rw [div_le_iff₀ ht0]; nlinarith
    nlinarith

lemma Q_lower (F S E0 E1 Y0 Y1 : ℝ) :
    (E0 * Y1 - E1 * Y0) ^ 2 * (F ^ 2 + S ^ 2) ≤
      ((F * E0 - S * Y0) ^ 2 + (F * E1 - S * Y1) ^ 2) * ((E0 ^ 2 + E1 ^ 2) + (Y0 ^ 2 + Y1 ^ 2)) := by
  set v0 := F * E0 - S * Y0
  set v1 := F * E1 - S * Y1
  have e1 : (S * (E0 * Y1 - E1 * Y0)) ^ 2 = (E0 ^ 2 + E1 ^ 2) * (v0 ^ 2 + v1 ^ 2) -
      (E0 * v0 + E1 * v1) ^ 2 := by simp only [v0, v1]; ring
  have e2 : (F * (E0 * Y1 - E1 * Y0)) ^ 2 = (Y0 ^ 2 + Y1 ^ 2) * (v0 ^ 2 + v1 ^ 2) -
      (Y0 * v0 + Y1 * v1) ^ 2 := by simp only [v0, v1]; ring
  nlinarith [sq_nonneg (E0 * v0 + E1 * v1), sq_nonneg (Y0 * v0 + Y1 * v1)]

lemma cK3_mono {b b' K2 K2' : ℝ} (hb : |b| ≤ b') (hK : 0 ≤ K2) (hK' : K2 ≤ K2') :
    cK3 b K2 ≤ cK3 b' K2' := by
  unfold cK3
  have : K2 ^ 3 ≤ K2' ^ 3 := pow_le_pow_left₀ hK hK' 3
  have : |b| ≤ |b'| := hb.trans (le_abs_self _)
  nlinarith

lemma rpow_upper {x U p : ℝ} (hx : 0 < x) (hxU : x ≤ U) (hU : 1 ≤ U) (hp0 : 0 ≤ p)
    (hp2 : p ≤ 2) : x ^ p ≤ U ^ 2 := by
  calc x ^ p ≤ U ^ p := Real.rpow_le_rpow hx.le hxU hp0
    _ ≤ U ^ (2:ℝ) := Real.rpow_le_rpow_of_exponent_le hU hp2
    _ = U ^ 2 := by norm_num

lemma rpow_lower {x m p : ℝ} (hm : 0 < m) (hmx : m ≤ x) (hp0 : 0 ≤ p) (hp2 : p ≤ 2) :
    (min 1 m) ^ 2 ≤ x ^ p := by
  have h1 : 0 < min 1 m := lt_min one_pos hm
  calc (min 1 m) ^ 2 = (min 1 m) ^ (2:ℝ) := by norm_num
    _ ≤ (min 1 m) ^ p := Real.rpow_le_rpow_of_exponent_ge h1 (min_le_left _ _) hp2
    _ ≤ x ^ p := Real.rpow_le_rpow h1.le ((min_le_right _ _).trans hmx) hp0

lemma abs_log_le_of {x L U : ℝ} (hL : 0 < L) (hLx : L ≤ x) (hxU : x ≤ U) :
    |Real.log x| ≤ |Real.log L| + |Real.log U| := by
  have h1 := Real.log_le_log hL hLx
  have h2 := Real.log_le_log (hL.trans_le hLx) hxU
  rw [abs_le]; constructor
  · linarith [neg_abs_le (Real.log L), abs_nonneg (Real.log U)]
  · linarith [le_abs_self (Real.log U), abs_nonneg (Real.log L)]

lemma tpow_le {t p : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) (hp : 3 / 4 ≤ p) :
    t ^ (p - 1) ≤ (t * (1 - t)) ^ (-(1/4:ℝ)) := by
  have h0 := ht.1
  have h1 : 0 < 1 - t := by linarith [ht.2]
  calc t ^ (p - 1) ≤ t ^ (-(1/4:ℝ)) :=
        Real.rpow_le_rpow_of_exponent_ge h0 ht.2.le (by linarith)
    _ ≤ (t * (1 - t)) ^ (-(1/4:ℝ)) :=
        Real.rpow_le_rpow_of_nonpos (mul_pos h0 h1) (by nlinarith [ht.2]) (by norm_num)

lemma one_le_tt {t c : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) (hc : c ≤ 0) : 1 ≤ (t * (1 - t)) ^ c := by
  have h0 := ht.1
  have h1 : 0 < 1 - t := by linarith [ht.2]
  exact Real.one_le_rpow_of_pos_of_le_one_of_nonpos (mul_pos h0 h1) (by nlinarith [ht.2]) hc

lemma cube_le {x : ℝ} (hx : 0 ≤ x) : (1 + x) ^ 3 ≤ 4 * (1 + x ^ 3) := by
  nlinarith [sq_nonneg (x - 1), mul_nonneg hx (sq_nonneg (x - 1))]

lemma dot_abs_le (E0 E1 Y0 Y1 : ℝ) :
    |E0 * Y0 + E1 * Y1| ≤ ((E0 ^ 2 + E1 ^ 2) + (Y0 ^ 2 + Y1 ^ 2)) / 2 := by
  rw [abs_le]; constructor <;> nlinarith [sq_nonneg (E0 - Y0), sq_nonneg (E1 - Y1),
    sq_nonneg (E0 + Y0), sq_nonneg (E1 + Y1)]

lemma sq_sum_le (S c F : ℝ) (hc : 1 ≤ c) :
    (|S| + c * F) ^ 2 ≤ 2 * c ^ 2 * (F ^ 2 + S ^ 2) := by
  have h1 : |S| ^ 2 = S ^ 2 := sq_abs S
  have h2 : 1 ≤ c ^ 2 := one_le_pow₀ hc
  nlinarith [sq_nonneg (|S| - c * F), sq_nonneg S]

lemma cK3_one (b : ℝ) {K2 : ℝ} (hK : 0 ≤ K2) : 1 ≤ cK3 b K2 := by
  unfold cK3; have := pow_nonneg hK 3; nlinarith [abs_nonneg b]

lemma gen_bound (w : Fin 3 → ℝ) (hw : ∀ j, 0 < w j) {a b m2 M2 c2 d2 CY CE : ℝ}
    (ha : 0 < a) (hb : 0 < b) (hm2 : 0 < m2) (hc2 : 0 ≤ c2) (hd2 : 0 < d2) (hCY : 0 ≤ CY)
    (hCE : 0 ≤ CE) :
    ∃ K, ∀ (t p : ℝ) (α β : Fin 3 → ℝ) (E0 E1 Y0 Y1 : ℝ), t ∈ Ioo (0:ℝ) 1 →
      p ∈ Icc (3/4:ℝ) (5/4) → α 0 = a * t → α 1 = b * (1 - t) → m2 ≤ α 2 → α 2 ≤ M2 →
      β 0 = a → β 1 = -b → |β 2| ≤ c2 * α 2 → d2 ≤ (E0 * Y1 - E1 * Y0) ^ 2 →
      Y0 ^ 2 + Y1 ^ 2 ≤ CY → E0 ^ 2 + E1 ^ 2 ≤ CE →
      ∀ n ≤ 3, nd n (cG w α β E0 E1 Y0 Y1) p ≤ K * (t * (1 - t)) ^ (-(3/4:ℝ)) := by
  have hw0 := hw 0; have hw1 := hw 1; have hw2 := hw 2
  obtain ⟨U, hU⟩ : ∃ U, U = 1 + a + b + |M2| := ⟨_, rfl⟩
  have hU1 : 1 ≤ U := by have := abs_nonneg M2; linarith
  obtain ⟨L0, hL0d⟩ : ∃ L0, L0 = w 2 * (min 1 m2) ^ 2 := ⟨_, rfl⟩
  have hL0 : 0 < L0 := by rw [hL0d]; have : 0 < min 1 m2 := lt_min one_pos hm2; positivity
  obtain ⟨U0, hU0⟩ : ∃ U0, U0 = (w 0 + w 1 + w 2) * U ^ 2 := ⟨_, rfl⟩
  have hU00 : 0 ≤ U0 := by rw [hU0]; positivity
  obtain ⟨CF, hCF⟩ : ∃ CF, CF = |Real.log L0| + |Real.log U0| := ⟨_, rfl⟩
  have hCF0 : 0 ≤ CF := by rw [hCF]; positivity
  obtain ⟨Kd, hKdd⟩ : ∃ Kd, Kd = 1 + |Real.log a| + |Real.log b| + |Real.log m2| +
    |Real.log M2| := ⟨_, rfl⟩
  have hKd : 1 ≤ Kd := by
    rw [hKdd]; have := abs_nonneg (Real.log a); have := abs_nonneg (Real.log b)
    have := abs_nonneg (Real.log m2); have := abs_nonneg (Real.log M2); linarith
  obtain ⟨cc, hccd⟩ : ∃ cc, cc = 4 + 2 * c2 := ⟨_, rfl⟩
  have hcc : 1 ≤ cc + 1 := by rw [hccd]; linarith
  obtain ⟨K0, hK0d⟩ : ∃ K0, K0 = 2 * (cc + 1) ^ 2 * (CE + CY) / d2 := ⟨_, rfl⟩
  have hK0 : 0 ≤ K0 := by rw [hK0d]; positivity
  obtain ⟨CM, hCMd⟩ : ∃ CM, CM = 2 * (CE + CY) := ⟨_, rfl⟩
  have hCM : 0 ≤ CM := by rw [hCMd]; positivity
  obtain ⟨K3, hK3d⟩ : ∃ K3, K3 = cK3 CF (2 * (CM * K0 + 1)) := ⟨_, rfl⟩
  have hK3 : 1 ≤ K3 := by rw [hK3d]; exact cK3_one _ (by positivity)
  obtain ⟨CS, hCSd⟩ : ∃ CS, CS = (w 0 + w 1 + w 2 * c2) * U ^ 2 := ⟨_, rfl⟩
  have hCS : 0 ≤ CS := by rw [hCSd]; positivity
  obtain ⟨CG, hCGd⟩ : ∃ CG, CG = Real.exp (3 * CF) * Real.sqrt CM * (U0 + CS) := ⟨_, rfl⟩
  have hCG : 0 ≤ CG := by rw [hCGd]; positivity
  refine ⟨6 * CG * (K3 * Kd) ^ 3 * (4 * 4097), ?_⟩
  intro t p α β E0 E1 Y0 Y1 ht hp hα0 hα1 hm hM hβ0 hβ1 hβ2 hdd hY hE n hn
  have hp0 : 0 < p := by linarith [hp.1]
  have ht0 := ht.1
  have ht1 : 0 < 1 - t := by linarith [ht.2]
  have hpos : ∀ j, 0 < α j := by
    intro j; fin_cases j
    · simp only [Fin.zero_eta]; rw [hα0]; positivity
    · simp only [Fin.mk_one]; rw [hα1]; positivity
    · exact hm2.trans_le hm
  have hαU : ∀ j, α j ≤ U := by
    have := abs_nonneg M2
    intro j; fin_cases j
    · simp only [Fin.zero_eta]; rw [hα0]
      calc a * t ≤ a * 1 := mul_le_mul_of_nonneg_left ht.2.le ha.le
        _ ≤ U := by rw [hU]; linarith
    · simp only [Fin.mk_one]; rw [hα1]
      calc b * (1 - t) ≤ b * 1 := mul_le_mul_of_nonneg_left (by linarith) hb.le
        _ ≤ U := by rw [hU]; linarith
    · exact hM.trans ((le_abs_self M2).trans (by linarith))
  have hex : ∀ j, Real.exp (Real.log (α j) * p) = α j ^ p := fun j =>
    (Real.rpow_def_of_pos (hpos j) p).symm
  obtain ⟨F, hFdef⟩ : ∃ F, F = cF w α p := ⟨_, rfl⟩
  have hF3 : F = w 0 * α 0 ^ p + w 1 * α 1 ^ p + w 2 * α 2 ^ p := by
    simp only [hFdef, cF, Fin.sum_univ_three, hex]
  have hpow : ∀ j, α j ^ p ≤ U ^ 2 := fun j =>
    rpow_upper (hpos j) (hαU j) hU1 hp0.le (by linarith [hp.2])
  have hpow0 : ∀ j, 0 ≤ α j ^ p := fun j => (Real.rpow_pos_of_pos (hpos j) p).le
  have hFU : F ≤ U0 := by
    rw [hF3, hU0]
    have := mul_le_mul_of_nonneg_left (hpow 0) hw0.le
    have := mul_le_mul_of_nonneg_left (hpow 1) hw1.le
    have := mul_le_mul_of_nonneg_left (hpow 2) hw2.le
    have : (w 0 + w 1 + w 2) * U ^ 2 = w 0 * U ^ 2 + w 1 * U ^ 2 + w 2 * U ^ 2 := by ring
    linarith
  have hFL : L0 ≤ F := by
    rw [hF3, hL0d]
    have := mul_le_mul_of_nonneg_left (rpow_lower hm2 hm hp0.le (by linarith [hp.2])) hw2.le
    have := mul_nonneg hw0.le (hpow0 0); have := mul_nonneg hw1.le (hpow0 1)
    linarith
  have hlogF : |Real.log F| ≤ CF := hCF ▸ abs_log_le_of hL0 hFL hFU
  -- `D`
  obtain ⟨D, hDd⟩ : ∃ D, D = Kd * (1 + ls t) := ⟨_, rfl⟩
  have hls := ls_nonneg t
  have hD1 : 1 ≤ D := by rw [hDd]; exact one_le_mul_of_one_le_of_one_le hKd (by linarith)
  have hKdls : Kd + ls t ≤ D := by
    rw [hDd]; have := mul_le_mul_of_nonneg_right hKd hls; linarith
  have hD : ∀ j, |Real.log (α j)| ≤ D := by
    have := abs_nonneg (Real.log a); have := abs_nonneg (Real.log b)
    have := abs_nonneg (Real.log m2); have := abs_nonneg (Real.log M2)
    intro j; fin_cases j
    · simp only [Fin.zero_eta]; rw [hα0, Real.log_mul ha.ne' ht0.ne']
      refine (abs_add_le _ _).trans ?_
      have : |Real.log t| ≤ ls t := abs_log_le_ls t
      rw [hKdd] at hKdls; linarith
    · simp only [Fin.mk_one]; rw [hα1, Real.log_mul hb.ne' ht1.ne']
      refine (abs_add_le _ _).trans ?_
      have : |Real.log (1 - t)| ≤ ls t := abs_log1_le_ls t
      rw [hKdd] at hKdls; linarith
    · show |Real.log (α 2)| ≤ D
      have := abs_log_le_of hm2 hm hM
      rw [hKdd] at hKdls; linarith
  -- `S`, `St`
  obtain ⟨S, hSd⟩ : ∃ S, S = cS w α β p := ⟨_, rfl⟩
  obtain ⟨St, hStd⟩ : ∃ St, St = cSt w α β p := ⟨_, rfl⟩
  have hx0 : 0 ≤ w 0 * α 0 ^ p := mul_nonneg hw0.le (hpow0 0)
  have hx1 : 0 ≤ w 1 * α 1 ^ p := mul_nonneg hw1.le (hpow0 1)
  have hx2 : 0 ≤ w 2 * α 2 ^ p := mul_nonneg hw2.le (hpow0 2)
  have hS3 : S = (w 0 * α 0 ^ p) / t - (w 1 * α 1 ^ p) / (1 - t) +
      (β 2 / α 2) * (w 2 * α 2 ^ p) := by
    rw [hSd]; simp only [cS, Fin.sum_univ_three, hex, hβ0, hβ1]
    rw [hα0, hα1]
    field_simp
    ring
  have hSt3 : St = (w 0 * α 0 ^ p) / t + (w 1 * α 1 ^ p) / (1 - t) +
      |β 2 / α 2| * (w 2 * α 2 ^ p) := by
    rw [hStd]; simp only [cSt, Fin.sum_univ_three, hex, hβ0, hβ1]
    rw [hα0, hα1]
    have q0 : w 0 * a / (a * t) = w 0 / t := by field_simp
    have q1 : w 1 * -b / (b * (1 - t)) = -(w 1 / (1 - t)) := by field_simp
    have q2 : |w 2 * β 2 / α 2| = w 2 * |β 2 / α 2| := by rw [mul_div_assoc, abs_mul, abs_of_pos hw2]
    rw [q0, q1, q2, abs_neg, abs_of_pos (div_pos hw0 ht0), abs_of_pos (div_pos hw1 ht1)]; ring
  have hq2 : |β 2 / α 2| ≤ c2 := by
    rw [abs_div, abs_of_pos (hpos 2), div_le_iff₀ (hpos 2)]; exact hβ2
  have hStle : St ≤ |S| + cc * F := by
    rw [hSt3, hS3, hF3, hccd]
    exact St_le hx0 hx1 hx2 ht hq2
  have hSt0 : 0 ≤ St := hStd ▸ cSt_nonneg p
  have hFpos : 0 < F := hL0.trans_le hFL
  obtain ⟨Q, hQd⟩ : ∃ Q, Q = cQ w α β E0 E1 Y0 Y1 p := ⟨_, rfl⟩
  have hdd0 : E0 * Y1 - E1 * Y0 ≠ 0 := by
    intro h; rw [h] at hdd; norm_num at hdd; linarith
  have hQpos : 0 < Q := hQd ▸ cQ_pos hw hdd0 p
  have hQlow := Q_lower F S E0 E1 Y0 Y1
  have hQe : Q = (F * E0 - S * Y0) ^ 2 + (F * E1 - S * Y1) ^ 2 := by rw [hQd, hFdef, hSd]; rfl
  rw [← hQe] at hQlow
  have hFS : F ^ 2 + S ^ 2 ≤ Q * (CE + CY) / d2 := by
    rw [le_div_iff₀ hd2]
    have h1 : d2 * (F ^ 2 + S ^ 2) ≤ (E0 * Y1 - E1 * Y0) ^ 2 * (F ^ 2 + S ^ 2) :=
      mul_le_mul_of_nonneg_right hdd (by positivity)
    have h2 : Q * ((E0 ^ 2 + E1 ^ 2) + (Y0 ^ 2 + Y1 ^ 2)) ≤ Q * (CE + CY) :=
      mul_le_mul_of_nonneg_left (by linarith) hQpos.le
    linarith
  have hQ : (F + St) ^ 2 ≤ K0 * Q := by
    have h1 : F + St ≤ |S| + (cc + 1) * F := by linarith
    have h2 : (F + St) ^ 2 ≤ (|S| + (cc + 1) * F) ^ 2 :=
      pow_le_pow_left₀ (by linarith) h1 2
    have h3 : (|S| + (cc + 1) * F) ^ 2 ≤ 2 * (cc + 1) ^ 2 * (F ^ 2 + S ^ 2) :=
      sq_sum_le S (cc + 1) F hcc
    have h4 : 2 * (cc + 1) ^ 2 * (F ^ 2 + S ^ 2) ≤ 2 * (cc + 1) ^ 2 * (Q * (CE + CY) / d2) :=
      mul_le_mul_of_nonneg_left hFS (by positivity)
    have : 2 * (cc + 1) ^ 2 * (Q * (CE + CY) / d2) = K0 * Q := by rw [hK0d]; ring
    linarith
  have hcM : cM E0 E1 Y0 Y1 ≤ CM := by
    unfold cM
    rw [abs_of_nonneg (by positivity : (0:ℝ) ≤ E0 ^ 2 + E1 ^ 2),
      abs_of_nonneg (by positivity : (0:ℝ) ≤ Y0 ^ 2 + Y1 ^ 2)]
    have := dot_abs_le E0 E1 Y0 Y1
    rw [hCMd]; linarith
  have hcM0 : 0 ≤ cM E0 E1 Y0 Y1 := cM_nonneg
  -- the derivative bound
  rw [hFdef, hStd, hQd] at hQ
  have hnd := nd_cG hw hdd0 hD1 hD (by linarith [hp.1]) hK0 hQ hn
  rw [← hFdef] at hnd
  have hk3 : cK3 (Real.log F) (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) ≤ K3 := by
    rw [hK3d]; refine cK3_mono hlogF (by positivity) ?_
    have := mul_le_mul_of_nonneg_right hcM hK0; linarith
  have hk3' : 1 ≤ cK3 (Real.log F) (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) := cK3_one _ (by positivity)
  have hpowK : (cK3 (Real.log F) (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) * D) ^ n ≤ (K3 * D) ^ 3 := by
    calc (cK3 (Real.log F) (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) * D) ^ n ≤ (K3 * D) ^ n :=
          pow_le_pow_left₀ (by positivity) (mul_le_mul_of_nonneg_right hk3 (by linarith)) n
      _ ≤ (K3 * D) ^ 3 := pow_le_pow_right₀ (one_le_mul_of_one_le_of_one_le hK3 hD1) hn
  -- the size of `G`
  set T := (t * (1 - t)) ^ (-(1/4:ℝ)) with hT
  have hT1 : 1 ≤ T := one_le_tt ht (by norm_num)
  have hStT : St ≤ CS * T := by
    rw [hSt3, hα0, hα1, Real.mul_rpow ha.le ht0.le, Real.mul_rpow hb.le ht1.le, hCSd]
    have e0 : w 0 * (a ^ p * t ^ p) / t = w 0 * a ^ p * t ^ (p - 1) := by
      rw [Real.rpow_sub_one ht0.ne']; ring
    have e1 : w 1 * (b ^ p * (1 - t) ^ p) / (1 - t) = w 1 * b ^ p * (1 - t) ^ (p - 1) := by
      rw [Real.rpow_sub_one ht1.ne']; ring
    rw [e0, e1]
    have ha' : a ^ p ≤ U ^ 2 := rpow_upper ha (by rw [hU]; have := abs_nonneg M2; linarith) hU1
      hp0.le (by linarith [hp.2])
    have hb' : b ^ p ≤ U ^ 2 := rpow_upper hb (by rw [hU]; have := abs_nonneg M2; linarith) hU1
      hp0.le (by linarith [hp.2])
    have t0 : t ^ (p - 1) ≤ T := tpow_le ht hp.1
    have t1 : (1 - t) ^ (p - 1) ≤ T := by
      have := tpow_le (t := 1 - t) ⟨ht1, by linarith⟩ hp.1
      rwa [show (1 - t) * (1 - (1 - t)) = t * (1 - t) by ring] at this
    have e2 : |β 2 / α 2| * (w 2 * α 2 ^ p) ≤ c2 * (w 2 * U ^ 2) * T := by
      have := hpow 2
      calc |β 2 / α 2| * (w 2 * α 2 ^ p) ≤ c2 * (w 2 * U ^ 2) :=
            mul_le_mul hq2 (mul_le_mul_of_nonneg_left this hw2.le) (by positivity) hc2
        _ ≤ c2 * (w 2 * U ^ 2) * T := le_mul_of_one_le_right (by positivity) hT1
    have := Real.rpow_pos_of_pos ha p
    have := Real.rpow_pos_of_pos hb p
    have := Real.rpow_pos_of_pos ht0 (p - 1)
    have := Real.rpow_pos_of_pos ht1 (p - 1)
    have f0 : w 0 * a ^ p * t ^ (p - 1) ≤ w 0 * U ^ 2 * T := by gcongr
    have f1 : w 1 * b ^ p * (1 - t) ^ (p - 1) ≤ w 1 * U ^ 2 * T := by gcongr
    have : (w 0 + w 1 + w 2 * c2) * U ^ 2 * T =
      w 0 * U ^ 2 * T + w 1 * U ^ 2 * T + c2 * (w 2 * U ^ 2) * T := by ring
    linarith
  have hQup : Q ≤ CM * (F + St) ^ 2 := by
    have := nd_cQ (α := α) (β := β) (E0 := E0) (E1 := E1) (Y0 := Y0) (Y1 := Y1) hw hD hp0
      (n := 0) (by norm_num)
    rw [nd_zero, pow_zero, pow_zero, one_mul, mul_one, ← hQd, ← hFdef, ← hStd,
      abs_of_pos hQpos] at this
    refine this.trans (mul_le_mul_of_nonneg_right ?_ (sq_nonneg _))
    exact hcM
  have hGle : cG w α β E0 E1 Y0 Y1 p ≤ CG * T := by
    unfold cG ch
    rw [← hFdef, ← hQd, Real.exp_add]
    have hinv : p⁻¹ ≤ 2 := by rw [inv_le_comm₀ hp0 (by norm_num)]; linarith [hp.1]
    have h1 : Real.exp (-(1 + p⁻¹) * Real.log F) ≤ Real.exp (3 * CF) := by
      apply Real.exp_le_exp.2
      have hpi : 0 < p⁻¹ := inv_pos.2 hp0
      have h1 := mul_le_mul_of_nonneg_left (neg_abs_le (Real.log F)) (by linarith : (0:ℝ) ≤ 1 + p⁻¹)
      have h2 := mul_le_mul (by linarith : 1 + p⁻¹ ≤ 3) hlogF (abs_nonneg _) (by norm_num)
      have : -(1 + p⁻¹) * Real.log F = -((1 + p⁻¹) * Real.log F) := by ring
      linarith
    have h2 : Real.exp (1 / 2 * Real.log Q) = Real.sqrt Q := by
      rw [Real.sqrt_eq_rpow, Real.rpow_def_of_pos hQpos, mul_comm]
    rw [h2]
    have h3 : Real.sqrt Q ≤ Real.sqrt CM * (F + St) := by
      rw [← Real.sqrt_sq (by linarith : 0 ≤ F + St), ← Real.sqrt_mul hCM]
      exact Real.sqrt_le_sqrt hQup
    have h4 : F + St ≤ (U0 + CS) * T := by
      have := le_mul_of_one_le_right hU00 hT1
      have : (U0 + CS) * T = U0 * T + CS * T := by ring
      linarith
    rw [hCGd]
    calc Real.exp (-(1 + p⁻¹) * Real.log F) * Real.sqrt Q
        ≤ Real.exp (3 * CF) * (Real.sqrt CM * ((U0 + CS) * T)) :=
          mul_le_mul h1 (h3.trans (mul_le_mul_of_nonneg_left h4 (Real.sqrt_nonneg _)))
            (Real.sqrt_nonneg _) (exp_pos _).le
      _ = _ := by ring
  -- the log factor
  have hls3 : (1 + ls t) ^ 3 ≤ 4 * 4097 * (t * (1 - t)) ^ (-(1/2:ℝ)) := by
    have h1 := ls_pow_le ht 3
    norm_num at h1
    have h2 : 1 ≤ (t * (1 - t)) ^ (-(1/2:ℝ)) := one_le_tt ht (by norm_num)
    have h3 : (1 + ls t) ^ 3 ≤ 4 * (1 + ls t ^ 3) := cube_le hls
    linarith
  have hTT : T * (t * (1 - t)) ^ (-(1/2:ℝ)) = (t * (1 - t)) ^ (-(3/4:ℝ)) := by
    rw [hT, ← Real.rpow_add (mul_pos ht0 ht1)]; norm_num
  have hG0 : 0 ≤ cG w α β E0 E1 Y0 Y1 p := (exp_pos _).le
  calc nd n (cG w α β E0 E1 Y0 Y1) p
      ≤ 6 * cG w α β E0 E1 Y0 Y1 p *
          (cK3 (Real.log F) (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) * D) ^ n := hnd
    _ ≤ 6 * (CG * T) * ((K3 * Kd) ^ 3 * (1 + ls t) ^ 3) := by
        rw [show (K3 * Kd) ^ 3 * (1 + ls t) ^ 3 = (K3 * D) ^ 3 by rw [hDd]; ring]
        gcongr
    _ ≤ 6 * (CG * T) * ((K3 * Kd) ^ 3 * (4 * 4097 * (t * (1 - t)) ^ (-(1/2:ℝ)))) := by
        gcongr
    _ = 6 * CG * (K3 * Kd) ^ 3 * (4 * 4097) * (T * (t * (1 - t)) ^ (-(1/2:ℝ))) := by ring
    _ = _ := by rw [hTT]

lemma min_le_aff {r t : ℝ} (ht0 : 0 ≤ t) (ht1 : t ≤ 1) : min 1 r ≤ 1 + (r - 1) * t := by
  rcases le_total 1 r with h | h
  · rw [min_eq_left h]; nlinarith
  · rw [min_eq_right h]; nlinarith

lemma aff_le_max {r t : ℝ} (ht0 : 0 ≤ t) (ht1 : t ≤ 1) : 1 + (r - 1) * t ≤ max 1 r := by
  rcases le_total 1 r with h | h
  · rw [max_eq_right h]; nlinarith
  · rw [max_eq_left h]; nlinarith

lemma edge_bound (H : Hex) (hH : IsHex H) (k : Fin 3) :
    ∃ K, ∀ p ∈ Icc (3/4:ℝ) (5/4), ∀ t ∈ Ioo (0:ℝ) 1, ∀ n ≤ 3,
      nd n (fun q => Gf H k q t) p ≤ K * (t * (1 - t)) ^ (-(3/4:ℝ)) := by
  obtain ⟨-, -, -, hB, hB1, hA0, hA1, hr⟩ := S05.proof_8 H hH k
  have hw := hH.1
  have hmr : 0 < min 1 (rr H k) := lt_min one_pos hr
  set V0 := eV0 H k
  set E := eE H k
  have hdet := edge_det_pos H hH k
  obtain ⟨K, hK⟩ := gen_bound (ew H k) (ew_pos H hH k) (a := A0 H k / H.w k)
    (b := A1 H k / H.w (k + 1)) (m2 := Bv H k * min 1 (rr H k) / H.w (k + 2))
    (M2 := Bv H k * max 1 (rr H k) / H.w (k + 2)) (c2 := |rr H k - 1| / min 1 (rr H k))
    (d2 := (det2 (eV0 H k) (eV1 H k)) ^ 2) (CY := 2 * (V0 0 ^ 2 + V0 1 ^ 2) + 2 * (E 0 ^ 2 + E 1 ^ 2))
    (CE := E 0 ^ 2 + E 1 ^ 2) (div_pos hA0 (hw _)) (div_pos hA1 (hw _))
    (div_pos (mul_pos hB hmr) (hw _)) (by positivity) (by positivity) (by positivity)
    (by positivity)
  refine ⟨K, fun p hp t ht n hn => ?_⟩
  have ht0 := ht.1.le
  have ht1 := ht.2.le
  have hml := min_le_aff (r := rr H k) ht0 ht1
  have hw2 := hw (k + 2)
  refine hK t p (eα H k t) (eβ H k) (E 0) (E 1) (ept H k t 0) (ept H k t 1) ht hp
    ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_ n hn
  · simp [eα]; ring
  · simp [eα]; ring
  · simp only [eα, Matrix.cons_val_two, Matrix.tail_cons, Matrix.head_cons]
    gcongr
  · simp only [eα, Matrix.cons_val_two, Matrix.tail_cons, Matrix.head_cons]
    gcongr
    exact aff_le_max ht0 ht1
  · simp [eα, eβ]
  · simp [eα, eβ]
  · simp only [eα, eβ, Matrix.cons_val_two, Matrix.tail_cons, Matrix.head_cons]
    have : 0 ≤ |rr H k - 1| := abs_nonneg _
    calc |Bv H k * (rr H k - 1) / H.w (k + 2)|
        = |rr H k - 1| * Bv H k * min 1 (rr H k) / (min 1 (rr H k) * H.w (k + 2)) := by
          rw [abs_div, abs_mul, abs_of_pos hB, abs_of_pos hw2]; field_simp
      _ ≤ |rr H k - 1| * Bv H k * (1 + (rr H k - 1) * t) / (min 1 (rr H k) * H.w (k + 2)) := by
          gcongr
      _ = _ := by field_simp
  · have e : E 0 * ept H k t 1 - E 1 * ept H k t 0 = -det2 (eV0 H k) (eV1 H k) := by
      simp only [E, eE, ept, det2, PiLp.add_apply, PiLp.sub_apply, PiLp.smul_apply, smul_eq_mul]
      ring
    rw [e, neg_sq]
  · have e0 : ept H k t 0 = V0 0 + t * E 0 := by
      simp [ept, V0, E, eE]
    have e1 : ept H k t 1 = V0 1 + t * E 1 := by
      simp [ept, V0, E, eE]
    rw [e0, e1]
    have : t ^ 2 ≤ 1 := by nlinarith
    nlinarith [sq_nonneg (V0 0 - t * E 0), sq_nonneg (V0 1 - t * E 1), sq_nonneg (E 0),
      sq_nonneg (E 1), mul_le_mul_of_nonneg_right this (sq_nonneg (E 0)),
      mul_le_mul_of_nonneg_right this (sq_nonneg (E 1))]
  · exact le_rfl

/-! ### differentiating under the integral sign -/

lemma chain3 {f0 f1 f2 f3 : ℝ → ℝ} {U : Set ℝ} (hU : IsOpen U) {x : ℝ} (hx : x ∈ U)
    (h0 : ∀ y ∈ U, HasDerivAt f0 (f1 y) y) (h1 : ∀ y ∈ U, HasDerivAt f1 (f2 y) y)
    (h2 : ∀ y ∈ U, HasDerivAt f2 (f3 y) y) (h3 : ContinuousOn f3 U) :
    ContDiffAt ℝ 3 f0 x ∧ iteratedDeriv 3 f0 x = f3 x := by
  have d0 : EqOn (deriv f0) f1 U := fun y hy => (h0 y hy).deriv
  have d1 : EqOn (deriv f1) f2 U := fun y hy => (h1 y hy).deriv
  have d2 : EqOn (deriv f2) f3 U := fun y hy => (h2 y hy).deriv
  have c2 : ContDiffOn ℝ 1 f2 U := by
    rw [show (1 : WithTop ℕ∞) = 0 + 1 from rfl, contDiffOn_succ_iff_deriv_of_isOpen hU]
    refine ⟨fun y hy => (h2 y hy).differentiableAt.differentiableWithinAt, by simp, ?_⟩
    exact (contDiffOn_zero.2 h3).congr d2
  have c1 : ContDiffOn ℝ 2 f1 U := by
    rw [show (2 : WithTop ℕ∞) = 1 + 1 from rfl, contDiffOn_succ_iff_deriv_of_isOpen hU]
    exact ⟨fun y hy => (h1 y hy).differentiableAt.differentiableWithinAt, by simp, c2.congr d1⟩
  have c0 : ContDiffOn ℝ 3 f0 U := by
    rw [show (3 : WithTop ℕ∞) = 2 + 1 from rfl, contDiffOn_succ_iff_deriv_of_isOpen hU]
    exact ⟨fun y hy => (h0 y hy).differentiableAt.differentiableWithinAt, by simp, c1.congr d0⟩
  refine ⟨c0.contDiffAt (hU.mem_nhds hx), ?_⟩
  have e1 : iteratedDeriv 1 f0 =ᶠ[𝓝 x] f1 := by
    filter_upwards [hU.mem_nhds hx] with y hy
    rw [iteratedDeriv_one]; exact d0 hy
  have e2 : iteratedDeriv 2 f0 =ᶠ[𝓝 x] f2 := by
    filter_upwards [hU.mem_nhds hx] with y hy
    have : iteratedDeriv 1 f0 =ᶠ[𝓝 y] f1 := by
      filter_upwards [hU.mem_nhds hy] with z hz
      rw [iteratedDeriv_one]; exact d0 hz
    rw [iteratedDeriv_succ, this.deriv_eq]; exact d1 hy
  rw [iteratedDeriv_succ, e2.deriv_eq]; exact d2 hx

/-- the `j`-th `p`-derivative of the length element (within `p > 0`). -/
noncomputable def Gw (H : Hex) (k : Fin 3) (j : ℕ) (p t : ℝ) : ℝ :=
  iteratedDerivWithin j (fun q => Gf H k q t) (Ioi 0) p

lemma Gf_cd (H : Hex) (hH : IsHex H) (k : Fin 3) (t : ℝ) :
    ContDiffOn ℝ 3 (fun q => Gf H k q t) (Ioi 0) := by
  have := ch_cd (α := eα H k t) (β := eβ H k) (ew_pos H hH k) (edge_det_ne H hH k t)
  exact this.exp

lemma Gw_hasDeriv (H : Hex) (hH : IsHex H) (k : Fin 3) {j : ℕ} (hj : j < 3) {p : ℝ}
    (hp : 0 < p) (t : ℝ) : HasDerivAt (fun q => Gw H k j q t) (Gw H k (j + 1) p t) p := by
  have hd := (Gf_cd H hH k t).differentiableOn_iteratedDerivWithin (m := j)
    (by exact_mod_cast hj) uS
  have hda := (hd p hp).differentiableAt (isOpen_Ioi.mem_nhds hp)
  have := hda.hasDerivAt
  unfold Gw
  rw [iteratedDerivWithin_succ, derivWithin_of_isOpen isOpen_Ioi hp]
  exact this

lemma Gw_contOn (H : Hex) (hH : IsHex H) (k : Fin 3) (t : ℝ) :
    ContinuousOn (fun q => Gw H k 3 q t) (Ioi 0) :=
  (Gf_cd H hH k t).continuousOn_iteratedDerivWithin le_rfl uS

lemma Gw_zero (H : Hex) (k : Fin 3) (p t : ℝ) : Gw H k 0 p t = Gf H k p t := by
  simp [Gw]

lemma Gw_bound (H : Hex) (hH : IsHex H) (k : Fin 3) :
    ∃ K, ∀ p ∈ Icc (3/4:ℝ) (5/4), ∀ t ∈ Ioo (0:ℝ) 1, ∀ j ≤ 3,
      |Gw H k j p t| ≤ K * (t * (1 - t)) ^ (-(3/4:ℝ)) := by
  obtain ⟨K, hK⟩ := edge_bound H hH k
  exact ⟨K, fun p hp t ht j hj => hK p hp t ht j hj⟩

lemma Gl_meas (H : Hex) (k : Fin 3) (p : ℝ) : Measurable (fun t => Gl H k p t) :=
  (measurable_deriv (arc H k p)).norm

lemma Gw_meas (H : Hex) (hH : IsHex H) (k : Fin 3) :
    ∀ j ≤ 3, ∀ q > 0, AEStronglyMeasurable (fun t => Gw H k j q t)
      (volume.restrict (Ι (0:ℝ) 1)) := by
  intro j
  induction j with
  | zero =>
    intro _ q hq
    refine (Gl_meas H k q).aestronglyMeasurable.congr ?_
    filter_upwards [ae_Ioo01] with t ht
    rw [Gw_zero]; exact (arc_deriv H hH k hq ht).2
  | succ j ih =>
    intro hj q hq
    have hq' : ∀ m : ℕ, 0 < q + 1 / ((m : ℝ) + 1) := fun m => by positivity
    refine aestronglyMeasurable_of_tendsto_ae atTop
      (f := fun (m : ℕ) t => ((m : ℝ) + 1) * (Gw H k j (q + 1 / ((m : ℝ) + 1)) t - Gw H k j q t))
      (fun m => ((ih (by omega) _ (hq' m)).sub (ih (by omega) q hq)).const_mul _) ?_
    refine Filter.Eventually.of_forall fun t => ?_
    have hD := Gw_hasDeriv H hH k (j := j) (by omega) hq t
    have hs := hD.tendsto_slope_zero_right
    have hseq : Tendsto (fun m : ℕ => 1 / ((m : ℝ) + 1)) atTop (𝓝[>] 0) := by
      refine tendsto_nhdsWithin_iff.2 ⟨tendsto_one_div_add_atTop_nhds_zero_nat, ?_⟩
      exact Filter.Eventually.of_forall fun m => mem_Ioi.2 (by positivity)
    refine (hs.comp hseq).congr fun m => ?_
    simp only [Function.comp_apply, smul_eq_mul, one_div, inv_inv]

/-- `∫_0^1 ∂_p^j G(p,t) dt`. -/
noncomputable def Phi (H : Hex) (k : Fin 3) (j : ℕ) (p : ℝ) : ℝ := ∫ t in (0:ℝ)..1, Gw H k j p t

lemma bnd_int (K : ℝ) :
    IntervalIntegrable (fun t : ℝ => K * (t * (1 - t)) ^ (-(3/4:ℝ))) volume 0 1 :=
  (ii_tt _ (by norm_num)).const_mul K

lemma ae_io {P : ℝ → Prop} (h : ∀ t ∈ Ioo (0:ℝ) 1, P t) : ∀ᵐ t ∂volume, t ∈ Ι (0:ℝ) 1 → P t :=
  (ae_restrict_iff' measurableSet_uIoc).1 (by filter_upwards [ae_Ioo01] with t ht; exact h t ht)

lemma Gw_int (H : Hex) (hH : IsHex H) (k : Fin 3) {j : ℕ} (hj : j ≤ 3) {p : ℝ}
    (hp : p ∈ Icc (3/4:ℝ) (5/4)) : IntervalIntegrable (fun t => Gw H k j p t) volume 0 1 := by
  obtain ⟨K, hK⟩ := Gw_bound H hH k
  refine (bnd_int K).mono_fun' (Gw_meas H hH k j hj p (by linarith [hp.1])) ?_
  filter_upwards [ae_Ioo01] with t ht
  rw [Real.norm_eq_abs]; exact hK p hp t ht j hj

lemma Phi_deriv (H : Hex) (hH : IsHex H) (k : Fin 3) {j : ℕ} (hj : j < 3) {p0 : ℝ}
    (hp0 : p0 ∈ Ioo (3/4:ℝ) (5/4)) : HasDerivAt (Phi H k j) (Phi H k (j + 1) p0) p0 := by
  obtain ⟨K, hK⟩ := Gw_bound H hH k
  have := intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le (μ := volume)
    (a := 0) (b := 1) (F := fun p t => Gw H k j p t) (F' := fun p t => Gw H k (j + 1) p t)
    (x₀ := p0) (s := Ioo (3/4) (5/4)) (bound := fun t => K * (t * (1 - t)) ^ (-(3/4:ℝ)))
    (isOpen_Ioo.mem_nhds hp0)
    (by filter_upwards [isOpen_Ioo.mem_nhds hp0] with q hq
        exact Gw_meas H hH k j hj.le q (by linarith [hq.1]))
    (Gw_int H hH k hj.le (Ioo_subset_Icc_self hp0))
    (Gw_meas H hH k (j + 1) hj p0 (by linarith [hp0.1]))
    (ae_io fun t ht x hx => by
      rw [Real.norm_eq_abs]; exact hK x (Ioo_subset_Icc_self hx) t ht (j + 1) hj)
    (bnd_int K)
    (ae_io fun t _ x hx => Gw_hasDeriv H hH k hj (by linarith [hx.1]) t)
  exact this.2

lemma Phi_cont (H : Hex) (hH : IsHex H) (k : Fin 3) :
    ContinuousOn (Phi H k 3) (Ioo (3/4:ℝ) (5/4)) := by
  obtain ⟨K, hK⟩ := Gw_bound H hH k
  intro p0 hp0
  refine (intervalIntegral.continuousAt_of_dominated_interval (μ := volume)
    (F := fun p t => Gw H k 3 p t) (bound := fun t => K * (t * (1 - t)) ^ (-(3/4:ℝ)))
    ?_ ?_ (bnd_int K) ?_).continuousWithinAt
  · filter_upwards [isOpen_Ioo.mem_nhds hp0] with q hq
    exact Gw_meas H hH k 3 le_rfl q (by linarith [hq.1])
  · filter_upwards [isOpen_Ioo.mem_nhds hp0] with q hq
    exact ae_io fun t ht => by
      rw [Real.norm_eq_abs]; exact hK q (Ioo_subset_Icc_self hq) t ht 3 le_rfl
  · exact ae_io fun t _ => ((Gw_contOn H hH k t) p0 (by simp; linarith [hp0.1])).continuousAt
      (isOpen_Ioi.mem_nhds (by simp; linarith [hp0.1]))

/-! ### one edge: length = integral of the length element -/

lemma eα_cont (H : Hex) (k : Fin 3) (j : Fin 3) : Continuous (fun s => eα H k s j) := by
  fin_cases j <;> simp [eα] <;> fun_prop

lemma ept_cont (H : Hex) (k : Fin 3) : Continuous (ept H k) := by
  unfold ept; fun_prop

lemma ept_ne (H : Hex) (hH : IsHex H) (k : Fin 3) (t : ℝ) (ht : t ∈ Icc (0:ℝ) 1) :
    ept H k t ≠ 0 := by
  intro h0
  have hd := edge_det_pos H hH k
  have e : det2 (eV0 H k) (ept H k t) = t * det2 (eV0 H k) (eV1 H k) := by
    simp only [det2, ept, PiLp.add_apply, PiLp.smul_apply, PiLp.sub_apply, smul_eq_mul]; ring
  rw [h0] at e
  simp only [det2, PiLp.zero_apply, mul_zero, sub_zero] at e
  rcases mul_eq_zero.1 e.symm with h | h
  · have : ept H k 0 = eV0 H k := by simp [ept]
    rw [h] at h0; rw [h0] at this
    rw [← this] at hd; simp [det2] at hd
  · simp only [det2] at hd; linarith

lemma cexp {f : ℝ → ℝ} {t : ℝ} (h : ContinuousAt f t) :
    ContinuousAt (fun s => Real.exp (f s)) t := Real.continuous_exp.continuousAt.comp h

lemma Gf_contAt (H : Hex) (hH : IsHex H) (k : Fin 3) (p : ℝ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    ContinuousAt (fun s => Gf H k p s) t := by
  have hpos := eα_pos H hH k ht
  have hterm : ∀ (c : Fin 3 → ℝ) (j : Fin 3),
      ContinuousAt (fun s => c j * Real.exp (Real.log (eα H k s j) * p)) t := fun c j =>
    continuousAt_const.mul (cexp (((eα_cont H k j).continuousAt.log (hpos j).ne').mul
      continuousAt_const))
  have hF : ContinuousAt (fun s => cF (ew H k) (eα H k s) p) t := by
    simp only [cF, Fin.sum_univ_three]
    exact ((hterm _ 0).add (hterm _ 1)).add (hterm _ 2)
  have hterm' : ∀ j : Fin 3, ContinuousAt (fun s => (ew H k j * eβ H k j / eα H k s j) *
      Real.exp (Real.log (eα H k s j) * p)) t := fun j =>
    (continuousAt_const.div (eα_cont H k j).continuousAt (hpos j).ne').mul
      (cexp (((eα_cont H k j).continuousAt.log (hpos j).ne').mul continuousAt_const))
  have hS : ContinuousAt (fun s => cS (ew H k) (eα H k s) (eβ H k) p) t := by
    simp only [cS, Fin.sum_univ_three]
    exact ((hterm' 0).add (hterm' 1)).add (hterm' 2)
  have hy0 : ContinuousAt (fun s => ept H k s 0) t :=
    ((PiLp.continuous_apply 2 _ 0).comp (ept_cont H k)).continuousAt
  have hy1 : ContinuousAt (fun s => ept H k s 1) t :=
    ((PiLp.continuous_apply 2 _ 1).comp (ept_cont H k)).continuousAt
  have hQ : ContinuousAt (fun s => cQ (ew H k) (eα H k s) (eβ H k) (eE H k 0) (eE H k 1)
      (ept H k s 0) (ept H k s 1) p) t := by
    unfold cQ
    exact ((hF.mul continuousAt_const).sub (hS.mul hy0)).pow 2 |>.add
      (((hF.mul continuousAt_const).sub (hS.mul hy1)).pow 2)
  unfold Gf cG ch
  exact cexp ((continuousAt_const.mul (hF.log (cF_pos (ew_pos H hH k) p).ne')).add
    (continuousAt_const.mul (hQ.log (cQ_pos (ew_pos H hH k) (edge_det_ne H hH k t) p).ne')))

lemma Gl_eq_Gw (H : Hex) (hH : IsHex H) (k : Fin 3) {p : ℝ} (hp : 0 < p) :
    ∀ᵐ t ∂(volume.restrict (Ι (0:ℝ) 1)), Gl H k p t = Gw H k 0 p t := by
  filter_upwards [ae_Ioo01] with t ht
  rw [Gw_zero]; exact (arc_deriv H hH k hp ht).2

lemma speedData (H : Hex) (hH : IsHex H) (k : Fin 3) {p : ℝ} (hp : p ∈ Icc (3/4:ℝ) (5/4)) :
    X04.SpeedData (arc H k p) (Gl H k p) := by
  have hp0 : 0 < p := by linarith [hp.1]
  refine ⟨?_, ?_, ?_, ?_, ?_, Gl_meas H k p⟩
  · intro t ht
    have hF : 0 < Fw H.u H.w p (ept H k t) := Fw_pos H hH (ept_ne H hH k t ht)
    have hc : ContinuousAt (fun s => Fw H.u H.w p (ept H k s)) t :=
      ((Fw_cont H.u H.w hp0).comp (ept_cont H k)).continuousAt
    exact ((hc.rpow_const (Or.inl hF.ne')).smul (ept_cont H k).continuousAt).continuousWithinAt
  · intro t ht
    exact ⟨(arc_deriv H hH k hp0 ht).1, rfl⟩
  · exact (Gw_int H hH k (j := 0) (by norm_num) hp).congr_ae
      (by filter_upwards [Gl_eq_Gw H hH k hp0] with t ht; exact ht.symm)
  · intro t; exact norm_nonneg _
  · intro t ht
    have hev : (fun s => Gf H k p s) =ᶠ[𝓝 t] Gl H k p := by
      filter_upwards [isOpen_Ioo.mem_nhds ht] with s hs
      exact ((arc_deriv H hH k hp0 hs).2).symm
    exact (Gf_contAt H hH k p ht).congr hev

lemma arc_var (H : Hex) (hH : IsHex H) (k : Fin 3) {p : ℝ} (hp : p ∈ Icc (3/4:ℝ) (5/4)) :
    eVariationOn (arc H k p) (Icc 0 1) = ENNReal.ofReal (Phi H k 0 p) := by
  rw [(speedData H hH k hp).var_eq, Phi]
  congr 1
  exact intervalIntegral.integral_congr_ae (ae_io fun t ht => (Gw_zero H k p t).symm ▸
    (arc_deriv H hH k (by linarith [hp.1]) ht).2)

/-! ### the polar curve, sector by sector -/

lemma ray0 (θ : ℝ) : ray θ 0 = Real.cos θ := by simp [ray]
lemma ray1 (θ : ℝ) : ray θ 1 = Real.sin θ := by simp [ray]

lemma plane_ext {v w : X04.Plane} (h0 : v 0 = w 0) (h1 : v 1 = w 1) : v = w := by
  ext i; fin_cases i <;> simp [h0, h1]

lemma ray_add_pi (θ : ℝ) : ray (θ + π) = -ray θ := by
  apply plane_ext <;> simp [ray0, ray1, Real.cos_add_pi, Real.sin_add_pi]

lemma ray_add_two_pi (θ : ℝ) : ray (θ + 2 * π) = ray θ := by
  apply plane_ext <;> simp [ray0, ray1]

lemma det_ray (a b : ℝ) : det2 (ray a) (ray b) = Real.sin (b - a) := by
  simp only [det2, ray0, ray1, Real.sin_sub]; ring

lemma Fw_negp (u : Fin 3 → X04.Plane) (w : Fin 3 → ℝ) (p : ℝ) (x : X04.Plane) :
    Fw u w p (-x) = Fw u w p x := by
  rw [show -x = (-1 : ℝ) • x by simp, Fw_smul]; simp

lemma polar_add_pi (H : Hex) (p θ : ℝ) : polar H p (θ + π) = -polar H p θ := by
  simp only [polar, ray_add_pi, Fw_negp, smul_neg]

lemma polar_add_two_pi (H : Hex) (p θ : ℝ) : polar H p (θ + 2 * π) = polar H p θ := by
  simp only [polar, ray_add_two_pi]

lemma var_trans {E : Type*} [PseudoEMetricSpace E] (g : ℝ → E) (a b c : ℝ) :
    eVariationOn (fun s => g (s + c)) (Icc a b) = eVariationOn g (Icc (a + c) (b + c)) := by
  rw [show (fun s => g (s + c)) = g ∘ fun t => t + c from rfl,
    eVariationOn.comp_eq_of_monotoneOn g (fun t => t + c) (fun x _ y _ h => by simpa using h)]
  congr 1
  ext x
  simp only [mem_image, mem_Icc]
  constructor
  · rintro ⟨y, ⟨h1, h2⟩, rfl⟩; constructor <;> linarith
  · rintro ⟨h1, h2⟩; exact ⟨x - c, ⟨by linarith, by linarith⟩, by ring⟩

lemma var_neg (f : ℝ → X04.Plane) (s : Set ℝ) :
    eVariationOn (fun x => -f x) s = eVariationOn f s := by
  unfold eVariationOn; simp only [edist_neg_neg]

lemma var_split {E : Type*} [PseudoEMetricSpace E] (f : ℝ → E) {a b c : ℝ} (hab : a ≤ b)
    (hbc : b ≤ c) : eVariationOn f (Icc a c) = eVariationOn f (Icc a b) + eVariationOn f (Icc b c) := by
  have := eVariationOn.Icc_add_Icc f (s := univ) hab hbc (mem_univ _)
  simp only [univ_inter] at this
  exact this.symm

lemma var_periodic {E : Type*} [PseudoEMetricSpace E] (f : ℝ → E) {T a : ℝ}
    (hf : ∀ x, f (x + T) = f x) (ha0 : 0 ≤ a) (haT : a ≤ T) :
    eVariationOn f (Icc 0 T) = eVariationOn f (Icc a (a + T)) := by
  rw [var_split f ha0 haT, var_split f haT (by linarith : T ≤ a + T)]
  have : eVariationOn f (Icc T (a + T)) = eVariationOn f (Icc 0 a) := by
    have h := var_trans f 0 a T
    simp only [hf, zero_add] at h
    rw [h]
  rw [this, add_comm]

lemma sin_le_zero_of {d : ℝ} (h1 : π ≤ d) (h2 : d ≤ 2 * π) : Real.sin d ≤ 0 := by
  have := Real.sin_nonneg_of_nonneg_of_le_pi (x := d - π) (by linarith) (by linarith)
  rw [Real.sin_sub_pi] at this; linarith

lemma lift_angle {θ0 θ : ℝ} (h0 : 0 ≤ θ0) (h0' : θ0 < 2 * π) (h1 : 0 ≤ θ) (h1' : θ < 2 * π)
    (hs : 0 < Real.sin (θ - θ0)) : ∃ φ, ray φ = ray θ ∧ θ0 < φ ∧ φ < θ0 + π := by
  rcases lt_or_ge θ0 θ with h | h
  · refine ⟨θ, rfl, h, ?_⟩
    by_contra hc
    have := sin_le_zero_of (d := θ - θ0) (by linarith) (by linarith)
    linarith
  · refine ⟨θ + 2 * π, ray_add_two_pi θ, by linarith, ?_⟩
    by_contra hc
    have := sin_le_zero_of (d := θ + 2 * π - θ0) (by linarith) (by linarith)
    rw [show θ + 2 * π - θ0 = (θ - θ0) + 2 * π by ring, Real.sin_add_two_pi] at this
    linarith

lemma det_smul_ray (r s a b : ℝ) : det2 (r • ray a) (s • ray b) = r * s * Real.sin (b - a) := by
  rw [← det_ray]; simp only [det2, PiLp.smul_apply, smul_eq_mul]; ring

lemma exists_polar_pos {x : X04.Plane} (hx : x ≠ 0) :
    ∃ r θ : ℝ, 0 < r ∧ 0 ≤ θ ∧ θ < 2 * π ∧ x = r • ray θ := by
  obtain ⟨r, θ, hr, h0, h1, hxr⟩ := polar_coords x
  refine ⟨r, θ, lt_of_le_of_ne hr ?_, h0, h1, hxr⟩
  rintro rfl; apply hx; rw [hxr, zero_smul]

lemma ne_of_det {x y : X04.Plane} (h : 0 < det2 x y) : x ≠ 0 ∧ y ≠ 0 := by
  constructor <;> rintro rfl <;> simp [det2] at h

lemma corner_angles (H : Hex) (hH : IsHex H) :
    ∃ φ0 φ1 φ2 r0 r1 r2 : ℝ, 0 ≤ φ0 ∧ φ0 < 2 * π ∧ 0 < r0 ∧ 0 < r1 ∧ 0 < r2 ∧
      H.C 0 = r0 • ray φ0 ∧ H.C 1 = r1 • ray φ1 ∧ H.C 2 = r2 • ray φ2 ∧
      φ0 < φ1 ∧ φ1 < φ2 ∧ φ2 < φ0 + π := by
  obtain ⟨-, -, -, -, h01, h12, h02⟩ := hH
  obtain ⟨r0, θ0, hr0, a0, b0, e0⟩ := exists_polar_pos (ne_of_det h01).1
  obtain ⟨r1, θ1, hr1, a1, b1, e1⟩ := exists_polar_pos (ne_of_det h01).2
  obtain ⟨r2, θ2, hr2, a2, b2, e2⟩ := exists_polar_pos (ne_of_det h12).2
  rw [e0, e1, det_smul_ray] at h01
  rw [e0, e2, det_smul_ray] at h02
  have s1 : 0 < Real.sin (θ1 - θ0) := pos_of_mul_pos_right h01 (by positivity)
  have s2 : 0 < Real.sin (θ2 - θ0) := pos_of_mul_pos_right h02 (by positivity)
  obtain ⟨φ1, q1, c1, d1⟩ := lift_angle a0 b0 a1 b1 s1
  obtain ⟨φ2, q2, c2, d2⟩ := lift_angle a0 b0 a2 b2 s2
  have e1' : H.C 1 = r1 • ray φ1 := by rw [e1, q1]
  have e2' : H.C 2 = r2 • ray φ2 := by rw [e2, q2]
  refine ⟨θ0, φ1, φ2, r0, r1, r2, a0, b0, hr0, hr1, hr2, e0, e1', e2', c1, ?_, d2⟩
  have h12' : 0 < r1 * r2 * Real.sin (φ2 - φ1) := by rw [← det_smul_ray, ← e1', ← e2']; exact h12
  have s3 : 0 < Real.sin (φ2 - φ1) := pos_of_mul_pos_right h12' (by positivity)
  by_contra hc
  have := Real.sin_nonpos_of_nonpos_of_neg_pi_le (x := φ2 - φ1) (by linarith) (by linarith)
  linarith

lemma sector_var (H : Hex) (hH : IsHex H) (k : Fin 3) {p : ℝ} (hp : 0 < p) {a b r0 r1 : ℝ}
    (hr0 : 0 < r0) (hr1 : 0 < r1) (h0 : eV0 H k = r0 • ray a) (h1 : eV1 H k = r1 • ray b)
    (hab : a < b) (hab' : b < a + π) :
    eVariationOn (polar H p) (Icc a b) = eVariationOn (arc H k p) (Icc 0 1) := by
  have hD := edge_det_pos H hH k
  set V0 := eV0 H k
  set V1 := eV1 H k
  set D := det2 V0 V1
  set A : ℝ → ℝ := fun θ => det2 V0 (ray θ) with hAd
  set B : ℝ → ℝ := fun θ => det2 (ray θ) V1 with hBd
  have hA : ∀ θ, A θ = r0 * Real.sin (θ - a) := fun θ => by
    rw [← det_ray]; simp only [hAd, h0, det2, PiLp.smul_apply, smul_eq_mul]; ring
  have hB : ∀ θ, B θ = r1 * Real.sin (b - θ) := fun θ => by
    rw [← det_ray]; simp only [hBd, h1, det2, PiLp.smul_apply, smul_eq_mul]; ring
  have hA0 : ∀ θ ∈ Icc a b, 0 ≤ A θ := fun θ hθ => by
    rw [hA]; exact mul_nonneg hr0.le (Real.sin_nonneg_of_nonneg_of_le_pi (by linarith [hθ.1])
      (by linarith [hθ.2]))
  have hB0 : ∀ θ ∈ Icc a b, 0 ≤ B θ := fun θ hθ => by
    rw [hB]; exact mul_nonneg hr1.le (Real.sin_nonneg_of_nonneg_of_le_pi (by linarith [hθ.2])
      (by linarith [hθ.1]))
  have hAB : ∀ θ ∈ Icc a b, 0 < A θ + B θ := by
    intro θ hθ
    rcases lt_or_eq_of_le hθ.2 with h | h
    · have : 0 < B θ := by
        rw [hB]; exact mul_pos hr1 (Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith [hθ.1]))
      linarith [hA0 θ hθ]
    · have : 0 < A θ := by
        rw [hA, h]; exact mul_pos hr0 (Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith))
      linarith [hB0 θ hθ]
  set τ : ℝ → ℝ := fun θ => A θ / (A θ + B θ) with hτ
  have key : ∀ θ, B θ * V0 0 + A θ * V1 0 = D * ray θ 0 ∧
      B θ * V0 1 + A θ * V1 1 = D * ray θ 1 := fun θ => by
    simp only [hAd, hBd, D, det2]; constructor <;> ring
  have hept : ∀ θ ∈ Icc a b, ept H k (τ θ) = (D / (A θ + B θ)) • ray θ := by
    intro θ hθ
    have hne := (hAB θ hθ).ne'
    obtain ⟨k0, k1⟩ := key θ
    have g : ∀ x y : ℝ, x + A θ / (A θ + B θ) * (y - x) = (B θ * x + A θ * y) / (A θ + B θ) :=
      fun x y => by field_simp; ring
    apply plane_ext
    · simp only [ept, hτ, PiLp.add_apply, PiLp.smul_apply, PiLp.sub_apply, smul_eq_mul]
      rw [g, k0]; ring
    · simp only [ept, hτ, PiLp.add_apply, PiLp.smul_apply, PiLp.sub_apply, smul_eq_mul]
      rw [g, k1]; ring
  have harc : ∀ θ ∈ Icc a b, arc H k p (τ θ) = polar H p θ := by
    intro θ hθ
    have hc : 0 < D / (A θ + B θ) := div_pos hD (hAB θ hθ)
    have hF := Phi_pos H hH (p := p) θ
    rw [arc, hept θ hθ, Fw_smul, abs_of_pos hc, Real.mul_rpow (by positivity) hF.le,
      ← Real.rpow_mul hc.le, show p * -(1 / p) = -1 by field_simp, Real.rpow_neg_one, smul_smul,
      polar]
    congr 1
    have := (hAB θ hθ).ne'
    have := hD.ne'
    field_simp
  have hmono : MonotoneOn τ (Icc a b) := by
    intro θ hθ θ' hθ' hle
    simp only [hτ]
    rw [div_le_div_iff₀ (hAB θ hθ) (hAB θ' hθ')]
    have idn : A θ' * B θ - A θ * B θ' = D * det2 (ray θ) (ray θ') := by
      simp only [hAd, hBd, D, det2]; ring
    rw [det_ray] at idn
    have : 0 ≤ Real.sin (θ' - θ) :=
      Real.sin_nonneg_of_nonneg_of_le_pi (by linarith) (by linarith [hθ.1, hθ'.2])
    nlinarith [mul_nonneg hD.le this]
  have hAc : Continuous A := by
    simp only [hAd, det2]; exact (continuous_const.mul ((PiLp.continuous_apply 2 _ 1).comp
      continuous_ray)).sub (continuous_const.mul ((PiLp.continuous_apply 2 _ 0).comp continuous_ray))
  have hBc : Continuous B := by
    simp only [hBd, det2]; exact (((PiLp.continuous_apply 2 _ 0).comp continuous_ray).mul
      continuous_const).sub (((PiLp.continuous_apply 2 _ 1).comp continuous_ray).mul continuous_const)
  have hcont : ContinuousOn τ (Icc a b) :=
    hAc.continuousOn.div (hAc.add hBc).continuousOn fun θ hθ => (hAB θ hθ).ne'
  have ha : a ∈ Icc a b := ⟨le_rfl, hab.le⟩
  have hb : b ∈ Icc a b := ⟨hab.le, le_rfl⟩
  have hτa : τ a = 0 := by simp only [hτ, hA a, sub_self, Real.sin_zero, mul_zero, zero_div]
  have hτb : τ b = 1 := by
    have h1 : B b = 0 := by rw [hB, sub_self, Real.sin_zero, mul_zero]
    have h2 : 0 < A b := by linarith [hAB b hb]
    simp only [hτ, h1, add_zero, div_self h2.ne']
  have himg : τ '' Icc a b = Icc 0 1 := by
    apply Subset.antisymm
    · rintro _ ⟨θ, hθ, rfl⟩
      exact ⟨hτa ▸ hmono ha hθ hθ.1, hτb ▸ hmono hθ hb hθ.2⟩
    · have := intermediate_value_Icc hab.le hcont
      rwa [hτa, hτb] at this
  calc eVariationOn (polar H p) (Icc a b) = eVariationOn (arc H k p ∘ τ) (Icc a b) :=
        eVariationOn.eq_of_eqOn fun θ hθ => (harc θ hθ).symm
    _ = eVariationOn (arc H k p) (τ '' Icc a b) := eVariationOn.comp_eq_of_monotoneOn _ _ hmono
    _ = _ := by rw [himg]

lemma Phi_nonneg (H : Hex) (k : Fin 3) (p : ℝ) : 0 ≤ Phi H k 0 p :=
  intervalIntegral.integral_nonneg zero_le_one fun t _ => by
    rw [Gw_zero]; exact (exp_pos _).le

lemma var_shift_pi (H : Hex) (p a b : ℝ) :
    eVariationOn (polar H p) (Icc (a + π) (b + π)) = eVariationOn (polar H p) (Icc a b) := by
  rw [← var_trans]
  simp only [polar_add_pi]
  exact var_neg _ _

lemma perim_eq (H : Hex) (hH : IsHex H) {p : ℝ} (hp : p ∈ Icc (3/4:ℝ) (5/4)) :
    Perim H p = 2 * ∑ k, Phi H k 0 p := by
  have hp0 : 0 < p := by linarith [hp.1]
  obtain ⟨φ0, φ1, φ2, r0, r1, r2, h0a, h0b, hr0, hr1, hr2, e0, e1, e2, c01, c12, c20⟩ :=
    corner_angles H hH
  set P := polar H p
  have v0 : eVariationOn P (Icc φ0 φ1) = ENNReal.ofReal (Phi H 0 0 p) := by
    rw [sector_var H hH 0 hp0 hr0 hr1 (by simp [eV0, e0]) (by simp [eV1, e1]) c01 (by linarith),
      arc_var H hH 0 hp]
  have v1 : eVariationOn P (Icc φ1 φ2) = ENNReal.ofReal (Phi H 1 0 p) := by
    rw [sector_var H hH 1 hp0 hr1 hr2 (by simp [eV0, e1]) (by simp [eV1, e2]) c12 (by linarith),
      arc_var H hH 1 hp]
  have v2 : eVariationOn P (Icc φ2 (φ0 + π)) = ENNReal.ofReal (Phi H 2 0 p) := by
    rw [sector_var H hH 2 hp0 hr2 hr0 (by simp [eV0, e2])
      (by simp [eV1, e0, ray_add_pi]) (by linarith) (by linarith), arc_var H hH 2 hp]
  have hper : eVariationOn P (Icc 0 (2 * π)) = eVariationOn P (Icc φ0 (φ0 + 2 * π)) :=
    var_periodic P (fun x => polar_add_two_pi H p x) h0a h0b.le
  have hπ := Real.pi_pos
  rw [Perim, hper, var_split P c01.le (by linarith : φ1 ≤ φ0 + 2 * π),
    var_split P c12.le (by linarith : φ2 ≤ φ0 + 2 * π),
    var_split P c20.le (by linarith : φ0 + π ≤ φ0 + 2 * π),
    var_split P (by linarith : φ0 + π ≤ φ1 + π) (by linarith : φ1 + π ≤ φ0 + 2 * π),
    var_split P (by linarith : φ1 + π ≤ φ2 + π) (by linarith : φ2 + π ≤ φ0 + 2 * π),
    show φ0 + 2 * π = (φ0 + π) + π by ring, var_shift_pi, var_shift_pi, var_shift_pi, v0, v1, v2]
  have n0 := Phi_nonneg H 0 p
  have n1 := Phi_nonneg H 1 p
  have n2 := Phi_nonneg H 2 p
  rw [Fin.sum_univ_three]
  simp only [ENNReal.toReal_add, ENNReal.ofReal_ne_top, ne_eq, not_false_eq_true,
    ENNReal.add_eq_top, or_self, ENNReal.toReal_ofReal n0, ENNReal.toReal_ofReal n1,
    ENNReal.toReal_ofReal n2]
  ring

lemma iter_Gl (H : Hex) (hH : IsHex H) (k : Fin 3) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    iteratedDeriv 3 (fun p => Gl H k p t) 1 = Gw H k 3 1 t := by
  have hev : (fun p => Gl H k p t) =ᶠ[𝓝 1] (fun p => Gf H k p t) := by
    filter_upwards [isOpen_Ioi.mem_nhds (show (1:ℝ) ∈ Ioi 0 by norm_num)] with p hp
    exact (arc_deriv H hH k hp ht).2
  rw [hev.iteratedDeriv_eq 3, Gw, iteratedDerivWithin_of_isOpen isOpen_Ioi (by norm_num : (1:ℝ) ∈ Ioi 0)]

lemma Ee_eq (H : Hex) (hH : IsHex H) (k : Fin 3) : Ee H k = Phi H k 3 1 :=
  intervalIntegral.integral_congr_ae (ae_io fun t ht => iter_Gl H hH k ht)

lemma one_mem_U : (1:ℝ) ∈ Ioo (3/4:ℝ) (5/4) := ⟨by norm_num, by norm_num⟩

end S05

open Set in
theorem S05.proof_9 : S05.stmt_9 := by
  intro H hH
  have hd : ∀ j < 3, ∀ y ∈ Ioo (3/4:ℝ) (5/4), HasDerivAt (fun p => 2 * ∑ k, S05.Phi H k j p)
      (2 * ∑ k, S05.Phi H k (j + 1) y) y := fun j hj y hy =>
    (HasDerivAt.fun_sum fun k _ => S05.Phi_deriv H hH k hj hy).const_mul 2
  have hc : ContinuousOn (fun p => 2 * ∑ k, S05.Phi H k 3 p) (Ioo (3/4:ℝ) (5/4)) :=
    continuousOn_const.mul (continuousOn_finsetSum (f := fun k p => S05.Phi H k 3 p) Finset.univ
      fun k _ => S05.Phi_cont H hH k)
  obtain ⟨hcd, hit⟩ := S05.chain3 isOpen_Ioo S05.one_mem_U (hd 0 (by norm_num))
    (hd 1 (by norm_num)) (hd 2 (by norm_num)) hc
  have hev : S05.Perim H =ᶠ[nhds 1] fun p => 2 * ∑ k, S05.Phi H k 0 p := by
    filter_upwards [isOpen_Ioo.mem_nhds S05.one_mem_U] with p hp
    exact S05.perim_eq H hH (Ioo_subset_Icc_self hp)
  refine ⟨hcd.congr_of_eventuallyEq hev, fun k => ?_, ?_⟩
  · exact (S05.Gw_int H hH k (j := 3) le_rfl (Ioo_subset_Icc_self S05.one_mem_U)).congr_ae
      (by filter_upwards [S05.ae_Ioo01] with t ht; exact (S05.iter_Gl H hH k ht).symm)
  · rw [hev.iteratedDeriv_eq 3, hit]
    simp only [S05.Ee_eq H hH]

/-! ## stmt_13: the `r = 1`, `w = 1` case -/

namespace S05
open MeasureTheory Set Filter Real Interval Topology

lemma tele1 : HasSum (fun q : ℕ => 1 / ((q:ℝ) + 1) - 1 / ((q:ℝ) + 2)) 1 := by
  have hf : ∀ q : ℕ, 0 ≤ 1 / ((q:ℝ) + 1) - 1 / ((q:ℝ) + 2) := fun q => by
    rw [sub_nonneg]; gcongr; linarith
  rw [hasSum_iff_tendsto_nat_of_nonneg hf]
  have e : ∀ n : ℕ, ∑ i ∈ Finset.range n, (1 / ((i:ℝ) + 1) - 1 / ((i:ℝ) + 2)) = 1 - 1 / ((n:ℝ) + 1) := by
    intro n
    have := Finset.sum_range_sub' (fun i : ℕ => 1 / ((i:ℝ) + 1)) n
    simp only [Nat.cast_add, Nat.cast_one, CharP.cast_eq_zero, zero_add, div_one] at this
    rw [← this]; refine Finset.sum_congr rfl fun i _ => ?_; ring
  simp_rw [e]
  have := (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)).const_sub 1
  simpa using this

lemma tele2 : HasSum (fun q : ℕ => 1 / ((q:ℝ) + 1) - 1 / ((q:ℝ) + 3)) (3 / 2) := by
  have h2 := (hasSum_nat_add_iff' 1).2 tele1
  simp only [Finset.range_one, Finset.sum_singleton] at h2
  have := tele1.add h2
  convert this using 1
  · funext q; push_cast; ring
  · norm_num

lemma zeta2_s1 : HasSum (fun q : ℕ => 1 / ((q:ℝ) + 2) ^ 2) (π ^ 2 / 6 - 1) := by
  have h0 : HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 2) (π ^ 2 / 6) := by
    have := D00.hasSum_Li2 (x := 1) (by norm_num)
    simpa using this
  have := (hasSum_nat_add_iff' 1).2 h0
  convert this using 1
  · funext q; push_cast; ring_nf
  · simp

lemma zeta2_s2 : HasSum (fun q : ℕ => 1 / ((q:ℝ) + 3) ^ 2) (π ^ 2 / 6 - 5 / 4) := by
  have h0 : HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 2) (π ^ 2 / 6) := by
    have := D00.hasSum_Li2 (x := 1) (by norm_num)
    simpa using this
  have := (hasSum_nat_add_iff' 2).2 h0
  convert this using 1
  · funext q; push_cast; ring_nf
  · simp [Finset.sum_range_succ]; norm_num

lemma zeta3_s1 : HasSum (fun q : ℕ => 1 / ((q:ℝ) + 2) ^ 3) (Z03.z3 - 1) := by
  have := (hasSum_nat_add_iff' 1).2 Z03.hasSum_z3
  convert this using 1
  · funext q; push_cast; ring_nf
  · simp

lemma zeta3_s2 : HasSum (fun q : ℕ => 1 / ((q:ℝ) + 3) ^ 3) (Z03.z3 - 9 / 8) := by
  have := (hasSum_nat_add_iff' 2).2 Z03.hasSum_z3
  convert this using 1
  · funext q; push_cast; ring_nf
  · simp [Finset.sum_range_succ]; norm_num

lemma sumA : HasSum (fun q : ℕ => 1 / (((q:ℝ) + 1) * ((q:ℝ) + 2) ^ 2)) (2 - π ^ 2 / 6) := by
  have := tele1.sub zeta2_s1
  convert this using 1
  · funext q; field_simp; ring
  · ring

lemma sumB : HasSum (fun q : ℕ => 1 / (((q:ℝ) + 1) * ((q:ℝ) + 2) ^ 3))
    (3 - π ^ 2 / 6 - Z03.z3) := by
  have := (tele1.sub zeta2_s1).sub zeta3_s1
  convert this using 1
  · funext q; field_simp; ring
  · ring

lemma sumC : HasSum (fun q : ℕ => 1 / (((q:ℝ) + 1) * ((q:ℝ) + 3) ^ 3))
    (17 / 16 - π ^ 2 / 24 - Z03.z3 / 2) := by
  have := ((tele2.mul_left (1/8)).sub (zeta2_s2.mul_left (1/4))).sub (zeta3_s2.mul_left (1/2))
  convert this using 1
  · funext q; field_simp; ring
  · ring

end S05

namespace S05
open MeasureTheory Set Filter Real Interval Topology

lemma Iint20 (N : ℕ) : Iint N 2 0 = 2 / ((N : ℝ) + 1) ^ 3 := by
  have hN : (0:ℝ) < (N : ℝ) + 1 := by positivity
  let n : ℝ := (N : ℝ) + 1
  have hc : ContinuousOn (fun t : ℝ => t ^ (N + 1) * Real.log t ^ 2 / n
      - 2 * (t ^ (N + 1) * Real.log t ^ 1) / n ^ 2 + 2 * t ^ (N + 1) / n ^ 3) (Icc 0 1) := by
    have h2 := G_cont N 2 two_ne_zero
    have h1 := G_cont N 1 one_ne_zero
    have a := h2.div_const n
    have b := (h1.const_smul (2:ℝ)).div_const (n ^ 2)
    have c := ((continuous_pow (M := ℝ) (N+1)).continuousOn (s := Icc 0 1)).const_smul (2:ℝ)
      |>.div_const (n ^ 3)
    exact (a.sub b).add c
  have hd : ∀ t ∈ Ioo (0:ℝ) 1, HasDerivAt (fun t : ℝ => t ^ (N + 1) * Real.log t ^ 2 / n
      - 2 * (t ^ (N + 1) * Real.log t ^ 1) / n ^ 2 + 2 * t ^ (N + 1) / n ^ 3)
      (t ^ N * Real.log t ^ 2) t := by
    intro t ht
    have hl := Real.hasDerivAt_log ht.1.ne'
    have h1 := ((((hasDerivAt_pow (N + 1) t).mul (hl.pow 2)).div_const n).sub
      ((((hasDerivAt_pow (N + 1) t).mul (hl.pow 1)).const_mul 2).div_const (n ^ 2))).add
      (((hasDerivAt_pow (N + 1) t).const_mul 2).div_const (n ^ 3))
    convert h1 using 1
    have := ht.1.ne'
    simp only [n, Pi.pow_apply, Nat.add_sub_cancel, Nat.cast_add, Nat.cast_one]
    field_simp
    ring
  have hi : IntervalIntegrable (fun t : ℝ => t ^ N * Real.log t ^ 2) volume 0 1 := by
    simpa using ii_g N 2 0
  have := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one hc hd hi
  unfold Iint
  simp only [pow_zero, mul_one]
  rw [this]
  simp only [n]
  simp <;> field_simp

lemma Iint30 (N : ℕ) : Iint N 3 0 = -6 / ((N : ℝ) + 1) ^ 4 := by
  have hN : (0:ℝ) < (N : ℝ) + 1 := by positivity
  let n : ℝ := (N : ℝ) + 1
  have hc : ContinuousOn (fun t : ℝ => t ^ (N + 1) * Real.log t ^ 3 / n
      - 3 * (t ^ (N + 1) * Real.log t ^ 2) / n ^ 2 + 6 * (t ^ (N + 1) * Real.log t ^ 1) / n ^ 3
      - 6 * t ^ (N + 1) / n ^ 4) (Icc 0 1) := by
    have h3 := G_cont N 3 (by norm_num)
    have h2 := G_cont N 2 two_ne_zero
    have h1 := G_cont N 1 one_ne_zero
    have a := h3.div_const n
    have b := (h2.const_smul (3:ℝ)).div_const (n ^ 2)
    have b' := (h1.const_smul (6:ℝ)).div_const (n ^ 3)
    have c := ((continuous_pow (M := ℝ) (N+1)).continuousOn (s := Icc 0 1)).const_smul (6:ℝ)
      |>.div_const (n ^ 4)
    exact ((a.sub b).add b').sub c
  have hd : ∀ t ∈ Ioo (0:ℝ) 1, HasDerivAt (fun t : ℝ => t ^ (N + 1) * Real.log t ^ 3 / n
      - 3 * (t ^ (N + 1) * Real.log t ^ 2) / n ^ 2 + 6 * (t ^ (N + 1) * Real.log t ^ 1) / n ^ 3
      - 6 * t ^ (N + 1) / n ^ 4)
      (t ^ N * Real.log t ^ 3) t := by
    intro t ht
    have hl := Real.hasDerivAt_log ht.1.ne'
    have hp := hasDerivAt_pow (N + 1) t
    have h1 := ((((hp.mul (hl.pow 3)).div_const n).sub
      (((hp.mul (hl.pow 2)).const_mul 3).div_const (n ^ 2))).add
      (((hp.mul (hl.pow 1)).const_mul 6).div_const (n ^ 3))).sub
      ((hp.const_mul 6).div_const (n ^ 4))
    convert h1 using 1
    have := ht.1.ne'
    simp only [n, Pi.pow_apply, Nat.add_sub_cancel, Nat.cast_add, Nat.cast_one]
    field_simp
    ring
  have hi : IntervalIntegrable (fun t : ℝ => t ^ N * Real.log t ^ 3) volume 0 1 := by
    simpa using ii_g N 3 0
  have := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one hc hd hi
  unfold Iint
  simp only [pow_zero, mul_one]
  rw [this]
  simp only [n]
  simp <;> field_simp

lemma Iint_refl0 (a b : ℕ) : Iint 0 a b = Iint 0 b a := by
  unfold Iint
  have := intervalIntegral.integral_comp_sub_left
    (fun s : ℝ => s ^ 0 * Real.log s ^ a * Real.log (1 - s) ^ b) (a := 0) (b := 1) 1
  simp only [sub_zero, sub_self, sub_sub_cancel] at this
  rw [← this]
  congr 1; funext s; ring

lemma Iint_refl1 (a b : ℕ) : Iint 1 a b = Iint 0 a b - Iint 1 b a := by
  have h := intervalIntegral.integral_comp_sub_left
    (fun s : ℝ => s ^ 1 * Real.log s ^ b * Real.log (1 - s) ^ a) (a := 0) (b := 1) 1
  simp only [sub_zero, sub_self, sub_sub_cancel] at h
  unfold Iint
  rw [eq_sub_iff_add_eq, ← h, ← intervalIntegral.integral_add]
  · congr 1; funext s; ring
  · exact (ii_g 1 a b)
  · have := ((ii_g 1 b a).comp_sub_left 1).symm
    simpa using this

lemma Iint_a1 (N a : ℕ) :
    HasSum (fun q : ℕ => -(Iint (N + q + 1) a 0 / ((q:ℝ) + 1))) (Iint N a 1) := by
  set g := fun t : ℝ => t ^ N * Real.log t ^ a with hg
  have hgi : IntervalIntegrable (fun t => |g t| * -Real.log (1 - t)) volume 0 1 := by
    refine ii_of_ls (by fun_prop : Measurable (fun t => |g t| * -Real.log (1 - t))).aestronglyMeasurable
      1 (a + 1) (fun t ht => ?_)
    have := g_bound (n := N) (a := a) (b := 1) ht
    simp only [pow_one] at this
    rw [abs_mul, abs_neg, abs_abs, ← abs_mul]; simpa [hg] using this
  have key := intervalIntegral.hasSum_integral_of_dominated_convergence
    (μ := volume) (a := 0) (b := 1)
    (F := fun (q : ℕ) (t : ℝ) => g t * (-(t ^ (q + 1) / ((q:ℝ) + 1))))
    (f := fun t => g t * Real.log (1 - t))
    (fun q t => |g t| * (t ^ (q + 1) / ((q:ℝ) + 1)))
    (fun q => (by fun_prop : Measurable (fun t : ℝ =>
      g t * (-(t ^ (q + 1) / ((q:ℝ) + 1))))).aestronglyMeasurable)
    (fun q => by
      refine Filter.Eventually.of_forall (fun t ht => ?_)
      rw [Set.uIoc_of_le zero_le_one] at ht
      rw [Real.norm_eq_abs, abs_mul, abs_neg,
        abs_of_nonneg (div_nonneg (pow_nonneg ht.1.le _) (by positivity))])
    (ae_io fun t ht => by
      have := (Real.hasSum_pow_div_log_of_abs_lt_one (x := t)
        (by rw [abs_of_pos ht.1]; exact ht.2)).summable
      exact this.mul_left _)
    (by
      refine hgi.congr_ae ?_
      filter_upwards [ae_Ioo01] with t ht
      have := (Real.hasSum_pow_div_log_of_abs_lt_one (x := t)
        (by rw [abs_of_pos ht.1]; exact ht.2)).mul_left (|g t|)
      exact this.tsum_eq.symm)
    (ae_io fun t ht => by
      have := (Real.hasSum_pow_div_log_of_abs_lt_one (x := t)
        (by rw [abs_of_pos ht.1]; exact ht.2)).neg.mul_left (g t)
      simpa using this)
  convert key using 1
  · funext q
    unfold Iint
    have : (fun t => g t * -(t ^ (q + 1) / ((q:ℝ) + 1))) = fun t => (-(1 / ((q:ℝ) + 1))) *
        (t ^ (N + q + 1) * Real.log t ^ a * Real.log (1 - t) ^ 0) := by
      funext t; simp only [hg]; rw [show N + q + 1 = N + (q + 1) by ring, pow_add]; ring
    rw [this, intervalIntegral.integral_const_mul]; ring
  · unfold Iint; congr 1; funext t; simp [hg]

end S05

namespace S05
open MeasureTheory Set Filter Real Interval Topology

lemma I011 : Iint 0 1 1 = 2 - π ^ 2 / 6 := by
  refine (Iint_a1 0 1).unique ?_
  convert sumA using 1
  funext q; rw [Iint10]; push_cast; field_simp; ring

lemma I021 : Iint 0 2 1 = -6 + π ^ 2 / 3 + 2 * Z03.z3 := by
  refine (Iint_a1 0 2).unique ?_
  convert sumB.mul_left (-2) using 1
  · funext q; rw [Iint20]; push_cast; field_simp; ring
  · ring

lemma I121 : Iint 1 2 1 = -17 / 8 + π ^ 2 / 12 + Z03.z3 := by
  refine (Iint_a1 1 2).unique ?_
  convert sumC.mul_left (-2) using 1
  · funext q; rw [Iint20]; push_cast; field_simp; ring
  · ring

lemma I020 : Iint 0 2 0 = 2 := by rw [Iint20]; norm_num
lemma I002 : Iint 0 0 2 = 2 := by rw [Iint_refl0, I020]
lemma I030 : Iint 0 3 0 = -6 := by rw [Iint30]; norm_num
lemma I003 : Iint 0 0 3 = -6 := by rw [Iint_refl0, I030]
lemma I012 : Iint 0 1 2 = -6 + π ^ 2 / 3 + 2 * Z03.z3 := by rw [Iint_refl0, I021]
lemma I130 : Iint 1 3 0 = -3 / 8 := by rw [Iint30]; norm_num
lemma I103 : Iint 1 0 3 = -45 / 8 := by rw [Iint_refl1, I003, I130]; norm_num
lemma I112 : Iint 1 1 2 = -31 / 8 + π ^ 2 / 4 + Z03.z3 := by
  rw [Iint_refl1, I012, I121]; ring

end S05

namespace S05
open Set Filter Topology Real

/-- a 3-jet of `f` on `U`: `f1, f2, f3` are the successive derivatives. -/
structure J3 (U : Set ℝ) (f f1 f2 f3 : ℝ → ℝ) : Prop where
  d0 : ∀ y ∈ U, HasDerivAt f (f1 y) y
  d1 : ∀ y ∈ U, HasDerivAt f1 (f2 y) y
  d2 : ∀ y ∈ U, HasDerivAt f2 (f3 y) y
  c3 : ContinuousOn f3 U

variable {U : Set ℝ} {f f1 f2 f3 g g1 g2 g3 : ℝ → ℝ}

lemma J3.c0 (h : J3 U f f1 f2 f3) : ContinuousOn f U :=
  fun y hy => (h.d0 y hy).continuousAt.continuousWithinAt
lemma J3.c1 (h : J3 U f f1 f2 f3) : ContinuousOn f1 U :=
  fun y hy => (h.d1 y hy).continuousAt.continuousWithinAt
lemma J3.c2 (h : J3 U f f1 f2 f3) : ContinuousOn f2 U :=
  fun y hy => (h.d2 y hy).continuousAt.continuousWithinAt

lemma J3.const (U : Set ℝ) (c : ℝ) : J3 U (fun _ => c) (fun _ => 0) (fun _ => 0) (fun _ => 0) :=
  ⟨fun _ _ => hasDerivAt_const _ _, fun _ _ => hasDerivAt_const _ _,
    fun _ _ => hasDerivAt_const _ _, continuousOn_const⟩

lemma J3.add (hf : J3 U f f1 f2 f3) (hg : J3 U g g1 g2 g3) :
    J3 U (fun y => f y + g y) (fun y => f1 y + g1 y) (fun y => f2 y + g2 y)
      (fun y => f3 y + g3 y) :=
  ⟨fun y hy => (hf.d0 y hy).add (hg.d0 y hy), fun y hy => (hf.d1 y hy).add (hg.d1 y hy),
    fun y hy => (hf.d2 y hy).add (hg.d2 y hy), hf.c3.add hg.c3⟩

lemma J3.sub (hf : J3 U f f1 f2 f3) (hg : J3 U g g1 g2 g3) :
    J3 U (fun y => f y - g y) (fun y => f1 y - g1 y) (fun y => f2 y - g2 y)
      (fun y => f3 y - g3 y) :=
  ⟨fun y hy => (hf.d0 y hy).sub (hg.d0 y hy), fun y hy => (hf.d1 y hy).sub (hg.d1 y hy),
    fun y hy => (hf.d2 y hy).sub (hg.d2 y hy), hf.c3.sub hg.c3⟩

lemma J3.cmul (c : ℝ) (hf : J3 U f f1 f2 f3) :
    J3 U (fun y => c * f y) (fun y => c * f1 y) (fun y => c * f2 y) (fun y => c * f3 y) :=
  ⟨fun y hy => (hf.d0 y hy).const_mul c, fun y hy => (hf.d1 y hy).const_mul c,
    fun y hy => (hf.d2 y hy).const_mul c, hf.c3.const_smul c⟩

lemma J3.mul (hf : J3 U f f1 f2 f3) (hg : J3 U g g1 g2 g3) :
    J3 U (fun y => f y * g y) (fun y => f1 y * g y + f y * g1 y)
      (fun y => f2 y * g y + 2 * (f1 y * g1 y) + f y * g2 y)
      (fun y => f3 y * g y + 3 * (f2 y * g1 y) + 3 * (f1 y * g2 y) + f y * g3 y) := by
  refine ⟨fun y hy => (hf.d0 y hy).mul (hg.d0 y hy), fun y hy => ?_, fun y hy => ?_, ?_⟩
  · have := ((hf.d1 y hy).mul (hg.d0 y hy)).add ((hf.d0 y hy).mul (hg.d1 y hy))
    convert this using 1; ring
  · have := (((hf.d2 y hy).mul (hg.d0 y hy)).add (((hf.d1 y hy).mul (hg.d1 y hy)).const_mul 2)).add
      ((hf.d0 y hy).mul (hg.d2 y hy))
    convert this using 1; ring
  · have a := hf.c0; have b := hf.c1; have c := hf.c2; have d := hf.c3
    have a' := hg.c0; have b' := hg.c1; have c' := hg.c2; have d' := hg.c3
    exact (((d.mul a').add (continuousOn_const.mul (c.mul b'))).add
      (continuousOn_const.mul (b.mul c'))).add (a.mul d')

lemma J3.exp (hf : J3 U f f1 f2 f3) :
    J3 U (fun y => Real.exp (f y)) (fun y => Real.exp (f y) * f1 y)
      (fun y => Real.exp (f y) * (f2 y + f1 y ^ 2))
      (fun y => Real.exp (f y) * (f3 y + 3 * (f1 y * f2 y) + f1 y ^ 3)) := by
  refine ⟨fun y hy => (hf.d0 y hy).exp, fun y hy => ?_, fun y hy => ?_, ?_⟩
  · refine (((hf.d0 y hy).exp).mul (hf.d1 y hy)).congr_deriv ?_
    ring
  · have := ((hf.d0 y hy).exp).mul (((hf.d2 y hy)).add ((hf.d1 y hy).fun_pow 2))
    refine this.congr_deriv ?_
    simp; ring
  · exact hf.c0.rexp.mul ((hf.c3.add (continuousOn_const.mul (hf.c1.mul hf.c2))).add (hf.c1.pow 3))

lemma J3.log (hf : J3 U f f1 f2 f3) (hp : ∀ y ∈ U, 0 < f y) :
    J3 U (fun y => Real.log (f y)) (fun y => f1 y / f y)
      (fun y => f2 y / f y - f1 y ^ 2 / f y ^ 2)
      (fun y => f3 y / f y - 3 * (f1 y * f2 y) / f y ^ 2 + 2 * f1 y ^ 3 / f y ^ 3) := by
  refine ⟨fun y hy => (hf.d0 y hy).log (hp y hy).ne', fun y hy => ?_, fun y hy => ?_, ?_⟩
  · have h0 := (hp y hy).ne'
    refine ((hf.d1 y hy).div (hf.d0 y hy) h0).congr_deriv ?_
    field_simp
  · have h0 := (hp y hy).ne'
    have := ((hf.d2 y hy).div (hf.d0 y hy) h0).sub (((hf.d1 y hy).fun_pow 2).div
      ((hf.d0 y hy).fun_pow 2) (pow_ne_zero 2 h0))
    refine this.congr_deriv ?_
    field_simp; ring
  · have hne : ∀ y ∈ U, f y ≠ 0 := fun y hy => (hp y hy).ne'
    refine ((hf.c3.div hf.c0 hne).sub ?_).add ?_
    · exact (continuousOn_const.mul (hf.c1.mul hf.c2)).div (hf.c0.pow 2)
        fun y hy => pow_ne_zero 2 (hne y hy)
    · exact (continuousOn_const.mul (hf.c1.pow 3)).div (hf.c0.pow 3)
        fun y hy => pow_ne_zero 3 (hne y hy)

lemma J3.explin (U : Set ℝ) (l : ℝ) :
    J3 U (fun y => Real.exp (l * y)) (fun y => l * Real.exp (l * y))
      (fun y => l ^ 2 * Real.exp (l * y)) (fun y => l ^ 3 * Real.exp (l * y)) := by
  have h : ∀ y, HasDerivAt (fun y => Real.exp (l * y)) (l * Real.exp (l * y)) y := fun y => by
    have := ((hasDerivAt_id' y).const_mul l).exp
    refine this.congr_deriv ?_
    ring
  refine ⟨fun y _ => h y, fun y _ => ?_, fun y _ => ?_, by fun_prop⟩
  · convert (h y).const_mul l using 1; ring
  · convert (h y).const_mul (l ^ 2) using 1; ring

lemma J3.inv : J3 (Ioi (0:ℝ)) (fun y => y⁻¹) (fun y => -(y ^ 2)⁻¹) (fun y => 2 * (y ^ 3)⁻¹)
    (fun y => -6 * (y ^ 4)⁻¹) := by
  refine ⟨fun y hy => hasDerivAt_inv (ne_of_gt hy), fun y hy => ?_, fun y hy => ?_, ?_⟩
  · have h0 : y ≠ 0 := ne_of_gt hy
    have := ((hasDerivAt_pow 2 y).inv (pow_ne_zero 2 h0)).neg
    convert this using 1; field_simp; ring
  · have h0 : y ≠ 0 := ne_of_gt hy
    have := ((hasDerivAt_pow 3 y).inv (pow_ne_zero 3 h0)).const_mul 2
    convert this using 1; field_simp; ring
  · intro y hy
    have h0 : y ≠ 0 := ne_of_gt hy
    exact ((continuousAt_const.mul ((continuous_pow 4).continuousAt.inv₀ (pow_ne_zero 4 h0)))).continuousWithinAt

lemma J3.iter (h : J3 U f f1 f2 f3) (hU : IsOpen U) {x : ℝ} (hx : x ∈ U) :
    iteratedDeriv 3 f x = f3 x :=
  (S05.chain3 hU hx h.d0 h.d1 h.d2 h.c3).2


lemma J3.congr (h : J3 U f f1 f2 f3) {g : ℝ → ℝ} (e : ∀ y, g y = f y) : J3 U g f1 f2 f3 := by
  have : g = f := funext e
  subst this; exact h


def kp1 (F1 : ℝ) : ℝ := -F1
def kp2 (F1 F2 : ℝ) : ℝ := -(F2 - F1 ^ 2) + 2 * F1
def kp3 (F1 F2 F3 : ℝ) : ℝ := -(F3 - 3 * F1 * F2 + 2 * F1 ^ 3) + 3 * (F2 - F1 ^ 2) - 6 * F1
def gg1 (F1 : ℝ) : ℝ := kp1 F1
def gg2 (F1 F2 : ℝ) : ℝ := kp2 F1 F2 + kp1 F1 ^ 2
def gg3 (F1 F2 F3 : ℝ) : ℝ := kp3 F1 F2 F3 + 3 * kp1 F1 * kp2 F1 F2 + kp1 F1 ^ 3
def uu2 (F1 S1 S2 : ℝ) : ℝ := S2 - 2 * S1 * F1
def uu3 (F1 F2 S1 S2 S3 : ℝ) : ℝ := S3 - 3 * S2 * F1 + 3 * S1 * (2 * F1 ^ 2 - F2)
def WW (F1 F2 S1 S2 S3 : ℝ) : ℝ :=
  3 * gg2 F1 F2 * S1 + 3 * gg1 F1 * uu2 F1 S1 S2 + uu3 F1 F2 S1 S2 S3
def HH (a F1 S1 S2 : ℝ) : ℝ := 3 * gg1 F1 * S1 ^ 2 + 3 * S1 * uu2 F1 S1 S2 + 3 * a * S1 ^ 3

lemma G_iter {w α β : Fin 3 → ℝ} {E0 E1 Y0 Y1 : ℝ} {F1 F2 F3 S1 S2 S3 : ℝ → ℝ}
    (hF : J3 (Ioi 0) (cF w α) F1 F2 F3) (hS : J3 (Ioi 0) (cS w α β) S1 S2 S3)
    (hF0 : cF w α 1 = 1) (hS0 : cS w α β 1 = 0) (hFp : ∀ q, 0 < cF w α q)
    (hQp : ∀ q, 0 < cQ w α β E0 E1 Y0 Y1 q) (hL : 0 < E0 ^ 2 + E1 ^ 2) :
    iteratedDeriv 3 (cG w α β E0 E1 Y0 Y1) 1 =
      Real.sqrt (E0 ^ 2 + E1 ^ 2) *
        (gg3 (F1 1) (F2 1) (F3 1) - (E0 * Y0 + E1 * Y1) / (E0 ^ 2 + E1 ^ 2) *
            WW (F1 1) (F2 1) (S1 1) (S2 1) (S3 1) +
          ((Y0 ^ 2 + Y1 ^ 2) * (E0 ^ 2 + E1 ^ 2) - (E0 * Y0 + E1 * Y1) ^ 2) / (E0 ^ 2 + E1 ^ 2) ^ 2 *
            HH ((E0 * Y0 + E1 * Y1) / (E0 ^ 2 + E1 ^ 2)) (F1 1) (S1 1) (S2 1)) := by
  have hA0 := (hF.cmul E0).sub (hS.cmul Y0)
  have hA1 := (hF.cmul E1).sub (hS.cmul Y1)
  have hQ := ((hA0.mul hA0).add (hA1.mul hA1)).congr (g := cQ w α β E0 E1 Y0 Y1)
    (fun y => by simp only [cQ]; ring)
  have hlF := hF.log (fun y _ => hFp y)
  have hlQ := hQ.log (fun y _ => hQp y)
  have hc := (((J3.const (Ioi (0:ℝ)) 1).add J3.inv).cmul (-1)).mul hlF |>.add (hlQ.cmul (1/2))
  have hch := hc.congr (g := ch w α β E0 E1 Y0 Y1) (fun y => by simp only [ch]; ring)
  have hG := hch.exp.congr (g := cG w α β E0 E1 Y0 Y1) (fun y => rfl)
  rw [hG.iter isOpen_Ioi (by norm_num : (1:ℝ) ∈ Ioi 0)]
  have hQ1 : cQ w α β E0 E1 Y0 Y1 1 = E0 ^ 2 + E1 ^ 2 := by simp [cQ, hF0, hS0]
  have hch1 : Real.exp (ch w α β E0 E1 Y0 Y1 1) = Real.sqrt (E0 ^ 2 + E1 ^ 2) := by
    rw [ch, hF0, hQ1, Real.log_one, mul_zero, zero_add, Real.sqrt_eq_rpow, Real.rpow_def_of_pos hL,
      mul_comm]
  rw [hch1, hQ1, hF0, hS0]
  simp only [Real.log_one, inv_one, one_pow]
  simp only [gg3, gg2, gg1, kp1, kp2, kp3, WW, HH, uu2, uu3]
  have hL0 : E0 ^ 2 + E1 ^ 2 ≠ 0 := hL.ne'
  field_simp
  ring

end S05

namespace S05
open MeasureTheory Set Filter Real Interval Topology

/-- `F_k(t)` at `p = 1` on an edge of the `r = 1`, `w = 1` hexagon. -/
noncomputable def hF (k : ℕ) (t : ℝ) : ℝ :=
  t / 2 * Real.log (t / 2) ^ k + (1 - t) / 2 * Real.log ((1 - t) / 2) ^ k + 1 / 2 * Real.log (1 / 2) ^ k
/-- `S_k(t)` at `p = 1`. -/
noncomputable def hS (k : ℕ) (t : ℝ) : ℝ :=
  1 / 2 * Real.log (t / 2) ^ k - 1 / 2 * Real.log ((1 - t) / 2) ^ k

lemma hasDerivAt_hF (k : ℕ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    HasDerivAt (hF k) (hS k t + k * hS (k - 1) t) t := by
  have h0 : t / 2 ≠ 0 := by linarith [ht.1]
  have h1 : (1 - t) / 2 ≠ 0 := by linarith [ht.2]
  have dP : HasDerivAt (fun s => Real.log (s / 2)) (1 / t) t := by
    have := ((hasDerivAt_id' t).div_const 2).log h0
    refine this.congr_deriv ?_; field_simp
  have dM : HasDerivAt (fun s => Real.log ((1 - s) / 2)) (-(1 / (1 - t))) t := by
    have := (((hasDerivAt_id' t).const_sub 1).div_const 2).log h1
    refine this.congr_deriv ?_
    have : 1 - t ≠ 0 := by linarith [ht.2]
    field_simp
  have := ((((hasDerivAt_id' t).div_const 2).mul (dP.pow k)).add
    ((((hasDerivAt_id' t).const_sub 1).div_const 2).mul (dM.pow k))).add_const
    (1 / 2 * Real.log (1 / 2) ^ k)
  unfold hF
  refine this.congr_deriv ?_
  have ht0 : t ≠ 0 := ht.1.ne'
  have ht1 : 1 - t ≠ 0 := by linarith [ht.2]
  unfold hS
  simp only [Pi.pow_apply, Pi.mul_apply]
  rcases k with _ | k
  · norm_num
  · simp only [Nat.add_sub_cancel, Nat.cast_add, Nat.cast_one, pow_succ]
    field_simp
    ring

lemma hF_cont (k : ℕ) (hk : k ≠ 0) : ContinuousOn (hF k) (Icc 0 1) := by
  have g := G_cont 0 k hk
  have a : ContinuousOn (fun t : ℝ => t / 2 * Real.log (t / 2) ^ k) (Icc 0 1) := by
    have := ContinuousOn.comp (s := Icc (0:ℝ) 1) (f := fun t : ℝ => t / 2) g
      ((continuous_id.div_const (2:ℝ)).continuousOn)
      (fun t (ht : t ∈ Icc (0:ℝ) 1) => ⟨by have := ht.1; simp; linarith, by have := ht.2; simp; linarith⟩)
    refine this.congr fun t _ => ?_
    simp
  have b : ContinuousOn (fun t : ℝ => (1 - t) / 2 * Real.log ((1 - t) / 2) ^ k) (Icc 0 1) := by
    have := ContinuousOn.comp (s := Icc (0:ℝ) 1) (f := fun t : ℝ => (1 - t) / 2) g
      (((continuous_const.sub continuous_id).div_const (2:ℝ)).continuousOn)
      (fun t (ht : t ∈ Icc (0:ℝ) 1) => ⟨by have := ht.2; simp; linarith, by have := ht.1; simp; linarith⟩)
    refine this.congr fun t _ => ?_
    simp
  exact (a.add b).add continuousOn_const

lemma gg3_eq (F1 F2 F3 : ℝ) :
    gg3 F1 F2 F3 = -F3 + 6 * (F1 * F2) - 6 * F1 ^ 3 + 3 * F2 - 9 * F1 ^ 2 - 6 * F1 := by
  simp only [gg3, kp1, kp2, kp3]; ring

lemma hF_one (k : ℕ) : hF k 1 = Real.log (1 / 2) ^ k := by simp [hF]; ring
lemma hF_zero (k : ℕ) : hF k 0 = Real.log (1 / 2) ^ k := by simp [hF]; ring

/-- The `L`-part: `∫_0^1 (g''' - a W) dt = 6c - 6c² + c³`, `c = log 2`. -/
lemma Lpart (a0 : ℝ)
    (hi : IntervalIntegrable (fun t => gg3 (hF 1 t) (hF 2 t) (hF 3 t) -
      (a0 + t) * WW (hF 1 t) (hF 2 t) (hS 1 t) (hS 2 t) (hS 3 t)) volume 0 1) :
    (∫ t in (0:ℝ)..1, (gg3 (hF 1 t) (hF 2 t) (hF 3 t) -
      (a0 + t) * WW (hF 1 t) (hF 2 t) (hS 1 t) (hS 2 t) (hS 3 t))) =
      6 * Real.log 2 - 6 * Real.log 2 ^ 2 + Real.log 2 ^ 3 := by
  set τ : ℝ → ℝ := fun t => (a0 + t) * (-hF 3 t + 6 * (hF 1 t * hF 2 t) - 6 * hF 1 t ^ 3 +
    3 * hF 2 t - 9 * hF 1 t ^ 2 - 6 * hF 1 t) with hτ
  have hc : ContinuousOn τ (Icc 0 1) := by
    have c1 := hF_cont 1 one_ne_zero
    have c2 := hF_cont 2 two_ne_zero
    have c3 := hF_cont 3 (by norm_num)
    have e1 : ContinuousOn (fun t => 6 * (hF 1 t * hF 2 t)) (Icc 0 1) := (c1.mul c2).const_smul (6:ℝ)
    have e2 : ContinuousOn (fun t => 6 * hF 1 t ^ 3) (Icc 0 1) := (c1.pow 3).const_smul (6:ℝ)
    have e3 : ContinuousOn (fun t => 3 * hF 2 t) (Icc 0 1) := c2.const_smul (3:ℝ)
    have e4 : ContinuousOn (fun t => 9 * hF 1 t ^ 2) (Icc 0 1) := (c1.pow 2).const_smul (9:ℝ)
    have e5 : ContinuousOn (fun t => 6 * hF 1 t) (Icc 0 1) := c1.const_smul (6:ℝ)
    have e6 : ContinuousOn (fun t : ℝ => a0 + t) (Icc 0 1) := (continuous_const.add continuous_id).continuousOn
    exact e6.mul (((((c3.neg.add e1).sub e2).add e3).sub e4).sub e5)
  have hd : ∀ t ∈ Ioo (0:ℝ) 1, HasDerivAt τ (gg3 (hF 1 t) (hF 2 t) (hF 3 t) -
      (a0 + t) * WW (hF 1 t) (hF 2 t) (hS 1 t) (hS 2 t) (hS 3 t)) t := by
    intro t ht
    have d1 := hasDerivAt_hF 1 ht
    have d2 := hasDerivAt_hF 2 ht
    have d3 := hasDerivAt_hF 3 ht
    have := ((hasDerivAt_id' t).const_add a0).mul
      ((((((d3.neg.add ((d1.mul d2).const_mul 6)).sub ((d1.pow 3).const_mul 6)).add
        (d2.const_mul 3)).sub ((d1.pow 2).const_mul 9)).sub (d1.const_mul 6)))
    refine this.congr_deriv ?_
    rw [gg3_eq]
    have e0 : hS 0 t = 0 := by simp [hS]
    simp only [WW, gg2, gg1, kp1, kp2, uu2, uu3, Nat.cast_ofNat, Nat.cast_one, e0]
    norm_num
    ring
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one hc hd hi]
  simp only [hτ, hF_one, hF_zero]
  have : Real.log (1 / 2) = -Real.log 2 := by rw [one_div, Real.log_inv]
  rw [this]; ring

end S05

namespace S05
open MeasureTheory Set Filter Real Interval Topology

lemma J3.congr4 {U : Set ℝ} {f f1 f2 f3 g g1 g2 g3 : ℝ → ℝ} (h : J3 U f f1 f2 f3)
    (e0 : ∀ y, g y = f y) (e1 : ∀ y, g1 y = f1 y) (e2 : ∀ y, g2 y = f2 y)
    (e3 : ∀ y, g3 y = f3 y) : J3 U g g1 g2 g3 := by
  have a : g = f := funext e0
  have b : g1 = f1 := funext e1
  have c : g2 = f2 := funext e2
  have d : g3 = f3 := funext e3
  subst a b c d; exact h

/-- `k`-th `p`-derivative of `Σ_j c_j exp(log α_j p)`. -/
noncomputable def Ed (c α : Fin 3 → ℝ) (k : ℕ) (q : ℝ) : ℝ :=
  ∑ j, c j * (Real.log (α j) ^ k * Real.exp (Real.log (α j) * q))

lemma Ed_J3 (c α : Fin 3 → ℝ) :
    J3 (Ioi 0) (fun q => ∑ j, c j * Real.exp (Real.log (α j) * q)) (Ed c α 1) (Ed c α 2)
      (Ed c α 3) := by
  have := (((J3.explin (Ioi (0:ℝ)) (Real.log (α 0))).cmul (c 0)).add
    ((J3.explin (Ioi (0:ℝ)) (Real.log (α 1))).cmul (c 1))).add
    ((J3.explin (Ioi (0:ℝ)) (Real.log (α 2))).cmul (c 2))
  refine this.congr4 ?_ ?_ ?_ ?_ <;> intro y <;> simp [Ed, Fin.sum_univ_three]

lemma cF_J3 (w α : Fin 3 → ℝ) : J3 (Ioi 0) (cF w α) (Ed w α 1) (Ed w α 2) (Ed w α 3) :=
  (Ed_J3 w α).congr (fun y => rfl)

lemma cS_J3 (w α β : Fin 3 → ℝ) :
    J3 (Ioi 0) (cS w α β) (Ed (fun j => w j * β j / α j) α 1) (Ed (fun j => w j * β j / α j) α 2)
      (Ed (fun j => w j * β j / α j) α 3) :=
  (Ed_J3 _ α).congr (fun y => rfl)

lemma hex_half (H : Hex) (hH : IsHex H) (hw : ∀ j, H.w j = 1) (hr : ∀ k, rr H k = 1)
    (k : Fin 3) : A0 H k = 1 / 2 ∧ A1 H k = 1 / 2 ∧ Bv H k = 1 / 2 := by
  have e := fun k => S05.proof_8 H hH k
  have r := fun k => hr k
  have tneg : ∀ j x, term H.u H.w j (-x) = term H.u H.w j x := by
    intro j x; simp [term, dot_neg']
  have rB : ∀ k, term H.u H.w (k + 2) (eV1 H k) = Bv H k := by
    intro k
    have hB := (e k).2.2.2.1
    have := r k
    unfold rr at this
    field_simp at this
    linarith
  obtain ⟨-, a0, b0, -⟩ := e 0
  obtain ⟨-, a1, b1, -⟩ := e 1
  obtain ⟨-, a2, b2, -⟩ := e 2
  have r0 := rB 0; have r1 := rB 1; have r2 := rB 2
  have h0 := r 0; have h1 := r 1; have h2 := r 2
  rw [h0] at b0; rw [h1] at b1; rw [h2] at b2
  simp only [A0, A1, Bv, eV0, eV1] at a0 b0 a1 b1 a2 b2 r0 r1 r2 ⊢
  norm_num at a0 b0 a1 b1 a2 b2 r0 r1 r2
  rw [tneg] at b2 r2
  fin_cases k <;> simp [tneg] <;> refine ⟨?_, ?_, ?_⟩ <;> linarith

end S05

namespace S05
open MeasureTheory Set Filter Real Interval Topology

lemma hex_vecs (H : Hex) (hH : IsHex H) (hw : ∀ j, H.w j = 1) (hr : ∀ k, rr H k = 1)
    (k : Fin 3) (t : ℝ) :
    ew H k = ![1, 1, 1] ∧ eα H k t = ![t / 2, (1 - t) / 2, 1 / 2] ∧
      eβ H k = ![1 / 2, -(1 / 2), 0] := by
  obtain ⟨a0, a1, b⟩ := hex_half H hH hw hr k
  refine ⟨?_, ?_, ?_⟩ <;> funext j <;> fin_cases j <;>
    simp [ew, eα, eβ, hw, hr, a0, a1, b] <;> ring

lemma Gw3_eq (H : Hex) (hH : IsHex H) (hw : ∀ j, H.w j = 1) (hr : ∀ k, rr H k = 1)
    (k : Fin 3) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    Gw H k 3 1 t = Len H k * (gg3 (hF 1 t) (hF 2 t) (hF 3 t) -
        (q0 H k / Lsq H k + t) * WW (hF 1 t) (hF 2 t) (hS 1 t) (hS 2 t) (hS 3 t)) +
      (V02 H k - q0 H k ^ 2 / Lsq H k) / Len H k *
        HH (q0 H k / Lsq H k + t) (hF 1 t) (hS 1 t) (hS 2 t) := by
  obtain ⟨hew, heα, heβ⟩ := hex_vecs H hH hw hr k t
  have hL := Lsq_pos H hH k
  have hLn : Len H k ^ 2 = Lsq H k := Real.sq_sqrt hL.le
  have hLp : 0 < Len H k := Real.sqrt_pos.2 hL
  have ht0 : 0 < t := ht.1
  have ht1 : 0 < 1 - t := by linarith [ht.2]
  have hQ := fun q => cQ_pos (α := eα H k t) (β := eβ H k) (Y0 := ept H k t 0) (Y1 := ept H k t 1)
    (E0 := eE H k 0) (E1 := eE H k 1) (ew_pos H hH k) (edge_det_ne H hH k t) q
  have hFp := fun q => cF_pos (α := eα H k t) (ew_pos H hH k) q
  simp only [hew, heα, heβ] at hQ hFp
  have hE : eE H k 0 ^ 2 + eE H k 1 ^ 2 = Lsq H k := by simp [Lsq, eE, X04.dot]; ring
  have hF0 : cF ![1, 1, 1] ![t / 2, (1 - t) / 2, 1 / 2] 1 = 1 := by
    simp only [cF, Fin.sum_univ_three, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, mul_one, one_mul]
    rw [Real.exp_log (by linarith), Real.exp_log (by linarith), Real.exp_log (by norm_num)]; ring
  have hS0 : cS ![1, 1, 1] ![t / 2, (1 - t) / 2, 1 / 2] ![1 / 2, -(1 / 2), 0] 1 = 0 := by
    simp only [cS, Fin.sum_univ_three, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, mul_one, one_mul]
    rw [Real.exp_log (by linarith), Real.exp_log (by linarith), Real.exp_log (by norm_num)]
    field_simp; ring
  have key := G_iter (cF_J3 _ _) (cS_J3 _ _ _) hF0 hS0 hFp hQ (by rw [hE]; exact hL)
  unfold Gw
  rw [iteratedDerivWithin_of_isOpen isOpen_Ioi (by norm_num : (1:ℝ) ∈ Ioi 0)]
  have hG : (fun q => Gf H k q t) = cG ![1, 1, 1] ![t / 2, (1 - t) / 2, 1 / 2]
      ![1 / 2, -(1 / 2), 0] (eE H k 0) (eE H k 1) (ept H k t 0) (ept H k t 1) := by
    funext q; simp only [Gf, hew, heα, heβ]
  rw [hG, key]
  have eF : ∀ m : ℕ, Ed ![1, 1, 1] ![t / 2, (1 - t) / 2, 1 / 2] m 1 = hF m t := by
    intro m
    simp only [Ed, Fin.sum_univ_three, hF]
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
      Matrix.tail_cons, mul_one, one_mul]
    rw [Real.exp_log (by linarith), Real.exp_log (by linarith), Real.exp_log (by norm_num)]
    ring
  have eS : ∀ m : ℕ, Ed (fun j => ![(1:ℝ), 1, 1] j * ![1 / 2, -(1 / 2), (0:ℝ)] j /
      ![t / 2, (1 - t) / 2, 1 / 2] j) ![t / 2, (1 - t) / 2, 1 / 2] m 1 = hS m t := by
    intro m
    simp only [Ed, Fin.sum_univ_three, hS]
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
      Matrix.tail_cons, mul_one, one_mul]
    rw [Real.exp_log (by linarith), Real.exp_log (by linarith), Real.exp_log (by norm_num)]
    field_simp
    ring
  rw [eF 1, eF 2, eF 3, eS 1, eS 2, eS 3, hE]
  have hq : eE H k 0 * ept H k t 0 + eE H k 1 * ept H k t 1 = q0 H k + t * Lsq H k := by
    simp [eE, ept, q0, Lsq, X04.dot]; ring
  have hY : ept H k t 0 ^ 2 + ept H k t 1 ^ 2 = V02 H k + 2 * t * q0 H k + t ^ 2 * Lsq H k := by
    simp [ept, q0, Lsq, V02, X04.dot]; ring
  rw [hq, hY]
  have hsq : Real.sqrt (Lsq H k) = Len H k := rfl
  rw [hsq]
  have hL0 : Lsq H k ≠ 0 := hL.ne'
  have hLn0 : Len H k ≠ 0 := hLp.ne'
  rw [← hLn]
  field_simp
  ring

end S05

namespace S05
open MeasureTheory Set Filter Real Interval Topology

/-- a finite linear combination of the monomials `t^n log^a t log^b(1-t)`. -/
noncomputable def mcomb {m : ℕ} (kv : Fin m → ℝ) (nv av bv : Fin m → ℕ) (t : ℝ) : ℝ :=
  ∑ i, kv i * (t ^ nv i * Real.log t ^ av i * Real.log (1 - t) ^ bv i)

lemma ii_mcomb {m : ℕ} (kv : Fin m → ℝ) (nv av bv : Fin m → ℕ) :
    IntervalIntegrable (mcomb kv nv av bv) volume 0 1 := by
  unfold mcomb
  have : IntervalIntegrable (∑ i, fun t : ℝ =>
      kv i * (t ^ nv i * Real.log t ^ av i * Real.log (1 - t) ^ bv i)) volume 0 1 :=
    IntervalIntegrable.sum _ fun i _ => (ii_g _ _ _).const_mul _
  convert this using 1
  funext t; simp [Finset.sum_apply]

lemma int_mcomb {m : ℕ} (kv : Fin m → ℝ) (nv av bv : Fin m → ℕ) :
    (∫ t in (0:ℝ)..1, mcomb kv nv av bv t) = ∑ i, kv i * Iint (nv i) (av i) (bv i) := by
  unfold mcomb
  rw [intervalIntegral.integral_finsetSum fun i _ => (ii_g _ _ _).const_mul _]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [intervalIntegral.integral_const_mul]; rfl

noncomputable def hk (a0 c : ℝ) : Fin 11 → ℝ :=
  ![2 + a0, -5 - 3 * a0, 4 + 3 * a0, -1 - a0, 2 * c, -4 * c, 2 * c, -2, 6, -6, 2]
def hn : Fin 11 → ℕ := ![0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]
def ha : Fin 11 → ℕ := ![3, 2, 1, 0, 2, 1, 0, 3, 2, 1, 0]
def hb : Fin 11 → ℕ := ![0, 1, 2, 3, 0, 1, 2, 0, 1, 2, 3]

/-- the `h²`-part integrand expanded into the monomials `t^n log^a t log^b(1-t)`. -/
lemma HH_expand (a0 : ℝ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    HH (a0 + t) (hF 1 t) (hS 1 t) (hS 2 t) = 3 / 8 * mcomb (hk a0 (Real.log 2)) hn ha hb t := by
  have e1 : Real.log (t / 2) = Real.log t - Real.log 2 := Real.log_div ht.1.ne' two_ne_zero
  have e2 : Real.log ((1 - t) / 2) = Real.log (1 - t) - Real.log 2 :=
    Real.log_div (by linarith [ht.2]) two_ne_zero
  have e3 : Real.log (1 / 2) = -Real.log 2 := by rw [one_div, Real.log_inv]
  simp only [HH, gg1, kp1, uu2, hF, hS, e1, e2, e3, mcomb, Fin.sum_univ_succ, Fin.sum_univ_zero,
    hk, hn, ha, hb]
  simp
  ring

lemma ii_HH (a0 : ℝ) :
    IntervalIntegrable (fun t => HH (a0 + t) (hF 1 t) (hS 1 t) (hS 2 t)) volume 0 1 := by
  refine ((ii_mcomb (hk a0 (Real.log 2)) hn ha hb).const_mul (3 / 8)).congr_ae ?_
  filter_upwards [ae_Ioo01] with t ht
  exact (HH_expand a0 ht).symm

lemma int_HH (a0 : ℝ) :
    (∫ t in (0:ℝ)..1, HH (a0 + t) (hF 1 t) (hS 1 t) (hS 2 t)) =
      Real.log 2 * π ^ 2 / 4 - π ^ 2 / 2 - 3 / 4 * Z03.z3 := by
  rw [intervalIntegral.integral_congr_ae (ae_io fun t ht => HH_expand a0 ht),
    intervalIntegral.integral_const_mul, int_mcomb]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, hk, hn, ha, hb]
  simp
  rw [I030, I021, I012, I003, I020, I011, I002, I130, I121, I112, I103]
  ring

end S05

namespace S05
open MeasureTheory Set Filter Real Interval Topology

lemma Gf_one_hex (H : Hex) (hH : IsHex H) (hw : ∀ j, H.w j = 1) (hr : ∀ k, rr H k = 1)
    (k : Fin 3) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) : Gf H k 1 t = Len H k := by
  obtain ⟨hew, heα, heβ⟩ := hex_vecs H hH hw hr k t
  have ht0 : 0 < t := ht.1
  have ht1 : 0 < 1 - t := by linarith [ht.2]
  have hL := Lsq_pos H hH k
  have hE : eE H k 0 ^ 2 + eE H k 1 ^ 2 = Lsq H k := by simp [Lsq, eE, X04.dot]; ring
  have hF0 : cF ![1, 1, 1] ![t / 2, (1 - t) / 2, 1 / 2] 1 = 1 := by
    simp only [cF, Fin.sum_univ_three, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, mul_one, one_mul]
    rw [Real.exp_log (by linarith), Real.exp_log (by linarith), Real.exp_log (by norm_num)]; ring
  have hS0 : cS ![1, 1, 1] ![t / 2, (1 - t) / 2, 1 / 2] ![1 / 2, -(1 / 2), 0] 1 = 0 := by
    simp only [cS, Fin.sum_univ_three, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, mul_one, one_mul]
    rw [Real.exp_log (by linarith), Real.exp_log (by linarith), Real.exp_log (by norm_num)]
    field_simp; ring
  simp only [Gf, hew, heα, heβ, cG, ch, cQ, hF0, hS0, Real.log_one, mul_zero, zero_add, one_mul,
    zero_mul, sub_zero]
  rw [show eE H k 0 ^ 2 + eE H k 1 ^ 2 = Lsq H k from hE, Len, Real.sqrt_eq_rpow,
    Real.rpow_def_of_pos hL, mul_comm]

lemma Phi0_hex (H : Hex) (hH : IsHex H) (hw : ∀ j, H.w j = 1) (hr : ∀ k, rr H k = 1)
    (k : Fin 3) : Phi H k 0 1 = Len H k := by
  unfold Phi
  rw [intervalIntegral.integral_congr_ae (g := fun _ => Len H k) (ae_io fun t ht => by
    rw [Gw_zero]; exact Gf_one_hex H hH hw hr k ht)]
  simp

lemma Ee_hex (H : Hex) (hH : IsHex H) (hw : ∀ j, H.w j = 1) (hr : ∀ k, rr H k = 1)
    (k : Fin 3) :
    Ee H k = Len H k * (6 * Real.log 2 - 6 * Real.log 2 ^ 2 + Real.log 2 ^ 3) +
      (V02 H k - q0 H k ^ 2 / Lsq H k) / Len H k *
        (Real.log 2 * π ^ 2 / 4 - π ^ 2 / 2 - 3 / 4 * Z03.z3) := by
  set a0 := q0 H k / Lsq H k
  set K := (V02 H k - q0 H k ^ 2 / Lsq H k) / Len H k
  set EL := fun t => gg3 (hF 1 t) (hF 2 t) (hF 3 t) -
    (a0 + t) * WW (hF 1 t) (hF 2 t) (hS 1 t) (hS 2 t) (hS 3 t) with hEL
  set HL := fun t => HH (a0 + t) (hF 1 t) (hS 1 t) (hS 2 t) with hHL
  have hLp : 0 < Len H k := Real.sqrt_pos.2 (Lsq_pos H hH k)
  have hG := Gw_int H hH k (j := 3) le_rfl (p := 1) ⟨by norm_num, by norm_num⟩
  have hae : ∀ᵐ t ∂volume, t ∈ Ι (0:ℝ) 1 → Gw H k 3 1 t = Len H k * EL t + K * HL t :=
    ae_io fun t ht => Gw3_eq H hH hw hr k ht
  have hsum : IntervalIntegrable (fun t => Len H k * EL t + K * HL t) volume 0 1 :=
    hG.congr_ae ((ae_restrict_iff' measurableSet_uIoc).2 hae)
  have hHi : IntervalIntegrable (fun t => K * HL t) volume 0 1 := (ii_HH a0).const_mul K
  have hLi : IntervalIntegrable (fun t => Len H k * EL t) volume 0 1 := by
    have := hsum.sub hHi
    refine this.congr_ae (Filter.Eventually.of_forall fun t => ?_)
    simp
  have hEi : IntervalIntegrable EL volume 0 1 := by
    have := hLi.const_mul (Len H k)⁻¹
    refine this.congr_ae (Filter.Eventually.of_forall fun t => ?_)
    simp only
    field_simp
  rw [Ee_eq H hH, Phi, intervalIntegral.integral_congr_ae hae,
    intervalIntegral.integral_add hLi hHi, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul, Lpart a0 hEi, int_HH a0]

/-- `iteratedDeriv 3` of the hexagon display at `0`. -/
lemma hexDisp_iter (P E z : ℝ) :
    iteratedDeriv 3 (fun d => T12.hexDisp P E z d) 0 =
      (6 * Real.log 2 - 6 * Real.log 2 ^ 2 + Real.log 2 ^ 3) * P +
        E * (Real.log 2 * π ^ 2 / 4 - 3 / 4 * z - π ^ 2 / 2) := by
  set U : Set ℝ := Ioi (-1)
  have hσ : J3 U (fun d => d / (1 + d)) (fun d => 1 / (1 + d) ^ 2) (fun d => -2 / (1 + d) ^ 3)
      (fun d => 6 / (1 + d) ^ 4) := by
    refine ⟨fun y hy => ?_, fun y hy => ?_, fun y hy => ?_, fun y hy => ?_⟩
    · have h0 : 1 + y ≠ 0 := by simp [U] at hy; linarith
      refine ((hasDerivAt_id' y).div ((hasDerivAt_id' y).const_add 1) h0).congr_deriv ?_
      field_simp; ring
    · have h0 : 1 + y ≠ 0 := by simp [U] at hy; linarith
      refine ((hasDerivAt_const y (1:ℝ)).div (((hasDerivAt_id' y).const_add 1).pow 2)
        (pow_ne_zero 2 h0)).congr_deriv ?_
      simp only [Pi.pow_apply]
      field_simp; ring
    · have h0 : 1 + y ≠ 0 := by simp [U] at hy; linarith
      refine ((hasDerivAt_const y (-2:ℝ)).div (((hasDerivAt_id' y).const_add 1).pow 3)
        (pow_ne_zero 3 h0)).congr_deriv ?_
      simp only [Pi.pow_apply]
      field_simp; ring
    · have h0 : 1 + y ≠ 0 := by simp [U] at hy; linarith
      exact (continuousAt_const.div ((continuousAt_const.add continuousAt_id).pow 4)
        (pow_ne_zero 4 h0)).continuousWithinAt
  set α := π ^ 2 * E / 24
  set β := E / 8 * (z + 2 * π ^ 2 / 3)
  have hp : J3 U (fun d => P + α * d ^ 2 - β * d ^ 3) (fun d => 2 * α * d - 3 * β * d ^ 2)
      (fun d => 2 * α - 6 * β * d) (fun _ => -6 * β) := by
    refine ⟨fun y _ => ?_, fun y _ => ?_, fun y _ => ?_, continuousOn_const⟩
    · refine (((hasDerivAt_pow 2 y).const_mul α).const_add P |>.sub
        ((hasDerivAt_pow 3 y).const_mul β)).congr_deriv ?_
      simp; ring
    · refine ((((hasDerivAt_id' y).const_mul (2 * α))).sub
        ((hasDerivAt_pow 2 y).const_mul (3 * β))).congr_deriv ?_
      simp; ring
    · refine ((hasDerivAt_const y (2 * α)).sub ((hasDerivAt_id' y).const_mul (6 * β))).congr_deriv ?_
      ring
  have hJ := ((hσ.cmul (Real.log 2)).exp).mul hp
  have hJ' := hJ.congr (g := fun d => T12.hexDisp P E z d) (fun d => by
    simp only [T12.hexDisp, α, β]
    rw [Real.rpow_def_of_pos (by norm_num : (0:ℝ) < 2)])
  rw [hJ'.iter isOpen_Ioi (by norm_num : (0:ℝ) ∈ Ioi (-1))]
  simp only [α, β]
  norm_num
  ring

end S05

open Set in
theorem S05.proof_13 : S05.stmt_13 := by
  intro H hH hw hr
  obtain ⟨hcd, -, hit⟩ := S05.proof_9 H hH
  refine ⟨hcd, ?_⟩
  rw [hit, S05.hexDisp_iter]
  have hP : S05.Perim H 1 = 2 * ∑ k, S05.Len H k := by
    rw [S05.perim_eq H hH ⟨by norm_num, by norm_num⟩]
    simp only [S05.Phi0_hex H hH hw hr]
  have hE : S05.calE H = 2 * ∑ k, (S05.V02 H k - S05.q0 H k ^ 2 / S05.Lsq H k) / S05.Len H k := by
    unfold S05.calE
    congr 1
    refine Finset.sum_congr rfl fun k _ => ?_
    unfold T12.hsqOverL S05.t12Edge
    simp only
    rw [S05.edgeH_sq H hH k, S05.edgeLen_eq]
  rw [hP, hE]
  simp only [S05.Ee_hex H hH hw hr, Fin.sum_univ_three]
  ring
