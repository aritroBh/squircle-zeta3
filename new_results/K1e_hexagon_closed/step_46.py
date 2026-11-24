
import time
WORK='research/08_theory/K1e_hexagon_closed'
md = r"""# K1e — P'''(1) of EVERY centrally symmetric hexagon in closed form (Prop. 7 for all ratios r)

**Bottom line (in progress).** For a general single-ratio edge we produce $E_{\rm arc}=P'''_{\rm arc}(1)$ as an EXACT
combination of the 11-element basis
$\{\mathrm{Li}_3(x),\ \mathrm{Li}_3(\tfrac{x}{x-1}),\ \mathrm{Li}_2(x)\ln(1-x),\ \ln^3(1-x),\ \mathrm{Li}_2(x),\ \ln^2(1-x),\ \ln(1-x),\ \zeta(3),\ \zeta(2)\ln(1-x),\ \zeta(2),\ 1\}$,
$x=1-r$, with coefficients = (polynomials in $A_1,B,r,q_0$ and the logs $lA_0,lA_1,lB$)/$L_2^{5/2}$ — no other atoms
are needed (the harmonic atom $A(x)=S_{1,2}(x)$ is eliminated by the classical Nielsen–Kölbig identity, and no other
monomial occurs). Verified: symbolically exact against the fixed engine pickle; numerically against all 6 verifier
arcs at $\ge 100$ digits (worst $|{\rm diff}|=5.6{\rm e}{-104}$). Reductions: Prop. 7 at $r=2,\tfrac12$ EXACT;
affine-regular local law at $r=1$ EXACT (symbolic). The $\zeta(3)$-coefficient for general $r$ is exactly the local
$Z_e$ (E1 law); the non-local $\zeta(3)$ at special $r$ is carried by the $\mathrm{Li}_3$ coefficients
($=\Phi_A,\ 2\Phi_A+\Phi_L$ of K1d, exact). Fresh hexagon with $r=3,1/3,1$: see §3.

## 0. Conventions (the weighted convention of K1c / verifier t3_hexagon.py)

Centrally symmetric hexagon with vertices $V_0,V_1,V_2,-V_0,-V_1,-V_2$ and the weighted norm
$N(x)=\sum_{j=0}^2 w_j\,|k_j\,\mathrm{cross}(V_j,x)|$ normalized by $N(V_q)=1$ ($q=0,1,2$).
Perimeter $P(d)=\sum_{\text{6 edges}}\int_0^1 \sqrt{G}\,dt$,
$G=F^{-2/(1+d)}\big(y'\cdot y'-2(y\cdot y')\Phi+(y\cdot y)\Phi^2\big)$, $\Phi=F'/((1+d)F)$,
$F=\sum_j a_j(t)e^{d\,LA_j}$, $LA_j=\log u_j$ with **unweighted** term values $u_j=a_j/w_j$.
On edge $V_0V_1$: term values $a_0=A_0t$, $a_1=A_1(1-t)$, $a_2=B(1+\kappa t)$, $\kappa=r-1$,
$r=\mathrm{cross}(V_{j_o},V_1)/\mathrm{cross}(V_{j_o},V_0)$ the unique non-local ratio.
Edge constraints (from $N(V_0)=N(V_1)=1$): $A_1+B=1$, $A_0=A_1-B\kappa=1-Br$.
Geometric data: $L_2=|V_1-V_0|^2$, $q_0=V_0\cdot(V_1-V_0)$, $V_{02}=|V_0|^2$, $\Delta=L_2V_{02}-q_0^2=L_2h^2>0$.
$P'''(1)=2\sum_{q=0}^{2}E_{\rm arc}^{(q)}$ (central symmetry doubles the three arcs).

## 1. The general per-edge closed form (Task 1)

$$E_{\rm arc}=\frac{1}{L_2^{5/2}}\Big[\;c_{L_3}\mathrm{Li}_3(x)+c_{L_{3L}}\mathrm{Li}_3\Big(\frac{x}{x-1}\Big)+c_{L_2l}\,\mathrm{Li}_2(x)\ln(1-x)+c_{l^3}\ln^3(1-x)+c_{L_2}\mathrm{Li}_2(x)+c_{l^2}\ln^2(1-x)+c_l\ln(1-x)+c_{\zeta_3}\zeta(3)+c_{\zeta_2l}\zeta(2)\ln(1-x)+c_{\zeta_2}\zeta(2)+c_1\Big]$$

with $x=1-r$, $lA_j=\log(a_j/w_j)$ unweighted logs, and (using $A_0=1-Br$, $A_1=1-B$; full machine-readable table in
`formula_table.txt` / `Earc_basis_clean.pkl`):

- $c_{\zeta_3}=-6(B-1)(Br-1)(q_0^2-\Delta')( -3BL_2+3Bq_0(r-1)+L_2)$, $\Delta'=L_2V_{02}-q_0^2$ —
  identically $=Z_e=6A_0A_1(L_2V_{02}-q_0^2)(3A_1L_2+3B\kappa q_0-2L_2)/L_2^{5/2}$ (E1 local law; symbolic diff 0, step_12).
- $c_{L_3}=-6B\,\Delta'\big(-3BL_2+3Bq_0(r-1)-L_2r+4L_2-3q_0(r-1)\big)$ $=2\Phi_A+\Phi_L$ (exact mod constraints, step_19).
- $c_{L_{3L}}=-6B\,\Delta'\big(-3BL_2r+3Bq_0r(r-1)+L_2r+2L_2-3q_0(r-1)\big)$ $=\Phi_A$ (exact).
- $c_{L_2l}=-6B(Br-1)\Delta'(-2L_2+3q_0(r-1))$.
- $c_{\zeta_2 l}=18Bq_0r(B-1)(Br-1)\Delta'$.
- $c_{\zeta_2}=6(B-1)(Br-1)\Delta'\big[(-3BL_2+3Bq_0(r-1)+2L_2-3q_0)lA_1\ …\big]$ (4 log-monomials, see table).
- $c_{L_2}$: 4 log-monomials; $c_{l^3}=B\cdot($degree-3 polynomial in $r$, 26 terms$)$; $c_{l^2},c_l,c_1$ explicit (table).

Key structural facts (all symbolic, `step_08–step_19`):
1. The 40 non-local monomials integrate (J-table) to atoms $\{A(x),\mathrm{Li}_3(x),\mathrm{Li}_2(x),\ln(1-x),\zeta(2)\}$;
   eliminating $A$ via $6A(x)=12\mathrm{Li}_3(x)+6\mathrm{Li}_3(\tfrac{x}{x-1})-6\mathrm{Li}_2(x)\ln(1-x)-\ln^3(1-x)$
   (Nielsen/Kölbig $S_{1,2}$ identity; re-verified 120 dps at 9 points incl. $x=-2,-3/2$, residuals $\le1.9{\rm e}{-121}$)
   closes in the 11-basis: **no other atom occurs** (Poly extraction: 0 extra monomials).
2. Raw coefficients have $\kappa$-poles up to $\kappa^{-4}$; under $A_0=A_1-B\kappa$ **all poles cancel exactly** and every
   coefficient becomes a polynomial in $\kappa$ over $L_2^{5/2}$ (only even powers of $\sqrt{L_2}$ in numerators).
3. All polylog/$\zeta$ coefficients carry the factor $\Delta'=L_2V_{02}-q_0^2=L_2h^2$.

### Numeric check vs verifier (130 dps; atoms at $x=-1,\tfrac12,0$ exact)
| arc | r | $|E_{\rm arc}-P3_{\rm arc}|$ |
|---|---|---|
| b=0 e1 | 2 | 9.0e-131 |
| b=0 e2 | 1/2 | 1.7e-129 |
| b=0 e3 | 1 | 2.5e-130 |
| b=1 e1 | 2 | 5.6e-104 |
| b=1 e2 | 1/2 | 5.2e-104 |
| b=1 e3 | 1 | 2.9e-105 |

## 2. Theorem (Task 2)

**Theorem (K1e).** Let $H$ be a centrally symmetric hexagon with weighted norm $N(x)=\sum w_j|k_j\mathrm{cross}(V_j,x)|$,
$N(V_q)=1$. Then
$$P'''(1)=2\sum_{q=0}^{2}E_{\rm arc}(A_1^{(q)},B^{(q)},r_q,L_2^{(q)},q_0^{(q)},V_{02}^{(q)};\,lA_0^{(q)},lA_1^{(q)},lB^{(q)})$$
with $E_{\rm arc}$ the closed form of §1 (in terms of the three ratios $r_q$ and edge data only).

- **Prop. 7 reduction** ($r=2,\tfrac12$): at $x=-1$ ($r=2$): $\mathrm{Li}_3(x)=-\tfrac34\zeta_3$,
  $\mathrm{Li}_3(\tfrac{x}{x-1})=\mathrm{Li}_3(\tfrac12)=\tfrac78\zeta_3-\tfrac{\pi^2}{12}\ln2+\tfrac16\ln^32$; at $x=\tfrac12$ the two swap.
  Then $c_{\zeta_3}+c_{L_3}\zeta_3[\mathrm{Li}_3]+c_{L_{3L}}\zeta_3[\mathrm{Li}_3(\tfrac{x}{x-1})]$
  $=Z_e+\Phi_A(-\tfrac58)+\Phi_L(-\tfrac34)$ ($r=2$) / $=Z_e+\Phi_A(1)+\Phi_L(\tfrac78)$ ($r=\tfrac12$) — matches K1d's
  $N_e$ on all 6 arcs EXACTLY (symbolic diffs 0, step_14).
- **Affine-regular reduction** ($r=1$, $\kappa=0$): all 40 non-local J-entries vanish at $x=0$ (series check) and
  $E_{\rm arc}|_{\kappa=0,\ \text{atoms}\to0}\equiv 6L_p|_{\kappa=0}$ identically in the parameters (symbolic diff 0,
  step_13) — the local/affine-regular law.
- **$\zeta(3)$-coefficient for general $r$:** in the 11-basis it is exactly $c_{\zeta_3}=Z_e$ (the local E1 value);
  the would-be non-local $\zeta(3)$ is carried by $c_{L_3},c_{L_{3L}}$ and materializes only at special $r$.
  **Canonicality:** the only classical trilog relations are Landen
  $\mathrm{Li}_3(x)+\mathrm{Li}_3(1-x)+\mathrm{Li}_3(\tfrac{x}{x-1})=\zeta(3)+\zeta(2)\ln(1-x)-\tfrac12\ln x\ln^2(1-x)+\tfrac16\ln^3(1-x)$
  (verified 110 dps, 7 points, worst 3.3e-111) and inversion
  $\mathrm{Li}_3(-x)-\mathrm{Li}_3(-1/x)=-\tfrac16\ln^3x-\tfrac{\pi^2}{6}\ln x$ (verified 80 dps, 5 points);
  EVERY such relation involves an argument outside the basis ($1-x$ resp. $1/x$), so no rewriting inside the basis can
  move $\zeta(3)$ between the explicit $\zeta(3)$ term and the $\mathrm{Li}_3$ coefficients: the coefficient is
  **canonical for this basis**.

## 3. Fresh hexagon test (Task 3) — PENDING

## 4. Paper paragraph (Task 4) — PENDING

## Files
`Earc_basis_clean.pkl` (11 coefficients, params $A_1,B,\kappa,L_2,q_0,V_{02},lA$), `Earc_basis_final.pkl`,
`formula_table.txt` (full human-readable table), `Earc_by_logmon.pkl` (20 log-monomial pieces),
`progress.log`, steps `step_06..step_19`.
"""
open(f"{WORK}/result.md","w").write(md)
with open(f"{WORK}/progress.log",'a') as fh:
    fh.write(f"[{time.strftime('%H:%M:%S')}] step_20: wrote result.md (v1): formula, structural facts, 6-arc check, theorem + reductions + canonicality; Tasks 3-4 pending.\n")
print("result.md written", len(md))
