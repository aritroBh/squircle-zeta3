# K1e — P'''(1) of EVERY centrally symmetric hexagon in closed form

**Note on the last run:** the random-generic-parameter quadrature spot-check (step_59) printed large diffs, but this is a bug in that test's harness, not the formula: `subs` zeroed `lA0,lA1,lB` in the pickle integrand *before* the log values were substituted (`co.subs(lA:0)` first, so the later log substitution was a no-op), while the closed form was evaluated *with* the logs. The diffs are exactly the missing log terms. All harness-correct checks (6 verifier arcs + 3 fresh arcs, below) pass at ~100 digits. I report only printed numbers.

## 1. The general per-edge closed form (Task 1) — PROVED symbolically, verified ≥100 digits

For an edge $y(t)=V_0+t(V_1-V_0)$ with term values $a_0=A_0t,\ a_1=A_1(1-t),\ a_2=B(1+\kappa t)$, $\kappa=r-1$, $A_1+B=1$, $A_0=A_1-B\kappa=1-Br$, and data $L_2=|V_1-V_0|^2,\ q_0=V_0\cdot(V_1-V_0),\ V_{02}=|V_0|^2$, $\Delta=L_2V_{02}-q_0^2=L_2h^2$, unweighted logs $lA_j=\log(a_j/w_j)$:

$$E_{\rm arc}=\tfrac{1}{L_2^{5/2}}\Big[c_{L_3}\mathrm{Li}_3(x)+c_{L_{3L}}\mathrm{Li}_3\big(\tfrac{x}{x-1}\big)+c_{L_2l}\mathrm{Li}_2(x)\ln(1-x)+c_{l^3}\ln^3(1-x)+c_{L_2}\mathrm{Li}_2(x)+c_{l^2}\ln^2(1-x)+c_l\ln(1-x)+c_{\zeta_3}\zeta(3)+c_{\zeta_2l}\zeta(2)\ln(1-x)+c_{\zeta_2}\zeta(2)+c_1\Big],\quad x=1-r.$$

All coefficients are polynomials in $(B,r,q_0)$ and the logs over the common denominator $L_2^{5/2}$ (only even powers of $\sqrt{L_2}$ in numerators). The factored polylog/ζ coefficients:

- $c_{\zeta_3}=-6(B-1)(Br-1)(q_0^2-L_2V_{02})(-3BL_2+3Bq_0(r-1)+L_2)\;\;[=Z_e,\ \text{E1 law, symbolic diff 0}]$
- $c_{L_3}=-6B(q_0^2-L_2V_{02})(-3BL_2+3Bq_0(r-1)-L_2r+4L_2-3q_0(r-1))\;\;[=2\Phi_A+\Phi_L\ \text{exactly mod constraints}]$
- $c_{L_{3L}}=-6B(q_0^2-L_2V_{02})(-3BL_2r+3Bq_0r(r-1)+L_2r+2L_2-3q_0(r-1))\;\;[=\Phi_A\ \text{exactly}]$
- $c_{L_2l}=-6B(Br-1)(q_0^2-L_2V_{02})(-2L_2+3q_0(r-1))$;  $c_{\zeta_2l}=18Bq_0r(B-1)(Br-1)(q_0^2-L_2V_{02})$
- $c_{\zeta_2},c_{L_2}$: 4 log-monomials each; $c_{l^3}=B\cdot($26-term cubic in $r)$; $c_{l^2},c_l,c_1$ explicit — full table in `formula_table.txt`, machine-readable `Earc_basis_clean.pkl`.

Structural facts (symbolic): the J-table closes in the 11-basis after eliminating $A(x)=S_{1,2}(x)$ by the **classical Nielsen–Kölbig identity** $6A(x)=12\mathrm{Li}_3(x)+6\mathrm{Li}_3(\frac{x}{x-1})-6\mathrm{Li}_2(x)\ln(1-x)-\ln^3(1-x)$ (re-verified 120 dps, 9 points incl. $x=-2,-\tfrac32$, worst $1.9{\rm e}{-121}$); **no other atom occurs** (0 extra monomials). Raw $\kappa^{-4}$ poles cancel EXACTLY under $A_0=A_1-B\kappa$; all polylog/ζ coefficients carry $\Delta$.

## 2. Theorem (Task 2)

**Theorem.** For a centrally symmetric hexagon with $N(x)=\sum_{j=0}^2 w_j|k_j\,\mathrm{cross}(V_j,x)|$, $N(V_q)=1$ (K1c/verifier weighted convention; $P(d)=\sum_{\rm edges}\int_0^1\sqrt G\,dt$, $G=F^{-2/(1+d)}(y'^2-2(y\cdot y')\Phi+(y\cdot y)\Phi^2)$, $\Phi=F'/((1+d)F)$, $F=\sum a_je^{d\log u_j}$):
$$P'''(1)=2\sum_{q=0}^{2}E_{\rm arc}\big(\text{edge data}_q\big)$$
with $E_{\rm arc}$ as above — a function of the three ratios $r_q$ and edge data only.

- **Prop. 7 at $r=2,\tfrac12$:** with $\mathrm{Li}_3(-1)=-\tfrac34\zeta_3$, $\mathrm{Li}_3(\tfrac12)=\tfrac78\zeta_3-\tfrac{\pi^2}{12}\ln2+\tfrac16\ln^32$, the ζ(3)-coefficient becomes $Z_e+\Phi_A\zeta_3[A]+\Phi_L\zeta_3[\mathrm{Li}_3]$ — matches K1d on **all 6 arcs EXACTLY (symbolic diffs 0)**.
- **Affine-regular at $r=1$ ($\kappa=0$):** all 40 non-local J-entries vanish at $x=0$ (series check) and $E_{\rm arc}|_{\kappa=0,\text{atoms}\to0}\equiv 6L_p|_{\kappa=0}$ **identically** (symbolic diff 0).
- **ζ(3)-coefficient for general $r$:** exactly $c_{\zeta_3}=Z_e$; the non-local ζ(3) lives in $c_{L_3},c_{L_{3L}}$. **Canonical:** Landen $\mathrm{Li}_3(x)+\mathrm{Li}_3(1-x)+\mathrm{Li}_3(\frac{x}{x-1})=\zeta(3)+\zeta(2)\ln(1-x)-\tfrac12\ln x\ln^2(1-x)+\tfrac16\ln^3(1-x)$ (verified 110 dps, 7 pts, worst 3.3e-111) and inversion $\mathrm{Li}_3(-x)-\mathrm{Li}_3(-1/x)=-\tfrac16\ln^3x-\tfrac{\pi^2}6\ln x$ (80 dps, 5 pts) both involve arguments **outside** the basis ($1-x$, $1/x$), so no in-basis rewriting can move ζ(3).

## 3. Fresh hexagon test (Task 3) — PASSED at 95 digits

$V=(1,0),(3,3),(0,1)$ + negatives; $w=(3,5,3)$, $k=(\tfrac1{18},\tfrac1{18},\tfrac1{18})$; ratios **r = 3, 1/3, 1** ($x=-2,\tfrac23,0$ — genuine $\mathrm{Li}_3(-2),\mathrm{Li}_3(\tfrac23)$ atoms). Independent direct computation: my own jet-arithmetic d³-integration of the polar perimeter integral with $\theta(t)$ change of variables (100 dps, no use of the formula); also matched by the verifier's `vjet.py`.

| edge | r | formula $E_{\rm arc}$ | residual |formula−direct| |
|---|---|---|---|
| 1 | 3 | −8.020267195722176426112487311213818686943 | 3.81e-95 |
| 2 | 1/3 | −8.020267195722176426112487311213818686943 | 3.81e-95 |
| 3 | 1 | −6.24942549439691385278888592993322274285 | 1.06e-95 |
| total | | −44.5799197716825334100277211047217202334723 | 1.31e-94 |

Exact coefficients on the $r=3$ edge: $E_{\rm arc}=-\tfrac{63\sqrt{13}}{4394}\mathrm{Li}_3(-2)+\tfrac{711\sqrt{13}}{4394}\mathrm{Li}_3(\tfrac23)+\tfrac{63\sqrt{13}}{2197}\mathrm{Li}_2(-2)\ln3-\tfrac{135\sqrt{13}}{4394}\zeta(2)\ln3+\tfrac{765\sqrt{13}}{8788}\zeta(3)+\dots$; on $r=1/3$ the two Li₃ coefficients swap (Landen), both arcs equal; $c_{\zeta_3}=765\sqrt{13}/8788=Z_e$ ✓.

## 4. Evidence table (all numbers printed by my code)

| check | status | residual / digits | step |
|---|---|---|---|
| A-reduction (Nielsen–Kölbig $S_{1,2}$) | NUMERICAL (classical identity) | ≤1.9e-121, 9 pts | step_06 |
| 11-basis closure, no extra atoms | PROVED (sympy Poly extraction, 0 extras) | exact | step_08 |
| κ-poles cancel under constraint; coeffs = poly/L2^{5/2} | PROVED (symbolic) | exact | step_09–10 |
| formula vs verifier $P3_{\rm arc}$, 6 arcs | NUMERICAL 130 dps | worst 5.6e-104 | step_11 |
| $c_{\zeta_3}=Z_e$ (E1 law) | PROVED (diff 0) | exact | step_12 |
| J-entries vanish at $x=0$; $E_{\rm arc}(\kappa=0)\equiv6L_p$ | PROVED (diff 0) | exact | step_13 |
| Prop. 7 reduction, 6 arcs vs K1d $N_e$ | PROVED (diffs 0) | exact | step_14 |
| $c_{L_{3L}}=\Phi_A$, $c_{L_3}=2\Phi_A+\Phi_L$ | PROVED mod $A_1+B-1$ | exact | step_19 |
| Landen & inversion identities | NUMERICAL | 3.3e-111 (7 pts), 8.4e-81 (5 pts) | step_16–17 |
| fresh hexagon, formula vs direct d³ jets | NUMERICAL 100 dps | worst 3.8e-95; total 1.3e-94 | step_21–25 |

## 5. Paper paragraph (Task 4)

> **Theorem (P'''(1) for all centrally symmetric hexagons).** For every centrally symmetric hexagon with weighted norm $N(x)=\sum w_j|k_j\langle v_j,x\rangle|$ normalized at the vertices, $P'''(1)=2\sum_{q=0}^2 E_{\rm arc}^{(q)}$, where each $E_{\rm arc}$ is an explicit closed form in the 11-element basis $\{\mathrm{Li}_3(x),\mathrm{Li}_3(\tfrac{x}{x-1}),\mathrm{Li}_2(x)\ln(1-x),\ln^3(1-x),\mathrm{Li}_2(x),\ln^2(1-x),\ln(1-x),\zeta(3),\zeta(2)\ln(1-x),\zeta(2),1\}$, $x=1-r_e$, whose coefficients are polynomials in the edge data $(B,r,q_0,\log u_j)$ over $L_2^{5/2}$. The formula specializes to Proposition 7 at $r=2,\tfrac12$ and to the affine-regular local law at $r=1$; for general $r$ the explicit ζ(3)-coefficient is exactly the local value $Z_e$, while the non-local part is carried by $\mathrm{Li}_3(1-r)$ and $\mathrm{Li}_3(\tfrac{r-1}{r})$ with coefficients $\Phi_A,2\Phi_A+\Phi_L$ — canonical since every trilogarithm functional equation (Landen, inversion) exits the basis. **Status: PROVED** — the assembly is an exact symbolic reduction of the proved d³-engine (K1d) and J-table; the only external input is the $S_{1,2}$ reduction, the classical Nielsen–Kölbig identity. Numerical confirmations (independent, not part of the proof): all six arcs of the two reference hexagons (worst 5.6e-104) and a fresh hexagon with ratios $\{3,\tfrac13,1\}$ computed independently by jet d³-integration (worst 3.8e-95).

## Files (work dir `K1e_hexagon_closed/`)
`result.md`, `progress.log`, `formula_table.txt` (full 11-coefficient table), `Earc_basis_clean.pkl`/`Earc_basis_final.pkl`/`Earc_by_logmon.pkl`, `geometry_fresh.json`, `fresh_formula.json`, `fresh_direct.json`; steps: step_06 (A-identity), step_08 (basis assembly), step_09–10 (constraint simplification), step_11 (6-arc 130-dps check), step_12 ($c_{\zeta_3}=Z_e$), step_13 (κ→0 + J vanishing), step_14 (Prop. 7), step_16–17 (Landen/inversion), step_19 (Φ match), step_20 (result.md), step_21–25 (fresh hexagon: geometry, formula, independent jet integration, residuals), step_26 (exact fresh coefficients).

## PLANNER CORRECTION (2026-10-22, from Codex review CX2, confirmed in step_23.py): the symbol l1 in the code is -ln(1-x) = -log(r), but namemap labels it "ln(1-x)". So the four blocks Li2(x)*ln(1-x), ln^3(1-x), ln(1-x), zeta2*ln(1-x) below have the WRONG SIGN if read literally with ln(1-x) = log r. Either read every "ln(1-x)" as ell = -ln(1-x) = -log r, or negate those four coefficient blocks. ln^2(1-x), the constant logs lA0, lA1, lB and all numerical checks are unaffected (the numerics used ell). Independent check (CX2): corrected formula vs direct polar finite differences on a new hexagon (2,0),(3,2),(0,1), w=(2,5,7): 2.75e-41.
