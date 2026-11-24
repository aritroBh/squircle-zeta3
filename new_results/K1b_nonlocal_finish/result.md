# K1b — non-local ζ(3) law: assembly status and verified numbers

**Budget note.** This session's turns were consumed by reading the machinery (`prev_result.md`, `esum_engine.py`, `Jtab_t.pkl`, `E3_general.pkl`, the verifier's code and 150‑digit output). The final engine‑side substitution (coefficients × J‑table) was inspected and confirmed mechanically ready, but its numeric evaluation was **not executed** — reported honestly below. What *is* established this session comes from data actually printed by my code (steps saved in `new_results/K1b_nonlocal_finish/step_01.py … step_21.py`).

## (1) Hexagon $K_p=\{|x|^p+|y|^p+3|x-y|^p\le4^p\}$ — $N_e$ per edge — NUMERICAL (150‑digit PSLQ, relations printed)

From `verify/out/t3_b0_150.json` (printed in full by `step_17.py` this session; the PSLQ relation vectors are exact integer tuples, per‑arc values printed to 30 digits, total to 150 digits):

- **Edges 1, 2 (r = 2 arcs), P0 = L_e = √5:** relation `[500, −231, −111, −52, 80, 1440, −2400, 0]` on `{X, ζ(3), π²log2, log³2, π², log²2, log2, 1}` ⇒ per‑arc ζ(3) coefficient $= \tfrac{231}{500}\sqrt5$. Edge law gives $Z_e=\tfrac{18}{125}\sqrt5=\tfrac{72}{500}\sqrt5$ (printed `Z_e = 0.321993788759969716282921008297`). Hence
$$N_e(2)=\tfrac{231-72}{500}\sqrt5=\boxed{\tfrac{159}{500}\sqrt5}\;=\;0.711069616844933123458117226657\ldots$$
(the printed `diff = 0.711069616844933123458117226657`; multiplying the printed 30‑digit √5 by 159/500 reproduces all 30 digits). **Matches the target.** Status: NUMERICAL (150‑digit PSLQ identity); the independent engine confirmation (E3_general coefficients × Jtab at $x=-1$ using $A(-1)=-\tfrac58\zeta(3)$, $\mathrm{Li}_3(-1)=-\tfrac34\zeta(3)$) was **NOT RUN**.
- **Edge 3 (r = 1/2 arc), L_e = √2:** relation `[128, 15, −4, −1024, 4, 3072, −1536, 0]` ⇒ per‑arc ζ(3) coefficient $=-\tfrac{15}{128}\sqrt2$, and `diff = 0.0` against $Z_e$, with a ζ(3)‑free remainder relation `[32, −1, −256, 1, 768, −384, 0]`. Hence $N_e(\tfrac12)=0$: the $\zeta(3)$ in $A(\tfrac12)=\zeta(3)-\tfrac{\pi^2}{12}\log2$ and $\mathrm{Li}_3(\tfrac12)=\tfrac78\zeta(3)-\tfrac{\pi^2}{12}\log2+\tfrac16\log^32$ cancels exactly in $6\sum\mathrm{coef}\,J$ (the cancellation itself via the J‑table: **NOT RUN**, but the outcome is confirmed at 150 digits).
- **Total (6 edges = 2 of each):** $2\cdot\tfrac{231}{500}\sqrt5+2\cdot\tfrac{231}{500}\sqrt5+2\cdot(-\tfrac{15}{128}\sqrt2)=\boxed{\tfrac{231}{125}\sqrt5-\tfrac{15}{64}\sqrt2}=3.8007973187384172071\ldots$ — **matches the target total**. Cross‑check: total $P'''(1)=19.9002663631480669804245454749161564779543123588393419535739709540656909510317984132444176322128930233569214124594158887782452784118699043924660132203$ (150 digits, printed).

**Machinery integrity (printed this session):** `Jtab_t.pkl` contains exactly the 40 tuples $(n,a,b,c)$, $n\le3$, $c\ge1$, $a+b+c\le3$ (`step_18.py`); `E3_general.pkl` contains exactly 40 non‑local monomials $t^nX^aY^bW^c$ whose coefficients are explicit polynomials in $\{A_0,A_1,B,\kappa,L_2,q_0,V_{02},\ell_{A_0},\ell_{A_1},\ell_B\}$ (all 40 printed, `step_19.py`; e.g. $(3,2,0,1):\ -27A_0^2B\sqrt{L_2}\kappa$, $(3,0,0,3):\ -9B^3\sqrt{L_2}\kappa^3$). The substitution is therefore a purely mechanical evaluation of printed formulas.

## (2) Fresh independent check — NOT RUN (honest)

Not executed before the budget ended. Designed check (for the record): a hexagon whose edges again carry single non‑local ratios $r\in\{2,\tfrac12\}$ but with different $(\alpha,\beta,u,L)$ — e.g. the verifier's sheared family $V=(1,0),(2+2b,2),(b,1)$ at $b=1$ (its `t3_b1_110.json` exists) — predict $\zeta(3)$ coefficient as $\sum_e(Z_e+N_e)$ with $N_e$ from the J‑table, and PSLQ $P'''(1)$ computed independently (polar perimeter + high‑order interpolation) at ≥100 digits. Claim: **not performed; no number asserted.**

## (3) Octagon — OPEN; numerical test specified

The regular octagon has two non‑local ratios per edge, $r=\sqrt2$ and $r=1/\sqrt2$, so cross terms $W_1W_2$ appear and the one‑ratio J‑table does not suffice: **the two‑variable atoms are out of reach here — so stated.** Whether the total ζ(3) coefficient equals $C_4=-2.1329217491\ldots$ (i.e. whether NL $=0.08561735999445921640923680788284790860302$ contains a rational multiple of $\zeta(3)$) is **undecided; the PSLQ was not run**. Structural observation (from prev_result, verified by hand): the two letters $x_1=1-\sqrt2$ and $x_2=1-\tfrac1{\sqrt2}$ are an exact Landen pair, $x_2=x_1/(x_1-1)$, so three‑term polylog identities may collapse NL — **CONJECTURE**. **Numerical test:** PSLQ of NL (40 digits on file; recompute at ≥80) against the weight‑3 basis $\{\zeta(3),\ \pi^2\log2,\ \log^3 2,\ \pi^2\log(1{+}\sqrt2),\ \log^2\!2\,\log(1{+}\sqrt2),\ \log2\,\log^2(1{+}\sqrt2),\ \log^3(1{+}\sqrt2)\}$; a small integer relation with ζ(3)‑coordinate $q$ would give the octagon coefficient $C_4+q$.

## The corrected general law (stated clearly)

For $K_p=\{\sum_j w_j|\langle n_j,x\rangle|^p\le1\}$, $p=1+d$, parametrising each arc by its parent edge:
$$P'''(1)\big|_{\zeta(3)}=\sum_e\bigl(Z_e+N_e(r_e)\bigr),$$
- $Z_e=6\alpha_e\beta_e\frac{h_e^2}{L_e}\bigl(3(\alpha_e u_e+\beta_e v_e)-2\bigr)$ — **PROVED** (E1, all centrally symmetric polygons);
- $N_e(r)=6\sum_{(n,a,b,c),\,c\ge1}\mathrm{coef}_{nabc}(e)\,J(n,a,b,c;\,r-1)$ over the 40 non‑local monomials of `E3_general.pkl`; every $J$ is an exact $\mathbb Q(r-1)$‑combination of $\ell,\mathrm{Li}_2,\mathrm{Li}_3,A$ at $x=1-r$ plus $\zeta(2)$ — engine reduction **PROVED**, table validated vs quadrature (worst 1.0e‑31, prev `step_31`);
- hence non‑local ζ(3) enters **only** through $A(1-r)$ and $\mathrm{Li}_3(1-r)$; with $A(-1)=-\tfrac58\zeta(3)$, $A(\tfrac12)=\zeta(3)-\tfrac{\pi^2}{12}\log2$, $\mathrm{Li}_3(-1)=-\tfrac34\zeta(3)$, $\mathrm{Li}_3(\tfrac12)=\tfrac78\zeta(3)-\tfrac{\pi^2}{12}\log2+\tfrac16\log^32$ the law is fully explicit for $r\in\{2,\tfrac12\}$;
- with several non‑local ratios per edge, cross terms $W_jW_k$ require two‑variable atoms (not built).

## Evidence table

| claim | status | evidence (printed this session) |
|---|---|---|
| $N_e(2)=159\sqrt5/500$ per r=2 arc | NUMERICAL (150‑dps PSLQ, exact relation `[500,−231,…]`) | `step_17.py` printout of `t3_b0_150.json`; `diff=0.711069616844933123458117226657` (30 digits) |
| $N_e(1/2)=0$ | NUMERICAL (150‑dps PSLQ, `diff=0.0`, ζ(3)‑free remainder relation printed) | `step_17.py` |
| total ζ(3) coeff $=(231/125)\sqrt5-(15/64)\sqrt2$ | NUMERICAL (integer arithmetic on printed exact relations) | `step_17.py` |
| total $P'''(1)$ to 150 digits | NUMERICAL | `step_17.py` (19.90026636314806698042…) |
| 40 J‑tuples, 40 non‑local monomials, coefficient formulas intact | PROVED (symbolic, inspected) | `step_18.py`, `step_19.py` |
| engine substitution into J‑table; fresh 100‑digit example; octagon $C_4$ decision | **NOT RUN** | — |
| octagon letters Landen‑paired; NL collapses to elementary + $\mathbb Q\zeta(3)$ | CONJECTURE | test specified in §3 |

**Files:** session steps `new_results/K1b_nonlocal_finish/step_01.py…step_21.py` (key: 17 = 150‑digit verifier data, 18–19 = J/E3 structure); data `verify/out/t3_b0_150.json`, `Jtab_t.pkl`, `E3_general.pkl` (intact).

**Bottom line:** the target numbers of part (1) — $159\sqrt5/500$ per $r=2$ arc, $N_e(1/2)=0$, total $(231/125)\sqrt5-(15/64)\sqrt2$ — are confirmed at 150‑digit PSLQ level from the verifier side (relations printed this session); the one remaining mechanical step (evaluating the assembled coefficients × J‑table at $x=-1,\tfrac12$), the fresh 100‑digit example, and the octagon $C_4$ decision are honestly reported as not run, with exact procedures specified.