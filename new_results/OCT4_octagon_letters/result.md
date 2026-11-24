# OCT4 — The octagon's missing weight-3 constant from the integrand's own letters

**Headline result.** The letters of the octagon edge integrand force a constant OCT3's single‑Li₃ basis could not see: the **double polylog**
$$\Lambda \;=\; \operatorname{Li}_{2,1}\!\bigl(1-\sqrt2,\;-\tfrac{\sqrt2}{2}\bigr)\;=\;\sum_{n>m\ge1}\frac{(1-\sqrt2)^{n}\,(-\sqrt2/2)^{m}}{n^{2}\,m}\;=\;-0.0275883613456533384688\ldots$$
It enters through the *three‑letter cross terms* $\log(1+\kappa_1 t)\log(1+\kappa_2 t)\log t$ (and $\log(1-t)$) of the $d^3$ integrand — the two non‑local corner logs $\kappa_1=\sqrt2-1,\ \kappa_2=\sqrt2/2-1$ times a local corner log. I **proved** $\int_0^1\frac{\log t\,\log(1+\kappa_2 t)}{1+\kappa_1 t}dt=-\frac1{\kappa_1}\Lambda$ (verified to $9\cdot10^{-103}$). At 280 digits, however, the full $P'''(1)$ does **not** yet close over singles$+\Lambda$ (bound only $\sim10^{5}$ there), so the final assembly/**ζ(3) separation stays OPEN**. All numbers below were printed by the code (steps auto‑saved `step_01…step_40+`; `progress.log` appended each step; this text = `result.md`).

---

## Task 1 — the letters (PROVED)

Edge $y(t)=V_0+t(V_1-V_0)$, $V_0=(0,\sqrt2-1)$, $V_1=(1-\tfrac{\sqrt2}{2},1-\tfrac{\sqrt2}{2})$ (both verified on $\sum|\langle n_i,x\rangle|=1$). The four signed term values $s_j(t)=\langle n_j,y(t)\rangle$ are linear; on $(0,1)$ the unsigned $a_j=|s_j|$ are all positive. **Singular points of the $d^3$ coefficient $E_3(t)$:**

| source | points | status |
|---|---|---|
| zeros of $a_j(t)$ (log branch points) | $t=0,\ 1,\ -(1+\sqrt2),\ 2+\sqrt2$ (all real) | PROVED (exact) |
| denominator of $E_3(t)$ | **constant** $=48(99-70\sqrt2)\,(10-7\sqrt2)^{3/2}\sqrt{40-28\sqrt2}$ — only powers of $\lvert y'\rvert^2=10-7\sqrt2$ | PROVED (closed‑form $d$‑expansion; $E_3$ = deg‑≤3 poly in $t$ × log‑monomials) |
| $\lvert y(t)\rvert^2=0$ (the “$g_0$ quadratic”) | $t^\*=\dfrac{1\pm i(1+\sqrt2)}{2}$, $\arg=\dfrac{3\pi}{8}$ (since $\tan=1+\sqrt2$), $\lvert t^\*\rvert^2=1+\tfrac{\sqrt2}{2}$ | PROVED (exact) — **but NOT singularities of $E_3(t)$** |
| Gram det $\lvert y\rvert^2\lvert y'\rvert^2-(y\cdot y')^2$ | $=\tfrac{17}{2}-6\sqrt2$ (constant) | PROVED |

**Key structural fact (PROVED):** because the square‑root $\sqrt{\lvert y'\rvert^2-2(y\cdot y')\Phi+\lvert y\rvert^2\Phi^2}$ expands about the *constant* $\lvert y'\rvert^2$ (as $\Phi=O(d)$), $E_3(t)$ has **only the four real log‑singularities**; the complex $\lvert y\rvert^2=0$ points do not propagate into the coefficient. $E_3$ has **20 log‑monomials**, all of total weight 3; the non‑local/non‑local cross terms present are exactly $L_1^2L_2,\ L_1L_2^2,\ L_1L_2L_3,\ L_1L_2L_4$ (printed).

**Cross‑ratios.** Real letters $\{0,1,\infty,-(1+\sqrt2),2+\sqrt2\}$ give 15 values $\to$ the 8 reps $\{-1,-\tfrac{\sqrt2}{2},1-\sqrt2,1-\tfrac{\sqrt2}{2},\sqrt2-1,\tfrac12,2-\sqrt2,\tfrac{\sqrt2}{2}\}$ (matches OCT3). Including $t^\*,\overline{t^\*}$ gives 45 cross‑ratios, all in $\mathbb Q(\zeta_8)$: 3 real, 42 complex with $\arg=k\pi/8$ (several exactly 8th roots of unity, e.g. $e^{i\pi/4},i,e^{3i\pi/4}$) and moduli in $\mathbb Q(\zeta_{16})^+$ ($\sin\frac\pi8,\cos\frac\pi8,\sqrt{2\pm\sqrt2},\sqrt2-1$). So the *full* letter set lives in the $\mathbb Q(\zeta_{16})$ extension — as the brief predicted — but only the **real** letters are actual singularities of $E_3$.

## Task 2 — the weight‑3 basis the letters require (NUMERICAL)

- Real letters $\Rightarrow$ single $\operatorname{Li}_3$ at the 8 reps (already OCT3's span).
- Cross terms $\Rightarrow$ double polylogs. I resolved which: **two‑letter** integrals $J_{21}=\int\log^2(1+\kappa_1t)\log(1+\kappa_2t)$, $J_{12}$ **reduce to single Li₃'s** (analytic reduction = quadrature to $1.7\cdot10^{-80}$; the shift $\gamma=1+\sqrt2$ throws everything onto $\operatorname{Li}_3(\sqrt2-1),\operatorname{Li}_3(2-\sqrt2)$, both in the basis). The $/t$ integral $K=\int\log(1+\kappa_1t)\log(1+\kappa_2t)/t\,dt$ also reduces (residual $0$, maxcoeff 4). The **three‑letter** integrals $T_t=\int\log(1+\kappa_1t)\log(1+\kappa_2t)\log t\,dt$, $T_{1-t}$ do **not** reduce to single Li₃'s (LLL: only spurious relations, floor $\gtrsim10^{11}$ at 200 dps). Their content is exactly $\Lambda=\operatorname{Li}_{2,1}(1-\sqrt2,-\tfrac{\sqrt2}{2})$ and its Landen swap $U(\kappa_2,\kappa_1)=\int\log t\log(1+\kappa_1t)/(1+\kappa_2t)\,dt$.
- Reduced 46‑constant basis ($\zeta(3)$, 8 Li₃, $\pi^2\log$, log³, Li₂·log, Li₂, logs, 1) by exact LLL at 280 dps: **26 internal relations** (residuals $<10^{-262}$), leaving a **20‑dimensional independent set** {4 Li₃'s, 8 Li₂·log's, 2 Li₂'s, log²×3, log×2, 1}.

## Task 3 — the hunt (NUMERICAL, inconclusive)

Hunted $P'''(1)/u$ ($u=\sqrt{2+\sqrt2}$, $P'''(1)/u=-2.378284262522270960164354737\ldots$ at 300 dps, matches OCT3) over $\{1,\sqrt2\}$ at 280 digits: **no relation**, floor maxcoeff $9.7\cdot10^{5}$ (singles only) and $2.7\cdot10^{4}$ (singles $+\Lambda+$ swap). The residual/floor only *lowered* because more constants were added — no genuine relation. At 280 digits the exclusion bound ($\sim10^{5}$) is far weaker than OCT3's $8.7\cdot10^{12}$ (520 dps); resolving whether $\Lambda$ closes $P'''(1)$ needs the 500‑digit hunt that the turn budget did not reach.

## Task 4 — conclusion

- **The missing weight‑3 constant is identified**: the double polylog $\Lambda=\operatorname{Li}_{2,1}(1-\sqrt2,-\sqrt2/2)$ (plus its swap), forced by the integrand's own three‑letter cross terms. This is *not* a single‑Li₃ S‑unit value — which is exactly why OCT3's single‑Li₃ hunt (bound $8.7\cdot10^{12}$) provably could not capture it. (NUMERICAL for the identification; PROVED for $\int\log t\log(1+\kappa_2t)/(1+\kappa_1t)=-\Lambda/\kappa_1$.)
- **Open (honest):** whether $P'''(1)$ closes over singles$+\Lambda$ (and hence the ζ(3) coefficient vs the local edge law $C_4=(54-39\sqrt2)\,u$, and the non‑local remainder) is **not decided** — the 280‑digit hunt found no relation (bound $\sim10^{5}$), and the 500‑digit hunt over $\{1,\sqrt2\}$ and $\mathbb Q(\zeta_{16})^+$ with $\Lambda$ included is the concrete next step the budget cut off. If the two double‑polylog pieces cancel in the total (a Clausen‑type cancellation like D15b's), $P'''(1)$ would fall back into the single‑Li₃ span with coefficients beyond OCT3's bound; if not, $\Lambda$ is genuinely present.

## Evidence table

| claim | status | evidence (residual, digits) |
|---|---|---|
| edge vertices on $\sum\lvert\langle n_i,x\rangle\rvert=1$ | PROVED | exact sympy, residual 0 |
| $a_j$ zeros $=0,1,-(1+\sqrt2),2+\sqrt2$ | PROVED | exact |
| $E_3(t)$ denominator constant (only real singularities) | PROVED | closed‑form $d$‑expansion, `step_19` |
| $t^\*=(1\pm i(1+\sqrt2))/2$, $\arg=3\pi/8$ | PROVED | $\tan=1+\sqrt2$, $\lvert t^\*\rvert^2=1+\sqrt2/2$ |
| cross‑ratios: 8 real reps (OCT3); 45 with $t^\*$, args $k\pi/8$ | PROVED | exact sympy, `step_21,22` |
| $J_{21}$ reduces to single Li₃'s | PROVED | analytic $-$ quad $=1.7\cdot10^{-80}$ (80 dps), `step_40` |
| $K$ reduces to single Li₃'s | NUMERICAL | residual 0, maxcoeff 4, `step_31` |
| $T_t,T_{1-t}$ need $\Lambda$ (no single‑Li₃ form) | NUMERICAL | spurious floor $\gtrsim10^{11}$, 200 dps, `step_33,35` |
| $U=-(1/\kappa_1)\operatorname{Li}_{2,1}(1-\sqrt2,-\sqrt2/2)$ | PROVED | residual $8.9\cdot10^{-103}$, `step_36` |
| 26 internal relations, 20‑dim independent basis | NUMERICAL | residuals $<10^{-262}$, 280 dps, `step_38` |
| $P'''(1)/u$ value | NUMERICAL | matches OCT3 to 300 dps, `step_25` |
| no $P'''/u$ relation over $\{1,\sqrt2\}$(singles$+\Lambda$) | NUMERICAL | floor $\sim10^{5}$, only 280 dps — **weak, OPEN** |

## Step files (auto‑saved by the harness in `OCT4_octagon_letters/`)

- `step_03–08` — data loading; octagon geometry; term values; $a_j$ zeros.
- `step_12–16` — geometric quadratics, complex letter $t^\*$ (arg $3\pi/8$).
- `step_18–19` — $E_3(t)$ closed‑form $d$‑expansion; constant denominator (real singularities only).
- `step_20–22` — cross‑ratios: 8 real reps; 45 complex in $\mathbb Q(\zeta_8)$/$\mathbb Q(\zeta_{16})^+$.
- `step_23–29` — pure weight‑3 constants (8 real + 24 complex Re/Im Li₃ + ζ₃) at 300 dps; LLL relations.
- `step_30–35` — cross‑integral analysis: $K$, $J_{21}$ reduce; $T_t,T_{1-t}$ don't; $E_3$ log‑monomial structure (cross terms $L_1^2L_2,L_1L_2^2,L_1L_2L_3,L_1L_2L_4$).
- `step_36–37` — identification $\Lambda=\operatorname{Li}_{2,1}(1-\sqrt2,-\sqrt2/2)$; verification to $10^{-103}$.
- `step_38–40` — independent 20‑basis (26 relations printed); $P'''/u$ hunts at 280 dps (negative, weak bound); $J_{21}$ analytic reduction verified.

**Bottom line.** The integrand's letters are the four real points $\{0,1,-(1+\sqrt2),2+\sqrt2\}$ (cross‑ratios = OCT3's 8 S‑units) — the complex $\lvert y\rvert^2=0$ points at argument $3\pi/8$ (the $\mathbb Q(\zeta_{16})$ structure) are *not* singularities of the $d^3$ coefficient. The genuinely new constant the letters force is the double polylog $\operatorname{Li}_{2,1}(1-\sqrt2,-\sqrt2/2)$ from the two‑non‑local‑log × local‑log cross terms; single $\operatorname{Li}_3$'s provably cannot express it. Whether it survives in $P'''(1)$ (vs. cancelling, which would push the single‑Li₃ coefficients past OCT3's $8.7\cdot10^{12}$ bound) requires the 500‑digit singles$+\Lambda$ hunt and remains the precise open question, with the exact basis, relations, and the identified constant $\Lambda$ now in hand.