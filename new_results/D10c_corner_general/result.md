# D10c — Corner asymptotics at general angle: M1(α) in closed form for all α, the strengthened rule proved at weight 2, and the weight-3 structure

Working directory: `research/08_theory/D10c_corner_general/` (scripts `step_12.py`–`step_33.py`, plus `corner_series.py`, `ir_utils.py`, `M_new_340.txt` reused from D10b).

## 0. Labels

- **PROVED** = verified exactly by sympy with `c = cos(α/2)` symbolic (printed simplification to 0).
- **PROVED (formal-log)** = sympy identity after applying formal log laws (`expand_log/logcombine(force=True)`); the algebraic core is exact, branch bookkeeping is then checked numerically.
- **NUMERICAL** = high-precision agreement with `corner_series.py` (fresh evaluations), residuals and digit counts in the evidence table.

---

## 1. The general-angle theorem (Part 1)

**Setup.** The corner function is $M(u;\alpha)=\int_0^1\!\big[(1+g)^{-1-u}\sqrt{1+g^{2-2u}-2\cos\alpha\, g^{1-u}}-1\big]g^{u-1}dg=\sum_k M_k(\alpha)u^k$. With $c=\cos(\alpha/2)$, $s=\sin(\alpha/2)$, $c_a=\cos\alpha=2c^2-1$, $r_1=-1+2c$, $r_2=-1-2c$, $P(t)=t^2-2c_at+1$, $D=(t^2-1)(t-c_a)(t-r_1)(t-r_2)$, the second Euler substitution $t=g+\sqrt{1-2c_a g+g^2}$ gives $g=(t^2-1)/(2(t-c_a))$, path $t:1\to t_1=1+2s$.

**(a) Substitution — PROVED** (`step_12.py` [1], all simplify to 0 with c free):
$1+g=\frac{(t-r_1)(t-r_2)}{2(t-c_a)}$, $S_0=\frac{P^2}{4(t-c_a)^2}$, $\frac{dg}{dt}=\frac{P}{2(t-c_a)^2}$, $1-c_ag=\frac{-c_at^2+2t-c_a}{2(t-c_a)}$, $c_ag-g^2=-\frac{(t^2-1)(t^2-2c_at+\cos2\alpha)}{4(t-c_a)^2}$, $\cos2\alpha=8c^4-8c^2+1$. Hence the whole integrand is rational over $\mathbb Q(c)$ — the surd $\sin\alpha$ never appears.

**(b) t-integrand — PROVED** (`step_13.py`, both coefficients simplify to 0):
$$M_1(\alpha)=\int_1^{t_1}\big[A^*(t)\log g+B^*(t)\log(1+g)\big]dt,\quad A^*=\frac{P^2}{D}-\frac{P}{(t^2-1)(t-c_a)}-\frac{t^2-2c_at+\cos2\alpha}{(t-c_a)(t-r_1)(t-r_2)},\quad B^*=-\frac{P^2}{D},$$
with $A^*=N/D$, $N=-(t-1)\big((t+1)^3-4c^2(t+1)-8c^4(t-1)\big)$ — the $(t-1)$ factor forces $\operatorname{res}_{t=1}A^*=0$.

**(c) Universal residues — PROVED** (`step_12.py` [3]: each residue minus claimed value simplifies to exactly 0; [4]: full partial fractions simplify to 0). At $\rho=(1,-1,c_a,r_1,r_2)$:
$$a(\rho)=(0,-2,1,-c,+c),\qquad b(\rho)=(-1,1,-1,+2c,-2c),\qquad A^*=\sum_\rho \tfrac{a(\rho)}{t-\rho},\quad B^*=\sum_\rho \tfrac{b(\rho)}{t-\rho}.$$

**(d) Antiderivative — PROVED (formal-log)** (`step_14.py` [1,2]): with $I(t;\rho,\sigma)=-\operatorname{Li}_2\!\big(\frac{t-\rho}{\sigma-\rho}\big)+C(\rho,\sigma)\log(t-\rho)$, $I(t;\rho,\rho)=\tfrac12\log^2(t-\rho)$, the pair identity $dI/dt=\log(t-\sigma)/(t-\rho)$ holds because $(\sigma-t)(\rho-\sigma)=(t-\sigma)(\sigma-\rho)$ identically (sympy: ratio = 1); assembled $dF/dt=A^*\log g+B^*\log(1+g)$ simplifies to 0.

**(e) Endpoint — PROVED** (`step_14.py` [3]): the $\log(t-1)$ divergence at $t\to1^+$ has coefficient $\log\!\big(\frac{(2-2c)(2+2c)}{2(1-c_a)\cdot1}\big)=\log 1=0$; the $\rho=1$ row has limit 0.

**(f) The closed form.** Define the weight matrix $W(\rho,1)=W(\rho,-1)=a(\rho)$, $W(\rho,c_a)=-(a(\rho)+b(\rho))$, $W(\rho,r_{1,2})=b(\rho)$:
$$W=\begin{pmatrix}0&0&1&-1&-1\\-2&-2&1&1&1\\1&1&0&-1&-1\\-c&-c&-c&2c&2c\\c&c&c&-2c&-2c\end{pmatrix}\quad(\text{rows, cols in order }1,-1,\cos\alpha,-1+2c,-1-2c).$$
Then for **all** $\alpha\in(0,\pi)$:
$$\boxed{\;M_1(\alpha)=\sum_{\rho,\sigma}W(\rho,\sigma)\big[I(t_1;\rho,\sigma)-I(1^+;\rho,\sigma)\big]-\sum_\rho(a(\rho)+b(\rho))\log2\,\big[\log(t_1-\rho)-\log(1-\rho)\big]\;}$$
with the $\rho=1$ row at the lower endpoint regularized to 0; $\operatorname{Li}_2(1)=\zeta(2)$ occurs with total coefficient $\gamma=\sum_{\rho\ne1}a(\rho)=\mathbf{-1}$ (proved; numerically confirmed). **Branch conditions:** $c\in(0,1)$ keeps the five poles distinct (collisions only at $\alpha=0,\pi$); all poles are $<1<t$ on the path, so every $t-\rho>0$; $\operatorname{Li}_2$ is principal (cut $[1,\infty)$); with $C(\rho,\sigma)=\log|\rho-\sigma|-i\pi$ for $\rho<\sigma$ the total is real and equal to $M_1$ (imaginary part $\sim10^{-150}$, NUMERICAL); with the principal convention take $\operatorname{Re}$. Expanded, the formula has **31 Li₂ terms** at the explicit arguments (functions of $c,s$, from `step_18.py`):
$$\tfrac{-1}{s},\ \tfrac{-s}{1\mp c},\ 1+s,\ \tfrac{1+s}{c^2},\ \pm\tfrac{1+s}{c},\ \tfrac{1+s}{s},\ \tfrac{\pm s(1+s)}{c^2},\ \tfrac{\pm s(1+s)}{c(1\pm c)},\ \tfrac{1+s\mp c}{1\mp c},\ \tfrac{\pm(1+s\mp c)}{c},\ \tfrac{\pm(1+s\mp c)}{c(1\mp c)},\ \tfrac{\pm(1+s\mp c)}{2c},\ \text{and lower-endpoint}\ \{1,\tfrac{\pm1}{c},\tfrac{1}{c^2},\tfrac{\pm s^2}{c^2},\tfrac{\pm(1\pm c)}{2c},\ldots\}$$
plus 19 log-products, 3 diagonal $\tfrac12\log^2$ terms, 5 $\log2\cdot\log$ terms, and $(-1)\zeta(2)$.

**(g) Numeric verification — NUMERICAL** (`step_16/17/19/20.py`): vs `corner_series.py` at the 8 required angles (2π/3, π/2, π/3, 3π/4, π/6, 5π/6, π/12, 1.3 rad): **67–68 digits** (step_16), **~90 digits** (step_19), and at dps 150 **147.8–150.3 digits** including near-endpoint α = 0.05, 0.2, 2.8, 3.0 (step_20); imaginary part $\le 2.1\!\times\!10^{-150}$ with the $-i\pi$ convention.

---

## 2. The strengthened rule (Part 2)

**Why all coefficients lie in $\mathbb Q\oplus\mathbb Q c$ — mechanism (PROVED from (c)–(f)).** At weight 2 the antiderivative assigns to each $(\rho,\sigma)$ pair **exactly one universal residue** as coefficient (diagonal pairs: residue/2; the $-i\pi$ pieces cancel globally). The residues lie in $\{0,\pm1,\pm2,\pm c,\pm2c\}\subset\mathbb Q\oplus\mathbb Qc$. Reduction of the Li₂ arguments to $(0,1]$ (reflection/inversion/Landen) only *permutes* coefficients among terms and appends ζ(2)- and log-product terms carrying the *same* single coefficient — at weight 2 no two coefficients are ever multiplied (that first happens at weight 3). Hence every coefficient of every Li₂ term, log-product, log² term and of ζ(2) lies in $\mathbb Z[\tfrac12]\{1,c\}=\mathbb Q\oplus\mathbb Qc$. *Parity remark:* under the involution $\tau:c\mapsto-c$ of the function field, $r_1\leftrightarrow r_2$ while $1,-1,c_a$ are fixed, and the residues are τ-covariant ($a(r_1)=-c\mapsto+c=a(r_2)$ etc.) — so the coefficient space is spanned by the τ-even $1$ and τ-odd $c$.

**Degrees (PROVED, sympy minpoly, `step_21.py`):** $K=\mathbb Q(\cos(\alpha/2))$: 5π/6 → degree **4** ($16x^4-16x^2+1$); π/12 → degree **8** ($256x^8-512x^6+320x^4-64x^2+1$); (also 3π/4 → 4, $8x^4-8x^2+1$; π/6 → 4).

**Confirmation (`step_22/23/24.py`, NUMERICAL):** reduced closed forms (canonical Li₂ reps in $(0,\tfrac12]$) at 5π/6 and π/12 agree with **fresh** `corner_coeffs` M1 to **220.1 / 219.3 digits**; *every* Li₂ coefficient, the ζ(2) coefficient, and all 12 distinct log-product coefficients are LLL-recognized in $\mathbb Q\oplus\mathbb Qc$ with maxcoeff ≤ 2 — in fact only $\{\pm\tfrac12,\pm1,\pm2\}\cup\{\pm\tfrac c2,\pm c,\pm2c\}$ occur, exactly the predicted span.

---

## 3. Weight 3 at the octagon (Part 3)

**(a) M2 integrand — PROVED** (`step_25/27/28.py`): $M_2(\alpha)=\int_1^{t_1}[R_{A_2}\log^2 g+R_{B_2}\log g\log(1+g)+R_{C_2}\log^2(1+g)]\,dt$ (no single-log terms), with $R_{C_2}=\frac{P^2}{2D}=-\frac{B^*}{2}$, $R_{B_2}=\frac{2(2c^2(t^2+1)-(t+1)^2)}{(t^2-1)(t-r_1)(t-r_2)}$, $R_{A_2}$ rational with denominator $2(t+1)(t-r_1)(t-r_2)(t-c_a)P^2$. **NUMERICAL** (`step_26.py`): t-integral vs `corner_series` M2 at 3π/4, π/2, 1.3: **75–76 digits**.

**(b) Residues — PROVED (exact, identity checks = 0):**
- $R_{B_2}$: $(-1,+1,+c,-c)$ at $(1,-1,r_1,r_2)$ — universal, in $\mathbb Q\oplus\mathbb Qc$; no pole at $c_a$.
- $R_{C_2}$: $(\tfrac12,-\tfrac12,-c,+c,\tfrac12)$ at $(1,-1,r_1,r_2,c_a)$ — universal, in $\mathbb Q\oplus\mathbb Qc$.
- $R_{A_2}$: $(0,-1,-\frac{c^2+1}{4c},+\frac{c^2+1}{4c},\frac12)$ at $(1,-1,r_1,r_2,c_a)$, **plus a quadratic-factor part** $\frac{s^2}{P}+\frac{4c^2s^2(t-1)}{P^2}$.

**(c) Structure — the letter set grows at weight 3:**
1. **Complex letters appear:** $R_{A_2}$ has *double poles at the roots $e^{\pm i\alpha}$ of $P$* — the first non-real letters (at weight 2, $P$ entered numerators only). Since $\frac1P$ splits over $e^{\pm i\alpha}$, the log² primitives force $\operatorname{Li}_3$ at $\frac{t-e^{\pm i\alpha}}{\sigma-e^{\pm i\alpha}}$; conjugate pairs recombine into real Li₃'s at special angles (at $\alpha=\pi/2$: $\operatorname{Li}_3(ib)+\operatorname{Li}_3(-ib)=\tfrac14\operatorname{Li}_3(-b^2)$, consistent with the squircle M₂'s real letters $\{b,\pm,1/\sqrt2\}$).
2. **The strengthened rule breaks at the residue level:** $\operatorname{res}_{r_{1,2}}R_{A_2}=\mp\frac{c^2+1}{4c}=\mp\tfrac14(c+c^{-1})\notin\mathbb Q\oplus\mathbb Qc$. Whether cancellations restore $\mathbb Q\oplus\mathbb Qc$ in the *total* M₂ is open.
3. **ζ(3) is forced** (PROVED for the log² parts): the primitives $G(t;\rho,\sigma)=L^2\log z-2L\operatorname{Li}_2(z)+\log z\log^2(1-z)+2\log(1-z)\operatorname{Li}_2(1-z)-2\operatorname{Li}_3(1-z)$ ($L=\log(\sigma-\rho)$, $z=\frac{t-\rho}{\sigma-\rho}$) give $\operatorname{Li}_3(1)=\zeta(3)$ at the lower endpoint whenever $\rho=1$: $R_{C_2}$ contributes $-2\,\mathrm{res}_{1}R_{C_2}\cdot(1+1-1)\zeta(3)=-\zeta(3)$ to $F(1^+)$ (i.e. $+\zeta(3)$ to $M_2$); $R_{B_2}$'s diagonal $(c_a,c_a)$ pair contributes $-2\zeta(3)$ to $M_2$; $R_{A_2}$ is regular at $t=1$ (contributes 0); the cross pairs $(\sigma,\tau)$, $\sigma\neq\tau$, contribute via the ladder $H(z;\xi)=\int\frac{\operatorname{Li}_2(z)}{z-\xi}dz$ (letters $\xi=\frac{\tau-\rho}{\sigma-\rho}$, cross-ratios inside $\mathbb Q(c)$) whose ζ(3)-content is ξ-independent rational — each contribution carries a **single** $R_{B_2}/R_{C_2}$ residue $\in\{\pm1,\pm\tfrac12,\pm c\}$, so **the ζ(3) coefficient is expected in $\mathbb Q\oplus\mathbb Qc$** (PROVED for the non-cross parts; CONJECTURE for the total pending the completed ladder; the $\rho=\tau$ pieces already cancel: $+c-c=0$).
4. **Ladder status** (`step_29–33.py`, small LLL — allowed): $H(z;\xi)$ closes on $\operatorname{Li}_3$ at $\{\frac z\xi,\frac{\xi-z}\xi,\frac z{z-\xi},\frac\xi{\xi-z}\}$ plus log-products and $\zeta(3)$ (two independent candidate relations with maxcoeff ≤ 6, residuals ~$10^{-102}$); clean extraction/verification of the H-relation was **not finished** — the full M₂(3π/4) closed form vs the 330-digit value remains open. No large-basis PSLQ was attempted.

---

## 4. Evidence table

| check | label | residual / digits |
|---|---|---|
| step_12 [1]: 6 substitution identities, c free | PROVED | all simplify to 0 |
| step_13: t-integrand = A* log g + B* log(1+g) | PROVED | X−A* = 0, Y−B* = 0 |
| step_12 [3]: 10 universal residues (5 poles × 2) | PROVED | res − claimed = 0 exactly |
| step_12 [4]: partial fractions A*, B* | PROVED | both = 0 |
| step_14 [1]: pair antiderivative identity | PROVED (formal-log) | ratio = 1 |
| step_14 [2]: dF/dt = integrand | PROVED (formal-log) | 0 |
| step_14 [3]: endpoint log(t−1) coefficient | PROVED | log 1 = 0 |
| step_16/17: compact formula vs series, 8 angles | NUMERICAL | 2.2e-68…4.6e-67 (67–68 dig); Im ≤ 2e-90 (D10b conv.); ζ(2) coeff = −1.0 |
| step_19: expanded formula (31 Li₂), 8 angles | NUMERICAL | 89.5–90.4 dig; Im ≤ 3e-90 |
| step_20: dps 150, 8 angles incl. 0.05, 0.2, 2.8, 3.0 | NUMERICAL | 147.8–150.3 dig; Im ≤ 2.1e-150 |
| step_21: degrees of K at 5π/6, π/12 | PROVED | minpolies deg 4, 8 |
| step_22/24: reduced forms 5π/6, π/12 vs fresh M1 | NUMERICAL | 220.1 / 219.3 dig; all coeffs in ℚ⊕ℚc, maxcoeff ≤ 2 |
| step_26: M2 t-integral vs series M2, 3 angles | NUMERICAL | 75–76 dig |
| step_27/28: RB2, RC2, RA2 partial fractions | PROVED | identity checks = 0 |
| step_29–33: ladder H candidate relations | NUMERICAL | maxcoeff ≤ 6, resid ~1e-101; extraction incomplete |

## 5. Files
`step_12.py` (substitution + integrand + residues + PF, PROVED), `step_13.py` (integrand, PROVED), `step_14.py` (antiderivative + endpoint, PROVED formal-log), `step_16/17.py` (compact form, conventions, γ=−1), `step_18/19.py` (expanded form + argument table), `step_20.py` (dps-150 robustness), `step_21.py` (degrees), `step_22/23/24.py` (strengthened rule at 5π/6, π/12), `step_25.py` (+`step_25_out.json`, M2 integrand), `step_26.py` (M2 numeric check), `step_27/28.py` (M2 partial fractions, PROVED), `step_29–33.py` (ladder attempt); `corner_series.py`, `ir_utils.py`, `M_new_340.txt` from D10b.

## 6. Open questions
1. Finish the cross-term ladder $H(z;\xi)$ symbolically and assemble the full M₂(3π/4) closed form; test against the 330-digit value.
2. Exact ζ(3) coefficient of M₂(α): identified contributions $\{+1\ (R_{C_2}),\ -2\ (R_{B_2}\text{ diag}),\ \text{cross via }H\}$, each $\in\mathbb Q\oplus\mathbb Qc$ — total membership open.
3. Do the $\mp(c+c^{-1})/4$ RA2-residues cancel in the total M₂ (rule restored) or does the D10b strengthening genuinely fail at weight 3?
4. Mechanism by which conjugate-letter pairs $e^{\pm i\alpha}$ collapse to real classical Li₃'s at cyclotomic angles (π/2 must reproduce the squircle's real-letter M₂).
5. Symbolic certificate for the $-i\pi$ branch convention making $F$ exactly real (verified to 150 digits at 16 angle-checks).