# verify_num: independent numerics for D10c Parts 1 and 2

**Verdict: confirmed, with corrections to how it is stated** (see §4). Part 1's boxed $M_1(\alpha)$ and D10c's evaluator `M1_expanded` agree with **rigorous** arb enclosures of $M_1$ at 10 angles (radius $2.5\times10^{-113}$, actual agreement about $10^{-131}$). At 4 angles the enclosures have 229 digits, and the agreement there is $10^{-240}$. The method is independent: values of $M(u;\alpha)$ at complex $u$, then a Cauchy/FFT sum. It uses no $u$-series and no $d/du$. A separate derivation by a different substitution gives the same exact symbol as the D10c formula. Part 2 holds at the new angles $2\pi/5$ (degree 2) and $\pi/5$ (degree 4). Every coefficient of my reduced closed form is an exact element of $\mathbb Q+\mathbb Qc$, and the reduced form matches 229-digit rigorous $M_1$ to $\le 6.6\times10^{-246}$.

All scripts and outputs are in `<scratch>/d10c_num/` (scratch; listed in §5). The only file written in the project is this one.

---

## 0. Method for (1): rigorous Cauchy/FFT extraction

With $g=e^{-x}$ and $w=e^{-(1-u)x}$,
$$M(u;\alpha)=\int_0^\infty \operatorname{expm1}\!\Big(-(1+u)\operatorname{log1p}(e^{-x})+\tfrac12\big[\operatorname{log1p}(-we^{i\alpha})+\operatorname{log1p}(-we^{-i\alpha})\big]\Big)\,e^{-ux}\,dx .$$
For $|w|<1$ both factors $1-we^{\pm i\alpha}$ have positive real part. So $\sqrt{1-2\cos\alpha\,w+w^2}=\sqrt{1-we^{i\alpha}}\sqrt{1-we^{-i\alpha}}$ with principal branches is the analytic continuation from real $u$. Using expm1/log1p avoids the catastrophic cancellation as $x\to\infty$.

- **Integration:** `acb.integral` (python-flint 0.9.0, arb) on $[0,X]$. The analytic flag is honoured: the integrand returns NaN unless $\operatorname{Re}(1+y)>0$ is certain for every log1p argument $y$.
- **Tail on $[X,\infty)$, for $|u|=r$:** $|\text{integrand}|\le 4(1+r)e^{-(1-r)x}+4e^{-(1-2r)x}$. This follows from $|\operatorname{log1p}y|\le2|y|$ for $|y|\le\tfrac12$ and $|\operatorname{expm1}E|\le2|E|$ for $|E|\le\log2$. The bound is added to each value as a ball radius.
- **Aliasing bound:** $\frac1{Nr}\sum_j M(r\omega^j)\omega^{-j}=M_1+\sum_{m\ge1}M_{1+mN}r^{mN}$. By Cauchy's estimate on $|u|=R=\tfrac12$, the error is at most $\frac{B}{R}\frac{q^N}{1-q^N}$ with $q=r/R$.
- **Hand bound $B$:** for $|u|\le\tfrac12$:
  - On $[0,2\log2]$ the integrand is at most $3e^{x/2}\le6$, contributing at most 8.32.
  - On $[2\log2,\infty)$ it is at most $5.76\,(3e^{-x/2}+2e^{-x})$, contributing at most 20.1.
  - So $|M|\le28.5$, and I used $B=50$. A sampled check gave $\max|M|=1.65$ on $|u|=\tfrac12$.
- **Symmetry:** $M(\bar u)=\overline{M(u)}$ halves the work.

**Configurations:**
- A: $r=1/32$, $N=96$, 440 bits.
- B: $r=1/16$, $N=128$, 440 bits.
- Both give aliasing radius $2.5\times10^{-114}$. A and B are independent quadrature node sets, and they agree to about $10^{-131}$.
- HI: $r=1/32$, $N=192$, 820 bits, $X=600$, radius $6.4\times10^{-230}$; run at $2\pi/5$, $\pi/5$, $5\pi/6$, $1.3$.

Runtime is about 3 s per angle (HI: about 35 s).

## 1. $M_1(\alpha)$: Cauchy/FFT versus `step_19.M1_expanded`

`M1_expanded` was evaluated at dps 140. The rigorous radius for every angle in the table is $2.5\times10^{-114}$.

| α | $M_1$ (Cauchy, 60 digits) | \|A−B\| | M1_expanded − Cauchy |
|---|---|---|---|
| π/2 | 0.252619096207044794618363520965654860533627150724382125576718 | 1.4e-131 | −7.6e-132 |
| 2π/3 | 0.125795678937973233753868197683859494817894210686283298857014 | 3.6e-132 | −2.5e-132 |
| 5π/6 | 0.0335548827346891916144662392772242926156842529212447204735174 | 1.1e-132 | 1.4e-132 |
| π/12 | 0.478809166914067907859335435642066402377323165111799706578438 | 7.4e-132 | −3.4e-133 |
| 0.7 | 0.440557223084364221529397307070738528367876979234506326820079 | 8.8e-133 | 1.6e-132 |
| 2.5 | 0.0498476854483588892982044872228417570715521412083192701525906 | 1.8e-133 | 1.8e-132 |
| *2π/5* | 0.330642100401566231485909044010145969676876803771592253336617 | 1.8e-131 | −7.8e-132 |
| *π/5* | 0.450199065302689521030793735254681650671404440143937580233009 | 1.8e-131 | −1.0e-131 |
| *0.05* | 0.480699710155960639256153815833853526233206220347470881464511 | 1.6e-131 | −1.1e-131 |
| *3.0* | 0.00250222317317432384821017382095704722005568660037595158359571 | 3.8e-133 | 2.1e-132 |

`M1_expanded` lies inside the rigorous ball at all 10 angles, which is at least 113 certified digits. The actual agreement is about 131 digits; the Cauchy midpoints are good to their 440-bit working precision. Im(M1_expanded) is at most $2.5\times10^{-140}$.

At the HI angles, D10c's boxed formula agrees with the 229-digit balls to about $2\times10^{-241}$ (`extra_checks.py`). Here I evaluated it with **real parts only** ($\operatorname{Re}\operatorname{Li}_2$, $\log|\rho-\sigma|$, no $-i\pi$), with ζ(2) coefficient $-1$.

## 2. α = π/2 and the squircle: an exact identity

| k | corner $M_k(\pi/2)$ (Cauchy) | squircle $M_k$ (mk_deg9.tsv) | diff |
|---|---|---|---|
| 0 | −0.553303299720515717370808039043… | −0.553303299720515717370808039043… | 4.7e-81 |
| 1 | 0.252619096207044794618363520966… | 0.252619096207044794618363520966… | 2.1e-81 |
| 2 | −0.0221619110872693994237549377398… | same | −3.1e-82 |
| 3 | 0.0227705214917865774567305097423… | same | 2.9e-82 |

The differences equal the TSV's own radii (about 80 digits), so the agreement is exact to the precision available.

**Normalisation link (proved here).** The squircle is $P(p)/8=1+M_{\rm sq}(q)$ with $q=e=1/p$ and
$$M_{\rm sq}(q)=-\int_0^{1/2}\big[t^{q-1}+(1-t)^{q-1}-\sqrt{t^{2q-2}+(1-t)^{2q-2}}\big]dt .$$
Substitute $t=g/(1+g)$, so that $dt=dg/(1+g)^2$ and $1-t=1/(1+g)$:
$$M_{\rm sq}(u)=-\int_0^1(1+g)^{-1-u}\big[g^{u-1}+1-\sqrt{1+g^{2u-2}}\big]dg .$$
Since $g^{u-1}\sqrt{1+g^{2-2u}}=\sqrt{1+g^{2u-2}}$, we have, for $0<u<1$,
$$M(u;\tfrac\pi2)-M_{\rm sq}(u)=\int_0^1\big[(1+g)^{-1-u}(1+g^{u-1})-g^{u-1}\big]dg=\frac{1-2^{-u}}u+\frac{2^{-u}}u-\frac1u=0 .$$
The middle term uses $g\mapsto1/g$: $\int_0^1(1+g)^{-1-u}g^{u-1}dg=\int_1^\infty(1+h)^{-1-u}dh$.

So $M(u;\pi/2)\equiv M_{\rm sq}(u)$ as analytic functions, with $u=e=1/p$ and no rescaling. Hence **$M_k(\pi/2)=M_k^{\rm sq}$ for every $k$**. The corner's "$-1$" is the squircle's "$1+$" in $P/8$: one unit of straight side per octant.

## 3. Q + Qc at new angles, by an independent derivation

**Field degrees (sympy minpoly).**
- α = 2π/5: $c=\cos(\pi/5)$ has minimal polynomial $4x^2-2x-1$, **degree 2**. So $\mathbb Q+\mathbb Qc=\mathbb Q(c)$ and the test is weak; it only checks that $s=\sin(\pi/5)\notin\mathbb Q(c)$ never enters a coefficient.
- α = π/5: $c=\cos(\pi/10)$ has minimal polynomial $16x^4-20x^2+5$, **degree 4**. This is a genuine test.

**3a. My own route** (`own_closed_form.py`). It differs from D10c's in two ways:
- An Euler substitution of the *first* kind, $\sqrt{S_0}=1+\tau g$. This gives $g=2(\tau+c_a)/(1-\tau^2)$ with $\tau$ running from $-c_a$ to $2s-1$.
- The rational parametrisation $v=\tan(\alpha/4)$, with $c=\frac{1-v^2}{1+v^2}$ and $s=\frac{2v}{1+v^2}$. Everything becomes rational in $(\tau,v)$, and every identity is an exact sympy `cancel`. The checks are: $\sqrt{S_0}^2=S_0$, the endpoints, both partial fractions, and the factorisations of $g$ and $1+g$. All of them are 0.

Simple poles only, at $\tau\in\{1,-1,-c_a,1+2c,1-2c\}$ ($c_a=\cos\alpha$, the lower endpoint is $\tau=-c_a$). Residues:

| pole τ | 1 | −1 | −c_a | 1+2c | 1−2c |
|---|---|---|---|---|---|
| coefficient of log g | 1 | 1 | 0 | +c | −c |
| coefficient of log(1+g) | 1 | −1 | −1 | −2c | +2c |

All residues are exactly $p+qc$. They correspond to D10c's universal residues under the Möbius map $\tau=(t-1-2c_a)/(t+1)$, which sends $t\in\{1,-1,c_a,r_1,r_2,\infty\}$ to $\tau\in\{-c_a,\infty,-1,1-2c,1+2c,1\}$. The residue at ∞ is minus the sum of the others. So the residue set is confirmed, but it is the same 1-form in other coordinates, not a structurally new mechanism.

The proof uses only real variables. The antiderivative is
$$J_{\rho\sigma}=\log|\sigma-\rho|\log|\tau-\rho|-\operatorname{Re}\operatorname{Li}_2\!\big(\tfrac{\tau-\rho}{\sigma-\rho}\big),\qquad J_{\rho\rho}=\tfrac12\log^2|\tau-\rho|,$$
which needs only $\frac{d}{dx}\operatorname{Re}\operatorname{Li}_2(x)=-\log|1-x|/x$. On the path $-1<1-2c<-c_a\le\tau\le2s-1<1<1+2c$ for every $c\in(0,1)$. At $\tau\to-c_a$ the coefficients of $\log(\tau+c_a)$ and of $\log^2(\tau+c_a)$ vanish exactly. This gives a complete proof of a closed form ($M_1=F(\tau_1)-F(\tau_0^+)$: 30 ReLi₂ terms, 48 log products, $2\zeta(2)$) with no formal logs. It matches Cauchy to about $10^{-131}$ at all 10 angles of §1 (`own_reduce.py`).

**Exact cross-check against D10c** (`symbol.py`). Both closed forms have **exactly the same weight-2 symbol** (sympy). It is $\sum\kappa_{\ell m}\,\ell\otimes m$ over the letters $\{2,v,1-v,1+v,1+v^2\}$:

| ⊗ | 2 | v | 1−v | 1+v | 1+v² |
|---|---|---|---|---|---|
| 2 | 2 | 2−2c | 2 | 2 | −4 |
| v | 2−2c | 0 | 2−2c | 2−2c | −4 |
| 1−v, 1+v | 0 | 0 | 0 | 0 | 0 |
| 1+v² | −2 | −2+2c | −2 | −2 | 4 |

Every entry is in $\mathbb Q+\mathbb Qc$. The antisymmetric part is nonzero and has a $c$-part: $v\wedge(1\mp v)$ and $v\wedge(1+v^2)$ carry $-c$. So the $c$-component genuinely needs Li₂ and cannot be absorbed into $\mathbb Q$. Li₂ functional equations do not change the symbol, so this form of Part 2 does not depend on which Li₂ arguments are chosen.

**3b. Reduced forms** (`own_reduce.reduce_at`, `run_part3.py`). Each ReLi₂ argument is reduced to $(0,\tfrac12]$ with these rules:
- reflection for $x\in(\tfrac12,1)$;
- Landen for $x<0$;
- inversion for $x>1$: $\operatorname{Re}\operatorname{Li}_2x=-\operatorname{Li}_2(1/x)+2\zeta(2)-\tfrac12\log^2x$.

Equal arguments are merged numerically at 250 digits. Log products are collected over the letters. Coefficients are tracked **exactly** as pairs $(p,q)$ meaning $p+qc$, so a coefficient outside $\mathbb Q+\mathbb Qc$ would show up as a failed exact solve. This is the "exact linear algebra" route. I did not run a separate PSLQ fit: with exact coefficients and a 245-digit match against a rigorous ball, it would only repeat the check.

| α | reduced form − rigorous $M_1$ (229-digit ball) | ReLi₂ terms | log-pair terms | ζ(2) coefficient | coefficient outside Q+Qc? |
|---|---|---|---|---|---|
| **2π/5** (deg 2) | **4.6e-246** | 11 | 13 | −2 | none |
| **π/5** (deg 4) | **6.6e-246** | 11 | 13 | 0 | **none** |
| 5π/6 (control, deg 4) | 1.3e-246 | 11 | 9 | 0 | none |
| 1.3 (generic) | 1.6e-246 | 11 | 13 | −1+c | none |

Li₂ part at **α = π/5**, with $v=\tan(\pi/20)$, $c=\cos(\pi/10)$:
$$4c\,\mathrm{Li}_2\tfrac{v^2}{1+v^2}-(2+2c)\mathrm{Li}_2\tfrac{2v^2}{1+v^2}+2\,\mathrm{Li}_2\tfrac{4v^2}{(1+v^2)^2}+(1+c)\mathrm{Li}_2\tfrac{v(1-v)}{1+v}-4c\,\mathrm{Li}_2\tfrac{v}{1+v}+2c\,\mathrm{Li}_2v$$
$$+(1-c)\mathrm{Li}_2\tfrac{v(1+v)}{1+v^2}+(1+c)\mathrm{Li}_2\tfrac{2v}{1+v}-2\,\mathrm{Li}_2\tfrac{2v}{1+v^2}+(1-c)\mathrm{Li}_2\tfrac{1-v}2+(2c-2)\mathrm{Li}_2\tfrac{1-v^2}2 .$$
The arguments are 0.02447, 0.04894, 0.09549, 0.11507, 0.13673, 0.15838, 0.17898, 0.27346, 0.30902, 0.42081 and 0.48746.

Log part at π/5, as products of letters:

| letter pair | coefficient |
|---|---|
| (1+v)² | 1−c |
| (1+v)(1+v²) | −6+2c |
| (1+v)(1−v) | −2c |
| (1+v)·2 | 4−2c |
| (1+v)·v | 2−2c |
| (1+v²)² | 9/2+c/2 |
| (1+v²)(1−v) | −3+3c |
| (1+v²)·2 | −4 |
| (1+v²)·v | −4 |
| (1−v)·2 | 3−c |
| (1−v)·v | 2−2c |
| 2·2 | 1/2+c/2 |
| 2·v | 2−2c |

The ζ(2) coefficient is 0.

At 2π/5 the Li₂ coefficient pattern is the same as a function of $v$. One difference: 1−s = 0.41221 < ½ stays as $\mathrm{Li}_2\frac{(1-v)^2}{1+v^2}$ with coefficient +2, where π/5 has $-2\,\mathrm{Li}_2(s)$. The ζ(2) coefficient there is −2. Full lists are in `run_part3.out` and `part3_report.json`.

**Conclusion for (3):** at the new degree-4 angle π/5, every Li₂, log-product and ζ(2) coefficient lies in $\mathbb Q+\mathbb Qc$. At 2π/5 every coefficient lies in $\mathbb Q(c)$; there, $\mathbb Q+\mathbb Qc$ is the same set.

## 4. Corrections and gaps for quoting D10c as a theorem

1. **Drop the formal-log / $-i\pi$ caveat. It is unnecessary.** For real $x\ne0,1$, $\frac{d}{dx}\operatorname{Re}\operatorname{Li}_2x=-\log|1-x|/x$. So the real-part version of the boxed formula ($\operatorname{Re}\operatorname{Li}_2$, $C=\log|\rho-\sigma|$) is a genuine real-variable antiderivative, because no argument crosses 0 or 1 along the path. Numerically it equals $M_1$ to $10^{-240}$. Consequences:
   - D10c 1(d) can be upgraded from "PROVED (formal-log)" to PROVED.
   - D10c open question 5 (certifying that the imaginary part vanishes) goes away.
   - The theorem should be stated with Re (or with the real reflection and inversion forms), not with the $-i\pi$ convention.
2. **Part 2 is an existence statement and must be worded as one.** The coefficient list depends on the representation. D10c reports "only $\{\pm\frac12,\pm1,\pm2\}\cup\{\pm\frac c2,\pm c,\pm2c\}$ occur, maxcoeff ≤ 2". That holds only for D10c's own unmerged t-route bookkeeping. In my equally valid reduced form, coefficients include $-2-2c$, $\pm4c$, $\tfrac{13}2+\tfrac c2$ and $-7+3c$. Do not put "maxcoeff ≤ 2" or "pure rational or pure c-multiple" into the theorem. The robust formulations are:
   - (a) "there is a representation with all coefficients in $\mathbb Z[\tfrac12]+\mathbb Z[\tfrac12]c$"; or
   - (b) the intrinsic one: "the weight-2 symbol of $M_1$ has all entries in $\mathbb Q+\mathbb Qc$", as in the table in §3a, verified exactly.
3. **A degree-2 angle such as 2π/5 cannot test Q + Qc against Q(c).** D10c's 5π/6 (degree 4) and π/12 (degree 8), and my π/5 (degree 4), are the real tests.
4. **Rigour of my numerics.** The Cauchy values are rigorous balls, given arb's correctness and the hand bounds in §0 (tail bound, $|M|\le28.5$ on $|u|\le\frac12$). The closed-form identities are exact rational-function identities in sympy.
5. **Note on D10c 1(g) (not an error).** The claim of "147.8–150.3 digits" at dps 150 used `corner_series` as the reference. My reference is independent and rigorous, and it gives at least 113 certified digits at 10 angles and at least 229 at 4.

## 5. Files (scratch, `<scratch>/d10c_num/`)

| file | contents |
|---|---|
| `cauchy_m.py` | arb integrand, tail and aliasing bounds, `Mk_cauchy` |
| `run_part1.py` → `part1_cauchy.json` | 10 angles × configs A, B; $k=0..3$ |
| `compare_part1.py` → `compare_part1.out` | §1 table, squircle comparison |
| `run_cauchy_hi.py` → `cauchy_hi.json` | 229-digit balls at 2π/5, π/5, 5π/6, 1.3 |
| `own_closed_form.py` → `own_cf.pkl`, `own_closed_form.out` | τ-route exact derivation, residues |
| `own_reduce.py` → `own_reduce.out` | assembly, numeric check, exact reduction to (0,½] |
| `run_part3.py` → `run_part3.out`, `part3_report.json` | §3b tables |
| `symbol.py` → `symbol.out` | exact symbol, own vs D10c |
| `extra_checks.py` → `extra_checks.out` | real-part D10c formula; sampled check of the bound on $\lvert M\rvert$ |
