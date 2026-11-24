# D13b. Exact constant in the large‑k asymptotics of the squircle coefficients $M_k$

## 1. Headline result

$$\boxed{\;M_k \;=\; -\,\frac{\sqrt{\pi}}{k}\,e^{-\sqrt{\pi k}}\Bigl[\cos\!\sqrt{\pi k}\;+\;O\!\left(k^{-1/2}\right)\Bigr]\;}$$

i.e. in $M_k \sim C k^{\beta} e^{-A\sqrt k}\cos(B\sqrt k+\phi_0)$: $A=B=\sqrt\pi$, $\beta=-1$, **$C=\sqrt\pi$** (amplitude constant $A_0=C/\sqrt\pi=1$), **$\phi_0=\pi$**. The planner is right: the previous agent's $C=\sqrt\pi/2$ was off by exactly the factor 2 = sum over the **two conjugate branch points** $z_\pm=\pm i\pi/2$. A sharper form (with $N=k+1$, $\theta=\sqrt{\pi N}$):

$$M_k=-\frac{\sqrt\pi}{N}\,e^{-\theta}\Bigl[\cos\theta\bigl(1+\tfrac{c_r}{\sqrt N}+\dots\bigr)+\tfrac{c_i}{\sqrt N}\sin\theta+\dots\Bigr],\quad c_r\approx0.84\text{–}0.86,\ c_i\approx0.045\text{–}0.050 .$$

---

## 2. Step 1 — $M_k$, $k\le 400$, certified $\ge 12$ digits (NUMERICAL, achieved with 52 digits)

Method: previous agent's stable representation $M(q)=-\int_0^\infty e^{-qL(u)}B(u(q-1))/(e^u+1)\,du$, $L=\log(1+e^{-u})$, evaluated as FLINT `acb_series` in $q$ at each tanh–sinh node in $u$. Parameters used (as instructed): **udeg = 8 and 9, mdps = 160, bits = 800, umax = 1500** ($> k+40\sqrt k=1200$). Printed runtimes: 207 s (1809 nodes) and 704 s (3618 nodes).

Printed control checks: $M_0=-0.553303299720515717370808039042959437225890490$ (45 digits correct); $M_{50},M_{100},M_{200},M_{250}$ match all previous control digits.

**Printed table (deg 9 value; last column = significant digits on which deg 8 and deg 9 agree):**

| $k$ | $M_k$ (degree 9) | deg8=deg9 sig digits |
|---:|---|---:|
| 50 | `-1.257059900405641281369e-7` | 100 (capped) |
| 100 | `-1.772316807359054100089e-10` | 100 |
| 150 | `4.281902399189632130486e-12` | 100 |
| 200 | `-1.145441719208160174755e-13` | 100 |
| 250 | `4.66656107539387895708e-15` | 76.6 |
| 300 | `-2.146317091805711938025e-16` | 66.3 |
| 350 | `4.318804365699294281484e-18` | 58.2 |
| 400 | `1.057861107585034098121e-18` | 52.4 |

Worst agreement over **all** $k\in[1,400]$: **52.375 significant digits** (at $k=400$). Requirement ($\ge12$) exceeded by 40 digits. Full arrays in `mk_deg8.tsv`, `mk_deg9.tsv` (produced by `step_01.py`, compared by `step_02.py`).

## 3. Step 2 — Fits (NUMERICAL)

With $R_k := M_k\,k\,e^{\sqrt{\pi k}}/\sqrt\pi = A_k\cos(\sqrt{\pi k}+\phi_k)$, the two-consecutive-$k$ phasor estimator (printed, `step_03b.py`):

| $k$ | $A_k$ | $\phi_k/\pi$ |
|---:|---:|---:|
| 50 | 0.98731 | 1.03813 |
| 100 | 0.98742 | 1.02654 |
| 200 | 0.99447 | 1.01847 |
| 300 | 0.99534 | 1.01543 |
| 398 | 0.99579 | 1.01305 |

Global least-squares fits $R_k=\mathrm{Re}[e^{i\psi_k}\sum_{j=0}^{J}c_j k^{-j/2}]$ (`step_04.py`), $A_0=|c_0|,\ \phi_0=\arg c_0$: across 12 windows ($J=2,3,4$; $k_{\min}=60..200$): $A_0\in[0.9940,1.0056]$, central $\approx1.000$; $\phi_0\in[3.1300,3.1638]$, with the best-conditioned windows ($J=2$) giving $\phi_0\in[3.13999,3.14689]=\pi(0.99949\text{–}1.00169)$.

$$\boxed{A_0 = 1.000\pm0.003\;(\text{spread }\pm0.005),\qquad \phi_0 = \pi\;(1.000\pm0.002)\quad\text{NUMERICAL}}$$

consistent with exactly $A_0=1,\ \phi_0=\pi$ (planner), excluding the old $\sqrt\pi/2$ claim decisively (it would give $A_0=0.5$).

First corrections in the fit convention ($A_k=A_0+A_1/\sqrt k+\dots$, $\phi_k=\phi_0+\phi_1/\sqrt k+\dots$): $A_1\approx-0.05\pm0.05$, $\phi_1\approx0.85\pm0.06$ (J=2 windows: $A_1\in[-0.072,-0.004]$, $\phi_1\in[0.845,0.898]$).

## 4. Step 3 — Saddle-point derivation with both branch points

### 4.1 Mechanism and the exact local constant (PROVED — symbolically exact, verified numerically)

Writing $M(q)=-\int_0^\infty e^{-qL(u)}B(u(q-1))/(e^u+1)\,du$ with $B(z)=1-\dfrac{e^z}{1+\sqrt{1+e^{2z}}}$: for $u\to\infty$ the weight reduces to $e^{-u}(1+O(e^{-u}))$ and $e^{-qL}=1+O(e^{-u})$, so beyond-all-orders-exactly $M_k\approx-\int_0^\infty e^{-u}[q^k]B(u(q-1))\,du$. $[q^k]B(u(q-1))=\frac1{2\pi i}\oint B(z)(1+z/u)^{-N}dz/u$ ($N=k+1$) is dominated by the nearest branch points $z_s=s\,i\pi/2$, $s=\pm1$. Sympy series (`step_07.py`), confirmed numerically (`step_06.py`, $10^{-40}$-offset discontinuity):

$$B(z_s+\zeta)=(1-si)+s\sqrt2\,\zeta^{1/2}+i\zeta-\tfrac{s\sqrt2}{2}\zeta^{3/2}-\tfrac{i}{2}\zeta^2+\tfrac{5s\sqrt2}{24}\zeta^{5/2}+\dots$$

**The claimed local constant $s\sqrt2\,\zeta^{1/2}$ is CORRECT** (above-cut branch; sign flips below). Discontinuity $\mathrm{Disc}=2s\sqrt2\,\zeta^{1/2}(1-\tfrac12\zeta+\tfrac5{24}\zeta^2-\dots)$, i.e. $c_3=-1/2$ (PROVED).

### 4.2 Saddle (PROVED — sympy-exact, `step_08.py`, `step_11.py`)

Hankel transfer ($\zeta^{1/2}\to\Gamma(\tfrac32)\lambda^{-3/2}$, $\lambda=N/(u+z_s)$) then $u$-saddle of $\phi(u)=-u-N\log(1+z_s/u)$: exact saddle $u_s(u_s+z_s)=Nz_s$, and printed sympy series:
$\phi(u_s)=-2\sqrt{Nz_s}+\tfrac{z_s}{2}-\tfrac{z_s^{3/2}}{12\sqrt N}+\dots$, $\phi''(u_s)=-2/\sqrt{Nz_s}+\dots$, prefactor $P(u_s)=(u_s+z_s)^{3/2}/u_s=z_s^{1/4}N^{1/4}(1+\tfrac54\sqrt{z_s/N}+\dots)$, Gaussian-moment loop correction $=\dfrac{\varepsilon}{2\sqrt{z_s}}+O(\varepsilon^2)$, and $2\sqrt{Nz_s}=\sqrt{\pi N}(1+si)$.

### 4.3 The constant — factor 2 from the two branch points

Per branch point: $M_k^{(s)}\approx\kappa_s N^{-1}e^{-2\sqrt{Nz_s}}(1+\beta_s/\sqrt N)$ with
$$\kappa_s=\frac{is}{\sqrt2}\,z_s^{1/2}e^{z_s/2}=\frac{is}{\sqrt2}\cdot\sqrt{\tfrac\pi2}e^{si\pi/4}\cdot e^{si\pi/4}=-\frac{\sqrt\pi}{2}\quad(s=\pm1,\ \textbf{both equal, real, negative}).$$
Summing **both**: $M_k\approx-\dfrac{\sqrt\pi}{2N}\,2\,\mathrm{Re}\,e^{-\sqrt{\pi N}(1+i)}=-\dfrac{\sqrt\pi}{N}e^{-\sqrt{\pi N}}\cos\sqrt{\pi N}$. Hence $C=\sqrt\pi$, and the overall minus sign = $\phi_0=\pi$ in the $\sqrt{\pi k}$-convention. The first-order coefficient assembled from the verified pieces:
$$\beta_s=\frac{1}{2\sqrt{z_s}}+\frac{\sqrt{z_s}}{2}-\frac{z_s^{3/2}}{12},\qquad c_r=\mathrm{Re}\,\beta_+=\tfrac{1}{2\sqrt\pi}+\tfrac{\sqrt\pi}{4}+\tfrac{\pi^{3/2}}{48}=0.8412151,\quad c_i=\mathrm{Im}\,\beta_+=-\tfrac{1}{2\sqrt\pi}+\tfrac{\sqrt\pi}{4}-\tfrac{\pi^{3/2}}{48}=0.0450118 .$$

### 4.4 Comparison with numerics (evidence table)

| test | prediction | numerical result | status |
|---|---|---|---|
| $\gamma_s$ in $B(z_s+\zeta)$ | $s\sqrt2=1.41421356$ | `1.41414285454` (r=$10^{-4}$, both s, both sides) | PROVED+NUM ✓ |
| $c_3$ | $-1/2$ | sympy exact; Disc series vs numeric `0.61697354766` vs `0.61696866035` @r=0.05 | PROVED ✓ |
| ζ-integral asymptotics at saddle (`step_15.py`) | ratio 1 | `1.0046453 - 0.00033i` (O($\lambda^{-3}$) size) | NUM ✓ |
| leading formula vs $M_k$ (`step_10.py`) | $1+(c_r\cos+c_i\sin)/\sqrt N$ | ratios 1.1413(k=50)→1.0481(k=400), matching $1+0.85/\sqrt N$ trend within O(1/N) | NUM ✓ |
| residual with analytic $\beta$ (`step_12.py`) | $N\!\cdot\!R_k$ bounded O(1) | $N R_k\in[-0.96,1.17]$, envelope $\sim1/\sqrt N$ decay ✓ | NUM ✓ |
| $c_r,c_i$ free fit (`step_09/13.py`) | 0.8412, 0.0450 | $0.860\pm0.006,\ 0.050\pm0.003$ (stable over 12 windows) | ~2% tension, see below |
| $A_1,\phi_1$ fit-convention | $-0.0450,\ 0.8412$ | $-0.05\pm0.05,\ 0.85\pm0.06$ | NUM ✓ consistent |

The analytic $A_0=1,\phi_0=\pi$ matches step 2 exactly. The first-order coefficients $c_r,c_i$: analytic derivation gives (0.8412, 0.04501); the stable free fits give (0.860, 0.050); my independent saddle cross-checks confirmed the ζ-part and all sympy-verifiable pieces, but a direct steepest-path numeric integration of the u-saddle (`step_16/17.py`) contained integrand bugs and is excluded as inconclusive. The residual ~2% in $c_r$ (0.019, i.e. one missing O($N^{-1/2}$) sub-piece) is unresolved; the first-order analytic coefficients are therefore reported as **CONJECTURE**, the leading constant as derived + numerically confirmed.

## 5. Statements with labels

1. **NUMERICAL (52 digits):** $M_k$ for $k\le400$; degrees 8/9 of tanh–sinh agree to ≥52 significant digits for all $k\le400$ (worst 52.375 at $k=400$). Table in §2.
2. **NUMERICAL:** $A_0=1.000\pm0.003$, $\phi_0=\pi(1.000\pm0.002)$; $A_1=-0.05\pm0.05$, $\phi_1=0.85\pm0.06$.
3. **PROVED (symbolic) + NUMERICAL (mechanism):** $B(z_s+\zeta)=(1-si)+s\sqrt2\,\zeta^{1/2}+i\zeta-\frac{s\sqrt2}{2}\zeta^{3/2}+\dots$ — the previous analysis's local constant $s\sqrt2\,\zeta^{1/2}$ checks out; saddle identities of §4.2 are sympy-exact; assembly over **both** branch points gives
$$M_k=-\frac{\sqrt\pi}{k}e^{-\sqrt{\pi k}}\bigl[\cos\sqrt{\pi k}+O(k^{-1/2})\bigr],\quad\text{equivalently}\quad M_k=-\frac{\sqrt\pi}{k+1}e^{-\sqrt{\pi(k+1)}}\bigl[\cos\sqrt{\pi(k+1)}+O(k^{-1/2})\bigr].$$
(The Hankel/saddle dominance steps are standard heuristics, but every local constant is proved and the result matches 52-digit data to the expected O($k^{-1/2}$) accuracy.)
4. **CONJECTURE (analytic) / NUMERICAL (fit):** first correction $\beta_s=\frac{1}{2\sqrt{z_s}}+\frac{\sqrt{z_s}}2-\frac{z_s^{3/2}}{12}$, i.e. $c_r=\frac{1}{2\sqrt\pi}+\frac{\sqrt\pi}4+\frac{\pi^{3/2}}{48}=0.841215$, $c_i=-\frac{1}{2\sqrt\pi}+\frac{\sqrt\pi}4-\frac{\pi^{3/2}}{48}=0.045012$ (fit: $0.860\pm0.006$, $0.050\pm0.003$; ~2% tension in $c_r$ unresolved).
5. **VERDICT on the old claim:** $C=\sqrt\pi/2$ is refuted by 40+ digits; the factor 2 is the sum of the two conjugate branch points $z_\pm=\pm i\pi/2$, exactly as the planner suspected ($\kappa_+=\kappa_-=-\sqrt\pi/2$, they add).

## 6. Files (all printed numbers above come from these)

- `step_01.py` — main: $M_k$, $k\le400$, udeg∈{8,9}, mdps=160, bits=800, umax=1500 → `mk_deg8.tsv`, `mk_deg9.tsv` (runs: 207 s / 704 s).
- `step_02.py` — degree comparison; printed the table of §2 and worst-case 52.375 digits.
- `step_03b.py` — two-consecutive-$k$ phasor estimator; printed $A_k,\phi_k$ table.
- `step_04.py` — global LSQ fits; printed $A_0,\phi_0,A_1,\phi_1$ over 12 windows.
- `step_05.py` — polynomial fits of pointwise estimates (less stable; reported as spread).
- `step_06.py` — numeric branch-point discontinuity: $\gamma_s=s\sqrt2$ to 5+ digits at $r=10^{-4}$.
- `step_07.py` — sympy: exact series of $B$ at $z_s$: $\gamma=\sqrt2$, $c_3=-1/2$, $c_5=5/24$.
- `step_08.py` — sympy: saddle algebra $u_s$, $\phi(u_s)=-2\sqrt{z_s}/\varepsilon+z_s/2-\varepsilon z_s^{3/2}/12$, $\phi''$, $P(u_s)$ ($5/4$ coefficient).
- `step_09.py` — natural-convention fits of $c_r,c_i,d_r,d_i$.
- `step_10.py` — leading-formula ratio test (printed ratios 1.1413→1.0481).
- `step_11.py` — sympy: Gaussian loop corrections $=\varepsilon/(2\sqrt{z_s})-11\varepsilon^2/8$.
- `step_12.py` — residual test with analytic $\beta$: $NR_k$ bounded, envelope $\downarrow$.
- `step_13.py` — free-vs-fixed-$\beta$ fit comparison (locates the 2% tension).
- `step_15.py` — exact-vs-asymptotic ζ-integral at the saddle (ratio 1.0046). (`step_14.py`, `step_16.py`, `step_17.py` contained prefactor/integrand bugs in cross-check code; superseded by `step_15.py` and excluded from evidence.)