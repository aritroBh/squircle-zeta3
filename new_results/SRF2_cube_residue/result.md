(no FINAL answer; see transcript.json)

Last reply:
FINAL:

# SRF2: the superellipsoid's cube-end series — the pole at p = −1 and its residue

**Work dir:** `new_results/SRF2_cube_residue/` (all runs auto-saved as `step_NN.py`, artifacts `EV_coeffs_K40.json`, `srf2_lib.py`). Every number quoted below was printed by my code.

---

## 0. Setup

With $\varepsilon=1/p$,
$$S(p)=24+48\,\varepsilon E(\varepsilon)+48\,\varepsilon^2 V(\varepsilon),$$
$$E(\varepsilon)=\int_0^1 a^{2\varepsilon-1}\big[(1+a)^{-1-2\varepsilon}\sqrt{1+a^{2-2\varepsilon}}-1\big]\,da,\qquad
V(\varepsilon)=\iint_{[0,1]^2} a^{2\varepsilon-1}c^{\,\varepsilon-1}\big[g-g_0\big]\,da\,dc,$$
$g=N^{-1-2\varepsilon}S$, $N=1+a+ac$, $S=\sqrt{1+a^{2-2\varepsilon}+(ac)^{2-2\varepsilon}}$, $g_0=g(a,0)$, $S_0=\sqrt{1+a^{2-2\varepsilon}}$.

**Notation for this report:** $q=1+2\varepsilon$, $u=ac/(1+a)$, $B_1(a)=\partial_c(g-g_0)|_{c=0}=-q\,a(1+a)^{-2-2\varepsilon}S_0$, $I=\int_0^1(\sqrt{1+a^4}-1)\,da/a^2$.

---

## 1. Task 1 — the pole and its residue (PROVED)

### 1.1 Theorem (pole of $V$ at $\varepsilon=-1$; residue exact) — PROVED, and confirmed numerically to 7 digits

$$V(\varepsilon)=\frac{R}{1+\varepsilon}+H(\varepsilon),\qquad H\text{ holomorphic near }-1,\qquad
R=I-2=-1-\frac{2\pi^{3/2}}{\Gamma(1/4)^2}=-1.847213084793979086606499123482191636481\ldots$$

**Proof.** (i) *Decomposition.* For $\operatorname{Re}\varepsilon>-\tfrac12$ the integral defining $V$ converges, and subtracting the $c$-linear term $cB_1$ gives, identically,
$$V=W+P_c,\qquad W=\iint a^{2\varepsilon-1}c^{\,\varepsilon-1}\big[g-g_0-cB_1\big]\,da\,dc,\qquad
P_c=-\frac{q}{\varepsilon+1}\,J(\varepsilon),$$
$$J(\varepsilon)=\int_0^1 a^{2\varepsilon}(1+a)^{-2-2\varepsilon}\sqrt{1+a^{2-2\varepsilon}}\,da .$$
With $g-g_0-cB_1=(1+a)^{-q}\big[\Phi(u)S+(1-qu)(S-S_0)\big]$, $\Phi(u)=(1+u)^{-q}-1+qu=O(u^2)$, the bracket is $O(a^2c^2+(ac)^{2-2\varepsilon})$ at the corner; at $\varepsilon=-1$ it collapses to $N(ac)^4/(S+S_0)$, so $W$ is given by a convergent integral on $-1<\varepsilon\le0$ and is continuous through $-1$ (bounds uniform as in the paper's 2D lemma `lem:holo`). Hence $(1+\varepsilon)W(\varepsilon)\to0$.

(ii) *The 1D integral $J$ and the boundary layer.* Write $\eta=\varepsilon+1$ and
$J=\frac{1}{2\varepsilon+1}+J_1$, $J_1=\int_0^1 a^{-2+2\eta}f_1\,da$, $f_1=(1+a)^{-2\eta}\sqrt{1+a^{4-2\eta}}-1=-2\eta\,a+O(a^2)$ uniformly in $\eta\in[0,\eta_0]$. Splitting $f_1=-2\eta a+(f_1+2\eta a)$: the linear piece contributes $-2\eta\int_0^1 a^{-1+2\eta}da=-1$ *exactly for every $\eta>0$* (a Dirac-delta mass at $a=0$: $-2\eta a^{-1+2\eta}\rightharpoonup-\delta_0$); the remainder is $O(a^{-2+2\eta}\cdot a^2)$ and dominated convergence gives $\to\int_0^1 a^{-2}(\sqrt{1+a^4}-1)\,da=I$. Therefore
$$J(-1^+)=-1+\bigl(-1+I\bigr)=I-2 .$$
Note the subtlety this resolves: the *face value* $J_1(0)=I$ differs from the *limit* $J_1(0^+)=-1+I$ by exactly $-1$ — the planner's "second $-1$ from the subtracted linear term $-(1+2\varepsilon)ac$". (Numerically confirmed: reliable split-integration gives $J(-1)=-0.8472130847939790866064991$ ($=-1+I$, 25 digits) while $-qJ(-1+\eta)\to-1.847213$.)

(iii) *Residue.* $R=\lim_{\varepsilon\to-1^+}(\varepsilon+1)V=\lim\big[-qJ(\varepsilon)\big]=-(-1)\,J(-1^+)=I-2$.

(iv) *Closed form of $I$ (PROVED).* Integration by parts: $I=1-\sqrt2+2K$, $K=\int_0^1\frac{a^2}{\sqrt{1+a^4}}da$. The substitution $t=\frac{\sqrt2\,a}{\sqrt{1+a^4}}$ (monotone $[0,1]\to[0,1]$, verified symbolically: integrand transforms to $2\sqrt2\,a^2(1+a^4)^{-3/2}da$) gives
$$L':=\int_0^1\frac{t^2}{\sqrt{1-t^4}}dt=2\sqrt2\int_0^1\frac{a^2}{(1+a^4)^{3/2}}da .$$
The algebraic identity $\dfrac{a^2}{(1+a^4)^{3/2}}=\dfrac12\dfrac{d}{da}\dfrac{a^3}{\sqrt{1+a^4}}-\dfrac12\dfrac{a^2}{\sqrt{1+a^4}}$ (sympy: difference simplifies to 0) yields $\int_0^1\frac{a^2}{(1+a^4)^{3/2}}da=\frac{1}{2\sqrt2}-\frac K2$, i.e. $L'=1-\sqrt2\,K$. But $L'=\frac14B\!\left(\frac34,\frac12\right)=\frac{\sqrt2\,\pi^{3/2}}{\Gamma(1/4)^2}$ (sympy-verified chain: $\Gamma(5/4)=\frac14\Gamma(1/4)$, $\Gamma(3/4)\Gamma(1/4)=\pi\sqrt2$). Hence $K=\frac{1}{\sqrt2}-\frac{\pi^{3/2}}{\Gamma(1/4)^2}$ and
$$I=1-\frac{2\pi^{3/2}}{\Gamma(1/4)^2}=0.1527869152060209133935008765178083635186\ldots\quad\therefore\quad R=-1-\frac{2\pi^{3/2}}{\Gamma(1/4)^2}.\ \blacksquare$$

Every proof step was also verified numerically to 39 printed digits (direct quadrature = closed form for $I$, $K$, $L'$, both substitutions/IBP identities).

### 1.2 Full singularity inventory in $|\varepsilon|\le1+\delta$

**$E$ — no singularity at $\varepsilon=-1$ (PROVED).** The continuation subtracts terms $c_l(\varepsilon)a^l$ from the bracket; the pole candidate at $\varepsilon=-l/2$ has residue $\propto\binom{-(1+2\varepsilon)}{l}\big|_{\varepsilon=-l/2}=\binom{0}{l}=0$ for all $l\ge1$ — because at $\varepsilon=-l/2$, $(1+a)^{-1-2\varepsilon}=(1+a)^{l-1}$ is a polynomial of degree $l-1$, so the bracket has no $a^l$ term. Hence $E$ is holomorphic on $\operatorname{Re}\varepsilon<1$; its only singularities are the poles at $\varepsilon=\frac{l+2j}{2(j-1)}$ ($j\ge2$, $l\ge0$), all $>1$, accumulating at $1^+$. In particular $E$ is regular at $-1$ (finite subtracted-integral value $-0.6136330724710304038\ldots$; the analytic-continuation *limit* $E(-1^+)=-\frac32+I_E$ differs by a boundary-layer constant $-\frac12$, same phenomenon as (ii), with no pole to multiply). So $E_k$ decays like the 2D $M_k$, exactly as SRF1 observed.

**$V$ — poles at the negative integers, and the $1^+$ families (PROVED structure; residues as shown).** From the $c$-side terms $c^mB_m$ (weight $c^{\varepsilon-1}$): simple poles at $\varepsilon=-m$, $m=1,2,3,\ldots$, with residues $\binom{2m-1}{m}K^{(m)}(-m^+)\neq0$, $K^{(m)}(\varepsilon)=\int_0^1 a^{2\varepsilon-1+m}(1+a)^{-q-m}S_0\,da$ evaluated by the same Hadamard + boundary-layer rule. There are **no** negative half-integer poles (at $\varepsilon=-l/2$, $N^{-q}=N^{l-1}$ is a polynomial in $a$ of degree $l-1$: no $a^l$ term). The remaining singularities are the accumulating families $\varepsilon=\frac{l+2j}{2(j-1)}$ ($a$-side) and $\varepsilon=\frac{j+2m}{2m-1}$ ($c$-side), all $>1$, accumulating at $1^+$. So on $|\varepsilon|\le1$ the only singularities are $\varepsilon=-1$ (simple pole, residue $R$) and $\varepsilon=+1$ (accumulation point).

Computed residues (formula PROVED, values NUMERICAL, confirmed by independent reliable continuation):
$$\rho_{-1}=R=I-2=-1.8472130847939790866,$$
$$\rho_{-2}=3\bigl(-2+I_2'\bigr)=-5.3806643388446858024,\qquad I_2'=\int_0^1\frac{(1+a)\sqrt{1+a^6}-1-a}{a^3}\,da,$$
$$\rho_{-3}=10\bigl(-\tfrac83+I_3\bigr)=-23.564129294356354798,\qquad I_3=\int_0^1\frac{(1+a)^2\sqrt{1+a^8}-1-2a-a^2}{a^4}\,da .$$

**Numerical confirmation of the residue (NUMERICAL, 7 digits).** Boundary-layer-aware split-integration of $J$ gives $-qJ(-1+\eta)=-1.4423,\,-1.6318,\,-1.7360,\,-1.8019,\,-1.8244,\,-1.8358,\,-1.8426$ at $\eta=0.2,0.1,0.05,0.02,0.01,0.005,0.002$; polynomial extrapolation to $\eta=0$: $-1.84721279$ (deg-2, residual $3.0\times10^{-7}$), $-1.84721297$ (deg-3, residual $1.2\times10^{-7}$) vs exact $-1.84721308$. And $(-1)^{40}V_{40}=-1.84721308989$ from the K40 coefficients (residual $5.1\times10^{-9}$, the predicted $\delta_{40}$).

---

## 2. Task 2 — exact asymptotics of $s_k$ (NUMERICAL)

$$s_k=48\bigl(E_{k-1}+V_{k-2}\bigr)=48R(-1)^k+r_k,\qquad r_k=48\bigl(E_{k-1}+\delta_{k-2}\bigr),\quad \delta_k:=V_k-R(-1)^k .$$

The correction decomposes as $\delta_k=(-1)^k\bigl[\frac{\rho_{-2}}{2^{k+1}}+\frac{\rho_{-3}}{3^{k+1}}+\cdots\bigr]+\delta_k^{(+1)}$:

- **$k\lesssim17$: dominated by the pole at $\varepsilon=-2$.** Printed $\delta_k\cdot2^{k+1}$ (alternating sign removed): $-4.75,-5.65,-5.89,-5.95,-5.80,-5.67,-5.64,-5.43,-5.58,-5.36,-5.34$ for $k=5\ldots15$, against $\rho_{-2}=-5.380664$; ratios of successive $|\delta_k|\to2$.
- **$k\gtrsim20$: dominated by the $\varepsilon\to1^+$ accumulation** (the radius-zero mechanism at the octahedron end): non-alternating, decaying in sign blocks, magnitudes $\approx4\times10^{-6}$ ($k=20$), $4.5\times10^{-7}$ ($k=30$), $5\times10^{-9}$ ($k=40$) — the 3D analogue of the 2D $M_k\sim e^{-\sqrt{\pi k}}$ tail. The $-3$ pole ($\rho_{-3}=-23.564129$; crude data fit $-23.9$ at $k=9$) sits between the two regimes.
- $E_k$ also decays (2D-like) and is smaller than $\delta_k$ here ($|E_{10}|=1.45\times10^{-5}$ vs $|\delta_{10}|=2.8\times10^{-3}$).

Checks against SRF1 (all printed): $|s_k|\to48|R|=88.666228\ldots$; $48(E_4+V_3)=90.7414583$ = SRF1's $s_5=90.7414583339$, decomposed as $88.6662281+2.0752303$; $48(E_5+V_4)=-93.2703788$ = SRF1's $s_6=-93.2703787987$, decomposed as $-88.6662281-4.6041507$. Root tests of SRF1 reproduced: radius exactly 1, dominant boundary point $\varepsilon=-1$.

---

## 3. Task 3 — $V_k$ to order 40, ≥30 digits where it matters (NUMERICAL)

Recomputed with the SRF1 tanh–sinh/`arb_series` engine at dps = 110 ($h=1/64$, $T=4.6$, 118 s; control $h=1/32$, 30 s). Relative $h$-differences: $1.2\times10^{-51}$ ($k=10$), $2.0\times10^{-38}$ ($k=20$), $2.0\times10^{-27}$ ($k=30$), $7.2\times10^{-18}$ ($k=40$) — i.e. ≥30 digits for $k\lesssim28$, 27 digits at $k=30$, 17 at $k=40$ (saved in `EV_coeffs_K40.json`; $E_k$ to order 41 alongside). Values through $V_7$ match the K30 file and were independently reproduced (h-convergence $2\times10^{-36}\ldots4\times10^{-28}$):

$V_0=-0.4404654980318348925$, $V_1=1.1567597697434974306$, $V_2=-1.6757604153964917033$, $V_3=1.9041854807791881529$, $V_4=-1.9460518042914523510$, $V_5=1.9213770634099109876$, $V_6=-1.8913219386639152842$, $V_7=1.8702381647127794442$.

Remainder behaves as predicted: $(-1)^kV_k-R\to0$, and after subtracting $R(-1)^k$ the residual is the $-2$ pole for $k\lesssim17$ and the $1^+$ stretched-exponential accumulation beyond (values tabulated in §2).

---

## 4. Paragraph for the paper (plain words; status of each claim)

> At the cube end, the series $S(p)=\sum_k s_k p^{-k}$ converges for every $p>1$ and its radius of convergence is exactly one (*radius $\ge1$: PROVED modulo the $V$-bounds lemma that follows the paper's 2D `lem:holo` technique; radius $\le1$: NUMERICAL*). Its coefficients settle to alternating constants: $s_k\approx48R\,(-1)^k$ with $|s_k|\to48|R|=88.666228\ldots$ (*NUMERICAL*, 40 coefficients). The nearest singularity that sets their size is at $p=-1$: the vertex term $V$ has a simple pole at $\varepsilon=1/p=-1$ with residue $R=-1-2\pi^{3/2}/\Gamma(1/4)^2=-1.8472130847939790866\ldots$ (*PROVED*, and confirmed to 7 digits both by direct continuation and by the coefficient limit). The residue is a regularised corner integral, $R=I-2$ with $I=\int_0^1(\sqrt{1+a^4}-1)\,da/a^2=1-2\pi^{3/2}/\Gamma(1/4)^2$ (*closed form PROVED* by Beta functions); of the "$-2$", one unit is the Hadamard finite part of $\int_0^1 a^{-2}\sqrt{1+a^4}\,da$ and the other is a boundary-layer delta mass contributed by the subtracted linear corner term $-(1+2\varepsilon)ac$ — which is why the naive regularised integral misses it. This is a vertex effect absent in 2D, where $\varepsilon=+1$ is the only boundary point of the squircle's $M(\varepsilon)$: in 3D the edge term $E$ is provably regular at $-1$ (*PROVED*: all candidate residues vanish) and its coefficients decay like the 2D $M_k$, while the vertex term $V$ carries the pole. The subleading corrections to $s_k$ come first from the next vertex pole at $\varepsilon=-2$ (residue $-5.380664\ldots$; formula $3(-2+I_2')$ PROVED, value NUMERICAL to 4 digits) and are then governed by the pole accumulation at $\varepsilon\to1^+$ — the octahedron-end radius-zero mechanism made visible — whose stretched-exponential tail overtakes the $2^{-k}$ terms around $k\approx20$ (*NUMERICAL*).

---

## 5. Evidence table

| # | claim | status | evidence (residual / agreement) | digits | script |
|---|-------|--------|---------------------------------|--------|--------|
| 1 | $I=1-2\pi^{3/2}/\Gamma(1/4)^2$ | PROVED | direct quad vs closed form identical to 39 digits; IBP + substitution $t=\sqrt2a/\sqrt{1+a^4}$ + $\frac14B(\frac34,\frac12)$ each verified symbolically (sympy: True/0) and numerically | 39 | step_38–40 |
| 2 | $R=\operatorname{Res}_{-1}V=I-2=-1.8472130847939790866$ | PROVED (derivation); NUMERICAL (confirmations) | $-qJ$ extrapolation $-1.84721279$ vs exact, residual $3\times10^{-7}$; $(-1)^{40}V_{40}=-1.84721308989$, residual $5\times10^{-9}$ | 7 / 8 | step_36–39, 49–50 |
| 3 | boundary-layer jump $J(-1^+)=J(-1)-1$ | PROVED (analysis) + NUMERICAL | face value $J(-1)=-0.8472130847939790866064991$ matches $-1+I$ | 25 | step_36–37 |
| 4 | decomposition $V=W+P_c$ | NUMERICAL | vs raw vlib $V(0.1)$: residual $5\times10^{-18}$; vs 31-term Taylor at $-0.3$: $1.6\times10^{-10}$ | 10–18 | step_25, 35 |
| 5 | $E$ holomorphic at $-1$ | PROVED (vanishing residues $\binom{0}{l}=0$) | finite value $-0.6136330724710304038$; $E(0.1)$ matches vlib to 15 digits | 15+ | step_45–46 |
| 6 | $\rho_{-2}=3(-2+I_2')=-5.3806643388$ | formula PROVED; value NUMERICAL | continuation $3K^{(2)}(-2+\eta)\to-5.3728$ ($\eta=0.005$) $\to-5.3807$; $\delta_k2^{k+1}\in[-5.95,-5.34]$ for $k=5..15$ | 4 | step_47, 50–51 |
| 7 | $\rho_{-3}=10(-8/3+I_3)=-23.564129$ | formula PROVED; NUMERICAL crude | pair-fit at $k=9$: $-23.9$ | 1–2 | step_52 |
| 8 | $V_k$, $k\le40$ | NUMERICAL | $h$-convergence $1.2{\rm e}{-51}$ ($k=10$), $2{\rm e}{-38}$ ($k=20$), $2{\rm e}{-27}$ ($k=30$), $7.2{\rm e}{-18}$ ($k=40$) | 17–51 | step_49–50 (`EV_coeffs_K40.json`) |
| 9 | $\delta_k$: $-2$ pole then $1^+$ accumulation | NUMERICAL | §2 table: $\delta_k2^{k+1}\to-5.4$ ($k\le16$); accumulation $4\times10^{-6}$ ($k=20$) $\to 5\times10^{-9}$ ($k=40$) | — | step_50–51 |
| 10 | $s_k=48R(-1)^k+r_k$ | NUMERICAL | $s_5$: $88.6662281+2.0752303=90.7414584$ vs SRF1 $90.7414583339$; $s_6=-93.2703788$ vs $-93.2703787987$ | 7 | step_50 + SRF1 |

Open item: rigorous uniform bounds for $W$ on $|\varepsilon+1|\le r<1$ (same technique as the paper's 2D `lem:holo`; not written out in full), and a closed form for $I_2', I_3$ (not attempted beyond the integral forms above).

---

## 6. `verify_srf2.py` (deliverable; runtime ≈ 2 min < 5 min)

All parts are verbatim from the executed steps indicated; outputs quoted above.

```python
# verify_srf2.py -- SRF2: pole of the vertex term V at eps=-1 and its residue
# R = -1 - 2*pi^{3/2}/Gamma(1/4)^2.  Runs in ~2 min.  Work dir: new_results/SRF2_cube_residue.
import sys, json, time
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
sys.path.insert(0, '<project>/new_results/SRF2_cube_residue')
import mpmath as mp
import flint
from flint import arb, ctx
from vlib import ts_nodes, V_series, tompf
from srf2_lib import expm1_s, expm1m, u_minus_log1p, Phi, tomp   # step_24, step_26 (stable integrands)

t0 = time.time(); ctx.dps = 60; mp.mp.dps = 60
G4 = mp.gamma(mp.mpf('0.25'))
R  = -1 - 2*mp.pi**mp.mpf('1.5')/G4**2
I  = 1 - 2*mp.pi**mp.mpf('1.5')/G4**2

# PART 1 (step_38): closed form of I and R, every proof step numerically
I_dir = mp.quad(lambda a: a*a/(mp.sqrt(1+a**4)+1), [0,1])
K     = mp.quad(lambda a: a*a/mp.sqrt(1+a**4), [0,1])
Lp    = mp.quad(lambda t: t*t/mp.sqrt(1-t**4), [0,1])
assert abs(I_dir - I) < mp.mpf('1e-50')
assert abs(1-mp.sqrt(2)+2*K - I) < mp.mpf('1e-50')                      # integration by parts
assert abs(Lp - mp.sqrt(2)*mp.pi**mp.mpf('1.5')/G4**2) < mp.mpf('1e-50') # (1/4)B(3/4,1/2)
assert abs(Lp - 2*mp.sqrt(2)*mp.quad(lambda a: a*a/(1+a**4)**mp.mpf('1.5'), [0,1])) < mp.mpf('1e-50')  # substitution
assert abs(mp.quad(lambda a: a*a/(1+a**4)**mp.mpf('1.5'), [0,1]) - (1/(2*mp.sqrt(2))-K/2)) < mp.mpf('1e-50')
print("PART 1 ok: I =", mp.nstr(I_dir, 30), " R =", mp.nstr(R, 30))

# PART 2 (step_36/37): residue via boundary-layer-aware split integration of J
def J_small(eps, a0=mp.mpf('0.5'), Lam=mp.mpf('160')):
    s = -2-2*eps; Lmax = 200
    g1 = [mp.mpf(1)]
    for l in range(1, Lmax+1): g1.append(g1[-1]*(s-l+1)/l)
    bh = [mp.mpf(1)]
    for j in range(1, 60): bh.append(bh[-1]*(mp.mpf('0.5')-j+1)/j)
    two2e = 2-2*eps; tot = mp.mpf(0)
    for j in range(0, 60):
        for l in range(0, Lmax+1):
            if l == 0 and j == 0: continue
            pw = l + two2e*j
            if pw > Lam: break
            c = g1[l]*bh[j]
            if c == 0: continue
            mu = 2*eps + pw
            tot += c * a0**(mu+1)/(mu+1)
    return tot
def J_big(eps, h=arb(1)/64, T=5.0):
    a0 = arb('0.5'); s = arb(0)
    for (x, xm, w, par) in ts_nodes(h, T):
        a = a0 + (1-a0)*x; la = a.log(); l1a = (1+a).log()
        A2 = a*a*(arb(-2*eps)*la).exp(); S0 = (1+A2).sqrt()
        f = S0*expm1_s(arb(-(2+2*eps))*l1a) + A2/(S0+1)
        s += (arb(2*eps)*la).exp()*f*w*(1-a0)
    return s
def J_rel(e):
    return 1/(2*e+1) + J_small(e) + tomp(J_big(arb(mp.nstr(e, 30))))
assert abs(J_rel(mp.mpf(-1)) - (-1 + I)) < mp.mpf('1e-24')              # face value = -1 + I
pts = []
for eta in ['0.05', '0.02', '0.01', '0.005']:
    e = -1 + mp.mpf(eta); q = 1+2*e
    pts.append((mp.mpf(eta), -q*J_rel(e)))
A = mp.matrix([[x**c for c in range(3)] for x, y in pts]); b = mp.matrix([y for x, y in pts])
R_est = mp.lu_solve(A, b)[0]
print("PART 2 ok: -qJ extrapolated R =", mp.nstr(R_est, 15), " residual:", mp.nstr(abs(R_est-R), 3))
assert abs(R_est - R) < mp.mpf('2e-6')

# PART 3 (step_47): rho_{-2} = 3(-2 + I2')
I2 = mp.quad(lambda a: (1+a)*a**3/(mp.sqrt(1+a**6)+1), [0, mp.mpf('0.3'), 1])
rho2 = 3*(-2+I2)
print("PART 3 ok: rho_{-2} =", mp.nstr(rho2, 20))

# PART 4 (step_49/50): V_k to order 40 (h=1/32) vs stored high-precision values + delta_k behaviour
flint.ctx.cap = 64; ctx.dps = 100
Vf, _ = V_series(0, 41, h=arb(1)/32, T=4.6)
V = [tompf(x, 110) for x in Vf]
stored = [mp.mpf(x) for x in json.load(open('EV_coeffs_K40.json'))['V']]
worst = max(abs(V[k]-stored[k])/abs(stored[k]) for k in range(1, 41))
print("PART 4 ok: worst rel diff vs K40 store =", mp.nstr(worst, 3))
d40 = (-1)**40*V[40] - R
d10 = ((-1)**10*V[10] - R)*2**11
print("  (-1)^40 V_40 - R =", mp.nstr(d40, 5), " (accumulation);  2^11 delta_10 =", mp.nstr(d10, 8), " vs rho_-2 =", mp.nstr(rho2, 8))
assert abs(d40) < mp.mpf('1e-7') and abs(d10 - rho2)/abs(rho2) < mp.mpf('0.06')

# PART 5 (step_45/46): E regular at -1
def E_int(eps, h=arb(1)/32, T=4.6):
    q = 1+2*eps; s = arb(0)
    for (a, am, wa, pa) in ts_nodes(h, T):
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(arb(-2*eps)*la).exp(); S0 = (1+A2).sqrt()
        f = Phi(a, q)*S0 + (1-q*a)*A2/(S0+1)
        s += (arb(2*eps)*la).exp()/a*f*wa
    return -arb(1+2*eps)/(2*eps+1) + s
E_m1 = tomp(E_int(arb(-1)))
print("PART 5 ok: E(-1) (subtracted-integral value) =", mp.nstr(E_m1, 20), " (finite: no pole)")

# PART 6 (step_50): s_k = 48(E_{k-1}+V_{k-2}) vs SRF1
E3 = [mp.mpf(x) for x in json.load(open('EV_coeffs_K40.json'))['E3']]
s = [mp.mpf(24)] + [48*(E3[k-1] + (stored[k-2] if k >= 2 else 0)) for k in range(1, 42)]
print("PART 6 ok: s5 =", mp.nstr(s[5], 12), " s6 =", mp.nstr(s[6], 12),
      " 48R(-1)^5 =", mp.nstr(48*R*(-1)**5, 12))
assert abs(s[5] - mp.mpf('90.7414583339')) < mp.mpf('1e-8') and abs(s[6] - mp.mpf('-93.2703787987')) < mp.mpf('1e-8')
print(f"ALL CHECKS PASSED in {time.time()-t0:.0f}s")
```

### Step files that produced the evidence (auto-saved in the work dir)
- `step_20–22.py`: load K30 coefficients; coefficient limit $(-1)^kV_k\to-1.84721308479$; reproduce vlib $V_0..V_7$ (h-diff $2{\rm e}{-36}\ldots4{\rm e}{-28}$).
- `step_23–31.py`: stable continuation library `srf2_lib.py`; bug hunt (expm1m→expm1_s in $J$); decomposition $V=W+P_c$ validated vs raw vlib at $0.1$ (18 digits) and vs Taylor at $-0.3,-0.6$.
- `step_32–37.py`: diagnosis of the mass-hiding/boundary layer; boundary-layer-aware split integration; $J(-1)=-0.84721308479$ (25 digits) vs $-qJ(-1+\eta)\to-1.847213$ — the jump of $-1$.
- `step_38–40.py`: closed-form proof of $I$ (39-digit checks + sympy: substitution, IBP identity = 0, Beta chain).
- `step_39.py`: Richardson extrapolation of $-qJ$: $R=-1.84721279$ (residual $3\times10^{-7}$).
- `step_44–46.py`: $E$ stable at $-1$ ($\Phi$-form); $E(-1)=-0.6136330724710304$; $E(0.1)$ vs vlib (15 digits).
- `step_47–48.py`: $\rho_{-2}=3(-2+I_2')=-5.3806643388$; continuation $3K^{(2)}(-2+\eta)\to-5.3807$.
- `step_49–50.py`: $V_k$, $E_k$ to order 40/41 at dps 110 (h=1/64, 118 s; control h=1/32, 30 s) → `EV_coeffs_K40.json`; $\delta_k$ analysis.
- `step_51–52.py`: remainder decomposition ($-2$ pole vs $1^+$ accumulation); $\rho_{-3}=10(-8/3+I_3)=-23.564129$.

**Summary of statuses:** the pole at $\varepsilon=-1$ with $R=-1-2\pi^{3/2}/\Gamma(1/4)^2$ and the closed form $I=1-2\pi^{3/2}/\Gamma(1/4)^2$ are PROVED (with 7-digit independent numerical confirmation); the regularity of $E$ at $-1$ is PROVED; the full singularity inventory is PROVED in structure; the residue formulas at $-2,-3$ are PROVED as formulas with NUMERICAL values confirmed to 4 and 1–2 digits; the coefficient asymptotics, the next-singularity picture, and all high-order coefficients are NUMERICAL.