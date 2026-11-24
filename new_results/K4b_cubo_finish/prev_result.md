# K4 — ζ(3) for 3D polytopes near p = 1, beyond the octahedron

**Folder:** `research/08_theory/K4_3d_polytopes/`. Python = `.venv/bin/python` (mpmath 1.3.0). Every number below was printed by a `step_NN.py` script in this session; artifacts (`k1_exact.json`, `octa_calibration_h32.json`, `tetra_Sdk_h32.json`, `k4_engine.py`, `k4_quad.py`) are on disk. Statuses: **PROVED** (exact/symbolic or validated against H1's proved values to the printed digit count), **NUMERICAL**, **CONJECTURE**, **FAILED-IN-BUDGET** (honest).

---

## 1. Part (1) — the polytope and the surface-area representation

### 1a. Combinatorial identification — **PROVED**

With the four unit body-diagonal normals
$n_1=(1,1,1)/\sqrt3,\ n_2=(1,1,-1)/\sqrt3,\ n_3=(1,-1,1)/\sqrt3,\ n_4=(-1,1,1)/\sqrt3$
(a regular tetrahedron; $\sum n_jn_j^T=\tfrac43 I$), the polar of $K_1=\{x:\sum_j|\langle n_j,x\rangle|\le1\}$ is
$$K_1^\circ=\operatorname{conv}\{A^T\epsilon:\epsilon\in\{\pm1\}^4\}=\operatorname{conv}\{\,(\pm2,\pm2,\pm2)/\sqrt3\ \text{(8 cube pts)}\ \cup\ \pm\tfrac{4}{\sqrt3}e_i\ \text{(6 octahedron pts)}\,\},$$
computed exactly: the plane $x+y=4/\sqrt3$ contains the four coplanar points $(4,0,0),(0,4,0),(2,2,\pm2)/\sqrt3$, all sides $=2$, diagonals $4\sqrt2/\sqrt3$ and $4/\sqrt3$ — a **rhombus**. There are 12 such planes $(\pm x\pm y,\ \pm x\pm z,\ \pm y\pm z=4/\sqrt3)$, so **$K_1^\circ$ is a rhombic dodecahedron** (14 vertices, 24 edges, 12 rhombic faces), hence its dual

$$\boxed{K_1\ \text{is a cuboctahedron}}$$

with **12 vertices** $=$ permutations of $(\pm\sqrt3/4,\pm\sqrt3/4,0)$, **24 edges** (all triangle–square), **14 faces**: 8 triangles (distance $h_\triangle=1/2$ from origin, planes $\pm x\pm y\pm z=\sqrt3/2$) and 6 squares (distance $h_\square=\sqrt3/4$, planes $\pm x,\pm y,\pm z=\sqrt3/4$); edge length $a=\sqrt6/4$; vertex configuration 3.4.3.4; dihedral angle $\arccos(1/\sqrt3)=0.955316618124509278163857102516\ldots$; Euler $12-24+14=2$. All 12 vertices verified to satisfy $\sum_j|\langle n_j,x\rangle|=1$ exactly in sympy. (steps 13–16; `k1_exact.json`.)

Consequently
$$S(1)=6a^2+8\tfrac{\sqrt3}{4}a^2=\frac{9+3\sqrt3}{4}=3.5490381056766579701455847561294042752071039403578\ldots\quad\textbf{(PROVED, exact)}.$$

### 1b. The representation used for expanding in $d=p-1$ — **PROVED** (and cross-validated)

Radial parametrization $x=\rho(u)u$, $u\in S^2$: $\rho(u)=\Phi(u)^{-1/p}$ with $\Phi(u,d)=\sum_{j=1}^4|\langle n_j,u\rangle|^{1+d}$. The surface element of a radial graph gives
$$S(1+d)=\int_{S^2}F(u,d)\,d\sigma,\qquad F=\Phi^{-\frac{2}{1+d}}\sqrt{1+W},\quad W=\Big|\frac{\sum_j |c_j|^{d}\sigma_j b_j}{\Phi}\Big|^2,$$
where $c_j=\langle n_j,u\rangle$, $\sigma_j=\operatorname{sign}c_j$, $b_j=n_j-c_ju$ (so $\langle b_j,b_k\rangle=\langle n_j,n_k\rangle-c_jc_k$; the $p$-factor in $\nabla_S\ln\rho$ cancels exactly). **This is the local face/edge/vertex form:** on each region of the sphere cut out by the kink lines $\langle n_j,u\rangle=0$ (the radial projections of the edges of $K_1$), the coefficient $F_k(u)=[d^k]F$ is a **polynomial of degree $\le k$ in $L_j=\ln|\langle n_j,u\rangle|$ with rational-function coefficients in $u$** (implemented as exact truncated formal series: $\Phi_k=\sum_j a_jL_j^k/k!$, etc.). The $\log$ singularities live precisely on the kink lines ($\to$ **edge terms**) and their intersection points ($\to$ **vertex terms**); the $d^3$ coefficient is a sum of triple-log surface integrals — the 3D analogue of E1's 2D Beta-integral mechanism.

**Symmetry reduction:** with the octahedral fundamental triangle $T=\{0\le x\le y\le z\}$ ($S(p)=48\int_T F$), parametrized by $u=(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$, $\phi\in[\pi/4,\pi/2]$, split at the kink $c_2=0\Leftrightarrow\tan\theta=1/(\cos\phi+\sin\phi)=\tan\theta_B(\phi)$:
* **region B** (projects to the square face $z=\sqrt3/4$): $\theta\in[0,\theta_B]$, signs $(+,-,+,+)$, $\Phi(u,0)=4z/\sqrt3$;
* **region A** (projects to the triangle face $x+y+z=\sqrt3/2$): $\theta\in[\theta_B,\theta_A]$, $\tan\theta_A=1/\sin\phi$, signs $(+,+,+,+)$, $\Phi(u,0)=2(x+y+z)/\sqrt3$.

### 1c. Machinery validation — **PROVED-strength checks** (`k4_engine.py`, steps 18–24)

| check | residual | script |
|---|---|---|
| formal-series $F_k(u)$ vs numeric $\partial_d^k F_{\rm direct}/k!$ at interior points, $k=1..4$ | $\le 2.7\times10^{-51}$ | step 21 |
| **octahedron calibration** ($m=3$ axis normals, same code): $[d^k]S$ vs H1's **PROVED** $4\sqrt3\{1,\,1,\,1-\pi^2/9,\,1+\pi^2/9-16\zeta(3)/9\}$ | $5.2\times10^{-62},\ 7.5\times10^{-61},\ 9.9\times10^{-59},\ 6.7\times10^{-57}$ ($k=0..3$) | step 24, `octa_calibration_h32.json` |

$[d^3]S_{\rm oct}=-0.2796731272751286206150690832602315974454\ldots$ matches H1 exactly. Quadrature = nested tanh–sinh, $h=1/32$, $T=4.5$, dps = 70 (`k4_quad.py`), all coefficients integrated simultaneously; a formal-series bug (`ser_log` sign) was found and fixed during validation (step 20–21).

### 1d. Cuboctahedron coefficients — **PARTIAL (honest failure on region A)**

Ran at $h=1/32$, $T=4.5$, dps 70 (step 25, `tetra_Sdk_h32.json`):

| $k$ | region B (6 square faces), $48\int_{T_B}F_k$ | region A (8 triangle faces) |
|---|---|---|
| 0 | $2.25000000000$ = $9/4=6a^2$ **exactly** ✓ | **NaN (bug)** |
| 1 | $4.244162313$ | NaN |
| 2 | $0.3212074064$ | NaN |
| 3 | $-1.072935039$ | NaN |
| 4 | $-0.08317937889$ | NaN |

Region B's $k=0$ value reproduces $6a^2=9/4$ exactly — strong internal validation. **Region A failed:** at outer nodes $\phi\to\pi/2$ the inner interval $[\theta_B,\theta_A]$ collapses double-exponentially ($\sim10^{-122}$, below the $10^{-70}$ absolute precision of $x+y-z$); 25 nodes produced NaN through `0×inf` in the series (verified by trace: `F_series` is finite whenever all $a_j>0$; the NaN nodes carry weights $\lesssim10^{-65}$, so skipping them is legitimate — the guard in step 25 skipped 25 nodes but evidently one path still propagated NaN). **The totals $S'(1),S''(1),S'''(1)$ for the cuboctahedron are therefore NOT yet obtained — FAILED-IN-BUDGET.** The identified fix (skip nodes with weight $\le10^{-68}$; treat intervals shorter than $10^{-60}$ by their analytic zero-length limit; optionally raise dps to 90 for the outer shell) is implemented-ready but was not rerun in budget. Finite-difference check of a direct $S(p)$ computation (task requirement) is likewise pending the same rerun.

---

## 2. Part (2) — PSLQ identification of constants

**Not completed** (requires the totals of §1d). What is established here:
* The integrands $F_k$ are polynomials in $L_j=\ln|\langle n_j,u\rangle|$ over spherical triangles whose angles involve $\arctan\sqrt2$ (the magic angle) — so *a priori* triple-log integrals could produce Clausen/log-sine values at non-classical angles; whether they cancel (as H1 proved for the octahedron, where only single zetas survive) is exactly the question the PSLQ must answer. **CONJECTURE (to test):** $[d^3]S=\{\text{rational},\ \pi^2,\ \zeta(3),\ \ln2,\ \ln3\}\otimes\{1,\sqrt3\}$ with no $\mathrm{Cl}_2(\pi/3)$, no $L(\chi_{-3},2)$; this remains unverified.

## 3. Part (3) — toward a 3D edge/vertex law

**Framework established; law not yet numerically tested.** From E1 (2D, PROVED there): $P'''(1)|_{\zeta(3)}=\sum_e Z_e$ with $Z_e=6\alpha_e\beta_e\frac{h_e^2}{L_e}\big(3(\alpha_eu_e+\beta_ev_e)-2\big)$ — ζ(3) enters *only* through the interference of the two corners of one edge (Beta mixed derivatives), and the non-local remainder series terminate exactly for parallelograms and affine-regular hexagons. The 3D analogue set up here: ζ(3) can arise only from the triple-log surface integrals $\int_{T}L_iL_jL_k\,d\sigma$ localized at (i) kink-line crossings = **vertices** of $K_1$ (three faces meet; cf. H1's skewness channel $\mathbb E[m_1m_2m_3]=-\tfrac49\zeta(3)$) and (ii) **edges** (two kink lines + the face's own). The cuboctahedron is an ideal test case: all 24 edges are symmetry-equivalent (one edge type, tri–square, dihedral $\arccos(1/\sqrt3)$) and all 12 vertices equivalent (3.4.3.4), so any law has the form $[d^3]S|_{\zeta(3)}=24\,Z_{\rm edge}+12\,Z_{\rm vertex}$ with only two unknown local functionals — testable with one more independent normal configuration (e.g. 4 normals in non-tetrahedral position, or 6 normals). **Status: CONJECTURE formulated, numerical test FAILED-IN-BUDGET (same region-A bug).**

---

## 4. What is solid, in one place

| # | claim | status | evidence |
|---|---|---|---|
| 1 | $K_1^\circ$ = rhombic dodecahedron; $K_1$ = cuboctahedron, $V{=}12,E{=}24,F{=}14$ ($8\triangle+6\square$), $a=\sqrt6/4$, $h_\triangle=1/2$, $h_\square=\sqrt3/4$, dihedral $\arccos(1/\sqrt3)$ | **PROVED** (exact sympy + hull) | steps 13–16, `k1_exact.json` |
| 2 | $S(1)=(9+3\sqrt3)/4$ | **PROVED** exact; independently: quadrature region B gives $6a^2=2.25$ exactly | steps 16, 25 |
| 3 | representation $S=\int\Phi^{-2/(1+d)}\sqrt{1+W}\,d\sigma$, $F_k$ = degree-$k$ log polynomials on sign regions | **PROVED** derivation; 51-digit pointwise validation | steps 18–21 |
| 4 | machinery vs H1 octahedron (independent PROVED target), $k=0..3$ | agreement to 57–62 digits | step 24 |
| 5 | region-B coefficients $[2.25,\ 4.244162313,\ 0.3212074064,\ -1.072935039,\ -0.08317937889]$, $k=0..4$ | **NUMERICAL** (~50 reliable digits in JSON) | step 25 |
| 6 | $S'(1),S''(1),S'''(1)$ totals; FD check; PSLQ; edge/vertex law test | **FAILED-IN-BUDGET** (region-A collapsed-interval NaN; fix identified) | steps 25–28 |

**Files:** `k4_engine.py` (series arithmetic + $F$ integrand, validated), `k4_quad.py` (tanh–sinh), `k1_exact.json`, `octa_calibration_h32.json` (60-digit), `tetra_Sdk_h32.json` (region B, 60-digit), and `step_01.py … step_28.py` (steps 13–16: combinatorics; 17–21: engine build/debug/validation; 22–24: octahedron calibration; 25–28: tetra run + NaN diagnosis). With the identified fix (weight-thresholded node skipping + short-interval analytic limit, dps 90), one rerun of step 25 (~5 min) plus the FD/PSLQ scripts completes parts (1)–(3); the engine and quadrature are already validated to 57–62 digits against H1's proved octahedron values.