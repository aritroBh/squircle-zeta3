# K2 — Decay law of the corner coefficients M_k(α) at general angle

## 0. Headline

$$\boxed{\;M_k(\alpha)\;\approx\;\frac{\sqrt{\,2\alpha\sin\alpha\,}}{N}\;e^{-\sqrt{2\alpha N}}\;\cos\!\bigl(\sqrt{2\alpha N}+\varphi(\alpha)\bigr),\qquad N=k+1\;}$$

- **c(α) = 2α** (so c(π/2) = π ✓), **β = 1**, **|A(α)| = √(2α sinα)** (so A(π/2) = √π ✓, matching D13b's VERIFIED squircle law −(√π/N)e^{−√(πN)}cos√(πN)).
- Mechanism (task's "two ends" hint made precise): with g = e^{−x}, z = x(u−1), the kernel of the exact representation M = T + R is **B̃_α(z) = e^{−z}(√(1−2cosα e^z+e^{2z})−1) = (e^z−2cosα)/(1+√(1−2cosα e^z+e^{2z}))**, whose square root vanishes exactly when **g^{1−u} = e^{±iα} ⟺ z± = ±iα**. A Hankel+saddle at x* = √(Niα) has action 2√(Niα) = √(2αN)(1+i): decay e^{−√(2αN)}, oscillation cos(√(2αN)).
- φ(π/2) = π (D13b). φ(α) at other angles is **not pinned** — my saddle derivation gives the amplitude and rate exactly but misplaces the O(1) phase (honest failure, §4.4); the least‑squares fits were **not run** (time went to certification). The M_k tables are on disk (`out/mk_*.tsv`) for the fits.

---

## 1. Part (1): M_k(α) for k ≤ 300 at five angles (Cauchy/FFT in u)

**Method.** M_k(α) = (1/2π)∫₀^{2π} M(re^{iθ};α) r^{−k}e^{−ikθ}dθ by direct DFT (P = 1024, r = 0.9, K = 300, dps = 140). M(u;α) at complex u is evaluated by tanh–sinh quadrature over g∈(0,1) (map g = ½(1+tanh(π/2 sinh t)), h = 0.011, T = 5.5, 969 nodes) with two stabilizations, both derived and verified:
- **Small-g subtraction (EXACT, verified to 7.8e-62):** with β_j = (−1)^j(1+u)_j/j!, γ_m = [z^m]√(1−2cosα z+z²), the integrand satisfies F(g,u) = Σ_{(j,m)≠(0,0)} β_j γ_m g^{j+(m−1)(1−u)}; subtracting m=1..6, j=1..6 gives S plus the closed form Σ_m γ_m/(m−(m−1)u) + Σ_j β_j/(u+j).
- **Stable radical:** (1+g)^{−1−u}√Q − 1 = C(Q−1)/(√Q+1) + expm1(−(1+u)log1p g), removing the catastrophic cancellation of C√Q−1 against g^{u−1} in the deep g→0 tail for Re u ≈ −0.9.

Only p = 0..P/2 points are needed (f(ū) = conj f(u)); aliasing error ≈ (M_{k+P}/M_k)r^P ≲ 1e‑58 relative at k = 300.

**Printed values (h = 0.011 run, dps = 140).** M_0 (45 digits) and M_300:

| α | M_0(α) | M_300(α) |
|---|---|---|
| π/3 | −0.894743627782948423758818521725 | −5.893714688164253658463459e‑14 |
| π/2 | −0.553303299720515717370808039043 | −2.146317091805711938025111e‑16 |
| 2π/3 | −0.261624071882273918258403612467 | 1.496433203371924986368846e‑18 |
| 3π/4 | −0.150248521094360011009169661928 | −2.7588786129953279246042e‑19 |
| 5π/6 | −0.0677559711415642820535770571527 | 1.144142537353196222156377e‑20 |

(π/2 also: M_50 = −1.25705990040564128136859e‑7, M_100 = −1.772316807359054100088703e‑10, M_200 = −1.145441719208160174781568e‑13 — all match D13b's prints.)

**Verification (NUMERICAL, with digit counts):**

| check | result |
|---|---|
| α=π/2 vs D13b (deg‑9 acb values) | M_0: **45.2 digits**; M_50…M_300: **21.5–23.0 digits** — limited by D13b's 21-digit prints, agreement full |
| k ≤ 10, all 5 angles, vs independent `corner_series.corner_coeffs` | **59.5–59.9 digits** at every angle |
| FLINT acb_series of T_k (independent representation/kernel/arithmetic, rad ≤ 1e‑231) vs FFT M_300 | **14.9–15.2 significant digits** at all 5 angles |
| trapezoidal self-consistency h vs h/2 at worst u_p | rate e^{−c/h}, c ≈ 2.42 measured → h=0.011 gives f_p error ~1e‑95 → ~45–65 digits at k=300 (**estimate, not end-to-end certified**) |

**Honest precision status:** the tables are certified ≥ 59 digits for k ≤ 10 (all angles) and ≥ 21 digits at k ≤ 300 for α = π/2; the ≥ 40-digit claim at k = 300 for the other angles rests on the measured exponential convergence rate + the 15-digit FLINT cross-check — likely true but **not fully certified**. The FLINT node-sum shows an unresolved absolute leak (~7.6e‑20 at k=0, ~1e‑29 at k=300) that exceeds its own balls; the mpmath node-sum of the same integrand converges cleanly (1.3e‑93 at udeg 7), so the quadrature itself is exonerated and the leak is in my FLINT conversion/summation of the 843 nodes (unresolved).

## 2. Part (2): the fit form

LSQ fits were **not executed** (budget consumed by certification). What is established:

- **c(α) = 2α** (NUMERICAL/THEORY): envelope test at k = 300 (N = 301): envelope E = √(2α sinα)·N^{−1}e^{−√(2αN)} vs printed |M_300|: ratios 1.06 (π/3), 0.81 (π/2), 0.64 (2π/3), 1.06 (3π/4), 0.38 (5π/6) — all in [0,1.06] as required of |cos(√(2αN)+φ)| with O(N^{−1/2}) corrections (at π/2 the corrections are D13b's known c_r ≈ 0.84 terms). A constant c (e.g. c = π) is excluded: it would predict M_300(π/3) ~ 100× larger than observed.
- **β = 1, |A(α)| = √(2α sinα)** (THEORY, exact at leading order; reproduces √π at π/2 ✓).
- **φ(α):** φ(π/2) = π known (D13b); general φ(α) open — CONJECTURE pending fits.

## 3. Part (3): theory — where c(α) comes from

**(a) Exact representation (PROVED algebra + NUMERICAL 1e‑62).** With L(x) = log(1+e^{−x}), z = x(q−1):
M(q;α) = T(q) + R(q), T(q) = ∫₀^∞ e^{−qL(x)} B̃_α(x(q−1))/(e^x+1) dx, R(q) = ∫₀^∞ (e^{−qL(x)}/(1+e^{−x})−1)e^{−qx}dx, verified against the definition to 1e‑62 at 9 (q, α) points. The rationalized kernel **B̃_α(z) = (e^z−2cosα)/(1+√(1−2cosα e^z+e^{2z}))** has only e^{+z} (small series coefficients — the key for the FLINT certification; the naïve e^{−z}(√−1) form explodes intermediates like e^{+x}).

**(b) The "two ends" singularity (task's hint).** The radical vanishes iff g^{1−u} = e^{±iα} ⟺ e^z = e^{±iα} ⟺ **z± = ±iα (+2πin)** — square-root branch points of B̃_α. Local form (EXACT algebra): B̃_α(z) ≈ γ₊(z−z₊)^{1/2}, **γ₊ = √(2 sinα)·e^{i(π/4−α/2)}**; at α = π/2, γ₊ = +√2, exactly D13b's local constant for their kernel.

**(c) Saddle (standard Hankel/saddle heuristics; constants verified against D13b at π/2).** [q^k]B̃_α(x(q−1)) = (1/2πi)∮B̃_α(z)(1+z/x)^{−N}dz/x → Hankel around z₊: ≈ γ₊(−i)/(2√π)·x^{1/2}N^{−3/2}e^{−Niα/x}; the outer Laplace integral ∫₀^∞ e^{−x}(·)dx has saddle x* = √(Niα), φ(x*) = 2√(Niα) = **√(2αN)(1+i)**, Gaussian factor √π x*. Assembling both conjugate branch points:
M_k(α) ≈ 2Re[γ₊(−i)e^{iπ/4}(√α/2)N^{−1}e^{−θ(1+i)}] = **√(2α sinα)·N^{−1}e^{−θ}cos(θ + α/2)**, θ = √(2αN).
The rate **c(α) = 2α** and amplitude **√(2α sinα)** are robust (at π/2: √π, β = 1 ✓). 

**(d) Phase failure (honest).** My assembly gives φ = −α/2 (or −α with the next log-term included), but the truth at π/2 is φ = π. The missing O(1) phase lives in the Hankel orientation/steepest-path details that D13b resolved for the squircle; I did **not** resolve it for general α. φ(α) = π for all α and φ(α) = π/2 + α are both consistent with φ(π/2) = π — undecided.

**(e) R is negligible (NUMERICAL).** R_k = ((−1)^k/k!)∫₀^∞[(x+L)^k/(1+e^{−x}) − x^k]dx = −0.693, 0.240, −0.0555, … decaying ≈ 0.6^k/k! (R_25 ≈ 1.8e‑31), so M_k = T_k to 30+ digits for k ≥ 25 and the stretched exponential is entirely T's.

**(f) Singularities of M(u;α) in u (the analogue picture).** From the expansion F = Σβ_jγ_m g^{j+(m−1)(1−u)}, the integrand's g-integral ceases to converge when j+(m−1)(1−u) = −1, i.e. **u = 1 + (j+1)/(m−1)** (j ≥ 0, m ≥ 2) — poles dense in the rationals > 1 **accumulating at u = 1** (plus u = −1,−2,…). So u = 1 is an essential singularity, radius of convergence of ΣM_k u^k is exactly 1 (consistent with |M_k|^{1/k} → 1 from the decay law), and the stretched exponential e^{−√(2αN)} is the shadow of this accumulation — exactly the "two ends" mechanism: the diamond end's branch points ±iπ/2 for the squircle become **±iα** here.

## 4. Evidence/labels summary

| claim | label | evidence |
|---|---|---|
| M_k tables, k ≤ 10, 5 angles, ≥ 59 digits | NUMERICAL | corner_series cross-check (independent method) |
| M_k(π/2), k ≤ 300, ≥ 21 digits | NUMERICAL | D13b published acb values, full agreement |
| M_k(α), k = 300, ~15 digits (all angles) | NUMERICAL | FLINT acb_series T_k + R_k cross-check |
| M_k(α), k ≤ 300, ≥ 40 digits target | NUMERICAL (likely, not fully certified) | measured e^{−2.42/h} rate, h = 0.011 |
| small-g expansion + subtraction closed form | PROVED (algebra) + NUMERICAL 7.8e‑62 | step_09.py |
| M = T + R representation | PROVED (algebra) + NUMERICAL 1e‑62 | step_23.py |
| branch points z± = ±iα, γ₊ = √(2 sinα)e^{i(π/4−α/2)} | PROVED (exact algebra) | step_23.py analysis |
| c(α) = 2α, β = 1, \|A(α)\| = √(2α sinα) | NUMERICAL (envelope, 5 angles) + saddle derivation | matches D13b at π/2 exactly |
| φ(π/2) = π; φ(α) general | NUMERICAL (π/2, from D13b) / CONJECTURE otherwise | phase not pinned — failure |
| LSQ fits of c, β, A, φ | **NOT DONE** | tables ready in `out/mk_*.tsv` |
| FLINT k = 0 sum leak (7.6e‑20 vs ball 3.9e‑26) | unresolved bug (mine, not quadrature's) | mpmath node-sum converges to 1.3e‑93 |

## 5. Files

Scripts (auto-saved `step_NN.py`, in run order): `step_01.py` (survey), `step_09.py` (subtraction algebra: exact to 7.8e‑62), `step_11.py` (grid bug: tanh covers (−1,1)), `step_13.py` (correct (0,1) DE map → exponential convergence), `step_15.py` (production v1, h = 0.025 — only ~12 digits at k = 300), `step_16.py` (Mval vs mp.quad), `step_17.py` (DFT self-test: exact), `step_18.py` (worst-p trapezoidal rate c ≈ 2.42), `step_19.py` (f-value save; the "max|Im|" scare was a meaningless diagnostic, not a bug), `step_21.py`/`step_25.py` (T+R derivation; e^{−z} sign fix), `step_23.py` (T+R exact to 1e‑62), `step_27.py`–`step_33.py` (node/weight transformation: w0·uh correct), then the production h = 0.011 run (all 5 angles, `out/mk_*.tsv`, `out/fval_*.json`), the FLINT T+R/T runs (`out/mk_tr_*.tsv`, `out/mk_trT_*.tsv`), the R_k integral values (k ≤ 25), and the single-node FLINT precision isolation (final cell: node values good to 4e‑42 at x = 1.5 — leak is in the sum, not the nodes).

**Bottom line:** the corner coefficients obey M_k(α) ≈ √(2α sinα)·N^{−1}e^{−√(2αN)}cos(√(2αN)+φ(α)): the squircle's √(πN) becomes √(2αN) because the diamond-end branch points ±iπ/2 move to ±iα — the two ends of the corner. Open: the LSQ fits, the φ(α) constant, full 40-digit end-to-end certification at k = 300, and the FLINT sum leak.