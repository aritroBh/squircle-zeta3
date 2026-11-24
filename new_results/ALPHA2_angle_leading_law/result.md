# ALPHA2: the leading-coefficient law for the corner-angle family M_k(alpha)

**Workdir:** `new_results/ALPHA2_angle_leading_law/`. Steps run under `nice -n 10`, one process at a time; `progress.log` records finished steps.

## Verdict

- **Task A (top table and angle factors): DONE.** The corner top table is a **single factor** times the squircle pattern at every depth checked: `beta_d^{(ij)}(rho; c) = F_d(c) * beta_d^{(ij)}(rho; squircle)` — identical ratio across all keys and poles (EXACT at 8 rational angles + squircle, `d <= 8`). Closed forms `F_d in Q(t)`, `t = tan^2(alpha/2)`, found for `d <= 8` (EXACT on 9 points each, with spare points); the pattern in `d` is CONJECTURE (sub-structure recorded).
- **Task B (assembly): DONE.** `c(w;alpha) = sum beta N` assembled at 50 digits at 8 angles, `w = 2..8`; matches the law below to 44-51 digits. (Note: the correct normalization is `c(w;alpha) = sum beta N`, **without** the `1/(d-1)!` of the task sheet — pinned exactly by `c(2) = T(alpha)`.)
- **Task C (generating function): DONE.** Corner measures derived (PROVED symbolically); corner Master Lemma established; the Mellin kernel is `(1+w)/(1-2w cos a + w^2)^{3/2}` and the generating function is an associated Legendre function:
  `E_a(s) = -pi P^1_{-1-s}(-cos a)/(sin a sin(pi s)) = Gamma(1-s) Gamma(2+s) {}_2F_1(1+s,-s;2;c^2)/(1-cos a)` (NUMERICAL 25-52 digits at 6 angles; reduces to `E(s)` at `a = pi/2` PROVED).
  **The law:** `c(w;alpha) = eps_{w-2} A_{w-1} F_{w-1}(c)/(2(w-2)!) * [y^w](E_a(y) - cos a * E_a(-y))` — NUMERICAL 44-51 digits, 8 angles x w=2..8; PROVED conditional on the top table (the LEAD2 proof goes through verbatim with the corner measures).

## 1. Setup and corner measures — PROVED (`step_26/27`, sympy)

With `c = cos(a/2)`, `r = (1-c)/(1+c)`, letters `{0,+-1,+-r}`, `t = (w-r)(w+1)/(2(1-r)w)`, set
`u := 2t-1 = (w^2-r)/((1-r)w) = kappa(z - 1/z)`, `z = w/sqrt(r)`, `kappa = sqrt(r)/(1-r) = tan(a/2)/2`.
Then `t = (1+u)/2` exactly, and with `S_k = sqrt(u^2+4k^2)`, `D = sqrt(1+4k^2) = 1/c`, `dz/z = du/S_k`,
the four pole measures `mu_rho = dz/(z-z0)` (PROVED symbolically, `step_27`):
```
mu_{w=1}  = (1/2)[1/S_k + 1/(u-1) + D/(S_k(u-1))] du
mu_{w=-1} = (1/2)[1/S_k + 1/(u+1) - D/(S_k(u+1))] du
mu_{w=-r} = (1/2)[1/S_k + 1/(u-1) - D/(S_k(u-1))] du
mu_{w=0}  = du/S_k
```
The LEAD2 pole combinations give the corner measures (PROVED, `step_28`):
```
m_odd^k  = mu_{-1} - mu_1   = [1 + D u/S_k]/(1-u^2)  = 1/(1-u) - u m_e^k
m_e^k    = mu_{-1} - mu_{-r} = [1 - D/S_k]/(1-u^2)   (even in u)
m_diag^k = mu_{-1} - mu_0/2 = (1-u)/2 m_e^k
m_oddlast^k = m_odd^k - 1/(2(1-u))
```
At `kappa = 1/2` they are exactly LEAD2's `m_odd, m_e, m_diag, m_oddlast`. The two structural identities
`m_odd = 1/(1-u) - u m_e` and `m_diag - m_e/2 = -u m_e/2` hold (sympy residual 0), so the LEAD2 collapse
(even d) and decomposition (odd d) apply verbatim.

## 2. Task A: the angle factors F_d — EXACT at 8 angles + squircle, d <= 8

Exact kernels computed with the ALPHA1 reducer for `c in {1/2,3/5,1/3,4/5,2/3,5/7,1/4,6/7}`, `k <= 18`
(`kern_*.pkl, kernx_*.pkl, kern18_*.pkl`; ~250 s per angle to k=16, ~200 s more to k=18). Residue
sequences polynomial in k of degree d-1 with spares (asserted). Ratio corner/squircle of the top
coefficients is a **single value F_d(c) per (c,d)** — no pole-dependent factors (`step_33/34/52/62`).

Closed forms (`t = tan^2(a/2) = (1-c^2)/c^2`), each verified on 9 points (8 angles + squircle `F_d(1)=1`)
with 2-4 spares; `f_d(t) := F_d * f_d(1)`:
```
F_2 = t                         f_2 = t                          f_2(1)=1
F_3 = (1+3t)/4                  f_3 = 1+3t                       f_3(1)=4
F_4 = t(2+3t)/5                 f_4 = 2t+3t^2                    f_4(1)=5
F_5 = (2+15t+15t^2)/32          f_5 = 2+15t+15t^2                f_5(1)=32
F_6 = t(17+60t+45t^2)/122       f_6 = 17t+60t^2+45t^3            f_6(1)=122
F_7 = (17+231t+525t^2+315t^3)/1088   f_7 = 17+231t+525t^2+315t^3  f_7(1)=1088
F_8 = t(62+378t+630t^2+315t^3)/1385  f_8 = 62t+378t^2+630t^3+315t^4 f_8(1)=1385
```
Equivalent `ca = cos a` form: **`F_d(ca) = P_d(ca)/[A_d (1+ca)^{floor(d/2)}]`** with `A_d` zigzag and
`P_d in Z[ca]`, `P_d(0) = A_d`, leading coeff `(-1)^{floor(d/2)}`:
```
P_2=1-x; P_3=2-x; P_4=5-6x+x^2; P_5=16-13x+x^2; P_6=61-89x+29x^2-x^3;
P_7=272-297x+60x^2-x^3; P_8=1385-2396x+1134x^2-124x^3+x^4
```
CONJECTURE-level sub-patterns: monic `h_d = f_d/L_d`, `L_d` = leading coeff; coefficient of `t^{m-j}`
(`m = floor(d/2)`) is `(d-2j) g_j(d)`: `g_1 = 1/3`, `g_2 = (5d-13)/90`, `g_3 = (14d-81)/315`; the
`x^{m-1}` coefficients of `P_{2m}` and `P_{2m+1}` are `(-1)^{m-1}(2^{2m-1}-m)` and `(-1)^{m-1}(4^m-m-1)`;
`P_d(-1) = d!/2^{floor((d-1)/2)}` (d<=8); `a(2m,1) = a(2m+1,0)` in `f_d = sum a(d,j) t^j`. No low-order
recurrence in d found (order-2/3/4 with low-degree polynomial coefficients all excluded).

## 3. Task C: corner Master Lemma and the Mellin kernel

Under `u = (1-w)/(1+w)` the corner measures map to `h_k(w) = [(1+w)/Q - 1]/w`, `Q = sqrt(1-2w cos a + w^2)`
(the Gegenbauer kernel `(1-2xw+w^2)^{-1} = sum C_n(x) w^n` enters here). PROVED by sympy (`step_64/65`):
`q_k(1/w) = w^2 q_k(w)` and `h_k(1/w) = w^2 h_k(w)` for `q_k = (1+w)/Q^3`; the primitives
`d/dw[1-(1+w)/Q] = (1+ca)(w-1)/Q^3` and `d/dw[(w-1)/((1-ca)Q)] = (1+w)/Q^3`.

Define on the strip `-1 < Re s < 1`: `E_a(s) := int_0^1 (w^{1+s}+w^{-s}) Q^{-3} dw`,
`WW(s) = int_0^inf w^{-s} q_k dw`, `FF(s) = int_0^inf w^{-s} h_k dw`, `GG(s) = int_0^inf w^{-s} ((1-w)/(1+w)) h_k dw`.

**Corner Master Lemma** (NUMERICAL 16-51 digits at 5-6 angles x 4-5 strip points, `step_29/30/41`):
```
(i)  WW(s) = E_a(s) + E_a(-s)            [PROVED: q-symmetry, LEAD2 Step-2 argument]
(ii) FF(s) = (1+ca)(E_a(s)-E_a(-s))/s    [PROVED: IBP with the primitive above, LEAD2 Step-3 argument]
(iii) GG(s) = 2pi/sin(pi s) - (1-ca) WW(s)/s   [LEAD2 Step-4 + Jfull formula below]
```
For (iii): `Jfull(s) := int_0^inf w^{s-1} Q^{-1} dw = pi P_{s-1}(-cos a)/sin(pi s)` (classical
Laplace-Mehler/GR 3.252.10 form; NUMERICAL 25+ digits at interior strip points, `step_44/46/47`; the
trig-linearized form is `Jfull(s) = int_0^{pi-a} sin^{s-1} u sin^{-s}(a+u) du`, reducing at `a=pi/2` to
LEAD2's `(1/2) B(s/2,(1-s)/2)`), then `GG = Jfull(-s) - Jfull(1-s) + 2pi/sin(pi s)` as in LEAD2 Step 4.

**The generating function** (`step_38/40/42`): with `I(mu) = int_0^inf w^{mu-1} Q^{-3} dw = (2/sin a) B(mu,3-mu) P^{-1}_{mu-2}(-cos a)`
(classical; NUMERICAL 44-51 digits, `step_37/39`), one gets `WW = I(1-s)+I(2-s)`, hence
```
E_a(s) = Gamma(1-s) Gamma(2+s) P^{-1}_{-1-s}(-cos a)/sin a
       = -pi P^1_{-1-s}(-cos a)/(sin a sin(pi s))
       = Gamma(1-s) Gamma(2+s) 2F1(1+s,-s;2;c^2)/(1-cos a)
```
NUMERICAL: matches the integral definition to 25-52 digits at 6 angles x 4 s-values. At `a = pi/2`:
`E_a(s) = Gamma(1+s/2)Gamma(1/2-s/2)/sqrt(pi) = E(s)` (PROVED by Legendre duplication, `step_38/45`).

## 4. The leading law

**Theorem (conditional on the top table `D_d^{(ij)}(rho;c) = F_d(c) * [squircle pattern]`).**
For every `w >= 2`, with `eps_m = (-1)^floor((m+1)/2)`, `A_n` zigzag, `F_d` as above,
```
c(w;alpha) = eps_{w-2} A_{w-1} F_{w-1}(c) / (2 (w-2)!)  *  [y^w](E_a(y) - cos a * E_a(-y))
```
equivalently `c(w;alpha) = eps_{w-2} A_{w-1} F_{w-1}(c) trig_w e^a_w/(w-2)!` with
`trig_w = sin^2(a/2)` (w even), `cos^2(a/2)` (w odd), `e^a_w = [y^w]E_a(y)`.

*Proof sketch (LEAD2 goes through verbatim).* The D-table factorizes as `F_d(c)` times the squircle
pattern, so `K_d(a) := c(d+1;a)(d-1)!/(delta_d F_d)` is given by the same measure combinations with the
corner measures. Even d: parity collapse (m_e^k even, m_diag^k - m_e^k/2 = -u m_e^k/2 odd) gives
`K_d/d! = (1/4)[z^d] FF_k(z) = c^2 e^a_{d+1}` by (ii). Odd d: the decomposition `m_odd^k = 1/(1-u) - u m_e^k`
gives `K_d = T_d - (u m_e^k part)` with `T_d/d! = (1-2^{-d}) zeta(d+1)` (kappa-independent, the A-bracket
involution of LEAD2); by (iii) `(u m_e^k part)/d! = (1/4)[z^d] GG_k(z) = (1/2)[z^d](pi/sin pi z) - (1-ca)e^a_{d+1}/2`;
the `pi/sin` terms cancel exactly (the same obstruction as in LEAD2), leaving `K_d/d! = s^2 e^a_{d+1}`.
Assembly `c(d+1;a) = delta_d F_d d [trig] e^a_{d+1}` gives the law. []

**Verification:** 8 angles x `w = 2..7` at 44-51 digits (`step_49/50`), `w = 8` at 2 angles, 44-45 digits
(`step_72`). Squircle `a = pi/2` recovers LEAD2 exactly (`trig = 1/2`, `E_a = E`, `F = 1`).
Cross-checks independent of the F-table: `e^a_2 = 2 T(a)/(1-ca)` (T = pi^2/12 - Li_2(c^2)/2) and
`e^a_3 = -J(c)/(2c^2)` with ALPHA1's exact `J(c)`, both 61 digits at 8 angles (`step_51/53/54`);
`e^a_1 = -ln sin(a/2)/cos^2(a/2)` 61 digits (`step_45`).

**Task B normalization note.** `c(w;alpha) = sum_{keys,rho} beta_d^{(ij)}(rho) N_{(i,j),rho}` with
`beta = [k^{d-1}] res` (no `1/(d-1)!`); `c(2;a) = T(a)` exactly (50 digits, `step_48`) pins this.

## 5. Open pieces

- Closed form/recurrence for `F_d` in d (data through d=8 above; sub-patterns in section 2). This is a
  tower-level question (top coefficients of the k-polynomials), sibling of TOW2 for the squircle.
- Analytic proof of the `Jfull` Legendre evaluation along LEAD2 lines (currently classical + 25-digit
  numerics); likewise a from-scratch proof of the `I(mu)` 3/2-power Mellin formula.
- `d >= 9` top tables (k <= 20 kernels) for sharper `F_d` pattern fitting.

## Evidence table

| Claim | Type | Evidence |
|---|---|---|
| 4 pole measures + combinations (m_e, m_odd, m_diag, m_oddlast) | PROVED | sympy residual 0 (`step_27/28`) |
| `m_odd = 1/(1-u) - u m_e`; `m_diag - m_e/2 = -u m_e/2` | PROVED | sympy residual 0 (`step_28`) |
| `q,h` w<->1/w symmetries; 2 IBP primitives | PROVED | sympy residual 0 (`step_64/65`) |
| Single factor F_d per (c,d) | EXACT, d<=8 | 8 angles+squircle, 0 ratio failures (`step_33/52/62`) |
| F_d closed forms in t, d<=8 | EXACT on 9 pts each | 2-4 spare points (`step_35/52/62`) |
| `F_d = P_d/[A_d(1+ca)^f]` | EXACT on data, CONJECTURE in d | d<=8 (`step_63/75`) |
| Master Lemma (i)(ii)(iii) | (i)(ii) PROVED, (iii) NUMERICAL 16-51 d | 5-6 angles x 4-5 strip pts (`step_29/30/41`) |
| `E_a` Legendre/2F1 closed form | NUMERICAL 25-52 d | 6 angles x 4 s (`step_38/40/42`) |
| `E_a` reduces to squircle E at pi/2 | PROVED + NUMERICAL exact | duplication + Taylor (`step_45`) |
| `Jfull = pi P_{s-1}(-ca)/sin(ps)` | NUMERICAL 25+ d | interior strip pts (`step_44/46`) |
| Leading law `c(w;a)` | NUMERICAL 44-51 d | 8 angles x w=2..7; w=8 x2 angles (`step_49/50/72`) |
| `e^a_2 = 2T/(1-ca)`, `e^a_3 = -J/(2c^2)` | NUMERICAL 61 d | 8 angles (`step_51/53/54`) |
| `c(2;a) = T(a)` assembly normalization | NUMERICAL 50 d | `step_48` |
| F_d d-pattern (g_j, x^{m-1} forms) | CONJECTURE | d<=8 data (`step_67/75`) |

## Files
- `assemble.py` (beta extraction + N-integrals + assembly), `Fd_data.pkl`, `Fd_fits.pkl`, `lawcheck.pkl`.
- `kern_*.pkl` (k<=16, new angles), `kernx_*.pkl` (k<=16, old), `kern18_*.pkl` (k<=18, all 8).
- Reducer (copied from ALPHA1): `cwfield.py, cwintegrand.py, cwred.py, cwbuild.py, lllkernel.py, hpl_tools.py, corner_series.py`.
- `verify_alpha2.py` (PASS/FAIL runner, < 15 min).
