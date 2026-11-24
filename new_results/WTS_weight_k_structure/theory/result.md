# WTS/theory: the weight-(d+1) kernels of M_k are polynomial in k (degree d-1) -- result

Folder: `new_results/WTS_weight_k_structure/theory/`.  Python: project venv.  Everything below that is called
"exact" was printed by a script in this folder (outputs `*.out`, data `*.pkl`/`*.json`).  Nothing was written outside this folder.

## VERDICT

**EXACT DATA ONLY for the general law, plus a proved reduction to a sharper, still unproved, statement about the tower at v = +-i.
No all-d proof (and no new proof for any single d >= 2) was found.**

| item | status |
|---|---|
| Law K(d): degree-d kernels are polynomials in k of degree <= d-1 for k >= d+1 | **exact data**: every key, every pole, d = 1..8, k <= 18 (Sec. 1). Degree is *exactly* d-1 for every key; k = d is *never* on the polynomial. d = 1 is proved for all k (WT2b/CX4/CX5, not my work). |
| Equivalent forms (generating function = eps^d P(x), x = eps/(1-eps), deg P <= d, etc.) | **proved in writing** (Lemma 1) |
| Kernel is a k-independent linear image of the degree-(d+1) tower coefficients | **proved by reading the algorithm + exact check k <= 14** (Lemma 2) |
| NEW tower law T: the coefficient of depth delta (pole order r = 2n-2-delta at v = +i, degree n) is a polynomial in k of degree exactly floor(delta/2) | **exact data**: 64 independent (n, delta) cells, n = 2..11, no exception; delta = 0 is the already proved Lemma D (WT2b/CX5) (Sec. 2.3) |
| T(n = d+1) => K(d) (for k >= d+2; the single extra point k = d+1 is only checked) | **proved in writing** (Lemma 3) |
| Why: x = 1/(p-1) is the natural variable; the 1/(1-eps) is the factor of (d/dsigma + eps d/dtheta) on functions of sigma - theta | **explanation, partly machine-checked, not a proof** (Sec. 2.4) |
| Proof of T for delta >= 1 (hence of K(d) for d >= 2) | **NOT obtained**. Precise obstruction in Sec. 3. |

k-range: K(d) holds for k >= d+1.  At level k = d the kernel is the "initial" one and is off the polynomial for every tested d
(d = 1: k = 1 is the exceptional R10; d = 2: k = 2 exceptional, as in WT3).

## 0. Scripts (all run with `nice -n 10`, one at a time, each < 20 min) and outputs

| file | purpose |
|---|---|
| `run_kernel_k.py K`, `run_all_kernels.sh`, `run_more_k.sh` | exact fast-reducer kernels of ALL keys at level K -> `kern_K.pkl` (K = 1..18; logs `kernels.log`, `kernels_more.log`) |
| `analyze_diffs.py` -> `analyze_diffs.out` | exact finite-difference test of K(d) |
| `gen_functions.py` -> `gen_functions_d8.out`, `gen_functions_d3.out` | Newton-form generating functions, re-expanded and compared with the data for every k |
| `make_tables.py` -> `genfun.json`, `genfun_d1_d3.md` | exact Newton data for all keys d <= 8; readable tables d = 1,2,3 |
| `tower.py` -> `tower_K.pkl` (K = 2..14), `tower.log` | same reducer, recording the +i higher part (the "tower") of every work item |
| `analyze_tower.py` -> `analyze_tower.out`, `summarize_tower.py` -> `summarize_tower.out` | depth law T |
| `check_tower_to_kernel.py` -> `check_tower_to_kernel.out` | Lemma 2 |
| `check_top.py` -> `check_top.out` | leading coefficient = zigzag numbers |
| `unfolded_test.py`, `analyze_unfolded.py` -> `unfolded_test.out`, `analyze_unfolded.out` | side experiment: reduction WITHOUT the v -> 1/v fold |
| `verify_dual.py` -> `verify_dual.out` | 30-digit check of the dual-exponent (Cauchy) forms |
| `verify_single_log.py` -> `verify_single_log.out` | exact check of the single-log form of the integrand |
| `verify_branch_form.py` -> `verify_branch_form.out` | 40-digit check of the sigma-coordinate algebra of Sec. 2.4 |
| `check_d1.py` -> `check_d1.out` | d = 1 against WT2 |
| `fastred.py`, `genintegrand.py`, `buildintegrand.py`, `binteg_K.pkl` | copied unchanged from WT2 (binteg_16..18 built here) |

## 1. Task 1 -- exact data

Setting: `K^{(k)}_{ij}(v) = sum_r c_r(k)/(v-r)`, r in {0, a, -a, b, -b}; a = 1+sqrt2, b = sqrt2-1 (a pole at v = +b never occurs).
`c_r(k)` lies in Q(sqrt2).  Levels k = 1..18 were computed (k = 16, 17, 18 took 174 s, 273 s, 336 s of reduction).  In every
kernel no pole at +-i and no constant part survives (asserted in `analyze_diffs.py`).

**Test of the conjecture** (`analyze_diffs.out`): for every key (i,j), d = i+j, and every pole, both Q(sqrt2)-components,
the d-th forward difference of k -> c_r(k), k = d+1..18, is identically 0, and the minimal vanishing order shows the degree is
exactly d-1.  Number of independent d-th differences available (per sequence):

| d | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| independent tests | 16 | 14 | 12 | 10 | 8 | 6 | 4 | 2 |
| keys | (1,0) | (2,0),(1,1) | (3,0),(2,1) | (4,0),(3,1),(2,2) | 3 keys | 4 | 4 | 5 |

All keys of d <= 8 PASS (full list in `analyze_diffs.out`); for d >= 9 there are not enough levels (d = 9 would need k = 19).
(d = 1..7 had already passed with k <= 15; levels 16..18 were run afterwards and everything still passes, and d = 8 was added.)

**Generating functions.**  Put x = eps/(1-eps).  Writing `D_m = Delta^{m-1} c_r (k = d+1)` (Newton coefficients) one has, exactly,
`sum_k c_r(k) eps^k = c_r(d) eps^d + sum_{m=1}^d D_m eps^{d+m}/(1-eps)^m = eps^d [ c_r(d) + D_1 x + ... + D_d x^d ]`.
`gen_functions.py` re-expands these rational functions and checks them against ALL levels k = 1..18 (including k = d and the
vanishing below d): "ALL GENERATING FUNCTIONS REPRODUCE THE DATA EXACTLY ... True" for d <= 8 (`gen_functions_d8.out`).
`genfun.json` holds `c_d` and `D_1..D_d` for every key/pole, d <= 8.

**d = 1 (check against WT2):** `check_d1.py` -> `check_d1.out`: c(1) = residues of R10 and D_1 = residues of Rw2 at all poles, so the generating function is exactly `eps*R10 + eps^2/(1-eps)*Rw2` (table below).

**d = 2 and d = 3, all keys and poles (machine-generated; `numerator/(1-eps)^m`, expanded):**


**key (1,0), d=1**  (poles v = 0, a, -a, b, -b with a=1+sqrt2, b=sqrt2-1; no pole at +b ever occurs)

| pole | G_r(eps) = numerator/(1-eps)^m (expanded; sqrt2 explicit) | eps^d[c_d + D_1 x + ... + D_d x^d] data (c_d; D_1..D_d) |
|---|---|---|
| -a | [eps**2 - 2*eps] / (1-eps)^1 | c_d=-2; D=['-1'] |
| 0 | [eps**2*(-1/2 + sqrt(2)/2) + eps*(1 - sqrt(2)/2)] / (1-eps)^1 | c_d=1 - sqrt(2)/2; D=['1/2'] |
| -b | [-eps**2/2] / (1-eps)^1 | c_d=0; D=['-1/2'] |
| a | [eps**2/2] / (1-eps)^1 | c_d=0; D=['1/2'] |

**key (2,0), d=2**  (poles v = 0, a, -a, b, -b with a=1+sqrt2, b=sqrt2-1; no pole at +b ever occurs)

| pole | G_r(eps) = numerator/(1-eps)^m (expanded; sqrt2 explicit) | eps^d[c_d + D_1 x + ... + D_d x^d] data (c_d; D_1..D_d) |
|---|---|---|
| -a | [-eps**4/4 + 7*eps**3/4 - eps**2] / (1-eps)^2 | c_d=-1; D=['-1/4', '1/2'] |
| 0 | [eps**3*(-1/4 + sqrt(2)/8) + eps**2*(1/2 - 3*sqrt(2)/8)] / (1-eps)^1 | c_d=1/2 - 3*sqrt(2)/8; D=['1/4 - sqrt(2)/4', '0'] |
| -b | [-eps**3/2] / (1-eps)^2 | c_d=0; D=['-1/2', '-1/2'] |
| a | [eps**3/4] / (1-eps)^1 | c_d=0; D=['1/4', '0'] |

**key (1,1), d=2**  (poles v = 0, a, -a, b, -b with a=1+sqrt2, b=sqrt2-1; no pole at +b ever occurs)

| pole | G_r(eps) = numerator/(1-eps)^m (expanded; sqrt2 explicit) | eps^d[c_d + D_1 x + ... + D_d x^d] data (c_d; D_1..D_d) |
|---|---|---|
| -a | [-eps**4/2 - eps**3/2] / (1-eps)^2 | c_d=0; D=['-1/2', '-1'] |
| 0 | [eps**4*(1/4 - sqrt(2)/8) + eps**3/4 + sqrt(2)*eps**2/8] / (1-eps)^2 | c_d=sqrt(2)/8; D=['1/4 + sqrt(2)/4', '1/2'] |

**key (3,0), d=3**  (poles v = 0, a, -a, b, -b with a=1+sqrt2, b=sqrt2-1; no pole at +b ever occurs)

| pole | G_r(eps) = numerator/(1-eps)^m (expanded; sqrt2 explicit) | eps^d[c_d + D_1 x + ... + D_d x^d] data (c_d; D_1..D_d) |
|---|---|---|
| -a | [eps**6/24 - 7*eps**5/12 + 29*eps**4/24 - eps**3/3] / (1-eps)^3 | c_d=-1/3; D=['5/24', '5/6', '1/3'] |
| 0 | [eps**5*(1/24 - sqrt(2)/48) + eps**4*(-7/24 + sqrt(2)/6) + eps**3*(1/6 - 7*sqrt(2)/48)] / (1-eps)^2 | c_d=1/6 - 7*sqrt(2)/48; D=['1/24 - sqrt(2)/8', '-1/12', '0'] |
| -b | [-eps**4/6] / (1-eps)^2 | c_d=0; D=['-1/6', '-1/6', '0'] |
| a | [eps**6/24 - eps**5/4 - eps**4/8] / (1-eps)^3 | c_d=0; D=['-1/8', '-1/2', '-1/3'] |

**key (2,1), d=3**  (poles v = 0, a, -a, b, -b with a=1+sqrt2, b=sqrt2-1; no pole at +b ever occurs)

| pole | G_r(eps) = numerator/(1-eps)^m (expanded; sqrt2 explicit) | eps^d[c_d + D_1 x + ... + D_d x^d] data (c_d; D_1..D_d) |
|---|---|---|
| -a | [eps**6/8 - eps**5/4 - 7*eps**4/8] / (1-eps)^3 | c_d=0; D=['-7/8', '-2', '-1'] |
| 0 | [eps**6*(-1/16 + sqrt(2)/16) + eps**5*(1/8 - sqrt(2)/16) + eps**4*(7/16 - sqrt(2)/16) + sqrt(2)*eps**3/16] / (1-eps)^3 | c_d=sqrt(2)/16; D=['sqrt(2)/8 + 7/16', '1', '1/2'] |
| -b | [-eps**6/16 - eps**5/8 - 5*eps**4/16] / (1-eps)^3 | c_d=0; D=['-5/16', '-3/4', '-1/2'] |
| a | [eps**6/16 + eps**5/8 + 5*eps**4/16] / (1-eps)^3 | c_d=0; D=['5/16', '3/4', '1/2'] |

Per-pole degrees: not every pole carries the full degree d-1 (e.g. d = 4, key (4,0): the poles v = 0 and v = a have degree 2, the
poles -a, -b degree 3; see the D-lists in `genfun.json`); the *key* always has degree exactly d-1.

**Leading coefficient (new exact regularity, `check_top.out`).**  Let `D_d = Delta^{d-1} c(d+1) = (d-1)! [k^{d-1}] c(k)`.  At the pole
v = -a, for every key (i,j) with i+j = d, d = 2..8 (exact):

`D_d(i,j; -a) = (-1)^{j+1} binom(d,j) C_d`, with `C_d = [x^d]( sech x + tanh x ) = (-1)^{floor(d/2)} A_d / d!`,

A_d = Euler zigzag (up/down) numbers 1, 2, 5, 16, 61, 272, 1385 for d = 2..8.  (`check_top.py` computes the sympy series of
(1+sinh x)/cosh x and prints "equal: True" for each d.)  The other poles carry the same C_d with simple binomial/parity patterns
(printed in `check_top.out`).  Any proof must reproduce this p -> 1 (eps -> 1) leading behaviour.

## 2. Task 2 -- structure

### 2.1 Lemma 1 (equivalent forms of the law).  PROVED (elementary)

Let `c(k)`, k >= d, be a sequence in a Q-vector space and `S(eps) = sum_{k>=d} c(k) eps^k`, x = eps/(1-eps).  The following are equivalent:
(a) c(k) = p(k) for k >= d+1 with p a polynomial of degree <= d-1;
(b) S = eps^d [ c(d) + sum_{m=1}^d D_m x^m ], D_m = Delta^{m-1} c(d+1);
(c) S = eps^d R(eps)/(1-eps)^d with R a polynomial of degree <= d.
Proof.  [eps^k] eps^d x^m = binom(k-d-1, m-1) for k >= d+1 (m >= 1), a polynomial in k of degree m-1, while m = 0 only contributes at
k = d; Newton's forward-difference formula gives (a)<=>(b).  eps^d x^m = eps^d (1-eps)^{-d} [eps^m (1-eps)^{d-m}] and the Bernstein
polynomials eps^m(1-eps)^{d-m}, 0 <= m <= d, span the polynomials of degree <= d: (b)<=>(c).  QED.
(The formulation "rational with only pole at eps = 1, of order <= d" gives the polynomial law only for k >= k0, k0 depending on the numerator degree;
the sharp range k >= d+1 is exactly the numerator bound deg R <= d, i.e. numerator of S of degree <= 2d and divisible by eps^d.  For the data both hold.)

### 2.2 Lemma 2 (tower -> kernel).  PROVED by inspection of the algorithm, checked exactly for k <= 14

Notation (WT2/CX5): `W_ij^{(k)}` work items, `H` = principal part of order >= 2 at +-i (the "tower"; the -i part is the complex
conjugate), `Gm_ij = (1 - eval_b) o Int H(W_ij)`, `h_ij^{(k)}(r)` = coefficient of (v-i)^{-r}, 2 <= r <= 2n-2, n = i+j (Lemma A).
For k > d = i+j,
`W_ij^{(k)} = -(i+1) Gm_{i+1,j} DL1 - [i>j](j+1) Gm_{i,j+1} DL2 - [j=i-1] i F(Gm_ii DL1)`  (no initial term: it lives only at k = d),
the kernel `K^{(k)}_{ij}` is the simple part of `W_ij^{(k)}` and the tower of key (i,j) is its higher part.  Since `Gm(r) = [sum_{m>=1} e_m (r-i)^{-m} + conj] - [same at b]`, `e_m = -h(m+1)/m`, and the residue of `Gm*DL` at a real letter r is `Gm(r) * Res_r(DL)`:
**every real residue of K^{(k)}_{ij} (and the +-i simple-pole coefficient, which is 0) is a Q(sqrt2,i)-linear combination, with
k-independent coefficients, of the coefficients h_{i+1,j}^{(k)}(.), h_{i,j+1}^{(k)}(.), h_{ii}^{(k)}(.) of degree d+1; and the tower of degree d is the
same kind of k-independent linear image (the operator Lambda_d; it is R-linear, not C-linear, because the fold and the -i conjugate part enter).**  Order bookkeeping: order r' at degree d+1 feeds order r = r'-1-j', j' >= 1
(integrate: -1; DL has a simple zero at +-i: at least another -1), i.e. depth `delta = 2n-2-r` is **non-decreasing** under Lambda
(delta_new = delta_old + (j'-1)); the depth-preserving part is the Lemma-B recursion of WT2b/CX4.
`check_tower_to_kernel.py` verifies, for every level k and every key of degree 1..k-1, that the real-pole part of the expression above, built from the
recorded tower of the parents, equals the saved kernel, that its +-i simple-pole part vanishes, and that its higher part equals the recorded tower of the key
(output below).

### 2.3 Tower law T (new exact data) and Lemma 3

**Conjecture T (tower depth law).**  For n >= 2, every key (a,b) with a+b = n, and every 2 <= r <= 2n-2 (depth delta = 2n-2-r), the coefficient
`h_{ab}^{(k)}(r)` in Q(sqrt2,i) is, for k >= n (k >= n+1 for the diagonal key (n/2,n/2)), a polynomial in k of degree exactly floor(delta/2).

Evidence (`analyze_tower.out`, `summarize_tower.out`; the towers were recorded for k <= 14, a cell counts as a genuine test only if it
has at least degree+2 levels): 64 (n,delta)-cells tested, **0 exceptions**; for n = 3..7 every depth 0..2n-4 is tested, for n = 8..11 depths up to 9, 7, 5, 3.
depth 0 (stabilisation of the top coefficient) and the diagonal shift at k = 2p+1 are the proved Lemma D; depth 1 also has degree 0,
depth 2,3 degree 1, depth 4,5 degree 2, ... exactly as predicted.  (The diagonal exceptions at k = n for even n are visible as "[False, True]" in the
column "k=n on polynomial".)

**Lemma 3 (T => K).  PROVED (given Lemma 2).**  Let d >= 1, n = d+1.  If T holds for the keys of degree n for all k >= n+1, then for k >= d+2
every real residue of every degree-d kernel is a polynomial in k of degree <= d-1 (and the +-i simple-pole coefficient, also a polynomial of degree <= d-1,
vanishes identically as soon as it vanishes at d values of k >= d+2, which the data show for d <= 8).
Proof.  By Lemma 2 the residues are fixed Q(sqrt2,i)-linear combinations of h^{(k)}(r), 2 <= r <= 2d, of degree n = d+1, i.e. depths
delta = 2d - r in [0, 2d-2]; by T their k-degrees are <= floor((2d-2)/2) = d-1.  Reality of the residues makes the coefficients real.  QED.
The exact degree d-1 of the data is exactly the maximal depth 2d-2 of the degree-(d+1) tower: this is the quantitative content of "lowering the
log degree by one costs at most one factor 1/(1-eps)" (one factor per pair of pole orders).
The extra level k = d+1 is outside the stated range of T only for the diagonal key of even n; it is covered by the data (Sec. 1), not proved.

**Lemma 4 (T is a finite-level statement).  PROVED (reformulation).**  Put `u_n^{(k)} := Lambda_n ... Lambda_{k-1} rho_k` (the tower of degree n at level k; rho = initial data), and
`nu^{(0)}_m = rho_m`, `nu^{(j+1)}_m = Lambda_m nu^{(j)}_{m+1} - nu^{(j)}_m`.  Then `Delta_k^j u_n^{(k)} = Lambda_n ... Lambda_{k-1} nu^{(j)}_k`.  Since Lambda is depth
non-decreasing, T (all k >= n, degrees >= n, depths <= D; for the diagonal keys of even degree the range is shifted by one level) follows from: for every m >= n, the depth-delta component of nu^{(j)}_m vanishes for
all j > delta/2, delta <= D.  Each such component involves only the levels m..m+j.  So T for fixed depth is a statement about the first
floor(delta/2)+1 steps of the iteration, for every n; a proof for all n needs closed forms of rho_m up to depth delta (Lemma C of WT2b gives
delta = 0 only).  This is what a proof must supply (Sec. 3).

### 2.4 Why x = eps/(1-eps): what is verified and what is only a sketch

(i) **Dual exponent (verified, `verify_dual.out`).**  With p = 1/eps and p' = p/(p-1) = 1/(1-eps) = 1+x, x = 1/(p-1):
Cauchy's formula (perimeter = integral of the support function, which for the l_p ball is the l_{p'} norm) gives
`P(p)/8 = int_0^1 (1+u^{p'})^{1/p'} (1+u^2)^{-3/2} du`; and
`M(eps) = (1+x) int_0^infty y^{x-1} (1+y^{1+x})^{-1-eps} f(y) dy`, `f(y) = (sqrt(1+y^2) - 1 - y)/2`, with f independent of eps.
Numerically (30 digits, p = 3, 4, 5/2): direct arc length vs Cauchy vs F-form agree to 8e-31, and M(eps) vs the y-form to <= 5e-32.
In the y-form the sqrt branch points y = +-i are eps-independent and the only eps-dependent singularity is the zero y_* = -e^{-i pi eps}
of 1+y^{1+x}.  So x (not eps) is the natural variable of the whole family.
(ii) **Single log (verified exactly, `verify_single_log.out`, k = 1..7).**  Writing w = (1-t)/t, s = l2-l1, y = w^{-1}e^{eps s}:
`integrand = (1/2) e^{eps l1} (1/t) [ sqrt(1+y^2) - 1 - y ]`.  The two blocks A, B and the sqrt depend on the formal logs only through the
pure exponential e^{eps l1} (never needs integration by parts) and the single log s.
(iii) **Moving branch point (exact algebra).**  In the shifted family l_i -> L_i + lambda_i(v) the sqrt radicand is
`e^{2X1}(1-t)^2 [ 1 - e^{-2((1-eps) sigma + eps D)} ]`, sigma = log(i w) (it vanishes to second order at v = i, because v = i is a critical point of t(v): sigma ~ c0 (v-i)^2), D = L1-L2.  The branch points solve `sigma = -x D`: they depend on (eps, D) only through the product x D
(for the true logs D = i pi/2, so sigma_* = -i pi x/2).
(iv) **Sketch of the mechanism (hand derivation, NOT machine-checked, not a proof).**  In the coordinate sigma the sqrt block is
`(1/2) e^{eps l1} l1_sigma Phi0(sigma - theta) d sigma`, `Phi0(z) = sqrt(1 - e^{-2z})`, `theta = eps s`, `l1 = -log(1 - i e^sigma)`
(so `l1_sigma = i e^sigma/(1 - i e^sigma)`, whose Taylor coefficients are the polylog values Li_{-n}(i)/n! -- Euler/tangent-number-like, cf. the exact
C_d above; not verified here).  The two coordinate facts t = 1/(1 - i e^sigma), l1_sigma and the radicand identity 1 + w^{-2} e^{2 eps s} = 1 - exp(-2((1-eps) sigma + eps D)) are
checked numerically (40 digits, `verify_branch_form.out`).  The reduction equation for the generating tower G(sigma,theta) (exponents in Z+1/2, "polar part" = exponents <= -3/2) is (hand derivation)
`(d_sigma + eps (d_theta + l1_sigma)) G = Polar[ (1/2) l1_sigma Phi0(sigma - theta) ]`.
Parity in rho = sqrt(sigma) kills every simple pole automatically (explains "no pole at +-i survives").  On functions of z = sigma - theta,
`(d_sigma + eps d_theta) f(z) = (1-eps) f'(z)`: **each depth-preserving integration costs 1/(1-eps) = 1+x, each d_theta costs eps, i.e. x per
theta-derivative, and each multiplication by l1_sigma costs eps and moves one unit of sigma = two orders of the tower**, which is exactly the T
pattern "one more power of x per two orders of depth".  In (z,theta) the equation reads `(d_z + x(d_theta + l1_sigma)) G = (1+x) Polar'[...]`.
The step I could not carry out: the polar-part projection is defined at fixed sigma and does not commute with the change z = sigma - theta,
so the (z,theta) form only holds modulo a "regular" subspace that itself depends on theta; controlling that correction at all depths is the missing lemma.

### 2.5 Side results

* **The law does not depend on the fold.**  `unfolded_test.py` reduces the UNFOLDED integrand (all keys (i,j), no v -> 1/v, no symmetrisation; +i and -i towers then
  decouple).  Its kernels are *different* from the folded ones (after folding they disagree at every degree < k, `unfolded_test.out`), yet
  `analyze_unfolded.out` shows the same law: degree exactly d-1 for k >= d+1, d <= 5 (levels k <= 12; d = 5 has 2 tests).  So the mechanism is the local tower at
  v = +i, not the v -> 1/v coupling; a proof can be done in the decoupled (simpler) local problem.  (The folded kernels are the ones of the stated conjecture; whether the folded
  and unfolded kernels are related by an exact form was not investigated.)
* **No pole at +-i** survives in the unfolded reduction either (asserted in the script for k <= 12).
* k = d is exceptional for all d <= 8 (not on the polynomial) -- the "initial" kernel.

## 3. What is proved, what is checked, what is missing

| statement | status |
|---|---|
| Lemma 1 (equivalent forms), Lemma 3 (T => K for k >= d+2), Lemma 4 (T finite-level) | proved in writing |
| Lemma 2 (tower -> kernel is k-independent linear) | proved by reading the algorithm; exact check k <= 14 |
| K(d), d = 1..8, k <= 18, all keys/poles, degree exactly d-1, k = d exceptional | exact data (d = 1 also proved, WT2b/CX5) |
| T: degree floor(delta/2) at depth delta | exact data, 64 cells, no exception; delta = 0 proved (Lemma D) |
| top coefficient = (-1)^{j+1} binom(d,j) [x^d](sech + tanh) at v = -a | exact data d = 2..8 |
| law also holds for the unfolded reduction | exact data d <= 5, k <= 12 |
| T for delta >= 1; therefore K(d), d >= 2 | **OPEN** |

**Obstruction (precise).**  The proof of K(1) used only the top Laurent coefficient at v = i (depth 0) because degree-2 towers have order <= 2.  K(2) needs depth <= 2
of the degree-3 tower, hence (Lambda is depth non-decreasing) depth <= 2 of the towers of ALL degrees n >= 3 and of the initial data rho_n for all n; K(d) needs depth <= 2d-2.
The required inputs are (1) closed forms of the initial Laurent coefficients rho_n at depth 1, 2, ... (they involve the Taylor coefficients of t(v), z(v) at v = i beyond the leading
one used in Lemma C; in the sigma-coordinate of Sec. 2.4 they are products of the Taylor coefficients of l1_sigma and of the expansion of Phi0, i.e. explicit), (2) the matrix of Lambda
on depth <= 2d-2 (also explicit: Taylor coefficients of DL1, DL2 at i), (3) the identities of Lemma 4 for these data (the analogues of (I1)-(I4) of WT2b), for all n.
Everything is explicit but I did not find the closed forms; the identities in (3) (equivalently T) hold in the tested cells (n <= 11, depths as listed).
A fully rigorous alternative would be to prove the missing commutation statement in Sec. 2.4 (iv).

## 4. Recommended next steps

1. Prove depth 1 and 2 of T (this gives K(2) = the weight-3 law proved for all k, independent of the WT3 job), using the sigma-model: needs only Taylor coefficients kappa_0, kappa_1, kappa_2 of l1_sigma and the first three terms of Phi0.
2. Verify T for large n (n up to 40) at small depth by computing only the top few orders of the towers (cheap, Lemma 4: only levels m..m+floor(delta/2)+1 are needed) -- this
   would push the exact evidence for K(2), K(3) to all k <= 40 without running the full reducer.
3. Prove the top-coefficient formula (zigzag numbers) from the sigma-model; it is the leading singularity at eps = 1 (p -> 1).

## 5. Raw outputs (verbatim)

### analyze_diffs.out (tested keys; (key) (verdict, minimal degree, d-1, #levels k>=d+1, #independent tests, k=d on polynomial?))
```
(1, 0) ('PASS', 0, 0, 17, 16, 'k=d on polynomial: False')
(2, 0) ('PASS', 1, 1, 16, 14, 'k=d on polynomial: False')
(1, 1) ('PASS', 1, 1, 16, 14, 'k=d on polynomial: False')
(3, 0) ('PASS', 2, 2, 15, 12, 'k=d on polynomial: False')
(2, 1) ('PASS', 2, 2, 15, 12, 'k=d on polynomial: False')
(4, 0) ('PASS', 3, 3, 14, 10, 'k=d on polynomial: False')
(3, 1) ('PASS', 3, 3, 14, 10, 'k=d on polynomial: False')
(2, 2) ('PASS', 3, 3, 14, 10, 'k=d on polynomial: False')
(5, 0) ('PASS', 4, 4, 13, 8, 'k=d on polynomial: False')
(4, 1) ('PASS', 4, 4, 13, 8, 'k=d on polynomial: False')
(3, 2) ('PASS', 4, 4, 13, 8, 'k=d on polynomial: False')
(6, 0) ('PASS', 5, 5, 12, 6, 'k=d on polynomial: False')
(5, 1) ('PASS', 5, 5, 12, 6, 'k=d on polynomial: False')
(4, 2) ('PASS', 5, 5, 12, 6, 'k=d on polynomial: False')
(3, 3) ('PASS', 5, 5, 12, 6, 'k=d on polynomial: False')
(7, 0) ('PASS', 6, 6, 11, 4, 'k=d on polynomial: False')
(6, 1) ('PASS', 6, 6, 11, 4, 'k=d on polynomial: False')
(5, 2) ('PASS', 6, 6, 11, 4, 'k=d on polynomial: False')
(4, 3) ('PASS', 6, 6, 11, 4, 'k=d on polynomial: False')
(8, 0) ('PASS', 7, 7, 10, 2, 'k=d on polynomial: False')
(7, 1) ('PASS', 7, 7, 10, 2, 'k=d on polynomial: False')
(6, 2) ('PASS', 7, 7, 10, 2, 'k=d on polynomial: False')
(5, 3) ('PASS', 7, 7, 10, 2, 'k=d on polynomial: False')
(4, 4) ('PASS', 7, 7, 10, 2, 'k=d on polynomial: False')
```

### summarize_tower.out
```
cells (n, depth) in analyze_tower.out: 169
cells with a genuine exact test (npts >= degree+2): 64
cells where minimal polynomial degree != floor(depth/2): []
  degree n=2: depths 0..0 tested (of 0..0); all have degree floor(depth/2)
  degree n=3: depths 0..2 tested (of 0..2); all have degree floor(depth/2)
  degree n=4: depths 0..4 tested (of 0..4); all have degree floor(depth/2)
  degree n=5: depths 0..6 tested (of 0..6); all have degree floor(depth/2)
  degree n=6: depths 0..8 tested (of 0..8); all have degree floor(depth/2)
  degree n=7: depths 0..10 tested (of 0..10); all have degree floor(depth/2)
  degree n=8: depths 0..9 tested (of 0..12); all have degree floor(depth/2)
  degree n=9: depths 0..7 tested (of 0..14); all have degree floor(depth/2)
  degree n=10: depths 0..5 tested (of 0..16); all have degree floor(depth/2)
  degree n=11: depths 0..3 tested (of 0..18); all have degree floor(depth/2)
```

### check_tower_to_kernel.out (Lemma 2)
```
k=2: kernels of degree 1..1 (keys checked: 1) equal tower->kernel linear image, no +-i simple poles: True  bad=[]
k=3: kernels of degree 1..2 (keys checked: 3) equal tower->kernel linear image, no +-i simple poles: True  bad=[]
k=4: kernels of degree 1..3 (keys checked: 5) equal tower->kernel linear image, no +-i simple poles: True  bad=[]
k=5: kernels of degree 1..4 (keys checked: 8) equal tower->kernel linear image, no +-i simple poles: True  bad=[]
k=6: kernels of degree 1..5 (keys checked: 11) equal tower->kernel linear image, no +-i simple poles: True  bad=[]
k=7: kernels of degree 1..6 (keys checked: 15) equal tower->kernel linear image, no +-i simple poles: True  bad=[]
k=8: kernels of degree 1..7 (keys checked: 19) equal tower->kernel linear image, no +-i simple poles: True  bad=[]
k=9: kernels of degree 1..8 (keys checked: 24) equal tower->kernel linear image, no +-i simple poles: True  bad=[]
k=10: kernels of degree 1..9 (keys checked: 29) equal tower->kernel linear image, no +-i simple poles: True  bad=[]
k=11: kernels of degree 1..10 (keys checked: 35) equal tower->kernel linear image, no +-i simple poles: True  bad=[]
k=12: kernels of degree 1..11 (keys checked: 41) equal tower->kernel linear image, no +-i simple poles: True  bad=[]
k=13: kernels of degree 1..12 (keys checked: 48) equal tower->kernel linear image, no +-i simple poles: True  bad=[]
k=14: kernels of degree 1..13 (keys checked: 55) equal tower->kernel linear image, no +-i simple poles: True  bad=[]
TOWER->KERNEL LEMMA verified for all levels above: True
```

### analyze_unfolded.out
```
unfolded levels: [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
key (0, 1) d=1: minimal poly degree in k (k>=d+1) = 0;  conjectured bound d-1 = 0;  points 11
key (1, 0) d=1: minimal poly degree in k (k>=d+1) = 0;  conjectured bound d-1 = 0;  points 11
key (0, 2) d=2: minimal poly degree in k (k>=d+1) = 1;  conjectured bound d-1 = 1;  points 10
key (1, 1) d=2: minimal poly degree in k (k>=d+1) = 1;  conjectured bound d-1 = 1;  points 10
key (2, 0) d=2: minimal poly degree in k (k>=d+1) = 1;  conjectured bound d-1 = 1;  points 10
key (0, 3) d=3: minimal poly degree in k (k>=d+1) = 2;  conjectured bound d-1 = 2;  points 9
key (1, 2) d=3: minimal poly degree in k (k>=d+1) = 2;  conjectured bound d-1 = 2;  points 9
key (2, 1) d=3: minimal poly degree in k (k>=d+1) = 2;  conjectured bound d-1 = 2;  points 9
key (3, 0) d=3: minimal poly degree in k (k>=d+1) = 2;  conjectured bound d-1 = 2;  points 9
key (0, 4) d=4: minimal poly degree in k (k>=d+1) = 3;  conjectured bound d-1 = 3;  points 8
key (1, 3) d=4: minimal poly degree in k (k>=d+1) = 3;  conjectured bound d-1 = 3;  points 8
key (2, 2) d=4: minimal poly degree in k (k>=d+1) = 3;  conjectured bound d-1 = 3;  points 8
key (3, 1) d=4: minimal poly degree in k (k>=d+1) = 3;  conjectured bound d-1 = 3;  points 8
key (4, 0) d=4: minimal poly degree in k (k>=d+1) = 3;  conjectured bound d-1 = 3;  points 8
key (0, 5) d=5: minimal poly degree in k (k>=d+1) = 4;  conjectured bound d-1 = 4;  points 7
key (1, 4) d=5: minimal poly degree in k (k>=d+1) = 4;  conjectured bound d-1 = 4;  points 7
key (2, 3) d=5: minimal poly degree in k (k>=d+1) = 4;  conjectured bound d-1 = 4;  points 7
key (3, 2) d=5: minimal poly degree in k (k>=d+1) = 4;  conjectured bound d-1 = 4;  points 7
key (4, 1) d=5: minimal poly degree in k (k>=d+1) = 4;  conjectured bound d-1 = 4;  points 7
key (5, 0) d=5: minimal poly degree in k (k>=d+1) = 4;  conjectured bound d-1 = 4;  points 7
```

### verify_dual.out / verify_single_log.out / verify_branch_form.out
```
mp.mp.dps = 30
p=3.0: direct 6.7449931401263401361 | Cauchy 6.7449931401263401361*(P/8 form: x8 included) | F-form 6.7449931401263401361
   |direct-Cauchy| = 7.89e-31  |direct-Fform| = 7.89e-31  |M - Y-form| = 4.93e-32
p=4.0: direct 7.0176979435640416471 | Cauchy 7.0176979435640416471*(P/8 form: x8 included) | F-form 7.0176979435640416471
   |direct-Cauchy| = 0.0  |direct-Fform| = 0.0  |M - Y-form| = 4.93e-32
p=2.5: direct 6.5467047229647988524 | Cauchy 6.5467047229647988524*(P/8 form: x8 included) | F-form 6.5467047229647988524
   |direct-Cauchy| = 7.89e-31  |direct-Fform| = 7.89e-31  |M - Y-form| = 0.0
nice -n 10 timeout 900 $PY -B verify_dual.py  0.13s user 0.02s system 25% cpu 0.558 total
closed form Psi = (1/2) e^{eps l1}/t [sqrt(1+y^2)-1-y] reproduces genintegrand coefficients for k=1..7, two exact t values: True
nice -n 10 timeout 1000 $PY -B verify_single_log.py  3.23s user 0.11s system 17% cpu 18.952 total
mp.mp.dps = 40
v= (0.1499999999999999944488848768742172978818 + 0.9699999999999999733546474089962430298328j)  |(a) diff| = 1.62e-41  |s-(sigma-i pi/2)| = 0.0  |(b) t diff| = 5.74e-42  |(b) l1_sigma diff| = 9.34e-28
v= (-0.1000000000000000055511151231257827021182 + 1.040000000000000035527136788005009293556j)  |(a) diff| = 1.15e-41  |s-(sigma-i pi/2)| = 0.0  |(b) t diff| = 5.74e-42  |(b) l1_sigma diff| = 3.95e-27
v= (0.0 + 0.9000000000000000222044604925031308084726j)  |(a) diff| = 0.0  |s-(sigma-i pi/2)| = 0.0  |(b) t diff| = 0.0  |(b) l1_sigma diff| = 1.75e-27
```
