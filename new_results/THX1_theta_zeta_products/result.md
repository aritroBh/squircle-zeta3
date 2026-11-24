# THX1 — the diamond-end Fourier coefficients are polynomials in single zeta values: Q5..Q9 PROVED exactly

**Verdict.** For every `n = 5,6,7,8,9` and every `m = 0..⌊n/2⌋−1`, the identity `D̃_{n,m} = (zeta-product polynomial)` stated in the paper (and hence `Q5..Q9` exactly) is **PROVED**, in two independent ways:

1. **Direct certificates.** Each `V = D̃_{n,m} − (zeta products)` is written as an explicit exact rational combination of **proved relation instances** among level-2 MZVs (finite double shuffle from absolute convergence, and the duality `t ↦ (1−t)/(1+t)`; *no regularization, no conjecture anywhere*), and the combination is verified by exact rational arithmetic. 16 certificates, all stored.
2. **Assembly from certified atoms.** Every atom `I_{m,j} = ∫₀^∞ s^m log^j(1+e^{−s}) ds` with `m+j+1 ≤ 9` is reduced to an explicit alternating-MZV basis, each with its own exact certificate; assembling `D̃_{n,m} = 2 Σ_j e_{m,n−1−j} I_{n−1−j,j}/j!` with the rational e-table, **every non-zeta basis element (log 2, Li_r(1/2), and the new alternating primitives) cancels identically** and the zeta coefficients match the paper's `Q_n(u)` exactly.

**Status of the general statement (all n).** OPEN, but reformulated precisely (§5): the master generating function is proved here, and the exact data through `n = 9` show that the zeta-landing combinations of the atoms are *exactly* the `D̃_{n,m}` combinations (plus the trivial singleton `(n−2)!η(n)` at odd `n`). The general theorem is equivalent to the dimension conjecture `dim V_n = ⌈n/2⌉` (§5.3), verified exactly for `n ≤ 9`.

---

## 1. Setup and the proved relation engine

Words over `{0,1,−1}` (letter `a` = form `dt/(t−a)`), value `G(w;1)`; `ζ(S) = (−1)^d G(w(S);1)`. Relation families used (each a *proved theorem* applied with rational arithmetic):

- **DS** — finite double shuffle `u ш v = st(u)⋆st(v)` for convergent words (absolute convergence).
- **DUAL** — the involution `t ↦ (1−t)/(1+t)` on convergent words (change of variables).

Codimension of the relation space at weight `w` (convergent words `4·3^{w−2}`), computed with FLINT over `GF(p)` and cross-checked exactly at `w ≤ 6`:

| w | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|
| codim | 1 | 2 | 3 | 5 | 8 | **13** | **21** | **34** | **55** |

= the Deligne–Goncharov dimension (Fibonacci) at every weight ≤ 9, so **DS+DUAL suffice through weight 9**. (Caution: the regularized-DS relations as implemented in RUL1 are *false* — all 162 weight-6 instances are numerically nonzero; the naive `T=0` comparison is not a theorem. They were excluded; DS+DUAL suffice.)

## 2. The atoms `I_{m,j}` as alternating MZVs

`log^j(1+e^{−s}) = j! Σ_{l≥j} s(l,j) e^{−ls}/l!` (signed Stirling) gives, exactly,

> **`I_{m,j} = (−1)^j j! m! ζ(\overline{m+2}, 1^{j−1}) = j! m! G(0^{m+1}(−1)^j; 1)`**,

an alternating MZV of weight `m+j+1` and depth `j` (PROVED identity; verified numerically vs quadrature: worst residual **2.06e-90** over 16 pairs at 100 dps).

Every `I_{m,j}` with `m+j+1 ≤ 9` is reduced EXACTLY to an explicit basis with a certificate (w ≤ 5: already in RUL1; w = 6..9: this work). Basis primitives beyond `log2, ζ2, ζ3, ζ5, ζ7, ζ9, Li_4..Li_9(1/2)`:

- `X6 = ζ(\bar5,\bar1)` (weight 6);
- `X7a = ζ(\bar5,1,1)`, `X7b = ζ(\bar4,2,1)` (weight 7);
- `ζ(\bar6,1,1)`, `ζ(\bar5,2,1)`, `ζ(\bar5,1,1,1)`, `ζ(\bar4,2,1,1)` (weight 8);
- `ζ(\bar7,1,1)`, `ζ(\bar6,2,1)`, `ζ(\bar5,3,1)`, `ζ(\bar5,1,1,1,1)`, `ζ(\bar4,2,1,1,1)`, and `X9fp = G((−1,−1,−1,0,0,0,−1,0,−1);1)` (weight 9; the 55th direction is the class of `−ζ(2,1,1,1,1,1,\bar2)`, represented by an E1-computable word).

Reductions (certified; numerically checked to ≥ 107 digits — full data in `I_reductions_*.pkl`):

```
I_{m,1} = m! η(m+2)                                  (pure zeta, all m)
I_{0,2} = ζ3/4,  I_{0,3} = -6Li_4 - (21/4)ζ3·log2 + (3/2)ζ2·log2² - log2⁴/4 + (12/5)ζ2², ...
I_{2,2} = -(29/8)ζ5 + 2ζ2ζ3
I_{3,2} = -(234/35)ζ2³ + 9ζ3² + (93/4)ζ5·log2 - 12 X6
I_{2,3} = (321/35)ζ2³ - (21/2)ζ3² - (279/8)ζ5·log2 + 18 X6
I_{4,2} = (84/5)ζ3ζ2² + 24ζ5ζ2 - (753/8)ζ7
I_{2,4} = -24ζ3ζ2² - 24ζ5ζ2 + (945/8)ζ7 + 72 X7a
I_{3,3} = -36 X7a   I_{3,4} = 144 ζ(\bar5,1,1,1)   I_{4,3} = -144 ζ(\bar6,1,1)
I_{5,3} = -720 ζ(\bar7,1,1)   I_{3,5} = -720 ζ(\bar5,1,1,1,1)
I_{6,2} = (2232/7)ζ2³ζ3 + 504ζ2²ζ5 + 720ζ2ζ7 - (68805/16)ζ9
I_{1,6} = ... - 288 ζ(\bar6,1,1) + 360 ζ(\bar5,1,1,1) + ...  (16 terms)
I_{2,6} = 1008ζ2³ζ3 + 1728ζ2²ζ5 + 2160ζ2ζ7 - 240ζ3³ - (214095/16)ζ9 - 3600ζ(\bar7,1,1) + 2160ζ(\bar5,1,1,1,1)
I_{0,5}..I_{0,8}: the full log2/Li towers (printed in the run logs, coefficients in I_reductions_*.pkl)
```

## 3. Q5..Q9: PROVED (16 certificates)

`D̃_{n,m} = 2 Σ_{j=1}^{n−1−2m} e_{m,n−1−j} I_{n−1−j,j}/j!` with the exact rational e-table (`e_{r,m} = 2^{3/2}[x^m](gd x)^{2r} H(x)`, `H(x) = 2^{−3/2}e^{x/2}sech^{3/2}x`; K3, PROVED). Each `V = D̃_{n,m} − (stated zeta products)` is proved to be a combination of DS/DUAL relation instances:

| target | certificate | size | method |
|---|---|---|---|
| D̃_{5,0}, D̃_{5,1} | EXACT-VERIFIED | 2 rels | CRT(3 primes)+reconstruction |
| D̃_{6,0..2} | EXACT-VERIFIED | ≤6 rels | CRT(3 primes) |
| D̃_{7,0..2} | EXACT-VERIFIED | 869–871 rels, height ≤ 289 bits | exact sparse solve (sympy DomainMatrix) |
| D̃_{8,0..3} | EXACT-VERIFIED | 2882-rel basis | Dixon p-adic lifting (pk ≈ 2^4340) |
| D̃_{9,0..3} | EXACT-VERIFIED | 8693-rel basis | Dixon lifting to 2^16120 |

and then `Q_n(u) = Σ_m D̃_{n,m} (−4u)^m/(2m)!` matches `Qpolys.pkl` symbolically for all `n ≤ 9`:

```
Q5 = -(3/4 + 4u)ζ5 - (7/4 + 4u)ζ2ζ3
Q6 = (4039/256 + 733u/8 + 31u²)ζ6 - (7/8 + 2u)ζ3²
Q7 = (129/16 + 46u + 16u²)ζ7 + (141/16 + 50u + 16u²)ζ2ζ5 + (1071/160 + 181u/5 + 56u²/5)ζ2²ζ3
Q8 = -(2891923/12288 + 495965u/256 + 60889u²/48 + 127u³)ζ8 + (141/16 + 50u + 16u²)ζ3ζ5 + (169/32 + 27u + 8u²)ζ2ζ3²
Q9 = -(22805/192 + 3907u/4 + 1916u²/3 + 64u³)ζ9 - (8343/64 + 4163u/4 + 660u² + 64u³)ζ2ζ7
     - (60093/640 + 29741u/40 + 2334u²/5 + 224u³/5)ζ2²ζ5 - (263713/4480 + 130973u/280 + 2062u²/7 + 992u³/35)ζ2³ζ3
     + (169/96 + 9u + 8u²/3)ζ3³
```

## 4. The cancellation, displayed (second proof route)

Assembling the certified atom reductions into `D̃_{7,0}` (e-weights `721/1920, 121/384, −17/144, −5/96, 1/120, 1/360`):

| element | I₅₁ | I₄₂ | I₃₃ | I₂₄ | I₁₅ | I₀₆ | sum |
|---|---|---|---|---|---|---|---|
| log2⁷ | 0 | 0 | 0 | 0 | 1/126 | −1/126 | **0** |
| ζ2·log2⁵ | 0 | 0 | 0 | 0 | −1/30 | 1/30 | **0** |
| ζ3·log2⁴ | 0 | 0 | 0 | 0 | 7/96 | −7/96 | **0** |
| Li₄·log2³ | 0 | 0 | 0 | 0 | 1/3 | −1/3 | **0** |
| Li₅·log2² | 0 | 0 | 0 | 0 | 1 | −1 | **0** |
| Li₆·log2 | 0 | 0 | 0 | 0 | 2 | −2 | **0** |
| Li₇ | 0 | 0 | 0 | 0 | 2 | −2 | **0** |
| X7a | 0 | 0 | 17/4 | −15/4 | −1/2 | 0 | **0** |
| ζ7 | 45423/1024 | −30371/1024 | 0 | −1575/256 | −159/64 | 2 | **129/16** |
| ζ2ζ5 | 0 | 121/16 | 0 | 5/4 | 0 | 0 | **141/16** |
| ζ2²ζ3 | 0 | 847/160 | 0 | 5/4 | 3/20 | 0 | **1071/160** |

Identical cancellation occurs for every `(n,m)`, `n = 5..9` (checked exactly).

## 5. The mechanism (task 3)

### 5.1 Master generating function (PROVED here)

Summing K3's coefficient formula over `n` (binomial resummation, verified numerically — differences equal the first omitted term exactly: 6.1e-7 at τ = 0.1, 6.2e-10 at τ = 0.05):

> **`D̃(τ;m) := Σ_{n≥2m+2} D̃_{n,m} τⁿ = 2^{5/2} ∫₀^∞ (gd x)^{2m} H(x) [(1+e^{−x/τ})^τ − 1] dx`**, `H(x) = e^{−x}(1+e^{−2x})^{−3/2}`.

### 5.2 A clean derivation of the coefficient formula (PROVED)

Insert `(1+e^{−s})^τ = Σ_{l≥1} C(τ,l) e^{−ls}` and the Laplace transform `∫₀^∞ e^{−ls}(gd τs)^{2m}H(τs) τ ds = 2^{−3/2} Σ_r e_{m,r} r! (τ/l)^{r+1}`:

> `D̃(τ;m) = 2 Σ_{r≥2m} e_{m,r} r! τ^{r+1} S_r(τ)`, `S_r(τ) = Σ_{l≥1} C(τ,l)/l^{r+1} = Σ_{k≥1} τ^k I_{r,k}/(k! r!)`,

i.e. `D̃_{n,m} = 2 Σ_{r+k+1=n} e_{m,r} I_{r,k}/k!` exactly. Also `F_r(t) := Σ_k t^k I_{r,k}/k! = ∫₀¹ log^r(1/u) ((1+u)^t − 1) du/u = r! Σ_{a≥1} C(t,a)/a^{r+1}` — the atoms are the Taylor coefficients of one parameterized beta-type integral per `r`.

### 5.3 The saturation statement (PROVED for n ≤ 9; conjectured for all n)

Let `V_n` = {rational combinations of the atoms `I_{n−1−k,k}` (`k = 1..n−1`) whose non-zeta part vanishes}. From the certified reductions:

| n | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|
| dim V_n | 3 | 3 | 4 | 4 | 5 |
| #e_m (rank) | 2 | 3 | 3 | 4 | 4 |

- **Even n: V_n = span(e_m) exactly (saturation).**
- **Odd n: V_n = span(e_m) ⊕ span{I_{n−2,1}}, and I_{n−2,1} = (n−2)! η(n) is trivially zeta.**

So through weight 9 the `D̃`-combinations are *all* the zeta-landing combinations (nothing hides), and the e-vectors are independent for all `n` (`e_{m,r} = 0` for `r < 2m`, `e_{m,2m} = 1` — triangular). **Conjecture: `dim V_n = ⌈n/2⌉` for all `n`; equivalently the general Q_n reduction.** The numerical D̃ tables through `n = 12` (K3) support it.

### 5.4 What is proved for all n (status)

- Top row `D̃_{2n,n−1} = 2(2n−2)! η(2n)` and the odd top (K3, PROVED) — zeta products.
- `m = 0` row for `n ≤ 10` (D11's exact perimeter coefficients, PROVED) — zeta products.
- All `(n,m)` with `n ≤ 9` (this work, PROVED).
- Beyond: OPEN (§5.3 conjecture; the cancellation engine is now fully explicit).

## 6. Evidence table (all numbers printed by the listed steps)

| item | statement | evidence | digits |
|---|---|---|---|
| I_{m,j} = j!m!G(0^{m+1}(−1)^j) | PROVED identity | max residual 2.06e-90 over 16 pairs vs quadrature | 100 |
| codim = DG at w ≤ 9 | 1,2,3,5,8,13,21,34,55 | FLINT rank GF(p); exact RREF at w ≤ 6 (13✓,21✓) | exact |
| D̃ formula = zeta targets | numeric anchor | ≤ 1.6e-98 (n ≤ 9) vs K3 PSLQ; D11 ≤ 1.8e-99 | 100 |
| master GF (5.1) | PROVED | τ=0.1: 6.1e-7, τ=0.05: 6.2e-10 (= first omitted term) | 50 |
| individual reductions | certificates + numerics | residuals ≤ 3.2e-149 (w≤7), ≤ 1.3e-143 (w8), ≤ 3.7e-107 (w9) | 107–150 |
| D̃ certificates | 16 exact certificates | verify_thx1.py: **43 PASS, 0 FAIL in 16 s** | exact |
| I-reduction certificates | 26 exact certificates (w6,7,8,9) | same verifier + run logs | exact |
| assembly cancellation | all (n,m), n ≤ 9 | every non-zeta coefficient exactly 0 (step logs) | exact |
| V_n saturation | dim = 3,3,4,4,5 | exact nullspaces from certified reductions | exact |

## 7. Honest caveats

- The w9 D̃ certificates were verified exactly in the main runs; `verify_thx1.py` re-verifies D̃_{9,0} exactly and D̃_{9,1..3} modulo a 62-bit prime (to stay under the 3-minute budget). All other certificates (n ≤ 8) are re-verified *exactly* by the verifier.
- Two PSLQ candidates (I_{0,7}, I_{1,6}) produced by PSLQ at 170 digits were **wrong**; they were corrected by the modular route (free-coordinate solve over GF(p) + rational reconstruction) and then certified exactly. Moral: 34+ constants needs ~800+ PSLQ digits; the modular route is the safe one.
- RUL1's RDS relations are false as implemented (all 162 at w6); unused here.

## 8. Files

- Engine: `rul1.py, rul1la.py, rul1_named.py, rul1_reduce.py, rul1num.py` (from RUL1), `thx1_elim.py`, `thx1_cert.py` (certificates), `thx1_dixon.py` (p-adic lifting).
- Data: `rels{5,6,7,8,9}_dsdual.pkl`, `cols{6,7,8,9}.pkl`, `struct_w{5,6,7,8,9}.pkl` (pivot structure), `w9_Ainv.npy` (mod-p inverse), `e_rm_ext.pkl`, `Dtilde_relations.pkl`, `Qpolys.pkl`.
- Certificates: `cert_D{5,6,7,8,9}_{m}.pkl` (16), `cert_I_{m}_{j}.pkl` (26).
- Reductions: `I_reductions_w2_7.pkl`, `I_reductions_w89.pkl`, `I_reductions_w9.pkl`, `extras_val.pkl`, `basis_vals_89.pkl`.
- Verifier: **`verify_thx1.py` (43 PASS, 0 FAIL, 16 s)**; logs: `progress.log`, `summary.md`; `step_*.py` = the full run (~145 steps: conventions/master GF (steps 14–28), relation engine + codims (steps 29–46), Q5–Q7 certificates (47–62), Dixon lifting + Q8/Q9 (63–86), bases + individual reductions (89–131), w9 certificates + saturation + assembly (132–146)).