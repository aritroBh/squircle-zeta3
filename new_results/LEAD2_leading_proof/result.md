# LEAD2 — the leading-coefficient law is a theorem (given the top-coefficient table)

## VERDICT

**PROVED (conditional on the D-table).** Given the top-coefficient table `D_d` of the kernel
generating functions (binomial x zigzag table; exact rational data for `d <= 9`, and for `d >= 10`
the tower law `T` of WTS/theory, which remains the one open input — attacked separately by TOW2),
the law holds for **all** `w >= 2`:

    c(w)  =  eps_{w-2} * A_{w-1} / (2 (w-2)!) * [y^w] E(y),
    E(y)  =  Gamma(1+y/2) Gamma(1/2-y/2) / sqrt(pi),

`eps_m = (-1)^floor((m+1)/2)`, `A_n` Euler zigzag numbers.  The **even-weight case (d = w-1 odd) is
proved here for the first time**; the odd-weight case (d even) is the LEAD1 proof, refereed and
re-derived: verdict **OK** (one documentation typo, no mathematical gap — see below).  Generating
function:

    sum_{w>=2} c(w) u^w  =  (1/2) [ u^2 (1-sinh u)/cosh^2 u  ⊙  E(u) ]        (Hadamard product)

(equivalent to the law coefficientwise; the reduced GF is `E(u) - 1 - u log 2`, a pure Gamma product).

## The new ingredient (odd d / even weights)

The LEAD1 obstruction was a "finite part": the key-sum generating function
`G_odd(x,y) = B(x,y)/2 + sqrt2 2^{-(x+y)} J(x,y)` has separately divergent pieces. The cure is to
keep every integral **convergent on the strip `-1 < Re s < 1`** and to split the key sum via the exact
identity `m_odd = 1/(1-u) - u m_e` (Lemma 2.1 of proof_lead.md):

1. **`u m_e` part**: the measure `u m_e` is odd, so the half-binomial sum collapses exactly
   (Lemma 5.1); its generating function is `GG(z)/4` with
   `GG(s) = int_0^inf w^{-s} (1-w)/(1+w) h(w) dw`, evaluated in closed form on the strip:
   `GG(s) = 2 pi/sin(pi s) - (E(s)+E(-s))/s` (Master Lemma (iii), via Mellin pieces
   `(1/2)B(s/2,(1-s)/2)` and `pi/sin(pi s)`).
2. **`1/(1-u)` part + boundary** (`T_d`): with `A_{p,q} = (-1)^p int r^p log^q(1-e^{-r})/(e^r-1) dr`,
   the two-parameter GF is exactly `H(x,y) = B(1+x,y)` (Lemma 5.2). The functional equation
   `y B(1+x,y) = x B(1+y,x)` gives the coefficient recursion `aa_{p,q} = aa_{q+1,p-1}`, whose
   involution `(p,q) -> (q+1,p-1)` folds the half-sum into a full diagonal (Theorem 5.4):
       T_d/d! = (1/2)[z^d] B(1+z,-z) = (1/2)[z^d](pi/sin pi z) = (1-2^{-d}) zeta(d+1).
3. **Assembly**: the `pi/sin(pi z)` finite parts (the `(1-2^{-d}) zeta(d+1)` terms) **cancel exactly**
   between `T_d` and the `u m_e` part, leaving
   `K_d/d! = (1/4)[z^{d+1}] WW(z)` with `WW(s) = int_0^inf w^{-s}(1+w)/(1+w^2)^{3/2} dw = E(s)+E(-s)`
   (Master Lemma (i)), i.e. `K_d = (d!/2) e_{d+1}` for all odd `d` (Theorem 5.5).

Both sides of the cancellation are the *same* `(1-2^{-d}) zeta(d+1) d!` obstruction seen from the two
pieces — this explains the LEAD1 numerics `(piece1 - boundary) - piece2 = (d!/2) e_{d+1}`.

## Master Lemma (the three strip generating functions; proof_lead.md Section 3)

On `-1 < Re s < 1`, with `h = [(1+w)/sqrt(1+w^2)-1]/w`, `q = (1+w)/(1+w^2)^{3/2}`, `g1 = (1-w)h/(1+w)`:

    int_0^inf w^{-s} q  dw  =  E(s) + E(-s)
    int_0^inf w^{-s} h  dw  =  (E(s) - E(-s))/s
    int_0^inf w^{-s} g1 dw  =  2 pi/sin(pi s) - (E(s)+E(-s))/s

All three are convergent integrals (no finite parts): the even-d law uses the second, the odd-d law
the first and third. Proofs: `E = (1/2)B(1+s/2,1/2-s/2)`, split the Beta integral at `t = 1/2`,
substitute `t = w^2/(1+w^2)`, integrate by parts (boundary vanishes exactly on the strip), and use
the inversion symmetries `h(1/w) = w^2 h(w)`, `q(1/w) = w^2 q(w)`.

## Task 1 — referee report on the even-d (odd-weight) proof: **OK**

Re-derived every step (proof_lead.md Sections 2-4): the u-substitution and pole-combo measures
(symbolic, exact); the parity collapse `K_d = -(1/2) int (L1-L2)^d m_e` (pairing `j <-> d-j`,
middle-key correction killed by oddness); the reduction `int m_e -> -(1/2) int h`; the
partial-fraction series of `F(s)` with its convergence (`Re s < 1` for the integral, absolute for the
series off poles; uniform convergence of the alternating binomial series justifies termwise
integration); the identity `(F+F-)/2 = (E-E-)/(2s)` via Beta + IBP with boundary terms vanishing
exactly for `-1 < Re s < 1`; coefficient extraction at `s = 0` (inside the strip). **Erratum (docs
only):** LEAD1's genfun.md prints `F(y) = sum_{n>=1} (-1)^n binom(2n,n)/4^n [1/(2n-y)+1/(2n+1-y)]`;
the second sum must start at `n = 0` (the term `1/(1-y)` was dropped in the writeup; LEAD1's numerics
used quadrature and are unaffected — verified: residual of the printed form is exactly `-1/(1-y)`).
No other issues: **OK**.

## Evidence table (all numbers printed by code; verify_lead2.py reruns everything)

| check | range | residual |
|---|---|---|
| D-table == binomial x zigzag, exact rationals | d = 1..9, all keys/poles | 0 (exact); polynomial law in k exact (d=9: determined fit) |
| u-measures + t=(1+u)/2 + m_odd = 1/(1-u) - u m_e | symbolic | sympy radsimp == 0 |
| WW = E+E-, FF = (E-E-)/s, GG = 2pi/sin - (E+E-)/s | s in {0.3,-0.4,0.73,0.2i,...} | <= 8.7e-46 (interior), 1.4e-17 (s=0.73) |
| Beta substitution (4); IBP E-E- = s(F+F-) | same points | <= 8.7e-46 |
| F-series (n>=0) == quadrature | y = 0.3,-0.45 | <= 7.8e-62 |
| recursion aa_{p,q} = aa_{q+1,p-1} | p=1..5, q=0..4 | <= 4.4e-58 |
| odd-d chain key sum = T_d - piece2 = (d!/2) e_{d+1} | d = 1,3,5,7,9,11 (70 dps) | <= 3.9e-51 |
| closed forms T_d = d!(1-2^{-d}) zeta(d+1), piece2 = T_d - (d!/2)e_{d+1} | d = 1,3,5,7,9,11 | <= 6.1e-64 / 3.7e-51 (70 dps); d=5,7: 2.5e-119 / 5.5e-114 (120 dps) |
| A_{p,q}: M-form == Bose-Einstein form (eq. 8) | (p,q) in {(1,0),(2,1),(3,2),(5,3),(4,4)} (120 dps) | 0.0 |
| even-d collapse + evaluation | d = 2,4,6,8,10 (70 dps) | <= 3.8e-53 |
| full pipeline kernels -> c(w) vs law | w = 2..10 (dps 100) | <= 8.6e-91 |

## What is proved vs assumed

| statement | status |
|---|---|
| D-table for d <= 8 | EXACT (rational data, spare Newton points; `Dtop_lead2.pkl`) |
| D-table for d = 9 | EXACT as determined fit (zero spare points; matches table everywhere) |
| D-table for d >= 10 | HYPOTHESIS (= tower law T of WTS/theory; TOW2) |
| reduction (1) -> key sum (3), u-measures | PROVED (exact algebra) |
| Master Lemma (WW, FF, GG closed forms) | PROVED (Beta + IBP + Mellin; convergence on the strip shown) |
| even-d collapse + evaluation | PROVED (all even d) |
| odd-d decomposition, involution, cancellation | PROVED (all odd d) — **new** |
| the law c(w) for all w, given the D-table | **PROVED** |
| GF form Hadamard `u^2(1-sinh u)/cosh^2 u ⊙ E(u)` | PROVED (zigzag sign lemma 6.1) |

## Deliverables (this folder)

- `proof_lead.md` — the complete proof (setup, Master Lemma, even-d refereeing, odd-d proof, assembly)
- `verify_lead2.py` — reruns every identity A-G, prints PASS/FAIL (< 15 min)
- `lead2_lib.py` — shared library; `cw_dps100.json`, `Dtop_lead2.pkl` — data
- `progress.log` — step log; `step_*.py` — producing steps
