# LEAD1 — the generating function of the leading coefficients c(w)

`c(w)` = coefficient of `k^{w-2}` in the weight-`w` part `W_w(k)` of the squircle coefficients `M_k`
(`P(p) = 8(1 + M_0/p + M_1/p^2 + ...)`, `M_k = sum_w W_w(k)`).  `lam = log 2`, `z_n = zeta(n)`.
Everything numeric below was printed by a script of this folder (`step_*.py`, driver `verify_lead1.py`).

## 1. The closed form

Let `A_n` be the Euler zigzag numbers (`sec x + tan x = sum A_n x^n/n!`: 1,1,1,2,5,16,61,272,1385,...),
`eps_m = +1,-1,-1,+1` of period 4 (`= (-1)^floor((m+1)/2)`), and

    E(y)  =  Gamma(1 + y/2) Gamma(1/2 - y/2) / sqrt(pi)
          =  2^y Gamma(1-y) Gamma(1+y/2) / Gamma(1-y/2)                       [duplication formula]
          =  exp( y log 2 + sum_{n>=2} s_n y^n/n ),
             s_n = zeta(n) (n even),   s_n = (1 - 2^{1-n}) zeta(n) = eta(n) (n odd),

with Taylor coefficients `e_w = [y^w] E(y)` (polynomials in lam, pi^2, odd zetas; e_0 = 1, e_1 = lam).  Then

    c(w)  =  eps_{w-2} * A_{w-1} / (2 (w-2)!)  *  e_w          (w >= 2)

Equivalently, since `eps_m A_{m+1}/m! = [z^m] (1 - sinh z)/cosh^2 z` (proof: `(1-sinh z)/cosh^2 z
= sech^2 z - sech z tanh z = d/dz(sech z + tanh z)`, and `[x^d](sech x + tanh x) = (-1)^floor(d/2) A_d/d!
= eps_{d-1} A_d/d!`, both checked symbolically for d <= 10 in step_54),

    F(u) = sum_{w>=2} c(w) u^w  =  (1/2) * [ u^2 (1 - sinh u)/cosh^2 u  \odot  E(u) ]          (Hadamard product),

and the "reduced" generating function is a pure Gamma product:

    sum_{w>=2}  [ 2 (w-2)! / (eps_{w-2} A_{w-1}) ]  c(w)  u^w  =  E(u) - 1 - u log 2.

Verified against all c(w), w = 2..9, at 406 digits (step_21) and rerun at dps 240 (verify_lead1.py, part D).

## 2. First instances (from the e_w, expanded; all confirmed by independent LLL, two precisions)

    c(2) = lam^2/4 + pi^2/24
    c(3) = -(lam^3/12 + pi^2 lam/24 + z3/8)
    c(4) = -(lam^4/48 + pi^2 lam^2/48 + lam z3/8 + pi^4/320)
    c(5) = (12 lam^5 + 20 pi^2 lam^3 + 180 lam^2 z3 + 9 pi^4 lam + 30 pi^2 z3 + 270 z5)/3456
    c(6) = (168 lam^6 + 420 pi^2 lam^4 + 5040 lam^3 z3 + 378 pi^4 lam^2 + 2520 pi^2 lam z3 + 22680 lam z5 + 61 pi^6 + 3780 z3^2)/362880
    c(7), c(8), c(9): see leading_coeffs.json / lll_relations_dps400.json (12, 14, 19 monomials).

Note `c(w) = A_{w-1} * (rational polynomial in lam, pi^2, zeta odd)/2(w-2)!` except for the top pi-power,
which also carries A-factors of lower index (e.g. `61 pi^6` in c(6) with 61 = A_6, `21437 = 17*1261 pi^8` in c(8)).

## 3. Derivation

### 3.1 The top kernel coefficients (input; EXACT for d <= 8)

`c(w) = (1/(w-3)!) * sum_{keys (i,j), i+j=w-1, i>=j} sum_{poles r} D_{w-1}^{(i,j)}(r) N_{(i,j),r}`,
`N_{(i,j),r} = int_b^a dv/(v-r) log^i t log^j(1-t)`, `t = (v-b)(v+a)/(4v)`, `a = 1+sqrt2`, `b = sqrt2-1`,
poles `r in {-a, 0, -b, a}` (pole `+b` never occurs).  The top Newton coefficients `D_d` were extracted
from the exact kernels `kern_k.pkl` (k = d+1..18) of WTS/theory (step_08, Dtop.json).  They obey
(exact check for all keys and poles, d = 1..8, step_36):

    delta_d := [x^d](sech x + tanh x)  (= eps_{d-1} A_d/d!)
    pole -a (all keys):                 D_d = (-1)^{j+1} binom(d,j) delta_d
    pole  a (odd d only):               D_d = (-1)^{j}   binom(d,j) delta_d * (1 if j<(d-1)/2 else 1/2)
    pole -b (even d only, j < d/2):     D_d = (-1)^{j}   binom(d,j) delta_d
    pole  0 (diagonals only):           D_d = (-1)^{j}   binom(d,j) delta_d / 2
       (even d, j=d/2; odd d, j=(d-1)/2, which also carries pole -b with the -(-1)^j.../2 value)
    all other (key,pole):               D_d = 0

This table is exact data for d <= 9 (d = 9: determined fit on k = 10..18, zero spare points, step_57);
for d >= 10 it is the content of the (unproved) tower law T of WTS/theory (the pole -a column alone is
check_top.py, d = 2..8).  All sqrt2-parts vanish.

### 3.2 The substitution  u = (v - 1/v)/2   (t = (1+u)/2)

With `u = (v-1/v)/2` one has EXACTLY `t = (1+u)/2`, `1-t = (1-u)/2`, `[b,a] -> [-1,1]`,
`dv = v/sqrt(1+u^2) du`, and the pole-combination measures become (exact algebra, quotient ring
mod s^2 = 1+u^2, step_38; `S = sqrt(1+u^2)`):

    odd  (j<(d-1)/2):  1/(v+a)-1/(v-a)            ->  [1 + sqrt2 u/S]/(1-u^2) du
    even (j<d/2):      1/(v+a)-1/(v+b)            ->  [1 - sqrt2/S]/(1-u^2)   du
    diag (even):       1/(v+a)-1/(2v)             ->  [1 - sqrt2/S]/(2(1+u))  du = (1-u)/2 * (even measure)
    oddlast:           1/(v+a)-1/(2(v-a))+1/(2(v+b))-1/(2v)
                                                    ->  1/(2(1+u)) + sqrt2 u/(S(1-u^2))  du = (odd measure) - 1/(2(1-u))

Hence, writing `L1 = log((1+u)/2)`, `L2 = log((1-u)/2)`, `M(i,j; m) = int_{-1}^1 L1^i L2^j m(u) du`
(all convergent for the occurring (i,j)):

    c(d+1) (d-1)!/delta_d  =  sum_{j=0}^{floor((d-1)/2)} (-1)^{j+1} binom(d,j) M(d-j, j; m_{d,j}).   (*)

### 3.3 Even d  (odd weights w = d+1): complete proof

The even measure `m_e(u) = [1-sqrt2/S]/(1-u^2)` is even in u and extends to `-1/4` at u = +-1.
For the middle key, `m_diag = (1-u)/2 m_e`.  Pairing term j with term d-j under u -> -u gives
`sum_{j<d/2} T_j = (F_d - T_{d/2})/2` with `F_d = sum_{j=0}^d T_j = -int log^d((1+u)/(1-u)) m_e du`
(the full binomial sum is exactly `-(L1-L2)^d`), and the middle-key correction is
`(1/2) T_{d/2}(u m_e) = 0` (odd integrand).  Therefore, convergent integrals throughout,

    c(d+1)(d-1)!/delta_d = -(1/2) int_{-1}^1 log^d((1+u)/(1-u)) m_e(u) du
                          = (1/2) int_0^1 (-log w)^d [(1+w)/sqrt(1+w^2) - 1] dw/w      (w = e^{-s}, u = tanh(s/2)).

This integral is evaluated in closed form for ALL even d by its generating function: with
`F(y) = int_0^1 w^{-y}[(1+w)/sqrt(1+w^2)-1] dw/w = sum_{n>=1} (-1)^n binom(2n,n)/4^n [1/(2n-y) + 1/(2n+1-y)]`
(termwise integration, valid |y|<1), one has

    (F(y) + F(-y))/2  =  (E(y) - E(-y))/(2y)   (|y|<1),

PROOF: `E(y) = (1/2) B(1+y/2, 1/2-y/2)`; split the Beta integral at `t = 1/2` and substitute
`t = w^2/(1+w^2)` (resp. `1-t`) to get `(E(y)-E(-y))/(2y) = (1/(2y)) int_0^1 (w^y - w^{-y})(w-1)/(1+w^2)^{3/2} dw`,
then integrate by parts with `W = 1 - (1+w)/sqrt(1+w^2)` (boundary terms vanish for |y|<1).  Both steps
are verified numerically in verify_lead1.py (part E).  Comparing Taylor coefficients at y = 0,

    (1/2) int_0^1 (-log w)^d [(1+w)/sqrt(1+w^2) - 1] dw/w  =  (d!/2) e_{d+1}   (even d; checked d <= 16, 60-80 digits).

Hence for even d, GIVEN the D-table at d:  `c(d+1) = delta_d/(d-1)! * (d!/2) e_{d+1} = eps_{d-1} A_d/(2(d-1)!) e_{d+1}`.  QED.

### 3.4 Odd d (even weights w = d+1): reduction + numerical identification

Write `m_o(u) = [1+sqrt2 u/S]/(1-u^2) = 1/(1-u) - u m_e(u)` (exact).  Then (*) splits:

* `u m_e` part: antisymmetric under x<->y, so the half-binomial sum collapses to the full sum
  `-(1/2) int_{-1}^1 log^d((1+u)/(1-u)) u m_e(u) du = (1/2) int_0^1 (-log w)^d (1-w)/(1+w)[(1+w)/sqrt(1+w^2)-1] dw/w`
  (PROVED, convergent).
* `1/(1-u)` part plus the oddlast boundary `-1/(2(1-u))`: with `A_{p,q} = int L1^p L2^q/(1-u) du
  = (-1)^p int_0^inf r^p log^q(1-e^{-r})/(e^r-1) dr` (Bose-Einstein/log integrals, MZV world;
  A_{p,0} = (-1)^p p! zeta(p+1)), the bracket
  `sum_{j<=mid} (-1)^{j+1} binom(d,j) A_{d-j,j} - (1/2)(-1)^{mid+1} binom(d,mid) A_{mid+1,mid}`, mid=(d-1)/2.

The total equals (d!/2) e_{d+1} for every odd d tested (d = 1..11, 60-80 digits, step_47/50).  The
single-integral target

    c(d+1)(d-1)!/delta_d  =  (1/(2(d+1))) int_0^1 (log w)^{d+1} (1+w)/(1+w^2)^{3/2} dw   (odd d)

is PROVED to equal `(d!/2) e_{d+1}` for all odd d: from `[(E(y)+E(-y))/2 - 1]/(2y)
= (1/(4y)) int_0^1 (w^y + w^{-y} - 2)(1+w)/(1+w^2)^{3/2} dw` (same Beta substitution + value at y = 0)
on extracting `[y^d]`.  The remaining algebraic step — proving the A-bracket + piece-2 sum equals this
integral for all odd d — was NOT completed (the naive parity trick fails for odd d: the B(x,y)/2 and
sqrt2 J parts of `G_odd = B/2 + sqrt2 2^{-(x+y)} J` are separately divergent at the endpoints and only
their sum is finite; see Section 4 notes).

### 3.5 Assembly

For both parities the key-sum evaluates to `(d!/2) e_{d+1}`, so with d = w-1

    c(w) = delta_d/(d-1)! * (d!/2) e_{d+1} = delta_d * (d/2) * e_{d+1},

and since `delta_d = eps_{d-1} A_d/d!`, the prefactor is `delta_d d/2 = eps_{d-1} A_d/(2 (d-1)!)`.
With `d = w-1` this is exactly `c(w) = eps_{w-2} A_{w-1}/(2 (w-2)!) * e_w`, the stated law.

## 4. Status: what is proved vs checked

| statement | status |
|---|---|
| D-table (Sec 3.1) for d <= 9 | EXACT (rational arithmetic on kern_1..18; steps 08, 36, 57) |
| D-table for all d | conjecture (= tower law T of WTS/theory, unproved) |
| reduction of c(w) to combo integrals (*), given the D-table | PROVED (exact algebra, step_38) |
| even-d collapse + evaluation => law for odd w | PROVED given the D-table at d = w-1 (Sec 3.3) |
| odd-d key-sum == single integral (=> law for even w) | CHECKED 37-80 digits for d = 1..9 with the exact kernel data (steps 33, 50, 61) and d = 11 on the integral side; not proved for all d |
| single integrals == (d!/2) e_{d+1}, both parities | PROVED via generating functions (Sec 3.3, 3.4); checked to d = 16 |
| c(w) values w = 2..9 (405 digits), w = 10 (250 digits) | controls c(2),c(3),c(4),c(6) to ~410 digits; two independent runs (420/250 dps) agree; law vs c(10) to 237 digits |
| identifications of c(5)..c(9) | LLL relations IDENTICAL at dps 420 (D=390) and dps 250 (D=225), residuals 1e-217 or better |
| the law vs all c(w), w <= 9 | 406 digits (step_21); rerun 239 digits (verify_lead1.py) |

In particular c(3), c(5), c(7), c(9) are covered by the PROVED even-d branch (given the D-table),
and c(2), c(4) were already proved independently (LEDGER/WT4).
