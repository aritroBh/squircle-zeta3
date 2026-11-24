import os, time
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
part1 = r"""# LEAD2 — proof of the leading-coefficient law (given the top-coefficient table)

**Theorem.** Let `c(w)` be the coefficient of `k^{w-2}` in the weight-`w` part `W_w(k)` of the squircle
moment `M_k`, and let `d = w-1`. Assume the top-coefficient table `D_d` of the kernel generating
functions (Section 1; exact rational data for `d <= 9`, the tower law `T` of WTS/theory for `d >= 10`).
Then for **every** `w >= 2`,

    c(w)  =  eps_{w-2} * A_{w-1} / (2 (w-2)!) * [y^w] E(y),
    E(y)  =  Gamma(1+y/2) Gamma(1/2-y/2) / sqrt(pi),

with `eps_m = (-1)^floor((m+1)/2)` and `A_n` the Euler zigzag numbers (`sec x + tan x = sum A_n x^n/n!`).
Equivalently `sum_{w>=2} c(w) u^w = (1/2)[ u^2 (1-sinh u)/cosh^2 u  ⊙  E(u) ]` (Hadamard product).

The even-weight case (`d` odd) is **new**; the odd-weight case (`d` even) is the LEAD1 proof, refereed
and re-derived here (Section 6).  Every step is either proved in writing for all `d` or verified by
exact/high-precision code (`verify_lead2.py`; ranges stated).

**Notation.**  `lam = log 2`; `e_w := [y^w]E(y)` (`e_0 = 1, e_1 = lam`); zigzag sign
`delta_d := [x^d](sech x + tanh x) = eps_{d-1} A_d/d!` (Lemma 5.2);
`S = sqrt(1+u^2)`; `a = 1+sqrt2`, `b = sqrt2 - 1 = 1/a`.

## 1. The input: kernel data and the D-table (HYPOTHESIS for all d)

From WTS/theory (`kern_k.pkl`, exact): for each key `(i,j)`, `i+j = d`, `i >= j >= 0`, and each pole
`r in {-a, 0, -b, a}`, the coefficient of `1/(v-r)` in the level-`k` kernel is a polynomial in `k` of
degree `d` (Newton differences over `k = d+1..18` kill the `(d+1)`-st difference exactly; verified
`d <= 9`, `verify_lead2.py` part A). Its top coefficient `D_d^{(i,j)}(r)` obeys the table

    pole -a (all keys):                 D = (-1)^{j+1} binom(d,j) delta_d
    pole  a (odd d only):               D = (-1)^j    binom(d,j) delta_d  * (1 if j < (d-1)/2 else 1/2)
    pole -b (even d, j < d/2):          D = (-1)^j    binom(d,j) delta_d
    pole  0 (diagonals):                D = (-1)^j    binom(d,j) delta_d / 2
       (even d, j = d/2; odd d, j = (d-1)/2, which also carries pole -b with value (-1)^{j+1} binom(d,j) delta_d/2)
    everything else:                    D = 0           (all sqrt2-parts vanish)

**Status:** exact rational data for `d = 1..9` (part A; LEAD1 step_36/step_57). For `d >= 10` this is
the tower law `T`, unproved (attacked by TOW2).  **Everything below is conditional on this table.**

The leading coefficient is

    c(w)  =  (1/(d-1)!)  sum_{keys} sum_{poles r} D_d^{(i,j)}(r) N_{(i,j),r},
    N_{(i,j),r} = int_b^a dv/(v-r) log^i t log^j(1-t),
    t = (v-b)(v+a)/(4v),  1-t = (a-v)(v+b)/(4v),  d = w-1.                              (1)

We set `K_d := c(d+1)(d-1)!/delta_d`; the goal is `K_d = (d!/2) e_{d+1}` for all `d >= 1`.

## 2. The substitution u = (v - 1/v)/2 and the pole-combination measures

**Lemma 2.1 (exact algebra).** With `u = (v-1/v)/2`, `v = u + S` (positive branch, `v > 0` for
`u in (-1,1)`... `v` runs from `b` to `a` as `u` runs from `-1` to `1`):
  (a) `t = (1+u)/2`, `1-t = (1-u)/2`  (exactly);
  (b) `dv = v/S du`, and the pole combinations of the D-table become

      j < (d-1)/2 (odd d):  1/(v+a) - 1/(v-a)                 -> m_odd(u) du,  m_odd = [1 + sqrt2 u/S]/(1-u^2)
      j < d/2     (even d): 1/(v+a) - 1/(v+b)                 -> m_e(u) du,    m_e   = [1 - sqrt2/S]/(1-u^2)
      diag (even, j = d/2): 1/(v+a) - 1/(2v)                  -> m_diag du,    m_diag = (1-u)/2 m_e
      oddlast (j = (d-1)/2): 1/(v+a) - 1/(2(v-a)) + 1/(2(v+b)) - 1/(2v)
                                                              -> m_oddlast du, m_oddlast = m_odd - 1/(2(1-u))

  (c) `m_e` is even in `u` and extends to `-1/4` at `u = ±1`; `m_odd(u) + m_odd(-u) = 2/(1-u^2)` and

      m_odd(u) = 1/(1-u) - u m_e(u).                                                (2)

  (d) the two branches of `u = (v-1/v)/2` are `v_± = u ± S` with `v_+ v_- = -1` (deck transformation
      `v -> -1/v`); this is the origin of the pole pairing `{±a} <-> {∓b}`.

**Proof.** Rational-function identities in `Q(sqrt2)(u)[S]/(S^2 = 1+u^2)`, reduced symbolically
(sympy `radsimp` gives 0 for each difference; `verify_lead2.py` part B). For (c): `(1-u^2) m_e` even,
and `(1-u^2) m_e(u) = 1 - sqrt2/S -> 0` at `u = ±1` with rate `(u∓1)(±1)/4·...`, giving `-1/4`;
(2) is `(1+u - u + u sqrt2/S)/(1-u^2) = m_odd`.  ∎

Since `L1 := log t = log((1+u)/2)`, `L2 := log(1-t) = log((1-u)/2)`, equation (1) and Lemma 2.1 give,
with `M(i,j; m) := int_{-1}^1 L1^i L2^j m(u) du` (convergent for every occurring pair: near `u = 1`,
`m_odd ~ 1/(1-u)` but `L1 ~ -(1-u)·...` supplies convergence since `i >= 1` throughout; near `u = -1`
all measures are finite except `m_e, m_oddlast` which are still `O(1)`, and `L1^i` is integrable):

    K_d = sum_{j=0}^{floor((d-1)/2)} (-1)^{j+1} binom(d,j) M(d-j, j; m_{d,j})            (3)

with `m_{d,j} = m_e` (even d, `j < d/2`), `m_diag` (`j = d/2`), `m_odd` (odd d, `j < (d-1)/2`),
`m_oddlast` (`j = (d-1)/2`).  [Verified numerically from the raw kernel integrals (1), d = 1..9,
part F: the assembled `c(w)` matches the law to 91-102 digits.]

## 3. The three strip generating functions (Master Lemma)

Define on the open strip `-1 < Re s < 1` three absolutely convergent integrals

    h(w) := [(1+w)/sqrt(1+w^2) - 1]/w,   q(w) := (1+w)/(1+w^2)^{3/2},   g1(w) := (1-w)/(1+w) h(w),

    FF(s) := int_0^inf w^{-s} h(w) dw,   WW(s) := int_0^inf w^{-s} q(w) dw,   GG(s) := int_0^inf w^{-s} g1(w) dw.

Convergence: `h(w) -> 1` (`w->0`), `h(w) ~ 1/w^2` (`w->inf`); `q ~ 1` (`w->0`), `q ~ 1/w^2` (`w->inf`);
`g1 ~ 1` (`w->0`), `g1 ~ -1/w^2` (`w->inf`). Hence all three are analytic on the strip.

**Master Lemma.**
  (i)  `WW(s) = E(s) + E(-s)`.
  (ii) `FF(s) = (E(s) - E(-s))/s`.
  (iii) `GG(s) = 2 pi/sin(pi s) - (E(s) + E(-s))/s`.  In particular `GG` is odd.

**Proof.**
Step 1 (Beta form). `E(s) = Gamma(1+s/2)Gamma(1/2-s/2)/sqrt(pi) = (1/2) B(1+s/2, 1/2-s/2)`
since `Gamma(3/2) = sqrt(pi)/2`. For `-1 < Re s < 1`, `B(1+s/2,1/2-s/2) = int_0^1 t^{s/2}(1-t)^{-1/2-s/2} dt`.
Split at `t = 1/2`. On `(0,1/2)` substitute `t = w^2/(1+w^2)` (so `1-t = 1/(1+w^2)`, `dt = 2w/(1+w^2)^2 dw`,
`w : 0 -> 1`): the half integral is `2 int_0^1 w^{1+s}/(1+w^2)^{3/2} dw`. On `(1/2,1)` substitute
`1-t = w^2/(1+w^2)`: the half is `2 int_0^1 w^{-s}/(1+w^2)^{3/2} dw`. Hence

    E(s) = int_0^1 (w^{1+s} + w^{-s})/(1+w^2)^{3/2} dw.                                (4)

[Verified to 45 digits at several strip points, part E.]

Step 2 (proof of (i)). From (4), `E(s) + E(-s) = int_0^1 (w^s + w^{-s}) q(w) dw`. But
`q(1/w) = w^2 q(w)`, i.e. the measure `q(w) dw` is invariant under `w -> 1/w`; therefore
`int_1^inf w^{-s} q dw = int_0^1 w^s q dw`, and `WW(s) = int_0^1 w^{-s} q dw + int_0^1 w^s q dw = E(s)+E(-s)`. ∎

Step 3 (proof of (ii)). From (4), `E(s) - E(-s) = int_0^1 (w^s - w^{-s})(w-1)/(1+w^2)^{3/2} dw`.
The function `W(w) := 1 - (1+w)/sqrt(1+w^2) = -w h(w)` is a primitive of `(w-1)/(1+w^2)^{3/2}`
(differentiate). Integration by parts on `[eps, 1]`:

    int (w^s - w^{-s}) W' dw = [(w^s - w^{-s}) W] - s int (w^{s-1} + w^{-s-1}) W dw.

Boundary: at `w = 1`, `w^s - w^{-s} = 0`; at `w -> 0`, `W(w) ~ -w` so `(w^s-w^{-s})W = O(w^{1-|Re s|}) -> 0`
for `-1 < Re s < 1`. Hence `E(s)-E(-s) = s int_0^1 (w^s + w^{-s}) h(w) dw`. Since `h(1/w) = w^2 h(w)`
(same check as for `q`), `int_1^inf w^{-s} h dw = int_0^1 w^s h dw`, so the right side is
`s (F(s) + F(-s)) = s FF(s)` with `F(s) := int_0^1 w^{-s} h dw`. ∎
[IBP identity verified to 45 digits, part E.]

Step 4 (proof of (iii)). Write `GG(s) = int_0^inf w^{-s-1}(1-w) g(w) dw` with
`g(w) := 1/sqrt(1+w^2) - 1/(1+w)`  (`g(w) = O(w)` at `0`, `g(w) = O(w^{-2})` at `inf`; hence
`M_g(s) := int_0^inf w^{s-1} g(w) dw` converges exactly on `-1 < Re s < 1`). On `0 < Re s < 1` the two
pieces converge separately:

    int_0^inf w^{s-1}/sqrt(1+w^2) dw = (1/2) B(s/2, (1-s)/2) = Gamma(s/2)Gamma((1-s)/2)/(2 sqrt(pi)),
    int_0^inf w^{s-1}/(1+w) dw       = pi/sin(pi s),

(standard: substitute `w^2 = t/(1-t)`, resp. the classical Beta integral), so
`M_g(s) = Gamma(s/2)Gamma((1-s)/2)/(2 sqrt(pi)) - pi/sin(pi s)` on `0 < Re s < 1`, and both sides are
analytic on the full strip `-1 < Re s < 1` (at `s = 0` the `1/s` poles cancel; elsewhere clear), hence
equal there. On `-1 < Re s < 0`, both `M_g(s)` and `M_g(s+1)` converge and
`GG(-s) = M_g(s) - M_g(s+1)`; by analytic continuation this holds on the whole strip. Put `s -> -s`:

    GG(s) = M_g(-s) - M_g(1-s)
          = [Gamma(-s/2)Gamma((1+s)/2) - Gamma((1-s)/2)Gamma(s/2)]/(2 sqrt(pi)) + 2 pi/sin(pi s).

Using `Gamma(-s/2) = -(2/s) Gamma(1-s/2)`, `Gamma(s/2) = (2/s) Gamma(1+s/2)` and
`Gamma(1-s/2)Gamma((1+s)/2) = sqrt(pi) E(-s)`, `Gamma((1-s)/2)Gamma(1+s/2) = sqrt(pi) E(s)`
(the defining formula of `E`), the bracket is `-(2 sqrt(pi)/s)(E(s)+E(-s))`, giving (iii).
`GG` odd: `pi/sin` is odd and `(E(s)+E(-s))/s` is odd. ∎  [Verified to 45 digits, part E.]

**Remark (erratum to LEAD1's genfun.md).** The partial-fraction expansion of `F` is

    F(s) = sum_{n>=1} (-1)^n binom(2n,n)/4^n * 1/(2n-s) + sum_{n>=0} (-1)^n binom(2n,n)/4^n * 1/(2n+1-s),

i.e. the second sum starts at `n = 0` (LEAD1 printed `n >= 1` in both sums, dropping `1/(1-s)`;
the LEAD1 numerics used quadrature and are unaffected). Proof: expand
`h(w) = sum_{n>=1}(-1)^n c_n w^{2n-1} + sum_{n>=0}(-1)^n c_n w^{2n}` (`c_n = binom(2n,n)/4^n`); both
series converge uniformly on `[0,1]` (alternating, `c_n` decreasing to 0), so termwise integration
against `w^{-s}` (`Re s < 1`) is legitimate; convergence of the result for all `s` off the poles is
absolute (`c_n ~ 1/sqrt(pi n)`). Verified via `nsum` to ~50 digits (part E).
"""
open('proof_part1.md','w').write(part1)
print('wrote proof_part1.md', len(part1))
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: wrote proof_part1.md (setup, D-table, u-measures, Master Lemma)\n")
