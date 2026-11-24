# LEAD2 — proof of the leading-coefficient law (given the top-coefficient table)

**Theorem.** Let `c(w)` be the coefficient of `k^{w-2}` in the weight-`w` part `W_w(k)` of the squircle
moment `M_k`, and let `d = w-1`. Assume the top-coefficient table `D_d` of the kernel generating
functions (Section 1; exact rational data for `d <= 8`, and for `d = 9` a determined Newton fit on
`k = 10..18` with zero spare points that still matches the table everywhere; for `d >= 10` the tower
law `T` of WTS/theory, the one open input — attacked separately by TOW2).
Then for **every** `w >= 2`,

    c(w)  =  eps_{w-2} * A_{w-1} / (2 (w-2)!) * [y^w] E(y),
    E(y)  =  Gamma(1+y/2) Gamma(1/2-y/2) / sqrt(pi),

with `eps_m = (-1)^floor((m+1)/2)` and `A_n` the Euler zigzag numbers (`sec x + tan x = sum A_n x^n/n!`).
Equivalently `sum_{w>=2} c(w) u^w = (1/2)[ u^2 (1-sinh u)/cosh^2 u  ⊙  E(u) ]` (Hadamard product).

The even-weight case (`d` odd) is **new**; the odd-weight case (`d` even) is the LEAD1 proof, refereed
and re-derived here (Section 4).  Every step is either proved in writing for all `d` or verified by
exact/high-precision code (`verify_lead2.py`; ranges stated).

**Notation.**  `lam = log 2`; `e_w := [y^w]E(y)` (`e_0 = 1, e_1 = lam`); zigzag sign
`delta_d := [x^d](sech x + tanh x) = eps_{d-1} A_d/d!` (Lemma 6.1);
`S = sqrt(1+u^2)`; `a = 1+sqrt2`, `b = sqrt2 - 1 = 1/a`.

## 1. The input: kernel data and the D-table (HYPOTHESIS for all d)

From WTS/theory (`kern_k.pkl`, exact): for each key `(i,j)`, `i+j = d`, `i >= j >= 0`, and each pole
`r in {-a, 0, -b, a}`, the coefficient of `1/(v-r)` in the level-`k` kernel is a polynomial in `k` of
degree `d` (Newton differences over `k = d+1..18`: the `d`-th differences are constant and the
`(d+1)`-st vanish exactly; verified for `d <= 8` with spare points, while `d = 9` is a determined fit
on `k = 10..18` with zero spare points; `verify_lead2.py` part A). Its top coefficient
`D_d^{(i,j)}(r)` obeys the table

    pole -a (all keys):                 D = (-1)^{j+1} binom(d,j) delta_d
    pole  a (odd d only):               D = (-1)^j    binom(d,j) delta_d  * (1 if j < (d-1)/2 else 1/2)
    pole -b (even d, j < d/2):          D = (-1)^j    binom(d,j) delta_d
    pole  0 (diagonals):                D = (-1)^j    binom(d,j) delta_d / 2
       (even d, j = d/2; odd d, j = (d-1)/2, which also carries pole -b with value (-1)^{j+1} binom(d,j) delta_d/2)
    everything else:                    D = 0           (all sqrt2-parts vanish)

**Status:** exact rational data for `d = 1..8` (with spare Newton points) and `d = 9` (determined
fit, zero spare points, matches everywhere), part A; LEAD1 step_36/step_57. For `d >= 10` this is
the tower law `T`, unproved (attacked by TOW2).  **Everything below is conditional on this table.**

The leading coefficient is

    c(w)  =  (1/(d-1)!)  sum_{keys} sum_{poles r} D_d^{(i,j)}(r) N_{(i,j),r},
    N_{(i,j),r} = int_b^a dv/(v-r) log^i t log^j(1-t),
    t = (v-b)(v+a)/(4v),  1-t = (a-v)(v+b)/(4v),  d = w-1.                              (1)

We set `K_d := c(d+1)(d-1)!/delta_d`; the goal is `K_d = (d!/2) e_{d+1}` for all `d >= 1`.

## 2. The substitution u = (v - 1/v)/2 and the pole-combination measures

**Lemma 2.1 (exact algebra).** With `u = (v-1/v)/2`, `v = u + S` (the branch with `v > 0`;
as `u` runs from `-1` to `1`, `v` runs from `b` to `a`):
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
(sympy `radsimp` gives 0 for each difference; `verify_lead2.py` part B). For (c): `(1-u^2) m_e` is
even in `u`, and `1 - sqrt2/S = (S^2-2)/(S(S+sqrt2)) = (u^2-1)/(S(S+sqrt2))`, hence
`m_e(u) = -1/(S(S+sqrt2)) -> -1/4` at `u = ±1` (in particular `m_e` is continuous on the closed
interval); and (2) is the direct simplification `1/(1-u) - u(1 - sqrt2/S)/(1-u^2)
= ((1+u) - u + u sqrt2/S)/(1-u^2) = (1 + u sqrt2/S)/(1-u^2) = m_odd`.  ∎

Since `L1 := log t = log((1+u)/2)`, `L2 := log(1-t) = log((1-u)/2)`, equation (1) and Lemma 2.1 give,
with `M(i,j; m) := int_{-1}^1 L1^i L2^j m(u) du` (convergent for every occurring pair: at `u = -1`
all four measures are finite and `L1^i` has at worst an integrable log singularity; at `u = 1`,
`m_e, m_diag` are finite by (c), while `m_odd ~ 1/(1-u)`, `m_oddlast ~ 1/(2(1-u))` — but then
`L1 = log((1+u)/2) = O(1-u)`, and every occurring key has `i = d-j >= 1`, so
`L1^i L2^j m_odd = O((1-u)^{i-1} L2^j)` is integrable):

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
`M_g(s) := int_0^inf w^{s-1} g(w) dw` converges on `-1 < Re s < 2`). On `0 < Re s < 1` the two
pieces converge separately:

    int_0^inf w^{s-1}/sqrt(1+w^2) dw = (1/2) B(s/2, (1-s)/2) = Gamma(s/2)Gamma((1-s)/2)/(2 sqrt(pi)),
    int_0^inf w^{s-1}/(1+w) dw       = pi/sin(pi s),

(standard: substitute `w^2 = t/(1-t)`, resp. the classical Beta integral), so
`M_g(s) = Gamma(s/2)Gamma((1-s)/2)/(2 sqrt(pi)) - pi/sin(pi s)` on `0 < Re s < 1`, and both sides are
analytic on the full strip `-1 < Re s < 1` (at `s = 0` the `1/s` poles cancel; elsewhere clear),
hence equal there. On `-1 < Re s < 0`, both `M_g(s)` and `M_g(s+1)` converge and
`GG(-s) = M_g(s) - M_g(s+1)`; both sides are analytic on `-1 < Re s < 1` (`GG(-s)` because `GG` is
analytic on the strip; `M_g(s) - M_g(s+1)` by the convergence ranges), so the identity extends to
the whole strip. Put `s -> -s`:

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

## 4. Even d (odd weights): the LEAD1 argument, refereed

**Lemma 4.1 (parity collapse).** For even `d`,  `K_d = -(1/2) int_{-1}^1 (L1-L2)^d m_e(u) du`.

**Proof.** Write `T_j := (-1)^{j+1} binom(d,j) M(d-j,j; m_e)`. Since `m_e` is even, under `u -> -u`
(`L1 <-> L2`) `M(d-j,j;m_e) = M(j,d-j;m_e)`; and since `d` is even, `(-1)^{d-j+1} binom(d,d-j) =
(-1)^{j+1} binom(d,j)`. Hence `T_j = T_{d-j}` and `sum_{j<d/2} T_j = (sum_{j=0}^d T_j - T_{d/2})/2`.
The full sum is `sum_j (-1)^{j+1} binom(d,j) M(d-j,j;m_e) = -int sum_j binom(d,j) L1^{d-j}(-L2)^j m_e du
= -int (L1-L2)^d m_e du`. In (3) the middle key carries `m_diag`, not `m_e`, so

    K_d = sum_{j<d/2} T_j + (-1)^{d/2+1} binom(d,d/2) M(d/2,d/2; m_diag)
        = -(1/2) int (L1-L2)^d m_e du + (-1)^{d/2+1} binom(d,d/2) M(d/2,d/2; m_diag - m_e/2).

Since `m_diag - m_e/2 = -u m_e/2` is odd in `u` while `L1^{d/2} L2^{d/2}` is even, the correction
integral vanishes.  ∎  [Collapse checked to 60-70 digits for d = 2,4,6,8, part C.]

**Theorem 4.2.** For even `d`, `K_d = (d!/2) e_{d+1}`.

**Proof.** Lemma 4.1 gives `K_d = -(1/2) int (L1-L2)^d m_e du
= -(d!/2) [z^d] int_{-1}^1 e^{z(L1-L2)} m_e(u) du`. With `u = (1-w)/(1+w)` (so `du = -2/(1+w)^2 dw`,
`1-u^2 = 4w/(1+w)^2`, `sqrt2/S = (1+w)/sqrt(1+w^2)`) one gets, for any `f`,

    int_{-1}^1 f(u) m_e(u) du = -(1/2) int_0^inf f((1-w)/(1+w)) h(w) dw,                    (5)

hence `int_{-1}^1 ((1+u)/(1-u))^z m_e du = -(1/2) FF(z)` on the strip. Therefore for even `d`

    K_d/d! = -(1/2) [z^d] (-(1/2) FF(z)) = (1/4) [z^d] FF(z)
           = (1/4) [z^d] (E(z)-E(-z))/z        [Master Lemma (ii)]
           = (1/4) * 2 e_{d+1} = e_{d+1}/2,

because `(E(z)-E(-z))/z = 2 sum_{m odd} e_m z^{m-1}` and `m = d+1` is odd.  ∎

**Referee notes on the LEAD1 version (genfun.md 3.3).** Re-derived independently:
(a) the collapse (Lemma 4.1) — OK;
(b) the `w`-substitution and the inversion symmetry `h(1/w) = w^2 h(w)` reducing `int_0^inf` to
`2 int_0^1` — OK (subsumed in (5) + Step 3);
(c) `F(s) = int_0^1 w^{-s} h dw` converges for `Re s < 1`, its partial-fraction series converges
absolutely off the poles, termwise integration justified by uniform convergence of the alternating
binomial series — OK, with the printed erratum (second sum starts at `n = 0`), Remark in Section 3;
(d) `(F(s)+F(-s))/2 = (E(s)-E(-s))/(2s)` for `-1 < Re s < 1` — OK (this is Master Lemma (ii);
the proof via `E = (1/2)B`, split at `t = 1/2`, substitutions `t = w^2/(1+w^2)`, and IBP with
boundary vanishing exactly on `-1 < Re s < 1` re-checked line by line);
(e) coefficient extraction `[y^d]` at `y = 0` — OK (`F` analytic at 0, strip contains 0).
Verdict: **OK** (one documentation typo in the printed series; no mathematical gap).

## 5. Odd d (even weights): the new proof

Throughout, `d` is odd, `mid = (d-1)/2`.

### 5.1 Decomposition of the key sum

By (2), `m_odd = 1/(1-u) - u m_e`, and `m_oddlast = m_odd - 1/(2(1-u))`. Hence (3) splits as

    K_d = [ sum_{j<=mid} (-1)^{j+1} binom(d,j) M(d-j,j; 1/(1-u))
            - (1/2)(-1)^{mid+1} binom(d,mid) M(mid+1,mid; 1/(1-u)) ]                       (6a)
          - sum_{j<=mid} (-1)^{j+1} binom(d,j) M(d-j,j; u m_e).                             (6b)

Call the bracket `T_d` and the last line `-(u m_e part)`.

**Lemma 5.1 (the `u m_e` part collapses; PROVED).**
`sum_{j<=mid} (-1)^{j+1} binom(d,j) M(d-j,j; u m_e) = -(1/2) int_{-1}^1 (L1-L2)^d u m_e(u) du`,
and its (odd-`d`) generating function is `sum_{d odd} (...) z^d/d! = GG(z)/4`. Explicitly,

    (u m_e part)/d! = (1/4)[z^d] GG(z) = (1-2^{-d}) zeta(d+1) - e_{d+1}/2.                (7)

**Proof.** `u m_e(u)` is odd, so `M(d-j,j; u m_e) = -M(j,d-j; u m_e)`; and for `d` odd the
coefficient satisfies `(-1)^{d-j+1} binom(d,d-j) = (-1)^j binom(d,j)`. The two sign flips cancel:
`b_j := (-1)^{j+1} binom(d,j) M(d-j,j;u m_e)` obeys `b_{d-j} = b_j`. Since `d` is odd there is no
fixed point, so `sum_{j<=mid} b_j = (1/2) sum_{j=0}^d b_j
= -(1/2) int sum_j binom(d,j) L1^{d-j}(-L2)^j u m_e du = -(1/2) int (L1-L2)^d u m_e du`. Applying (5) with the extra factor
`u = (1-w)/(1+w)`: `int_{-1}^1 f(u) u m_e du = -(1/2) int_0^inf f((1-w)/(1+w)) g1(w) dw`, hence
`int_{-1}^1 ((1+u)/(1-u))^z u m_e du = -(1/2) GG(z)`. Extracting the odd powers of `z`
(`-(1/2) int sinh(z(L1-L2)) u m_e du = -(1/4)[-(1/2)GG(z) + (1/2)GG(-z)] = GG(z)/4`, using that
`GG` is odd by Master Lemma (iii)) gives the GF `GG(z)/4`. The explicit form is
`[z^d] GG(z) = 2 [z^d](pi/sin pi z) - [z^d](E(z)+E(-z))/z = 4(1-2^{-d}) zeta(d+1) - 2 e_{d+1}`
(odd `d`), using `[z^d] pi/sin(pi z) = 2(1-2^{-d}) zeta(d+1)` (from
`pi/sin(pi z) = 1/z + sum_{n>=1} (-1)^{n+1} (2^{2n}-2) B_{2n} pi^{2n}/(2n)! z^{2n-1}` and
`zeta(2n) = (2 pi)^{2n} |B_{2n}|/(2 (2n)!)`). ∎  [(7) verified to 51-72 digits for d = 1..11, part D.]

### 5.2 The Bose-Einstein-log bracket T_d

Define the convergent integrals (`p >= 1`)

    A_{p,q} := M(p,q; 1/(1-u)) = int_{-1}^1 L1^p L2^q/(1-u) du
             = (-1)^p int_0^inf r^p log^q(1-e^{-r})/(e^r-1) dr                            (8)

(the second form via `(1+u)/2 = e^{-r}`). Note `M(mid+1,mid; 1/(2(1-u))) = A_{mid+1,mid}/2`, so

    T_d = sum_{j=0}^{mid} (-1)^{j+1} binom(d,j) A_{d-j,j} - (1/2)(-1)^{mid+1} binom(d,mid) A_{mid+1,mid}.   (9)

**Lemma 5.2 (two-parameter generating function; PROVED).** For `Re y > 0`, `Re x > -1`,

    H(x,y) := sum_{p,q} (A_{p,q}/(p! q!)) x^p y^q = int_{-1}^1 e^{x L1 + y L2}/(1-u) du = B(1+x, y),

where the coefficient identity `A_{p,q}/(p!q!) = [x^p y^q] B(1+x,y)` holds rigorously for all `p >= 1`,
`q >= 0` (Laurent coefficient in `y`; `B(1+x,y)` has a simple pole `1/y`).

**Proof.** The integral equals `B(1+x,y)` by `u = 2s-1` (Section 3, Step 1 pattern) or by
`s = e^{-r}` in (8): `int_0^inf e^{-xr}(1-e^{-r})^y/(e^r-1) dr = int_0^1 s^x (1-s)^{y-1} ds = B(1+x,y)`.
For `p >= 1`, `a_p(y) := (-1)^p int_0^inf r^p (1-e^{-r})^{y-1} e^{-r} dr` is analytic on `Re y > -1`
(integrand `O(r^{p+Re y-1})` at 0, exponentially small at infinity) and `A_{p,q} = a_p^{(q)}(0)`;
while `[x^p] B(1+x,y) = (1/p!) int_0^1 (log s)^p (1-s)^{y-1} ds = a_p(y)/p!` by `s = e^{-r}`. So the
Taylor coefficients in `x` agree as analytic functions of `y` near `0`, and the `y`-Laurent
coefficients agree.  ∎

**Lemma 5.3 (functional equation and recursion; PROVED).** `y B(1+x,y) = x B(1+y,x)` (both sides
equal `P(x,y) := Gamma(1+x)Gamma(1+y)/Gamma(1+x+y)`, symmetric, analytic at 0, `P(0,0) = 1`).
Hence with `aa_{p,q} := A_{p,q}/(p!q!)`,

    aa_{p,q} = aa_{q+1, p-1}     (p >= 1, q >= 0).                                        (10)

**Proof.** `aa_{p,q} = [x^p y^q] B(1+x,y) = [x^p y^{q+1}] P`; `aa_{q+1,p-1} = [x^{q+1} y^p] P`;
equal by symmetry of `P`. ∎  [Recursion checked to 58 digits for p = 1..5, q = 0..4, part E.]

**Theorem 5.4 (the A-bracket; PROVED).** For odd `d`,

    T_d/d!  =  (1/2) [z^d] (pi/sin(pi z))  =  (1 - 2^{-d}) zeta(d+1).                     (11)

**Proof.** Work on the finite set `S_d = {(p,q): p >= 1, q >= 0, p+q = d}` with
`f(p,q) := (-1)^{q+1} aa_{p,q}`. The map `sigma(p,q) = (q+1, p-1)` is a bijection of `S_d`
(it is an involution: `sigma^2 = id`), preserving `f`: `aa` is preserved by (10), and
`(-1)^{q'+1} = (-1)^{p-1+1} = (-1)^p = (-1)^{q+1}` because `p+q = d` is odd. Its unique fixed point
is `D_* = ((d+1)/2, (d-1)/2)`. With `A_* = {(p,q) in S_d: p > q}` one checks
`sigma(A_*) = (S_d \ A_*) ∪ D_*` (cardinality check: `|S_d| = d`, `|A_*| = (d+1)/2`, and
`sigma(p,q) = D_*` iff `(p,q) = D_*`). Hence `2 sum_{A_*} f = sum_{S_d} f + f(D_*)`, and

    T_d/d! = sum_{A_*} f - (1/2) f(D_*) = (1/2) sum_{S_d} f.

Since `aa_{0,d} = [y^d] B(1,y) = [y^d] (1/y) = 0` for `d >= 1`, the sum extends to `p >= 0`:
`sum_{S_d} f = sum_{p+q = d} (-1)^{q+1} aa_{p,q} = -[z^d] B(1+z, -z)`. Finally

    B(1+z, -z) = Gamma(1+z) Gamma(-z) / Gamma(1) = z Gamma(z) * Gamma(1-z)/(-z) = -pi/sin(pi z)

(reflection formula), so `T_d/d! = (1/2)[z^d](pi/sin pi z)`. The Bernoulli expansion
`pi/sin(pi z) = 1/z + sum_{n>=1} (-1)^{n+1} (2^{2n}-2) B_{2n} pi^{2n}/(2n)! z^{2n-1}` with
`pi^{2n} |B_{2n}|/(2n)! = 2^{1-2n} zeta(2n)` gives, for `d = 2n-1`,
`[z^d] pi/sin(pi z) = 2(1-2^{-d}) zeta(d+1)`. ∎  [(11) verified to 63-72 digits for d = 1..11, part D.]

### 5.3 Assembly for odd d

**Theorem 5.5.** For odd `d`, `K_d = (d!/2) e_{d+1}`.

**Proof.** By (6), `K_d = T_d - (u m_e part)`. By Theorem 5.4, `T_d/d! = (1/2)[z^d](pi/sin pi z)`.
By Lemma 5.1 and Master Lemma (iii),

    (u m_e part)/d! = (1/4)[z^d] GG(z) = (1/4)( 2[z^d](pi/sin pi z) - [z^{d+1}] WW(z) ).

The `pi/sin(pi z)` contributions cancel exactly:

    K_d/d! = T_d/d! - (u m_e part)/d!
           = (1/2)[z^d](pi/sin pi z) - (1/2)[z^d](pi/sin pi z) + (1/4)[z^{d+1}] WW(z)
           = (1/4)[z^{d+1}](E(z)+E(-z))
           = (1/4)(e_{d+1} + (-1)^{d+1} e_{d+1}) = e_{d+1}/2,

since `d+1` is even. ∎  [Full chain key sum = T_d - piece2 = (d!/2) e_{d+1} verified to 51-72 digits
for d = 1..11, part C; T and piece2 closed forms to d = 11, part D.]

**Remark.** The finite parts that blocked LEAD1 are the identical `(1-2^{-d}) zeta(d+1) d!` terms in
`T_d` and in the `u m_e` part; the involution (10) exhibits the first as a single diagonal Beta value,
and the strip representation of `GG` exhibits the second as the `2 pi/sin(pi z)` pole term. They are
literally the same obstruction seen from the two sides, and they cancel.

## 6. Assembly: the law and its generating-function forms

**Lemma 6.1 (zigzag signs).** `delta_d := [x^d](sech x + tanh x) = eps_{d-1} A_d/d!` for all `d`,
and `eps_m A_{m+1}/m! = [z^m] (1-sinh z)/cosh^2 z`.

**Proof.** `(sec + tan)(ix) = sec(ix) + tan(ix) = sech x + i tanh x` (using `sec(ix) = sech x`,
`tan(ix) = i tanh x`). Comparing coefficients in `sum_d A_d (ix)^d/d! = sech x + i tanh x`:
`[x^d] sech x = Re(A_d i^d)/d! = A_d cos(pi d/2)/d!` (even `d`), `[x^d] tanh x = Im(A_d i^d)/d!
= A_d sin(pi d/2)/d!` (odd `d`). Since `sech` is even and `tanh` odd, `delta_d =
(-1)^{d/2} A_d/d!` (d even), `(-1)^{(d-1)/2} A_d/d!` (d odd); both cases equal
`eps_{d-1} A_d/d!` with `eps_m = (-1)^floor((m+1)/2)` (check `d mod 4`). For the second identity,
`(1-sinh z)/cosh^2 z = sech^2 z - sech z tanh z = d/dz(sech z + tanh z)`, so
`[z^m](1-sinh z)/cosh^2 z = (m+1) delta_{m+1} = (m+1) eps_m A_{m+1}/(m+1)! = eps_m A_{m+1}/m!`. ∎
[Both identities checked symbolically for `d, m <= 20`, part G.]

**Theorem 6.2 (the law).** Given the D-table of Section 1 for all `d`, for every `w >= 2`:

    c(w) = eps_{w-2} A_{w-1}/(2 (w-2)!) * e_w,   e_w = [y^w] E(y).

**Proof.** `c(d+1) = delta_d K_d/(d-1)!` (definition of `K_d`) `= delta_d (d!/2) e_{d+1}/(d-1)!`
(Theorems 4.2 and 5.5) `= delta_d (d/2) e_{d+1}`. With `delta_d = eps_{d-1} A_d/d!` (Lemma 6.1) the
prefactor is `eps_{d-1} A_d/(2 (d-1)!)`, and `w = d+1` gives the statement. ∎

**Corollary 6.3 (generating functions).**
  (a) `F(u) := sum_{w>=2} c(w) u^w = (1/2)[ u^2 (1-sinh u)/cosh^2 u  ⊙  E(u) ]`
      (Hadamard = coefficientwise product), since `eps_{w-2} A_{w-1}/(w-2)! =
      [u^{w-2}] ((1-sinh u)/cosh^2 u) = [u^w] (u^2 (1-sinh u)/cosh^2 u)` by Lemma 6.1.
  (b) The reduced GF is a pure Gamma product:
      `sum_{w>=2} [2(w-2)!/(eps_{w-2} A_{w-1})] c(w) u^w = E(u) - 1 - u log 2`.
  (c) `E(y) = 2^y Gamma(1-y) Gamma(1+y/2)/Gamma(1-y/2)`: Legendre duplication
      `Gamma(z)Gamma(z+1/2) = 2^{1-2z} sqrt(pi) Gamma(2z)` at `z = 1/2 - y/2` gives
      `Gamma(1/2-y/2) = 2^y sqrt(pi) Gamma(1-y)/Gamma(1-y/2)`, substitute in the definition of `E`.
      Hence with `log Gamma(1+t) = -gamma t + sum_{n>=2} (-1)^n zeta(n) t^n/n`,
      `log E = y log 2 + sum_{n>=2} s_n y^n/n`, `s_n = zeta(n)` (n even), `s_n = (1-2^{1-n}) zeta(n)`
      (n odd).  [All checked symbolically/numerically, part G.]

## 7. Verification guide (verify_lead2.py, ALL PASS, < 15 min)

| part | statement | check | digits/range |
|---|---|---|---|
| A | D-table = binomial x zigzag, d=1..9 | exact rationals from kern_1..18 (Newton differences) | exact |
| B | u-measures m_odd, m_e, m_diag, m_oddlast; t=(1+u)/2; m_odd=1/(1-u)-u m_e | sympy radsimp == 0 | symbolic |
| C | odd-d chain: key sum = T_d - piece2 = (d!/2) e_{d+1}; even-d: key sum = -(1/2) int log^d m_e = (d!/2) e_{d+1} | quadrature dps 70 | d odd 1..11, even 2..10; 51-72 digits |
| D | closed forms T_d = d!(1-2^{-d}) zeta(d+1); piece2 = T_d - (d!/2)e_{d+1} | quadrature dps 70 | d = 1..11 odd |
| E | Master Lemma: WW = E+E-, FF = (E-E-)/s, GG = 2pi/sin - (E+E-)/s; Beta sub (4); IBP; F-series; recursion (10) | quad + nsum dps 60 | strip points; 26-62 digits |
| F | full pipeline kernels -> c(w), w = 2..10 | quadrature dps 100 | 91-102 digits vs law |
| G | E(y) Gamma identities, duplication, log E expansion, zigzag signs, Hadamard GF form | sympy + series | symbolic / 30 digits |

All numbers printed by the scripts; no identity in this document is asserted beyond the ranges above
except those proved in writing for all d (Lemmas 2.1, 4.1, 5.1, 5.2, 5.3, Theorems 4.2, 5.4, 5.5, 6.1,
Master Lemma).

## Referee errata (LEAD2R, Kimi K3, 2026-10-27 22:1x: no mathematical gap; recorded by the M1 session, the text above is unchanged)
1. Eq. (3): the printed range j = 0..⌊(d−1)/2⌋ omits the even-d diagonal key; read j = 0..⌊d/2⌋ with the measure dictionary
   (Lemma 4.1 already adds the diagonal back, so the proof is unaffected).
2. §1: the kernel residues are polynomials in k of degree d−1 (the (d−1)-st differences are constant, the d-th vanish), not "degree d".
3. Lemma 5.2: read the displayed series equality coefficientwise; the analytic argument (a_p(y) analytic on Re y > −1 for p ≥ 1) is
   the valid content.
Table status after LEAD2R: exact for d ≤ 10 (levels k ≤ 21; spare points 3 at d = 9, 1 at d = 10), so the law holds for w ≤ 11 for the
exactly fitted polynomials. End-to-end check c(w) = law for w = 2..11: LEAD2R_referee_leading/claude_check_d9_d10.out.
