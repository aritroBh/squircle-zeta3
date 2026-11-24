# SM5: the tau^5 coefficient of the soft-max perimeter

**Outcome.**  c5 = (107/32) zeta(5) - zeta(2) zeta(3) = 1.48992283096334707297341790385558242735871285792...

- Numerically, c5 agrees with Richardson extrapolation of high-precision quadrature to about 55 digits (difference 4.6e-56).
- Exactly, Theorem 1 of `../LAY1_layer_gluing/result.md` at N = 5 reduces c5 to four moments of g(s) = log(1+e^{-2s}). The reflection identity (Lemma 5) and one elementary integral then give the closed form. Every intermediate integral is checked at 60 digits.
- c5 is of pure weight 5 in ordinary zeta values: Q-linear in zeta(5) and zeta(2)zeta(3), with no alternating-zeta, log 2 or Li_n(1/2) content.
- The method also reproduces the p-norm (squircle) value -sqrt2 (3/4 zeta(5) + 7/4 zeta(2) zeta(3)) exactly.

Setting: I(tau) = int_0^{2pi} h dtheta with h = (sum_v exp(v.n/tau))^tau, and I = 8 + (pi^2/3) tau^2 + 2 zeta(3) tau^3 + (zeta(4)/4) tau^4 + c5 tau^5 + O(tau^6).
Files (all in this folder): `sm5_numerics.py`, `sm5_richardson.py`, `sm5_pslq.py`, `sm5_derive.py`, `sm5_squircle_test.py`, `sm5_verify_integrals.py`, `sm5_h2_check.py`, each with an `*_output.txt`. Raw quadrature values are in `sm5_I_values.pkl`. The scripts need mpmath 1.3 and sympy 1.14 (I used a throwaway venv in this folder and deleted it).

---

## 1. Numerics (item 1)

**Method** (`sm5_numerics.py`). I(tau) = 8 int_0^{pi/4} (cos+sin)(1 + e^{-2 sin/tau} + e^{-2 cos/tau} + e^{-2(cos+sin)/tau})^tau dtheta, by mpmath tanh-sinh at dps 170 (maxdegree 12). The breakpoints are 0, tau, 2tau, 4tau, ... doubling out to pi/4, so they cluster at the layer theta = 0. The quadrature error estimate returned by mpmath is at most 9e-176 to 9e-182 at every tau used (table below).
tau = 1/n with n = 20 * 2^i, i = 0..17, so n runs from 20 to 2,621,440.

**Extrapolation** (`sm5_richardson.py`). Y(tau) = (I - 8 - (pi^2/3) tau^2 - 2 zeta(3) tau^3 - (zeta(4)/4) tau^4)/tau^5 = c5 + c6 tau + c7 tau^2 + ...
The remainder is a power series in tau because Theorem 1 gives no log terms and all other contributions are O(e^{-const/tau}). I apply Neville elimination with ratio 2. Order k removes tau, ..., tau^k.

Raw Y(tau), before extrapolation (this is the "rough fit" level):

| n = 1/tau | Y(tau) | quad err est |
|---|---|---|
| 20 | 1.50635981895740487114683337364 | 8.9e-182 |
| 320 | 1.49011160648397066437440087243 | 8.0e-174 |
| 5120 | 1.48993143287640992986551612534 | 8.0e-175 |
| 81920 | 1.48992335609778559840886785366 | 8.0e-176 |
| 2621440 | 1.48992284734860047242701676379 | 8.0e-177 |

Extrapolated value at order k, using the smallest-tau end of the table:

| k | value (45 digits) | abs. change from order k-1 |
|---|---|---|
| 5 | 1.48992283096334707297341790385550419031875854 | 9.2e-31 |
| 7 | 1.4899228309633470729734179038555824273564037 | 5.5e-39 |
| 9 | 1.48992283096334707297341790385558242735871286 | 3.2e-45 |
| 11 | 1.48992283096334707297341790385558242735871286 | 4.7e-50 |
| 13 | 1.48992283096334707297341790385558242735871286 | 1.5e-53 |
| 15 | 1.48992283096334707297341790385558242735871286 | 9.4e-56 |
| 17 | 1.48992283096334707297341790385558242735871286 | 2.1e-56 |

The change between successive odd orders falls by roughly four to eight decades per step until it flattens near 1e-56, where the cancellation in Y limits the digits. A coarser first run (n = 40 * 2^i, i <= 10, dps 110) gave 27 stable digits, which is consistent with the table.

**Result.** c5 = 1.489922830963347072973417903855582427358712857924387 (60-digit value in `sm5_c5_numeric.txt`). I claim **55 digits**: the successive-order change is 2e-56 at the last step, and the value agrees with the closed form of item 2 to 4.6e-56. That is far beyond the 30 requested.

---

## 2. Exact derivation (item 2)

**Data for the soft-max** (result.md of LAY1, Section 4): four kinks, H = T = 1, S = 0, A = 1 (so a_k = 0, P_0 = 1 and P_b = 0 for b >= 1, B_n = 0 for n >= 1). For s >= 0 the exact window formula is L = tau log(1 + e^{-y} + b) with y = 2 sin(tau s)/tau and b exponentially small.

`sm5_derive.py` (sympy) gives:

- y = 2s - tau^2 s^3/3 + tau^4 s^5/60 - ...
- With F(y) = log(1+e^{-y}) and L = tau F(y): g1 = log(1+e^{-2s}), g2 = g4 = 0, g3 = s^3/(3(1+e^{2s})), and g5 = s^5 (10 s e^{2s} - 3 e^{2s} - 3)/(180 (1+e^{2s})^2). g5 is not needed at N = 5: Theorem 1 uses g_k only for k <= N-1.
- Note g3 = -(s^3/6) g1', since g1' = -2/(1+e^{2s}).
- Tent coefficients, s >= 0: M0 = 1, M1 = s, M2 = -s^2/2, M3 = -s^3/6.
- Layer coefficients (exp(tau g1 + tau^3 g3) - 1): G1 = g1, G2 = g1^2/2, G3 = g1^3/6 + g3, G4 = g1^4/24 + g1 g3.

Theorem 1 at N = 5 (with B_5 = 0, and four kinks each contributing a symmetric s <-> -s pair) gives c5 = 4 int_R Phi_4 ds = 8 int_0^inf Phi_4 ds, where

  Phi_4 = M0 G4 + M1 G3 + M2 G2 + M3 G1
        = g^4/24 + s g^3/6 - s^2 g^2/4 + g g3 - s^3 g/6 + s g3     (g = g1).

**Reduction to moments.** Substitute g3 = -(s^3/6) g'. Integration by parts, with vanishing boundary terms, gives int s^k g^j g' = -(k/(j+1)) int s^{k-1} g^{j+1} for k >= 1. With gamma_{ab} := (1/(a! b!)) int_0^inf s^a g^b ds this yields

  int_0^inf Phi_4 ds = gamma_04 + gamma_13 + 3 gamma_31.

**Evaluation.**
- Lemma 5 (reflection) at total degree 4 in (lambda, mu). The right side's degree-4 part is (zeta(5)/2) mu^4 - (1/4)(zeta(2)zeta(3) + zeta(5)) lambda mu^3 + (1/8)(zeta(2)zeta(3) + zeta(5)) lambda^2 mu^2. Matching coefficients gives two relations:
  - gamma_22 = -3 gamma_31 + zeta(2)zeta(3)/16 + zeta(5)/16
  - gamma_04 = -gamma_13 + 2 gamma_31 - zeta(2)zeta(3)/8 + zeta(5)/8
  - gamma_13 and gamma_31 stay free. The unknown gamma_13 cancels from the integral.
- Hence int_0^inf Phi_4 ds = 5 gamma_31 + (zeta(5) - zeta(2)zeta(3))/8.
- gamma_31 is elementary: int_0^inf s^3 log(1+e^{-2s}) ds = sum_k (-1)^{k+1} (1/k) 6/(2k)^4 = (3/8) eta(5) = (45/128) zeta(5), so gamma_31 = (15/256) zeta(5).
- Therefore 8 int_0^inf Phi_4 ds = 40 gamma_31 + zeta(5) - zeta(2)zeta(3) = (75/32 + 1) zeta(5) - zeta(2)zeta(3):

  **c5 = (107/32) zeta(5) - zeta(2) zeta(3).**

**Numerical checks of the integrals** (`sm5_verify_integrals.py`, dps 60, quad):

| check | residual |
|---|---|
| gamma_31 - 15 zeta(5)/256 | 2.9e-62 |
| gamma_22 + 3 gamma_31 - zeta(2)zeta(3)/16 - zeta(5)/16 | 9.9e-62 |
| gamma_04 + gamma_13 - 2 gamma_31 + zeta(2)zeta(3)/8 - zeta(5)/8 | -8.0e-62 |
| 8 int_0^inf Phi_4 (direct, explicit g, g3) - closed form | -1.6e-61 |
| Richardson c5 of item 1 - closed form | 4.6e-56 |

The four numerical values are gamma_04 = 0.0015675048..., gamma_13 = 0.0024003921..., gamma_22 = 0.0061170496..., gamma_31 = 0.0607574856... (45 digits in the output file).

**Method test: squircle** (`sm5_squircle_test.py`). Setting: H = T = 1/sqrt2, window L = tau F(2 artanh(tan(tau s))/tau), y = 2s + (4/3) tau^2 s^3 + (4/3) tau^4 s^5, g1 = log(1+e^{-2s}), g2 = g4 = 0, g3 = -4 s^3/(3(1+e^{2s})) = (2/3) s^3 g1'. The same reduction and the same two reflection relations give

  c5(squircle) = -(7 sqrt2/4) zeta(2)zeta(3) - (3 sqrt2/4) zeta(5),

which equals the p-norm value -sqrt2 (3/4 zeta(5) + 7/4 zeta(2)zeta(3)), symbolic difference 0 (here gamma_13 and gamma_04 also drop out). The direct integral 8 int Phi_4 = -5.99340657187070660107725941712303200938264339396 agrees with the closed form to 0.0 at 50 digits (`sm5_verify_output.txt`). It is also consistent with LAY1's numerical limit of about -5.99. So the method was validated before the soft-max result was trusted.

**Hypothesis (H2) of Theorem 1 at N = 5 for the soft-max** (the remainder bound). I did verify it, by the argument in the style of LAY1 Section 4:

- Exactly, with x = tau|s| <= pi/4: L = tau log(1 + e^{-y} + b), y = 2 sin x/tau = 2|s| - delta, where delta = 2(x - sin x)/tau and b <= 2 e^{-sqrt2/tau}.
- Taylor with Lagrange remainder for x - sin x gives delta = tau^2|s|^3/3 + O(tau^4 |s|^5), with 0 <= delta <= 0.206|s|. Hence y >= 1.79|s| on the whole window.
- Second-order Taylor of F with |F'(y)| <= e^{-y} and |F''(y)| <= e^{-y} gives |tau F(y) - tau g1 - tau^3 g3| <= C tau^5 (1+|s|)^6 e^{-1.79|s|}. Here g3 = -(1/3)|s|^3 F'(2|s|).
- The b term is at most 2 tau e^{-sqrt2/tau}. By the non-decaying-error remark in (H2) it is <= C tau^5 e^{-c'|s|} on |s| <= pi/(4tau), since sqrt2 > pi/4.
- So (H2) holds at N = 5 with g1, g3 as above, g2 = g4 = 0 (the exponent has only odd powers of tau), and decay rate c < 1.79 (the polynomial weight costs a little).
- (H1) is trivial (A = 1). (H3) is vacuous, because W = pi/4 means the four windows tile T. Hypothesis (H0) is unchanged.
- Grid check (`sm5_h2_check.py`, dps 250, tau from 1 down to 0.005, 401 points in s on each window, exact L including b): the largest ratio of |L - tau g1 - tau^3 g3| to tau^5 (1+s)^6 e^{-s} is 0.29 (at tau = 0.5, s = 0, where the difference is the b term). So the bound with constant C = 1 holds on the grid.

The proof is a sketch at the level of detail of LAY1 Section 4 (constants not written out), not a line-by-line verification of every inequality; it is supported by the grid check above.

---

## 3. Comparison and PSLQ (item 3)

Item 1 and item 2 agree: the numerical c5 and (107/32) zeta(5) - zeta(2) zeta(3) differ by 4.6e-56, which is at the noise floor of the extrapolation.

PSLQ on the numerical c5 alone (`sm5_pslq.py`, 48 to 50 working digits, maxcoeff 10^6, independent of the exact derivation):

| basis | relation found |
|---|---|
| {c5, zeta(5), zeta(2)zeta(3)} | 32 c5 - 107 zeta(5) + 32 zeta(2)zeta(3) = 0 |
| plus zeta(4) | same relation, coefficient of zeta(4) = 0 |
| wide: {zeta(5), zeta(2)zeta(3), zeta(4)log2, zeta(3)log^2 2, zeta(2)log^3 2, log^5 2, Li5(1/2), Li4(1/2)log 2} | same relation, every extra coefficient 0 |

**Is c5 pure weight 5 in ordinary zeta values only?** Yes. c5 = (107/32) zeta(5) - zeta(2) zeta(3), a Q-combination of zeta(5) and zeta(2)zeta(3) alone, with no alternating-zeta, log 2, Li5(1/2) or Li4(1/2) log 2 terms. The exact derivation shows the same thing: every integral reduces to gamma_31 = (15/256) zeta(5) and the zeta(5), zeta(2)zeta(3) terms of Lemma 5. The pure statement here is "in the span of these two constants, as derived"; the PSLQ negative results on the wide basis are numerical evidence consistent with, and now explained by, the derivation.

Observation, no proof attempted: the squircle gives -sqrt2 (3/4 zeta(5) + 7/4 zeta(2)zeta(3)) and the soft-max (107/32) zeta(5) - zeta(2)zeta(3). Both are combinations of the same two constants; the soft-max rational coefficients are not those of the squircle scaled by a simple factor.

---

## 4. Status

| item | status |
|---|---|
| Theorem 1 at N = 5 applies to the soft-max, i.e. (H0)-(H3) at N = 5 | **PROVED** (short argument in the style of LAY1 Section 4, plus grid check; constants not written out) |
| c5 = (107/32) zeta(5) - zeta(2) zeta(3) for the soft-max, remainder O(tau^6) | **PROVED** (Theorem 1 at N = 5 + Lemma 5 + gamma_31 = 15 zeta(5)/256) |
| Reduction int_0^inf Phi_4 = gamma_04 + gamma_13 + 3 gamma_31 and the two degree-4 reflection relations | **PROVED** (sympy; integration by parts; Lemma 5) |
| Each moment gamma_ab (a+b = 4) and 8 int Phi_4 against the closed forms | **NUMERICAL**, residuals <= 1.6e-61 at dps 60 |
| c5 from quadrature + Richardson | **NUMERICAL**, 1.48992283096334707297341790385558242735871285792 (55 digits) |
| Agreement of numerical and exact c5 | **NUMERICAL**, difference 4.6e-56 |
| Squircle c5 = -sqrt2 (3/4 zeta(5) + 7/4 zeta(2)zeta(3)) from Theorem 1 at N = 5 (method test) | **PROVED** for the formula given the squircle (H2) at N = 5 (the squircle (H2) at N = 5 was NOT re-derived here; only N = 3, 4 are in LAY1); numerically consistent (-5.9934065718707...) |
| Squircle (H2) at N = 5 | **OPEN** (not examined; LAY1 Section 3 covers N = 3, 4 only) |
| PSLQ: c5 in span{zeta(5), zeta(2)zeta(3)} | **NUMERICAL** (relation found, coefficients 107/32 and -1, at 48 digits) |
| PSLQ on the weight-5 alternating basis: no extra constants | **NUMERICAL** (all extra coefficients 0) |
| c5 is pure weight 5 in ordinary zeta values | **PROVED** (follows from the closed form) |
| Tau^6 and higher coefficients | **OPEN** (not attempted) |
