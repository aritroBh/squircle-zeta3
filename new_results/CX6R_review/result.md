# CX6R: independent referee review of CX6 "an explicit rest bound"

**Under review:** `new_results/CX6_decay_rest/proof.md`, `verify_rest.py`, `verify_rest.log`.
**Claim (1):** for every `N = k+1 >= 100`, `|M_k - M_k^main| <= 10,600,000 * e^{-2 sqrt(N)}`.
**Claim (2):** `|M_k + (sqrt(pi)/N) e^{-sqrt(pi N)} cos sqrt(pi N)| <= 1.21e9 N^{-3/2} e^{-sqrt(pi N)}`.

## OVERALL VERDICT: **PROVED**

Every section was re-derived independently, every contour deformation and exchange was checked by
hand, and every identity and bound was tested numerically to 20-50 digits (residuals below).
The two constants were recomputed with mpmath independently of the arb run, and `verify_rest.py`
was re-run (passed in 13.2 s, output identical to `verify_rest.log`). No gap and no error found.
Four harmless remarks are listed at the end; none affects validity.

## What the proof does (as I understand it)

1. `M(eps) = -F(eps)/2` (paper eq:F, lem:holo) is rewritten by `g = t/(1-t)`, `x = -log g` as
   `-int_0^inf e^{-x}(1+e^{-x})^{-1-eps} C((1-eps)x) dx`, `C(z) = e^z(1+e^{-z}-sqrt(1+e^{-2z}))`  (3).
2. Kernel identity `C(z) = 1 - (2/pi) Im int_Gamma H K e^{-sz} ds` on `Gamma: s = 1/2 + r e^{i pi/4}`,
   `H(s) = (1/2)B(s/2,3/2)`, `K = 1/(2 cos(pi s/2))`  (5), proved by indenting Gamma onto the real
   axis above the odd-integer poles; the residue sum is the binomial series of `1 - C(z)`  (6).
3. Inserting (5) into (3) and extracting `[eps^k]` gives the exact representation
   `M_k = -J_k + (2/pi) Im int_Gamma H K P_k ds`  (7), `J_k = (-1)^k (log 2)^{k+1}/(k+1)!`,
   `P_k = (-1)^k/k! int_0^1 t^s (s log t + log(1+t))^k/(1+t) dt`.
   With `Q_k = s^k/(1+s)^{k+1}`, `A_k = P_k - Q_k`, `E = e^{i pi s/2}`, and
   `M_k^main = (2/pi) Im int_ray H E Q_k ds`, the exact rest decomposition is
   `M_k - M_main = -J_k + (2/pi) Im { int_Gamma H K A_k + int_Gamma H(K-E)Q_k - int_0^{1/2} H E Q_k }` (8).
4. Uniform bounds: `|A_k| <= 2 rho^k/|s+2|^{k+1} + 2 sum_{i>=1} (2^i/i!) rho^{k-i}/|s+i+1|^{k-i+1}` (12)
   via `L = -log t`, `beta(L) = log(1+e^{-L})`, `|beta| <= 2 e^{-Re L}`, per-term contour rotation;
   `|H|<=2`, `|K|<=2e^{-ar}`, `|K-E|<=2e^{-3ar}` on Gamma (13), `a = pi/(2 sqrt 2)`;
   `F(j,q;lam) <= (4/lam) e^{lam/2 + lam q/sqrt2} e^{-sqrt(3 lam j q/sqrt2)}` (14, AM-GM).
5. Assembly for `N >= 100`: `|M_k - M_main| <= C_* e^{-2 sqrt N}`,
   `C_* = 2 + T_0 + K_0 (e^{pi/4} + e^x + e^{ex})`, `x = 2e^{pi/4}`,
   `K_0 = 64 e^{a/2+pi/4}/(pi a)`, `T_0 = 32 e^{3a/2+3 pi/4}/(3 pi a)`  (15)-(19).
6. arb: `C_* = 10578907.82... < 10,600,000` (1); with `delta = 2 - sqrt(pi) > 0` and
   `sup_{u>0} u^3 e^{-delta u} = (3/(e delta))^3 = 114.096375...`, plus the supplied main-term
   theorem (254.23, certified by D13s_fix_recertify / D13u_close_h3: 254.2275176 <= 254.23),
   total `1,209,421,834.09 < 1.21e9` (2).

## Independent re-derivation notes (hand, all confirmed)

- (3): re-derived from scratch. Branch matches lem:holo (`sqrt(a^2+b^2) = a sqrt(1+(b/a)^2)` for
  `t <= 1/2`). Both displayed forms equal; absolutely convergent for `Re eps < 1`.
- (5): re-did the deformation. The ray through 0 lies ABOVE Gamma (line `Im s = Re s - 1/2`), so the
  wedge between them contains no pole; indentations above the odd integers are traversed clockwise
  (contribution `-i pi Res`); `Res_{2n+1} K = (-1)^{n+1}/pi` (checked against `d/ds cos`);
  vertical segment `|K(2J+iy)| = 1/(2 cosh(pi y/2))`, `|int| <= 1/(4J)` (uses `int_0^inf sech(pi y/2)dy = 1`);
  series (6) absolutely convergent at `z = 0` by (4); binomial identification uses
  `H(2n+1) = sqrt(pi) Gamma(n+1/2)/(4 (n+1)!)` and `sum binom(-1/2,n) u^{n+1}/(n+1) = 2(sqrt{1+u}-1)`.
- (7): Fubini domination `Re(s(1-eps)) >= (1 - sqrt2/4) Re s > 0` for `|eps| <= 1/4`
  (since `|s| <= sqrt2 Re s` on Gamma); Im-to-holomorphic complexification for coefficient
  extraction is legitimate; `J(eps) = (1-2^{-eps})/eps` gives the stated `J_k`.
- (8): the sign of `-int_0^{1/2}` is the crux; re-derived from the positively oriented closed
  contour `0 -> 1/2 -> out Gamma -> connector -> back the ray`. Sign CONFIRMED (my first guess of
  the opposite sign was wrong because Gamma is the LOWER of the two parallel rays).
- (9): origin residue of `H K P_k` is `J_k/2`; clockwise pi/4 arc gives `-i pi J_k/8`, hence
  `M_k = -(3/4) J_k + (2/pi) Im int_ray H K P_k`. DCY1 step_50's missing factor `log 2/(k+1)` confirmed.
- (10)-(12): re-derived; `beta(L) = int_0^1 e^{-L}/(1+v e^{-L}) dv`; `|1+v e^{-L}| >= 1/2` on
  `|arg L| <= pi/4` (cases `Re L </>= log 2`, using `log 2 < 1 < pi/2`); rotations to
  `-arg(s+i+1) in [-pi/4,0]` stay inside the wedge where (11) holds (`1+e^{-L}` has no zeros there
  and stays in the right half-plane, so the principal log is analytic); large-arc decay
  `Re((s+1+i)L) >= R|s+i+1| cos(pi/4) > 0`; factorial cancellation yields exactly (12).
- (13): `|1+e^{i pi s}|^2 = (e^{-y} - sin y)^2 + cos^2 y`, `y = 2ar`; case `y >= log 2`:
  `>= (1-e^{-y})^2 >= 1/4`; case `y <= log 2`: `>= cos^2(log 2) > 1/4` since `log 2 < 1 < pi/3`.
- (14): re-derived. Note `|s+q| >= (rho + q/sqrt2)/sqrt2` is the sharp inequality; the proof's
  `/2` is weaker but valid (Remark 1). AM-GM direction `e^{-X-Y} <= e^{-2 sqrt(XY)}` checked;
  reserved exponent accounting `lam r/4 + 3 lam t/4 - (3lam/4)(1/2 + q/sqrt2) = lam r` checked.
- (15)-(19): all elementary comparisons re-checked: `(i+1)(k-i) - 2(k-1) = (i-1)(k-i-2) >= 0` for
  `1 <= i <= k/2`; `(3 pi/2)(N-2) >= 4 N` needs `N >= 18`; `(N-1)/2 >= 2 sqrt N` needs `N >= 18`;
  `3^{-k} <= e^{-k/2}` (`log 3 > 1/2`); `N! >= 2^{N-1}`; all comfortably true for `N >= 100`
  (verified on a dense grid, step_14).
- (2): `u^3 e^{-delta u} <= (3/(e delta))^3` with `u = sqrt N >= 10`; `u* = 3/delta = 13.18 > 10`,
  so the global sup applies; `10.6e6 * 114.096375... + 254.23 = 1,209,421,834.09 < 1.21e9`.

## Evidence table (all numbers printed by my step files; "digits" = agreement)

| check | where | result | residual / ratio | digits |
|---|---|---|---|---|
| re-run `verify_rest.py` (venv python) | step_02 | PASSED in 13.2 s; output identical to verify_rest.log | — | — |
| (3) logit form vs original `M = -F/2`, eps in {0, .1, -.2, .3, -.4, .9} | step_03 | holds | <= 1.0e-40 | 40 |
| `M_0` from (3) vs table | step_03 | matches | — | 30 |
| `H(s) = (1/2)B(s/2,3/2)` vs definition integral | step_04 | holds | 2.7e-51 | 50 |
| (4) `|H| <= 1/x`, `<= sqrt(pi/2) x^{-3/2}` on 200-pt grid | step_04 | holds | max 0.256, 0.166 | — |
| (5) `C(z) = 1 - (2/pi) Im int_Gamma`, z in {0,.3,1,2.5,10} | step_04 | holds | <= 6.7e-51 | 50 |
| (6) residue series vs `1 - C(z)` | step_04 | holds | <= 6.8e-51 | 50 |
| `Res_{2n+1} K = (-1)^{n+1}/pi`, n=0..3 | step_04b | exact match | < 1e-19 | 19 |
| (13) on Gamma, 4000-pt grid | step_04b | `max|K|e^{ar} = max|(K-E)|e^{3ar} = 1.393 < 2`; `min|1+e^{i pi s}| = 0.7178 > 1/2` | — | — |
| vertical-segment formula + `int sech = 1` | step_04b | exact | 0.0 | — |
| **(7)** vs certified table, k = 0..8, 20 | step_05 | holds | <= 1.2e-41 | 41 |
| **(9)** (ray form, `-3/4 J_k`) vs table, k = 3, 8, 10 | step_06 | holds; wrong-sign arc variant fails at 5e-8 | <= 2.8e-42 | 42 |
| **(8)** full decomposition vs table & M_main, k = 3, 10, 30 | step_07 | holds | <= 2.9e-41 | 41 |
| rotation real-axis = ray for M_main (nsum) | step_08b | consistent | ~6e-14 (nsum-limited) | 13 |
| wedge `int_Gamma = int_ray - int_0^{1/2}` for `H E Q_k` | step_08c | holds | <= 2.7e-44 | 44 |
| real-axis = `int_Gamma + int_0^{1/2}` (chain closes to ray) | step_08c | holds exactly (diff = segment value) | — | 44 |
| (10) L-form of `A_k` vs t-form | step_09/09b | holds (t-side quadrature-limited) | <= 3.9e-24 | 23 |
| (12) pointwise grid k <= 30, r in [0,60] | step_10 | holds, worst ratio 0.3154 | — | — |
| rotated finite-sum `A_k = sum T_i` vs direct, k = 5..30 | step_11 | holds | <= 7.3e-45 | 45-67 |
| (12) termwise at k = 99, r in {0,4,12}; `T_0` vs exact series | step_12b | holds, worst ratio 0.5009; `T_0` err 1.8e-135 | — | — |
| (14) worst case over k = 99 regime (`j = k-i`, `q = i+1`) | step_12b | holds, worst ratio 0.00414 at `(j,q) = (99,1)` | — | — |
| (14) for `(k,1,3a)`, k = 99/200/400 | step_12b | holds, ratios 1.7e-4 / 2.3e-5 / 1.6e-6 | — | — |
| (14) for `(k,2,a)` (the (15) first piece), k = 99/200/400 | step_17 | holds, ratios 6.4e-4 / 1.3e-4 / 1.5e-5 | — | — |
| **(15)** actual `(2/pi)|Im int H K A_k|` vs bound, k = 99/200/400 | step_13/17 | holds, ratios 7.9e-7 / 1.7e-7 / 7.2e-9 | — | — |
| **(17)** actual vs `T_0 e^{-sqrt(9 pi k/4)}`, k = 99/200/400 | step_13/17 | holds, ratios 2.0e-6 / 1.2e-7 / 1.2e-8 | — | — |
| (18) exact integral `= 3^{-k}/k` | step_14 | holds | <= 7.9e-150 | 149 |
| (18) bound ratio at k = 99 | step_13 | 0.6118 < 1 | — | — |
| (16) ratio at k = 99 | step_16 | 9.6e-17 < 1 | — | — |
| all-N elementary chains (16),(18),`|J_k|`, saddle exponents | step_14 | hold for all N = 100..3000 dense + spots to 1e7 | — | — |
| sec.3 generating-function identity, eps in {.05,-.1,.2,-.24,.001} | step_16d | holds | <= 9.9e-25 | 24-30 |
| constants `K_0, T_0, C_*, delta, conversion, total` (mpmath, indep. of flint) | step_15 | match arb log to 30 digits | — | 30 |
| `C_* < 10,600,000` | step_15 | margin 21,092.18 | — | — |
| `total < 1.21e9` | step_15 | margin 578,165.91 | — | — |
| `sup u^3 e^{-delta u}` grid vs `(3/(e delta))^3` | step_14 | 114.09635859 <= 114.09637546 | — | — |
| `H(1) = pi/4`, recurrence `(2n+1)/(2n+4)` | step_15 | holds | <= 1.5e-101 | 100 |
| (20) ingredient `r|H(r omega)| <= sqrt 2` | step_15 | max 0.899 < 1.4142 | — | — |
| `P_0` digamma identity (sec.8) | step_15b | holds | <= 2.9e-23 | 23+ |
| **claim (1) scan**, 78 k-values in 99..400: `|M_k-M_main| e^{2 sqrt N}/C_*` | step_14 | max 1.06e-11 at k = 105 | — | — |
| **claim (2) scan**, same k's: ratio to `1.21e9 N^{-3/2} e^{-sqrt(pi N)}` | step_14 | max 1.36e-9 at k = 113 | — | — |
| k = 100/200/400 rests vs arb log | step_17 | identical to all 18 printed digits | < 1e-33 | 18+ |
| main-term input 254.23 certified (`254.2275176 <= 254.23`) | D13s_fix_recertify, D13u_close_h3 result.md | confirmed as input | — | — |

## Section-by-section verdicts

| section | contents | verdict |
|---|---|---|
| §1 logit form (3), H, (4) | CORRECT | hand re-derivation + 40-digit tests (step_03, step_04) |
| §2 shifted-contour identity (5),(6) | CORRECT | poles/indentations/orientation/closing segment all re-derived; 50-digit tests (step_04, step_04b) |
| §3 extraction (7), main term, rotation | CORRECT | 41-digit (7) tests at 11 k's; rotation chain to 44 digits (steps 05, 08, 08b, 08c, 16d) |
| §4 DCY1 endpoint correction, (9) | CORRECT | arc sign re-derived independently; 42-digit tests; wrong-sign variant fails (step_06) |
| §5 uniform bounds (10)-(14) | CORRECT | rotation validated to 45-67 digits; worst ratios 0.315 (12), 0.5009 (12 termwise k=99), 0.00414 (14) |
| §6 assembly (15)-(19) | CORRECT | (15),(17) actual-vs-bound at k=99/200/400; all-N chains on dense grid (steps 13, 14, 16, 17) |
| §7 constants + conversion to (2) | CORRECT | arb run reproduced; constants recomputed independently with mpmath (step_15) |
| §8 validation & reproduction | CORRECT | TSV as described; arb enclosures reproduced; my rests match to all 18 printed digits (steps 02, 17) |

## Remarks (none affects validity)

1. In (14) the proof uses `|s+q| >= (rho + q/sqrt 2)/2`; the sharp inequality has `/sqrt 2`.
   The proof's version is weaker, hence valid; it only inflates `C_*` slightly.
2. The arb `main_arb` integrand `(s/(1+s))^{N-2}/(1+s)^2 ...` (removable-pole cancellation) is the
   right form for `N >= 2`; it is used only at k = 100, 200, 400. Worth one line in the text.
3. My own Gauss-Legendre checks needed >= 96 nodes once integrands contain degree-~99 monomial
   factors (an early 48-node run produced a spurious termwise ratio; fixed and re-verified with
   level convergence checks). This concerns my review code only, not the proof or its arb verifier.
4. §4's commentary on DCY1 (Abel summability of the cut series, `c_i^{(m)} = (-1)^m e_i(1,...,1/m)`
   vs reversed index) was verified only through the printed k = 8/10/12 elementary terms in
   `verify_rest.log` (which I reproduced); it plays no role in the proof of (1) or (2).

## Files

- `cx6r_common.py` — shared high-precision Gauss-Legendre engine + kernels (H, K, E, Q_k, P_k, A_k, J_k, table loader).
- `step_01.py` inventory; `step_02.py` re-run of verify_rest.py; `step_03.py` (3); `step_04.py`, `step_04b.py` (4)-(6),(13), residues;
  `step_05.py` (7); `step_06.py` (9); `step_07.py` (8); `step_08.py`, `step_08b.py`, `step_08c.py` rotation chain;
  `step_09.py`, `step_09b.py` (10); `step_10.py` (12) grid; `step_11.py` rotated-sum validity;
  `step_12.py`, `step_12b.py` (12),(14) at k = 99 (12b supersedes the node-count bug in 12);
  `step_13.py` (15),(17),(18) at k = 99; `step_14.py` full-table scan of (1),(2) + all-N chains;
  `step_15.py` constants + small checks; `step_16.py`, `step_16b-d.py` generating-function identity + branch checks;
  `step_17.py` more worst cases + k = 100/200/400 vs arb log.

**Overall verdict: PROVED.** Both (1) and (2) are correct as stated, with all constants verified.
