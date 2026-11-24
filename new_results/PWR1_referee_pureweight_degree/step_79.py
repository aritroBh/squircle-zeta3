result = open('result.md','w')
result.write(r'''# PWR1 referee report: thm:pw (pure weight) and prop:degree (degree of the Fourier coefficients)

**Referee:** independent re-derivation of every identity, exact symbolic rebuilds (formal-Z engines, my own code:
`pwr1_engine.py`, `pwr1_rhombus.py`, hand-rolled rational series for the e-table), and fresh high-precision
quadratures (70-105 dps). Every number below was printed by code in this folder (steps step_15-*, `verify_pwr1.py`,
output `verify_pwr1.out`; 21/21 checks PASS, runtime 106 s).

**Objects.**
- (A) thm:pw: log(P(p)/(4 sqrt 2)) = sum_{k>=2} w_k tau^k, w_k in Q[zeta(2),..,zeta(k)] homogeneous of weight k,
  tau = 1-1/p; same for every rhombus (w_k in Q[q^2, zeta(2..k)]); and S_n(p)/(S_n(1) Rat_n(tau)) = exp(sum w_k(n) tau^k)
  in every dimension n. Written proof reviewed: research/08_theory/D16_pure_weight/result.md (Proof (mechanism) M1-M3,
  single-appearance lemma, d16_engine.py) + H1 result.md sec.3 (T6) for general n.
- (B) prop:degree: [tau^n] (-1)^k (pi/sqrt2) a_k(h_p) = Theta_n(k^2), deg exactly floor(n/2)-1, alternating-MZV
  coefficients of weight n; leading coeff for n=2m: (-1)^{m-1}(2^{2m-1}-1) zeta(2m). Written proof reviewed:
  research/08_theory/K3_fourier_all_orders/result.md sec.1.1-1.3.

---

## A. thm:pw -- step-by-step verdict table

| step | verdict | what I did / printed evidence |
|---|---|---|
| H1 T1 representation S_n = 2^n eps^{n-1} Gamma(eps)^n/Gamma(1+(n-1)eps) E[(sum G_i^{2 eta})^{1/2}], and D16's B_n-form with Gamma(1+a), a=-(n-1)tau/n | CORRECT | re-derived (octant formula + Dirichlet/Lukacs); equivalence of the two forms proved exactly (size-bias shift: E_{1+a}[F(G) prod G_j^{-tau/n}] = (Gamma(eps)/Gamma(1+a))^n E_eps[F], and 1+a-tau/n = eps). Numeric: quadrature n=2 (P, 70-100 dps), n=3 (simplex integral, 50 dps) both match the resulting series (remainder tests below) |
| M1: explicit moment E[prod(e^{2 tau m_i}-1)] = sum_S (-1)^{j-s} prod_v Gamma(1+a+2 tau(nu_v-s/n))/Gamma(1+a)^n, counting (n)_r s!/(prod lambda_i! prod mult!), E[X^j] assembly | CORRECT | re-derived line by line (independence of G_i; equality-pattern count; S-subset expansion). My independent exact engine implementing ONLY this formula reproduces D16's tau_series_2d.json (P/sqrt2, logP, B2) through tau^10 -- difference exactly 0, formal Z's; reproduces H1 routeA jsons (n=2,3) and logform_general_n.out logB through tau^8 -- diff exactly 0; reproduces D11's exact f_n (independent third party, eta=-tau) through eta^8 -- diff exactly 0 |
| M2: log Gamma(1+x) = -gamma x + sum (-1)^k zeta(k) x^k/k; every x = tau x rational, so zeta(k) enters at exactly tau^k; products/exp/log preserve grading | CORRECT | verified exactly by engine: coefficient of tau^k is homogeneous of weight k (wt zeta(j)=j) for every k<=10 (n=2) and k<=8 (n=3..6) -- printed True at every order |
| M3a: gamma cancels inside each log-moment (sum_v (nu_v - s/n) = 0) | CORRECT | re-derived (sum_v nu_v = s); engine tracks a formal GAM symbol: GAM-coefficient exactly 0 at every order tau^1..tau^10 (n=2,3), tau^8 (n=4,5,6) |
| M3b: gamma cancels in the flat-face factor; integer shifts give only rationals | CORRECT | re-derived: n log Gamma(1+a) contributes -n gamma a = +gamma(n-1) tau; -log Gamma(1-(n-1)tau) contributes -gamma(n-1) tau; sum 0. Gamma(n-(n-1)tau) = Gamma(1-(n-1)tau) prod_{i=1}^{n-1}(i-(n-1)tau), log matches log Gamma(n) via sum log i; harmonic sums stay rational |
| Rat_n(tau) = (1-tau)^{n-1}/prod_{i=1}^{n-1}(1-(n-1)tau/i) = prod_{i=1}^{n-2} i/(i-(n-1-i)d) | CORRECT | re-derived both forms and their equivalence; engine: logRat(2) = 0 exactly (Rat_2 = 1 -- the perimeter needs NO prefactor), logRat(3) = log(1-tau)-log(1-2 tau), i.e. Rat_3 = (1-tau)/(1-2 tau) = 1/(2-p) (H1's prefactor) |
| Lambda_n = -sum_{k>=2} zeta(k)(n-1)^k(1-n^{1-k}) tau^k/k | CORRECT | re-derived from n log Gamma(1+a) - log Gamma(1-(n-1)tau); matches engine exactly |
| w_1 = 0 (log starts at tau^2) | CORRECT | E[X] has no tau^1 term (gamma-part vanishes identically, zeta-parts start at k=2); engine: w_1 = 0 for n=2..6 |
| B_n at Z=0 is 1 (no rational terms in log B_n) | CORRECT | re-derived: E[X^j] at Z=0 gives n^{-j} sum_s binom(j,s)(-1)^{j-s} n^{j-s} n^s = (1-1)^j = 0 for j>=1 (each moment M=1 at Z=0). This kills rational/weight-0 pollution of log B_n |
| Single-appearance lemma: zeta(k) occurs in log only at tau^k, rationally | CORRECT | re-derived: each log-moment's Z_k-part is (rational) Z_k tau^k; exp linearizes single Z's at the same order; since B_n-1 has no weight-0 terms, products in log B_n have weight <= degree with equality only via the linear term, so a single Z_k can appear only at tau^k. Engine: u_k(2) = 1/4, 1/4, 1/32, -3/16, 91/192, 129/64, -5039/1024, -22805/768, 561331/5120 (k=2..10), matching D16's list exactly |
| n-dimensional part (H1 T6): A_2(n) = -(n-1)(n^2-2n-1)/(2n), A_3(n) = -(n^2-1)((n-1)^3-2)/(3n^2); W_4 | CORRECT | engine n=3..6 matches: e.g. w_2(4) = -21/8 Z2 = A_2(4) zeta(2), w_3(4) = -125/16 Z3 = A_3(4) zeta(3); also D16's u_2(n), u_3(n) closed forms match (n=3,4 checked exactly) |
| Rhombus |x/a|^p+|y/b|^p=1 | CORRECT (with a re-derivation note) | D11c's integral rep P = 4r(1+eta) int_0^1 (t(1-t))^{eta/2} sqrt(cosh eta l + q sinh eta l) dt verified against direct arclength quadrature at 50 digits, 4 parameter sets (diffs 0, -2.1e-50). My own exact q-symbolic engine (operator/beta form, re-derived): log(P/(4r)) through tau^6 is exactly pure weight with coefficients in Q[q^2]; the grading mechanism re-proved (the (-tau)^{2m} prefactor + 2m mu-derivatives give total excess tau-degree - weight = 0; gamma cancels as in M3; prefactor (1-tau) cancels the rational -log(1-tau) from Gamma(2-tau), hence Rat=1). w_2 = (1-q^2) zeta(2)/4, w_3 = (1-q^2) zeta(3)/4 exactly; quadrature fits (a,b)=(1,2),(2,3),(1,3): w2 to 11 digits, w3 to 7.5-8 digits, w4 to 3.8-5.4 digits, w5 to 2-3.3 digits (asymptotics-limited). kappa_3 = 1-q^2 = 4a^2b^2/(a^2+b^2)^2 consistent with D11c |
| Paper's displayed series (P/sqrt2 = 4 + zeta2 tau^2 + zeta3 tau^3 - 33/16 zeta4 tau^4 - (3/4 zeta5 + 7/4 zeta2 zeta3) tau^5 + ...) | CORRECT | matches my engine's formal-Z series exactly through tau^10 (cross-checked against D16 json and D11) |
| Interchanges of limits/sums/integrals | CORRECT | all manipulations are at the level of formal Taylor series in tau: B_n expansion needs only j<=N moments for tau^N; each moment is a finite Gamma-product; log Gamma expansions used only at arguments tau x rational (|x|<1 for tau small); expectation/derivative interchange justified by E[e^{c|m_i|}] < inf for small c (log-moments of Gamma exist to all orders). Series is asymptotic only (radius 0), but the theorem is explicitly "to every finite order" |

### thm:pw numerical evidence (all printed this session; computations at 50-105 digits)

| test | result |
|---|---|
| exact engine vs D16 tau_series_2d.json (P/sqrt2, logP, B2), formal Z | difference exactly 0 through tau^10 |
| exact engine vs H1 routeA (n=2,3) and logB table | difference exactly 0 through tau^8 |
| exact engine (n=2, tau=-eta) vs D11 exact f_n | difference exactly 0 through eta^8 |
| gamma coefficient, n=2..6 | exactly 0 at every order (k<=10 for n=2,3; k<=8 for n=4,5,6) |
| weight homogeneity, n=2..6 | True at every order |
| quadrature remainder test, n=2, logP at 9 taus in [0.0005,0.006] (100 dps): R(tau)=logP-sum_{2..8} w_k tau^k | R/(w9 tau^9) = 0.99704 at tau=0.0005 -> 1; (R-w9 tau^9)/(w10 tau^10) = 1.00277 -> 1 (validates w2..w10 collectively to ~30 digits: R at tau=0.0005 is -4.56e-28, matching prediction to 0.3%) |
| coefficient extraction (LSQ, 13 points, k=2..12) vs exact | w2: 20.5 digits; w3: 16.8; w4: 13.9; w5: 11.2; w6: 9.0; w7: 6.9; w8: 5.2; w9: 3.7; w10: 2.4 (asymptotics-limited, as expected for a divergent series) |
| purity residuals (fitted coeffs minus exact pure-weight values) | w2 - zeta2/4 = -1.17e-21; w4 + 19 pi^4/2880 = -8.2e-15; w5 + (zeta2 zeta3/2 + 3 zeta5/16) = +7.3e-12 -- excludes gamma/rational/lower-weight pollution at those levels |
| n=3 quadrature (50 dps, 5 taus) | remainder R/(w9(3) tau^9) = 1.0016 at tau=0.001 -> 1; LSQ fit w2..w8 matches exact |
| rhombus quadrature vs exact engine | see rhombus row above |

**Overall verdict, thm:pw: PROVED.**

Repairs (documentation level only; the proof is correct as written):
- (R1, MINOR FIX) In K3's result.md the displayed formula I_{m,j} = j! m! sum_{l>=j} (-1)^{l-j} s(l,j)/(l! l^{m+1}) is correct only if s(l,j) denotes the UNSIGNED Stirling numbers of the first kind (with signed s(l,j) the factor (-1)^{l-j} must be removed). State the convention. (Their code used quadrature, so no numerical consequence.)
- (R2, MINOR FIX) In D16's result.md, "4 u_k(2) = kappa_k exactly, k=3..11" holds for ODD k as stated/checked there, but a reader could over-generalize: for even k the eta- vs tau-sign convention flips the sign (kappa_2 = -1 while 4 u_2(2) = +1; kappa_4 = -1/8 while 4 u_4(2) = +1/8). State the sign convention explicitly.

---

## B. prop:degree -- step-by-step verdict table

| step | verdict | what I did / printed evidence |
|---|---|---|
| a_k - A_k = (-1)^k B(tau,2k), B(tau,kappa) = (8/pi) int_0^{pi/4} (h_p-h_1) cos(kappa(pi/2-2theta)) dtheta, B even entire in kappa | CORRECT | re-derived (h_p symmetric under theta -> pi/2 - theta; cos 4k theta = (-1)^k cos(2k(pi/2-2theta)) on [0,pi/4]); B even (cos) and entire (compact interval, entire integrand). A_k = -(-1)^k 4 sqrt2/(pi(16k^2-1)) verified numerically: |a_k(h_1) - A_k| = 7.4e-13 for k=0..5 |
| B(tau,kappa) = sum_m D(tau;m)(-kappa^2)^m/(2m)!, termwise in kappa | CORRECT | compact interval + uniform convergence of the cos series: termwise integration fine |
| tau-coefficient extraction commutes with the m-sum; uniformity in k | CORRECT (the flagged point resolves cleanly) | [tau^N] B(tau,2k) = sum_m [tau^N] D(tau;m) (-4k^2)^m/(2m)! and the sum is FINITE: D(tau;m) = O(tau^{2m+2}), so only m <= N/2 - 1 contribute. No uniformity-in-k issue arises at any fixed order N |
| layer identity h_p - h_1 = (1+e^{-2y})^{-1/2} (e^{tau L(y/tau)} - 1), y = -log tan theta, L(s)=log(1+e^{-s}) | CORRECT | re-derived exactly from h_p = cos theta (1+tan^{1/tau} theta)^tau and tan theta = e^{-y} |
| Watson expansion: D(tau;m) = (8/pi) tau int_0^inf H(tau s) gd(tau s)^{2m} sum_j tau^j L(s)^j/j! ds, H(x) = e^{-x}(1+e^{-2x})^{-3/2} = 2^{-3/2} e^{x/2} sech^{3/2} x; Taylor coefficients = asymptotic series | CORRECT | re-derived (theta -> y -> y = tau s; d theta = -e^{-y}/(1+e^{-2y}) dy; pi/2-2theta = gd(y) since gd(x)=2 arctan e^x - pi/2). Justification: standard Watson-lemma situation -- for y >= delta the integrand is flat in tau (e^{-y/tau} scale), in the layer everything is analytic in tau uniformly on compact s-sets, tails bounded by e^{-s}; moments int s^r L^j finite. Hence D(tau;m) is smooth and its Taylor coefficients are the stated ones |
| D~_{n,m} = 2 sum_{j=1}^{n-1-2m} e_{m,n-1-j} I_{n-1-j,j}/j!, e_{r,m} = 2^{3/2}[x^m] gd^{2r} H in Q | CORRECT | re-derived; my own exact rebuild of the e-table (hand-rolled rational series, gd' = sech via Euler numbers) matches all 50 entries of K3's e_rm.pkl exactly; e_{r,2r}=1 and e_{r,2r+1}=1/2 proved (gd ~ x, H(0)=2^{-3/2}, odd part of H at 0 = (x/2) 2^{-3/2}); e_{r,2r+2} = -5/8 - r/3 verified r<=5 |
| I_{m,j} = int_0^inf s^m log^j(1+e^{-s}) ds = alternating MZV of weight m+j+1 | CORRECT | re-derived via log^j(1+x) = j! sum_l s(l,j) x^l/l! (signed Stirling) and int s^m e^{-ls} = m!/l^{m+1}: I_{m,j} = j! m! sum_l s(l,j)/(l! l^{m+1}), the sum being Li_{m+2,{1}^{j-1}}(-1,..) of weight m+j+1. Quadrature (105 dps) = Stirling series (m>=1: agreement 1e-7..3e-15, truncation-limited; direct closed forms: I_{m,1} = m! eta(m+2) exact, I_{0,2} = zeta(3)/4 exact) |
| degree <= floor(n/2)-1 | CORRECT | j >= 1 in the D~ formula forces D_{n,m} = 0 for n < 2m+2 |
| degree attained, even n: D~_{2n,n-1} = 2(2n-2)! eta(2n) > 0 | CORRECT (PROVED) | only j=1 contributes; e_{n-1,2n-2}=1; I_{2n-2,1} = (2n-2)! eta(2n). Numerically: diffs 0, 0, 2.2e-103, 1.4e-101, 9.1e-100 for n=2..6 |
| degree attained, odd n: D~_{2n+1,n-1} = (2n-2)![zeta(2n+1) + 2 sum_{j=1}^{n-1} eta(2j) zeta(2n+1-2j)] > 0 | CORRECT (PROVED) | re-derived: D~ = (2n-1)! eta(2n+1) + 2(2n-2)![eta(2n+1) - A_{2n}] with I_{m,2} = 2m![eta(m+3)-A_{m+2}], A_s = sum (-1)^{l-1} H_l/l^s; classical evaluation 2 A_{2m} = (2m+1) eta(2m+1) - zeta(2m+1) - 2 sum eta(2j) zeta(..) verified to <=1.1e-106 for m=1..4; the closed form then matches D~ to <=8.7e-106 (n=2), 1.4e-104 (n=3), 8.9e-103 (n=4), 0 (n=1,5). Positivity manifest (all terms > 0), so the degree is attained |
| leading coefficient for n=2m: (-1)^{m-1}(2^{2m-1}-1) zeta(2m) | CORRECT | algebraic identity 2(-4)^{m-1}(1-2^{1-2m}) = (-1)^{m-1}(2^{2m-1}-1) verified exactly (sympy, m=1..8); numeric: matches D~-built top to <=2.2e-103, m=1..6 |
| coefficients are alternating MZVs of weight n | CORRECT | each I_{r,j} in D~_{n,m} has r+j = n-1, weight r+j+1 = n, alternating |
| Theta_n(k^2) polynomials, n<=8 | CORRECT | my independent D~-build (quadrature I's + my e-table) vs K3's PSLQ-identified zeta-product Qpolys: max |diff| over u in {0,1,2,5,7} = 1.6e-51 (n=2), 1.1e-51 (n=3), 4.5e-49 (n=4), 3.9e-49 (n=5), 9.2e-48 (n=6), 4.0e-47 (n=7), 2.4e-45 (n=8). Paper's displayed Theta_5 = -(3/4+4u) zeta5 - (7/4+4u) zeta2 zeta3 matches |
| a_k quadrature extraction (independent) | CORRECT | fresh a_k (70-90 dps) minus A_k, LSQ in tau; extracted (-1)^k (pi/sqrt2) [tau^n] a_k vs my Theta_n(k^2), k=0..6: rel.err n=2: 6.5e-12, n=3: 3.1e-8, n=4: 1.7e-7, n=5: 6.5e-5, n=6: 1.8e-4 (fresh 5-pt run: n=2: 8.9e-10, n=3: 2.9e-6, n=4: 2.4e-5). Degree test: overfit u^{floor(n/2)} coefficient negligible relative to scale at every n<=6 |

**Overall verdict, prop:degree: PROVED.**

Repairs: only (R1) above touches this proof (the Stirling sign convention in K3's displayed I_{m,j} formula).

---

## Verify script

`verify_pwr1.py` (runtime 106 s < 5 min): 21 checks A1-A11 (thm:pw) and B1-B9 (prop:degree), all PASS;
output in `verify_pwr1.out`. Referee engines: `pwr1_engine.py` (formal-Z exact series for S_n, moments,
Rat_n, Lambda_n, weight/gamma tracking), `pwr1_rhombus.py` (q-symbolic rhombus engine). Data: quad_logP_n2*.json,
quad_logS_n3.json, quad_rhombus.json, quad_ak.json, w_exact_n2_tau10.json, w_exact_n3_tau10.json, e_rm_mine.pkl,
I_mj_quad.json, fit_w2_w12.json, tops_check.json, theta_extracted.json.
''')
result.close()
open('progress.log','a').write(
"step_67-69: Theta_n D~-build vs K3 Qpolys n<=8: 1.6e-51..2.4e-45; leading coeff m=1..6 exact; e-table extended r<=5\n"
"step_70-72: a_k quadrature 90dps (k=0..6, 8 taus), A_k verified 7.4e-13; extraction Theta_n(k^2): rel.err 6.5e-12..1.8e-4 (n<=6); degree overfit negligible; n=7,8 1.2e-2,9.6e-2 (asymptotics)\n"
"step_73-75: verify_pwr1.py written; run: 21/21 PASS, 106s; result.md completed\n")
print("result.md written,", len(open('result.md').read()), "chars; progress.log updated")
import os
print(sorted(f for f in os.listdir('.') if not f.startswith('step_')))
