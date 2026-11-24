result = '''# WT2r: adversarial review of WT2b's claimed proof of the weight-two conjecture for all k

Reviewer: independent re-derivation + re-verification in this directory (rev_01..rev_16), using the
WT2b machinery (fastred/genintegrand/buildintegrand copied here) and the paper's own scripts
(paper/v2/code/mk_integrands.py). Every claim below was printed by code in this directory.

## VERDICT: PROOF COMPLETE
(modulo the *setting*: the paper's master integrand Phi and its reduction algorithm ibp_reduce,
which are inputs to WT2b, not part of its claim. The ε-expansion genintegrand.py is PROVED from
the master formula by exact algebra and was re-verified against mk_integrands.t_integrand for
k=4..7, all monomials, exact.)

One documentation bug found (not a math error): genintegrand.py's docstring writes the sqrt-part
coefficient with a stale minus sign ("-q/(2 t(1-t) q2^k) * S_m(t)"); the code computes the correct
+ sign, matching the paper's t_integrand exactly (rev_14). One cosmetic off-by-one: WT2b's (V2)
bound check would flag Gm_(1,0) (bound 2d-3 = -1 for d=1); the stored Gm is the zero function.

## Review question 1 — Lemma A (pole order <= 2d-2 at +-i): PROVED FOR ALL k
- DL1, DL2 have simple zeros at v=+-i (symbolic: values 0; DL1'(i)=(-1-i)/2=tau1_0, DL2'(i)=(-1+i)/2=tau2_0,
  DL1'(-i)=conj(tau1_0)=tau2_0) [rev_02].
- F (fold v->1/v/v^2) preserves top pole order with factor exactly -1, PROVED by the exact rational
  identity F(c/(v+i)^m) = c i^{-m} v^{m-2}/(v-i)^m and N(v)=i^{-m}v^{m-2}+1 with N(i)=0 (so the order
  stays <= m and the top coefficient is -c); verified m=1..80 [rev_02b].
- A/B exponential blocks: AB_part(v) = -(v^2+1)/(2 k! v (v-b)(v+a)) — denominator has only real poles,
  NO +-i poles ever (structural, all k; PF form verified k=1..6) [rev_02b].
- sqrt block initial data: denominator contains (v^2+1)^{2d-2} exactly (buildintegrand formula), so
  order <= 2d-2 at the base level; PF check d=2..9 all m [rev_02b].
- Induction: HigherPart(W_ij^{(k)}) for k>d is a sum of the three feeds; each feed has order <= 2d-2
  (Gm at degree d+1 has order <= 2d-1 by IH; DL zero lowers by 1; F preserves). Int lowers by 1, so
  Gm <= 2d-3. Descending-degree induction within each level; valid. The constant part of Gm
  (eval_b correction) times DL gives only simple real poles. Every term of the system is covered,
  including the F-branch (the folded A-feed from key (i,i) onto canonical key (i,i-1), coefficient i,
  using DL1) and the A/B blocks. Data: bounds hold on all 239 records k<=12 [rev_01].

## Review question 2 — Lemma B (top-order recursion): PROVED FOR ALL d (derived symbolically, not from data)
- (B1) [u^{-(2d-2)}](G*DL1) = tau1_0*[u^{-(2d-1)}]G, because DL1 = tau1_0 u + O(u^2) and, by Lemma A,
  G has no term above order 2d-1; the Cauchy product then leaves only the top term, and lower terms of
  G contribute nothing. Verified formally in exact fastred arithmetic, d=2..13 (also DL2/tau2_0, also
  "lower terms contribute nothing") [rev_03e].
- (B2) Int maps top at order 2d to top at order 2d-1 with factor -1/(2d-1): symbolic in d
  (int u^{-2d} du = -u^{1-2d}/(2d-1)) [rev_03]; exact fastred check d=2..13 [rev_03e].
- (B3) F-branch conjugation: F maps the -i top coefficient c at order m to +i top coefficient -c
  (rev_02b identity, all m); DL1'(-i) = tau2_0. Hence the adjacent-key feed is
  i*tau2_0*(h_ii - conj h_ii)/(2d-1) — exactly the formula in the code [rev_03, rev_03e].
- Assembled recursion re-verified against reducer data: k=3..12 ALL MATCH [rev_03].

## Review question 3 — Lemma C (rho closed form): PROVED FOR ALL d
rho_pq = K_d (-1)^p (1+[p>q]) C(d,p),  K_d = (-1)^d (2d-3)!! i^{d+1}/(2 sqrt2 d!).
- Finite-difference identity sum_s (-1)^s C(d,s) s^m (d-s)^{d-m} = d! (-1)^m: PROVED via
  sum_s (-1)^{d-s} C(d,s) s^n = d! S(n,d) (Stirling EGF (e^x-1)^d), expand (d-s)^{d-m}; only r=d-m
  survives, giving (-1)^{d-m} d!, times (-1)^d => (-1)^m d!. Verified exactly d=1..80, Stirling lemma
  d<=24,n<30 [rev_04].
- U_d^{(m)}(tau=(1+i)/2) = i^d (-1)^m C(d,m): PROVED by the generating function
  sum_m U_d^{(m)} x^m y^{d-m} = ((1-t)^2(e^{2x}-1)+t^2(e^{2y}-1))^d; at tau: tau^2=i/2, (1-tau)^2=-i/2,
  so the bracket is (i/2)(e^{2y}-e^{2x}); binomial expansion + the finite-difference identity.
  Verified exactly d=1..40 [rev_04].
- The analytic formula R~(v) = (1/(2 sqrt2)) sum_j binom(1/2,j) z^{2-2j} U_j(t)/(v t(1-t))
  (z=(v^2+1)/(2 sqrt2 v)) equals the v-converted sqrt part for ALL d: pure algebra
  (dt/dv=(v^2+1)/(4v^2), t(1-t) = -(v^2-a^2)(v^2-b^2)/(16v^2)) reduces it to buildintegrand's formula
  -(sqrt2/2) 8^d 4^{-D} v^{2d-1-D} Stilde (v^2+1)^{2-2d}/((v^2-a^2)(v^2-b^2)); also checked d<=6 in WT2b
  step_41. "Only j=d contributes to the top": z has a simple zero at i, so z^{2-2j} has pole order
  2j-2 < 2d-2 for j<d and U_j(t)/(v t(1-t)) is regular at v=i (t(i)(1-t(i))=1/2); verified d=2..15 [rev_11].
- Top coefficient: (1/(2 sqrt2)) binom(1/2,d) (z'(i))^{2-2d} * U_d(tau)/(i*1/2), z'(i)=sqrt2/2,
  (z'(i))^{2-2d}=2^{d-1}, 1/(i/... )=-2i.
- Fold doubling for p>q: folded top = -conj(top of (q,p) piece) (F factor -1 + reality), and
  (-i)^{d+1}(-1)^{q+1} = (-1)^{p+1} i^{d+1} since p+q=d: PROVED; verified d=2..26 [rev_07].
- Assembly: binom(1/2,d) = (-1)^{d-1} (2d-3)!!/(2^d d!) PROVED by the ratio recurrence
  (both ratios = (1/2-d)/(d+1), base d=1: 1/2) [rev_09, rev_11].
- Evidence: rho_pq == top coefficient of the ACTUAL level-d integrand, d=2..26 ALL MATCH [rev_07]
  (far beyond WT2b's d<=12 vs reducer data, also re-verified [rev_15 item 1]).
- Source audit: genintegrand == paper's mk_integrands.t_integrand, k=4..7, every monomial, exact
  [rev_14]. (Earlier failures in rev_12* were my comparison-routine bugs — symbol identity mismatch and
  a flaky splitter — not the claim.)

## Review question 4 — identities (I1)-(I4): PROVED FOR ALL d,p (WT2b had only checked d<=16, p<=10)
All four proved symbolically with sympy (d, p symbolic positive integers), and pushed numerically far:
- (I1) K_{d+1}(d+1) = -i(2d-1) K_d: symbolic, K_{d+1}/K_d = i(1-2d)/(d+1) [rev_09]; exact d=2..400 [rev_08].
- (I2) non-adjacent recursion (p>=q+2): assembled LHS-RHS = 0 symbolically [rev_10]; reduces to (I1),
  absorption (p+1)C(d+1,p+1)=(d+1)C(d,p), symmetry C(d,d-p)=C(d,p), i(tau1-tau2)=1 [rev_09];
  exact d=2..300 all q [rev_08].
- (I3) adjacent recursion: assembled LHS-RHS = 0 symbolically (after dividing (-1)^p) [rev_10]; reduces
  to (I1) at d=2p-1, C(2p,p)=2C(2p-1,p), C(2p,p+1)=(p/(p+1))C(2p,p), tau2-tau1=i, and rho_pp pure
  imaginary (K_{2p} = i(-1)^p (4p-3)!!/(2 sqrt2 (2p)!)) [rev_09, rev_10]; exact p=2..200 [rev_08].
- (I4) Im hat_p = Im rho_pp: symbolic, both sides = 2^{-2p-3/2} (4p)!/((4p-1) p!^2 (2p)!) [rev_10];
  exact p=1..200 [rev_08].

## Review question 5 — Lemma D logic + final step: PROVED, no hidden assumption
- Induction logic (level k ascending, degree d descending within each level) is valid: feeds come
  strictly from degree d+1; base d=k is the initial data. Adjacent key (p,p-1): at k=2p uses h_pp=rho_pp
  via (I3); for k>2p uses h_pp=hat_p, and (I4) (hat_p-conj hat_p = rho_pp-conj rho_pp) makes the same
  identity apply. Diagonal key (p,p) has only the A-feed, hence constant hat_p for k>2p.
  Re-verified against all data k<=12 [rev_15 item 2]; WT2b also checked k=13,14.
- No hidden assumption about lower-order principal parts: for d=2 keys the bound is 2d-2=2, so
  HigherPart IS the single top coefficient h u^{-2} + conj (verified: max order over k=2..12 is 2
  [rev_15 item 3]). Hence Gm_20, Gm_11 are determined exactly by rho_20=-sqrt2 i/4, rho_11=sqrt2 i/4,
  hat_1=sqrt2(1+i)/4: exact rational arithmetic gives Gm_20=g20, Gm_11^{(2)}=-g20, Gm_11^{(k>=3)}=g11
  (derived symbolically [rev_15 item 4]; data k=2..12 [item 5], k=13 spot check [rev_16]).
- The epsilon-identity (1-eps)(W_10 - eps R_10) = eps^2 Rw2 re-verified symbolically: c2=Rw2, c3=0,
  no other eps powers [rev_15 item 7]; no +-i poles survive (Rw2's denominator is 2v(v-b)(v+a)(v-a)
  [rev_16]); SimplePart(W_10^{(k)})=Rw2 for k=2..12 [item 6]. No hidden assumption: the F-branch enters
  W_10 with coefficient 1 as F(Gm_11 DL1), included in the identity.

## Evidence table (this review; all exact, residual 0)
| check | range | result | file |
|---|---|---|---|
| fixed-point system (*) vs reducer | 238 items, k<=12 | MATCH | rev_01 |
| Lemma A bounds on data | all 239 records | HOLD | rev_01 |
| DL zeros/tau values, F factor -1 identity | symbolic, m<=80 | PROVED | rev_02, rev_02b |
| A/B blocks: no +-i poles | structural + k<=6 | PROVED | rev_02b |
| sqrt block base order <= 2d-2 | structural + d<=9 | PROVED | rev_02b |
| (B1) product-top rule, (B2) Int factor, (B3) fold factor | symbolic d; exact d=2..13 | PROVED | rev_03, rev_03e |
| Lemma B recursion vs data | k=3..12 | MATCH | rev_03 |
| finite-difference identity | d=1..80 exact; proof via Stirling | PROVED | rev_04 |
| U_d(tau) evaluation | d=1..40 exact; proof via GF | PROVED | rev_04 |
| fold doubling of top | d=2..26 exact | MATCH | rev_07 |
| rho_pq vs actual level-d integrand | d=2..26 exact | MATCH | rev_07 |
| binom(1/2,d) closed form | symbolic recurrence + base | PROVED | rev_09, rev_11 |
| only j=d contributes to top | d=2..15 exact | MATCH | rev_11 |
| genintegrand vs paper t_integrand | k=4..7 all monomials | MATCH | rev_14 |
| (I1) | symbolic d; exact d<=400 | PROVED | rev_09, rev_08 |
| (I2) | symbolic d,p; exact d<=300 | PROVED | rev_10, rev_08 |
| (I3) | symbolic p; exact p<=200 | PROVED | rev_10, rev_08 |
| (I4) | symbolic p; exact p<=200 | PROVED | rev_10, rev_08 |
| Lemma D stabilization vs data | k<=12 (+13,14 WT2b; k=13 here) | AS PREDICTED | rev_15, rev_16 |
| Gm's from rho's (symbolic derivation) | exact | MATCH g20,-g20,g11 | rev_15 |
| SimplePart(W_10)=Rw2 | k=2..12 | MATCH | rev_15 |
| (1-eps)(W_10-eps R_10)=eps^2 Rw2 | symbolic identity | TRUE | rev_15 |
| no +-i poles in W_10(eps) | symbolic | NONE | rev_16 |

## One-paragraph status line for the paper
The WT2b proof of Conjecture 3 for all k>=2 is COMPLETE. Lemma A (pole-order bound) and Lemma B
(top-order recursion, including the F-branch conjugation formula) are proved for all k/d by symbolic
Laurent algebra around v=+-i; Lemma C's closed form rho_pq = K_d(-1)^p(1+[p>q]) binom(d,p) is proved
for all d from the master integrand via a generating-function evaluation at t=(1+i)/2 and the
finite-difference identity sum_s (-1)^s C(d,s) s^m (d-s)^{d-m} = d!(-1)^m (a Stirling-number
consequence); the closing identities (I1)-(I4), previously only checked for d<=16, are now PROVED for
general d/p (each reduces to the K_d ratio and elementary binomial identities; also verified exactly to
d<=400/p<=300). Lemma D's stabilization then follows by induction with no hidden assumptions, and the
exact epsilon-identity (1-eps)(W_10-eps R_10)=eps^2 Rw2 yields K_k=Rw2 for all k>=2. The only external
dependencies are the paper's master integrand and reduction algorithm (the setting, independently
re-verified here against mk_integrands for k<=7 and against the reducer data for k<=12,13).
'''
work = '<project>/research/08_theory/WT2r_review'
open(f'{work}/result.md','w').write(result)
import os
print("result.md written,", len(result), "chars")
print("rev files:", sorted(f for f in os.listdir(work) if f.startswith('rev_')))
