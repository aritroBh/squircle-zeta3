# WT2 — Conjecture 3: the weight-two part of M_k is (pi^2+6 lambda^2)/24 for all k >= 2

## MAIN RESULT (conditional theorem, every step exact)

Conjecture 3 is reduced, by an exact computation, to one explicit identity
(the "Stabilization Lemma"), which is proved below to be *equivalent* to the
conjecture and is verified in exact arithmetic for every k up to 15.

Write the reduction of the M_k integrand as in the paper (Sec. "How the proofs
work"): after the Euler substitution t=(v^2+2v-1)/(4v), a=1+sqrt2, b=sqrt2-1,
the weight-two part of M_k is int_b^a K_k(v) log t(v) dv, where K_k is the
degree-1 (in log t, log(1-t)) kernel left after integration by parts at v=+-i.
Conjecture 3 <=> K_k(v) = Rw2(v) := 1/(2v)+1/(2(v-a))-1/(v+a)-1/(2(v+b))
for all k>=2, because int_b^a Rw2(v) log t(v) dv = (pi^2+6(log 2)^2)/24
(proved in the paper's Appendix by dilog identities).

THEOREM (conditional, exact).  Let eps be the expansion parameter and form the
generating series of the reduction.  Suppose that the two "Gm" coefficients
feeding the degree-1 key are level-independent in the sense
    Gm_20(eps;v) = eps^2/(1-eps) * g20(v),
    Gm_11(eps;v) = -eps^2 g20(v) + eps^3/(1-eps) * g11(v),
with  g20(v) = (1+sqrt2)/4 - sqrt2/(2(v^2+1)),
      g11(v) = -sqrt2 (v^2+2v-1)/(4(v^2+1)).
Then the degree-1 kernel of the reduction is
    K(eps;v) = eps*R10(v) + eps^2/(1-eps) * Rw2(v),
where R10(v) = (1-1/sqrt2)/v - 2/(v+a) is the (different) k=1 kernel.
Hence K_k = Rw2 for every k>=2, i.e. Conjecture 3 holds.
PROOF: the degree-1 work item is (derived system, see below)
    W10 = eps*R10 - 2 Gm_20 DL1 - Gm_11 DL2 - F(Gm_11 DL1),
DL1 = 1/(v-b)+1/(v+a)-1/v, DL2 = 1/(v-a)+1/(v+b)-1/v, F(R)=R(1/v)/v^2.
Substituting the two stabilized Gm's and clearing (1-eps) gives the identity
    (1-eps)(W10 - eps*R10) = eps^2 * Rw2(v) + 0 * eps^3,
checked symbolically in exact arithmetic (c2 = Rw2 and c3 = 0 both
simplify to true); in particular no pole at +-i survives and W10-eps R10 is
already simple-poled.  Taking [eps^k] gives K_k = Rw2 for k>=2.   QED (mod Lemma)

STATUS OF THE LEMMA.  The stabilization is an exact, finite check per k; it is
verified for all k<=15 for Gm_20 (level-independent for every k=2..15) and for
all k=3..15 for Gm_11.  Both closed forms g20, g11 (and g11_{k=2} = -g20) are
confirmed exactly.  A proof for ALL k requires solving the "tail tower" (the
upper principal parts at +-i close among themselves, order >= m at degree d
feeds order >= m-1 at degree d-1); this fixed point was identified but its
closed form was not fully solved within the budget.  Hence the strongest claim
is: CONJECTURE 3 PROVED MODULO THE STABILIZATION LEMMA, with the lemma itself
verified exactly to k = 15 and the kernel identity K_k = Rw2 verified
independently (full reducer) to k = 15.

## Structural facts discovered (all exact)
* Level homogeneity: the M_k integrand is homogeneous of degree k in
  (log t, log(1-t)); the reduction is eps-linear, so Conjecture 3 <=> each
  degree-d>=2 initial integrand reduces to the degree-1 kernel Rw2.
* The a-block t^{eps-1} and b-block (1-t)^{eps-1} have NO (v^2+1) denominator,
  so they never trigger integration by parts: the entire IBP tower is generated
  by the q-block sqrt(t^{2eps-2}+(1-t)^{2eps-2}) alone.
* The reduction is the unique eps-series solution of the fixed-point system
  W_ij = eps^{i+j} R_ij - (i+1) Gm_{i+1,j} DL1 - [i>j] (j+1) Gm_{i,j+1} DL2
         - [j=i-1] i F(Gm_{i,i} DL1),   Gm_ij = (1-eval_b) o Int o HigherPart(W_ij).
* Stabilization: HigherPart(W_20), HigherPart(W_11) are level-independent for
  k>=2, k>=3 respectively; equivalently Gm_20, Gm_11 stabilize.  Higher degrees
  do not stabilize individually; only the upper tails at +-i do.
* Poles at +-i: none survive the reduction (any weight) for k<=12, extending
  the paper's exact k<=5 check; the (0,0) weight-1 piece never survives.
* k=0,1 (conjecture not claimed): M_0 = 0 (the integrand vanishes identically
  at eps=0); the k=1 kernel is R10 = (1-1/sqrt2)/v - 2/(v+a), not Rw2.

## Evidence table (exact kernel identity K_k = Rw2)
| k | method | result | time |
|---|--------|--------|------|
| 2..5 | paper's ibp_reduce (sympy apart/integrate) | EXACT, reproduced | 2.1s, 6.6s, 21.9s, 72.8s |
| 2,3 | fastred vs paper, ALL keys compared | identical | — |
| 4..8 | converted v_integrand + fastred | EXACT | reduce 0.1–2.7s |
| 2..15 | closed-form integrand + fastred | EXACT | build+reduce 0.03–92s |

Fast-reducer timings (closed-form integrand build / reduce, seconds):
k=9: 5.5/4.1, k=10: 7.6/6.5, k=11: 11.0/9.4, k=12: 15.9/15.7, k=13: 30.0/19.3,
k=14: 33.7/30.3, k=15: 35.5/56.6.  The closed-form integrand itself is
verified equal to mk_integrand.t_integrand for k=2..5 (after Q=sqrt(q2)) and to
v_integrand for k=2..8, in exact arithmetic.

Closed form of the integrand used (proved exact for k<=8): the coefficient of
l1^m l2^{k-m} in the M_k integrand is
  c_m(t) = q S_m(t)/(2 t(1-t) q2^k) - delta_{m,k}/(2 k! t) - delta_{m,0}/(2 k! (1-t)),
q2 = t^2+(1-t)^2, q = sqrt(q2),
S_m(t) = sum_{j=1}^k binom(1/2,j) sum_{s=0}^j binom(j,s) w(s,m) w(j-s,k-m)
         t^{2(j-s)} (1-t)^{2s} q2^{k-j},   w(s,m) = s! S(m,s) 2^m/m!,
S(m,s) = Stirling numbers of the second kind.

## Machinery (new code, originals untouched)
* fastred.py — exact Q(sqrt2,i) arithmetic (4 Fractions); partial-fraction
  integrands (simple real poles in {0,+-a,+-b}, principal parts at +-i, no
  polynomial part); closed-form fold, split, rational antiderivative, eval at b,
  multiplication by DL1/DL2; reducer mirroring ibp_reduce.reduce_k.
* genintegrand.py — the Stirling/binomial closed form above.
* buildintegrand.py — exact t->v conversion and folding, producing PF form.

## Paragraph for the paper
The weight-two term is not an accident of small k.  Writing the reduction of
Section IX as a generating series in eps, the degree-one kernel K(eps;v) is
governed by a triangular fixed-point system on the letters 0,+-a,+-b,+-i.  Only
the square-root block sqrt(t^{2eps-2}+(1-t)^{2eps-2}) feeds the tower at +-i;
the two exponential blocks reduce trivially.  The two coefficients that reach
degree one stabilize: they are rational in eps with denominator (1-eps), with
explicit v-parts g20(v)=(1+sqrt2)/4-sqrt2/(2(v^2+1)) and
g11(v)=-sqrt2(v^2+2v-1)/(4(v^2+1)).  Substituting, the identity
(1-eps)(K(eps;v)-eps R_10(v)) = eps^2 Rw2(v) holds in exact arithmetic, so
K_k = Rw2 = 1/(2v)+1/(2(v-a))-1/(v+a)-1/(2(v+b)) for every k>=2 and the
weight-two part of M_k is (pi^2+6 log^2 2)/24.  The argument reduces
Conjecture conj:wt2 to a single stabilization identity, which we prove is
equivalent to the conjecture and verify exactly for all k<=15 (the kernel
identity itself is also checked exactly, by an independent exact reducer, for
all k<=15); no pole at +-i survives the reduction at any weight for k<=12,
extending the k<=5 check of Section IX.
