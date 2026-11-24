# WT3 — the weight-three part of M_k for every k

## VERDICT: REDUCED TO LEMMAS (verified k <= 15); the reduction theorem is PROVED for all k

**Conjecture WT3.** For every k >= 3 the weight-three part of M_k is
    W3(k) = W3' - (k-3) * ( zeta(3)/8 + pi^2*lambda/24 + lambda^3/12 ),   lambda = log 2,
with W3' = weight-3 part of M_3 (paper, thm:M3).  Stronger kernel form: for k >= 3,
    R_20(k;v) = (1-sqrt2)/(4v) + 1/(4(v-a)) + (2k-7)/(4(v+a)) - (k-2)/(2(v+b)),
    R_11(k;v) = (2k-5+sqrt2)/(4v) - (2k-5)/(2(v+a)),    a = 1+sqrt2, b = sqrt2-1.
Matched EXACTLY for k = 3..15 (kernels and values), with k = 3,4 controlled
key-by-key against the paper's slow reference ibp_reduce.reduce_k, and every
value checked against 70-digit quadrature.  (k = 2 is exceptional; its kernel
gives exactly thm:M2's weight-three part — also verified.)

## Task 1 (kernels): DONE.  Degree-2 kernels for k=3..15 via the WT2 fast
reducer + closed-form integrand (wt3_kernels.json).  Control: fast reducer ==
slow ibp_reduce.reduce_k EXACTLY on ALL keys at k=3,4.  Poles only at
{0, a, -a, -b}; no +-i poles survive; only endpoint pole at v=a vs log^2 t.

## Task 2 (values): DONE.  W3(k) exact for k=3..15 in the paper's weight-3
Q(sqrt2)-basis {zeta(3), Li3(b), Li3(-b), Li3(1/sqrt2), Li2(b), lam, L, pi2}
(m3_tails_exact basis: Li3 kept at b,-b,1/sqrt2; Li2 reduced by classical
relations; Kummer relations included).  Controls: k=3 == thm:M3 W3' term by
term; k=4 == rem:k4 (W3' - z3/8 - pi2 lam/24 - lam^3/12) term by term; k=2
pipeline == thm:M2 weight-3 part exactly.  Every k: |exact - quadrature| <
3e-69 at dps=70 (printed).

## Task 3 (pattern): DONE.  (a) stabilization: FALSE.  (b) kernels ARE
polynomial in k: TRUE, degree 1 for k >= 3 (first differences constant,
exact for k=3..15).  (c) eps-generating series rational: YES, denominator
(1-eps)^2:
    W_20(eps;v) = eps^2 R_20(2) + eps^3 R_20(3)/(1-eps) + eps^4 S_20/(1-eps)^2,
    W_11(eps;v) = eps^2 R_11(2) + eps^3 R_11(3)/(1-eps) + eps^4 S_11/(1-eps)^2,
    S_20 = 1/(2(v+a)) - 1/(2(v+b)),  S_11 = 1/(2v) - 1/(v+a).
Hence W3(k) affine with slope DELTA = -zeta3/8 - pi2 lam/24 - lam^3/12
(independently confirmed value-by-value, k=3..15).

## Task 4 (all-k proof): the reduction theorem is PROVED for all k; one lemma remains
PROVED for all k (proof_wt3.md):
 (P1) the fixed-point system of the reduction (WT2/WT2b; verified k<=15);
 (P2) its degree-2 specialization W_20 = -3 Gm_30 DL1 - Gm_21 DL2,
      W_11 = -2 Gm_21 DL1 (k>=3);
 (P3) the REDUCTION THEOREM: the affine Gm law
      Gm_30^(k) = Gm_30^(3) + (k-3) s30,  Gm_21^(k) = Gm_21^(3) + (k-3) s21,
      s30 = -(1/12)(1 - 2 sqrt2 v/(v^2+1)),  s21 = -3 s30,
      implies the kernel law for all k>=3 — by the four exact identities
      (C1)-(C4) (SimplePart computations, residual 0, no +-i poles survive),
      and hence W3(k) = W3' + (k-3) DELTA for all k.
REMAINING (verified k=3..15): the affine Gm law itself.  Its proof is reduced
to the tail tower: (T1) the order-m recursion (exact; DL(i)=0 gives the
2-shift); (T2) the polynomial-growth invariant nu_x(d,m)=floor((2d-4-m)/2)
(verified d<=8, L<=15, incl. the diagonal one-level-late caveat); (T3)
initial-data closed forms: t=0 is WT2b's Lemma C (PROVED all d), t=1 is
**PROVED here for all d**: rho_ij^(d)[2d-3] = -i(d-2) rho_ij^(d)[2d-2], via
(log A)'(i) = -i(d-2) which follows from t'(i)=0 in Stilde(v)=(4v)^D S(t(v));
t>=2 explicit but messy.  What remains for a complete proof: the zone-t
identities (stabilization at t=1, constant-drift at t=2) from these closed
forms, by WT2b's Lemma-D method.

## Files
- result.md (this file), proof_wt3.md (full proof status), progress.log
- wt3_kernels.json (exact degree-2 kernels k=3..15), w3_values.json (exact W3(k) k=3..15)
- verify_wt3.py (reruns tasks 1-3 + W1/W2 + C1-C4), step_01.py, step_02.py, step_04.py, step_05.py
- fastred.py, genintegrand.py, buildintegrand.py (WT2 machinery, unchanged)
