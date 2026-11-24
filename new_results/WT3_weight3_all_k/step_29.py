
import os
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK)
result = r'''# WT3 — the weight-three part of M_k for every k

## VERDICT (current): PATTERN PROVED for k <= 15, all-k proof in progress

Main discovery: the degree-2 kernels are **affine in k** for all k >= 3, with an
explicit slope; equivalently the weight-three part W3(k) of M_k is affine in k
with slope DELTA = -zeta(3)/8 - pi^2*lambda/24 - lambda^3/12 (lambda = log 2).

**Conjecture WT3 (matched exactly for k = 3..15).** For every k >= 3,
  W3(k) = W3' - (k-3) * ( zeta(3)/8 + pi^2*lambda/24 + lambda^3/12 ),
where W3' is the weight-three part of M_3 (paper, thm:M3).  At the kernel level
the stronger statement holds: for every k >= 3,
  R_20(k;v) = (1-sqrt2)/(4v) + 1/(4(v-a)) + (2k-7)/(4(v+a)) - (k-2)/(2(v+b)),
  R_11(k;v) = (2k-5+sqrt2)/(4v) - (2k-5)/(2(v+a)),
with a = 1+sqrt2, b = sqrt2-1.  (k = 2 is exceptional, like k = 1 was for
weight two; its kernel is R_20(2) = (1/2-3*sqrt2/8)/v - 1/(v+a),
R_11(2) = (sqrt2/8)/v, giving exactly thm:M2's weight-three part.)

## Task 1: kernels (DONE).  Degree-2 kernels (after the v->1/v fold) computed
with the WT2 fast reducer + closed-form integrand for k = 3..15
(wt3_kernels.json).  CONTROL: at k = 3 and 4 the fast reducer agrees with the
slow sympy reference ibp_reduce.reduce_k EXACTLY, key by key, on ALL reduced
keys (not just degree 2).  Poles only at {0, a, -a, -b}; no pole survives at
+-i or at +b; the only endpoint pole is at v = a against log^2 t (as in the
paper's k=3,4).

## Task 2: values (DONE).  W3(k) evaluated exactly for k = 3..15 in the paper's
weight-3 Q(sqrt2)-basis (z3, Li3(b), Li3(-b), Li3(1/sqrt2), Li2(b), lam, L,
pi2; basis of m3_tails_exact.py: Li3 kept at b, -b, 1/sqrt2 only, Li2 reduced
by classical relations, Kummer relations included).  Controls: k=3 equals
thm:M3 W3' term by term; k=4 equals W3' - zeta3/8 - pi2 lam/24 - lam^3/12
(rem:k4) term by term.  Every k checked against direct mpmath quadrature of
the kernel integral at 70 dps: |exact - quadrature| < 3e-69 for all k
(w3_values.json, step_02.log).  Bonus control: the same pipeline at k=2
reproduces thm:M2's weight-three part EXACTLY (residual 3e-70).

## Task 3: pattern (DONE, exact for k<=15).
(a) kernels do NOT stabilize: FALSE (they drift linearly).
(b) kernels ARE polynomial in k: TRUE, degree 1, for k >= 3:
    R_20(k) = R_20(3) + (k-3) S_20,  S_20 = 1/(2(v+a)) - 1/(2(v+b));
    R_11(k) = R_11(3) + (k-3) S_11,  S_11 = 1/(2v) - 1/(v+a).
(c) eps-generating series rational: YES, denominator (1-eps)^2:
    W_20(eps;v) = eps^2 R_20(2;v) + eps^3 R_20(3;v)/(1-eps) + eps^4 S_20(v)/(1-eps)^2
    W_11(eps;v) = eps^2 R_11(2;v) + eps^3 R_11(3;v)/(1-eps) + eps^4 S_11(v)/(1-eps)^2
(first-difference constant for all k=3..15, exact; hence W3(k) affine with
slope DELTA = integral of the slope kernels = -zeta3/8 - pi2 lam/24 - lam^3/12,
independently confirmed value-by-value for k=3..15).

## Task 4: all-k proof (IN PROGRESS).  Route: WT2-style eps fixed-point
system; the slope must come from the drift of the degree-3 Gm terms
(Gm_30, Gm_21) feeding degree 2.  See proof_wt3.md.

## Files
- step_01.py / step_01.log: kernels k=3..15 + slow-reference control (PASS).
- step_02.py / step_02.log: exact values + controls + quadrature (PASS).
- wt3_kernels.json: exact degree-2 kernels (residues over Q(sqrt2)) k=3..15.
- w3_values.json: exact W3(k) in the basis, k=3..15.
- fastred.py, genintegrand.py, buildintegrand.py: WT2 machinery (unchanged).
'''
with open('result.md','w') as fh:
    fh.write(result)
print("result.md written")
with open('progress.log','a') as fh:
    fh.write("step_03: k=2 pipeline validated vs thm:M2 EXACT (residual 3e-70); k=2 exceptional, affine law is k>=3; result.md drafted\n")

# Now read WT2b stabilization results
WT2B = os.path.join(BASE, 'research/08_theory/WT2b_stabilization')
print(sorted(os.listdir(WT2B)))
