result = r'''# WT2b — PROOF of the Stabilization Lemma for ALL k (hence Conjecture 3 for all k >= 2)

**STATUS: PROVED for all k (exact, with machine-verified components).**
Route used: (c)+(a) — an induction on the level k with an explicit invariant of the
"tail tower", namely the *top-order stabilization* (Lemma D below); the invariant
closes because the leading principal parts at v = +i obey a scalar recursion whose
initial data has a closed form satisfying three algebraic identities (I1)-(I3).
Everything is checked against WT2's exact data (all keys, k <= 12) and the
Stabilization Lemma itself against exact reducer output for k <= 14/15 (WT2).

## 0. Setup (from WT2, all facts re-verified in code)

- Euler variable: t = (v^2+2v-1)/(4v), a = 1+sqrt2, b = sqrt2-1, ab = 1, a-b = 2.
  z := q(v) = sqrt(t^2+(1-t)^2) = (v^2+1)/(2 sqrt2 v),  t(1/v) = 1 - t(v),
  t(i) = tau := (1+i)/2, tau(1-tau) = 1/2, tau^2 = i/2, (1-tau)^2 = -i/2,
  z'(i) = sqrt2/2,  dt/dv = (v^2+1)/(4v^2) = z/(sqrt2 v).
- The M_k integrand family is the eps-expansion of the master (genintegrand.py)
      Phi(eps) = A e^{eps l1} + B e^{eps l2} - Q' sqrt(1+u),
  Q' = q/(t(1-t)), u = ((1-t)^2 X + t^2 Y)/q2, X = e^{2 eps l1}-1, Y = e^{2 eps l2}-1,
  l1 = log t, l2 = log(1-t).  Expanding (binomial series + (e^x-1)^s = s! sum_m S(m,s) x^m/m!)
  gives the level-d integrand; only the sqrt block has poles at v = +-i.  The
  coefficient of l1^m l2^{d-m} (v-space, with measure dt/dv) is
      R~_m(v) = (1/(2 sqrt2)) sum_{j=1}^d binom(1/2,j) z^{2-2j} U_j^{(m)}(t(v)) / (v t(1-t)),
      U_j^{(m)}(t) = sum_s binom(j,s) w(s,m) w(j-s,d-m) (1-t)^{2s} t^{2(j-s)},
      w(s,m) = s! S(m,s) 2^m/m!.
  [VERIFIED == buildintegrand.sqrt_part_PF for d <= 6, all m.]
- Reduction = unique eps-series solution of the triangular fixed-point system on
  canonical keys (i,j), i >= j >= 0, i >= 1 (d := i+j):
      W_ij = eps^d R_ij - (i+1) Gm_{i+1,j} DL1 - [i>j](j+1) Gm_{i,j+1} DL2
             - [j=i-1] i F(Gm_{ii} DL1),         (*)
      Gm_ij = (1 - eval_b) o Int o HigherPart(W_ij),
  where DL1 = (log t)' = (v^2+1)/(v(v-b)(v+a)), DL2 = (v^2+1)/(v(v-a)(v+b)),
  F(R) = R(1/v)/v^2, HigherPart = principal parts of order >= 2 at v = +-i,
  R_ij = level-d initial data (zero for k != d).  [VERIFIED == reducer output for
  all 238 (k,key) pairs, k <= 12.]  DL1, DL2 have a SIMPLE ZERO at v = +-i:
  tau1_0 := lim_{v->i} DL1/(v-i) = (-1-i)/2,  tau2_0 := DL2/(v-i)|_i = (-1+i)/2 = conj(tau1_0).

## Lemma A (pole-order bound).  For all k and every canonical key (i,j), d=i+j:
  HigherPart(W_ij^{(k)}) has pole order <= 2d-2 at +-i, and Gm_ij^{(k)} has <= 2d-3.
PROOF.  Induction on k, top-down on d.  Initial data: R~_m has order <= 2d-2
  (worst term j=d gives z^{2-2d}; denominators v,t,1-t are analytic nonzero at i;
  the A/B exponential parts have no +-i poles).  Each feed reduces order by >= 2:
  Int lowers order by 1 and DL1,DL2 vanish to order 1 at +-i; F preserves order
  (F maps u=v-i to w = u/(1-iu) up to the factor -1, see Lemma B).  Hence
  order(W_ij) <= max(2d-2, (2(d+1)-3) - 1) = 2d-2; order(Gm_ij) <= 2d-3.  QED
  [VERIFIED on all data k <= 12.]

## Lemma B (top-order recursion).  Let u = v-i and, for d >= 2,
      h_ij^{(k)} := [u^{-(2d-2)}] HigherPart(W_ij^{(k)})|_{v=i}   (top +i coefficient).
  Then with rho_ij := h_ij^{(d)} (initial-data value; feeds vanish at k = d):
      h_ij^{(k)} = delta_{k,d} rho_ij + ( (i+1) tau1_0 h_{i+1,j}^{(k)}
                   + [i>j] (j+1) tau2_0 h_{i,j+1}^{(k)}
                   - [j=i-1] i tau2_0 conj(h_{ii}^{(k)}) ) / (2d-1),
  where for j = i-1 the middle term is i tau2_0 h_{ii}^{(k)}, so the last two
  combine to i tau2_0 (h_{ii}^{(k)} - conj(h_{ii}^{(k)})).
PROOF.  By Lemma A, Gm_{i+1,j} has top order exactly <= 2d-1 at +i; since DL1 has
  a simple zero at i, [u^{-(2d-2)}](Gm DL1) = [u^{-(2d-1)}]Gm . tau1_0 and no lower
  coefficient of Gm contributes; Int gives [u^{-(2d-1)}]Gm = -h_{i+1,j}/(2d-1).
  Same for DL2.  For the F-branch (j = i-1): writing w = 1/v + i = u/(1-iu) and
  1/v^2 = -1/(1-iu)^2 gives PP+(F(R)) = -sum_{r>=2} e_r u^{-r} (1-iu)^{r-2} + e_1/u,
  whose top coefficient is -e_s, e_s = top of R = Gm_ii DL1 at -i = conj(-h_ii
  tau1_0/(4i-3)) = -conj(h_ii) tau2_0/(4i-3) (reality of Gm_ii, DL1; 4i-3 = 2d-1).
  Collecting terms gives the recursion.  [VERIFIED against all data k <= 12.]  QED

## Lemma C (initial data, closed form).  For p >= q >= 0, p >= 2? (d = p+q >= 2):
      rho_pq = K_d (-1)^p (1 + [p>q]) binom(d,p),
      K_d = (-1)^d (2d-3)!! i^{d+1} / (2 sqrt2 d!).
PROOF.  Only the j = d term of R~_m contributes to order 2d-2:
      rho~(m) = (1/(2 sqrt2)) binom(1/2,d) (z'(i))^{2-2d} U_d^{(m)}(tau) /(i tau(1-tau)).
  Now (z'(i))^{2-2d} = (sqrt2/2)^{2-2d} = 2^{d-1}, 1/(i tau(1-tau)) = -2i, and
      U_d^{(m)}(tau) = sum_s binom(d,s) w(s,m) w(d-s,d-m) (-i/2)^s (i/2)^{d-s}
                     = (i/2)^d [x^m y^{d-m}] (e^{2y} - e^{2x})^d
                     = (i/2)^d . 2^d (-1)^m binom(d,m) = i^d (-1)^m binom(d,m),
  because [x^m y^{d-m}](e^{2y}-e^{2x})^d = (2^d/(m!(d-m)!)) sum_s (-1)^s binom(d,s)
  s^m (d-s)^{d-m} and the sum is (-1)^d d! times the leading coefficient
  (-1)^{d-m} of the degree-d polynomial s |-> s^m (d-s)^{d-m} (d-th finite
  difference), i.e. d! (-1)^m.  Assembling with binom(1/2,d) =
  (1/2)^d (-1)^{d-1} (2d-3)!!/d! gives rho~(m) = K_d (-1)^m binom(d,m).
  The folded piece (canonical keys i>=j fold the (q,p) term via t -> 1-t, which
  sends U to its t<->1-t image and yields the same value K_d (-1)^p binom(d,p);
  VERIFIED pf_fold(sqrt part) == sqrt part with S(1-t) for d <= 6): hence
  rho_pq = K_d (-1)^p binom(d,p) . (1 + [p>q]).  [Closed form VERIFIED against
  reducer initial data for ALL canonical keys, d = 2..12.]  QED

## Lemma D (top-order stabilization).  For every canonical key (p,q), d = p+q:
  (i)  off-diagonal p > q:  h_pq^{(k)} = rho_pq  for ALL k >= d;
  (ii) diagonal p = q:      h_pp^{(2p)} = rho_pp, and h_pp^{(k)} = hat_p for ALL
       k >= 2p+1, where hat_p := (p+1) tau1_0 rho_{p+1,p}/(4p-1).
PROOF.  Induction on the level k >= 2; within each level, downward induction on d.
  At k = d the feeds vanish and h_pq^{(d)} = rho_pq.  For k > d apply Lemma B;
  the feed keys have degree d+1 <= k, so by the induction hypothesis (already
  established at this level for degree d+1, and at all earlier levels) they take
  their stabilized values.  Constancy then reduces to the following identities,
  all proved from Lemma C (tau1_0 = (-1-i)/2, tau2_0-tau1_0 = i):
  (I1) K_{d+1} (d+1) = -i (2d-1) K_d.   [immediate from (2d-1)!!/(2d-3)!! = 2d-1]
  (I2) non-adjacent off-diagonal (p >= q+2):
       rho_pq (2d-1) = (p+1) tau1_0 rho_{p+1,q} + (q+1) tau2_0 rho_{p,q+1}.
       Proof: with rho = 2K(-1)^p binom, use (p+1)binom(d+1,p+1) = (d+1)binom(d,p),
       (q+1)binom(d+1,p) = (d+1)binom(d,q) = (d+1)binom(d,p); the RHS becomes
       2K_{d+1}(-1)^p (d+1) binom(d,p) (tau2_0-tau1_0) = 2K_{d+1}(-1)^p (d+1) i
       binom(d,p) = 2K_d (-1)^p (2d-1) binom(d,p) by (I1) = LHS.
  (I3) adjacent keys (p, p-1), p >= 2:
       rho_{p,p-1} (4p-3) = (p+1) tau1_0 rho_{p+1,p-1} + p tau2_0 (rho_pp - conj rho_pp).
       Proof: K_{2p} = (4p-3)!! i (-1)^p/(2 sqrt2 (2p)!) (pure imaginary) and
       K_{2p-1}(-1)^p = -(4p-5)!!/(2 sqrt2 (2p-1)!), so rho_pp = i(4p-3)!!/(2 sqrt2
       p! p!) and LHS = -2p (4p-3)!!/(2 sqrt2 p! p!).  The RHS assembles to
       (4p-3)!!/(2 sqrt2 p! p!) . p[(i-1)-(1+i)] = -2p (4p-3)!!/(2 sqrt2 p! p!).
  (I4) Im(hat_p) = Im(rho_pp): hat_p is (p+1)(-1-i)/2 . rho_{p+1,p}/(4p-1) with
       rho_{p+1,p} real, so Im(hat_p) = -(p+1) rho_{p+1,p}/(2(4p-1)); using
       K_{2p+1} = -i(4p-1)/(2p+1) K_{2p} (by (I1)) and (p+1)binom(2p+1,p+1) =
       (2p+1)binom(2p,p), this equals Im(rho_pp).  Hence
       hat_p - conj(hat_p) = rho_pp - conj(rho_pp), so (I3) covers both k = 2p
       (where h_pp^{(2p)} = rho_pp) and k >= 2p+1 (where h_pp^{(k)} = hat_p).
  With (I2),(I3) the off-diagonal values stay at rho forever (base cases k = d and
  k = d+1 respectively; propagation by the level induction), and (ii) follows for
  diagonal keys with no further identity.  [All identities VERIFIED in exact
  arithmetic for d <= 14; the whole statement VERIFIED against reducer data for
  ALL keys, k <= 12.]  QED

## THEOREM (Stabilization Lemma, all k).
      Gm_20(eps;v) = eps^2/(1-eps) g20(v),      g20 = (1+sqrt2)/4 - sqrt2/(2(v^2+1)),
      Gm_11(eps;v) = -eps^2 g20(v) + eps^3/(1-eps) g11(v),
      g11 = -sqrt2 (v^2+2v-1)/(4(v^2+1)).
PROOF.  By Lemma A the degree-2 keys (2,0) and (1,1) have +-i pole order <= 2, so
  their HigherPart is exactly the top coefficient (order-1 terms are excluded from
  HigherPart by definition):
  - (2,0), off-diagonal: h_20^{(k)} = rho_20 = -sqrt2 i/4 for all k >= 2 (Lemma D(i)).
  - (1,1), diagonal: h_11^{(2)} = rho_11 = sqrt2 i/4; h_11^{(k)} = hat_1 =
    2 tau1_0 rho_21/3 = sqrt2 (1+i)/4 for all k >= 3 (Lemma D(ii); uses h_21 = rho_21
    = -3 sqrt2/4 for k >= 3, i.e. (I3) at p=2).
  Hence HigherPart(W_20)(eps) = eps^2/(1-eps) (-sqrt2 i/(4u^2) + conj) and
  HigherPart(W_11)(eps) = eps^2 (sqrt2 i/(4u^2)+conj) + eps^3/(1-eps)(sqrt2(1+i)/(4u^2)+conj).
  Applying (1-eval_b) o Int (eps-linear, level-preserving):
      Int(+i part) of c u^{-2} is -c u^{-1}; with the conjugate part:
      rho_20:  sqrt2 i/(4(v-i)) - sqrt2 i/(4(v+i)) = -sqrt2/(2(v^2+1));
               value at b is -(1+sqrt2)/4, so Gm_20^{(k)} = g20 for all k >= 2.
      rho_11:  sqrt2/(2(v^2+1)); value at b is (1+sqrt2)/4, so Gm_11^{(2)} = -g20.
      hat_1 :  -sqrt2(1+i)/(4(v-i)) - sqrt2(1-i)/(4(v+i)) = -sqrt2(v-1)/(2(v^2+1));
               value at b is sqrt2/4, so Gm_11^{(k)} = g11 for all k >= 3.
  [All three antiderivative/eval identities computed exactly; the resulting Gm's
  match the exact reducer output for every k <= 12.]  QED

## COROLLARY (Conjecture 3, all k >= 2).
  The degree-1 work item is W_10 = eps R_10 - 2 Gm_20 DL1 - Gm_11 DL2 - F(Gm_11 DL1)
  (system (*) at key (1,0)).  Substituting the stabilized Gm's, the exact
  eps-series identity
      (1-eps)(W_10 - eps R_10) = eps^2 Rw2,  Rw2 = 1/(2v) + 1/(2(v-a)) - 1/(v+a) - 1/(2(v+b))
  holds in Q(sqrt2)(eps)(v) (verified symbolically: the left side is eps^2 c2 + eps^3 c3
  with c2 = Rw2, c3 = 0, no other terms; in particular no poles at +-i survive in W_10
  at any level).  Hence [eps^k] SimplePart(W_10) = Rw2 for all k >= 2 (and R_10 for k=1).
  By the paper's reduction theorem the degree-1 kernel of M_k is K_k = Rw2, and by the
  Appendix dilog identity int_b^a Rw2(v) log t(v) dv = (pi^2 + 6 (log 2)^2)/24.  QED

## Evidence table (all exact)
| component | statement | check | result |
|---|---|---|---|
| fixed-pt system (*) | == reducer, all W_ij | k<=12, 238 items | MATCH |
| Lemma A | pole orders <=2d-2 / <=2d-3 | all data k<=12 | HOLD |
| Lemma B | top-order recursion | all keys k=3..12 | MATCH |
| Lemma C | rho closed form | all keys d=2..12 | MATCH |
| R~ formula | == buildintegrand | d<=6, all m | MATCH |
| fold identity | pf_fold == S(1-t) build | d<=6 | MATCH |
| U_d(tau) formula | i^d (-1)^m binom(d,m) | d<=8 | HOLD |
| (I1)-(I4) | stabilization identities | exact, d<=14 / p<=7 | HOLD |
| Lemma D | top-order stabilization | all keys k<=12 | AS PREDICTED |
| Theorem | Gm_20=g20 (k>=2), Gm_11 | all k<=12 | MATCH |
| Corollary | (1-eps)Wcas = eps^2 Rw2 | symbolic eps, exact | TRUE |
| Stabilization | Gm_20 (k=2..14), Gm_11 (k=3..14) | WT2 exact reducer | MATCH |
| kernel identity | K_k = Rw2 | WT2 exact reducer, k<=15 | MATCH |

## Paragraph for the paper
We prove Conjecture 3 for all k.  In the eps-generating form of the reduction, the
integration-by-parts tower on the letters +-i closes on principal parts, and its
leading (top pole-order) coefficients at v = +i satisfy a scalar recursion
h_ij = delta rho_ij + ((i+1) tau_1 h_{i+1,j} + (j+1) tau_2 h_{i,j+1} [+ fold])/(2d-1)
with tau_{1,2} = (-1 +- i)/2, whose initial data has the closed form rho_pq =
K_d (-1)^p (1+[p>q]) binom(d,p), K_d = (-1)^d (2d-3)!! i^{d+1}/(2 sqrt2 d!).  Three
algebraic identities (a recurrence for K_d and two binomial evaluations, one of
which reflects the v -> 1/v fold of the diagonal keys) imply that the top
coefficients are level-independent: off-diagonal keys keep their initial value, and
diagonal keys stabilize one level later to an explicit twisted value.  Since the
keys (2,0) and (1,1) that feed the weight-two kernel have pole order <= 2, their
higher parts are exactly these stable top coefficients, which yields the
Stabilization Lemma Gm_20 = eps^2/(1-eps) g20, Gm_11 = -eps^2 g20 + eps^3/(1-eps) g11
with g20 = (1+sqrt2)/4 - sqrt2/(2(v^2+1)), g11 = -sqrt2(v^2+2v-1)/(4(v^2+1)).
Substituting into the degree-one work item gives (1-eps)(W_10 - eps R_10) = eps^2 Rw2
identically in eps, so the weight-two kernel of M_k equals Rw2 for every k >= 2 and
the weight-two part of M_k is (pi^2 + 6 log^2 2)/24.

## Files
- step_31..step_41: tail-tower data extraction, pattern discovery (polynomiality in k
  with degree d-1-ceil(r/2)), top-order recursion discovery, closed form for rho,
  verification of (I1)-(I4), stabilization vs reducer, final identity.
- fastred.py / genintegrand.py / buildintegrand.py: WT2 machinery (copied).
- record_W2b.pkl: exact W, HigherPart, Gm for all keys, k <= 12.
'''
with open('result.md','w') as f:
    f.write(result)
with open('progress.log','a') as f:
    f.write('\nWT2b: PROOF COMPLETE (all k). Top-order stabilization (Lemma D) via closed-form rho_pq + identities (I1)-(I4); Stabilization Lemma proved; Conjecture 3 proved for all k>=2. See result.md\n')
print("result.md written,", len(result), "chars")
