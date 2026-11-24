
import os
Z3DIR = 'Zeta' + chr(40) + '3' + chr(41)
BASE = os.path.join('<home>/Downloads/Research', Z3DIR)
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK)
proof = r'''# WT3 — proof document: the weight-three part of M_k for all k

Status: **REDUCED TO LEMMAS (every step either proved for all k or verified in exact
arithmetic for k <= 15).** The single remaining lemma is the **affine Gm law**;
its proof is reduced further to explicit tail-tower lemmas (polynomial growth +
initial-data closed forms + annihilation identities), each verified k <= 15.

## Notation
a = 1+sqrt2, b = sqrt2-1, L = log a, lam = log 2, DL1 = 1/(v-b)+1/(v+a)-1/v,
DL2 = 1/(v-a)+1/(v+b)-1/v, F(R)(v) = R(1/v)/v^2.  a=1+sqrt2, b=sqrt2-1.
The reduction of the M_k integrand (paper, Sec. IX; ibp_reduce.py) writes
M_k = sum_{i+j<=k} int_b^a R_ij(v) (log t)^i (log(1-t))^j dv with simple poles only;
the piece of log-degree d is a pure iterated integral of weight d+1.
The **weight-3 part** of M_k is the degree-2 piece:
    W3(k) = int_b^a [R_20(v) log^2 t + R_11(v) log t log(1-t)] dv,
where R_20, R_11 are the degree-2 kernels (after the v->1/v fold, keys (2,0),(1,1)).

## The main theorem (value form)
**Conjecture WT3.** For every k >= 3,
    W3(k) = W3' - (k-3) * ( zeta(3)/8 + pi^2 lam/24 + lam^3/12 ),
with W3' the weight-3 part of M_3 (paper thm:M3).
Kernel form (stronger): for every k >= 3,
    R_20(k;v) = (1-sqrt2)/(4v) + 1/(4(v-a)) + (2k-7)/(4(v+a)) - (k-2)/(2(v+b)),
    R_11(k;v) = (2k-5+sqrt2)/(4v) - (2k-5)/(2(v+a)).
(k=2 is exceptional; its kernels give exactly thm:M2's weight-3 part.)

## What is PROVED for all k

**(P1) The fixed-point system of the reduction (WT2, WT2b; verified k<=15).**
With eps the expansion parameter, the reduction is the unique eps-series solution of
the triangular system on canonical keys (i,j), d = i+j:
    W_ij = eps^d R_ij - (i+1) Gm_{i+1,j} DL1 - [i>j](j+1) Gm_{i,j+1} DL2 - [j=i-1] i F(Gm_ii DL1),
    Gm_ij = (1 - eval_b) o Int o HigherPart(W_ij).
Verified against the exact reducer for every key, k = 2..15 (step_05(a): ALL MATCH).

**(P2) The degree-2 specialization (verified k<=15).**
    W_20^(k) = delta_{k,2} R_20 - 3 Gm_30 DL1 - Gm_21 DL2,
    W_11^(k) = delta_{k,2} R_11 - 2 Gm_21 DL1.
(step_04: PASS for k = 2..15.)

**(P3) THE REDUCTION THEOREM (proved for all k).**
*Theorem.* Suppose the **affine Gm law** holds:
    Gm_30^(k) = Gm_30^(3) + (k-3) s30(v),   Gm_21^(k) = Gm_21^(3) + (k-3) s21(v)   for all k >= 3,
with s30(v) = -(1/12)(1 - 2 sqrt2 v/(v^2+1)),  s21(v) = -3 s30(v) = (1/4)(1 - 2 sqrt2 v/(v^2+1)).
Then the kernel form of Conjecture WT3 holds for all k >= 3, hence so does the value form.
*Proof.* By (P2), for k >= 3 the kernel is SimplePart(W_20^(k)), SimplePart(W_11^(k)),
and by the affine law
    SimplePart(W_20^(k)) = SimplePart(-3 Gm_30^(3) DL1 - Gm_21^(3) DL2)
                          + (k-3) SimplePart(-3 s30 DL1 - s21 DL2),
    SimplePart(W_11^(k)) = SimplePart(-2 Gm_21^(3) DL1) + (k-3) SimplePart(-2 s21 DL1).
The four exact partial-fraction identities (exact arithmetic over Q(sqrt2,i), residual 0):
    (C1) SimplePart(-3 Gm_30^(3) DL1 - Gm_21^(3) DL2) = (1-sqrt2)/(4v) + 1/(4(v-a)) - 1/(4(v+a)) - 1/(2(v+b)) = R_20(3),
    (C2) SimplePart(-3 s30 DL1 - s21 DL2) = 1/(2(v+a)) - 1/(2(v+b)) = S_20,
    (C3) SimplePart(-2 Gm_21^(3) DL1) = (1+sqrt2)/(4v) - 1/(2(v+a)) = R_11(3),
    (C4) SimplePart(-2 s21 DL1) = 1/(2v) - 1/(v+a) = S_11,
with **no pole at v = +-i surviving in any of the four** (checked: order-1 +-i
coefficients vanish identically).  Hence R_20(k) = R_20(3) + (k-3) S_20 and
R_11(k) = R_11(3) + (k-3) S_11 for all k >= 3, i.e. the kernel form.  The value
form follows by exact integration (paper's m3_tails_exact machinery; the slope
integral evaluates to -zeta(3)/8 - pi^2 lam/24 - lam^3/12; done once, exactly,
and rechecked term-by-term for every k = 3..15 in step_02).  QED.
(The inputs Gm_30^(3), Gm_21^(3) are the exact k=3 values, cross-checked against
the independent WT2b record and against the paper's ibp_reduce at k=3,4.)
Since (C1)-(C4) are exact rational-function identities verified in exact
partial-fraction arithmetic, the implication is proved for all k.

## The remaining lemma and its reduction (all verified k <= 15)

**LEMMA (affine Gm law).** For all k >= 3, with s30, s21 as above,
    Gm_30^(k) = Gm_30^(3) + (k-3) s30(v),  Gm_21^(k) = Gm_21^(3) + (k-3) s21(v).
Verified in exact arithmetic for k = 3..15 (step_04: PASS; Gm's cross-checked
vs the independent WT2b record at k = 3..12).

Reduction of the lemma to the tail tower (all verified k <= 15):

**(T1) Tail-tower recursion (verified k<=15).**  The +i principal-part
coefficients obey the order-m recursion: with h_ij^(k)[m] the order-m (m >= 2)
coefficient of HigherPart(W_ij^(k)) at v = +i, g[m'] = -h[m'+1]/m', and
tau1_n, tau2_n the Taylor coefficients of DL1, DL2 at v = i
(tau1_0 = (-1-i)/2, tau2_0 = (-1+i)/2),
    h_ij^(k)[m] = delta_{k,d} rho_ij[m]
       - (i+1) sum_{m'>=m} g_{i+1,j}[m'] tau1_{m'-m}
       - [i>j](j+1) sum_{m'>=m} g_{i,j+1}[m'] tau2_{m'-m}
       - [j=i-1] i * phi_m(Gm_ii DL1),
phi the v->1/v fold contribution (CX4-corrected diagonal term).  This is the
+i-projection of (P1), hence exact; the full system (P1) is verified k<=15.

**(T2) Polynomial-growth invariant (verified k<=15).**  For the level
differences X_ij^(L) = W_ij^(L) - W_ij^(L-1): for every key (i,j) with
d = i+j >= 2 and every order m in 2..2d-2, the order-m coefficient of
HigherPart(X_ij^(L)) at +i is
  - for i > j (off-diagonal): polynomial in L of degree <= floor((2d-4-m)/2)
    for all L >= d+1 (vanishing when the bound is negative);
  - for i = j (diagonal): the same for L >= 2p+2 (one level later, mirroring
    WT2b Lemma D(ii); at L = 2p+1 there is the exceptional "first-propagation"
    value).
Verified for all keys with d <= 8 and all orders, L <= 15 (step_05(c), plus
the diagonal tabulation).  The degree-3 instance: HigherPart(X_30^(L)),
HigherPart(X_21^(L)) have order 2 constant in L for L >= 4 and orders 3,4
zero — which, integrated, IS the affine Gm law (the constant order-2 value,
together with the constant orders 3,4 and the eval_b constant, gives s30, s21;
the values are fixed by the tower, cf. T3).

**(T3) Initial-data closed forms.**  The level-d initial integrand's +i
principal part, rho_ij^(d)[2d-2-t] (t = 0,1,...,2d-3), satisfies:
  t = 0 (WT2b Lemma C, PROVED for all d): rho_ij^(d)[2d-2] = K_d (-1)^i (1+[i>j]) binom(d,i),
      K_d = (-1)^d (2d-3)!! i^{d+1}/(2 sqrt2 d!).
  t = 1 (verified d = 2..15): rho_ij^(d)[2d-3] = -i(d-2) * rho_ij^(d)[2d-2].
  Higher t: rho[2d-2-t]/rho_top = beta_t is polynomial in d of degree t
      (verified d <= 12, t <= 6); e.g. beta_2 = -(2d^2-9d+11)/4 + i * (parity term).
Via the analytic-factor formula rho[2d-2-t] = A_ij^{(t)}(i)/t!,
A_ij(v) = -(sqrt2/2) 8^d 4^{-D} v^{2d-1-D} Stilde_ij(v) (v+i)^{2-2d}/((v^2-a^2)(v^2-b^2)),
the t = 1 identity is equivalent to the logarithmic-derivative identity
    (log A_ij)'(i) = -i(d-2)   <=>   Stilde_ij'(i)/Stilde_ij(i) = -i D_ij,
a single exact identity per (d, key) about the Stirling-number polynomial
Stilde — the natural next target for a WT2b-Lemma-C-style proof.

**(T4) Annihilation identities (verified k<=15 in data).**  At odd orders m,
the feeds in (T1) have polynomial degree one higher than the invariant (T2)
allows; their top-degree parts cancel exactly among the DL1-, DL2- and
fold-feeds (verified for all keys d <= 8, L <= 15).  Proving these identities
for all d from the closed forms (T3) is what remains for a complete proof.

## Evidence table (all exact arithmetic; residuals 0 or as stated)
| check | range | result | file |
|---|---|---|---|
| kernels R_20, R_11 via fast reducer | k=3..15 | exact, affine in k | step_01, wt3_kernels.json |
| control vs slow ibp_reduce, ALL keys | k=3,4 | MATCH key by key | step_01.log |
| W3(k) exact in the paper basis | k=3..15 | exact | step_02, w3_values.json |
| W3(k) vs 70-dps quadrature | k=3..15 | residual < 3e-69 (dps=70) | step_02.log |
| W3(3) == thm:M3 W3' term by term | k=3 | PASS | step_02.log |
| W3(4) == rem:k4 (W3' - z3/8 - pi2 lam/24 - lam^3/12) | k=4 | PASS | step_02.log |
| pipeline at k=2 == thm:M2 weight-3 part | k=2 | PASS (residual 3e-70) | step_03 |
| fixed-point system (P1) all keys | k=2..15 | ALL MATCH | step_05.log |
| degree-2 equations (P2) | k=2..15 | PASS | step_04.log |
| reduction identities (C1)-(C4) | all k (exact) | PASS, no +-i poles | step_C |
| affine Gm law (the Lemma) | k=3..15 | PASS | step_04.log |
| polynomial-growth invariant (T2) | d<=8, L<=15 | matches incl. diagonal caveat | step_05 |
| initial data t=1 = -i(d-2) rho_top | d=2..15 | ALL MATCH | step_52 |
'''
with open('proof_wt3.md','w') as fh:
    fh.write(proof)
print("proof_wt3.md written", len(proof))
with open('progress.log','a') as fh:
    fh.write("step_05: full fixed-point system ALL MATCH k<=15; poly-growth invariant confirmed (off-diag + diagonal caveat)\n")
    fh.write("step_52: initial-data t1 = -i(d-2) rho_top verified d=2..15; beta_2 closed form (real part) d<=12; proof_wt3.md drafted\n")
