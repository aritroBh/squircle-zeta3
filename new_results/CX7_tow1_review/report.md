# CX7 independent referee report

**Overall verdict: PROOF COMPLETE, with the explicit write-up corrections in item 7.**

The operative T1 recursion and the all-degree symbolic certificates prove the depth<=2 tower law. At degree three this exhausts the higher part, yielding the affine Gm law for every k>=3. The exact C1–C4 identities and base/slope integrations then prove both claimed degree-two kernels and the weight-three value law for every k>=3. No unproved all-depth conjecture or finite-range extrapolation is needed. A fresh independent k=16 reduction agrees exactly.

Two submitted displayed statements are false: the combined B_d(delta) block double-counts the adjacent feed, and the Taylor-coefficient conjugation statement misses (-1)^n. The explicit corrections below agree with the existing implementation and certificates, so they leave no mathematical gap. This verdict applies to the proof with those corrections; it does not endorse the literal erroneous displays. Only checks actually run are described as fresh verification.

## Audit items

1. **OK, with wording corrections below:** tower law to affine Gm law to C1–C4 to the weight-three theorem. See the completed item 1 below.
2. **OK:** genuinely symbolic all-degree identities, matching the full T1 recursion (not the erroneous displayed B block).
3. **OK:** initial data through depth two and S(t*), S'(t*), re-derived and numerically checked.
4. **OK:** fold identity, E3, adjacent keys, p=1, and first diagonal levels.
5. **OK in the full recursion and certificates; ERROR in one displayed block:** corrected Lemma B is used in the actual argument; the B_d display requires the recorded indicator fix.
6. **OK:** fresh k=16 integrand and fast reduction give both claimed kernels exactly.
7. **OK (errata specified):** necessary corrections and scope/status updates are listed below; source files remain read-only.

## Execution record

- Created this report before detailed proof inspection or mathematical computation.
- All files outside this review folder remain read-only. Resumed computations use the requested interpreter directly with `-B`; the user's resumed instruction supersedes the original request for `nice`.

## Findings recorded during inspection

- TOW1's full T1 recursion retains the separate direct and folded diagonal feeds and is structurally correct. However, its displayed `B_d(delta)` block immediately afterwards repeats CX4's double-counting error: the second indicator must be `[p >= q+2]` when the final adjacent term already contains `p tau2_0 (u-conj(u))`. The full T1 implementation does not make this mistake. This distinction must be kept explicit when importing E0.
- The sentence `tau2_n = conj(tau1_n)` in TOW1 section 1 is false for odd n. The correct relation is `conj(tau1_n) = (-1)^n tau2_n`; `proof_common.py` implements this relation.
- The relevant depth range is `0 <= delta <= min(2,2d-4)`: degree two has only depth zero. The adjacent symbolic identities require p >= 2. At p=1 their child has degree one and no higher part; apparent denominators `4p-4` must not be evaluated there.
- Degree three has pole order at most four, hence exactly depths 0,1,2 exhaust its higher part. A proved tower law at these depths proves affine Gm functions for every k >= 3 by linear integration and normalization at b; finite checks are then only slope/base calibration.
- Read the previous CX4 review and CX5 corrected write-up directly from their supplied files. No prior-memory claim is being substituted for current inspection.

## 1. Logic chain and unrestricted k-range — OK

Carry forward the recorded pole-order argument: for degree three, the only higher orders are 4, 3, 2, precisely depths 0, 1, 2. Both relevant keys, (3,0) and (2,1), are off-diagonal, so their tower laws begin at k=d=3, including the initial level. Depths 0 and 1 are constant and depth 2 is affine. Complex conjugation supplies the -i part; principal-part integration and subtraction at b are linear. Thus the entire Gm functions are affine for every k>=3, with no unexamined higher-depth contribution and no diagonal delay at degree three.

The slopes need not be fitted: K_3=-sqrt(2)/8 gives rho_30(0)=sqrt(2)/4 and rho_21(0)=-3sqrt(2)/4. Hence e_30(2,1)=-sqrt(2)/12 and e_21(2,1)=sqrt(2)/4. Integrating the paired double poles and normalizing at b gives

    s30(v) = -1/12 + sqrt(2)*v/[6(v^2+1)],    s21(v) = -3*s30(v).

Here b/(1+b^2)=1/(2sqrt(2)). The base functions are obtained from the exact degree-three initial data. C1–C4 are four finite rational-function identities with these fixed inputs; their truth is independent of k, so substitution of the proved affine law proves the displayed kernels for all k>=3. At degree two the initial source occurs only at k=2 and is absent throughout this range. Linearity of the convergent degree-two integrals then gives the value law once the base value and the single slope integral are evaluated.

Freshly ran WT3's `claude_check_C1C4.py`: the four identities with the written inputs all returned True, as did its k=7 affine consistency checks. The deliberately retained diagnostic using Gm_21 from the wrong level returned False, confirming why the old verifier needed its correction. `check_value.py` reuses only the definition prefix of WT3's exact integration script, excluding its batch and all file-writing code. It freshly evaluates the single slope integral and the k=3 base integral, reduces them by the supplied exact polylogarithm identities, and obtains zero formal residual against, respectively,

    -zeta(3)/8 - pi^2*log(2)/24 - log(2)^3/12

and the paper's W3'. Independent 70-digit quadrature errors were 0 at the working precision for the slope and 3.62227e-71 for the base (`value.out`). Thus no further extrapolation from finite-k values is used to fix either constant.

There is no circularity: initial integrand -> triangular IBP system and pole bound -> initial tower data -> E0/E1/E2a/E2b/E3 -> tower law -> degree-three Gm -> C1–C4 -> kernels -> integral. At fixed arbitrary k, descending induction in degree suffices; the initial diagonal parent at k=2p is handled by E3, not by assuming its stabilized value. The finite k<=14/15 tests corroborate this argument, rather than supply its unrestricted premise.

## 2. Symbolic all-degree certificates — OK

Read all four `prove_*.py` scripts and `proof_common.py` and matched their terms to T1 and the Newton conversion. Freshly executed all four using the project Python with `-B`; `symbolic.out` records six exact zero remainders (E1 adjacent; E2b non-adjacent; E2b adjacent; E3 at each of depths 0,1,2). These are rational identities in unevaluated symbols d,p or p, not loops over numerical degrees.

For non-adjacent keys the normalized parent top coefficients are i(2d-1)/(p+1) and -i(2d-1)/(q+1), q=d-p. Their depth-one and depth-two coefficients are exactly the stated initial-data ratios. For adjacent keys d=2p-1, the two off-diagonal ancestors used to construct the stable diagonal give ratios i(2d-1)/(p+1) and (2d-1)(2d+1)/[p(p+1)]. The scripts retain both the direct diagonal feed and the separate fold, including the three fold orders at depth two. They correctly use ehat(2,0)-ehat(2,1) in the constant Newton coefficient. The conjugates are legitimate: the adjacent normalizing rho is real, whereas diagonal rho_pp is purely imaginary; E3's extra minus sign in `Dc` accounts for the latter. `tc1c` has the corrected odd-Taylor sign.

The non-adjacent E1 calculation in the text reduces to

    -i*(d-1)*(2d-1)/(2d-2) + 3*i/2 = -i*(d-2).

E2a is the E0 identity rescaled by the child/parent denominators 2d-3 and 2d-1; it is not an additional empirically assumed slope relation. Diagonal identities follow from the definition of ehat, because a diagonal has only its first off-diagonal parent. Thus the symbolic cases cover precisely the identities needed in the induction.

Domain restrictions matter: d>=2 and p>=q+2>=2 for non-adjacent keys; adjacent identities have p>=2; depth delta must satisfy delta<=min(2,2d-4). In these domains all denominators used are nonzero, and the binomial(4p-6,2) replacement by (4p-6)(4p-7)/2 is exact. Normalization divides only by nonzero rho(0). A zero rational remainder therefore establishes every admissible integer degree, not just those in the recorded data.

Execution note: the first `runpy` invocation lacked the sibling import path and exited before any identity was evaluated; after explicitly adding that read-only source directory, the run succeeded. No source files were changed.

## 3. Initial data and S-values — OK

Re-derivation: q2(t*)=0, so only j=d survives in S(t*). The Stirling support forces s=m in that sum, giving binom(1/2,d)*i^d*(-1)^m*binom(d,m). Differentiation leaves only j=d and j=d-1. The former contributes B_d'(t*)=B_d(t*)*2[(d-2m)-di]; the latter has s=m-1,m, with w(n-1,n)=(n-1)2^(n-1). Their sum gives the displayed S' formula, including the endpoint cases m=0,d when out-of-range binomials are zero. Dividing and cancelling gives the particularly simple form

    (S_m'/S_m)(t*) = -2*(d-2m)/(2d-3) - 2*d*i.

The factors (4v)^D cancel between Stilde and its prefactor. One can therefore use the analytic factor directly as

    A(v) = -(sqrt(2)/2)*8^d*v^(2d-1)*S_p(t(v))*(v+i)^(2-2d)
           /[(v^2-a^2)(v^2-b^2)].

Since t'(i)=0 and t''(i)=-i/2, its logarithmic derivatives are

    (log A)'(i) = -i*(d-2),
    (log A)''(i) = (3d-3)/2 - (i/2)*(S_p'/S_p)(t*).

The resulting rho(0) is Lemma C; rho(1)/rho(0)=-i(d-2), and equivalently

    rho(2)/rho(0) = -(2d^2-9d+11)/4 + i*(d-2p)/[2*(2d-3)].

The last two formulas are asserted only at admissible depths. The off-diagonal multiplicity is 2 because S_q(1-t)=S_p(t), while a diagonal is counted once.

`check_details.py` freshly verifies the logarithmic derivatives symbolically. It also constructs S directly from the finite Stirling sum (not from the closed-form value) and uses 80-digit numerical differentiation to check S, S', and the Taylor coefficients of A. Samples (d,m) are (2,2),(3,0),(3,2),(3,3),(4,2),(7,4),(11,8); rho is checked for the canonical samples at every admissible depth <=2. The largest absolute error divided by max(1,abs(expected)) was 2.2564e-81. These numerical checks corroborate the all-degree derivation; they are not its proof. See `details.out`.

The first local assertion compared an unsimplified SymPy expression structurally with zero; simplifying it corrected the check. The completed script exits successfully.

## 4. Fold, E3, and boundary cases — OK

The rational identity DL1(1/v)/v^2=-DL2(v) gives, for any diagonal Gm,

    -p*Gm(v)*DL2(v) - p*F(Gm*DL1)(v)
      = -p*DL2(v)*[Gm(v)-Gm(1/v)].

This is the sum of two distinct feeds. The normalized primitive's additive constant cancels in the bracket. `check_details.py` freshly verifies the DL identity by symbolic simplification. For E3, the relevant discrepancy is the initial diagonal tower minus the stable polynomial extrapolated one level back: Delta=rho-sum_s(-1)^s*ehat(.,s). It is not rho-ehat(.,0) at depth two. Since rho_pp(0) is imaginary, conjugation of a normalized discrepancy requires the minus sign used in `prove_E3.py`. Substitution in the complete T1 feed yields the three symbolic zeros in `symbolic.out`. Depth triangularity ensures that unexamined deeper discrepancies cannot enter these coefficients. E3 establishes blindness at the admissible depths <=2; no all-depth blindness is inferred.

The one exceptional level for the adjacent child (p,p-1) is k=2p. Its diagonal parent is then initial, while the child's Newton basis expects ehat extrapolated to k=2p. E3 equates exactly these feeds. At k>=2p+1 the parent follows its stable polynomial and the alternating Newton conversion applies directly. At k=2p-1 the child is itself initial, so no parent recursion is imposed.

Fresh exact boundary checks in `details.out`: p=2,3,4; depths 0,1,2; diagonal and adjacent tower values at k=2p,2p+1,2p+2; and equality of the initial/extrapolated diagonal feeds at k=2p. The first instance, (2,1) at k=4 fed by (2,2), is included.

At p=1 the child is (1,0), degree one, and has no higher part by the pole bound. There is no admissible child depth, so formulas containing 4p-4 must not be evaluated there. The diagonal (1,1) itself has only depth zero: rho_11(0)=sqrt(2)*i/4 at k=2, and ehat_11(0,0)=sqrt(2)*(1+i)/4 for k>=3. The local check verifies these at k=2,3,4, constructs the paired-pole primitive, and confirms that its combined direct/fold feed has no higher part at degree one. This handles the boundary without extending the symbolic formulas across singular denominators.

## 5. Corrected Lemma B — OK in the operative proof; displayed ERROR to fix

Carry forward the previous inspection rather than repeat CX4's counterexample computation. The full T1 formula and `Lambda_ij`/`LambdaE_ij` retain the direct DL2 feed with indicator [p>q] and a separate fold. At fixed depth the combined diagonal contribution is p*tau2_0*(u_pp-conj(u_pp))/(2d-delta-1). Consequently the already-combined block must read

    B_d(delta)(u) = [(p+1)*tau1_0*u_(p+1,q)
                    + [p>=q+2]*(q+1)*tau2_0*u_(p,q+1)
                    + [q=p-1]*p*tau2_0*(u_pp-conj(u_pp))]
                   /(2d-delta-1).

At delta=0 this is exactly CX5's corrected Lemma B. Keeping [p>q] in this combined display adds the direct diagonal term a second time. The scripts audited in item 2 implement the separate-feed version, which is correct, and their symbolic adjacent identities therefore certify the corrected block.

E0 is available as a proved result: CX4/CX5 give general binomial proofs of I1–I4 and descending induction at arbitrary fixed k, including preservation of the imaginary part of the initial/stable diagonal coefficient. The E2a rescaling imports this corrected E0, so it has no dependency on the false combined display. The printed error is substantive as an equation, but is repaired by the single indicator replacement and does not require changing the implementation or finding an additional identity.

## 6. Independent point beyond the recorded k-range — OK

Ran `check_k16.py` with the specified project Python, directly and with `-B`. It imports `buildintegrand.py` and `fastred.py` from WT2, freshly builds the level-16 integrand, and runs the full fast reducer. It does not load saved integrands, tower data, kernel JSON, or candidate tower values. The expected kernels are constructed separately from the formulas in the brief.

The exact results were

    R20(16;v) = (1-sqrt(2))/(4v) + 1/[4(v-a)] + 25/[4(v+a)] - 7/(v+b),
    R11(16;v) = (27+sqrt(2))/(4v) - 27/[2(v+a)].

Both differences from the claimed formulas are identically zero in exact Q(sqrt(2),i) partial-fraction arithmetic; neither output has any pole at +/-i or a polynomial part. Fresh integrand construction took 77.36 s and the whole run 197.15 s on this machine. `k16.out` records the successful result; `check_k16.py` is the runnable independent check. This supplies the requested new point beyond k<=15, while the proof of unrestricted k is the symbolic argument in items 1–5.

## 7. Required write-up corrections — specified

1. **TOW1 §2, B_d(delta) display:** change the second indicator from [p>q] to [p>=q+2] when retaining the already-combined adjacent term. The full T1 formula's separate direct term keeps [p>q]. This repairs a false displayed equation, not merely its prose.
2. **TOW1 §1, conjugation:** replace tau2_n=conj(tau1_n) by conj(tau1_n)=(-1)^n*tau2_n. For example tau2_1=-conj(tau1_1). The implementation already has the correct relation.
3. **Domains:** state d>=2 and 0<=delta<=min(2,2d-4) throughout the depth-two theorem; degree two has only depth zero. Adjacent E1/E2b/E3 statements about these child depths require p>=2. Handle p=1 separately by the degree-one pole bound, without evaluating the formal denominators.
4. **TOW1 §4 hypotheses:** explicitly list E2a, the depth-two slope equation, together with E0/E1/E2b/E3, and immediately point out its proved reduction to corrected E0. Specify that diagonal entries in the fixed-point family mean stable ehat values; E3 handles their exceptional initial levels.
5. **TOW1 §7 and result table:** replace the wording that the slopes are only verified through k<=14 with their closed-form derivation from e(2,1)=-rho(0)/(2d-3), as in item 1. Explain that orders 2,3,4 exhaust the degree-three higher part, so the affine law starts at k=3 for both keys.
6. **Status text:** TOW1's opening finite-range status is stale: E1/E2b/E3 now have all-degree symbolic certificates. WT3's opening status, “Conjecture WT3,” and “remaining lemma” passages should be updated to cite the completed depth<=2 tower proof. P1/P2 follow from the exact triangular IBP system; finite reducer comparisons are validation, not the reason those equations hold for every k.
7. **Keep verification and proof distinct:** list the symbolic rational certificates, admissible domains, and their zero remainders separately from V1–V8's finite checks. C1 must use both Gm functions at level 3. Do not reintroduce any of the corrected verifier's old control-flow/input mistakes.
8. **Clean the initial-data derivation:** replace S3's unfinished cancellation line by B_(d-1)(t*)/B_d(t*)=-i*(d-1)*(d-2m)/d; state t''(i)=-i/2 without an ellipsis. In `towlib.r_i`, the docstring should say S'/S, not S''/S; its executable formula is correct.
9. **All-depth outlook:** retain it as conjectural. The title of TOW1's result must not imply an established polynomial law at every weight. The speculative factor 3 at depth four in §7 conflicts with the later observed factor 2 in §7b; reconcile that text. Also replace the all-depth assertion that S_m^(t)(t*)/S_m(t*) is polynomial in (d,m) by the appropriately qualified rational-function statement: already S'/S has denominator 2d-3. None of these all-depth claims is needed for the theorem reviewed here.

All corrections above are recommendations for the coordinating write-up. Per the brief, no file outside CX7 was edited.

## Verification provenance and reproducibility

- Fresh successful checks in this resumed audit: all four symbolic scripts (`symbolic.out`); `check_details.py` (`details.out`); `check_k16.py` (`k16.out`); `check_value.py` (`value.out`); and the existing C1–C4 diagnostic described in item 1.
- The local `verify_tower.out` currently contains only V1–V3 and the earlier `nice` warning. The full 12/12 PASS output is present in `../TOW1_tower_law/claude_rerun_verify.out`, which I inspected. I do not describe the truncated local file as a fresh complete run. The WT3 coordinated 10/10 PASS output was also inspected at `../WT3_weight3_all_k/claude_rerun_verify4.out`; that is prior evidence, distinct from the fresh checks above.
- Re-run local checks using `<project>/.venv/bin/python -B` followed by the local script name (quote the executable path in a shell). They require only the project's installed packages and read-only source/data directories. No new dependency was installed.
