# REFEREE_CODEX — G9 adversarial review

2026-11-23. Supplied files have not been edited. All code, logs and scratch products from this review are new files in this directory. A rerun of supplied code is identified as such; it is not an independent derivation. The fresh S1–S6 chain, exact coefficient comparison and r=2 fold have completed successfully.

## Verdicts

| Result or sub-claim | Verdict | Qualification |
|---|---|---|
| Proposition 9: the two local coefficient laws | **PROVED** | Computer-assisted; independently recomputed exactly, with both differences zero. |
| Proposition 9: Beta rules, fully local case, parallelogram and regular-polygon substitutions | **PROVED** | Hand arguments are correct; independent algebra agrees. |
| Diamond, regular hexagon and sheared-parallelogram displayed forms | **PROVED** | Supplied exact local computation rerun successfully. Regular hexagon also independently checked from the polar definition. |
| Theorem 26: exact HPL representation and continuation | **PROVED AFTER REPAIRS** | Exact symbolic construction, no numerical fitting in S1–S6. Correct the normalization error and auxiliary statements below. Fresh exact rerun and complete table comparison pass. |
| Theorem 26: explicit zeta(5) and zeta(2)zeta(3) coefficients are local | **PROVED AFTER REPAIRS** | True for the specified formal expression; proof has a factor-120 error. No-new-constants argument survives scrutiny of the actual summation algorithm. |
| Proposition 10: nonzero antisymmetrised (3,2) symbol component and its pair | **PROVED AFTER REPAIRS** | Independent exact symbol calculation reproduces the pair. Define the weight projection, quotient and regularization precisely. |
| Proposition 10: no classical Li5-plus-products symbol | **PROVED AFTER REPAIRS** | S-unit case is proved in the text. An elementary extension argument supplied below proves the arbitrary-argument version over a specified field. |
| Proposition 10: non-proportionality to the quoted M4 pair | **PROVED** | Elementary exact calculation, independently repeated. The M4 calibration itself was not independently recomputed here. |
| An unconditional nonclassicality assertion about the numerical period | **NOT PROVED** | The manuscript explicitly makes this conditional on a period conjecture. |
| Independence/minimality of the 307 functions or 20 numerical constants | **NOT PROVED** | Neither is established by a symbolic expression or high-precision agreement. |
| Full classical-plus-C symbol identity and complete octagon closed form | **NOT PROVED in this review** | Matching two pairings is not a proof of the full identity. The complete numerical closed forms remain numerical. |

No fatal numerical or algebraic counterexample has been found. The substantive repairs concern what the octagon object means and what quotient is being tested; the hexagon coefficient proof has a genuine but readily repaired normalization mistake.

## Numbered findings

1. **Must fix — factor 120 in Lemma z5's proof.** `sec_g9_hex5.tex:127,131`. The lemma states that the coefficients in the unmultiplied sum are Z/120. Its proof instead calls `5!*24*integral` the coefficient in that sum and concludes it is Z. The unmultiplied coefficient is `24*integral`; multiplication by 120 gives the edge contribution. Delete the factor 5! there and conclude Z/120, or explicitly extract the coefficient of `E_e^(5)=120*sum`. This does not invalidate the stated coefficient law or the table convention.

2. **Must fix — W5 is defined twice, with different meanings.** `sec_g9_oct5.tex:33-51,55,58`; printed excerpt `main_tex_statements_l1715-1778.tex:51`. At :37 W5 is the entire integral of the degree-five log part; at :47–51 it is the projection retaining only terms of iterated-integral weight five. Those are different: integration by parts also produces lower weights. Define

   `W5 = gr_5(960*(3-2*sqrt(2))*integral(E5_top/L))`,

   with constant logarithms assigned weight one. The discarded terms have weight at most four. Say **“The symbol of W5 has nonzero antisymmetrised deconcatenation component of bidegree (3,2), modulo shuffle products.”** The pair detects this associated-graded object; it is not an unconditional statement about the number P^(5)(1).

3. **Must fix — divergent endpoint evaluations in the integration-by-parts lemma.** `sec_g9_oct5.tex:40-44`. For words starting with 1, G(w;1) can diverge, so `0*G(w;1)` is not a literal value. Use the limit of `(t^(n+1)-1)*G(w;t)` as t tends to 1, which is zero for logarithmic divergences, or specify tangential regularization. The integration-by-parts identity is valid with that convention. Its remainder has weight **at most** k−1, not necessarily exactly k−1. The lower-endpoint argument is otherwise sound.

4. **Must fix — coincident letters require a stated symbol convention.** `sec_g9_oct5.tex:22-26`; implementation `scripts/g9_oct5_symbol.py:111-121`. The displayed recursion uses a ratio that can have zero numerator or denominator. “Entries equal to 0 omitted” does not specify how to treat those ratios. The code actually uses `[a_i-a_(i+1)] - [a_i-a_(i-1)]` and drops each zero difference separately. State this split-entry rule and the shuffle/tangential regularization realizing it. For G(p,p;1), simply omitting a whole zero-numerator ratio would lose its negative log(p) contribution. My recomputation uses the implemented convention and finds the printed pair; this is a definition repair, not a detected numerical failure.

5. **Must fix — finish the arbitrary-argument classical-symbol argument.** `sec_g9_oct5.tex:63,106`. The hand lemma proves vanishing only when both entries lie in V=QX+QY; the text invokes references for the general case without identifying the pairing with the cited cobracket. A shorter repair avoids that dependency. Set U=F* tensor Q, torsion removed, for the specified field F, and extend the two coordinate functionals of V to U. Project tensor entries to V and apply the same pairing. Every classical Li5 symbol remains `-a tensor b tensor b tensor b tensor b` after projection and pairs to zero by the existing hand calculation. Shuffle products remain shuffle products and also pair to zero. A putative equality with W5 over F is therefore impossible. Specify algebraic scalar coefficients, the field of arguments, and branches modulo torsion. The same argument works in a specified larger field containing F. No converse to a depth conjecture is needed for this negative result.

6. **Minor — harmonic-number range is wrong in one general assertion.** `sec_g9_hex5.tex:74`. The local moments for a+b<=5 can involve H_(N+1)^(5), as the (0,5) moment does. S3 includes H5 correctly. Change i<=4 to i<=5 for the whole local table, and retain i<=4 for the non-local summands, where a+b<=4. The no-new-zeta argument is unaffected.

7. **Minor — repair two hand/machine labels and the continuation wording.** `sec_g9_hex5.tex:27,84,118`. HPL generation in (d) uses an exact matrix inversion, so it is hand framework plus machine step S2, not wholly “by hand.” The continuation's right side is initially analytic on the **punctured** slit plane; equality there with the analytic integral establishes removability at zero. The current proof calls it analytic on the unpunctured domain before proving pole cancellation. The identity-theorem argument itself is sound once phrased in that order.

8. **Must fix if intended unconditionally — occurrence is not numerical minimality.** Printed excerpt `:35,48`; supplement `sec_g9_hex5.tex:156,177`. A nonzero coefficient of Li4(1/2), Li5(1/2), H00101 or H01011 in a chosen representation does not prove that no other representation eliminates it. The manuscript does acknowledge numerical independence assumptions and lack of minimality, but “what does not depend on [the choice] is that ... occur” and “genuine multiple polylogarithms” should not be read as an unconditional necessity theorem. Replace with: **“These functions occur in the displayed representation. Under the stated independence assumption, the value lies outside the span of the specified local constants.”** For the generic HPL theorem, say only that its explicit expression uses the thirteen listed functions; no minimal generating-set claim is proved here.

9. **Minor — Table II should name the symbol object.** Live paper `paper/ieee/main.tex:2175` labels the row “octagon P^(5)(1) needs depth two,” though its status column says symbol level. Replace the row label with **“Octagon: nonzero (3,2) component of the weight-five symbol.”** This makes the status unambiguous without requiring the reader to recover the conditional period statement.

10. **Minor — the r=2 certificate justification overstates absolute convergence.** `sec_g9_hex5.tex:175` says the double-shuffle relations follow from absolute convergence. The external relation generator includes the one-letter word (-1,), whose Euler-sum series is alternating harmonic and only conditionally convergent. The relations can instead be justified by finite-cutoff stuffle identities followed by convergence, together with absolutely convergent iterated-integral shuffle identities. Give that justification for the actual DS family. This does not require numerical fitting, and no invalid reduction was detected.

## Hand arguments checked

- **Beta moments (`sec_g9_order5.tex:44-63`): correct.** In log Beta, the z5 term is `(epsilon^4 delta + 2 epsilon^3 delta^2 + 2 epsilon^2 delta^3 + epsilon delta^4)*z5`. The z2 and z3 terms are `-epsilon*delta*z2` and `(epsilon^2*delta+epsilon*delta^2)*z3`. Their product gives the stated z2*z3 coefficient. Multiplying by a!b!B gives 24B and −12B in exactly the stated cases. Weight counting excludes all other cases.
- **Local extraction (`:69-105`): correct.** A weight-five zeta monomial forces full log degree five and no corner-log shifts. The leading-log two-variable series is the correct p=1 radial length element. The finite polynomial identity has been independently checked below.
- **Parallelogram and regular-polygon substitutions (`:116-134`): correct.** Signed L0,L1 and the substitution u0=−1/2 have been handled consistently.
- **Shape (`sec_g9_hex5.tex:47-52`): correct.** Phi starts at order d and has no L denominator; each square-root term contributes L^(1−2j), j<=5, times a polynomial in L² and the other edge data. Hence denominator L^9 and only odd L powers. This proves an allowed denominator, not a claim that every geometry attains the worst pole.
- **HPL generation (`:76-84`): correct for the defined words ending in 1.** The shuffle algebra closes on these words. My independent matrices are invertible through weight five; no independence of evaluated functions is needed for this spanning assertion.
- **No new constants (`:97-103`): correct for the actual S4 reduction.** S4 carries z2,z3 as formal coefficients; its partial fractions, shifts, finite initial sums, atom recognition and rational shuffle inversion never introduce a zeta symbol. The source Beta moment has a+b<=4, so its zeta monomial weight is <=4. This proves absence of explicit z5 and z2*z3 in this particular formal J representation. It does **not** imply their absence after special-value reduction of the HPLs.
- **Local coefficients for every r (`:126-134`): correct after finding 1.** Constant corner logarithms reduce the remaining log degree and cannot supply a new explicit weight-five zeta monomial. Only c=0, a+b=5 contributes.
- **Continuation (`:107-120`): correct after finding 7.** Holding the corner-log symbols independent while varying kappa is legitimate for this finite polynomial/integral identity. Equality on the small real interval extends on the punctured slit plane; the integral gives the removable value at zero. No unjustified power-series convergence at r>=2 is assumed.
- **Octagon integration by parts (`sec_g9_oct5.tex:39-44`): correct with findings 2–3.** Polynomial division gives the printed leading coefficient. Lower-weight terms cannot contribute to the length-five tensor being paired.
- **Li5 and product controls (`:67-74`): correct.** The repeated tensor entries cancel against the degree-two or degree-three Lie polynomial. Lie elements annihilate proper shuffles, and the two surviving 3+2 block orders cancel under antisymmetrization. The extension in finding 5 completes the field issue.
- **Non-proportionality (`:108`): correct.** The two candidate multipliers have distinct rational parts; I printed the exact difference below.

## Precise replacement for Proposition 10's scope

Let F=Q(sqrt(2)), V be the rational span of the multiplicative classes of 2 and 1+sqrt(2), with torsion removed, and W5 be the associated weight-five component defined in finding 2. In the shuffle Hopf algebra on V, take indecomposables modulo proper shuffle products, then the antisymmetrised reduced deconcatenation of bidegree (3,2). The symbol of W5 has the two printed nonzero coordinates against the specified Lie elements. Therefore it is not the symbol of an F-linear combination of classical Li5 values and lower-weight products with arguments in F. For the numerical period itself, a corresponding motivic/period-injectivity assumption is required.

This states exactly what the computation obstructs. It does not assert a unique depth-two formula or prove the full classical-plus-C decomposition. The extension argument in finding 5 permits a specified larger field of arguments as well.

## Exactness audit of S1–S6

**The J table is derived, not fitted and then checked. No PSLQ, LLL, numerical linear solve, rational reconstruction from decimal data, or coefficient fitting occurs in S1–S6.** I read the complete summation engine and the construction portions of all six scripts.

| Step | Actual construction and ring | Numerical role |
|---|---|---|
| S1, `g9_hex5_s1_engine.py:17-77` | Truncated d-series through degree five. F, D, inversion, logarithm, exponential and square root use rational coefficients. Formal edge data and log variables; coefficients in Q[A0,A1,B,kappa,q0,V02,lA0,lA1,lB,L²,(L²)^−1] times sqrt(L²), subject to the edge relations. | None in construction. External E3 pickle is only compared afterward. |
| S2, `g9_hex5_s2_atoms.py:19-186` | Enumerates Lyndon words, forms and inverts rational shuffle matrices; ordered set partitions expand products of finite harmonic sums, followed by H_n=H_(n−1)+1/n^j. Produces 26 atom expressions in 13 formal HPL generators. | Numerical checks are separate. |
| S3, `g9_hex5_s3_jloc.py:9-33` | Formal log-Gamma/Beta expansion through degree five, coefficients in Q[DEN,H1,...,H5,z2,z3,z5], z4=2*z2²/5. Stores 21 pairs (a,b). | Quadrature begins only after the exact table is saved. |
| S4, `g9_hex5_s4_jtab.py:18-46`; `g9_esum5.py:26-107` | Builds the exact log^c series coefficient, substitutes the symbolic Beta moment at N=n+q, aligns H_(q+d) with H_q, applies rational partial fractions, shifts denominators and subtracts finite initial terms. Recognizes sums x^q H_q^e/q^a as S2 atoms. Output lies in Q[x,x^−1,HPL13,z2,z3]. | All 210 expressions are saved before the 40-key/four-kappa quadrature block at :49-109. Those 160 checks are not derivations or fits. |
| S5, `g9_hex5_s5_blocks.py:7-34` | Reads S1/S3/S4, substitutes exact integer harmonic values for local terms, packages exact building blocks. | None. |
| S6, `g9_hex5_s6_assembly.py:7-43` | Multiplies every E5 coefficient by its exact local or non-local integral, decomposes formal HPL/zeta/x monomials and accumulates edge-data coefficients. | None. |

The chain's mathematical credibility depends on these symbolic algorithms and their implementation, not on decimal agreement. It is a computer-assisted proof, not a hand derivation of every coefficient. I independently checked M1, the shuffle matrices, one whole table evaluation and the octagon symbol, but did not independently reimplement all 210 general J reductions or the symbolic S1/S6 engine.

PSLQ in the research history is separate: regular hexagon discovery, the r=2 numerical relation, and initial discovery of x=−1 special values. The generic theorem uses none of those special values. The r=2 corollary additionally requires exact special-value certificates; a numerical check of those values alone would not suffice.

## Completed fresh-chain results

All six steps completed. Their reported per-step times sum to 901 seconds (`referee_exact_types.out`). S4 saved its 210 expressions after a reported 563 seconds; including its numerical check block, that shell step took 628 seconds. Its 160 quadrature comparisons had worst difference `1.3e-76`.

The comparison script was copied into `scripts/referee_work/` before execution, so its table writer created a **new** table there and did not overwrite the supplied one. `referee_chain_compare.out` and `referee_r2fold_rerun.out` print:

```
monomials: rebuilt 307  notes 307  same key set: True
coefficients differing (exact, sympy expand of the difference): 0
x-exponents: {-6: 39, -5: 44, -4: 44, -3: 44, -2: 44, -1: 44, 0: 48}
largest weight of a term of the table (atoms + zeta weights + degree in lA0,lA1,lB): 5
120 x [zeta5] - Z_e^(5) = 0 ;  120 x [zeta2 zeta3] - Z_e^(2,3) = 0  (after A1 = 1-B, A0 = 1-B-B kappa)
wrote g9_hex5_table.txt: 307 lines, 5348394 bytes, sha256 f0fd1f7a771301273e14f4f8f94a5e109e50f78c10c3265a23dc080bf9c2030a
13 harmonic polylogarithms at x = -1: worst |series (CVZ, 170 terms) - closed form in the constants| = 9.15e-111
all 20 coefficients equal the integer relation of the referee (weighted convention): True
[zeta5] = -20187/2000 sqrt5 ; [zeta2 zeta3] = -858/125 sqrt5 ; remainder outside the span of the 20 constants: 0
E_e (r=2 edge, weighted convention) = -83.5391363475376327634345072141820934358
```

The original and regenerated table bytes are identical. Their SHA-256 is `f0fd1f7a771301273e14f4f8f94a5e109e50f78c10c3265a23dc080bf9c2030a`. Exact-type inspection found **zero SymPy Float expressions among all 307 coefficients and all 210 J expressions**. These observations support the source-level exactness audit; they do not turn the quadrature checks into proofs.

The optional r=2 folding script names a dictionary `fitted` and calls `nsimplify`, but its input coefficients are already exact rational coefficients extracted from a symbolic polynomial. That is not a fit to decimal data. This script is outside the generic S1–S6 theorem construction.

## Independent recomputations

Runnable code: `referee_checks.py`. It imports no author modules. Commands use `<project>/.venv/bin/python`.

### M1 and the shuffle matrices

`python referee_checks.py algebra`, output `referee_algebra.out`.

I first integrate F^−2 in t, obtaining `1/((1+xi0*U)*(1+xi1*V))`, since the square-root numerator is independent of t. I recover its square root from S²=R coefficient by coefficient, divide by that denominator, and convert powers of exp(u)−1 using finite differences of powers. This is a different computational route from the supplied binomial expansion followed by t-integration.

```
Independent M1 differences: 0 0
Parallelogram differences: 0 0
Regular edge differences (after dividing out c): 0 0
Shuffle matrix weight 1 shape (1, 1) determinant 1
Shuffle matrix weight 2 shape (2, 2) determinant 2
Shuffle matrix weight 3 shape (4, 4) determinant 6
Shuffle matrix weight 4 shape (8, 8) determinant -96
Shuffle matrix weight 5 shape (16, 16) determinant 5760
```

The determinants use my enumeration order. Shuffle coefficients were obtained by choosing letter positions, independently of the supplied recursive shuffle implementation.

### Regular hexagon from the definition

`python referee_checks.py numeric regular`, output `referee_regular.out`. Polar radius `R(theta)=F_p(theta)^(-1/p)` and speed `R*sqrt(1+(F_theta/(p*F))²)`; no edge-jet engine or closed form enters the finite differences. Precision 110 digits; eight Richardson levels from h=0.001. Printed results:

```
regular closed form: 656.41536196762475447664760905831120964826426464898548604185355096
regular polar FD:    656.41536196762475447664760905831120964826425381407381657416382461
regular difference: -1.0834912e-41
```

This is substantially more than the requested 25-digit numerical check, but is numerical evidence, not interval-certified accuracy or my own exact derivation of the entire regular-hexagon formula.

### Entire text table at new weights

`python referee_checks.py numeric table`, output `referee_table.out`. Corners `(1,1), (-27/10,33/10), (-2,1)`; **new weights (2,3,4)**. I read the plain-text table, not its pickle. Coefficients are evaluated with multiprecision constants; the 13 HPLs are evaluated by Taylor-coefficient recurrences from their differential definitions, using 650 coefficients at 110 digits, rather than the supplied nested-sum evaluator.

```
Read text-table rows: 307
table edge 0 r 13/10 value -487.11543574317237386639146448654440113644058927775351831934706976
table edge 1 r 1/2 value -638.96596377980454670733374823918392677968041194283714156591838095
table edge 2 r 20/13 value -117.08459762409190723517854484044642269929704609502668149992284827
table perimeter derivative: -2486.331994294137655617807515132349501230836094631234682770376597964044
new weights (2,3,4) polar FD: -2486.331994294137655617807515132349501230836100788855509334546525115222
table minus polar FD: 6.157620827e-42
```

The comparison is against my independent polar finite differences. All three edge ratios lie within the Taylor disk.

### Octagon symbol from the 252 rows

`python referee_checks.py symbol`, output `referee_symbol.out`. I used SymPy's exact Q(sqrt(2)) field, exact finite S-unit factorization (no floating search), my own expanded-tensor recursion, multiset permutations and position-enumerated shuffles. The row input was freshly regenerated by supplied P1. Thus this is an independent **symbol computation from the rows**, not an independent second derivation of all row coefficients.

```
exhaustive binary shuffle controls: 128
independent Li5 symbol controls: passed
independent block local [0, 0]
independent block single-p [-53955/8 + 4770*sqrt(2), -10215/2 + 14445*sqrt(2)/4]
independent block single-q [-53955/8 + 4770*sqrt(2), -10215/2 + 14445*sqrt(2)/4]
independent block mixed [48015/4 - 8505*sqrt(2), 11115 - 15525*sqrt(2)/2]
INDEPENDENT TOTAL: [-1485 + 1035*sqrt(2), 900 - 540*sqrt(2)]
```

The shuffle controls cover every binary-word product of weights k and 5−k for k=1,2,3,4, not random samples. The two Li5 controls compare complete tensors for G(0,0,0,0,p;1) and G(0,0,0,0,q;1) with the classical formula. The pairing is nonzero in Q(sqrt(2)), not merely nonzero as a floating approximation.

Additional algebra printed in `referee_extra.out`:

```
Diamond displayed form - 120*a5 (after pi2=6*z2): 0
Octagon non-proportionality, two possible multipliers:
  -495/28 + 345*sqrt(2)/28
  270/19 - 225*sqrt(2)/19
Their difference: -16965/532 + 12855*sqrt(2)/532
Local moment (a,b)=(0,5), coefficient of H5: -24*DEN
```

Here a5 is the displayed expression at live `paper/ieee/main.tex:1159`; this checks its algebraic agreement with the diamond formula, without re-proving Proposition 1.

### Independent exact spot certificates for every J entry

`python referee_jtab_check.py`, output `referee_jtab_check.out`. This reads the freshly regenerated J expressions, independently expands HPLs from their differential recurrences, and computes the expected series coefficients with a separate log-Beta exponential recurrence. Every J entry is checked from x^−6 through x^3, including cancellation of all apparent principal parts. It also checks the explicit zeta-monomial weight bound entry by entry.

```
Fresh J entries: 210
Independent exact coefficients checked per entry: x^-6 through x^3
Coefficient equalities checked: 2100 mismatches: 0
Largest explicit zeta-monomial weight: 4
```

These are finite exact coefficient checks, not a proof of equality of the complete analytic functions. They independently supplement the deterministic summation derivation and quadrature checks; they do not replace that derivation.

## Supplied-code reruns, distinguished from independent work

- `referee_beta_rerun.out`: 504 moment entries, zero rule failures; worst quadrature difference `4.8e-35` at 40 digits. These are the supplied script's checks.
- `referee_local_rerun.out`: exact differences zero for diamond, regular hexagon and shear; 91 regular-hexagon integrand monomials. Printed regular coefficients `2835*sqrt(3)/32` and `1395*sqrt(3)/32`.
- `referee_rows_rerun.out`: supplied P1 regenerated 252 rows and compared 1303 coefficients with the notes, zero mismatches.
- `referee_symbols_rerun.out`: supplied P2 reproduced both pairings and the pairing of C. This does not prove the entire classical-plus-C symbol identity.
- `referee_rul1_rerun.out`: external read-only `verify_rul1.py`, run with bytecode writes disabled, rebuilt the rational DS/duality relation system and passed **41/41** checks in its printed 55 seconds. It also confirms the unit-normal-convention fold against the historical data. I inspected its exact linear-algebra construction but did not independently prove every relation generator or rebuild the named-constant substitutions. This is a verifier rerun, not an independent certificate verifier.
- Fresh S1–S6 completed in `scripts/referee_work/`; the complete exact comparison and fresh r=2 fold passed. See the terminal results above. Existing author logs were not substituted for this run.

## Coverage of the printed claims

Locations below refer to `main_tex_statements_l1715-1778.tex`, except Table II locations, which refer to the live external `paper/ieee/main.tex`. The excerpt was used for the statements; the external file was read for Table II, verification-status wording and the diamond a5 formula.

### Proposition 9 and its adjacent text

| Printed sentence or clause | Evidence and status |
|---|---|
| :12–19: local Beta evaluation and the two printed laws | Hand Beta/extraction plus exact finite identity M1; **PROVED**, independently recomputed. |
| :21: fully local polygons contain no further explicit z5 or z2*z3 in this form | Hand log-degree argument; **PROVED** as a representation statement. Independence of actual numbers is an assumption, not a consequence. |
| :23–24: both regular 2m-gon sums | Hand substitution; independent exact check; **PROVED**. The sums are local coefficients for general m, not a formula for the complete non-local derivative. |
| :28, “proof ... at weight five” and raw-integral difference zero | Supplied M1 and my integration-first computation establish the stated finite identity. |
| :28, every parallelogram and both short sums | Fully local argument and substitution; **PROVED**. |
| :28, “second exact derivation, independent of the first” | Provenance claim about the notes. Not established by merely rerunning M2. This review supplies its own exact M1 route, but not an independent exact full-length computation. |
| :28, diamond closed form | Supplied M2 rerun exact zero. Agreement with the displayed a5 also checked algebraically in this review. Proposition 1 itself was not re-proved. |
| :28, shear and the quoted quadrature precision | Exact form confirmed by M2 rerun. The historical quadrature and its quoted residual were not rerun; retain as an attributed numerical check. |
| :28, regular-hexagon zeta coefficients | Exact M2 rerun and local law; **PROVED**. |
| :28, high-precision polar runs and earlier finite differences | Historical numerical evidence, not independently repeated at those settings. My independent polar check supplies its own printed result above. |
| :28–34, 91 monomials and the displayed regular-hexagon closed form | Exact M2 derivation, rerun; **PROVED**. Historical integer relations discovered it but are not ingredients of the proof. |
| :35, neither half-argument Li4 nor Li5 is needed for the regular hexagon | The displayed exact formula proves existence of such a representation. No numerical independence is required for existence. |
| :35, non-local hexagon edge at r=2 in twenty constants | Additional corollary requiring exact special-value reductions and a fold. RUL1 and the fold were rerun successfully; all twenty coefficients agree exactly. The generic theorem alone does not supply these reductions. |
| :35, repeated PSLQ relation, maximum precision and unit-normal conversion | Historical numerical/provenance claims. They are not ingredients in S1–S6. No fresh integer-relation search was run here. |
| :35, “far kink lines add ...” | Exact difference of coefficients in the specified constant presentation, once the exact fold is accepted. Not a canonical coefficient of a number. |
| :35, basis-independent occurrence of half-argument Li values | Needs the qualification in finding 8; the symbolic expression alone does not establish numerical minimality. |

### Theorem 26 and its adjacent text

| Printed sentence or clause | Evidence and status |
|---|---|
| :40, all nondegenerate hexagons, arbitrary positive weights, opposite-edge sum | Hand arc-length formula, central symmetry and dominated differentiation; **PROVED**. |
| :40–45, explicit polynomial in thirteen functions, weight <=5, denominator L^9 | Hand framework plus exact S1–S6; **PROVED AFTER REPAIRS** in the explanatory proof, not a numerical fit. |
| :45, cancellation at x=0 and the local limit | Analytic integral and identity theorem on the punctured slit plane; **PROVED** with finding 7's wording. A numerical limit check is unnecessary for this conclusion. |
| :45, both explicit zeta coefficients equal the local laws for every r | Correct no-new-constants argument plus local extraction; **PROVED AFTER REPAIRS** to the factor 120 in the supplement. |
| :48, table checked against quadrature | Numerical check of derived expressions. It does not replace the symbolic derivation; only forty of the 210 keys are checked at four kappa values in the stated 160-check block. |
| :48, power-series domain and continuation to all r>0 | Hand argument valid as discussed above; not an assumption of series convergence outside its disk. |
| :48, differentiation under the integral sign | Supplied smoothness lemma is stated for all k and covers k=5. Its bound and use are sound on inspection. |
| :48, independent referee's six hexagons and quoted digits | Statement about another review. This pass did not rerun those six tests or validate their historical settings; it independently checked one new-weight hexagon and the regular one. |
| :48, eleven Nielsen words and the two other words | Word classification is immediate from the list. Their occurrence in the table is exact. An intrinsic necessity/minimality claim is not shown by this fact; use finding 8's wording. |
| :48, non-local z5 carried by HPLs rather than an explicit coefficient | Correct **for this formal representation**. Special values of HPLs may of course contain z5 and z2*z3. |
| :48, exact r=2 coefficients | Requires RUL1 plus fold, beyond S1–S6. Certificate and fresh fold reruns passed; all twenty coefficients agree exactly. |
| :48, r=1 gives the local law exactly | Hand argument proves the entire local edge value, not just its two zeta coefficients. |
| :48, unit-normal edge comparison and golden-hexagon finite differences | Historical numerical checks; not rerun here. The writer's weaker numerical replication is not a counterexample to a stronger historical precision claim. Do not promote either historical check to proof. |

### Proposition 10 and the following paragraphs

| Printed sentence or clause | Evidence and status |
|---|---|
| :51, definition of the weight-five part | Needs the associated-graded projection in finding 2. Total log degree before integration is not literal homogeneous weight afterward. |
| :51, “It has a nonzero (3,2) part” | **PROVED AFTER REPAIRS** when “it” is the symbol of the projected W5 and “part” is the antisymmetrised deconcatenation class modulo shuffle products. |
| :51, the two exact pairings | Supplied P1/P2 and my independent symbol computation from all rows agree exactly; **PROVED** under the specified regularized symbol convention. |
| :51, no Li5-plus-products symbol | S-unit case hand-proved. Arbitrary field arguments require finding 5's extension or an explicit theorem and identification; **PROVED AFTER REPAIRS** with the elementary extension given here. |
| :51, not a multiple of M4's pair | Exact scalar comparison, independently repeated; **PROVED** relative to the quoted M4 calibration. |
| :51, obstruction for the number assumes a period conjecture | Conditional statement only. No unconditional numerical nonclassicality theorem is established. “Standard period conjecture” should specify the intended period realization/injectivity if this implication is used substantively. |
| :54, all 756 integral rewrites and M4 calibration | Described by external verifier and writer report, but not independently recomputed here. This review's pairing calculation uses the top-weight reduction directly from the 252 rows. |
| :54–58, full symbol is classical plus C | Not implied merely by our matching pairing of C. Full rank/symbol identity not checked in this pass; **NOT PROVED in this review**. |
| :59, octagon numerical value and high-precision history | Numerical/provenance claim; not recomputed from the perimeter in this pass. The independent symbol computation does not verify the numerical value. |
| :59, depth two enters through non-local octagon terms | Correct at the symbol-obstruction level: independent local block zero and total nonzero. Avoid an unconditional claim about the numerical period. |
| :61, decagon, 12-gon and complete nine-Li32 octagon formula | Outside the three statements' proofs. Not checked here; the manuscript itself calls the complete closed forms numerical. A vanishing symbol obstruction alone does not establish a full numerical closed form. |

### Table II rows

- Live `main.tex:2174` (order-five edge law and examples): “Proved” is justified, with computer assistance disclosed; quoted quadrature precision remains historical numerical evidence.
- `:2175` (octagon needs depth two): retain **symbol level**, and replace its label as in finding 9. The value column is a separate numerical claim.
- `:2196` (every hexagon, both orders): the order-five exact computer-assisted construction is justified; this pass does not recertify all claims of the earlier order-three theorem.
- `:2197` (regular hexagon closed form): exact local Beta derivation justifies “Proved”; integer relations are checks/discovery only.
- `:2204` (non-local order-five hexagon edge): exactness additionally rests on the x=−1 certificate system and fold. It is not derived from a numerical integer relation once those exact reductions are supplied.
- `:2205` (decagon no depth two): outside this review; symbol level must remain explicit.

“Refereed” is defined by the live paper at :2077 as a separate AI verification pass, not journal peer review. This review does not silently reinterpret that word.

## Not checked, limitations, and reproducibility

- I have not independently reimplemented the entire symbolic S1 engine, all 210 general J summations, or S6 assembly. A fresh exact rerun and independent numerical values are different kinds of evidence.
- I did not prove independence of numerical constants/functions, a period conjecture, the claimed minimal number of HPLs, or necessity of particular numerical constants without the stated assumptions.
- I did not rerun historical PSLQ/LLL searches, golden-hexagon settings, the six-hexagon review, the high-precision shear quadrature, or the full octagon/decagon/12-gon numerical closures.
- I did not verify the full classical-plus-C rank computation or independently derive M4's calibration.
- For the arbitrary-argument classical-symbol obstruction, I supplied an elementary tensor argument. I did not rely on an unverified claim that a particular theorem in the cited papers is exactly this pairing. The cited [Matveiakin–Rudenko paper](https://arxiv.org/abs/2208.01564) and [Charlton–Gangl–Radchenko–Rudenko paper](https://arxiv.org/abs/2210.11938) were located online; their full theorems were not audited here.
- The supplied shell chain uses process exit status; several comparison scripts print discrepancies without asserting failure. A zero shell exit status alone is not the verification criterion. I inspect the actual comparison outputs.
- No PDF rebuild or visual-layout review was performed; this is a mathematical and computational review.
- My finite-difference comparisons are high-precision numerical checks, not rigorous interval enclosures. Exact algebra checks and symbol checks use assertions; the script fails on a mismatch.

The full classical-plus-C sentence could truthfully be replaced in this supplement by: “A separate calculation in the notes claims the full symbol decomposition; here we verify that C has the same computed (3,2) pair as W5.” The latter assertion is all the P2 check establishes.

**Final verdict:** Proposition 9 is PROVED. Theorem 26 and Proposition 10 are PROVED AFTER REPAIRS specified above. Unconditional assertions about numerical independence or nonclassicality remain NOT PROVED. No supplied source file was edited.
