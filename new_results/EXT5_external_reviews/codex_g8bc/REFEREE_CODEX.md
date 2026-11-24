# Independent adversarial referee report

Review completed on the supplied copies. No supplied file was edited. New runnable checks and their outputs are listed below. **No fatal mathematical defect was found in the new function-level argument or the real-letter integral derivation.** There are mandatory wording/definition repairs and missing reproducibility inputs.

## Verdicts

“PROVED” below assesses the new argument using its explicitly cited earlier results. It does **not** mean that the earlier all-N machine proofs were independently reproduced: their essential files are absent. This distinction is material for both the M3 expansion and the exact identity with Y.

| Claim | Verdict | Qualification |
|---|---|---|
| Hyperlogarithm independence over C(c), including mixed lengths | **PROVED** | The minimal-relation and residue argument is valid. |
| Classical factors have finite-alphabet expansions with the stated components; every classical expression has delta22(S4^b)=0 | **PROVED** | Arbitrary weights, products, complex integration constants and rational prefactors are covered. |
| G8B(a): no interval of rational-argument classicality | **PROVED** | Uses the earlier M3 expansion, M5, as an imported result. No repair to the new obstruction argument is needed. |
| G8B(b): function invariant S4(M3), identification with Omega and displayed nonzero cobracket | **PROVED** | Same imported expansion. I independently verified the printed 56-term symbol and its cobracket. |
| W3(a): absolutely convergent real-letter integral equals M2-W2 | **PROVED** | N1 independently recomputed exactly; the analytic derivation is valid. |
| W3(b): length-three real-letter representation, including the “64 terms” wording and definition | **PROVED AFTER REPAIRS** | There are 64 product-logarithm integrals, or 105 ordered word terms after the stated shuffle expansion. Clarify regularisation. |
| W3(c): equality with Y(c)+Y(-c)-zeta(2)(1-c²)artanh(c)/(2c) | **PROVED** | Follows algebraically from the imported earlier M2 theorem. Independently confirmed numerically here, not re-proved for all c. |
| Main-text unqualified “an iterated integral” and “not classical as a function of the angle” | **PROVED AFTER REPAIRS** | Specify a linear combination and rational arguments in c; see exact wording below. |
| Standalone reproduction of the earlier M5 and M1/M2 machine proofs from this packet | **NOT PROVED** | Required data/libraries/scripts are missing. This is a reproducibility limitation, not a counterexample to the new reasoning. |

## Numbered findings and exact repairs

1. **Must fix — missing proof/check inputs.** `scripts/g8b_check.py:69-79`; `sec_g8_M3alpha.tex:77-84`; `sec_g8_M3func.tex:215,277`; `sec_g8_M2alpha.tex:113`.

   The supplied G8B checker actually fails with `FileNotFoundError` for `scripts/g8_terms_k3_part1.json`; `scripts/g8_terms_k3_partz2.json` is also absent. The complete 74-term R3 and 14-term Psi lists, the M5 all-N checker and its library are unavailable. Consequently I cannot independently inspect the asserted 18 shorter R3 terms or Psi's 6 length-one and 8 length-two terms. The saved `.out` reports these counts but is not a reproduced computation. The exact earlier M2 proof used for W3(c) also calls absent `g8_allN.py`/`g8_zsum.py`.

   **Repair:** include both exact JSON inputs and all transitive dependencies of the named all-N checks, with the invocation and recorded results. Alternatively explicitly identify these expansions/theorems as imported previously established premises and do not describe this packet as a standalone reproduction. The new theorem's logical dependence on M5 is already disclosed correctly; the files are needed to verify that premise independently.

2. **Must fix — distinguish regularised integrals from ordinary simplex integrals.** `sec_g8_W3int.tex:33-37,47,206`.

   The conventions define the words actually needed, but the universal recursive definition is incomplete: applied to the pure word (0,0,0), “otherwise” gives the divergent integral of log²(y)/(2y). Also the displayed nested-simplex formula is not an ordinary convergent integral for a word ending in zero. For example I(j,0,0;Y), j>=1, uses the regularised inner value log²(y)/2. The assertion that “every integral converges” and that only the pure word (0,0) needs a convention obscures this distinction. The expanded representation contains **21 words ending in zero**, including (j,0,0) for j=1,2,3,4.

   **Repair, replacing the definition/convergence qualification:** “Set I(empty;y)=1 and I(0^r;y)=log^r(y)/r! for every r>=1, using the positive real tangent at 0. For all other words use the displayed recursion. Interpret the simplex notation with this trailing-zero regularisation. For all words occurring in the formula, the remaining recursive improper integrals converge absolutely after these values have been assigned.” The stated length-two shuffle identities then hold exactly. No extra endpoint correction or lower-length term is introduced with this convention.

3. **Must fix — the count 64 is not a count of individual ordered words.** `sec_g8_W3int.tex:47,206`; `REPORT_G8C.md`, “The exact reading proved”.

   The nonzero quadratic monomial coefficients number **10,11,13,15,15** for outer indices j=0,...,4, totaling 64. Each off-diagonal product lambda_a lambda_b gives two words; each diagonal gives twice one word. There are **23 diagonal** and **41 off-diagonal** monomials, hence **23+2*41=105 distinct ordered words** with nonzero coefficients. There are 128 summands if the two identical diagonal words are left uncollected. The TeX does explicitly count product-logarithm integrals, but the final singular “an iterated integral” and the report/user shorthand conflate the two counts.

   **Repair:** “W3 is a sum of 64 integrals of products of two logarithms against a logarithmic differential. After shuffle expansion this gives a Q(c)-linear combination of 105 ordered, regularised length-three iterated integrals, with coefficients in Q+Qc+Q/c.” This is the count in this expansion, not a minimality assertion.

4. **Must fix — main-text singular integral.** `main_tex_statements_l1500-1580.tex:12`; also `sec_g8_M2alpha.tex:15`.

   The literal assertion that W3 is one word integral with a constant coefficient is not established. A single ordinary integral with a quadratic logarithmic integrand is established, but that is different terminology.

   **Exact replacement:** “For 0<alpha<pi, M2(alpha)=W2+W3, where W2=zeta(2)/2-Li2(c²)/2 and W3 is a Q(c)-linear combination of regularised iterated integrals of length three, with real letters 0,-1,-1/(1-c²),-1/(1+c),-1/(1-c) and endpoint (1-c²)^(-1/2).” The following explicit formula in terms of Y retains its existing dependence on the earlier M2 theorem.

5. **Must fix — qualify the main-text non-classicality/depth assertion.** `main_tex_statements_l1500-1580.tex:48`; compare the correct scope at `sec_g8_M3func.tex:38-41,271`.

   G8B excludes finite C(c)-linear combinations of products of constants, log g and Li_n(g), with g rational in c. It does not exclude algebraic arguments, arguments rational in alpha instead of c, arbitrary other functions of alpha, or classical values at individual angles. It does not prove that subtracting the two proposed Li22 functions leaves a classical function. “It needs depth-two functions” must not be read as a proved depth-two-plus-classical decomposition.

   **Exact replacement:** “On no nonempty open interval in c=cos(alpha/2) can M3 be expressed as a finite C(c)-linear combination of products of constants, logarithms and classical polylogarithms of rational functions of c, even allowing arbitrary polylogarithm weights and branches as in Definition G8B. The two displayed depth-two functions match its (2,2) invariant; this does not establish that their subtraction leaves a classical function.” If a shorter sentence is desired, use “Thus M3 is not classical with rational arguments in c, in the sense of Definition G8B.”

6. **Minor — false lower bound in an otherwise valid limit argument.** `sec_g8_W3int.tex:180`.

   The statement t-k>=2 for c<=1/2 is false: at t=1,c=1/2 it is 3/2.

   **Exact repair:** replace it by “t-k>=2s²>=3/2”. This supplies the required uniform bound and leaves the dominated-convergence proof intact.

7. **Minor — the exceptional set belongs to a presentation.** `sec_g8_M3func.tex:20,196,226`.

   E(C) is not intrinsic to the meromorphic function C: one may add cancelling factors with new exceptional points. The class itself is well defined, and the proof fixes a finite presentation, so there is no logical failure.

   **Exact repair:** write “For a chosen finite presentation of C, let E be the union of ...” and use that E in the alphabet enlargement. Holomorphic classical factors and meromorphic rational-prefactor combinations on U are otherwise consistent.

8. **Minor — remove stale gap statements when these sections are integrated.** `sec_g8_M3alpha.tex:44,62,126,234`; `sec_g8_M2alpha.tex:35`.

   These passages still say function invariance or the real-letter integral is not proved. Their present-tense global reading conflicts with the new sections, although their original section-local caution is understandable.

   **Repair:** replace the old gap descriptions by cross-references to G8B and the W3 integral theorem, preserving the remaining caveats on algebraic arguments, special values and the converse criterion.

## Detailed audit of A

- **Minimal relation, normalisation, residues — accepted.** `sec_g8_M3func.tex:77-97`. Lexicographic minimisation of (maximum length N, number k of nonzero coefficients at N) is well founded. Dividing by a nonzero rational coefficient creates only finitely many additional poles; this is allowed. Differentiation leaves the normalised top coefficient zero and introduces only shorter words from differentiating L. Unless the resulting relation is coefficientwise zero it contradicts minimality. Therefore all top coefficients are constants and the next-length equations are rational identities. A derivative of a rational function has zero residue at every finite point, including points outside U. Taking residues at distinct letters forces all top constants to vanish. Poles of coefficients at letters or elsewhere do not evade this argument. Mixed lengths are explicitly allowed. On simply connected U the chosen primitives are single-valued; multivalued analytic continuations do not create a relation among these fixed branches.
- **Finite alphabet and enlarged branch class — accepted.** `:9-23,130-138,177-210`. A finite expression uses finitely many rational arguments, each with finitely many finite zeros/poles of g and 1-g. Their union suffices. The recursive definition of Pn intentionally includes more than genuine analytic-continuation branches of Li_n, through arbitrary constants at successive integrations. It still satisfies the stated finite expansion. Constants at constant arguments are already arbitrary complex constants.
- **All component degrees — accepted.** `:184-193`. For 0<=j<n, S_j^b(Pn)=P_(n-j)(b)[g]^tensor j; for j=n it is -[1-g] tensor [g]^(n-1); for j>n it is zero. This includes n=1 and j=0. The formula is unchanged by the allowed integration constants except for the displayed base-point values.
- **All products, arbitrary weights and rational prefactors — accepted.** `:195-210`. In total degree four, two or more positive component degrees give a shuffle product killed by delta22. A sole degree-four factor is zero for n<4, has the repeated [g] slots for n=4, and is a constant times [g]^4 for n>4; every case is killed. Degree-zero factors contribute their values at b, including factors that themselves have large positive length. Rational prefactors are scalars for S_k^b and delta22, exactly as defined. Their derivatives are needed only in the independence proof. No hidden weight bound or purity assumption is used.
- **Base point and cancellation — accepted with the stated restriction.** `:59-68,99-128,158-171`. The constant-coefficient triangular base-point formula proves invariance of length and the top component. Thus S4(M3) itself, not merely its cobracket, is invariant because length(M3)<=4. S4^b of a general longer function need not be invariant, nor need its cobracket. A counterexample to that stronger, unclaimed statement is F=L^3_(2,1,1,0,0), alphabet {0,1,2}, with b=3 and b'=4 in a suitable domain avoiding the letters. S4^3(F)=0, whereas S4^4(F)=log(4/3) e0 tensor e1 tensor e1 tensor e2, with nonzero delta22. The proof correctly compares a classical competitor and M3 at the **same** base point; its classical vanishing lemma holds at every fixed base point. Cancellation of different lengths cannot violate uniqueness.
- **Interval and alphabet enlargement — accepted.** `:225-232`. Removing a finite exceptional set leaves a point of every nonempty interval. A small complex disc through that point avoids all coefficient poles and letters. Equality on its real diameter gives equality of holomorphic functions by the identity theorem. Independence on the enlarged finite alphabet and injectivity of the exterior-power embeddings preserve the obstruction.
- **Lower terms — accepted subject to the missing imported data.** `:215-223`. The displayed Phi formula has length <=1 even after rational multiplication. If the asserted Psi list has lengths <=2 and the remaining R3 list lengths <=3, multiplication by zeta values cannot raise length, since they are complex scalar constants. This proves S4(M3)=S(Omega). I verified the Phi formula and the full printed Omega table; I cannot independently certify the unprovided shorter lists. Weight and length must not be interchanged here: zeta(3)Phi may have weight four but length at most one.
- **Ordering, [2], two terms — accepted.** `:114-155,239-253`. Reversal places the outermost integration letter last; f1=-(c-1)^(-1) supplies the minus signs. In this function invariant [2]=0 because d log 2=0. This is a quotient of the letter space, not a claim that log 2 vanishes as a value or coefficient. Pulling back the printed Li21/Li22 differentials at x=1/2 gives top tensor [1-y/2] tensor [y] tensor ([y]-[1-y]) tensor [y]. It gives delta22(f_-)=-p wedge r and delta22(f_+)=q wedge r, independently reproducing the two-term identity. A formal symbol retaining [2] is a different object; the new section correctly does not treat its residual components as this invariant.

## Detailed audit of B

`sec_g8_W3int.tex:59-99`: the domination g^(-1/2)(1+|log g|)^(2j) is sufficient for real-u differentiation under the integral sign. R and each required derivative are bounded on [0,1] for each fixed angle in (0,pi). The u² coefficient, including the subtraction term, is correct.

`:103-153`: the Euler inverse, Jacobian and real domain are correct. One can justify monotonicity immediately from R>|g-k|, giving 1+(g-k)/R>0. The five real poles are distinct throughout 0<c<1. Independently recomputed N1 shows that every complex pole lies in V', with the stated real simple-pole residues and zero remainder.

`:157-180`: both integration-by-parts boundaries vanish, using L(t1)=0 and V L²=O((t-1)log²(t-1)) at the lower endpoint. The identity V L'=1/(2(t+1)) is exact. Differentiation in c gives W2'=log(1-c²)/c; the moving-endpoint contribution vanishes. The c->0 limit is zeta(2)/2, after correcting the bound in finding 6.

`:183-208`: the y change of variables, endpoint 1/s and logarithms have no missing additive log s. All nonzero letters are negative. The unsplit j=0 integrand is integrable because ell_m=O(y); after expansion each required recursive integral is also absolutely convergent with the regularisation in finding 2. For mixed inner words I(0,b), I(b,0), b!=0, the sizes are O(y) and O(y(1+|log y|)), respectively. The four I(j,0,0), j>=1, integrate log² y against a bounded kernel. The divergent pure word (0,0,0) does not occur.

## Main-paper sentences: support and truthful reading

1. `main_tex_statements_l1500-1580.tex:12`, the decomposition and explicit W2: **PROVED** by B(a). “W3 is an iterated integral ...”: **PROVED AFTER REPAIRS**, using the replacement in finding 4. The following explicit Y formula is supported through the earlier M2 theorem and independently checked numerically here.
2. `:40-45`, Omega's function-level (2,2) part, ordering and nonzero displayed formula: **PROVED**, using the imported M3 expansion and the definition of the length-four invariant. Suggested clarification: “Its length-four function invariant, with letters modulo nonzero constants and rational prefactors treated as scalars, has the following (2,2) part.” No claim about ordinary differentiation of rational prefactors is intended.
3. `:48`, “So M3(alpha) is not classical as a function of the angle: it needs depth-two functions”: **PROVED AFTER REPAIRS**, as in finding 5. The HPLs' identification with multiple polylogarithms is an earlier representation fact, not a new converse theorem.
4. `:48`, “One depth-two term is not enough ...”: **PROVED with the stated restriction**, made explicit by `sec_g8_M3func.tex:234-235`: no M3=hG+C with rational h, classical C and G a constant-coefficient finite hyperlogarithm combination on a common finite alphabet. The two nonzero components would have constant ratio, whereas (1-c)/(1+c) is nonconstant. Do not extend this to unrestricted G; trivially G=M3 would invalidate such a claim.
5. `:48`, “Two terms match it: modulo constants ...”: **PROVED for the function invariant**, with the exact Li22 convention of the supplement. I independently recomputed the constant-free equality. The parenthetical assertion of four remaining components when [2] is formally retained is an older formal-symbol computation, **not newly verified here** and not part of the function invariant proved by A.
6. `:48`, “We do not claim that Omega minus these two terms is classical”: **correct and necessary**. The next “This is proved at the level of symbols” should specify its antecedent: replace with “The displayed equality of (2,2) parts is proved for the constant-free length-four function invariant.” It must not be read as proving a classical remainder.
7. All subsequent sentences in that paragraph about a period conjecture, the square and hexagon specialisations, special-value Lie classes, S-unit equations, lattice reductions, CMZV proofs and numerical accuracies are **NOT PROVED by A or B**. Some are explicitly imported earlier claims; this review does not adjudicate their separate proofs. Truthful connective wording is: “The following special-angle statements require their separate arguments; function-level non-classicality with rational arguments does not decide them.” In particular no unconditional non-classicality of any single number, and no obstruction for algebraic polylogarithm arguments, follows from G8B.

## Independent recomputations and numbers

Runnable checks, requiring only the installed Python environment:

```
<project>/.venv/bin/python referee_exact.py
<project>/.venv/bin/python referee_numeric.py
```

Both scripts import no supplied checker and contain assertions that fail on disagreement.

**Exact (`referee_exact.py`, `referee_exact.out`):** direct u² radicand/exponential jet; N1 A/B/C identities with zero rational remainders; V formula, endpoint and IBP kernel; y identities and all letters; counts 64/105/21; printed 56-term Omega -> 24 tensors with the stated cobracket; 768 basis shuffle cases on four letters; formal Li_n chains n=1,...,8 in every degree 0,...,n+1; the pulled-back two-Li22 cobracket. A search for a rational relation among 7 words of lengths <=2 on alphabet {0,1}, with coefficients in span(1,c,1/c,1/(c-2)), has a 51-jet matrix of rank **28/28**. This is a finite counterexample search, not a substitute for the independence proof.

**Numerical (`referee_numeric.py`, `referee_numeric.json`, `referee_numeric.out`):** M2 is extracted from the defining M(u;alpha) integral after g=exp(-x), using rationalisation/expm1 for endpoint stability, 48 Fourier nodes, radius 1/8 and 70 working digits. This uses neither the differentiated fA/fB/fC formula nor the Euler decomposition to compute the reference M2. Independent real-letter quadrature and independent HPL coefficient summation are compared to it. This is numerical evidence, not a certified error bound.

| alpha | M2 from definition |
|---|---|
| pi/3 | -0.042871754976852882532306425052508863463709733192 |
| pi/2 | -0.022161911087269399423754937739791514977955046329 |
| 2pi/3 (obtuse) | -0.009567895652956384775719722385651527995008301633 |

| alpha | W3=M2-W2 | absolute error versus real-letter integral | absolute error versus Y expression |
|---|---|---|---|
| pi/3 | -0.376104091935813048896980675113240709192442546486 | 2.633e-50 | 2.633e-50 |
| pi/2 | -0.553508681278876364708634360982964055215330759529 | 5.637e-51 | 5.637e-51 |
| 2pi/3 | -0.698208609535703299552335391464758334694554718780 | 3.510e-52 | 3.510e-52 |

The W2 integral agrees with (zeta(2)-Li2(c²))/2 at working precision (reported difference 0 at each angle). HPL truncation orders were 1301, 552 and 286. The displayed values are rounded to 48 significant digits, within the observed agreement. The JSON preserves more digits, some beyond the demonstrated agreement.

## What I did not check

I did not independently reproduce the absent all-N M5 or M1/M2 programs, inspect the absent complete R3/Psi JSON lists, validate the other earlier machine steps or special-angle period/CMZV claims, verify external bibliographic attributions, prove a classical-remainder converse, handle algebraic arguments/prefactors, establish minimality of either representation, formally certify quadrature/software arithmetic, or rebuild/visually QA the PDFs. The TeX was the source of record. The supplied G8B check fails because of missing inputs; its historical output was not treated as a successful rerun.
