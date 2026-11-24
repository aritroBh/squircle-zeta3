# Independent adversarial referee report

Scope: the supplied copies only. Line numbers below are local to those copies. I edited no supplied file. I found no counterexample to the three explicit coefficient formulas and no fatal defect in the all-n computation. Several statements require repairs or exceed what these sections prove.

## Verdicts

| Result or sub-claim | Verdict | Reason |
|---|---|---|
| Theorem 20: displayed HPL formula for every 0<c<1 | **PROVED** | Sound all-n reduction, exact M1/M2, including odd and low coefficients. |
| Theorem 20: W_3 is a Q(c)-combination of weight-three HPLs and zeta(2) times weight-one HPLs, with real letters | **PROVED** | Direct subtraction and the sign rule. |
| Theorem 20: literal single weight-three iterated-integral assertion / notes' Euler-substitution representation | **NOT PROVED** | The supplement expressly supplies the preceding weaker formulation instead. |
| Theorem 20: ten classical reductions, eight arguments, displayed classical formula and zeta(3) coefficient | **PROVED** | Exact derivatives, basepoints, and substitution; independent checks below. |
| M_3: the 74-term R_3 formula and explicit Psi | **PROVED AFTER REPAIRS** | Add the standalone constant to the coefficient criterion; the implemented calculation already handles it correctly. |
| Proposition 5: displayed expression's constant-free tensor and nonzero (2,2) projection | **PROVED** | Independently reproduced from the printed table and JSON. |
| Proposition 5 as worded “of the function”; resulting function-level non-classicality | **NOT PROVED** by the submitted argument | Expression independence is explicitly left unproved. There is a standard repair, described below. |
| One term r(c)G cannot match the obstruction, with G having constant symbol coefficients | **PROVED** at the stated tensor level | Exhaustive rank/ratio argument, not a search over candidate MPLs. |
| The displayed two terms match the obstruction | **PROVED AFTER REPAIRS** | State explicitly that constants are removed and that only the cobrackets agree. |
| Two terms plus classical functions equal Omega as functions | **NOT PROVED** | Neither the tensor calculation nor equality of cobrackets proves this. |
| Theorem 21: displayed M_4 formula for every 0<c<1 | **PROVED AFTER REPAIRS** | Repair the criterion's treatment of constants; change the inaccurate 22-HPL description. |
| Theorem 21: zeta(3) ladder and elementary pi^4 contribution | **PROVED** | Accepted ladder plus the correctly computed rational coefficients and elementary sums. |
| Special-corner classifications and weight-five depth-two consequences in the surrounding paper | **NOT PROVED** by these sections | Separate computations and evidence are absent here. |

“NOT PROVED” identifies the reach of this submission, not a demonstrated false mathematical result.

## Numbered findings and repairs

1. **Must fix — the coefficient criterion omits constants used in its applications.** `sec_g8_M2alpha.tex:79-84`; applications `sec_g8_M3alpha.tex:74-81`, `sec_g8_M4alpha.tex:22,57`. The displayed L has HPL terms and a series starting at n=1, but no standalone constant. Allowing empty regular words does not resolve the issue: h_{w,delta} is defined only for nonempty words. Both Psi and Psi_2 contain +1/2. The script handles this correctly: `scripts/g8_allN.py:153-154,168,208-210` excludes constants from high-order coefficients and includes them in low orders. My boundary audit prints a constant coefficient -1/2 for each harmonic-only part and zero after adding +1/2. **Repair:** add a_0+a_1 c to L, require the listed HPL words to be nonempty, and make condition (ii) concern the entire finite Laurent expression, including that polynomial. Its coefficients in degrees at least two vanish, so the proof otherwise stays unchanged. This repairs the applications without changing any computed coefficient.

2. **Must fix — “22 terms of weights one to three” is inaccurate.** `sec_g8_M4alpha.tex:19`; `main_tex_statements_l1500-1580.tex:71`. The data contain 21 paired HPL terms, with weight counts 3,6,12, plus the constant 1/2. The constant is weight zero, not an HPL term of weight one to three. The explanation at `sec_g8_M4alpha.tex:22` does not make the theorem's wording true. **Replacement:** “Psi_2 is 1/2 plus a rational combination of 21 such terms of weights one to three.” The 109 count for R_4 is correct: weights 2,3,4,5 have counts 1,8,27,73. Counts refer to the paired He convention, not to individual HPL summands after expansion.

3. **Must fix — function-level wording exceeds the supplied symbol proof.** `sec_g8_M3alpha.tex:35,59,123,138,231`; `main_tex_statements_l1500-1580.tex:40,48`. The supplement defines S on an expression, treats r(c) as a scalar, and explicitly declines to prove expression independence. Consequently “of the function” and “So M_3 is not classical as a function” do not follow from the argument as submitted. **Minimal replacement:** “For the displayed 56-term expression, in the constant-free tensor space with rational-function scalars, its (2,2) projection is ... and is nonzero.” To retain the function conclusion, establish the canonical associated-graded symbol over Q(c). A useful precise repair is the injective hyperlogarithm realization in [Brown, Theorem 1.1 and its following discussion](https://www.ihes.fr/~brown/RHpaper5.pdf): choose a common finite alphabet containing zeros and poles of all rational arguments, expand the classical functions into iterated integrals there, and clear rational denominators. This supplies the needed uniqueness of word expansions. It is not “linear independence of HPLs and classical functions” literally: those families overlap, for example H_{0,0,1}=Li_3. The repair must also specify the weight-four classical class under consideration. I have not silently inserted this missing argument into the submitted proof's verdict.

4. **Minor — the stated denominator range includes a vanishing denominator.** `sec_g8_M2alpha.tex:93` allows n+t with t>=-1 and then says none vanishes at n>=1. The listed n-1 vanishes at n=1. The actual chi/rho construction uses n+t with t>=0 (`scripts/g8_allN.py:86-106`); only 2n+t needs t>=-1. **Repair:** distinguish these ranges. Also, `sec_g8_M2alpha.tex:59`'s assertion t>=j+1 is for the negative-shift sum; for positive shifts the bound is t>=1. The conclusion remains correct in both cases. My actual-range audit found no pole at a positive integer.

5. **Must fix — specify the quotient in the two-term claim.** `main_tex_statements_l1500-1580.tex:48`; compare the precise qualification in `sec_g8_M3alpha.tex:194-215`. The two-term equality is correct for the constant-free iterated-differential symbols, with the rational prefactors treated as scalars. Retaining [2] produces a nonzero residual, independently reproduced below. Thus the qualification is mathematically necessary if the sentence could otherwise mean the refined symbol. It is redundant only after a global convention has already fixed the constant-free quotient. **Replacement:** “In the constant-free symbol, with coefficients extended to Q(c), the displayed two Li_{2,2} terms have the same (2,2) cobracket as Omega.” This asserts neither equality of full symbols nor equality of functions. Distinguishing the constant-free and refined conventions is standard; see [Duhr–Gangl–Rhodes, Section 3.3 and Section 4](https://www.maths.dur.ac.uk/users/herbert.gangl/Polygons_symbols_article.pdf).

6. **Minor — M4's substitution checks 32 coefficient keys, not 31.** `sec_g8_M2alpha.tex:207`; `REPORT_G8.md:17`. Independently collecting the nonzero (function atom, constant atom) keys gives 32. The supplied comparison prints all their zero differences. **Repair:** replace 31 by 32. There is no resulting formula error.

7. **Must fix — distinguish the proved W_3 description from a single iterated integral.** `sec_g8_M2alpha.tex:15,31,35,123`; `main_tex_statements_l1500-1580.tex:12`. The proof establishes

   W_3=Y(c)+Y(-c)-zeta(2)(1-c^2)artanh(c)/(2c).

   This is a rational-function combination of weight-three HPLs and products of zeta(2) with weight-one HPLs. The rational prefactors and constant factor require a stated convention if this is called “an iterated integral of weight three.” The separate Euler-substitution integral is expressly not reproved. **Replacement:** “W_3 is a Q(c)-linear combination of weight-three HPLs and zeta(2) times weight-one HPLs, with letters 0,1,-1.” Alternatively supply the omitted integral construction and its definition of weight.

8. **Minor — make the symbol-space terminology unambiguous.** `sec_g8_M3alpha.tex:48-69,146`. V, the exterior powers, and the displayed tensor map are precise definitions. Here B_2 means the degree-two tensor quotient by shuffle products, identified with Lambda^2 V; it should be called B_2^{symb} to avoid confusing it with the usual dilogarithmic/Bloch group. The formal exterior product is not the wedge of one-variable differential forms. For arbitrary rational arguments, first enlarge the letter space: their letters need not lie in V. In the one-term lemma, projection along extra nonconstant letters is an additional projection, not merely “modulo constants.” It is legitimate because an equality in the larger space implies the projected equality. State that explicitly. The map's normalization is fixed by its displayed formula, so no hidden sign or factor convention is needed for these computations.

9. **Minor — the differential lemma needs a removable-singularity qualification.** `sec_g8_M3alpha.tex:156-175`. The series domain includes x=0, y=0, and y=1, whereas the separated dlog expressions are not defined there term by term. **Repair:** first state the formulas away from xy(1-y)=0, interpret dlog f as df/f, and extend the combined expressions through their removable singularities. The actual paths x=1/2, y=1±c, 0<c<1, avoid these exceptional points. For y=1+c, dlog(1-y)=dlog(-c)=dc/c; no real logarithm of a negative number is required. The four series-derived differential formulas and their signs are otherwise correct.

10. **Minor — the M_4 gap box contradicts its proof.** `sec_g8_M4alpha.tex:35` lists the zeta(3) and pi^4 consequences as untreated, although `:55,59-75` proves them. **Repair:** remove those two items from the untreated list, leaving the generic/special-corner depth-two assertions there. Also read “nonzero for every c” at `sec_g8_M3alpha.tex:121` as a nonzero coefficient vector in the function-letter basis: it is not a nonzero period specialization at every angle, as the square example immediately demonstrates.

## All-n audit: why the calculation is sufficient

The criterion's essential implication is correct: equality of formal word coefficients as rational functions implies equality of evaluated sums at every allowed integer. There is no need to prove that evaluation is injective.

The dictionary gives h_w(m)=gamma epsilon^m m^{-kappa} Z_v(m-1) for m>=1. For m=2n+delta with delta>=-1, its denominator is at least 1. Its shift is j=delta-1>=-2. In a negative shift the removed indices are 2n+t>=2n-1>=1; the remaining upper limit may be zero, where nonempty Z sums are zero and the empty sum is one. Recursion decreases word length, including at j=-2. Positive shifts use only t>=1. Thus the endpoints n=1, j=-2 and delta=-1 are valid, not an eventual-n argument.

Code-to-mathematics correspondence:

- `g8_allN.py:28-53` expands exp(uA)(1+y(exp(u(b-a))-1))^(1-2n) formally to order k. The iterative divisions by r and t produce the exponential factorials; truncation in u cannot alter retained coefficients.
- `:56-83` constructs exp(Delta_*) to total degree k and multiplies by i!m!, exactly as D-sharp requires.
- `:86-134` substitutes the three harmonic-number expressions, rho_j, and -1/[n(2n-1)]. Its stuffle product has the correct positive merge for strict nested sums.
- `g8_zsum.py:83-153` implements the shift and dictionary recursions. Negative powers of signs are reduced by parity, so they do not introduce floating-point arithmetic.
- `g8_allN.py:164-177,281-292` implements the stated Y and the paired HPL sign convention, with delta=parity-e. The extra -Li_2(c^2)/2 contributes exactly -1/(2n^2) to the even coefficient.
- `:202-221,249-275` checks the low Laurent coefficients, including standalone constants. With the smallest actual exponent -3 and nonempty HPLs starting at degree at least one, orders below -4 are automatically zero. No untested small n or low degree is concealed.

The arithmetic uses QQ polynomial rings, the exact rational-function field QQ(n), integers, and Fractions. The word dictionaries are normal forms for the defined formal objects: equality there is decided coefficientwise. This is not a demonstrated complete decision procedure for identities of evaluated sequences. A nonzero formal residual would be inconclusive about sequence equality. The submission neither proves the converse nor needs it; this audit does not assert injectivity or exhibit a kernel relation.

`CODEX_allN_audit.out` prints the actual shift range -2,-1,0,1,2,3 and the distinct monic denominator factors of the collected mu/HPL coefficients and recorded shift coefficients:

    n, n+1, n+1/2, n+3/2, n-1/2.

None has a positive integer root. Intermediate chi/rho expressions may also contain n+t with larger nonnegative t; these are safe as well. Thus cancellation in QQ(n) is not concealing an invalid original evaluation. The same audit prints

    P_{2,1} = -7/12 + s2/3
    P_{3,1} = 35/72 - s2/36 - s3/3
    P_{4,1} = -571/1080 + (103/432)s2 + s3/36 + (4/15)s4 - s2^2/15.

Fresh M1/M2/M5/M8/M9 runs pass every required word comparison, odd comparison, and low-order comparison. The supplied selftest also passes its printed 98 words and 4116 coefficient checks. The finite sample tests against g6_common are controls, not the reason the result holds for all n.

## Classical form: derivatives, constants, branches, substitution

My independent code imports no supplied verification program. It independently checks the five dilogarithm reductions and all ten weight-three derivative/basepoint identities. Three formulas were separately transcribed directly from TeX: H_{1,0,1}, H_{1,0,-1}, H_{-1,1,-1}. Those transcriptions also agree exactly with the data. The other seven were checked from the JSON, compared with the printed formulas by inspection. The fresh supplied M3 run additionally verifies all six weight-two HPL identities and duplication.

| Word w | Required derivative C_w' | Own exact residual | Own basepoint C_w(0+) |
|---|---|---|---|
| 0,0,1 | Li_2(x)/x | 0 | 0 |
| 0,0,-1 | -Li_2(-x)/x | 0 | 0 |
| 0,1,-1 | H_{1,-1}(x)/x | 0 | 0 |
| 0,-1,1 | H_{-1,1}(x)/x | 0 | 0 |
| 1,0,1 | Li_2(x)/(1-x) | 0 | 0 |
| 1,0,-1 | -Li_2(-x)/(1-x) | 0 | 0 |
| 1,-1,1 | H_{-1,1}(x)/(1-x) | 0 | 0 |
| -1,0,1 | Li_2(x)/(1+x) | 0 | 0 |
| -1,0,-1 | -Li_2(-x)/(1+x) | 0 | 0 |
| -1,1,-1 | H_{1,-1}(x)/(1+x) | 0 | 0 |

Every row has domain 0<x<1. In the order used in `sec_g8_M2alpha.tex:132`, the eight Li arguments have ranges (0,1), (-1,0), (0,1), (1/2,1), (1/2,1), (0,1), (-1,0), (0,1). All use the real restriction of the principal branch. All log arguments x,1-x,1+x,2 are positive. In the same order, 1-a(x) is

    1-x, 1+x, x, x/(1+x), (1-x)/2,
    2x/(1+x), 2/(1+x), (1-x)/(1+x).

Thus every logarithm splitting used in differentiation is valid without an imaginary correction. At x->0+, the Li arguments tend to the stated 0,1,-1,1/2 special points; continuity and the given special values suffice. Every log x factor is multiplied by a positive power of log(1±x)=O(x), so these products vanish. No rational pole multiplies a Li_3 in these ten identities, and no unjustified full Taylor series at argument 1 is needed to obtain the constants.

For substitution step M4, the supplied exact comparison checks all 32 coefficient keys against its independently typed classical expression. My separate collection confirms 32 keys and the coefficient 2x+7/4 of zeta(3). My independently typed classical formula also agrees with the defining integral at all four angles below. No uniqueness of this basis or transcendence assumption is needed.

## Proposition 5: independent symbol computation and exact scope

`CODEX_independent_algebra.py` uses determinants for the two successive exterior products. It imports no supplied symbol library. It independently parses the printed 56-term table and proves equality with the weight-four slice of the 74-term JSON. Its output gives 24 nonzero tensors and, for p=x∧y, q=x∧z, r=y∧z,

    delta(Omega) = (0, (c-1)/c, (-c-1)/c)

in the basis (p∧q,p∧r,q∧r). This is exactly the displayed formula. The signs of g(1)=-[1-c] and reversal of word order are correct. The hand shuffle cancellations and Li_4 annihilation lemma are also correct; they give the needed necessary obstruction without a converse theorem.

The one-term exclusion is exhaustive over the class actually named. A constant-coefficient symbol in any enlarged letter space projects to a constant vector (a,b,d). Multiplication by one r(c) yields r(c)(a,b,d). Matching Omega would force a=0 and b/d=(1-c)/(1+c), which is nonconstant. Additional letters, a larger depth, or an arbitrary constant linear combination within G cannot evade this argument. It excludes one common rational prefactor, not an arbitrarily defined G already carrying nonconstant symbol coefficients. It applies to cancellation as well as matching, since an overall minus sign can be absorbed in r.

Independently expanding the generic Li_{2,2} symbol gives the printed 15 tensors and the same five generic cobracket components. Pullback at x=1/2 gives

    delta Li_{2,2}(1/2,1-c) = -p∧r,
    delta Li_{2,2}(1/2,1+c) =  q∧r.

Consequently the displayed prefactors reproduce delta(Omega). With t=[2] retained, the independently computed difference, two-term expression minus Omega, is

    (1/c-1)(x∧y)∧(y∧t)
    +(1-1/c)(x∧y)∧(z∧t)
    +(1+1/c)(x∧z)∧(y∧t)
    -(1+1/c)(x∧z)∧(z∧t).

These four components do not vanish identically. This does not prove that no other corrections could eliminate them.

At the square corner the displayed three representatives do map into the span of [2] and [1+sqrt(2)], so their iterated exterior square collapses. This is a valid linear-algebra check of those representatives. It does not recover constant letters lost earlier or establish a symbol of a numerical period from the constant-free function symbol.

## Theorem 21: ladder, pi^4, and counts

The accepted ladder, specialized to Q_{zeta(3)}=xx'(x+x'), gives [s3]P_{4,n}=-[s2]P_{3,n}. The sign is correct; the n=1 audit prints +1/36 and -1/36 respectively. The fresh all-n run independently confirms this identity.

M9 freshly prints

    [s4]P_{4,n}=(n+3)/[(2n-1)(2n+1)(2n+3)],
    [s2^2]P_{4,n}=-n/[(2n-1)(2n+1)(2n+3)].

The partial fractions in `sec_g8_M4alpha.tex:63-64` are correct. The sums at `:68-70` follow directly by shifting the artanh series; the constant subtractions -1 and -1/c^2-1/3 are correct. They yield the displayed F_1 and F_2, and combining zeta(4)=pi^4/90 with zeta(2)^2=pi^4/36 yields exactly the factor 1/2880 and the displayed bracket. I checked this algebra by hand and reran the exact supplied partial-fraction/series verification. This proves the stated formal decomposition followed by specialization to real zeta values, without asserting uniqueness of a real-number split.

Independent data census:

| Expression | HPL counts by weight | Standalone constant |
|---|---|---|
| R_3 | weight 2: 2; weight 3: 16; weight 4: 56 | none |
| Psi | weight 1: 6; weight 2: 8 | 1/2 |
| R_4 | weight 2: 1; weight 3: 8; weight 4: 27; weight 5: 73 | none |
| Psi_2 | weight 1: 3; weight 2: 6; weight 3: 12 | 1/2 |

## Independent numerics from the definition

`CODEX_independent_numeric.py` imports no supplied code. I substituted g=exp(-t) in the defining integral, expanded its integrand analytically in u through degree four, and quadrature-integrated those coefficient functions. This avoids the supplied Cauchy-integral method. The implementation rationalizes the cancellation at large t. HPL values use a separately written coefficient recurrence from (1-ax)H_w'=H_tail for a=±1 and xH_w'=H_tail for a=0. Paired terms are evaluated at c and -c directly, independently testing the sign convention.

Working precision was 75 decimal digits. The HPL truncation orders were 1200,650,400,300 at the four angles; the transcript also prints comparisons with 50 fewer terms. The following definition values and errors were actually printed by the new code. The last two angles are obtuse.

| Angle | k | M_k from definition | Absolute HPL-form error |
|---|---:|---|---:|
| pi/3 | 2 | -0.04287175497685288253230642505250886346370973319237454663 | 5.1817011e-77 |
| pi/3 | 3 | 0.0343591536761547609039954312098991930921887551172688404 | 4.8362544e-76 |
| pi/3 | 4 | 0.01132305145366556066779367839642295140323962846973444683 | 1.2781529e-75 |
| pi/2 | 2 | -0.02216191108726939942375493773979151497795504632874345621 | 1.2954253e-77 |
| pi/2 | 3 | 0.02277052149178657745673050974232780025018369172921726138 | 5.3112437e-76 |
| pi/2 | 4 | 0.006027993112609122504516108250877605571476734148292802895 | 8.4310596e-76 |
| 2pi/3 | 2 | -0.009567895652956384775719722385651527995008301633351144017 | 1.5545103e-76 |
| 2pi/3 | 3 | 0.01071736293228946541840957629122489585571897573846433576 | 2.612441e-76 |
| 2pi/3 | 4 | 0.002101755152606312397514057423977901270444590652315008808 | 3.5570219e-76 |
| 2acos(1/3) | 2 | -0.004025944948505109500781252223609780955578436078109873169 | 3.0334542e-76 |
| 2acos(1/3) | 3 | 0.004571210718072147109736283222672342747424704096432675895 | 2.6556218e-76 |
| 2acos(1/3) | 4 | 0.0007246280460453527427004111039799214485600052387994256375 | 2.517848e-75 |

For my independently typed classical M_2 expression, the printed errors in the same angle order are 2.4181272e-76, 5.9157755e-76, 4.7498927e-76, 1.8028002e-76. These are numerical discrepancies, not certified enclosures or claims of rigorous 75-digit accuracy. They comfortably meet the requested 30-digit numerical comparison.

**Planted error caught:** adding 1 to the fourth R_4 data coefficient produces residual

    1.54316482491555992778976398406092382371133077

at pi/2. The new check detects it. The fresh supplied mutant suite also detects all its listed coefficient, letter, omitted-term, and exponent mutants.

## Coverage of the paper's claims

The following references are to `main_tex_statements_l1500-1580.tex`. “Not supplied” means these three sections do not establish the assertion; it does not reject separate work that was not provided.

| Paper sentence or clause | Coverage / status |
|---|---|
| Lines 3-7: Beta-integral proof idea and formal zeta-coefficient discussion | Accepted g6 foundation; formal/real split qualification is necessary. |
| Line 9: HPL definitions | Definitions agree with the audited calculus. |
| Lines 12-20: M_2 formula, W_2, and Y | PROVED by g6 + hand criterion + M1/M2. W_3's literal iterated-integral description needs finding 7. |
| Line 24: fitted discovery and proof for every n | Discovery history not checked; the all-n proof is now supplied and valid. |
| Line 27: no zeta(3), and square-corner W_2 | No explicit zeta(3) in the formal coefficient split; square substitution into W_2 is elementary. Identification with T is not checked because its defining equation is outside this excerpt. The zeta split is representation dependent, as the later classical expression illustrates. |
| Line 27: polynomial-in-k weight kernels at four angles, 19-angle numerical claim, six-angle slope J | Not supplied by these sections. |
| Lines 29-33: leading-weight law, Ferrers functions, fitted F_d, eight-angle and high-weight checks | Not supplied by these sections. |
| Line 35: ten classical identities, eight arguments, coefficient 7/4+2c, derivative and constant checks | PROVED by M3/M4 and the independent checks above. |
| Line 35: historical 120-digit comparison | Not reverified to 120 digits here; my independent comparison uses the printed precision above. |
| Line 37: M_3 decomposition and 74 HPL terms | PROVED AFTER the constant-criterion repair, by M5. |
| Line 37: at the square it gives W_4; relation to c_0=c_1=0 and the other theorem | The square coefficient value was numerically checked here; the separate W_4 identification and referenced theorem were not audited. |
| Line 37: (2,2) tensor map, vanishing on products and Li_4 | PROVED for the defined symbol quotient; keep the B_2 convention explicit. |
| Lines 40-45: Omega's 56 terms, letters, formula, nonzero projection | PROVED for the displayed expression. “Of the function” needs finding 3. |
| Line 48: M_3 is not classical as a function, hence needs depth two | NOT PROVED by the submitted expression-level argument; add the missing symbol-independence argument or restrict the claim. |
| Line 48: HPLs are MPLs | Accepted g6 dictionary/convention; no new obstruction here. |
| Line 48: one term r(c)G is insufficient | PROVED, exhaustive in the explicitly named constant-symbol-coefficient class. |
| Line 48: two displayed terms suffice | PROVED AFTER specifying the constant-free cobracket quotient. Full-symbol/function sufficiency is not proved. |
| Line 48: checked by independent programs / exact algebra | Fresh provided and independent implementations agree here; historical report provenance was not audited. |
| Line 48: period-conjectural conclusion for a number at one angle | Not established here; needs the numerical period symbol and a precise conjectural framework, not a function-letter substitution alone. |
| Line 48: square and hexagon letter spans make the chosen exterior test vanish | The linear-algebra observation is correct; vanishing is not a proof of classicality. |
| Line 48: hexagon Lie symbol, 21 S-unit solutions, lattice-reduction formula, level-6 conversion and exact relations, precision claims | Not supplied or audited here. |
| Lines 48,59: lists of non-classical polygon/rational corners and complete rational classification | Not supplied or audited here. |
| Lines 50-59: c=1/3 Lie class, W_2/W_3/W_4 closed forms, level-6 proof and precision claims | Not supplied or audited here, apart from the independently printed M_3 value. |
| Line 59: quadratic search, four extra corners, golden/silver formulas and precision claims | Not supplied or audited here. |
| Lines 64-70: M_4 formula and definition of Psi | PROVED AFTER repairing the coefficient criterion, by g6, M8/M9, ladder and elementary sums. |
| Line 71: 109 terms, weights 2-5, exponent range -3..1 | PROVED by the checked data; paired-term convention required. |
| Line 71: 22 HPL terms of weights 1-3 | Incorrect description of these data; replace by 21 plus constant 1/2. |
| Line 75: explicit formal coefficients, all-n harmonic identities, odd coefficients zero | PROVED by the audited exact computation. |
| Line 75: c^1280 check and historical six-angle/square high-precision checks | The finite coefficient equality follows from the all-n result; those historical numerical runs were not audited. |
| Line 78: zeta(3) ladder consequence | PROVED by the accepted ladder and checked directly in M9. |
| Line 78: pi^4 part is elementary | PROVED by M9 and the partial fractions/sums. |
| Line 78: generic weight-five (3,2) obstruction, explicit pairing, hexagon/1/3/square pairings | Not supplied or audited here; they do not follow merely from the coefficient formula. |
| Line 78: golden-corner obstructions, explicit Li_{3,2} corrections, hexagon classical remainder | Not supplied or audited here. The unidentified products/zeta(5) are correctly acknowledged as unidentified. |
| Line 80: broader ladder | Accepted g6 result; this audit checked its needed M_4 specialization. |

The M2 gap boxes honestly exclude the Euler integral, uniqueness, and independence claims. The M3 gap box honestly excludes the function-invariance step, converse, period and special-corner assertions, and a function-level two-term repair; those exclusions must also govern the paper's wording. The M4 gap box is accurate about missing depth-two consequences and uniqueness, after the contradiction in finding 10 is removed.

A truthful replacement for the paper's strongest unsupported transition is: “The displayed expression has a nonzero constant-free (2,2) tensor projection. No single rational prefactor times a constant-coefficient symbol reproduces it; the displayed two Li_{2,2} terms reproduce this projection. A function-level conclusion requires the expression-independence argument for the associated-graded symbol.”

## Reproducibility and limitations

New runnable checks and their actual transcripts:

- `CODEX_independent_algebra.py` / `.out`: independent tensor/exterior algebra, printed-table/data match, generic and specialized Li_{2,2}, constants retained/dropped, census, exact classical derivatives/basepoints, three direct transcriptions, zeta(3) coefficient and substitution-key count.
- `CODEX_independent_numeric.py` / `.out`: independent definition quadrature, HPL evaluation, classical evaluation, truncation comparison and planted error.
- `CODEX_allN_audit.py` / `.out`: instrumentation of the supplied implementation for denominator factors, shift ranges, n=1 and low coefficients. This check is not an independent reimplementation of the mu-side symbolic engine.
- Fresh supplied-run transcripts: `CODEX_allN2.out`, `CODEX_allN3.out`, `CODEX_allN4.out`, `CODEX_selftest.out`, `CODEX_cor3.out`, `CODEX_symbol.out`, `CODEX_pi4.out`, `CODEX_mutants.out`, `CODEX_numeric.out`. The original numerical run also reproduces its reported agreement.

Run the new checks from this directory using the supplied venv Python, with `PYTHONDONTWRITEBYTECODE=1`. Run the original programs from `scripts/`, directing any new output to new files rather than overwriting the supplied `.out` files.

I did not formally verify SymPy, Python, mpmath, or the programs; certify quadrature with interval arithmetic; prove injectivity of the evaluated nested-sum representation; rederive all accepted g6 sections; reconstruct the original fitting history; audit unavailable earlier referee reports; verify the period conjecture, special-corner classifications, or weight-five symbol computations; rebuild or visually review the PDF. The independent algebra code's seven remaining weight-three inputs are from the supplied JSON rather than a second mechanical TeX parser. These limitations do not replace the specific successful all-n and independent checks recorded above.
