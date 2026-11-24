**Overall verdict: PROOF COMPLETE, with the corrections below, for the depth-eight tower law and the weight-six law for every integer k>=6.** Task D's broader equivalence has a repairable gap and is excluded from this verdict. There is an incorrect printed depth-seven top formula; the computed family and all executed certificates use the correct formula.

CX10 independent referee report. Scope: admissible depths 0<=delta<=min(8,2d-4), off-diagonal levels k>=d, and stable diagonal levels k>=2p+1. The exceptional diagonal initial value at k=2p is retained separately. All writes are confined to this review folder. I read TOW4 `result.md` and `proof_depth8.md` fully, the seven requested Python files, and the relevant inherited CX9/WTS/WT6S arguments. Only checks actually executed are reported as fresh passes.

1. **OK — forced recursion and closing conditions.** The submitted `forced78.py` was executed with only its output pickle redirected into this folder. All nine closing conditions delta=0..8 pass, including F(7,4)=0 and F(8,5)=0 (`forced.log`). The generated 25-entry string cache equals the submitted cache exactly. Symbolic sanity checks confirm F equals TOW3's e at all 16 indices through depth six, and every top through depth eight equals the zigzag formula (`initial.log`). The recursion is triangular in s: its same-depth parent F(delta,s) is already known when constructing F(delta,s+1). Initial F(delta,0)=a_delta and the closing equation give every non-adjacent Newton index; above top both sides vanish. Adjacent equations and E3 remain separate requirements. Finite descending induction in degree at each fixed k then applies, retaining actual diagonal initial values separately.
2. **OK — fast adjacent proof.** All nine submitted fast cases rerun successfully: (7,0..3),(8,0..4); exact polynomial numerator zero, exit 0 (`remaining.log`, plus `adj_8_0.log`). Original TOW3 slow `prove_adj56.py` and the fast method both pass (5,1),(5,2),(6,2),(6,3), including a non-top depth-six case. The copied slow source is byte-identical to TOW3's original. `check_method.py` additionally proves equality of the lhs and every individual slow/fast feed term at those four cases in Q(i)(p), independently obtains zero field residuals, and checks that adding 1 is detected as nonzero (`method.log`). The denominator method is sound: cancel each fraction n_j/d_j, factor each nonzero denominator as c_j times monic factors, let D be their polynomial lcm, and compute sum n_j(D/d_j). D is nonzero in Q(i)[p], each division has an asserted zero remainder, and contents c_j are retained. Thus a zero polynomial is a rational identity, not multiplication by zero. No numerical evaluations replace the symbolic p variable.
3. **OK — initial data.** `check_initial.py` passes (`initial.log`): L-cache minus the proved tow2lib.L(r) cancels symbolically to zero for every r=1..8. Both tau definitions agree exactly with the tables through n=8 and obey conjugation symmetry; tau1(8)=(-169+41i)/128 and tau2(8)=(-169-41i)/128. An independently implemented Bell coefficient recurrence n*a_n=sum_{r=1}^n r*L_r*a_{n-r} agrees with the partition assembly and F(delta,0) through depth eight. Only R1..R4 are needed since t'(i)=0; their proofs are inherited from CX9. Fresh exact TOWN1 checks pass at all depths 0..8, 21 initial samples per depth (189 total), degrees 6,7,8,12,20,40,52 and p in {d,d-1,ceil(d/2)}. Each depth includes six diagonal initial values and one adjacent key. These checks use rho(0)*a_delta, with no extra off-diagonal factor. They are fresh spot checks, not a rerun of the claimed 717-per-depth total.
4. **OK, with domain correction — E3 and domains.** Both original E3 certificates pass freshly at delta=7 and delta=8: exact zero numerator, exit 0 (`remaining.log`). The direct DL2 and conjugated fold feeds are both included. e00_pp is anti-real, so physical conjugation in its units is Dc=-Delta(conj=True); the adjacent child's e00 is real. Admissible new depths require d>=6. Adjacent keys have d=2p-1, so their new depths require p>=4. Stable diagonal d=2p first permits p=3. The formal p-identities may extend outside these tower domains; that does not create a depth-7/8 coefficient for the degree-five key (3,2).
5. **OK — k = 6 and boundary kernels.** Independent `check_kernels.py` passes (`kernels.log`): all three full partial-fraction dictionaries obey K6=5K7-10K8+10K9-5K10+K11, including 48/48 real-residue QSI components. Every +/-i coefficient and constant vanishes separately for all three keys at each k=6..11 (18 full-kernel zero checks). The complete difference kernels Delta^j K6, j=0..4, equal WT6S's five coefficient kernels for each key (15 comparisons); integrating gives E_{5,j+1}. All 24 rational/sqrt2 component polynomials in `deg5_residue_polynomials.json` equal the Newton expansion from `e5m_binomial_coeffs.json` symbolically, and reproduce the full saved kernels at k=6..11. Thus the boundary and identification needed after WTS Lemma 3 are supplied, independently of TOW4's verifier.
6. **GAP as stated; repairable — Task D (not needed for the theorem).** The finite telescoping equivalence is correct: with P=floor(delta/2), U_dp=sum_s(-1)^s F(dp,s), and C_delta=Lambda_delta(P)-F(delta,P), one has C_delta=(-1)^P(feed(U)-a_delta). This was rerun exactly at delta=7,8 (`details.log`). Evaluating x=-1 is legitimate for these finite Newton polynomials, coefficient by coefficient in depth. However, the asserted equivalence to the full all-depth tower theorem omits the adjacent identities and E3 at all depths. Equivalence with a separately prescribed zigzag-descending value family also needs top matching, not merely the top equations. A fresh symbolic example changes only a8 by d^4/Pi4: its forced changes are Delta_d^s(d^4)/Pi4, so closing still holds (fifth difference zero), but the top changes by 24/Pi4. Top equations allow that change. Thus top equations alone do not justify the claimed ascending/descending identification. At actual depths 7,8 top matching, adjacency, and E3 are separately verified, so the theorem is unaffected. Repair Task D by stating the finite non-adjacent telescoping reduction, and explicitly retaining these additional premises for the full theorem. Also fix indexing: for (1+x)Uhat=rho+x Lambda(Uhat), x^0 says F0=rho; x^1 is the value/Newton s=0 equation; x^(s+1) is N(delta,s), with the top at x^(P+1). TOW4's claim that x^0 is the value equation reverses CX8's actual correction.
7. **ERRORS in the exposition, repairs specified.** Freshly rebuilt F(7,3) is **5i(d-5)/Pi_3**, not the displayed -2i(d-5)/Pi_3. The cache and all executed certificates use the correct value. `check_details.py` mechanically parses and checks all seven displayed F formulas: this is the sole incorrect one. Both saved a7/a8 closed forms and the printed a7 numerator equal the Bell formula exactly (`details.log`). For delta=7,8, adjacent admissibility 4p-6>=delta means **p>=4**, not the text's p>=3. The latter bound belongs to depths 5,6. Feed denominators 2d-delta'-1 are nonzero on the domain but are not all odd (delta' odd gives even denominators). Final wording recommendations appear below.

**Additional details supporting items 1–4.** All e/F formulas in the libraries are normalized by rho_pq(0), not absolute tower coefficients. The off-diagonal factor two is already inside that normalization. In the adjacent script the two-generation ratios are

```
rA0 = i(4p-3)/(p+1)
rB0 = (4p-1)(4p-3)/(p(p+1)).
```

Its ehat is normalized by the real adjacent child's rho(0). The E3 script instead normalizes by the diagonal rho_pp(0) and uses `rB0d=2i(4p-1)/(p+1)`. This accounts for its conjugation minus sign. The tau swap `conj(tau1_n)=(-1)^n tau2_n` has been checked exactly through n=8.

The stable diagonal basis is binom(k-2p-1,s). To feed the adjacent child's Newton equation in the basis binom(k-2p,s), the correct coefficient is

```
Ediag(delta,s) = sum_{t>=0} (-1)^t ehat(delta,s+t).
```

This follows by expanding binom(z-1,j) in binom(z,s); the sum terminates. The submitted fast and slow implementations both use it. In the original fold formula put mu=4p-4-t and order=4p-4-delta. Its factors reduce as

```
(-i)^mu = (-i)^(-4-t),   i^(order-2) = i^(-6-delta),
binom(mu-2,mu-order) = binom(4p-6-t,delta-t).
```

These are valid because p is an integer and i^(4p)=(-i)^(4p)=1. The reduced expressions then define rational functions in symbolic p; no arbitrary complex-p exponent identity is assumed. The falling-factorial binomial has fixed nonnegative lower index delta-t. Thus the fast method preserves the entire inherited fold, including its conjugation and signs.

There are no omitted Newton equations: depths seven and eight require s=0..3 and s=0..4 respectively; all larger s vanish by the support convention. For a fixed finite k, start at degree k with actual initial data and descend in degree. Non-adjacent and adjacent Newton identities give T1 at stable off-diagonal keys. The single-parent ehat construction gives it at stable diagonals. At a diagonal initial key retain actual rho, using E3 only to justify replacing its feed to its adjacent child by stable extrapolation. This finite triangular recursion is unique, so the candidate is the actual tower. No uniqueness of a fit or induction from infinite degree is used.

The degree bounds are exact at the new depths. The off-diagonal tops are `5i(d-5)/Pi3` and `16/Pi4`, nonzero for d>=6. Independently simplifying the stable diagonal top coefficients gives (`method.log`)

```
ehat(7,3)/rho_pp(0) = 5(1+i)(2p-5)/((4p-7)(4p-5)(4p-3))
ehat(8,4)/rho_pp(0) = 16(1-i)/((4p-9)(4p-7)(4p-5)(4p-3)).
```

Both are nonzero for diagonal p>=3. Every feed denominator 2d-delta'-1 is at least 3 on the admissible domain; Pi factors, p, and p+1 are likewise nonzero. Clearing denominators therefore loses no admissible key.

**Full exact tower confirmation — OK.** The original `check_full78.py` was freshly run through the local read-only wrapper, exit 0. It reports **49,953 cells, 0 bad**, including 1,596 stable diagonal cells, plus **222 separate diagonal initial cells, 0 bad** (`full.log`). This exhausts that script's stated range: d<=52, delta<=min(8,2d-4), k<=min(d+8,52), using the supplied TOWN1 pickle. The distinction between initial and stable diagonal data is preserved. This is confirmation, not the all-degree/all-k proof premise.

**Weight-six implication of item 5.** Degree-six towers exhaust depths 0..8. Their diagonal (3,3) is stable only from k=7, so WTS Lemma 3 supplies degree-five kernel residues of degree at most four in k for k>=7. Each imaginary simple-pole coefficient has the same degree bound; its five exact zeros at k=7..11 force identical vanishing. Constant parts are also a fixed linear coordinate of the parent-to-kernel map (and the rational expressions are proper); the saved constants vanish at all six checked levels. The complete k=6 equality places that level on the same polynomial. The fixed real-linear tower-to-kernel map, followed by the fixed integral map, preserves polynomial degree. The fifteen WT6S coefficient-kernel identifications then give

```
W6(k) = sum_{m=1}^5 binom(k-6,m-1) E_{5,m},   every integer k>=6.
```

Here E_{5,m} is the exact fixed-integral combination in WT6S's table. This verdict does not upgrade any separate numerical closed-form relation or symbol-level classicality claim. I did not check k=5's claimed 15/48 discrepancy; it is outside the theorem and not needed.

**Task D precision.** The perturbation in item 6 is a check on the proposed general logical implication, not a change to the actual proved Bell initial data. It shows why top equations alone cannot identify the ascending family with a separately prescribed zigzag-descending family. The finite-depth top matching is independently established here. One may evaluate x=-1 in the truncated polynomial and use the telescoping residual; one may not substitute x=-1 into an arbitrary infinite formal x-series without first proving suitable termination or another evaluation rule. An all-depth tower equivalence must additionally include adjacency and E3 at every depth.

**Required write-up corrections (item 7, final).**

- Replace the printed F(7,3) by **5i(d-5)/Pi3** in both `result.md` and `proof_depth8.md`; retain A4=5. The other six displayed F formulas pass mechanical symbolic comparison. Both a7/a8 pickle formulas and the printed a7 numerator also pass.
- State adjacent **p>=4** for the new depths. Retain diagonal p>=3 and the defining inequality delta<=2d-4. Replace “all feed/Pi denominators are nonzero odd integers” by “all relevant denominators are nonzero”; feed denominators may be even.
- Replace “at k=2p E3 gives the stable extrapolation” by “the actual initial diagonal value is retained, and E3 makes its adjacent child blind to its difference from stable extrapolation.” E3 does not put the initial diagonal coefficient on its stable polynomial.
- Qualify Task D and correct its x-indexing as in item 6. The phrase “closing is the sole missing condition” applies to the non-adjacent system, not the full folded system without further hypotheses.
- The written fast-method validation lists depth-seven/eight cases, not TOW3 cases. Cite the four actual depth-five/six slow/fast comparisons supplied here, especially (6,2).
- Make the k=6 full-kernel certificate, the five k=7..11 imaginary-pole zeros, and the WT6S difference-kernel identification explicit. A quartic obtained by fitting alone does not establish the all-k result.
- `verify_tow4.py` writes the source cache through its forced subprocess and infers several passes from stdout substrings without consistently requiring a successful return code; its full-law match is also weaker than checking both counts. For robust reproducibility, redirect generated output and require exit 0 plus each explicit result. The underlying certificates have been run and inspected here, so this is not a proof gap.
- The full a8 numerator is elided in the proof text despite `result.md` saying full numerators are printed there. Point accurately to the checked pickle or print the expression. The copied tow2lib docstrings still describe R3/R4 as pending although CX9 supplies their proofs; update those status comments if maintaining the copied library.
- `paper/squircle_note/main.tex`, `prop:w6k`, currently says 6<=k<=52 and labels CX10 pending. Its structural range may become every integer k>=6 with the theorem and boundary certificate above. Preserve the independent status of the numerical/symbol statements. No paper file was edited here.

**Fresh certificate inventory.** All arithmetic in these checks is exact; I did not rerun the optional 40-digit tau numeric check.

| Check | Fresh result | Retained evidence |
|---|---|---|
| Forced rebuild, depths 0..8 | 9/9 closings; identical 25-entry source cache | `forced.log`, local `forced_F78.pkl` |
| Adjacent depth 7 | (7,0..3), all exact numerator zero | `remaining.log` |
| Adjacent depth 8 | (8,0..4), all exact numerator zero | `adj_8_0.log`, `remaining.log` |
| Original TOW3 slow and TOW4 fast | (5,1),(5,2),(6,2),(6,3), all pass both ways | `remaining.log` |
| Individual slow/fast feeds in Q(i)(p) | all four cases equal term by term; zero residuals | `method.log` |
| E3 at 7,8 | both exact numerator zero | `remaining.log` |
| L-cache, tau, Bell, lower-depth/top bridges | all symbolic comparisons pass | `initial.log` |
| Initial tower samples, depths 0..8 | 189/189 exact matches | `initial.log` |
| Original full law checker | 49,953 + 222 cells, 0 bad | `full.log` |
| k=6 full PF boundary | 48/48 real QSI components; full kernels equal | `kernels.log` |
| Imaginary poles and constants k=6..11 | all three kernels at all six levels zero | `kernels.log` |
| WT6S identification | 15 difference kernels; 24 polynomial identities | `kernels.log` |
| Printed formulas, Task D residuals and premise example | expected unique printed top error; both master identities zero | `details.log` |

**Reproduction and limits.** Use the specified venv Python with `-B`, cwd this folder. `run_certificate.py` executes the original submitted source, redirecting only forced78's pickle output into this review directory; it adds no mathematical changes. Reproduce with:

```sh
PY='<project>/.venv/bin/python'
"$PY" -B run_certificate.py forced78.py
"$PY" -B run_certificate.py prove_adj78_fast.py 8 0
"$PY" -B run_remaining.py
"$PY" -B check_initial.py
"$PY" -B check_method.py
"$PY" -B check_kernels.py
"$PY" -B check_details.py
"$PY" -B run_certificate.py check_full78.py
```

The symbolic runs were performed sequentially, directly, without nice/ps. I did not run the monolithic `verify_tow4.py`, regenerate TOWN1/WTS data, rerun all CX8/CX9 base certificates, redo WT6S integration/symbol calculations, or prove any all-depth extension. The corrected CX9 depth<=6 theorem and written WTS tower-to-kernel reduction are inherited premises. No source, data, dependency, memory, or paper outside this folder was modified.
