# CX8 independent referee report

**Overall verdict: PROOF COMPLETE, with the write-up corrections in item 7.**

The tower law holds at every admissible depth `0 <= delta <= min(4,2d-4)`, for all degrees and keys: off-diagonal from `k=d`, and stable diagonal from `k=2p+1`. The exceptional diagonal initial level `k=2p` is not on that polynomial in general; E3 proves precisely the blindness needed by its adjacent child. The fitted `e(4,1)` is certified by symbolic identities, not extrapolation. Degree-four towers exhaust the higher part, so WTS Lemmas 2–3 give quadratic degree-three kernels for `k>=5`. The exact kernel checks below extend the same quadratic to `k=4` and establish the required vanishing of imaginary simple poles. With the existing fixed-level WT4 integral certificates, this proves

\[
W_4(k)=W_4''+(k-4)D_1+\binom{k-4}{2}D_2,\qquad k\ge4.
\]

This verdict concerns the depth-four theorem and its weight-four consequence. It does not certify the submission's additional all-depth/master-equation claims. Several of those descriptions need qualification, as listed below. No source file outside this review folder was edited.

## 1. Fitted e(4,1) and completeness of guess-and-verify — OK

Normalize off-diagonal coefficients by `e00=rho_pq(0)`. The candidate uses `e(delta,0)=rho(delta)`, the inherited depth-two coefficients, and

\[
\frac{e(3,1)}{e00}=\frac{i(d-3)}{2d-3},\quad
\frac{e(4,1)}{e00}=\frac{d^3-9d^2+\frac{113}{4}d-\frac i2d+im-29}{(2d-3)(2d-5)},\quad
\frac{e(4,2)}{e00}=\frac{2}{(2d-3)(2d-5)},\qquad m=p.
\]

Here is the complete new identity set; there is no unexamined Newton index at these depths.

| Depth | Newton indices s | Non-adjacent p>=q+2 | Adjacent (p,p-1) |
|---|---|---|---|
| 3 | 0,1 | both exact symbolic zeros | both exact symbolic zeros |
| 4 | 0,1,2 | all three exact symbolic zeros | all three exact symbolic zeros |

`check_details.py` independently verifies all nine non-adjacent identities through depth four as rational functions in the unevaluated symbols `d,m`. This includes the five new cases in the table. The five submitted `prove_adj34.py DELTA S` runs verify the adjacent column symbolically in `p`, retaining both the direct diagonal feed and the fold. The diagonal entries use

\[
\widehat e_{pp}(\delta,s)=(p+1)\sum_{\delta'\le\delta}
\frac{\tau_{1,\delta-\delta'}e_{p+1,p}(\delta',s)}{4p-\delta'-1}.
\]

To feed an adjacent key, the script correctly converts this to `sum_t (-1)^t ehat(delta',s+t)`. Its conjugation substitutes `i -> -i` in the normalized off-diagonal coefficients and uses `conj(tau1_n)=(-1)^n tau2_n`. The adjacent child's normalizing `e00` is real, and its two-generation parent ratio is real, so this conjugation is legitimate. Diagonals need no additional off-diagonal-style N identity with `ehat(delta,0)=rho_pp(delta)`: that equality would be false. Their stable polynomial follows directly from their sole off-diagonal parent, in the shifted Newton basis.

For each fixed integer k, start at degree k with the true initial data and descend in degree. At an off-diagonal initial key the candidate agrees with rho. At a stable key the Newton identities make it satisfy T1. At a diagonal initial key retain its actual rho; when this feeds an adjacent child at `k=2p`, E3 gives the same result as the stable extrapolation. The finite triangular recursion has a unique solution, so the candidate equals the actual tower. This proves equality for arbitrary k without any assumption about the uniqueness of a fit or a recursion at infinite degree.

The fitted `N1..N4` for S derivatives are likewise certified by rational identities and their equality to the library polynomials was explicitly checked. No other fitted all-depth e-family is needed or endorsed. The 32 fit points and checks at degrees 12,13 are discovery/validation evidence, not premises of this proof; I did not rerun those fits.

The assertion of **exact** polynomial degree also holds. Off-diagonal top coefficients at depths 3 and 4 are the nonzero expressions above for every admissible d>=4. On a stable diagonal the corresponding top Newton coefficients, normalized by `rho_pp(0)`, are

\[
\frac{\widehat e_{pp}(3,1)}{\rho_{pp}(0)}=
\frac{(1+i)(2p-3)}{4p-3},\qquad
\frac{\widehat e_{pp}(4,2)}{\rho_{pp}(0)}=
\frac{2(1-i)}{(4p-5)(4p-3)}.
\]

These were verified symbolically in `check_details.py` and are nonzero for p>=2. The lower depths are inherited from the corrected TOW1 proof.

## 2. Initial data rho(3), rho(4) and S derivatives — OK

Cancel the artificial homogenization degree D first. Apart from a constant independent of v, the analytic factor is

\[
A(v)=-\frac{\sqrt2}{2}\,8^d\,
\frac{v^{2d-1}S_m(t(v))(v+i)^{2-2d}}{(v^2-a^2)(v^2-b^2)},\quad
t(v)=v/4+1/2-1/(4v).
\]

Thus `t(i)=t*=(1+i)/2`, `t'(i)=0`, and

\[
t^{(r)}(i)=\frac{(-1)^{r-1}r!}{4i^{r+1}},\quad r\ge2;
\qquad t''(i)=-i/2,\ t'''(i)=3/2,\ t''''(i)=6i.
\]

The direct differentiation check verifies these through r=8. Let `R1=S'/S`, `R2=S''/S`, evaluated at t*, and let `ell_r=(log A)^(r)(i)`. Independent differentiation of the algebraic factor and the chain rule give

\[
\begin{aligned}
\ell_1&=-i(d-2),\\
\ell_2&=(3d-3)/2-iR_1/2,\\
\ell_3&=i(7d-5)/2+3R_1/2,\\
\ell_4&=(-45d+21)/4+6iR_1-3(R_2-R_1^2)/4.
\end{aligned}
\]

With `L_r=ell_r/r!`, the required normalized initial data are

\[
\begin{aligned}
a_3:=\rho(3)/\rho(0)&=L_3+L_1L_2+L_1^3/6,\\
a_4:=\rho(4)/\rho(0)&=L_4+L_1L_3+L_2^2/2+L_1^2L_2/2+L_1^4/24.
\end{aligned}
\]

`check_details.py` verifies that these equal the submitted cached formulas symbolically (also at depths 0,1,2). In particular, the adjacent certificates do not depend on an unchecked cache. Only R1 and R2 enter rho through depth four; R3 and R4 concern the later outlook, but were checked as requested.

For an independent derivation of R3 and R4, put `h=t-t*`, `z*=1-t*`, `j=d-n`, `s=m-r`. Then `q2=2ih+2h^2`. Write `P_r(x)=Stirling(x,x-r)` and `Pi_n(d)=product_{j=1}^n(2d-2j-1)`. The normalized finite S sum gives, modulo h^(l+1),

\[
\frac{S_m(t_*+h)}{S_m(t_*)}=
\sum_{n=0}^{l}\frac{(-2)^n(2ih+2h^2)^n}{\Pi_n(d)}
\sum_{r=0}^{n}P_r(m)P_{n-r}(d-m)t_*^{2(r-n)}z_*^{-2r}
(1+h/t_*)^{2(d-m-n+r)}(1-h/z_*)^{2(m-r)}.
\]

The cancellation used here is
`[binom(1/2,d-n)/binom(1/2,d)] / d_falling_n = (-2)^n/Pi_n(d)`.
The Stirling tails through order four are fixed by their two-term recurrence and zero anchors, including

\[
P_3(x)=x(x-1)(x-2)^2(x-3)^2/48,
\quad
P_4(x)=x(x-1)(x-2)(x-3)(x-4)(15x^3-150x^2+485x-502)/5760.
\]

Extracting `l! [h^l]` in this independent implementation gives `N_l/Pi_l` with exact zero rational residual for l=3,4 and symbolic d,m. The submitted `prove_S12.py` and both `prove_S34.py` runs independently return zero as well. Endpoints m=0,d are handled after polynomial cancellation, without evaluating uncancelled 0/0 ratios. The actual finite S polynomials were additionally checked exactly for every d=1..4, m=0..d, l=1..4, covering the small-degree cases where the generic sum reaches j<=0.

Numerical spot checks use the finite S polynomial, Horner evaluation and 80-digit mpmath differentiation, not the candidate R formulas as a definition of S. Samples `(d,m)` are `(4,0),(4,2),(4,3),(4,4),(5,3),(7,4),(11,8)`. All S derivatives through four agree with the closed forms. For the six canonical samples, rho(3) and rho(4) from direct derivatives of A agree with exact initial tower coefficients loaded through `towlib` from WTS data. Their normalized formulas were also compared to those tower coefficients in exact arithmetic. Maximum absolute error divided by `max(1,abs(expected))` across the numerical checks: **8.5690111e-79**. See `details.out`.

## 3. Symbolic coverage, scripts and domains — OK

The submitted S certificates use symbolic d,m; the adjacent and blindness certificates use symbolic p; the top proof uses symbolic d,s (with p in its parent-ratio brackets). Loops over derivative splits, fixed depths or the finite list of required identities do not specialize those degree/key parameters. The top proof cancels a general numerator and proves only the top equations; that is sufficient here because the value/slope equations also pass.

Domains must be imposed before invoking these identities:

- d>=2 and `0<=delta<=min(4,2d-4)`. Degree two has only depth zero; degree three only depths zero through two.
- New non-adjacent identities at depths 3,4 require d>=4 and p>=q+2>=2, q>=0.
- New adjacent identities and E3 at depths 3,4 require p>=3, since the child degree is 2p-1. At p=2 the child has only depths 0,1,2; at p=1 it has no higher part.
- The diagonal itself first admits depths 3,4 at p=2, with stable levels starting at k=5.

In these domains, all actual parent/feed denominators are positive integers and all Pi factors are nonzero odd integers. The initial normalizing rho(0) is nonzero by TOW1's binomial formula. The fold's generalized binomial is a falling factorial divided by a fixed factorial and equals the ordinary binomial in the admissible domain. Reducing powers of i modulo four is valid because p is an integer. Thus clearing denominators proves the needed identities everywhere in the domain, not just generically away from an overlooked integer root.

Fresh submitted-script execution is logged separately in `proofs_progress.out` and the corresponding `.out` files. The three value/slope scripts are **step_48.py, step_49.py, step_52.py**, not the pointers 51,52,55 in the write-up. All fourteen invocations completed successfully; their actual zero remainders/true zero-numerator results were inspected. The two value identities took 34.88 and 274.64 seconds including interpreter startup; the slope took 1.79 seconds. The slowest adjacent identity was `(4,0)`, 81.84 seconds. The independent simplified implementation checks the same identities without the large radical expressions.

## 4. Blindness and first diagonal levels — OK

`prove_E3b34.py 3` and `prove_E3b34.py 4` freshly return exact zero numerators. These use the actual discrepancy

\[
\Delta_{pp}(\delta)=\rho_{pp}(\delta)-\sum_s(-1)^s\widehat e_{pp}(\delta,s),
\]

which extrapolates the stable diagonal polynomial one level backwards. The code normalizes by the purely imaginary `rho_pp(0)`; conjugating the *absolute* discrepancy therefore needs the additional minus sign in `Dc`. The submitted code has that sign. Both the direct DL2 feed and the folded DL1 feed are retained, with no double counting of the corrected TOW1 B block.

The new referee harness checks **129 exact tower cells** for diagonal/adjacent keys with p=1..4, every admissible depth <=4, and the first levels through k=2p+3. It also compares the initial and extrapolated diagonal feeds at k=2p for p=2,3,4 and every admissible child depth <=4. This includes `(2,2)` at k=4,5,6,7 and `(3,2)` at k=5,6,7,8. All pass. At p=1 the diagonal has only depth zero and its child has no higher part; no singular formal depth-three/four expression is evaluated there.

## 5. Weight four: the k=4 boundary and imaginary poles — OK

For k>=5, WTS Lemmas 2–3 apply to the complete degree-four towers (orders 2 through 6, exactly depths 0 through 4). Integration of their higher parts, normalization at b, multiplication by DL1/DL2, folding and taking the simple part are fixed maps independent of k. Conjugation in these maps preserves polynomiality in the real integer k. Consequently every degree-three residue, including each imaginary simple-pole coefficient, has degree at most two.

`check_kernels.py` loads **kern_4.pkl, kern_5.pkl, kern_6.pkl, kern_7.pkl** and compares every exact partial-fraction component for both `(3,0)` and `(2,1)`. It verifies

\[
K^{(4)}_{pq}=3K^{(5)}_{pq}-3K^{(6)}_{pq}+K^{(7)}_{pq}.
\]

This is the quadratic at k=4 extrapolated from 5,6,7, including zero/missing residues and the constant component. Both imaginary lists are zero at k=4,5,6,7. In particular the three levels 5,6,7 force each imaginary simple-pole polynomial to vanish identically for k>=5. This supplies the extra finite check needed by Lemma 3, rather than silently assuming absence of those poles at all k.

For reproducibility the exact differences about k=4 are:

| Kernel | K4 | first difference | second difference |
|---|---|---|---|
| (3,0) | `(1/24-sqrt(2)/8)/v - 1/[8(v-a)] + 5/[24(v+a)] - 1/[6(v+b)]` | `-1/(12v) - 1/[2(v-a)] + 5/[6(v+a)] - 1/[6(v+b)]` | `-1/[3(v-a)] + 1/[3(v+a)]` |
| (2,1) | `(7/16+sqrt(2)/8)/v + 5/[16(v-a)] - 7/[8(v+a)] - 5/[16(v+b)]` | `1/v + 3/[4(v-a)] - 2/(v+a) - 3/[4(v+b)]` | `1/(2v) + 1/[2(v-a)] - 1/(v+a) - 1/[2(v+b)]` |

The first/second differences also freshly agree exactly with `WT4_weight4_all_k/deltas_deg3.json`, the targets used in the existing D1/D2 integral certificates. Thus the polynomial is the one with the stated constants, not merely some unspecified quadratic. Integrating it against the fixed log monomials yields the claimed W4 formula for every k>=4. The closed-form evaluations of those fixed integrals and W4(4) are inherited established WT4/M4 results; their large relation certificates were inspected but **not rerun** in this audit.

## 6. Fresh full-reducer point — GAP in optional validation (time limit)

A bounded attempt at k=22 used `WT2_weight2_all_k/buildintegrand.py` and `fastred.py`, attempting to build the integrand afresh without loading a saved integrand, tower, or level-22 kernel. Only the expected quadratic uses saved kernels at levels 5,6,7. The process hit its **360-second limit during integrand construction**, before the fresh-integrand completion message or any fresh-level comparison. `fresh22.out` records the partial output and explicit timeout; `fresh_progress.out` records the same status. No fresh k=22 kernel, zero residual, or independent new-level data point is claimed.

This is repairable by rerunning the existing local check with a longer compute allowance or a faster integrand builder; no mathematical repair to the theorem is needed. The user made this check conditional on completing within a few minutes, so I stopped the attempt at that limit. The overall PROOF COMPLETE verdict rests on the symbolic identities and exact boundary kernels, not on this optional validation.

## 7. Required write-up corrections — OK, specified

1. **Normalize rho correctly.** Replace the ambiguous `rho(delta)/(1+[p>q]) * e00-factor` by `a_delta=rho_pq(delta)/rho_pq(0)=A^(delta)(i)/(delta! A(i))`. The displayed a3,a4 formulas are correct with this normalization.
2. **Normalize Delta correctly.** The displayed factored Delta(3), Delta(4) are `Delta_pp(delta)/rho_pp(0)`, not the absolute discrepancy defined immediately above them. Both normalized formulas were verified symbolically. State the conjugation sign caused by imaginary rho_pp(0).
3. **Separate off-diagonal e from stable diagonal ehat.** The universal e formulas, including e(4,1)/e00 and the zigzag tops, are off-diagonal formulas. Stable diagonal coefficients come from the displayed ehat map and shifted basis. Do not set their constant Newton coefficient equal to their exceptional initial rho. List the domains in item 3.
4. **Make the finite certificate list explicit.** Cite steps 48,49,52 for the two non-adjacent values and the slope. Cite all five adjacent `(delta,s)` cases and both E3 arguments. Explain the finite descending-induction uniqueness step, so the fit cannot be mistaken for proof.
5. **Correct zigzag generating functions.** For positive Euler zigzag numbers, `sec x + tan x = sum A_n x^n/n!`. The coefficients of `sech x + tanh x` alternate signs; the implementation takes absolute values. Its current verbal equality without signs is false. In result.md's first odd-top formula, normalize the left side by e00 as well.
6. **Qualify all-depth claims.** Change the title's “top-coefficient law PROVED at all depths” to “non-adjacent top equations proved for an arbitrary degree-independent numerator.” Those equations cannot determine the zigzag constants, and adjacent top equations and blindness beyond the proved depths remain open. “Value identities are the sole missing ingredient” contradicts these other open conditions. None is needed for the reviewed theorem.
7. **Fix master-equation status and indexing.** The promised full generating-equation derivation is not in proof_all_depths.md section 5 as claimed. For the displayed `(1+x) Uhat = rho + x Lambda(Uhat)`, the x^0 coefficient is `Uhat_0=rho`; the value equation `e_0+e_1=Lambda(e_0)` is the x^1 coefficient. A generating identity alone does not prove the required Newton termination. Supply definitions, convergence/formal interpretation, diagonal shifts and blindness, or retain this as an unreviewed outlook. Do not use it as a premise for depth four.
8. **Reconcile status text.** Section 5 and library comments still say the S'''/S and S''''/S symbolic proof is pending although prove_S34.py supplies it. The result also alternates between “rho through depth eight proved” and “Bell assembly unfinished”; distinguish a finite symbolic formula supplied by proved ingredients from an assembled/checked closed form. The advertised single verify_tow2.py is absent; per-lemma certificates suffice mathematically, but the deliverable statement should say so.
9. **State the weight-four boundary argument explicitly.** Lemma 3 supplies k>=5 only. Add `K4=3K5-3K6+K7` for both degree-three keys and the imaginary-pole zeros at 5,6,7. With this report the paper's weight-four range may be changed to “every k>=4,” retaining the existing exact fixed-integral certificates.
10. **Clarify WTS linearity.** With the fold, the map on +i tower coordinates involves their complex conjugates. Call it a fixed real-linear map, or include conjugate coordinates explicitly, rather than calling it complex-linear. Its polynomial-degree conclusion is unaffected. Carry forward CX7's corrected combined B indicator and Taylor-conjugation identity.

## Reproducibility and verification limits

All commands use `<project>/.venv/bin/python -B` and all outputs are local. `run_proofs.py` runs each submitted certificate in a subprocess, records its output and return status, and removes only the historical external `os.chdir(work)` calls in the step scripts. No mathematical expression is replaced. The independent `check_details.py` and `check_kernels.py` leave runnable assertions for the added checks. No dependency was installed, and no bytecode is written outside this folder.

The first boundary harness run failed because Python divided two ordinary integers in its new e(2,1) helper and produced a float. Sympifying the inputs fixed this referee-only error; the complete corrected harness passes. No submitted symbolic identity was changed. The historical 1348-cell and 72-cell totals were not rerun as such and are not claimed as fresh results here. No fresh TOWN1 large pickle check or full WT4 relation-certificate rerun is claimed.
