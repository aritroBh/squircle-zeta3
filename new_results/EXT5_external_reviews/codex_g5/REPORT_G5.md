# REPORT_G5 (decay law group): Theorems 12 and 13 of the full version

Files (all in `paper/ieee/proofs/`): `sec_g5_decay.tex` (Thm 12), `sec_g5_varpi1.tex` (Thm 13), `macros_g5.tex`, `wrap_g5.tex`, `wrap_g5.pdf`,
`scripts/g5_consts.py`, `g5_varpi1_algebra.py`, `g5_checks.py` (each with `.out`). Nothing else written; the temp files I made were in the
scratchpad only (no tmp folder in the proofs directory).

## Build
`tectonic -C --keep-logs wrap_g5.tex`: 18 pages, "??" absent in the PDF text, no undefined reference in the final PDF, one overfull box of 0.25 pt
(sec_g5_varpi1.tex:50), a few underfull boxes. Pages: introduction 1, Thm 12 section about 9 (pp. 2-10), Thm 13 section about 8 (pp. 10-18).
(The log file also lists the first-pass "undefined" warnings of tectonic's multi-pass run; the final PDF has none.)

## What is proved by hand
Thm 12 section: exact representation of M_k (logit form, kernel identity by contour indentation, differentiation under the integrals, rotation
of the main term to the line s = 1/2 + r*omega); bound of the correction kernel W_k by rotating a finite sum of integrals; the integral lemma F(j,q;lambda);
the rest bound (Prop. rest) with the crude form C_* e^{-2 sqrt N} and a refined form c_rest N^{-2} e^{-theta}; assembly of 1.21e9 as in the paper
(10,600,000 * 114.096 + main term).
Thm 13 section: sigma-variable form of the main term, modulus lemma (convexity of mu, bound 0 <= Re R <= gamma_0 y e^{-3 sigma}), Watson bound for H with
explicit constant, region I crude bound (exact Bessel K_{3/2} integral), region II saddle point expansion with explicit remainder (Cauchy estimates on a strip,
Gaussian moments, tails by convexity), identification of the model with T (= -i sqrt(pi)/N e^{-z}(1+varpi_1/sqrt N)).

## By script
`g5_consts.py` (Arb enclosures, 200 bits, under 1 s): B_1 = 70.0985853351..., B_2 = 170.6577698580..., C_* = 10,578,907.8235... < 10,600,000;
(3/(e delta))^3 = 114.0963...; c_rest = 114,536.10 < 120,000; c_1 = 5.0886 < 5.1; c_I = 70,217.69 < 71,000; E_II(10^4) = 7912.19 < 8000;
C_13 = 71,000 + 120,000 = 191,000 <= 1,000,668; Thm 12 total 1,209,428,681.35 < 1.21e9 (margin 571,318).
`g5_varpi1_algebra.py` (exact SymPy): Lambda_N = -i sqrt(pi)/N; varpi_1 sqrt(2/pi) = eta_1 + (nu_0 - 1/8)(1+i)/(sqrt2 pi); closed forms of nu_0, eta_1; varpi_1 to 60 digits.
`g5_checks.py` (mpmath floating point, about 4 minutes; numerical only): J versus T, and M_k from its definition by a Cauchy integral, N = 100..420.
Quoted maxima: N^2 e^theta |M_main + X| = 1.7511724 (N=116), |M_k + X| = 1.743787 (N=158): equal to the notes' 1.751172383 and 1.743787030.
M_0, M_1 reproduce the paper's values to 20 digits.

## Gapnote boxes (all are "Computer-assisted step", finite Arb evaluations)
1. sec_g5_decay.tex:203 constants B_1, B_2, C_*, c_rest (reason: numerical evaluation of closed expressions).
2. sec_g5_varpi1.tex:99 c_1(100) and S_1(100) > 1.
3. sec_g5_varpi1.tex:244 E_II(10^4) and the intermediate constants.
No step is left unproved in the text; the monotonicity in N of every bound is argued by hand (log-derivatives) and then evaluated at the left endpoint.

## New arguments, need a referee
Essentially the whole main-term proof is ours and differs from the notes:
- sec_g5_varpi1.tex:40-95 sigma-variable representation, modulus lemma with convexity and the inequality -eps^3/6 <= log(1+eps+eps^2/2) - eps <= 0; Watson bound with C = 3/2 u^2 (lines 80-93).
- sec_g5_varpi1.tex:95-137 region I (100 <= N <= 10^4): |J| <= c_1/N e^{-theta}; this alone gives the claim there because 10^6 (here 71,000) is far above N * (true size).
- sec_g5_varpi1.tex:139-245 region II (N >= 10^4): factorisation f = e^{-z} A e^{-z(cosh sigma - 1)}, Cauchy bounds on the strip, window remainder, tails, identification; monotonicity claims at lines 237-243.
- sec_g5_decay.tex:228-232 monotonicity of the normalised rest pieces (log-derivative bracket) and the refined constant c_rest, used for Thm 13 instead of the notes' rest bound with exponent sqrt(2 pi).
- sec_g5_decay.tex:150-168 the notes' hand step for A_k, rewritten (contour rotation per term); it follows the notes, but the write-up is ours.
- sec_g5_decay.tex:(proof of Thm 12, last paragraph) the explicit O(k^{-1/2}) constant for formula (Mk).

## Disagreements between main.tex and the notes (no statement changed)
1. main.tex proves Thm 12 with a main-term bound 254.23 N^{-3/2}e^{-sqrt(pi N)} from 85,210 interval cells (notes D13s/D13u, "supplied input, not rerun"). Not reproduced; we use our Prop. main (71,000 N^{-2} e^{-theta}), which gives 7,100 N^{-3/2}.
2. main.tex proves Thm 13 via 10^6 + 667.357 (rest bound 1030 N^{-1/2} e^{-sqrt(2 pi N)} from the notes RST1). We do not prove the 1030 bound (it is outside the statements of Thms 12 and 13) and do not use it; our split is 71,000 + 120,000 = 191,000, which implies C = 1,000,668. The notes' own DCY2 corollary with the weaker rest bound is 10^6 + 1.85e10, which does NOT give 1,000,668; only the RST1 route does.
3. The 1.76 for 100 <= N <= 400 in main.tex is consistent with our check (1.7512).
4. Notes' DCY2 original text contained a sign error (|M - X| instead of |M + X|, "constant 1409"); repaired in DCY2F; irrelevant to our proof.
5. The statement of Thm 12 ends "In particular (Mk) holds"; the label is not available in the wrapper, so the formula is copied (noted in the text).
6. CX6 text says |s+q| >= (rho+q/sqrt2)/2 in one place; the sharp version (/1) holds and is what we prove (Lemma F).
7. Core paper does not state Thms 12, 13, so headings say so.

## Machine steps and runtimes
g5_consts.py: M1, under 1 s. g5_varpi1_algebra.py: M2, about 1 s. g5_checks.py: about 200-250 s (floating point, not a proof).

## What a referee should look at first
(a) region II, Steps 2-4 and 8 of sec_g5_varpi1.tex (strip constants, Cauchy bounds, window remainder, monotonicity in N), since it replaces the notes' 279-cell sweep;
(b) Lemma mu (convexity and the tangent-line tail bounds);
(c) the claim that 71,000 + 120,000 suffices (the true constant is about 1.75);
(d) sec_g5_decay.tex:228-232 (monotonicity of the rest pieces).
