# REPORT_G8C: Theorem 20, the sentence "W_3 is an iterated integral of weight three whose letters are all real"

Files (all new, in `paper/ieee/proofs/`): `sec_g8_W3int.tex` (section, 7 pages, pp. 32-38 of the build), `macros_g8c.tex`, `wrap_g8c.tex`,
`scripts/g8c_structure.py` + `.out` (exact, 4 s), `scripts/g8c_numeric.py` + `.out` (13 s), this report. No existing file was touched.

## The exact reading proved (Theorem "Reading proved here", `g8c:thm:W3`)
With c = cos(alpha/2), s = sin(alpha/2), Y* = 1/s, letters a_0=0, a_1=-1, a_2=-1/(1-c^2), a_3=-1/(1+c), a_4=-1/(1-c) (all real, in Q(c), <= 0),
lambda_j = log(1 - y/a_j) (lambda_0 = log y), l_g = lambda_0+lambda_2-lambda_1, l_m = lambda_3+lambda_4-lambda_1, q = (1+c^2)/(4c):
* (a) W_3 := M_2 - W_2 equals sum_j int_0^{Y*} [A_j l_g^2 + B_j l_g l_m + C_j l_m^2] d lambda_j, with (A,B,C)_j =
  (0,-1,1/2), (1/2,0,1/2), (-1,1,-1/2), (-q,c,-c), (q,-c,c) for j=0..4. Absolutely convergent.
* (b) Hence W_3 is a Q(c)-combination (coefficients in Q + Qc + Q/c; 64 nonzero terms) of weight-three iterated integrals I(j,a,b; Y*) with
  letters in {0,-1,-1/(1-c^2),-1/(1+c),-1/(1-c)} and real end point Y* = (1-c^2)^(-1/2). Equivalently one integral of a one-form with
  quadratic coefficients in weight-one functions.
* (c) W_2 + W_3 = M_2; W_3 = Y(c)+Y(-c) - zeta(2)(1-c^2)/(2c) artanh c (exact, through Theorem 20 of `sec_g8_M2alpha.tex`; the integral
  representation (a) itself does not use Theorem 20). The complex letters e^{+-i alpha} occur (double poles) only in an exact derivative V'
  that integrates by parts to W_2.

## Hand / machine / numerical
* Hand: differentiation under the integral and the u^2 coefficient (Lemma 6.1, with the domination g^(-1/2)(1+|log g|)^(2j)); Euler substitution
  (Lemma 6.2); closed forms of F_B, F_C; V(1)=0, V dlog g/dt = 1/(2(t+1)), g(t1)=1; the integration by parts (Lemma 6.3); the identification
  W_2 = zeta(2)/2 - Li2(c^2)/2 by differentiating in c (Lemma 6.4); the substitution t = 1+2 s^2 y and the letters (Lemma 6.5); the shuffle/product
  step giving I-integrals; the convergence checks.
* Machine (one named step, **N1**, exact sympy in Q(c)(t) and, for the y-variable, in Q(tau)(y) with c=(1-tau^2)/(1+tau^2), s=2tau/(1+tau^2)):
  (N1.0) the hand formula for f_A,f_B,f_C equals the u^2 coefficient of the u-expansion of the integrand, F_B, F_C closed forms, Euler formulas;
  (N1.1) F_A = V' + F_A^s, V = (t-1)(t-k)/(2P), F_A^s with the five simple poles and residues (0,-1,1/2,-q,q), no polynomial part;
  (N1.2) F_B, F_C partial fractions and residues. Script `g8c_structure.py`; also prints the 64 coefficients of part (b).
* Numerical (40+ digits, 48 working digits; angles pi/3, pi/2, 2pi/3, 1.8, 5pi/6; the last three obtuse): M_2 from the definition (Cauchy in u),
  g-integral, Euler t-integral (agree to 1e-50), W_2 integral vs zeta(2)/2-Li2(c^2)/2 (<= 4.3e-50), y-integral W_3^sharp vs M_2 - W_2, and vs the
  closed form Y(c)+Y(-c)-...: all agree to 8e-46..1.3e-45 (the aliasing limit of the Cauchy integral).

## New argument, needs referee (file:line in `sec_g8_W3int.tex`)
* 59-99: Lemma 6.1, differentiation of the corner integral under the sign twice and the closed form of the u^2 coefficient. The notes only
  re-derived it by sympy (`m2_structure.py` [1]); the domination argument (Faa di Bruno, h <= g^(3/4)) is ours.
* 105-113 and 125-132: Euler substitution with domain and monotonicity, hand derivation of F_B, F_C (the notes only print them).
* 155-181: Lemma 6.3 (IBP with boundary terms; the notes and verify.md 3.3 state them) and Lemma 6.4 (W_2). The proof of Lemma 6.4 is the
  derivative argument of `verify.md` 3.3, re-written with Leibniz rule (singular part log(t-1) independent of c) and the limit c -> 0+ by
  dominated convergence; the notes' own proof used three dilogarithm antiderivatives checked numerically.
* 183-195: Lemma 6.5. The change of variable y = x/s is ours: it makes all letters rational in c and removes log s from the logarithms; the
  notes use x = (t-1)/(2s) in [0,1] with letters 0, s, 1/s, tau, 1/tau (our letters times s) and log(x/s).
* 197-212: part (b), the shuffle/regularisation convention and the statement that exactly 64 terms with coefficients in Q+Qc+Q/c occur.

## Disagreements and readings
* Paper: "W_3 is an iterated integral of weight three whose letters are all real" (main.tex:1511, verbatim). Notes: a weight-3 integral on the
  real alphabet in x in [0,1], with coefficients (A,B,C) in Q(c) (the table of result.md section 1 item 4). Proved here: a Q(c)-linear
  combination of such integrals, so the paper's singular "an iterated integral" is true only in that reading; the section says so after the
  theorem. The end point is the algebraic number 1/sqrt(1-c^2) (or, after rescaling to 1, the letters carry sqrt(1-c^2) and the logs carry log s).
* The notes call W_3 "pure weight 3". Here that means length three; the coefficients depend on c. No zeta(2) factor appears in (a).
* No numerical disagreement between the notes (x-table) and our y-table: they are the same residues; the x-letters are -s*a_j.
* verify.md corrections (section 5) concerning Theorem C: only item 4 (the short W_2 proof): applied. The others concern Theorem A / B / the
  phrasing of the zeta(3) coefficient and are outside this section.
* The supplement's earlier sentence (gap box in `sec_g8_M2alpha.tex:35`) says the integral representation is "not reproved"; after this section
  is spliced that box can be rewritten to point to Section `g8c:sec:W3` (we did not edit it).

## Which report confirmed what
`G1/verify.md` 3.3 confirmed the NOTES' Theorem C (sympy re-run, IBP boundary terms, letter order, quadrature 5.7e-47 at 13 angles) and
supplied the W_2 derivative proof. It did not see this text. This text is checked only by `g8c_structure.py` and `g8c_numeric.py`; unrefereed.

## Gaps (`\gapnote` at the end of the section)
(i) N1 is a machine step, not formally verified. (ii) A single constant-coefficient iterated integral / minimality is not claimed. (iii) The
64 coefficients are printed, only their form is stated. (iv) Part (c)'s identity with Y(c)+Y(-c)-... is exact only via Theorem 20, itself
computer-assisted (M1, M2).

## Build
`tectonic -C --keep-logs wrap_g8c.tex`: 0 undefined references, 0 Overfull boxes, 0 "??" in pdftotext, no font errors; 54 pages in all, the
new section is Section 6, pages 32-38. `wrap_g8c.tex` prints sec_g6_*, `sec_g8_M2alpha.tex` (as at build time), the new section, and then
`sec_g8_M3alpha.tex`, `sec_g8_M4alpha.tex` (needed so that references from the M2 section resolve). `macros_g8c.tex` has one guarded
`\newtheorem*{gcmain}` besides `\providecommand` lines (needed for an unnumbered statement environment that does not disturb the theorem counter).
Temp files: the scratch folder was deleted; `wrap_g8c.pdf/.log` remain as build products.
