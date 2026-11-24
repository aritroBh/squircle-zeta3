# Referee report: lem:holo, thm:M1 (+ app:M1, lem:dilog), thm:one — paper/ieee/main.tex draft 36

All code here is my own (m1_checks.py, cauchy_Mk.py, growth.py); logs alongside. Nothing outside this folder was touched.

## 1. Lemma 1 (`lem:holo`) — PROVED (one wording repair)
- Definitions checked: P/4 = 2 - eps F, M = -F/2, P = 8(1 + sum M_k p^{-k-1}), M_k = -F^{(k)}(0)/(2 k!). Correct.
- Proof checked: for 0<t<=1/2 and Re eps<1, |z|<=1 with z=1 only at t=1/2, so 1+z is never in (-inf,0];
  |b|<=2^{1-Re eps}<=|a|, |1+sqrt(1+z)|>=1, integrand bounded by 2^{2-Re eps}, uniform on compacts. Holomorphy and
  differentiation under the integral follow. Radius >= 1 since dist(0, {Re eps = 1}) = 1.
- Eq. (M1t) re-derived by hand (differentiate, use t<->1-t); confirmed numerically (route A vs route A', 110 digits).
- Repair (wording only): main.tex says "sqrt(a^2+b^2) := a sqrt(1+z) is holomorphic" without naming the branch.
  Suggested: "... := a_eps sqrt(1+z), with the principal branch (1+z is never in (-inf,0]; Re sqrt(1+z) >= 0), is holomorphic
  in eps and is the positive root for real eps." (squircle_note already has this.)

## 2. Theorem 1 (`thm:M1`), Appendix A, Lemma `lem:dilog` — PROVED AFTER REPAIRS (one missing identity)
Numerics (m1_checks.log, growth.log):
- Closed form vs eq.(M1t) quadrature, vs -F'(0)/2 from the unsymmetrized derivative, vs the v-form integral: agree to 110 digits (|diff| ~ 1e-110).
- Closed form vs M_1 computed from the definition F (Cauchy integral of the Lemma-1 integrand on |eps|=1/4, N=160, dps 100): 89 digits.
- kappa = 8 M_1 = 2.020952769656358356946908... ; printed 2.02095276965635835695... correct; Ahmed's 2.0209527696 matches (truncated).
- 8 M_0 = -4K checked. Signs/factors kappa = 8 M_1 = -4F'(0) consistent in Sec. VIII, IX-B and squircle_note.
Line-by-line:
- Euler substitution: t, 1-t, q, dt, endpoints t(b)=0, t(a)=1 — exact (sympy). The M_1 kernel equals c*/v - 2/(v+a) exactly
  (sympy): the would-be pole at v=b cancels and there are no poles at +-i, so no integration by parts and no boundary term occurs for M_1.
- int dv/v = 2L, int dv/(v+a) = L - lambda/2, the constant -2 lambda^2 + 2 sqrt2 lambda L: correct.
- All six I(x,y) values correct (110 digits); I(0,0) and I(-a,-a) checked by hand.
- Lemma dilog: all six lines correct (110 digits) and re-derived by hand. Domains: chi_2-Landen at x=b in (0,1), fixed point; duplication at w=b;
  reflection at 1-b in (0,1); Landen at -sqrt2 and at (2+sqrt2)/4 (both <1); five-term at (-sqrt2, -b): x<1, y<1, x+y<1, arguments -1, -b^2, b^2 — OK;
  inversion at -a^2 (<0) OK.
- Assembly: c*-bracket = 10 Li2(b) - 19pi^2/24 + 5L^2 and rest = 3pi^2/4 + lambda^2/2 + 2sqrt2 lambda L - 6L^2 - 10 Li2(b): re-derived by hand; sum = Theorem 1.
REPAIR (gap in the printed derivation of the I-table):
  The appendix says "We use the real dilogarithm on (-inf,1]" and lists inversion only for z<0. But the stated substitution
  w=(v-y)/(x-y) gives, for three of the integrals, w-ranges inside [1,inf): I(b,0): [1, a^2]; I(b,-a): [1, a/sqrt2]; I(0,-a): [1+b^2, 2].
  Their values need Re Li2 on the cut (Li2(2) = pi^2/4 and) the inversion Re Li2(z) + Li2(1/z) = pi^2/3 - (1/2)log^2 z for z>1, which is not stated.
  The results printed are right (they are exactly what this inversion gives). Corrected text, either:
  (a) add to the identity list: "inversion for z>1: Re Li2(z) + Li2(1/z) = pi^2/3 - (1/2) log^2 z"; or
  (b) replace the substitution sentence by: "for x != y, the substitution w=(v-y)/(x-y) if x=-a, and w=(x-y)/(v-y) if x is 0 or b,
      turns I(x,y) into logarithms and dilogarithms at arguments in [-1,1]." Then the ranges are I(b,0): [b^2,1]; I(b,-a): [2-sqrt2,1];
      I(0,-a): [1/2,(2+sqrt2)/4]; I(-a,0): [-1,-b^2]. (With w=(x-y)/(v-y): log|v-x| = log|x-y| - log w + log|1-w|, dv/(v-y) = -dw/w.)
      Li2(1/2) then enters I(0,-a) and needs the classical value pi^2/12 - lambda^2/2.
- Minor: "Every logarithm below is that of an S-unit +-2^{j/2}a^m" should read "the logarithm of the absolute value of an S-unit".

## 3. Theorem 11 (`thm:one`) — PROVED, conditional on Theorem `thm:zero` (which is marked dagger); wording repair
- Radius >= 1 from Lemma 1; the only point of |eps|=1 not in Re eps<1 is eps=1: correct. limsup |M_k|^{1/k} = 1 <=> radius exactly 1: correct.
- Gap in wording: "M is not analytic at eps = 1" — M is not defined at eps = 1. The needed step is: if the Taylor series had radius >1, its sum S
  would be holomorphic near eps=1 and equal M on {|eps|<R, Re eps<1} (connected); then for p in [1,1+delta) P(p) = 8(1+S(1/p)/p) is real-analytic
  in p, and since P is C^inf at p=1 (Thm zero), its Taylor series at p=1 equals that of this analytic function and has positive radius,
  contradicting Thm zero. Only right-hand derivatives (p>1) are used, which is legitimate because Thm zero asserts P is C^inf at p=1.
  Suggested proof text: "Lemma 1 gives radius at least 1, and eps=1 is the only point of |eps|=1 outside {Re eps<1}. If the radius exceeded 1,
  the sum of the series would continue M holomorphically across eps=1, so p -> P(p) = 8(1+M(1/p)/p) would agree on p >= 1 with a function
  analytic at p=1. As P is C^inf at p=1, its Taylor series there would then converge, contradicting Theorem zero."
- Dependence: the theorem rests on Theorem zero, which the paper itself marks dagger (sketch here, full proof in notes). The statement of
  thm:one carries no dagger, while Table (line 2068) lists it "Proved dagger". Suggest adding the dagger to the statement, or "(given Theorem zero)".
  Remark: Theorem `thm:decay` (also dagger), if accepted, gives limsup |M_k|^{1/k} = 1 independently, since |cos sqrt(pi N)| >= 1/2 infinitely often.
- Numerics (cauchy_r095.log vs cauchy_r090.log, k <= 220; two circles agree to 1e-12 relative or better): |M_k|^{1/k} rises
  slowly (0.47 at k=10, 0.80 at k=100, 0.864 at k=220), tracking exp(-sqrt(pi/k)); M_k / [-(sqrt(pi)/N) e^{-sqrt(pi N)} cos sqrt(pi N)]
  = 1.05-1.15 for k = 20..220. Consistent with radius one and with eq.(Mk); numerics cannot by themselves prove a limsup.

## Sentences claiming more than proved
- thm:one statement (no dagger) vs Table "Proved dagger": add the dagger / "given Theorem zero".
- Appendix A: "We use the real dilogarithm on (-inf,1]" — false as printed (values at a^2, a/sqrt2, 2sqrt2 b, 2 are used implicitly). Fix as above.
- Nothing in protected Sections I-VII, the abstract or Conclusion bears on these three results; nothing to report there.
