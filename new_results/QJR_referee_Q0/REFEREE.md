# QJR: referee report on QJ (closed form of J; P_1, P_2a, P_2c; vertex constant Q_0 of thm:S)

Object refereed: `new_results/QJ_superellipsoid_Q0/result.md` (and its LaTeX block for
`paper/ieee/proofs/sec_g3_S.tex`). Earlier context: `paper/ieee/proofs/referee_g3/REFEREE.md`, items E5, G2.
All code below is mine, written from the definitions; I did not import or run the QJ scripts.

## Verdicts

| piece | verdict |
|---|---|
| Lemma J (closed form of J) | **CONFIRMED** |
| P_2c (incl. the complex chi_2 steps (a)-(e)) | **CONFIRMED** (the chi_2-Landen continuation is justified; see 2.2) |
| P_2a (Tonelli step, tau-substitution, F_+-, IBP) | **CONFIRMED** |
| P_1 (Euler substitution, letters, T1-T6, endpoint pole, Z_1-Z_3) | **CONFIRMED AFTER REPAIR** (R3: the 18 contributions must be printed, not only logged) |
| Table of 15 dilogarithm instances | **CONFIRMED** (all 15: identity, domain, branch, value) |
| Assembly 6(P_1-P_2a-P_2c) = Q_0 | **CONFIRMED** (exact 0) |
| LaTeX for the supplement | **CONFIRMED AFTER REPAIR** (R1-R5 below) |

Status word for the closed form of Q_0: **proved** (after R1-R5). Every step is a written hand
argument using only classical dilogarithm identities with stated domains; the computer-algebra checks
confirm it but are not needed for it.

## 1. Lemma J (`r1_J.py`, `r1_J.log`)

- Substitution u=r+R: I re-derived r=(u^2-q^2)/(2u), R=(u^2+q^2)/(2u), dr=(u^2+q^2)/(2u^2) du, and the
  transformed integrand equals g(u) exactly (sympy, `simplify` = 0). The factor u-q cancels, so g has no
  pole at u=q; g=O(u^{-3}).
- Partial fractions checked exactly modulo B^2=1+q^2; residues sum to 0 (so F(U)->0); F'=g with the true
  beta; the two residue simplifications ((1/beta -+ 1/(1+q))), (q+1)^2-beta^2=2q,
  (q+1+beta)^2=2(q+beta)(1+beta) and the log(2q) coefficient all exact.
- Positivity of u-u_k on [q,oo): q+1-beta = 2q/(q+1+beta) > 0. Correct.
- Original definition of J by quadrature (60 digits, 8 values q = 1/sqrt2, 0.72, 0.75, 0.85, 0.95,
  0.999, and also 0.3, 1.7 outside the range) vs closed form: all <= 2.3e-62 (working precision).

## 2. P_2c (`r2_P2.py`, `r2_P2.log`, `r6_claims.py`)

2.1 Substitutions. beta=(x^2+1)/(2cx): beta(1)=sqrt(3/2), beta(sqrt3)=sqrt2, 2beta^2-3=3(x^2-1)^2/(4x^2),
dbeta/sqrt(2beta^2-3)=dx/(sqrt2 x), (beta+1)/(beta-1)=(x^2+2cx+1)/(x^2-2cx+1): all exact. The x-form
integral equals G(sqrt3)-G(1) to 1.7e-111; G'(x) equals the integrand (numerical derivative, 1e-111).
Branch: Im(x omega)=x/sqrt3 != 0 for x>0, so x omega, -x omega never meet [1,oo); principal
log(1-z) is the right branch of d Li_2/dz there. Correct.

2.2 The complex steps.
- (a) Fourier value at theta=alpha_c and alpha_c+pi in [0,2pi]: correct (1.7e-111).
- (b) inversion at z and -z for non-real z; principal log z = (1/2)log3 + i alpha_c,
  log(-z) = (1/2)log3 + i(alpha_c-pi): correct (3.3e-111). 1/(sqrt2+i)=(sqrt2-i)/3 exact.
- (c) **chi_2-Landen continuation — justified.** I checked the domain myself: chi_2 is analytic off
  (-oo,-1]U[1,oo); m(t)=(1-t)/(1+t) is real exactly for real t, with m(t) in [1,oo) iff t in (-1,0] and
  m(t) in (-oo,-1] iff t < -1 (pole at t=-1), so chi_2(m(t)) is analytic on C\(-oo,0]; log t needs
  C\(-oo,0]; log((1+t)/(1-t)) needs t not in (-oo,-1]U[1,oo). The intersection
  C\((-oo,0]U[1,oo)) is open and connected and contains (0,1), where the real identity holds; the
  identity theorem gives it on the whole set, with principal branches. t=(sqrt2-i)/3 lies in it.
  Numerically the complex identity holds on 127 test points of that set, including points 1e-8 above
  and below both cuts (max residual 1.1e-39 at 40 digits). (1-t)/(1+t)=b e^{i pi/4},
  (1+t)/(1-t)=a e^{-i pi/4} exact. Step value correct (1.7e-111).
- (d) duplication at |y|=b<1, y^2=ib^2, Re Li_2(ix)=Li_2(-x^2)/4 (power series): correct (4e-112).
- (e) Lemma A at theta=pi/4: correct (4e-112).
- Solving (b)-(e) symbolically and inserting (a) gives exactly the printed P_2c (sympy, exact).

## 3. P_2a (`r2_P2.py`)

- q=sinh(phi): c^2cosh^2phi-1=(2sinh^2phi-1)/3 exact; G(phi_0)=arcosh 1=0, G(L)=arcosh(2/sqrt3)=(1/2)log3.
  G' ~ (phi-phi_0)^{-1/2} near phi_0: integrable; G absolutely continuous, so the integration by parts
  is valid (boundary term phi_0 G(phi_0)=0).
- Tonelli: the integrand cosh(phi)/sqrt(k^2cosh^2phi-1) is positive on the region
  {phi_0<=phi<=L, 1/cosh(phi)<=k<=c}; its sections are as stated (1/cosh(phi_0)=c, 1/cosh L=1/sqrt2,
  arcosh(1/k)>=phi_0 for k<=c), and the total is W < oo. Inner primitive (1/k)arcosh(sinh(phi)/mu)
  checked; mu<=1 iff k>=1/sqrt2. Correct.
- tau-substitution, dk/k=tau dtau/(2-tau^2), tau/(2-tau^2)=K(tau)/2, the reflection tau->-tau: exact.
- F_-, F_+ derivatives exact; arguments w/a in [b^2/sqrt2, 1/sqrt2] subset (0,1), -w/b <= -1/sqrt2 < 0:
  both real and continuous on the interval.
- W in its five forms agrees to 2e-112. My own exact evaluation of
  W=(1/4)([F_-]-[F_+]) (exact S-unit logarithms) equals the printed W exactly; with table lines 2, 3, 15
  the bracket is -12 Li_2(b)+Li_2(-b^4)-2L^2+5pi^2/6 and P_2a equals Lemma P exactly.

## 4. P_1 (`r4_P1.py`, `r4_P1.log`)

I rebuilt P_1 independently: own letter decomposition, own residues by limits, own primitive choice
(T1 for the pole a; T2/T3/T4 by the sign of d=c+p; T5; T6; plain log for lambda), own exact endpoint
evaluation with exact S-unit logarithms (each log decomposition is verified exactly in Q(sqrt2)).
- Substitution: q(1)=1/sqrt2, q(a)=1, q'(v)=sqrt2(v^2-1)/(4v^2) > 0, sqrt(2q^2-1)=(v^2-1)/(2v),
  the measure (v^2+1)/(2v^2) dv, and the three factorisations: all exact.
- The five coefficient functions are sums of simple poles only, and equal the printed C_l exactly.
- Endpoint pole: sum_l r_(l,a) log l(a) + r_(lambda,a) lambda = 0 exactly.
- Every primitive: F' = integrand (exact differentiation after Li_2'(z) = -log(1-z)/z, then evaluation
  at three interior points to 1e-50; T6 by numerical differentiation of -2 Re Li_2(Z(v)), 1e-40).
  Hypotheses hold on [1,a]: every pole c != a has c < 1; for T3, v-c > |d|.
- T6 branch: Z(v)=(c-v)/(c+i) has Im Z = -(c-v)/(c^2+1), real only at v=c where Z=0; Z never in [1,oo).
- Z_1, Z_2, Z_3 data (Re Z_1=1/2, 1-Z_1=conj Z_1, |Z_1|^2=b/sqrt2, 1-Z_2=a e^{-i pi/4},
  |Z_2|^2=sqrt2 a, arg Z_2=5pi/8, Z_3/(Z_3-1)=b e^{-i pi/4}, |1-Z_3|^2=a/sqrt2, arg(1-Z_3)=-pi/8,
  tan(pi/8)=b, tan(3pi/8)=a): all exact (`r6_claims.log`). I redid the three reductions by hand
  (reflection; reflection + inversion + Lemma A; Landen + Lemma A) and got the printed values; numerically
  <= 2.2e-112.
- My unreduced P_1 minus the printed 18-term expression = **0 exactly**; reduced with my own table values,
  P_1 minus Lemma P = **0 exactly**; unreduced value vs v-form quadrature 2.2e-112.
- **Gap (repair R3):** the LaTeX says "Adding the eighteen terms (listed one by one in 04_P1.log)". For a
  hand proof the eighteen contributions must be in the text; they are routine but not obvious.

## 5. The 15 identity instances (`r3_table.py`, `r3_table.log`)

Re-derived from scratch with my own implementations of reflection (0<x<1), inversion (x<0), Landen
(x<1), duplication (|x|<=1) and Abel's five-term relation (x<1, y<1, x+y<1, all five arguments < 1,
asserted in code). For each line: the identity named is the one that produces the stated value; the
argument is in the stated domain; every other dilogarithm on the line is already known (enforced in
code: Li_2(b), Li_2(-b^4), Li_2(1/2), Li_2(-1), or an earlier line); the value solved exactly equals
the printed value (15/15). Residual of each instance at 110 digits: <= 7.8e-121; each obtained value vs
direct mpmath Li_2 at 120 digits: <= 1.6e-120. The Abel side data in the text (arguments -1, -1/sqrt2,
1/sqrt2; -2a, b/(2sqrt2), -1/sqrt2; -a^4, b^2/sqrt2, -a^2/sqrt2; right sides; x+y=-5 in (15)) are
correct. All arguments are real, so there is no branch question beyond the stated real domains.

## 6. Direct quadrature of the original definitions (`r9_quad.py`, `r9b_k3.py`)

| quantity (original definition) | digits | |quad - closed form| |
|---|---|---|
| J(q), 8 values of q (sect. 1) | 60 | <= 2.3e-62 |
| P_1 = int_0^1 f(q(s)) ds | 60 | 1.5e-75 |
| P_2a = int_0^1 arsinh(q)/sqrt(1+q^2) ds | 60 | 1.2e-75 |
| P_2c = int_0^1 arsinh(1/q)/sqrt(1+q^2) ds | 60 | 6.9e-76 |
| 6(P_1-P_2a-P_2c) vs Q_0 closed form | 60 | 2.2e-75 |
| Q_0 = 6 int_0^1 J(q(s)) ds, J by inner quadrature (nested) | 45 | 6.2e-61 |
| Q_0 = -pi^{-1/2} int_0^oo k_0(y)^3 dy (nested; x-form, coarse) | 42 | 7.6e-33 (inner-quadrature limited) |
| Q_0 = -pi^{-1/2} int_0^oo k_0(y)^3 dy (nested; u-form, finer split) | 45 | 1.5e-62 |

Q_0 = -2.83293568349484776046205725151576921593768332376794122833171353128521157348857646879972043650351108...

## 7. Assembly (`r5_assembly.py`)

6(P_1 - P_2a - P_2c) - Q_0(theorem) = **0** exactly, with P_1 from my own reduction (r4), P_2a and P_2c
from my own exact reductions (r2).

## 8. The LaTeX: repairs

Order and content are right (Lemma J; then Step 3 after Lemma P, using Lemma A, which precedes it), but:

**R1 (undefined symbols / conventions).** chi_2 and the branch conventions are defined only in the prose
of result.md, not in the LaTeX block, yet chi_2 is used in the P_2c paragraph. Insert at the start of
Step 3, before "We use the following classical identities":

```latex
All logarithms and dilogarithms are principal: $\log$ is analytic on $\CC\setminus(-\infty,0]$,
$\Li_2$ is analytic on $\CC\setminus[1,\infty)$, $\frac{\dd}{\dd z}\Li_2(z)=-\frac{\log(1-z)}z$ there, and
$\Li_2(\bar z)=\overline{\Li_2(z)}$. We write $\chi_2(z)=\frac12\bigl(\Li_2(z)-\Li_2(-z)\bigr)$; it is analytic on
$\CC\setminus\bigl((-\infty,-1]\cup[1,\infty)\bigr)$.
```

and in (c) replace "both sides are analytic on the connected set
$\CC\setminus((-\infty,0]\cup[1,\infty))$, so it holds there." by

```latex
both sides are analytic on the connected set $\CC\setminus((-\infty,0]\cup[1,\infty))$
(for real $t$, $\frac{1-t}{1+t}\in[1,\infty)$ only if $-1<t\le0$, and $\frac{1-t}{1+t}\le-1$ only if $t<-1$;
$\log t$ and $\log\frac{1+t}{1-t}$ are analytic there), so it holds there.
```

**R2 (headings and the old J text).** The new "Step 2" lemma is inserted under the existing heading
"Steps 2 and 3: evaluation (computer-assisted)" and after the old sentence describing a different
method (rationalise, poles $0,-1,-s(1-s)$) with the $J$ display. Replace the heading by
`\subsubsection*{Steps 2 and 3: evaluation}`, replace the sentence "The inner integral $J$ is elementary
(rationalise ... $\frac{1+q^2}2$):" and its display by

```latex
The inner integral $J$ is elementary (Lemma~\ref{g3:lem:J}), and so $Q_0=6(P_1-P_{2a}-P_{2c})$.
```

drop the new `\subsubsection*{Step 2: ...}` heading (Lemma J follows directly), and wrap the Step 3
text from "The integrands of $P_1,P_{2a},P_{2c}$ are continuous" to "this proves Lemma~\ref{g3:lem:P}." in
`\begin{proof}[Proof of Lemma~\ref{g3:lem:P}] ... \end{proof}` (keeping the table outside the proof as a float
is fine).

**R3 (the eighteen terms of P_1).** Replace "Adding the eighteen terms (listed one by one in
\path{04_P1.log}) gives" by "The eighteen contributions are listed in Table~\ref{g3:tab:P1}; adding them gives",
and add:

```latex
\begin{table}[ht]
\centering\footnotesize\setlength{\tabcolsep}{3pt}
\caption{The eighteen contributions to $P_1$ (letter, pole $c$, residue, primitive, value of the primitive between $v=1$ and $v=a$).}\label{g3:tab:P1}
\begin{tabular}{@{}llrll@{}}
letter & $c$ & res. & prim. & contribution\\\hline
$\log(v+a)$ & $a$ & $-2$ & T1 & $-2\Li_2(\frac b{\sqrt2})$\\
$\log(v+a)$ & $-a$ & $2$ & T4 & $\lambda L+\frac34\lambda^2$\\
$\log(v+a)$ & $b$ & $2$ & T2 & $3\lambda L+\frac32\lambda^2-2\Li_2(-\frac1{\sqrt2})+2\Li_2(-\frac b2)$\\
$\log(v+a)$ & $-b$ & $-2$ & T2 & $-2\lambda^2+2\Li_2(-\sqrt2)-2\Li_2(-\frac1{\sqrt2})$\\
$\log(v+b)$ & $a$ & $-2$ & T1 & $-2\Li_2(\frac12)$\\
$\log(v+b)$ & $-a$ & $2$ & T3 & $\lambda L+\frac34\lambda^2+2\Li_2(b)-2\Li_2(\sqrt2\,b)$\\
$\log(v+b)$ & $b$ & $2$ & T2 & $-2L^2+\lambda L+\lambda^2+2\Li_2(-\frac1{\sqrt2})-2\Li_2(-a)$\\
$\log(v+b)$ & $-b$ & $-2$ & T4 & $-2\lambda^2$\\
$\log v$ & $0$ & $-\sqrt2$ & T4 & $-\frac{\sqrt2}2L^2$\\
$\log v$ & $-a$ & $-2$ & T3 & $-\lambda L-\frac34\lambda^2-2\Li_2(\frac12)+2\Li_2(\frac1{\sqrt2})$\\
$\log v$ & $-b$ & $2$ & T3 & $2\lambda^2+2\Li_2(\frac b{2\sqrt2})-2\Li_2(\frac b{\sqrt2})$\\
$\log(v^2+1)$ & $0$ & $\sqrt2$ & T5 & $\frac{\sqrt2}2\Li_2(-1)-\frac{\sqrt2}2\Li_2(-a^2)$\\
$\log(v^2+1)$ & $a$ & $2$ & T6 & $4\operatorname{Re}\Li_2(Z_1)$\\
$\log(v^2+1)$ & $b$ & $-2$ & T6 & $2L^2-2\lambda L-\frac32\lambda^2+4\operatorname{Re}\Li_2(Z_2)-4\operatorname{Re}\Li_2(Z_3)$\\
$\lambda$ & $0$ & $-\frac{\sqrt2}2$ & $\log$ & $-\frac{\sqrt2}2\lambda L$\\
$\lambda$ & $-a$ & $-3$ & $\log$ & $-\frac32\lambda^2$\\
$\lambda$ & $b$ & $-2$ & $\log$ & $-2\lambda L-\lambda^2$\\
$\lambda$ & $-b$ & $3$ & $\log$ & $3\lambda^2$
\end{tabular}
\end{table}
```

(For the pole $c=b$ of $\log(v^2+1)$ the entry includes $-2\log(b^2+1)\,[\log(v-b)]_1^a$; for the pole
$a$ the constants $\log\ell(a)$ have been removed, as explained in the text. I computed every entry
independently, `r4_P1.log`.)

**R4 (symbol clashes, optional but advisable).** In the same section $G$ already denotes
$\Gamma(1+\varepsilon)^2/\Gamma(1+2\varepsilon)$, and $c$ is both $\sqrt{2/3}$ (P_2c, P_2a) and the pole
variable (P_1, T2-T6). Rename the antiderivatives in P_2c and P_2a to $\mathcal G$, and $c=\sqrt{2/3}$ to
$c_0$.

**R5 (other places that still say "computer-assisted").** Besides the Plan paragraph and "End of the
proof" (both noted in result.md): in "Which report confirmed the proof", replace "The step not written
out by hand (Lemma~\ref{g3:lem:P}) is in the box above." by "Steps~2 and~3 are now written out by hand
(Lemma~\ref{g3:lem:J} and the proof of Lemma~\ref{g3:lem:P}); an independent referee re-derived every step
(\path{new_results/QJR_referee_Q0/REFEREE.md})." In `main.tex`, the summary-table row of
Thm.~\ref{thm:S} reads "Proved$^\dagger$ (computer-assisted; refereed)"; once nothing else in that row
rests on a computer step, it should read "Proved$^\dagger$ (refereed)". Also use full paths
(`new_results/QJ_superellipsoid_Q0/01_J.py`, ...) in the `\path{}` references.

## 9. Status and the sentence for main.tex

Status word: **proved**.

Current sentence (`paper/ieee/main.tex`, after `thm:S`): "The closed form of the vertex constant $Q_0$
is proved, with computer assistance, by an exact reduction of 54 dilogarithm values with the classical
functional equations; it brings in three new constants, $\log3$, the cube angle $\alpha_c$ and
$\Li_2(-(\sqrt2-1)^4)$."

Replacement:

```latex
We prove the closed form of the vertex constant $Q_0$ by hand. The inner integral is elementary, and the three integrals that remain reduce, by fifteen instances of the classical dilogarithm identities, to $\Li_2(\sqrt2-1)$ and $\Li_2(-(\sqrt2-1)^4)$. We also checked every step exactly by computer algebra. The closed form brings in three new constants, $\log3$, the cube angle $\alpha_c$ and $\Li_2(-(\sqrt2-1)^4)$.
```

## Files

`rcommon.py` (exact Q(sqrt2) arithmetic, exact S-unit logs), `r1_J.py`/`.log`, `r2_P2.py`/`.log`,
`r3_table.py`/`.log`/`r3_table_values.json`, `r4_P1.py`/`.log`/`r4_P1_unreduced.json`/`r4_P1_reduced.json`,
`r5_assembly.py`/`.log`, `r6_claims.py`/`.log`, `r9_quad.py` with `r9_quad_P.log`, `r9_quad_J.log`,
`r9_quad_K.log`, `r9b_k3.py`/`.log`. Python: `<project>/.venv/bin/python`.
