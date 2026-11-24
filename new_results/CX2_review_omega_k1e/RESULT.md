Claim A passes, with the range stated as **odd n >= 3** and one sentence making the infinite composition precise. Claim B **fails as displayed**: four logarithmic terms have the wrong signs. Its implementation uses a different logarithm convention and passes the new independent test after that convention is restored. Its assertion of canonicality is not proved.

All computations below were run with `<project>/.venv/bin/python`. Only this directory was modified. `review.py` contains the runnable checks; `EVIDENCE.txt` preserves the printed outputs. Numerical agreement here is high precision evidence, not a certified interval bound or a proof for all hexagons.

**Claim A: the proposed proof is complete after specifying the formal composition.**

Use the square-root branches with value 1 at zero. The factorization is correct:

\[
g(t)=\sqrt{(1+e^{4t})/2}=e^t\sqrt{\cosh(2t)}.
\]

Let \(\mathcal H=\{\sum_{m\ge0}h_m t^m/m!:h_m\in\mathbb Z\}\). Its product rule is

\[
(fg)^{(m)}(0)=\sum_{k=0}^m\binom mk f^{(k)}(0)g^{(m-k)}(0),
\]

so it is a ring. For \(z=(\cosh(2t)-1)/4\), the constant and odd derivatives vanish, and the derivative of even order \(m\ge2\) is \(2^{m-2}\). In particular \(z\in\mathcal H\) and its order of vanishing is two.

The ordinary Taylor series

\[
\sqrt{1+4y}=1+\sum_{j\ge1}2(-1)^{j-1}C_{j-1}y^j
\]

is correct. Its first printed coefficients are
`2, -2, 4, -10, 28, -84, 264, -858`.
Every nonconstant coefficient is an even integer.

The essential missing sentence is: **for a fixed derivative order m, only j <= floor(m/2) contributes after substituting y=z(t)**. Thus this infinite sum is coefficientwise finite in the formal power-series ring. Every contributing power \(z^j\) is Hurwitz, and multiplying by its even coefficient makes every derivative even. No interchange of an infinite numerical sum with differentiation is required. Equivalently, the analytic identity also holds in a sufficiently small neighborhood of zero.

Consequently, if \(c_m=(\sqrt{\cosh(2t)})^{(m)}(0)\), then \(c_0=1\) and \(c_m\in2\mathbb Z\) for every \(m>0\). The product with \(e^t\) gives

\[
V_n=\sum_{k=0}^n\binom nk c_k\equiv1\pmod2.
\]

For \(n\ge2\), \(2^{n-1}\) is even, so \(V_n\ne2^{n-1}\). The paper's displayed formula therefore implies \(\omega_n\ne0\) for every odd \(n\ge3\).

The qualification on n matters: the printed exact calculation gives \(V_1=1\) and the formal expression \(\omega_1=0\); \(\zeta(1)\) is not a finite zeta value. Do not write “every odd n” without specifying the intended range.

**Exact arithmetic and consistency with the paper.**

Two separate square identities were used:

\[
c_n=\frac{h_n-\sum_{k=1}^{n-1}\binom nkc_kc_{n-k}}2,
\qquad
V_n=\frac{2^{2n-1}-\sum_{k=1}^{n-1}\binom nkV_kV_{n-k}}2,
\]

where \(h_n=2^n\) for positive even n and zero for odd n, with \(c_0=V_0=1\). These are obtained from \((\sqrt{\cosh(2t)})^2=\cosh(2t)\) and \(g(t)^2=(1+e^{4t})/2\), respectively. Every division, parity assertion, and equality with the binomial-transform formula was checked exactly through n=61. All V_n, including even indices and n=0, were printed; their residues modulo two are all one. The full integers are in `EVIDENCE.txt` under `exact.log`.

| n | printed V_n | resulting omega_n |
|---|---:|---:|
| 3 | 7 | 1 |
| 5 | 1 | -3/4 |
| 7 | 967 | 129/16 |

These match the coefficients of the standalone \(\zeta(3)\), \(\zeta(5)\), and \(\zeta(7)\) monomials in `main.tex:1084`, `main.tex:1087`, and `main.tex:945`. The other \(\zeta(3)\) terms in a_5 do not affect the \(\zeta(5)\) check. Negative V_n also occur; the argument correctly uses parity, not positivity.

This verifies the proposed nonvanishing deduction and the requested low-order consistency checks. It is not an independent derivation of the paper's formula for omega_n at every order. Also, “coefficient” means the coefficient in the stated formal zeta-polynomial representation; the proof does not establish independence of numerical zeta values or their products.

**Claim B: an explicit counterexample to the formula as printed.**

I used the strictly convex hexagon with successive vertices

\[
(2,0),(3,2),(0,1),(-2,0),(-3,-2),(0,-1),
\qquad w=(2,5,7),\quad k=(3/32,1/24,1/84).
\]

Thus \(w_jk_j=(3/16,5/24,1/12)\). Exact checks printed \(N(V_q)=(1,1,1)\) and consecutive edge cross products `(5, 1, 3, 5, 1, 3)`, so this is a valid nondegenerate example. It differs from both reference hexagons and K1e's fresh hexagon. The weights remain outside the p-th powers, exactly as requested.

For the three successive half-boundary edges, the derived data are:

| edge | A0 | A1 | B | r | L2 | q0 | V02 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 3/4 | 5/6 | 1/6 | 3/2 | 5 | 2 | 4 |
| 2 | 5/8 | 1/4 | 3/4 | 1/2 | 10 | -11 | 13 |
| 3 | 1/6 | 3/8 | 5/8 | 4/3 | 5 | -1 | 1 |

The constant logarithms supplied to the coefficient polynomials were:

| edge | lA0 | lA1 | lB |
|---|---|---|---|
| 1 | log(3/8) | -log(6) | -log(42) |
| 2 | -log(8) | -log(28) | log(3/8) |
| 3 | -log(42) | log(3/16) | -log(8) |

All normalization and parameter constraints were checked exactly. I also parsed `formula_table.txt` independently and checked every numerator against `Earc_basis_clean.pkl` after A1=1-B and kap=r-1: **all eleven differences were exactly zero**. The discrepancy is therefore present in both supplied representations.

The cause is visible in K1e's own source. `step_23.py:14` defines `l1=-ln(1-x)`. Its `namemap` at lines 36-39 then labels l1 and its odd products as **positive** `ln(1-x)`. The validation code `step_48.py:26` continues to evaluate `l1=-mp.log(1-x)`. Since x=1-r, that is **-log(r)**, whereas the displayed formula in `result.md:9` uses **+log(r)**.

To use the existing numerical coefficients correctly, replace each occurrence of the logarithmic atom ln(1-x) in that formula by ell=-ln(1-x)=-log(r). Alternatively, retain ln(1-x)=log(r) and negate the four coefficient blocks labeled:

- `Li2(x)*ln(1-x)`;
- `ln^3(1-x)`;
- `ln(1-x)`;
- `zeta2*ln(1-x)`.

The squared-log block and the constant parameter logs lA0, lA1, lB retain their existing values.

**Independent polar computation, without K1e's derivative engine.**

For \(u=(\cos\theta,\sin\theta)\), set

\[
s_j=k_j\operatorname{cross}(V_j,u),\quad
H_p=\sum_jw_j|s_j|^p,\quad
T_p=\sum_jw_j|s_j|^{p-1}\operatorname{sgn}(s_j)
             k_j\operatorname{cross}(V_j,u_\theta).
\]

Then \(\partial_\theta H_p=pT_p\), \(\rho=H_p^{-1/p}\), and the polar perimeter integral is

\[
P(p)=\int_0^{2\pi} H_p^{-1/p}\sqrt{1+(T_p/H_p)^2}\,d\theta.
\]

I integrated this expression on each angular sector, using only the change of variable
\(y(t)=(1-t)V_q+tV_{q+1}\), \(\theta=\arg y(t)\),
\(d\theta/dt=\operatorname{cross}(V_q,V_{q+1})/|y(t)|^2\).
This removes endpoint angle cancellation and retains the polar definition. The code evaluates actual perimeters at neighboring p values, not Taylor jets or the supplied symbolic integrand. The values were combined only after integration:

\[
D(h)=\frac{P(1+2h)-2P(1+h)+2P(1-h)-P(1-2h)}{2h^3},
\]

\[
R_{k,0}=D(h_0/2^k),\qquad
R_{k,j}=\frac{4^jR_{k,j-1}-R_{k-1,j-1}}{4^j-1}.
\]

The integral was split at every kink direction, with tanh-sinh quadrature. The p=1 result was checked against the exact polygon perimeter \(2(2\sqrt5+\sqrt{10})\). Central samples below one are valid polar curves even when convexity is lost; the perimeter integral still applies. Finite-order Richardson extrapolation does not require a convergent infinite Taylor series.

| calculation | printed P'''(1) |
|---|---|
| Direct polar, 100 dps, h0=0.001, six levels | 11.17600084483001768696522886994373866088609545038151929536816322691467 |
| Direct polar, 120 dps, h0=0.0005, seven levels | 11.1760008448300176869652288699428381889157919133613739858079886723487 |
| K1e coefficients with ell=-log(r) | 11.17600084483001768696522886994283818891576440763410415183162821869985 |
| Literal displayed formula, log(r) | -18.88376858042687109242138455043862587994939790558884650147675379527205 |

The two direct runs differ by `9.00471970303537020145309560175e-31`, comfortably exceeding the requested twenty-digit check. The higher-resolution Richardson diagonal changed by `2.20292429005e-34` on its last refinement. Relative to its final value:

- corrected formula minus direct = `-2.75057272698339763604536488497e-41`;
- literal formula minus direct = `-30.05976942525688877938661342038146406886518981895`.

The per-edge comparison also passes after the sign correction:

| edge | corrected formula E_arc | formula minus direct, 120 dps |
|---|---|---|
| 1 | 23.9395250569344214178704457082666678109914363541563617501610847208138390638437157354347553 | -1.3752862705708576118e-41 |
| 2 | -5.93306015204506855108311981062137641003147251965939200654068759321083529856301036519968888 | -1.2025980217823685162e-48 |
| 3 | -12.418464482474344023304711462673872306502081630679917667704583018253078878292678566288527 | 2.7338960971967066349e-49 |

Every perimeter sample and Richardson row was printed and is retained in `EVIDENCE.txt`. These residuals are observed differences, not rigorous error estimates.

**The local zeta coefficient survives; its claimed canonicality does not.**

Using the supplied coefficient polynomial, the general symbolic difference between the **full coefficient** of zeta(3) and

\[
Z_e=\frac{6A_1(A_1-B\kappa)(L_2V_{02}-q_0^2)
       (3A_1L_2+3B\kappa q_0-2L_2)}{L_2^{5/2}}
\]

was exactly zero. This statement holds symbolically in the parameters, not merely at the sampled ratios. It also gives the printed exact edge coefficients \(36\sqrt5/25\), \(-27\sqrt{10}/25600\), and \(-3\sqrt5/50\) in the new example. The logarithm sign repair does not change this coefficient.

However, the coefficient field must be specified before the eleven functions can be called a basis.

- Over **R(r) or C(r)** they are immediately dependent. The three constant atoms 1, zeta(2), zeta(3) are scalar multiples of 1; zeta(2)log(r) is a scalar multiple of log(r). There are at most eight independent functions over these fields.
- Over **Q(r)** those scalar multipliers are not automatically allowed. Independence of the eleven functions would already imply Q-linear independence of **1, zeta(2), zeta(3)**. That is an open period-independence problem, also explicitly acknowledged by this paper in its discussion following Theorem `thm:flat`. It cannot be concluded from the listed functional equations.
- Inspecting Landen and inversion separately and finding an outside argument in each does not rule out eliminating the outside functions by combining identities. It is not an independence proof.
- Numerical real matrix rank cannot establish Q(r)-independence. Nor does independence over Q(r), if eventually proved, automatically address the larger coefficient class actually used here, which includes parameter logarithms depending on r.

For an explicit numerical check, I evaluated the eleven functions at
`r = 1/4, 1/3, 2/5, 1/2, 2/3, 3/4, 3/2, 2, 3, 4, 5, 6`
at 130 dps. The sample matrix has eight numerically resolved singular values and three at about `1e-131`; the exact constant-column dependencies explain these. Its eight resolved singular values run from `15.5580662633` down to `8.59755308708e-7`. This refutes real-linear independence, not Q(r)-independence.

I then retained the eight analytic functions consisting of the first seven atoms and 1. To search for rational-function relations, I tested polynomial coefficient degrees zero, one, and two by sampling the functions r^j*f_i(r), with respectively 8, 16, and 24 columns. Clearing a denominator turns any bounded-degree rational relation into this kind of polynomial relation.

| coefficient degree | normalized determinant | smallest/largest singular value |
|---|---|---|
| 0 | 6.1540381190276228431e-19 | 5.7573651391134872718e-8 |
| 1 | -1.5266361814082992176e-96 | 4.9983303443787660725e-19 |
| 2 | -7.9306308243626242784e-238 | 8.1152485597717818636e-30 |

No relation was detected in these bounded searches. These degree bounds apply after clearing a common denominator, not separately to the original rational numerators and denominators. This is numerical evidence for the analytic part; it neither proves arbitrary-degree independence nor resolves the constant-period obstruction. As another limited structural check, the four weight-three atoms have independent top symbols: in the word columns `baa,bab,bba,bbb`, their printed integer coefficient matrix has determinant -6. This excludes a constant-coefficient relation among their top symbols; it is not being promoted to a proof over Q(r).

A defensible replacement is: **“In the specified reduction convention, retaining the two trilogarithms unevaluated, the explicit zeta(3) coefficient is Z_e.”** A basis-independent or uniquely numerical coefficient has not been established. At special ratios, further evaluations can transfer zeta(3) between the trilogarithms and the explicit constant, as K1e itself notes.

**Additional defects and limits to fix before paper insertion.**

1. **The harmonic atom is mislabeled.** `result.md:19` says A(x)=S_{1,2}(x), but `esum_engine.py` and `step_20.py` define
   \(A(x)=\sum_{n\ge1}H_nx^n/n^2\). With the paper's own Nielsen convention (`main.tex:1307`),
   \[
   A(x)=\operatorname{Li}_3(x)+S_{1,2}(x),\qquad
   S_{1,2}(x)=\tfrac12\int_0^x\frac{\log^2(1-t)}t\,dt.
   \]
   The printed series are \(A=x+3x^2/8+11x^3/54+25x^4/192+\cdots\) and \(S_{1,2}=x^2/4+x^3/6+11x^4/96+\cdots\). The identity actually used for A is correct; its identification with S_{1,2} is not. I checked the identity against the integral definition at x=1/2 and x=-2. It can also be proved by differentiation: after the dilogarithm Landen identity, the derivative of the stated A-reduction is \((\operatorname{Li}_2(x)+\log^2(1-x)/2)/x\), and both expressions for A have value zero at x=0. The symbolic derivative difference printed zero.
2. **The proposed paper paragraph swaps the two trilogarithm coefficients.** `result.md:61` lists Phi_A and 2Phi_A+Phi_L in the order Li3(1-r), Li3((r-1)/r). The correct order is **2Phi_A+Phi_L, Phi_A**, as its own earlier formulas state. I independently substituted the stated Phi formulas and obtained two exact zero differences.
3. **Numerator versus full coefficient is inconsistent.** In `result.md:9`, each displayed c has a further factor L2^(-5/2) outside the bracket. Thus the exact statement is c_zeta3/L2^(5/2)=Z_e, not c_zeta3=Z_e if c denotes the displayed numerator. The pickle evaluation handles the denominator correctly.
4. **Define the constant parameter logs explicitly.** The table uses lA0=log(A0/w_j0), lA1=log(A1/w_j1), and lB=log(B/w_jo). They are not the full t-dependent log(a_j(t)/w_j), despite the wording in `result.md:7`; the log(t), log(1-t), and log(1+kap*t) parts have already been integrated.
5. **Specify geometry, positivity, and branches.** The natural theorem concerns nondegenerate, strictly convex centrally symmetric hexagons with positive weights and the specified normalized family. Then L2>0, Delta>0, 0<B<1, A0=1-Br>0, and r>0. Both trilogarithm arguments are below one, so the real branches can be stated directly. Degenerate hexagons require a separate limiting statement. The result also depends on the chosen weighted p-family, not just the shape at p=1; weights cannot be silently absorbed into k for p other than one.
6. **The discarded generic test is not proof evidence.** The zero-before-log substitution bug described in K1e's note is present in `step_59.py:27-35`. Its third test also gives A0=-1 and therefore is not a physical hexagon-edge parameter set. I did not count it as a validation. The independent example above satisfies all positivity constraints.
7. **The repaired general formula has not been fully re-proved in this review.** Exact table/pickle agreement, the general Z_e identity, and the new polar test pass. I did not re-derive every inherited J-table entry or the full general d^3 engine. Thus this review supplies a counterexample to the displayed theorem and strong independent evidence for the repaired implementation, not a fresh complete proof of its inherited symbolic machinery.

**Reproduction.**

From this directory, use the specified Python executable with:

```text
review.py exact
review.py formula
review.py basis --dps 130
review.py numeric --dps 100 --h 0.001 --levels 6
review.py numeric --dps 120 --h 0.0005 --levels 7
```

The numeric mode integrates the polar definition first and compares the supplied closed form afterward. Assertions check geometry, normalization, the exact parity calculations, table/pickle agreement, and the corrected numerical formula. The evidence records Python 3.14.7, SymPy 1.14.0, and mpmath 1.3.0.
