# Referee report: eight Table II rows of `paper/ieee/main.tex` (draft 37)

Rows as named in `paper/ieee/README_verification_map.md` (draft 36): R36, R47, R84, R96, R99, R100, R101, R106.
In draft 37 they are rows 36, 47, 83, 95, 98, 99, 100, 105 of Table II. From the octagon row on, each row moved up by one.
None of the eight is stated in the protected text (front matter, Sections I–VII, the first three Conclusion paragraphs, the first two Future Work paragraphs, the Acknowledgment). All are stated in Sections IX and XII.

All code here is my own. No project script was imported. Every script has a `.log` next to it.
Python: `.venv/bin/python` (mpmath, sympy). The status words are judged by the paper's own paragraph "The words in the table" (Sec. `sec:verify`).

## Summary

| row (map / draft 37) | verdict | main independent check |
|---|---|---|
| R36 / 36: $M_1(\frac{2\pi}3)$ | CONFIRMED AS STATED | printed form = Prop. `prop:M1short` at $c=\frac12$ to 121 digits; = corner integral to 82 digits; hand reduction exact |
| R47 / 47: $M_2(\alpha)$, classical $\Li_3$ | CONFIRMED AS STATED (one optional wording fix) | classical form = own HPL series of Thm. `thm:M2alpha` to 51 digits, = corner integral to 41 digits, 5 angles |
| R84 / 83: regular octagon $P'''(1)$ | STATUS WORD TOO STRONG: should be "Proved$^\dagger$ (computer-assisted)" | the 34-term table = own polar-integral $P'''(1)$ to 41 digits; finite differences to 7 digits; $C_4$, the $\zeta(3)$ numbers, the $A$-identity, $\mathcal D$ and its swap all confirmed |
| R96 / 95: Apostol's integral | CONFIRMED AS STATED | sum = quadrature; three closed values to 40 digits; 1.7330025902 confirmed by two representations (40 digits); 0.2003353, 0.2003409 reproduced; no $\varepsilon^4$ term |
| R99 / 98: squircle sine | CONFIRMED AFTER REPAIRS (one sentence) | $b_n'$ law, $b_1''$, $b_n''$ ($n\le17$), $N$, $N''(2)$, $N'''(2)$ to 40 digits |
| R100 / 99: squircle and polar | CONFIRMED AS STATED | duality to 50 digits; $c_1,\dots,c_6$ to 36–43 digits; $a_1,\dots,a_9$ to 44–50 digits |
| R101 / 100: even path $\lvert x\rvert^s$ | CONFIRMED AS STATED (cosmetic: $s>-1$) | asymptotics; $\sum\tilde b_n(2)$ three ways to 35 digits; cross-parity sum exact by hand |
| R106 / 105: Beukers' integral | CONFIRMED AS STATED | sum = nested quadrature; $B(1)$, $B(2)$ to 60 digits; $\varepsilon^2$ coefficient $-\pi^2/12$ to 61 digits; $\varepsilon^3,\varepsilon^4$ and the vanishing $\varepsilon^5$ term confirmed |

Lighter treatment, said once here: R84 (I did not re-derive the 730-term edge integrand, the IBP engine or $D_1$–$D_{10}$; I checked the final value instead). R47 (I did not re-run the ten HPL→$\Li_3$ derivative proofs one by one; I checked them together numerically). R106 (the "no rational approximations" sentence is Numerical and I did not test it). R99 and R100 (I read the proofs of the $b_n'$ law and of the weight bound but did not redo them).

---

## R36 (draft 37 row 36): $M_1(\frac{2\pi}3)$ in four constants

**(a) Claim.** Table row: "$M_1(\frac{2\pi}3)$ in four constants & Proved$^\dagger$ (Prop. `prop:M1short` at $c=\frac12$) & exact; 189 digits as a check".
Text (Sec. IX-B, after Thm. `thm:Malpha`): "at $\alpha=\frac{2\pi}3$ Proposition `prop:M1short` below and four classical dilogarithm identities (the value of $\Li_2(\frac12)$, duplication at $\frac12$, reflection at $\frac14$ and Landen at $-\frac12$) give exactly $24M_1(\tfrac{2\pi}{3})=\pi^2+48\lambda^2-6\log^23-24\lambda\log3-12\Li_2(\tfrac13)$".

**(b) Numbers** (`r36_m1.py`).
- I computed $M_1(\alpha)$ from the corner integral by differentiating under the integral sign. A central difference of $M(u;\alpha)$ itself confirms the derivative formula to $10^{-18}$.
- Prop. `prop:M1short` agrees with the integral to $\le2\cdot10^{-60}$ at $\alpha=0.3,\frac\pi3,\frac\pi2,\frac{2\pi}3,2.6,3.0$.
- The printed $2\pi/3$ form minus Prop. `prop:M1short` at $c=\frac12$ is $2.7\cdot10^{-121}$, and minus the integral (80 digits) it is $3\cdot10^{-82}$.
- $M_1(\frac{2\pi}3)=0.12579567893797323375386819768\ldots$
- At $c=1/\sqrt2$ the Proposition reproduces Thm. `thm:M1` to $7\cdot10^{-61}$, and $8M_1=2.020952769656358356947\ldots$ as printed.

**Hand reduction (exact).** At $c=\frac12$, $M_0=2\lambda-\frac32\log3$. Use $\Li_2(\frac12)=\frac{\zeta(2)}2-\frac{\lambda^2}2$, reflection at $\frac14$, and duplication $\Li_2(-\frac12)=\frac12\Li_2(\frac14)-\Li_2(\frac12)$. This gives $24M_1=60\lambda^2-36\lambda\log3+6\Li_2(\frac14)$, with the $\zeta(2)$ cancelling. Then duplication again and Landen, $\Li_2(-\frac12)=-\Li_2(\frac13)-\frac12\log^2\frac32$, give exactly the printed form.

**(c) Status.** The reduction from the Proposition is exact and done by hand. The Proposition itself is proved for all $n$ in `research/10_frontier/G1_zeta3_any_angle/verify.md` §3.4: a harmonic-sum identity on the Taylor coefficients, plus the proved $\mu_{1,n}$ formula. "Proved$^\dagger$" fits.

**Notes pointer (map, not paper).** The verification map gives `D10c_corner_general/verify_proof.md` as the proof of `prop:M1short` and R36. That file proves the 19+8-dilogarithm Theorem `thm:Malpha`, not the four-dilogarithm form. The proof of the four-dilogarithm form is G1 `verify.md` §3.4. That proof was written by the G1 verifier, so it has had no second, independent reading; my numerical check is the only one.

**(d)** No over-claim.

**Verdict: CONFIRMED AS STATED.**

---

## R47 (row 47): $M_2(\alpha)$ with classical $\Li_3$ at eight arguments

**(a) Claim.** Row: "$M_2(\alpha)$ with classical $\Li_3$ at eight arguments & Proved$^\dagger$ & exact (34 checks); 120 at 3 angles".
Text (after Thm. `thm:M2alpha` and the leading-coefficient paragraph): "each of the ten in $Y(c)+Y(-c)$ equals a combination of classical polylogarithms of weight three … $\Li_3$ at the eight arguments $c$, $-c$, $1-c$, $\frac1{1+c}$, $\frac{1+c}2$, $\pm\frac{1-c}{1+c}$ and $\frac{2c}{1+c}$ … that of $\zeta(3)$ being $\frac74+2c$ (in our notes) … So the classical form has the status of Theorem `thm:M2alpha`; it agrees with $M_2(\alpha)$ to 120 digits at $\alpha=\frac\pi3,\frac\pi2,\frac{2\pi}3$."

**(b) Numbers** (`r47_m2.py`).
- I computed $M_0,\dots,M_3$ from the corner integral by a Cauchy integral in $u$ (64 points, $|u|=0.3$), without differentiating anything.
- I wrote my own HPL series (letters $0,\pm1$) for the paper's displayed $Y(c)$.
- I typed the classical form B2′ by hand from `G1_zeta3_any_angle/result.md`. It has the eight $\Li_3$ arguments listed in the paper, $\Li_2$ only at $c,-c,\frac{1+c}2$, and the $\zeta(3)$ coefficient $\frac74+2c$.

At $\alpha=\frac\pi3,1.0,\frac\pi2,\frac{2\pi}3,2.6$:
- the classical form minus the HPL form is $\le7\cdot10^{-51}$;
- each of them minus the corner integral is $\le3\cdot10^{-41}$.

$M_0$ is reproduced to $\le4\cdot10^{-41}$. $M_2(\frac\pi2)=-0.0221619110872693994237549377\ldots$

**(c) Status.** `COR3k_classical_proof/result.md` proves the ten identities by exact derivative cancellation in $\mathbb Q(x)$, with exact constants at $0^+$ (34/34 items). Thm. `thm:M2alpha` is proved for all $n$; per the LEDGER this is in `research/10_frontier/W1_corner_family/B_taylor_pi4`, a location the map does not cite for this row. "Proved$^\dagger$, inheriting Thm. `thm:M2alpha`" fits.

**(d) Optional wording.** A coefficient of $\zeta(3)$ depends on which $\Li_3$ values are kept, as G1 `verify.md` §3.2 also notes. Suggest: "with coefficients in $\mathbb Q(c)$; in this form the coefficient of $\zeta(3)$ is $\frac74+2c$".

**Verdict: CONFIRMED AS STATED.**

---

## R84 (row 83): regular octagon, $P'''(1)$ closed form; $\mathcal D$ and its swap reduce

**(a) Claim.** Row: "Regular octagon: $P'''(1)$ closed form; $\mathcal D$, swap reduce (Sec. IX-C) & Proved$^\dagger$ & exact; about 410 (check)".
Text (Sec. IX-C): "$P'''(1)=48\,u\,S$, where $S$ is an explicit combination of 34 terms … So the closed form is proved$^\dagger$: the assembly is exact, and the $A$-atom identities and the reductions of $\mathcal D$ and its swap are proved. … its coefficient in $P'''(1)$ is $\bigl(\frac{87}2-\frac{117\sqrt2}4\bigr)u$. The local part is exactly $C_4$ … The rest, $\frac{39\sqrt2-42}4\,u$, is the non-local part".

**(b) Numbers** (`r84_octagon.py`, `r84_closed.py`, `r84_fd.py`).
- I wrote my own polar-integral engine with truncated Taylor arithmetic in $d=p-1$. As a test it reproduces the diamond: $P(1)=4\sqrt2$, $P'(1)=0$, $P''(1)=\pi^2\sqrt2/3$, and $P'''(1)=6\sqrt2(\zeta(3)-\pi^2/3)$, all to 30 digits.
- Octagon: $P(1)=8(\sqrt2-1)\sqrt{2-\sqrt2}$ and the printed $P'(1)$ both agree to 35 digits.
- $P'''(1)=-4.394496305276052159979110113806099913822$ (own value).
- The 34 nonzero terms of `OCT6_octagon_closed/result.md`, typed by hand, give the same number to $1.3\cdot10^{-41}$.
- A plain finite difference of $P(p)$, with no termwise differentiation ($h=0.004,0.002,0.001$, Richardson), gives $-4.3944963$.
- Exact arithmetic: $48(\frac{29}{32}-\frac{39\sqrt2}{64})=\frac{87}2-\frac{117\sqrt2}4$. $C_4$ from Thm. `thm:Cm` equals $(54-39\sqrt2)u$ to 25 digits. The difference is $\frac{39\sqrt2-42}4$.
- The printed identity for $A(x)-A(1-x)$ holds to $10^{-45}$ at $x=0.3,0.5,0.7,1-\frac{\sqrt2}2$. At $x=\frac12$ both sides vanish by antisymmetry, which pins the constant exactly.
- $\mathcal D=\Li_{2,1}(1-\sqrt2,-\frac{\sqrt2}2)$ and its swap, summed directly from the printed series definition, equal PRV1's single-$\Li_3$ forms to $\le7\cdot10^{-46}$.

**(c) Status.**
- Both depth-two reductions are proved twice: by PRV1 (exact derivative identities and exact endpoint pins) and by OCT_R34 (a level-8 certificate).
- The $A$-identity is proved by differentiation, and the constant is fixed by symmetry.
- The rest of the proof is a large exact computer computation, which I did not re-run: the 730-term edge integrand reduced to 83 monomials (OCT5/OCT6), the IBP/partial-fraction engine for the 20 mixed integrals, and the closed forms of $D_1$–$D_{10}$. The choice of the 34-term basis also used 18 lattice-found classical $\Li_2$/$\Li_3$ relations, which were re-checked numerically (745 digits) rather than derived symbolically. They are standard Landen and duplication instances.

By the paper's own definition, "Computer-assisted means that the proof rests on an exact computation that we did not reduce to steps one can check by hand", this proof is computer-assisted. The decagon row, proved by the same pipeline, already says "Proved (computer-assisted)". So the octagon's plain "Proved$^\dagger$" is too strong and inconsistent with its neighbour.

**(d)** Replace "So the closed form is proved$^\dagger$: the assembly is exact, …" with "So the closed form is proved$^\dagger$, computer-assisted: the assembly is an exact computer computation, and the $A$-atom identities and the reductions of $\mathcal D$ and its swap are proved." In Table II write "Proved$^\dagger$ (computer-assisted)".

**Verdict: STATUS WORD TOO STRONG. It should be "Proved$^\dagger$ (computer-assisted)". All printed formulas and decimals are correct.**

---

## R96 (row 95): Apostol's integral over the squircle

**(a) Claim.** Row: "the sum, $Z(\frac12)$, $Z(1)$, $Z(2)$; $\zeta(2)-Z(p)$ through $\varepsilon^4$ & Proved$^\dagger$ & over 100 (sum vs quadrature, 12 values of $p$); closed values: residual $\le4.7\cdot10^{-97}$; $\varepsilon^3$ coefficient $0.2003409$ at $p=100$".
Text (Sec. XII-D): the sum formula; $Z(\frac12)=\frac{\pi^2}9-4\log^2\varphi$, $Z(1)=\frac{\zeta(2)}3$, $Z(2)=\frac{7\zeta(2)}{12}$; $\zeta(2)-Z(p)=\varepsilon\int_0^\infty\frac{1-G}{x^2}-\frac{\zeta(2)}2\varepsilon^2+\frac{\zeta(3)}6\varepsilon^3+O(\varepsilon^5)$; "The integral is $1.7330025902\ldots$"; the checks $0.2003353$ ($p=50$) and $0.2003409$ ($p=100$).

**(b) Numbers** (`r96_apostol.py`, `r96_integral.py`).
- The sum equals a one-dimensional quadrature of the double integral at $p=\frac12,1,2,3,7.5$, to working precision (40 digits).
- The three closed values agree to $\le3\cdot10^{-42}$.
- $\int_0^\infty(1-G)/x^2=1.73300259020287278417012270424573195251$. It agrees with the different representation $\int_0^\infty\log\frac{u-\log(1-e^{-u})}u\,du$ of `AV1_apostol_squircle` to 30 digits.
  - Caution for anyone re-checking: a naive quadrature of $(1-G)/x^2$ is wrong in the fourth digit (I got 1.73333 at first), because $1-G$ cancels catastrophically at the tanh-sinh nodes near 0. Use the series of $\log G$ there, as `r96_integral.py` does.
- $(\zeta(2)-Z-\varepsilon I+\frac{\zeta(2)}2\varepsilon^2)/\varepsilon^3=0.200335262$ ($p=50$) and $0.2003409281$ ($p=100$), as printed.
- $(r-\frac{\zeta(3)}6)p^2\to-0.018891$, which is the predicted $\varepsilon^5$ coefficient $-[x^5]G\cdot\zeta(-3)=(\zeta(2)\zeta(3)-3\zeta(5))/60$. So the $\varepsilon^4$ term is 0 and the error is $O(\varepsilon^5)$, as stated.

**(c) Status.** I re-derived the Mellin argument. $f=1-G$ has fundamental strip $-2<\operatorname{Re}s<0$. Shifting left picks up the pole of $\zeta(2+s)$ at $s=-1$ (giving $F(-1)\varepsilon$; there is no double pole because $[x]G=0$) and the poles $-g_k\zeta(2-k)\varepsilon^k$. With $G=1-\zeta(2)x^2+2\zeta(3)x^3+\cdots$ this gives $-\frac{\zeta(2)}2$, $\frac{\zeta(3)}6$ and 0. The closed values are proved in AV1 by elementary steps (Lemma chain for $Z(2)$; golden-ratio dilogarithms for $Z(\frac12)$). "Proved$^\dagger$" fits.

**(d)** None.

**Verdict: CONFIRMED AS STATED.**

---

## R99 (row 98): the squircle sine near the circle

**(a) Claim.** Row: "$b_n'(2)$, $b_1''(2)$: Proved$^\dagger$; $N'''(2)$: Proved; $b_n''(2)$, $3\le n\le17$, and the $\frac3{64}\pi^2\zeta(3)$ in $b_1'''(2)$: Numerical".
Text (Sec. XII-F, "Other Deformations"):
- $b_1'(2)=\frac{\log2-1}2$ and $b_n'(2)=(-1)^{(n-1)/2}/((n^2-1)(2\lfloor\frac{n-1}4\rfloor+1))$;
- the norm $N(p)$ is pure Gamma, and "$N'''(2)$ contains $\frac{87}{16}\zeta(3)$ (proved)";
- "Beyond Gamma, $b_1'''(2)$ contains, besides the $\frac{87}{32}\zeta(3)$ from $N'''$, a term $\frac3{64}\pi^2\zeta(3)$ from the products $b_n'b_n''$ (numerical; one more series in $b_1'''$ we have not evaluated)."

**(b) Numbers** (`r99_sine.py`, `r99_more.py`, `r99_R.py`, `r99_R75.py`).
- My own engine uses $\Phi_p(y)=\frac\pi2\,y\,{}_2F_1(\frac1p,\frac1p;1+\frac1p;y^p)\,p\sin(\pi/p)/\pi$, with $p$-derivatives by `mpmath.taylor` and quadrature in $y=\sin\phi$.
- $b_n'(2)$ matches the law for $n=3,\dots,11$ to $\le10^{-28}$ (25-digit run) and $\le5\cdot10^{-44}$ (40-digit run).
- $b_1'(2)$ is exact to working precision.
- $b_1''(2)$ equals the notes' $\frac38-\log2+\frac{\log^22}4+\frac{5\pi^2}{96}-\frac{\pi^4}{768}$ to $6\cdot10^{-42}$.
- The notes' $b_n''$ forms for $n=3,\dots,17$ agree to $\le3\cdot10^{-42}$. Their $\pi^2$ and $\log2$ coefficients follow the $C_n$ and $b_n'$ laws exactly.
- I derived $N(p)$ independently: $\frac4\pi\int_0^{\pi/2}S_p^2=\frac{2\sin(\pi/p)}\pi\frac{\Gamma(3/p)\Gamma(1-1/p)}{\Gamma(1+2/p)}$, matching the paper's formula.
- The notes' $N''(2)$ and $N'''(2)=\frac1{16}[87\zeta(3)+\cdots]$ agree to $\le2\cdot10^{-40}$.
- $b_1'''(2)=-0.3026228510415431767747785\ldots$ as in the notes.
- Parseval at orders 1–3 holds up to the truncation of the sum at $n\le61$.
- Given the $C_n$ law, $\sum_{n\ge3}b_n'C_n=(1-\zeta(3))/64$ to 30 digits.

**(c) Status.**
- The $b_n'$ law is proved in AV8 (exact integrals, both residue classes symbolic).
- $b_1''$ is proved by Parseval from the proved $\sum(b_n')^2$.
- $N'''$ is exact Gamma calculus.
- These words fit.

**(d) Over-claim.** "$b_1'''(2)$ contains … a term $\frac3{64}\pi^2\zeta(3)$". What is established is only this: the $\pi^2$ part $C_n\pi^2$ of $b_n''$, which itself is known numerically for $n\le17$ only, contributes $\frac3{64}\pi^2\zeta(3)$ through $-3\sum b_n'C_n\pi^2$. The remaining series $R=\sum_{n\ge3}b_n'A_n$ is unknown. I isolated it: $R=-0.0092253434023847181165512739452\ldots$ (same 40 digits at 50 and 75 digits of working precision). PSLQ on $\{1,\log2,\log^22,\log^32,\pi^2,\pi^2\log2,\zeta(3),\pi^2\zeta(3),\pi^4,\pi^4\log2,G,\pi G\}$ finds only relations that change completely between 44 and 66 digits, so they are spurious. $R$ may therefore carry $\pi^2\zeta(3)$ and change the $\frac3{64}$. Proposed wording: "Beyond Gamma, the products $b_n'b_n''$ contribute to $b_1'''(2)$ a term $\frac3{64}\pi^2\zeta(3)$ through the $\pi^2$ parts of the $b_n''$ (numerical, from the second-order law checked for $n\le17$); one more series in $b_1'''$, which we have not evaluated, could change this coefficient." In the table: "the $\frac3{64}\pi^2\zeta(3)$ contributed to $b_1'''(2)$ by the $\pi^2$ parts: Numerical".

**Verdict: CONFIRMED AFTER REPAIRS (the one sentence above).**

---

## R100 (row 99): the squircle and its polar

**(a) Claim.** Row: "$V(\eta)=V(1/\eta)$; only powers of $\pi$ at the circle through $\delta^8$; $\zeta(2k+1)$ first at order $2k+2$ at the square end, $k\le5$ & Duality and odd orders: Proved; weight bound: Proved$^\dagger$; circle end and first orders: Numerical & $c_n$ ($n\le12$): 87–102; circle end: 38–48; duality: 50–69".
Text: the duality; evenness in $\delta$; "$c_n$ lie in $\mathbb Q[\log2,\zeta(2),\zeta(3),\zeta(5),\dots]$ with weight at most $n$ (proved…)"; "$c_4=\zeta(3)-\frac{\pi^2}3-\frac{\pi^2\log2}6+\cdots$".

**(b) Numbers** (`r100_polar.py`).
- I derived $V(K_p,K_{p'})=\frac12\int h_{K_{p'}}\,dS_{K_p}$, with $h_{K_{p'}}(u)=\lVert u\rVert_p$, as $V=\frac2p\int_0^1[t^{1-p}+(1-t)^{1-p}]^{1/p}dt$. For computing I use the smooth form $\frac4p\int_{\log2}^\infty e^{-x/p}(1+(e^x-1)^{1-p})^{1/p}dx$.
- $V(2)=\pi$ and $V(1)=4$.
- Duality holds to $\le5\cdot10^{-51}$ at $\eta=0.1,0.3,0.5$.
- Square end, Cauchy integral on $|\eta|=0.02$: $c_1,\dots,c_6$ equal the notes' exact forms (including the paper's $c_4$) to $10^{-43}$–$10^{-36}$.
  - Caution: at complex $\eta$ the integrand has branch points at $|\eta|\approx\pi/x$, so a large Cauchy radius gives wrong coefficients. The cut-off explained in the script handles this.
- Circle end, $|\delta|=0.15$: $a_1,a_3,\dots,a_9$ are $\le10^{-43}$, and $a_2,\dots,a_8$ equal the notes' $\pi\,\mathbb Q[\pi^2]$ forms to $\le6\cdot10^{-45}$.

**(c) Status.** Duality and odd orders are proved (a mixed area is symmetric). I read the weight bound (AV9 Q1/Q3: termwise Beta-derivative expansion, top-weight theorem) but did not redo it; "Proved$^\dagger$" is reasonable. In AV9 the circle-end coefficients through $\delta^{10}$ and the square-end coefficients through order 12 are exact symbolic computations, so "Exact (computation)" would also be honest. "Numerical" understates; it does not over-claim.

**(d)** None required. Optional: the circle end is exact through $\delta^{10}$ in the notes, while the text says $\delta^8$.

**Verdict: CONFIRMED AS STATED.**

---

## R101 (row 100): the even path $|x|^s$

**(a) Claim.** Row: "no $1/n^3$ Fourier mode at any real $s$; $\zeta(3)$ in odd-sector sums such as $\sum_n\tilde b_n(2)=2\pi\log2-\frac{7\zeta(3)}\pi$ & Proved$^\dagger$ & 60; …80 (three methods); $\tilde b_3(\frac52)$: 20 (independent)".
Text: the asymptotic formula for $a_n(s)$; the $1/n^3$ argument; $\tilde b_n(2)\ni-\frac4\pi(1-(-1)^n)n^{-3}$; $\sum a_n(1)\tilde b_n(1)=-\frac{7\zeta(3)}\pi$.

**(b) Numbers** (`r101_evenpath.py`).
- Residual of $a_n(s)$ minus the corner term minus three boundary terms ($c_1,c_2,c_3$ from integration by parts), times $n^8$: $O(1)$ at $s=2.5,1.3,2.0001$, $n=200,201$. So only even boundary powers occur. At $s=2$ the residual is $10^{-12}$.
- $a_7(2)=-4/49$ exactly.
- The closed form of $\tilde b_n(2)$ holds at $n=3,4,9$.
- $\sum\tilde b_n(2)$ by series, by $\frac1\pi\int_0^\pi t^2\cot\frac t2$, and by the closed form: all $1.6767860083890118652503058147891649$.
- By hand: $a_n(1)\tilde b_n(1)=-\frac8{\pi n^3}$ for odd $n$ and 0 for even $n$, so the sum is $-\frac{7\zeta(3)}\pi$ (also numerically, 35 digits).
  - Caution: `nsum` gives a wrong value, $-2.6721$, if fed the even/odd pattern directly.
- $\tilde b_3(\frac52)=3.523084206853680186260234$ (no printed value to compare).

**(c) Status.** Standard endpoint asymptotics, and the sums are exact by hand. "Proved$^\dagger$" fits.

**(d) Cosmetic.** "at any real $s$" should read "at any real $s>-1$" (for $s\le-1$ the integral diverges).

**Verdict: CONFIRMED AS STATED.**

---

## R106 (row 105): Beukers' integral over the squircle

**(a) Claim.** Row: "the sum for $B(p)$; $B(p)$ through $\varepsilon^4$; $B(1)=2\zeta(3)/3$; $B(2)=\zeta(3)/12+\pi G_{\mathrm C}/3$; no rational approximations at finite $p$ & Sum, $\varepsilon^2$, $\varepsilon^3$, $\varepsilon^4$ terms, $B(1)$, $B(2)$: Proved$^\dagger$; approximations: Numerical & sum vs quadrature: $7\cdot10^{-109}$; $\varepsilon^2$: 101; $B(1)$: 70; $B(2)$: 60; $\varepsilon^3$: $0.6005$ at $p=400$".
Text: $B=\sum\phi(n/p)/n^3$; $B=\zeta(3)-\frac{\pi^2}{12}\varepsilon^2+\frac{\zeta(3)}2\varepsilon^3-\frac{3\zeta(4)}{16}\varepsilon^4+\cdots$; $\phi=1-\zeta(3)w^3+\frac94\zeta(4)w^4+\cdots$; $\int_0^\infty(\phi-1)w^{-3}dw=-\pi^2/12$.

**(b) Numbers** (`r106_beukers.py`).
- The sum equals a nested quadrature of the double integral at $p=1,2,3$, to working precision (25 digits).
- $B(1)-\frac23\zeta(3)=8\cdot10^{-62}$ and $B(2)-(\frac{\zeta(3)}{12}+\frac{\pi G}3)=0$ at 60 digits.
- My own exact series for $\phi$ (from $\log\Gamma$ and $\psi$) gives $\phi_1=\phi_2=0$, $\phi_3=-\zeta(3)$ and $\phi_4=\frac94\zeta(4)$.
- $\int_0^\infty(\phi-1)w^{-3}dw=-\pi^2/12$ to $1.6\cdot10^{-61}$.
- The remainder after the $\varepsilon^4$ term, divided by $\varepsilon^5$, is $3.6\cdot10^{-4},1.8\cdot10^{-4},8.9\cdot10^{-5}$ at $p=100,200,400$. That is $\varepsilon\cdot\phi_6\zeta(-3)$, so the $\varepsilon^5$ term is 0 and the next one is $\phi_6\zeta(-3)\varepsilon^6$, as the Mellin argument predicts.
- $(B-\zeta(3)+\frac{\pi^2}{12}\varepsilon^2)/\varepsilon^3=0.60052$ at $p=400$, as printed ($0.6005$).

**(c) Status.** `AV2P_beukers_proofs/result.md` gives complete proofs:
- the Mellin expansion, with the reason for the missing $\varepsilon$ and $\varepsilon^2\log\varepsilon$ terms;
- $C_2=-\pi^2/12$ through a second (Watson) representation;
- $B(1)$ by the $T+U$ Clausen cancellation;
- $B(2)$ by the even/odd split and Clausen identities.

I read these and re-derived the Mellin step. "Proved$^\dagger$" fits. I did not test the "approximations" claim (Numerical).

**(d)** None required. Optional: "$+\cdots$" could be "$+O(\varepsilon^6)$", since $\zeta(-2)=0$ removes the $\varepsilon^5$ term. The text's "101 digits" for the $\varepsilon^2$ quadrature is conservative (the notes have 131).

**Verdict: CONFIRMED AS STATED.**

---

## Files

`r36_m1.py`, `r47_m2.py`, `r84_octagon.py` (engine), `r84_closed.py`, `r84_fd.py`, `r96_apostol.py`, `r96_integral.py`, `r99_sine.py` (engine), `r99_more.py`, `r99_R.py` and `r99_R75.py` (the unevaluated series $R$; spurious PSLQ at two precisions), `r100_polar.py`, `r101_evenpath.py`, `r106_beukers.py`. Each has a `.log`.
