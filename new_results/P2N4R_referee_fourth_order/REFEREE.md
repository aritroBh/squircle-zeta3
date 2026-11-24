# P2N4R: referee report on `08_theory/P2N4_fourth_order_term/result.md` (the $n^{-4}$ term of thm:P2asym)

**Verdict: CONFIRMED.** The proof of Lemma `g3:prop:n4` is correct as written, including its constants. It needs no mathematical repair. Three editorial repairs follow (E1–E3). They change notation and the surrounding text of the supplement, not the argument.

**The printed sentence of thm:P2asym** ("The terms of order $n^{-3}$ and $n^{-4}$, and the weight-four coefficient $a_4$, are explicit too (in our notes).") **is fully backed by written proofs once Lemma `g3:prop:n4` is pasted into `sec_g3_P2asym.tex` and the two places that still call the $n^{-4}$ term "only formal" are updated (E2).** The three parts are:
- the $n^{-3}$ term: Lemma `g3:prop:n3` (refereed in `referee_g3`). Its constant $C_3$ is now explicit as a corollary, $C_3\le3.1\cdot10^4$ (E3).
- $a_4$: Lemma `g3:lem:a4` (refereed in `referee_g3`, Q1). I re-checked its polynomial identity (`ref_extra.log`, difference 0) and its value numerically (below).
- the $n^{-4}$ term: this lemma.

Until the lemma is pasted in, the supplement file still says the $n^{-4}$ term is only formal, so the sentence is not yet backed by the supplement as it stands.

## Files (all mine, written without reusing P2N4 code)

| script | what | log |
|---|---|---|
| `ref_symbolic.py` | sympy. Starts from the exact formula of thm:P2, with Legendre polynomials from `sympy.legendre` and $g=\exp(\sum a_j\varepsilon^j)$ (the definition of the $a_j$). It computes the $x^0..x^5$ coefficients of $E_n$, compares them with $d$, with both forms of $Q_4$, with $p_4$ and with the Step-6 identity, and expands the pieces $P_1..P_6$ to $x^5$. | `ref_symbolic.log` |
| `ref_constants.py` | Every constant of Steps 0–7, in exact rationals and exact sympy reals. It also re-checks the 3.1 chain of Lemmas 2.2 and 2.6. | `ref_constants.log` |
| `ref_numeric.py` | mpmath at 110 digits. It evaluates the exact formula, checks it against arc-length quadrature, and computes $a_1..a_5$ both by Cauchy integral (256 nodes) and from the closed forms. It computes the residuals at 7 angles and 12 values of $n$ from 9 to $10^8$ (non-integer $n=9.5$ included), and the true $\sup|H|$. | `ref_numeric.log` |
| `ref_extra.py` | The $a_4$ identity of Lemma `g3:lem:a4`, and an explicit $C_3$ for Lemma `g3:prop:n3`. | `ref_extra.log` |

## Line-by-line check

**Reused lemmas, with their hypotheses and ranges.**
- Theorem 15 holds for all real $n>1$. For $n\ge9$ we have $\varepsilon\le\frac18$, so no removable points arise.
- Lemma `g3:lem:g` (b, c): $|F|\le2.1$ and $|\pi\varepsilon/\sin\pi\varepsilon|\le1.12$ on $|\varepsilon|\le\frac14$, uniformly in $z\in[0,1]$. Recomputed: $|F|\le2.0768$, $1.1186$, so $|g|\le2.3231\le2.4$.
- Lemma `g3:lem:taylor`: $|H|\le2.4e^{1/4}=3.0817<3.1$ on $|\varepsilon|\le\frac14$. This uses $|a_1|<1$, from Lemma `g3:lem:ranges`.
- Lemma `g3:lem:Lam`(c): $|\Lambda_m|\le1-c=2z$.
- Lemma `g3:lem:ranges`: $-1<a_1<0$ and $0<a_2<2.15$, for all $z\in(0,1)$.

All are quoted with the right hypotheses. $H$ is holomorphic on the open disc, continuous on the closed disc, and real on the real axis, so Lagrange's form of the remainder applies. As a sanity check, the true $\max|H|$ on $|\varepsilon|=\frac14$ over a $z$-grid is $1.111$.

**Step 0.**
- (i) The Lagrange constants are $\binom{k+J}{J+1}$, with factors $\beta^{k+J+1}$, $\beta=\frac98$. They are $4,6,4,5$, with $\beta^5,\beta^5,\beta^5,\beta^6$. Correct (`ref_symbolic.log`). The case $J=-1$ (empty sum) gives $\varepsilon^5\le\beta^5x^5$ and $\varepsilon^6\le\beta^6x^6$.
- (ii) Correct.
- (iii) $(\log n+1)/n$ and $(\log n+1)^2/n$ are decreasing for $n>e$; the derivatives were checked. Their values at $n=9$ are $0.35525<0.3553$ and $1.13581<1.1359$.
- (iv) $H=\exp(\sum_{j\ge2}a_j\varepsilon^j)$ has coefficients $1,0,a_2,a_3,a_4+\frac12a_2^2,\dots$ (sympy). The Cauchy estimates on $|\varepsilon|=\frac15<\frac14$ give $387.5$ and $1937.5$.
- (v) The discs $|\varepsilon-\varepsilon_0|\le\frac1{10}$, $\varepsilon_0\le\frac18$, lie in $|\varepsilon|\le0.225<\frac14$. So $|H^{(5)}|\le5!\cdot3.1\cdot10^5$, and since $H'(0)=0$ the degree-4 Taylor polynomial is the one written. Correct.

**Taylor expansion of $q$ to third order (Step 2).** The identity $q\varepsilon^2-x^2=(2-\lambda)x^3+(3-3\lambda+\frac12\lambda^2)x^4+r_2-\lambda r_3+\frac12\lambda^2r_4-\vartheta\frac{\lambda^3}6\varepsilon^5$ is exact. The splitting
$\vartheta\varepsilon^5\lambda^3-\ell^3x^5=(\varepsilon^5-x^5)\lambda^3+x^5(\lambda^3-\ell^3)-(1-\vartheta)\varepsilon^5\lambda^3$ is exact. The three bounds use $\lambda-\log n=-a_1\in(0,1)$, $1-\vartheta\le\varepsilon\lambda$ and (iii), and are correct. The extracted term $-\frac{za_2}3\log^3n\,x^5$ is exactly the $L^3$ part of the $x^5$ coefficient of $E_n$ (sympy, my own expansion). $P_1$ is the only piece that contains $L^3$.

**Steps 3–4.**
- $q\varepsilon^3=\varepsilon^3-\lambda\varepsilon^4+\eta\lambda^2\varepsilon^5$ with $0\le\eta\le\frac12$. Correct.
- $q\varepsilon^4-x^4=r_4-(1-q)\varepsilon^4$ with $0\le(1-q)\le\varepsilon\lambda$. Correct.
- $|P_4|\le2z\cdot3.1\cdot10^5\beta^5x^5$, since $0<q\le1$. Correct.
- $P_5$: $\frac1{1-2x}=1+2x+\frac{4x^2}{1-2x}$, so $|e_5|\le\frac{144}7zx^5$. Correct.

**Step 5 ($T_4$).**
- The $m=2$ remainder is $\frac{9/8}{1-3x/2}\le1.35$, and the $m=3$ remainder is $\frac{4}{3(3-4x)}\le\frac{12}{23}$. The identities were checked symbolically and the constants exactly.
- For $m\ge4$: $m-(m+1)x\ge\frac{8m-1}9\ge\frac{31}9$, and $\sum_{m\ge4}x^{m+1}\le\beta x^5$. So the bound is $2z\cdot\frac{81}{248}x^5$. Correct.
- $\Lambda_3=2z(1-z)(1-5z+5z^2)$, from `sympy.legendre`.

**Step 6.** $p_4=-8z(1-z)+\frac34\Lambda_2-\frac13\Lambda_3=-\frac16z(1-z)(20z^2-2z+43)$, difference 0. The pieces $P_1..P_6$, expanded formally, reproduce the full $x^3$, $x^4$ **and $x^5$** coefficients of $E_n$ (difference 0). So no term is lost in the splitting.

**$Q_4$.** My own expansion of the exact formula gives the following, all with difference **exactly 0**:
- $x^4$ coefficient minus $Q_4$ (compact form) $=0$;
- $x^4$ coefficient minus $Q_4$ (expanded form) $=0$;
- $x^3$ coefficient minus $(d-2za_2L)=0$.

$Q_4$ does not involve $a_5$ or $a_6$. For reference, the full $x^5$ coefficient is $-\frac{za_2}3L^3+z(a_1a_2+4a_2+a_3)L^2-z(a_1^2a_2+8a_1a_2+2a_1a_3+a_2^2+12a_2+8a_3+2a_4)L+(\dots)$. The constant term is in `ref_symbolic.log`.

**Step 7: constants, in exact rational arithmetic** (`ref_constants.log`).
- $C_0=1\,153\,627.05\le1.154\cdot10^6$, $C_1=12\,615.67\le1.262\cdot10^4$ and $C_2=720.17\le721$.
- $\kappa=1+1/2.1972=1.45512<1.4552$, and $\log9\ge2.1972$.
- $C'=C_0/L_0^2+C_1\kappa/L_0+C_2\kappa^2=248\,840.34<2.49\cdot10^5$. Using the printed roundings instead gives $248\,922.8$, still $<2.49\cdot10^5$.
- $C_4=\frac{2.15}3+C'/L_0=113\,254.1<1.14\cdot10^5\le1.2\cdot10^5$.

**From $\Delta'$ to (b).** $|R_n|\le\frac{za_2}3\log^3n\,x^5+C'z\log^2n\,x^5\le(\frac{2.15}3+\frac{C'}{\log9})z\log^3n\,x^5$. This needs $\log n\ge\log9$ and $a_2<2.15$, and both hold for every $\theta$.
- $a_2$ is bounded uniformly on $(0,\pi)$: $0<a_2<\zeta(2)+\frac12$, with limits $a_2\to\zeta(2)$ as $\theta\to0$ and $a_2\to\frac12$ as $\theta\to\pi$. My numbers: $1.6449339$ at $\theta=10^{-3}$, $0.5000041$ at $\pi-10^{-3}$.
- None of the bounds uses a lower bound on $z$ or on $1-z$. The bounds on $|a_3|$ and $|b_4|$ are Cauchy bounds independent of $z$.
- Every error term carries $2z$: through the prefactor in $P_1..P_4$, through $z(1-z)\le z$ in $P_5$, and through $|\Lambda_m|\le2z$ in $P_6$.

So the estimate is uniform on all of $(0,\pi)$. It vanishes like $z$ as $\theta\to0$ and stays bounded as $\theta\to\pi$.

**$\log^3n$ is the true power.** $a_2>0$ at every angle (Lemma `g3:lem:ranges`). By (a), $R_nn^5/\log^3n=-\frac{za_2}3+O(z/\log n)\to-\frac{za_2}3\ne0$. Correct.

## Numerical check (`ref_numeric.log`, 110 digits, own code)

**Inputs.**
- The exact formula agrees with arc-length quadrature to $5\cdot10^{-46}$ at $\theta=0.4$ and $2.9$, for $n=9$ and $13.5$.
- The closed forms of $a_1..a_4$ agree with Cauchy-integral Taylor coefficients of $\log g$ to $\le4\cdot10^{-105}$ at every angle.

**Angles.** $\theta\in\{10^{-3},0.3,1,\frac\pi2,2.5,\pi-0.01,\pi-10^{-3}\}$, so $z$ runs from $2.5\cdot10^{-7}$ to $1-2.5\cdot10^{-7}$.

**Values of $n$.** $9,9.5,10,20,50,10^2,\dots,10^8$.

**Columns.**
- col1 $=R_nn^5/(z\log^3n)$. Bound in (b): $1.2\cdot10^5$. Limit: $-a_2/3$.
- col2 $=(R_n+\frac{za_2}3\frac{\log^3n}{n^5})n^5/(z\log^2n)$. Bound in (a): $2.5\cdot10^5$.
- col3 $=(R_nn^5-Q_5(\log n))\,n/(z\log^4n)$, where $Q_5$ is the full $x^5$ coefficient, with $a_5$ from the Cauchy integral.

**Results.**
- $\sup|\text{col1}|=1.395$ and $\sup|\text{col2}|=5.344$, both attained at $\theta=10^{-3}$. These are far inside the proved bounds, which are pessimistic by about $10^5$ because of the Cauchy bound for $H^{(5)}$.
- col1 moves monotonically toward $-a_2/3$ for $n\ge10^3$. At $10^8$ it has reached about half of the limit ($-0.258$ against $-0.548$ at $\theta=10^{-3}$, and $-0.107$ against $-0.167$ near $\pi$). The gap shrinks like $\text{col2}/\log n$, which is what (a) predicts.
- col3 stays bounded ($\sup|\text{col3}|=1.25$, and $\le0.035$ for $n\ge10^3$). This confirms that $Q_5$, with its $-\frac{za_2}3\log^3n$ leading part, is the next term, so col2 tends to its formal limit $a_1a_2+4a_2+a_3$ at rate $1/\log n$.
- My values match the P2N4 table digit for digit wherever the grids overlap.

## Repairs (editorial; the proof itself stands)

**E1 (notation clash).** `sec_g3_P2asym.tex` already uses $\ell=\log(1-z)$ and $R_\theta(\varepsilon)$. The proof's $\ell_n=\log n$ and its $R_n$ invite confusion. Write $L_n$ for $\log n$, and $E^{(4)}_n$ for the residual:
```latex
% in Lemma g3:prop:n4 and its proof: replace \ell_n by L_n:=\log n, and R_n by E^{(4)}_n, e.g.
For $n\ge9$ and $0<\theta<\pi$ let $E^{(4)}_n=E_n-\frac{d(\theta)-2za_2\log n}{n^3}-\frac{Q_4(\theta,n)}{n^4}$ \dots
\[\text{(a)}\quad E^{(4)}_n=-\frac{z\,a_2}{3}\,\frac{\log^3n}{n^5}+\Delta'_n,\ |\Delta'_n|\le2.5\cdot10^5\,z\,\frac{\log^2n}{n^5};\qquad
\text{(b)}\quad |E^{(4)}_n|\le1.2\cdot10^5\,z\,\frac{\log^3n}{n^5}.\]
```

**E2 (the supplement must stop calling the $n^{-4}$ term formal in two places, not one).** `result.md` replaces only the gapnote. The statement section of `sec_g3_P2asym.tex` also says "the $n^{-4}$ term is only formal (Section~\ref{g3:sec:P2asym:a4})", and the Plan paragraph omits the lemma. Replace those as follows.
```latex
The last sentence of the theorem is about terms that go beyond the displayed expansion. We prove the $n^{-3}$ term (Lemma~\ref{g3:prop:n3}), $a_4$ (Lemma~\ref{g3:lem:a4}) and the $n^{-4}$ term (Lemma~\ref{g3:prop:n4}).
```
Add "and Lemma~\ref{g3:prop:n4} the $n^{-4}$ term" to the Plan. The gapnote should drop "needs a referee": the lemma was refereed once, by `08_theory/P2N4R_referee_fourth_order/REFEREE.md` (verdict CONFIRMED). Use:
```latex
\gapnote{The last sentence of Theorem~16 is covered: the $n^{-3}$ term is Lemma~\ref{g3:prop:n3}, $a_4$ is Lemma~\ref{g3:lem:a4}, and the $n^{-4}$ term is Lemma~\ref{g3:prop:n4}, with remainder $-\frac{za_2}3\frac{\log^3n}{n^5}+O\bigl(z\frac{\log^2n}{n^5}\bigr)$ uniformly in $\theta$. Lemma~\ref{g3:prop:n4} is a new argument, refereed once (\path{08_theory/P2N4R_referee_fourth_order/REFEREE.md}, confirmed).}
```

**E3 (optional; closes G4 of `referee_g3`, an explicit $C_3$).** Add after Lemma `g3:prop:n4`:
```latex
\begin{corollary}\label{g3:cor:C3}
Lemma~\ref{g3:prop:n3} holds with $C_3=3.1\cdot10^4$.
\end{corollary}
\begin{proof}
$E_n-\frac{d-2za_2\log n}{n^3}=\frac{Q_4}{n^4}+E^{(4)}_n$. By Step~0(iv), $|a_3|\le387.5$ and $|a_2^2+2a_4|=2|b_4|\le3875$; also $0<a_2<2.15$, $\lambda_n\in(\log n,\log n+1)$ and $|p_4(z)|\le\frac{43}6z$. Hence $|Q_4|\le z\bigl(2.15\max(3,(\lambda_n-3)^2-3)+775|\lambda_n-3|+3875+\frac{43}6\bigr)\le935\,z\log^2n$ for $\log n\ge\log9$. By Lemma~\ref{g3:prop:n4}(b) and $\frac{\log n}n\le\frac{\log9}9$, $|E^{(4)}_n|\le1.2\cdot10^5\cdot\frac{\log9}9\,z\frac{\log^2n}{n^4}\le29\,300\,z\frac{\log^2n}{n^4}$. Add.
\end{proof}
```
The constant 935 is a numerical supremum over a grid of $L=\log n\in[\log9,202]$ (`ref_extra.log`). The ratio decreases to $2.15$ as $L\to\infty$, so the grid captures the maximum, which is at $L=\log9$. At that point the bound is $(6.45+622.2+3875+7.17)/4.828=934.2$.

**Minor (no change needed).** "$J\ge-1$" in Step 0(i) uses the empty-sum convention. This is fine, but a reader may prefer the two cases $\varepsilon^5\le\beta^5x^5$ and $\varepsilon^6\le\beta^6x^6$ stated directly, as they already are in the display.
