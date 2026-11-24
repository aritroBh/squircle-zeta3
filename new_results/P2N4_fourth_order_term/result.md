# P2N4: the $n^{-4}$ term of Theorem 16 (`thm:P2asym`)

**Verdict: PROVED** (new argument, needs a referee). For $n\ge9$ and every $0<\theta<\pi$,
$E_n=\frac{d-2za_2\log n}{n^3}+\frac{Q_4}{n^4}+R_n$ with $R_n=-\frac{z a_2}{3}\frac{\log^3n}{n^5}+\Delta'_n$,
$|\Delta'_n|\le2.5\cdot10^5\,z\frac{\log^2n}{n^5}$, hence $|R_n|\le1.2\cdot10^5\,z\frac{\log^3n}{n^5}$.
The bound is uniform in $\theta\in(0,\pi)$ and carries the factor $z=\sin^2\frac\theta2$, so it vanishes as $\theta\to0$ and stays bounded as $\theta\to\pi$.
The power $\log^3n$ is the true one: $a_2>0$ at every angle, so $R_n$ is not $O(\log^2n/n^5)$ for any fixed $\theta$.

$$Q_4=z\Bigl[a_2\bigl((\lambda_n-3)^2-3\bigr)-2a_3(\lambda_n-3)+a_2^2+2a_4\Bigr]-\tfrac16z(1-z)(20z^2-2z+43),\qquad \lambda_n=\log n-a_1 .$$

Expanded in $\log n$, this is the notes'/referee's $Q_4$ (sympy difference $0$):
$Q_4=z\,a_2\log^2n-2z(a_1a_2+3a_2+a_3)\log n+z(a_1^2a_2+6a_1a_2+6a_2+2a_1a_3+6a_3+a_2^2+2a_4)+\frac{10}3z^4-\frac{11}3z^3+\frac{15}2z^2-\frac{43}6z$.

## Status, kept apart

- **PROVED** (below): the lemma, with explicit constants $2.5\cdot10^5$ and $1.2\cdot10^5$ for $n\ge9$. The proof uses only Theorem 15, Lemma 2.2 (`g3:lem:g`), Lemma `g3:lem:Lam`(c) of Section `g3:sec:P2`, Lemma 2.5 (`g3:lem:ranges`), the bound $|H|\le3.1$ from the proof of Lemma 2.6 (`g3:lem:taylor`), and the closed forms of $a_1,\dots,a_4$ (Lemmas `g3:lem:a123`, `g3:lem:a4`, needed only to make $Q_4$ explicit).
- **DERIVED** (exact computer algebra, `symbolic.py` / `symbolic.log`): $Q_4$ equals the formal $n^{-4}$ coefficient of $E_n$ computed by sympy from the exact formula of Theorem 15 (difference $0$). It also equals the referee's $Q_4$ (difference $0$). The pieces $P_1,\dots,P_6$ of the proof expand separately to the same $n^{-3}$ and $n^{-4}$ coefficients (difference $0$). The $n^{-5}$ coefficient has $\log^3n$-part $-\frac{za_2}3$, as part (a) of the lemma says. Every algebraic identity used for the remainders is checked there.
- **NUMERICAL** (`numeric.py` / `numeric.log`, mpmath, 90 digits): table below. The explicit constants are in `constants.py` / `constants.log` (exact rational arithmetic, rounding upward).

## LaTeX, to paste after Lemma 2.7 (`g3:prop:n3`)

The $\gapnote$ of Section `g3:sec:P2asym:a4` should then say that the $n^{-4}$ term is proved (Lemma `g3:prop:n4`); a replacement is at the end of this file.

```latex
\begin{lemma}[the $n^{-4}$ term]\label{g3:prop:n4}
Put $\lambda_n=\log n-a_1$, $p_4(z)=-\frac16z(1-z)(20z^2-2z+43)$ and
\[
Q_4(\theta,n)=z\Bigl[a_2\bigl((\lambda_n-3)^2-3\bigr)-2a_3(\lambda_n-3)+a_2^2+2a_4\Bigr]+p_4(z).
\]
For $n\ge9$ and $0<\theta<\pi$ let $R_n=E_n-\frac{d(\theta)-2za_2\log n}{n^3}-\frac{Q_4(\theta,n)}{n^4}$, with $d$ as in Lemma~\ref{g3:prop:n3}. Then
\[
\text{(a)}\quad R_n=-\frac{z\,a_2}{3}\,\frac{\log^3n}{n^5}+\Delta'_n,\quad|\Delta'_n|\le2.5\cdot10^5\,z\,\frac{\log^2n}{n^5};\qquad
\text{(b)}\quad |R_n|\le1.2\cdot10^5\,z\,\frac{\log^3n}{n^5}.
\]
Since $a_2>0$ (Lemma~\ref{g3:lem:ranges}), the power $\log^3n$ cannot be lowered at any angle.
\end{lemma}

In powers of $\log n$,
$Q_4=z\,a_2\log^2n-2z(a_1a_2+3a_2+a_3)\log n+z(a_1^2a_2+6a_1a_2+6a_2+2a_1a_3+6a_3+a_2^2+2a_4)+\frac{10}3z^4-\frac{11}3z^3+\frac{15}2z^2-\frac{43}6z$.
All of $a_1,\dots,a_4$ are explicit (Lemmas~\ref{g3:lem:a123} and~\ref{g3:lem:a4}), so $Q_4$ is explicit.

\begin{proof}
Let $n\ge9$ (real), $x=\frac1n\le\frac19$, $\ell_n=\log n\ge\log9>2.1972$, $\varepsilon=\frac x{1-x}=\frac1{n-1}\le\frac18$, and $q=e^{-\varepsilon\lambda_n}$ as in the proof of Lemma~\ref{g3:prop:unif}. We write $\lambda=\lambda_n$. Then $\ell_n<\lambda<\ell_n+1$, since $-1<a_1<0$.

\emph{Step 0: tools.}
(i) For $k\ge1$, $J\ge-1$ and $0\le x\le\frac19$, Taylor's formula with Lagrange remainder for $(1-x)^{-k}$ gives
\[
0\le(1-x)^{-k}-\sum_{j=0}^{J}\tbinom{k+j-1}{j}x^j\le\tbinom{k+J}{J+1}\bigl(\tfrac98\bigr)^{k+J+1}x^{J+1},
\]
because the $(J+1)$-st derivative is $(k)_{J+1}(1-x)^{-k-J-1}$, which is positive and at most $(k)_{J+1}(\frac98)^{k+J+1}$. Put $\beta=\frac98$, so $\beta^5<1.8021$ and $\beta^6<2.0273$. Since $\varepsilon^k=x^k(1-x)^{-k}$, this gives
\begin{gather*}
0\le r_2:=\varepsilon^2-x^2-2x^3-3x^4\le4\beta^5x^5,\qquad 0\le r_3:=\varepsilon^3-x^3-3x^4\le6\beta^5x^5,\\
0\le r_4:=\varepsilon^4-x^4\le4\beta^5x^5,\qquad 0\le\varepsilon^5-x^5\le5\beta^6x^6,\qquad \varepsilon^5\le\beta^5x^5,\quad \varepsilon^6\le\beta^6x^6 .
\end{gather*}
(ii) For $t\ge0$: $e^{-t}=1-t+\frac{t^2}2-\vartheta\frac{t^3}6$ with $\vartheta=e^{-\xi}$ for some $\xi\in(0,t)$, so $0<\vartheta\le1$ and $0\le1-\vartheta\le t$. Also $0\le e^{-t}-1+t\le\frac{t^2}2$ and $0\le1-e^{-t}\le t$. We use these with $t=\varepsilon\lambda$.
(iii) $x\lambda<\frac{\ell_n+1}n\le\frac{\log9+1}9<0.3553$ and $x\lambda^2<\frac{(\ell_n+1)^2}n\le\frac{(\log9+1)^2}9<1.1359$: both functions of $n$ decrease for $n\ge e$.
(iv) Coefficients of $H$. By the proof of Lemma~\ref{g3:lem:taylor}, $H$ is holomorphic in $|\varepsilon|<\frac14$ and $|H|\le3.1$ on $|\varepsilon|\le\frac14$. Since $\log H=\sum_{j\ge2}a_j\varepsilon^j$, we have $H=1+a_2\varepsilon^2+a_3\varepsilon^3+b_4\varepsilon^4+\cdots$ with $b_4=a_4+\frac12a_2^2$. Cauchy's estimate on the circle $|\varepsilon|=\frac15$ (strictly inside $|\varepsilon|<\frac14$) gives $|a_3|\le3.1\cdot5^3=387.5$ and $|b_4|\le3.1\cdot5^4=1937.5$, for every $z\in(0,1)$.
(v) Taylor to order $5$. For $\varepsilon_0\in[0,\frac18]$ the closed disc $|\varepsilon-\varepsilon_0|\le\frac1{10}$ lies in $|\varepsilon|\le0.225<\frac14$. Cauchy's estimate on its boundary gives $|H^{(5)}(\varepsilon_0)|\le5!\cdot3.1\cdot10^5$. $H$ is real and smooth on $[0,\frac18]$, so Taylor's formula with Lagrange remainder gives, for $0\le\varepsilon\le\frac18$,
\[
H(\varepsilon)=1+a_2\varepsilon^2+a_3\varepsilon^3+b_4\varepsilon^4+\rho_5(\varepsilon),\qquad|\rho_5(\varepsilon)|\le3.1\cdot10^5\,\varepsilon^5 .
\]

\emph{Step 1: splitting.} In the proof of Lemma~\ref{g3:prop:unif}, $E_n=T_1+T_2+T_3+T_4$ with $T_2=2zq\rho$ and $\rho=H-1-a_2\varepsilon^2=a_3\varepsilon^3+b_4\varepsilon^4+\rho_5$. So $E_n=P_1+\dots+P_6$ with
\begin{gather*}
P_1=T_1=2za_2(q\varepsilon^2-x^2),\quad P_2=2za_3q\varepsilon^3,\quad P_3=2zb_4q\varepsilon^4,\quad P_4=2zq\rho_5,\\
P_5=T_3=-\frac{4z(1-z)x^3}{1-2x},\qquad P_6=T_4=\sum_{m\ge2}(-1)^m\frac{x}{m-(m+1)x}\Lambda_m\,x^m .
\end{gather*}
Here we used $n^2(n-2)=n^3(1-2x)$ and $\frac{\varepsilon}{m-\varepsilon}=\frac{x}{m-(m+1)x}$. For each piece we write its part of order $x^3$ and $x^4$ and bound the rest.

\emph{Step 2: $P_1$.} By (ii), $q\varepsilon^2=\varepsilon^2-\lambda\varepsilon^3+\frac{\lambda^2}2\varepsilon^4-\vartheta\frac{\lambda^3}6\varepsilon^5$. With $r_2,r_3,r_4$ of (i),
\[
q\varepsilon^2-x^2=(2-\lambda)x^3+\bigl(3-3\lambda+\tfrac12\lambda^2\bigr)x^4+r_2-\lambda r_3+\tfrac12\lambda^2r_4-\vartheta\tfrac{\lambda^3}6\varepsilon^5 .
\]
We split the last term as $-\frac16\ell_n^3x^5-\frac16\bigl(\vartheta\varepsilon^5\lambda^3-\ell_n^3x^5\bigr)$. Now
\[
\bigl|\vartheta\varepsilon^5\lambda^3-\ell_n^3x^5\bigr|\le(\varepsilon^5-x^5)\lambda^3+x^5(\lambda^3-\ell_n^3)+(1-\vartheta)\varepsilon^5\lambda^3 .
\]
By (i) and (iii) the first term is at most $5\beta^6(x\lambda)\lambda^2x^5\le5\beta^6\cdot0.3553\,\lambda^2x^5$. The second is at most $3\lambda^2(\lambda-\ell_n)x^5\le3\lambda^2x^5$. The third is at most $\varepsilon\lambda\cdot\varepsilon^5\lambda^3\le\beta^6(x\lambda^2)\lambda^2x^5\le\beta^6\cdot1.1359\,\lambda^2x^5$. With $0<a_2<2.15$ (Lemma~\ref{g3:lem:ranges}),
\[
P_1=2za_2\Bigl[(2-\lambda)x^3+\bigl(3-3\lambda+\tfrac12\lambda^2\bigr)x^4\Bigr]-\frac{za_2}3\ell_n^3x^5+e_1,
\]
\[
|e_1|\le z\,x^5\Bigl[2\cdot2.15\,\beta^5\,(4+6\lambda+2\lambda^2)+\tfrac{2.15}3\bigl(5\cdot0.3553\,\beta^6+3+1.1359\,\beta^6\bigr)\lambda^2\Bigr].
\]

\emph{Step 3: $P_2$ and $P_3$.} By (ii), $q\varepsilon^3=\varepsilon^3-\lambda\varepsilon^4+\eta\varepsilon^5\lambda^2$ with $0\le\eta\le\frac12$. Hence
\[
P_2=2za_3\bigl[x^3+(3-\lambda)x^4\bigr]+e_2,\qquad|e_2|\le2z\cdot387.5\,x^5\bigl(6\beta^5+4\beta^5\lambda+\tfrac12\beta^5\lambda^2\bigr).
\]
Also $q\varepsilon^4-x^4=r_4-(1-q)\varepsilon^4$ with $0\le(1-q)\varepsilon^4\le\lambda\varepsilon^5$, so
\[
P_3=2zb_4x^4+e_3,\qquad|e_3|\le2z\cdot1937.5\,x^5\bigl(4\beta^5+\beta^5\lambda\bigr).
\]

\emph{Step 4: $P_4$ and $P_5$.} As $0<q\le1$, $|P_4|\le2z\cdot3.1\cdot10^5\,\beta^5x^5$. Next, $\frac1{1-2x}=1+2x+\frac{4x^2}{1-2x}$, so $P_5=-4z(1-z)(x^3+2x^4)+e_5$ with $|e_5|\le\frac{16z}{1-2/9}x^5=\frac{144}7z\,x^5$.

\emph{Step 5: $P_6$.} We use $|\Lambda_m|\le2z$ (Lemma~\ref{g3:lem:Lam}(c)). For $m=2$: $\frac x{2-3x}=\frac x2+\frac{3x^2}4+\frac{9x^3/8}{1-3x/2}$, so the term $m=2$ is $\frac12\Lambda_2x^3+\frac34\Lambda_2x^4$ plus at most $2z\cdot\frac98\cdot\frac65x^5=2z\cdot1.35\,x^5$. For $m=3$: $\frac1{3-4x}=\frac13+\frac{4x}{3(3-4x)}$, so the term $m=3$ is $-\frac13\Lambda_3x^4$ plus at most $2z\cdot\frac{4}{3(3-4/9)}x^5=2z\cdot\frac{12}{23}x^5$. For $m\ge4$: $m-(m+1)x\ge\frac{8m-1}9\ge\frac{31}9$, so these terms add up to at most $2z\cdot\frac9{31}\sum_{m\ge4}x^{m+1}\le2z\cdot\frac9{31}\beta\,x^5$. From $\Pleg_2=\frac{3c^2-1}2$, $\Pleg_4=\frac{35c^4-30c^2+3}8$ and $c=1-2z$, Lemma~\ref{g3:lem:Lam}(c) gives $\Lambda_3=\frac{\Pleg_2(c)-\Pleg_4(c)}7=2z(1-z)(1-5z+5z^2)$. So
\[
P_6=\tfrac12\Lambda_2x^3+\bigl(\tfrac34\Lambda_2-\tfrac13\Lambda_3\bigr)x^4+e_6,\qquad|e_6|\le2z\bigl(1.35+\tfrac{12}{23}+\tfrac{81}{248}\bigr)x^5 .
\]

\emph{Step 6: the coefficients.} The $x^3$ terms of Steps 2--5 add up to $2za_2(2-\lambda)+2za_3-4z(1-z)+\frac12\Lambda_2=d-2za_2\ell_n$, as in Lemma~\ref{g3:prop:n3}. The $x^4$ terms are
\[
2za_2\bigl(3-3\lambda+\tfrac12\lambda^2\bigr)+2za_3(3-\lambda)+2zb_4-8z(1-z)+\tfrac34\Lambda_2-\tfrac13\Lambda_3 .
\]
Here $2a_2(3-3\lambda+\frac12\lambda^2)=a_2\bigl((\lambda-3)^2-3\bigr)$, $2a_3(3-\lambda)=-2a_3(\lambda-3)$, $2b_4=2a_4+a_2^2$, and
\[
-8z(1-z)+\tfrac32z(1-z)(1-2z)-\tfrac23z(1-z)(1-5z+5z^2)=-\tfrac16z(1-z)(20z^2-2z+43)=p_4(z).
\]
So the $x^4$ terms add up to $Q_4$, and
\[
R_n=-\frac{za_2}3\ell_n^3x^5+\Delta'_n,\qquad\Delta'_n=e_1+e_2+e_3+P_4+e_5+e_6 .
\]

\emph{Step 7: the constants.} Collecting Steps 2--5, $|\Delta'_n|\le z\,x^5\,(C_0+C_1\lambda+C_2\lambda^2)$ with $C_0\le1.154\cdot10^6$, $C_1\le1.262\cdot10^4$ and $C_2\le721$. The term $P_4$ dominates: it gives $2\cdot3.1\cdot10^5\beta^5<1.118\cdot10^6$ of $C_0$. Now $\lambda<\ell_n+1\le\kappa\ell_n$ with $\kappa=1+\frac1{2.1972}<1.4552$, and $\ell_n\ge2.1972$. So $C_0\le C_0\frac{\ell_n^2}{2.1972^2}$, $C_1\lambda\le C_1\kappa\frac{\ell_n^2}{2.1972}$ and $C_2\lambda^2\le C_2\kappa^2\ell_n^2$. Hence $|\Delta'_n|\le C'z\,\ell_n^2x^5$ with $C'<2.49\cdot10^5$. This is (a). For (b), $|R_n|\le\bigl(\frac{2.15}3+\frac{C'}{2.1972}\bigr)z\,\ell_n^3x^5<1.14\cdot10^5\,z\,\ell_n^3x^5$.

Finally, (a) shows $R_n\,n^5/\log^3n\to-\frac{za_2}3\ne0$ as $n\to\infty$, for every fixed $\theta$. So no bound $O(\log^2n/n^5)$ holds.
\end{proof}
```

Remarks on the proof (not for pasting):
- The radii. Every Cauchy estimate uses a circle strictly inside $|\varepsilon|<\frac14$. The circle $|\varepsilon|=\frac15$ is used for the coefficients $a_3,b_4$, and the circles $|\varepsilon-\varepsilon_0|=\frac1{10}$ ($\varepsilon_0\le\frac18$, so they lie in $|\varepsilon|\le0.225$) for $H^{(5)}$, exactly as in Lemmas 2.6 and 2.7.
- Edge cases. $n$ may be any real $\ge9$; the removable points $n=1+\frac1j$ of Theorem 15 are $<2$. Every remainder carries the factor $z$: in $P_1$ to $P_4$ through the prefactor $2z$, in $P_5$ through $z(1-z)$, and in $P_6$ through $|\Lambda_m|\le1-c=2z$. The bounds on $a_2,a_3,b_4$ do not depend on $z$. So the estimate is uniform on $(0,\pi)$, it tends to $0$ like $\sin^2\frac\theta2$ as $\theta\to0$, and it stays bounded as $\theta\to\pi$ ($z\to1$). We do not need $Q_4$ to have a limit at the ends; it does, since $a_1,\dots,a_4$ extend continuously to $z\in[0,1]$.
- Constants. $P_4$ contributes $1.12\cdot10^6$ out of $1.15\cdot10^6$ in $C_0$, and it comes only from the Cauchy bound for $H^{(5)}$ on circles of radius $\frac1{10}$. The true constant is about $1.4$ (table), so the proven one is pessimistic by a factor of about $10^5$. Lemma 2.6 (factor 1400) behaves the same way. A Cauchy remainder on one circle $|\varepsilon|=0.24$ would lower $P_4$ by a factor of about $40$. We kept the supplement's circles so the proof matches Lemma 2.7 line by line.
- Lemma 2.7 states its $C_3$ without a value. Steps 0--5 above also give an explicit $C_3$ if one stops one order earlier; we did not write it out.

## Numerical table (NUMERICAL)

`numeric.py`, mpmath at 90 digits. $L_\theta(n)$ comes from the exact formula of Theorem 15. The series in $m$ is summed until the terms fall below $10^{-95}$. $a_1,\dots,a_4$ come from the closed forms; they agree with the Taylor coefficients of $\log g$ from a Cauchy integral to $1.6\cdot10^{-89}$. The exact formula agrees with direct quadrature of the arc length to $7\cdot10^{-48}$ at $n=50,400$, $\theta=0.7,2.4$.

col1 $=R_n\,n^5/(z\log^3n)$ (lemma (b): $|\cdot|\le1.2\cdot10^5$; limit $-a_2/3$).
col2 $=(R_n+\frac{za_2}3\frac{\log^3n}{n^5})\,n^5/(z\log^2n)$ (lemma (a): $|\cdot|\le2.5\cdot10^5$; formal limit $a_1a_2+4a_2+a_3$, from the $n^{-5}$ coefficient in `symbolic.log`).

| $\theta$ | $z$ | col1 at $n=9$ | $10^2$ | $10^4$ | $10^6$ | $10^8$ | limit $-a_2/3$ | col2 at $10^8$ | formal limit |
|---|---|---|---|---|---|---|---|---|---|
| 0.01 | 2.5e-5 | -1.395 | -0.0866 | -0.0931 | -0.1899 | -0.2582 | -0.5483 | 5.344 | 6.580 |
| 0.3 | 0.0223 | -1.374 | -0.0844 | -0.0929 | -0.1892 | -0.2570 | -0.5446 | 5.297 | 6.516 |
| 1 | 0.2298 | -1.144 | -0.0606 | -0.0897 | -0.1813 | -0.2441 | -0.5058 | 4.820 | 5.875 |
| $\pi/2$ | 0.5 | -0.788 | -0.0262 | -0.0839 | -0.1667 | -0.2212 | -0.4408 | 4.045 | 4.851 |
| 2 | 0.7081 | -0.479 | -0.0006 | -0.0775 | -0.1494 | -0.1948 | -0.3709 | 3.244 | 3.817 |
| 2.8 | 0.9711 | 0.0153 | 0.0254 | -0.0605 | -0.1029 | -0.1268 | -0.2091 | 1.518 | 1.668 |
| 3.13 | 0.99998 | 0.0896 | 0.0243 | -0.0552 | -0.0889 | -0.1071 | -0.1668 | 1.099 | 1.168 |

Over the 7 angles and $n\in\{9,10,20,50,10^2,\dots,10^8\}$: $\sup|\text{col1}|=1.40$ and $\sup|\text{col2}|=5.34$. Both are bounded, as the lemma says. col1 moves toward $-a_2/3$ like $\text{col2}/\log n$, which is slow ($\log10^8=18.4$). col2 moves toward $a_1a_2+4a_2+a_3$ with corrections of order $1/\log n$. This is consistent with $\log^3n$ being the right power. The angles $0.01$ and $3.13$ show the uniformity at both ends: the scaled residual (divided by $z$) changes little as $z\to0$ or $z\to1$.

## Files

- `symbolic.py`, `symbolic.log`: sympy. The exact formal expansion of $E_n$ to $x^5$ from Theorem 15, with $Q_4$ (compact form and referee form) and $d$ (difference $0$). It also checks the per-piece expansions of the proof, the $L^3$ and $L^2$ parts of the $x^5$ coefficient, $\Lambda_1,\Lambda_2,\Lambda_3$, $p_4$, and the identities of Steps 1, 4 and 5.
- `constants.py`, `constants.log`: the constants $C_0,C_1,C_2,C',C_4$ of Step 7 in exact rational arithmetic ($C'=248\,840$, $C_4=113\,254$).
- `numeric.py`, `numeric.log`: the table above, the check of $a_1,\dots,a_4$, and quadrature against the exact formula.

## Suggested replacement for the gapnote of Section `g3:sec:P2asym:a4`

```latex
\gapnote{The last sentence of Theorem~16 is now covered: the $n^{-3}$ term is Lemma~\ref{g3:prop:n3}, $a_4$ is Lemma~\ref{g3:lem:a4}, and the $n^{-4}$ term is Lemma~\ref{g3:prop:n4}, with remainder $-\frac{za_2}3\frac{\log^3n}{n^5}+O\bigl(z\frac{\log^2n}{n^5}\bigr)$ uniformly in $\theta$. Lemma~\ref{g3:prop:n4} is a new argument and needs a referee.}
```
