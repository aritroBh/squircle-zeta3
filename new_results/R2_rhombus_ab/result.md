# R2: perimeter derivatives at the rhombus

**PROVED.** Put $r=\sqrt{a^2+b^2}$, $q=(a^2-b^2)/(a^2+b^2)$ and
$\kappa=a^2b^2/(a^2+b^2)^2=(1-q^2)/4$. For all positive semi-axes,

$$
\boxed{\begin{aligned}
P(1;a,b)&=4\sqrt{a^2+b^2},\\
P'(1;a,b)&=0,\\
P''(1;a,b)&=\frac{4\pi^2a^2b^2}{3(a^2+b^2)^{3/2}},\\
P'''(1;a,b)&=\frac{24a^2b^2}{(a^2+b^2)^{3/2}}
\left(\zeta(3)-\frac{\pi^2}{3}\right).
\end{aligned}}
$$

Equivalently, for $d=p-1$,

$$P(p;a,b)=4r+\frac{2\pi^2a^2b^2}{3r^3}d^2
+\frac{4a^2b^2}{r^3}\left(\zeta(3)-\frac{\pi^2}{3}\right)d^3+O(d^4).$$

There are no $\log(b/a)$ terms in these derivatives. In fact, no finite-order
coefficient requires them: reflection makes every odd logit moment vanish.
The logarithm of the aspect ratio does enter the divergence scale below.

## Exact representation and differentiation

**PROVED.** With $e=1/p$, $\eta=e-1$, $a=r\cos\phi$, $b=r\sin\phi$,
$0<\phi<\pi/2$, parametrize one quadrant by
$x=a t^e$, $y=b(1-t)^e$. Then

$$\sqrt{(x')^2+(y')^2}
=e\sqrt{a^2t^{2\eta}+b^2(1-t)^{2\eta}}
=er[t(1-t)]^{\eta/2}g_\phi(\eta l(t)),$$

$$l(t)=\log\frac{t}{1-t},\qquad
g_\phi(z)=\sqrt{\cos^2\phi\,e^z+\sin^2\phi\,e^{-z}}
=\sqrt{\cosh z+q\sinh z}.$$

Consequently the exact formula, valid for real $p>0$, is

$$\boxed{P(p;a,b)=4r(1+\eta)\int_0^1
[t(1-t)]^{\eta/2}g_\phi(\eta l(t))\,dt.}$$

The exponent outside $g_\phi$ is $\eta/2$, rather than $\eta$.
The positive square root is used on the real axis. Differentiation to every
fixed order at $\eta=0$ is justified by an integrable bound consisting of
endpoint powers with exponents greater than $-1$, times powers of endpoint
logarithms. This proves smoothness and finite-order Taylor remainders; it
does not assert analyticity.

Set $h(t)=\tfrac12\log[t(1-t)]$. Reflection sends $l\mapsto-l$ and fixes $h$,
so $\int_0^1h^jl^{2k+1}\,dt=0$. The required moments are:

| Moment on the unit interval | Exact value |
| --- | --- |
| $\int_0^1 h\,dt$ | $-1$ |
| $\int_0^1 h^2\,dt$ | $2 - \frac{\pi^{2}}{12}$ |
| $\int_0^1 l^2\,dt$ | $\frac{\pi^{2}}{3}$ |
| $\int_0^1 h^3\,dt$ | $-6 + \frac{3 \zeta\left(3\right)}{2} + \frac{\pi^{2}}{4}$ |
| $\int_0^1 h*l^2\,dt$ | $- (2 \zeta\left(3\right) + \frac{\pi^{2}}{3})$ |

The beta generating function gives these moments without endpoint
approximations:

$$M(\eta,s)=\int_0^1e^{\eta h+s l}\,dt
=B(1+\eta/2+s,1+\eta/2-s),$$

$$A(\eta,s)=(1+\eta)M(\eta,s)
=\frac{\Gamma(1+\eta/2+s)\Gamma(1+\eta/2-s)}{\Gamma(1+\eta)}.$$

Writing $g_\phi(z)=\sum_{k\ge0}c_kz^k$, its beginning is

$$g_\phi(z)=1+\frac q2z+\frac{2-q^2}{8}z^2
+\frac{3q^3-2q}{48}z^3+O(z^4).$$

The gamma expansion yields, through total degree three,

$$\log A=-\frac{\pi^2}{24}\eta^2+\frac{\pi^2}{6}s^2
+\frac{\zeta(3)}4\eta^3-\zeta(3)\eta s^2
+O((|\eta|+|s|)^4).$$

As a finite-order formal identity,

$$\frac{P}{4r}=\left.\sum_{k\ge0}c_k\eta^k\partial_s^kA(\eta,s)\right|_{s=0}
=1+\frac{(1-q^2)\pi^2}{24}\eta^2
-\frac{(1-q^2)\zeta(3)}4\eta^3+O(\eta^4).$$

Finally $\eta=-d/(1+d)=-d+d^2-d^3+O(d^4)$ gives the stated derivatives.
For $a=b=1$ this reduces exactly to the supplied check:

$$\frac{P(p;1,1)}{\sqrt2}=4+\frac{\pi^2}{6}d^2
+\left(\zeta(3)-\frac{\pi^2}{3}\right)d^3+O(d^4).$$

## Independent numerical evidence

**NUMERICAL.** All computations used the specified Python interpreter,
mpmath at **110 decimal digits**, and exact rational SymPy stencil weights.
These are high-precision agreement checks, not interval-certified error bounds.
For a residual $E=|v-v_{\rm ref}|$, scaled digits means
$\lfloor-\log_{10}(E/\max(1,|v_{\rm ref}|))\rfloor$, capped at **100**.
A zero residual is zero at working precision and receives that cap.
Displayed values use **38 significant digits**; residuals use **6**.

The independent direct quadrature uses the trigonometric parametrization
$x=a\cos^{2/p}\theta$, $y=b\sin^{2/p}\theta$:

$$P_{\rm direct}=\frac8p\int_0^{\pi/2}
\sqrt{a^2\cos^{4/p-2}\theta\sin^2\theta
+b^2\sin^{4/p-2}\theta\cos^2\theta}\,d\theta.$$

Both quadratures are reflected onto half their domains to avoid subtracting
nearly equal numbers at the upper endpoints.

| Axes $(a,b)$ | $p$ | Direct perimeter | Absolute residual | Scaled digits |
| --- | --- | --- | --- | --- |
| (1,2) | 0.9 | 8.9710300363393093561158166620244843352 | 0.0 | 100 |
| (1,2) | 1.1 | 8.9650396312803044604346243745557869974 | 6.92951e-104 | 100 |
| (1,2) | 1.5 | 9.261790263805834271400204784101887598 | 5.46855e-76 | 76 |
| (3,1) | 0.9 | 1.2670680235214890801384306190279323256e+1 | 0.0 | 100 |
| (3,1) | 1.1 | 1.2665807325437913653774449396539653838e+1 | 4.72669e-104 | 100 |
| (3,1) | 1.5 | 1.2932468323322938929397965210903684526e+1 | 2.19854e-76 | 76 |

Derivatives below use only $P_{\rm direct}$, with the centered **7-point**
stencil at offsets $-3,-2,-1,0,1,2,3$, computed at
$h=10^{-12}$ and $h/2=5\times10^{-13}$.
Orders one and two have truncation error $O(h^6)$; order three has
$O(h^4)$ error. Order zero is the direct integral at the center.
Both step sizes independently passed the **30-digit** assertion for every
entry. The first derivative is checked by absolute residual because its
exact value is zero. The same scale definition covers the other entries.

| Axes $(a,b)$ | Derivative order | Closed form evaluated | Direct finite difference, $h/2$ | Absolute residual, $h$ | Absolute residual, $h/2$ | Scaled digits, $h/2$ |
| --- | --- | --- | --- | --- | --- | --- |
| (1,1) | 0 | 5.6568542494923801952067548968387923143 | 5.6568542494923801952067548968387923143 | 0.0 | 0.0 | 100 |
| (1,1) | 1 | 0.0 | -4.3394038349612650413978965456007522628e-71 | 2.77722e-69 | 4.3394e-71 | 70 |
| (1,1) | 2 | 4.6525761330925863561050406241129368599 | 4.6525761330925863561050406241129368599 | 1.6323e-68 | 2.55047e-70 | 70 |
| (1,1) | 3 | -1.7715665747400515929969780068967686733e+1 | -1.7715665747400515929969780068967686733e+1 | 2.26806e-44 | 1.41754e-45 | 46 |
| (1,2) | 0 | 8.9442719099991587856366946749251049418 | 8.9442719099991587856366946749251049418 | 0.0 | 0.0 | 100 |
| (1,2) | 1 | 0.0 | 3.1127172295320500070457690846600722465e-71 | 1.99214e-69 | 3.11272e-71 | 70 |
| (1,2) | 2 | 4.7080760217316064662066878384170935177 | 4.7080760217316064662066878384170935177 | 1.28311e-68 | 2.00486e-70 | 70 |
| (1,2) | 3 | -1.7926993288964738099929772159666923676e+1 | -1.7926993288964738099929772159666923676e+1 | 1.62691e-44 | 1.01682e-45 | 46 |
| (2,5) | 0 | 2.1540659228538016125002841966161318225e+1 | 2.1540659228538016125002841966161318225e+1 | 0.0 | 0.0 | 100 |
| (2,5) | 1 | 0.0 | 6.6894936469711421217769887545971013177e-71 | 4.28128e-69 | 6.68949e-71 | 70 |
| (2,5) | 2 | 8.4263886295816453547595562721260142711 | 8.4263886295816453547595562721260142711 | 2.29057e-69 | 3.57902e-71 | 71 |
| (2,5) | 3 | -3.2085253448639068563085677554332038479e+1 | -3.2085253448639068563085677554332038479e+1 | 3.49638e-44 | 2.18523e-45 | 46 |
| (1,10) | 0 | 4.0199502484483561080877059651038304748e+1 | 4.0199502484483561080877059651038304748e+1 | 5.32245e-110 | 5.32245e-110 | 100 |
| (1,10) | 1 | 0.0 | -8.6041920953556046515413337637008659956e-72 | 5.50668e-70 | 8.60419e-72 | 71 |
| (1,10) | 2 | 1.2964519381840393447564738272605085961 | 1.2964519381840393447564738272605085961 | 1.07699e-69 | 1.6828e-71 | 70 |
| (1,10) | 3 | -4.936514425003379729372492786738101566 | -4.936514425003379729372492786738101566 | 4.49712e-45 | 2.8107e-46 | 46 |

## Literature comparison, offline

**PROVED relative to the statement supplied in the task.** The result
$P'(1)=0$ agrees with the quoted Montilla–Rodriguez–Linares statement
(arXiv:2607.04592, Corollary 3.5).

**UNVERIFIED literature match.** The task identifies Rodriguez–Montilla,
arXiv:2607.09048, Corollary 5.3 / Equation (104), as a source for $P''(1)$,
but supplies no expression from that equation. No paper text is present in
this directory, and no network access was used. I therefore cannot claim
an algebraic match to Equation (104). The independently derived candidate
for comparison is
$4\pi^2a^2b^2/[3(a^2+b^2)^{3/2}]$; the table verifies this formula numerically.

## Every coefficient and divergence

**PROVED: coefficient structure.** The exact generating identity is

$$\log A(\eta,s)=\sum_{m\ge2}\frac{(-1)^m\zeta(m)}m
\left[(\eta/2+s)^m+(\eta/2-s)^m-\eta^m\right].$$

For any fixed order $n$, truncate this identity to total degree $n$,
exponentiate formally, and apply
$\sum_{k=0}^n c_k\eta^k\partial_s^k|_{s=0}$.
Only even $k$ survive. Thus each coefficient in $\eta$ is a finite
polynomial, with rational coefficients, in the $c_{2k}$ and the ordinary
zeta values $\zeta(2),\ldots,\zeta(n)$; zeta monomials have total weight $n$.
Each $c_{2k}$ is itself an even polynomial in $q$ with rational coefficients.
The change to $d$ mixes these polynomials from lower orders. Multiplying
by $4r$ restores the perimeter. No explicit aspect-ratio logarithms occur
in these finite coefficients. All coefficients above the constant vanish
at $q=\pm1$, as the degenerate-axis limit has constant perimeter.

The infinite operator expression is a **formal asymptotic expansion**;
termwise integration of an infinite convergent $g_\phi$ series over the
entire interval would be unjustified.

**PROVED: the Taylor series has zero radius of convergence** for every
fixed $a,b>0$, in both $\eta$ and $d$. Here is a factorial-growth argument
which also rules out cancellation by the weight. The nearest branch points
of $g_\phi$ are

$$z=\log(b/a)\pm i\pi/2,\qquad
R=\sqrt{\log^2(b/a)+\pi^2/4}.$$

Let $I(\eta)=P/[4r(1+\eta)]$. With the logit coordinate and reflection,

$$I(\eta)=\int_0^\infty e^{-u}H(\eta u)K(\eta,u)\,du,$$
$$H(z)=e^{-z/2}[g_\phi(z)+g_\phi(-z)],\qquad
K(\eta,u)=(1+e^{-u})^{-2-\eta}.$$

The Taylor series $H(z)=\sum H_nz^n$ has radius exactly $R$.
For unequal axes the nearest branch points of the two summands are
distinct; for equal axes the summands coincide and add. The nonzero
exponential multiplier cannot remove their square-root singularities.

Split $K=1+(K-1)$. The first part gives formal coefficients $n!H_n$,
so its order-one Borel transform is exactly $H$. To bound the remainder,
write $K-1=\sum_{j\ge0}k_j(u)\eta^j$. For $u\ge0$,

$$|k_0(u)|\le2e^{-u},\qquad |k_j(u)|\le e^{-ju}/j!\quad(j\ge1).$$

For any $0<\rho<R$, choose $C_\rho$ with $|H_k|\le C_\rho\rho^{-k}$.
If $E_n$ denotes the remainder's coefficient of $\eta^n$, integration of
these bounds gives

$$\frac{|E_n|}{n!}\le\frac{C_\rho}{(2\rho)^n}
\left[1+\frac12\sum_{j=1}^\infty\frac{(2\rho)^j}{(j!)^2}\right].$$

Here one uses $(n-j)!/n!\le1/j!$ and $(j+1)^{n-j+1}\ge2^{n-j+1}$.
Thus the remainder's Borel transform is analytic for $|z|<2R$ and cannot
cancel the branch points of $H$ at distance $R$.
Writing $I\sim\sum I_n\eta^n$ therefore proves

$$\limsup_{n\to\infty}\left|\frac{I_n}{n!}\right|^{1/n}=\frac1R.$$

Multiplication by $1+\eta$ changes the Borel transform to itself plus its
primitive. A primitive has a higher branch exponent and cannot cancel the
leading square-root branch. The perimeter's Taylor series in $\eta$
therefore also has zero radius. Since $\eta=-d/(1+d)$ is locally analytic
and invertible, a convergent Taylor series in $d$ would imply one in
$\eta$, a contradiction. Finite-order Taylor expansions remain valid.

**NUMERICAL.** The corresponding Borel radii for the derivative checks are:

| Axes $(a,b)$ | $R$ |
| --- | --- |
| (1,1) | 1.5707963267948966192313216916397514421 |
| (1,2) | 1.7169315985765248532176365624549724701 |
| (2,5) | 1.8185130754522537619469604195225491799 |
| (1,10) | 2.7873462667474125047784184568663680487 |

No conjectural claim is needed for these results.

## Reproduction

Run from this directory:

```sh
"<project>/.venv/bin/python" r2_analysis.py
```

`r2_analysis.py` performs the symbolic identities, the independent
quadratures, and the finite differences; it asserts every requested
numerical agreement and regenerates this report. It prints all numerical
entries and then the complete report. `rhombus_formula.py` is the supplied
equal-axis reference and was read without modification; its finite
operator truncation is not used as an exact quadrature in these checks.
