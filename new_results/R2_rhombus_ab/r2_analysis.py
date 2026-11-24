"""Offline symbolic derivation, independent quadrature, and result.md for R2.

Run with <project>/.venv/bin/python r2_analysis.py
"""

from pathlib import Path

import sympy as sp
from mpmath import mp


def symbolic_check():
    eta, s, z, d, q = sp.symbols("eta s z d q", real=True)
    g = sp.series(sp.sqrt(sp.cosh(z) + q * sp.sinh(z)), z, 0, 4).removeO()
    log_a = sp.expand(sum(
        (-1) ** m * sp.zeta(m) / m
        * ((eta / 2 + s) ** m + (eta / 2 - s) ** m - eta ** m)
        for m in range(2, 4)
    ))
    # Through total degree three, exp(log_a) = 1 + log_a.
    beta_a = 1 + log_a
    normalized = sp.expand(sum(
        g.coeff(z, k) * eta**k
        * sp.diff(beta_a, s, k).subs(s, 0)
        for k in range(4)
    ))
    expected_eta = 1 + (1 - q**2) * (
        sp.pi**2 * eta**2 / 24 - sp.zeta(3) * eta**3 / 4
    )
    assert sp.simplify(normalized - expected_eta) == 0
    in_d = sp.series(normalized.subs(eta, -d / (1 + d)), d, 0, 4).removeO()
    expected_d = 1 + (1 - q**2) * (
        sp.pi**2 * d**2 / 24
        + (sp.zeta(3) / 4 - sp.pi**2 / 12) * d**3
    )
    assert sp.simplify(in_d - expected_d) == 0

    beta_m = sp.Poly(sp.expand(beta_a * (1 - eta + eta**2 - eta**3)), eta, s)
    beta_m = sum(c * eta**i * s**j for (i, j), c in beta_m.terms() if i + j <= 3)
    moments = {
        "integral h": sp.diff(beta_m, eta).subs({eta: 0, s: 0}),
        "integral h^2": sp.diff(beta_m, eta, 2).subs({eta: 0, s: 0}),
        "integral l^2": sp.diff(beta_m, s, 2).subs({eta: 0, s: 0}),
        "integral h^3": sp.diff(beta_m, eta, 3).subs({eta: 0, s: 0}),
        "integral h*l^2": sp.diff(beta_m, eta, s, 2).subs({eta: 0, s: 0}),
    }
    a, b = sp.symbols("a b", positive=True)
    r = sp.sqrt(a*a + b*b)
    series = (4 * r * in_d).subs(q, (a*a - b*b) / (a*a + b*b))
    derivatives = [sp.factor(sp.diff(series, d, n).subs(d, 0)) for n in range(4)]
    expected = [
        4*r,
        sp.S.Zero,
        4*sp.pi**2*a*a*b*b/(3*r**3),
        24*a*a*b*b/r**3*(sp.zeta(3) - sp.pi**2/3),
    ]
    assert all(sp.simplify(x-y) == 0 for x, y in zip(derivatives, expected))
    equal_axes = sp.simplify(series.subs({a: 1, b: 1}) / sp.sqrt(2))
    assert sp.simplify(equal_axes - (
        4 + sp.pi**2*d**2/6 + (sp.zeta(3)-sp.pi**2/3)*d**3
    )) == 0
    print("PROVED: symbolic checks passed", flush=True)
    print("g(z) =", g)
    print("log A through total degree 3 =", log_a)
    print("P/(4r), eta series =", sp.factor(normalized))
    print("P/(4r), d series =", sp.factor(in_d))
    for name, value in moments.items():
        print(name, "=", value)
    for n, value in enumerate(derivatives):
        print(f"P derivative {n} at p=1 = {value}")
    print("a=b=1: P/sqrt(2) =", equal_axes, "+ O(d^4)", flush=True)
    return moments


def perimeter_formula(p, a, b):
    eta = 1/p - 1
    r = mp.sqrt(a*a + b*b)
    c, s = a*a/(r*r), b*b/(r*r)

    def f(t):
        u = mp.log(t) - mp.log1p(-t)
        z = eta*u
        g_plus = mp.sqrt(c*mp.exp(z) + s*mp.exp(-z))
        g_minus = mp.sqrt(c*mp.exp(-z) + s*mp.exp(z))
        return mp.exp(eta*(mp.log(t)+mp.log1p(-t))/2)*(g_plus+g_minus)

    # Reflection keeps every endpoint evaluation away from t=1.
    return 4*(1+eta)*r*mp.quad(f, [0, mp.mpf(1)/8, mp.mpf(1)/2])


def perimeter_direct(p, a, b):
    """Independent arc length: x=a*cos(theta)^(2/p), y=b*sin(theta)^(2/p)."""
    power = 4/p - 2

    def f(theta):
        sn, cs = mp.sin(theta), mp.cos(theta)
        first = mp.sqrt(a*a*cs**power*sn**2 + b*b*sn**power*cs**2)
        reflected = mp.sqrt(b*b*cs**power*sn**2 + a*a*sn**power*cs**2)
        return first + reflected

    return 8/p*mp.quad(f, [0, mp.pi/8, mp.pi/4])


def closed_derivatives(a, b):
    r = mp.sqrt(a*a + b*b)
    return [
        4*r,
        mp.mpf(0),
        4*mp.pi**2*a*a*b*b/(3*r**3),
        24*a*a*b*b/r**3*(mp.zeta(3)-mp.pi**2/3),
    ]


def error_metrics(value, expected):
    error = abs(value - expected)
    scaled = error/max(mp.mpf(1), abs(expected))
    # ponytail: numerical agreement, capped below working precision; use ball
    # integration if certified enclosures rather than numerical evidence are needed.
    digits = 100 if scaled == 0 else min(100, int(mp.floor(-mp.log10(scaled))))
    return error, digits


def fmt(value, digits=38):
    return mp.nstr(value, digits, min_fixed=0, max_fixed=0)


def numerical_check():
    mp.dps = 110
    print("NUMERICAL: working decimal precision =", mp.dps, flush=True)
    representation = []
    for aa, bb in [(1, 2), (3, 1)]:
        a, b = mp.mpf(aa), mp.mpf(bb)
        for pp in ["0.9", "1.1", "1.5"]:
            p = mp.mpf(pp)
            direct = perimeter_direct(p, a, b)
            transformed = perimeter_formula(p, a, b)
            error, digits = error_metrics(transformed, direct)
            assert digits >= 30, (aa, bb, pp, error, digits)
            row = [f"({aa},{bb})", pp, fmt(direct), fmt(error, 6), str(digits)]
            representation.append(row)
            print("integral", *row, flush=True)

    nodes = list(range(-3, 4))
    weights = [sp.finite_diff_weights(n, nodes, 0)[n][-1] for n in range(4)]
    print("finite difference nodes =", nodes)
    for n, row in enumerate(weights):
        print(f"derivative {n} weights =", row)
    step = mp.mpf("1e-12")
    print("finite difference steps =", fmt(step), fmt(step/2), flush=True)
    derivative_rows = []
    radius_rows = []
    for aa, bb in [(1, 1), (1, 2), (2, 5), (1, 10)]:
        a, b = mp.mpf(aa), mp.mpf(bb)
        closed = closed_derivatives(a, b)
        samples = {}
        approximations = []
        for h in [step, step/2]:
            for k in nodes:
                p = 1+k*h
                if p not in samples:
                    samples[p] = perimeter_direct(p, a, b)
            estimates = []
            for n in range(4):
                value = mp.fsum(
                    mp.mpf(int(w.p))/int(w.q)*samples[1+k*h]
                    for k, w in zip(nodes, weights[n])
                )/h**n
                error, digits = error_metrics(value, closed[n])
                assert digits >= 30, (aa, bb, n, h, error, digits)
                estimates.append((value, error, digits))
            approximations.append(estimates)
        for n in range(4):
            value, error, digits = approximations[1][n]
            row = [
                f"({aa},{bb})", str(n), fmt(closed[n]), fmt(value),
                fmt(approximations[0][n][1], 6), fmt(error, 6), str(digits),
            ]
            derivative_rows.append(row)
            print("derivative", *row, flush=True)
        radius = mp.sqrt(mp.log(b/a)**2 + mp.pi**2/4)
        radius_rows.append([f"({aa},{bb})", fmt(radius)])
        print("Borel radius", radius_rows[-1][0], radius_rows[-1][1], flush=True)
    return representation, derivative_rows, radius_rows


def table(headers, rows):
    return "\n".join([
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"]*len(headers)) + " |",
        *("| " + " | ".join(row) + " |" for row in rows),
    ])


def write_report(moments, representation, derivative_rows, radius_rows):
    moment_table = table(["Moment on the unit interval", "Exact value"], [
        [name.replace("integral ", r"$\int_0^1 ") + r"\,dt$", "$"+sp.latex(value)+"$"]
        for name, value in moments.items()
    ])
    integral_table = table(
        ["Axes $(a,b)$", "$p$", "Direct perimeter", "Absolute residual", "Scaled digits"],
        representation,
    )
    derivative_table = table(
        ["Axes $(a,b)$", "Derivative order", "Closed form evaluated", "Direct finite difference, $h/2$",
         "Absolute residual, $h$", "Absolute residual, $h/2$", "Scaled digits, $h/2$"],
        derivative_rows,
    )
    radius_table = table(["Axes $(a,b)$", "$R$"], radius_rows)
    report = r"""# R2: perimeter derivatives at the rhombus

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

""" + moment_table + r"""

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

""" + integral_table + r"""

Derivatives below use only $P_{\rm direct}$, with the centered **7-point**
stencil at offsets $-3,-2,-1,0,1,2,3$, computed at
$h=10^{-12}$ and $h/2=5\times10^{-13}$.
Orders one and two have truncation error $O(h^6)$; order three has
$O(h^4)$ error. Order zero is the direct integral at the center.
Both step sizes independently passed the **30-digit** assertion for every
entry. The first derivative is checked by absolute residual because its
exact value is zero. The same scale definition covers the other entries.

""" + derivative_table + r"""

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

""" + radius_table + r"""

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
"""
    Path(__file__).with_name("result.md").write_text(report, encoding="utf-8")
    print(report, flush=True)


if __name__ == "__main__":
    moment_values = symbolic_check()
    integral_evidence, derivative_evidence, borel_evidence = numerical_check()
    write_report(moment_values, integral_evidence, derivative_evidence, borel_evidence)
