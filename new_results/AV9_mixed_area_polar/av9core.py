
# av9core.py : exact series machinery for AV9 mixed-area expansions
import sympy as sp

eta, t, x, y = sp.symbols('eta t x y')
lv, lw = sp.symbols('lv lw')
L2 = sp.log(2)
zeta = sp.zeta

def series_exp_var(poly_in_var, var, deg):
    e = [sp.Integer(1)]
    pc = [sp.expand(poly_in_var).coeff(var, k) for k in range(1, deg+1)]
    for n in range(1, deg+1):
        en = sp.expand(sum(k*pc[k-1]*e[n-k] for k in range(1, n+1)))/n
        e.append(sp.expand(en))
    return e  # list of coefficients e[0..deg]

def series_log1p_var(poly_in_var, var, deg):
    # returns coefficients of log(1+u), u = poly_in_var with zero constant term
    res = [sp.Integer(0)]*(deg+1)
    um = poly_in_var  # u^m
    for m in range(1, deg+1):
        c = sp.Rational((-1)**(m+1), m)
        um_exp = sp.expand(um)
        for k in range(m, deg+1):
            res[k] = sp.expand(res[k] + c*um_exp.coeff(var, k))
        um = sp.expand(um*poly_in_var)
        # truncate um
        um = sum(um.coeff(var, k)*var**k for k in range(deg+1))
    return res

def beta_int_table(maxdeg):
    """int_0^1 lv^a lw^b dv = a! b! [x^a y^b] B(x+1,y+1), exact in zeta(k)."""
    logkern = sum((-1)**k * zeta(k)/k * (x**k + y**k - (x+y)**k) for k in range(2, maxdeg+1))
    logkern_t = sp.expand(logkern.subs({x: t*x, y: t*y}))
    ec = series_exp_var(logkern_t, t, maxdeg)
    exppart = sum(ec[k]*t**k for k in range(maxdeg+1))
    geom = sum((-1)**m * t**m*(x+y)**m for m in range(maxdeg+1))
    Bfull = sp.expand(exppart*geom)
    Bfull = sum(Bfull.coeff(t, k)*t**k for k in range(maxdeg+1)).subs(t, 1)
    Bpoly = sp.Poly(sp.expand(Bfull), x, y)
    tab = {}
    for (a,b), coeff in Bpoly.terms():
        tab[(a,b)] = coeff*sp.factorial(a)*sp.factorial(b)
    tab[(0,0)] = sp.Integer(1)
    return tab

def G_coeffs(N):
    """G(beta) = log(e^{beta lv}+e^{beta lw}) = ln2 + sum g_j beta^j."""
    u = sum((lv**k+lw**k)*sp.Rational(1,2)/sp.factorial(k) * eta**k for k in range(1, N+1))
    lc = series_log1p_var(u, eta, N)
    g = {0: L2}
    for j in range(1, N+1):
        g[j] = lc[j]
    return g

def integrate_poly(expr, tab):
    poly = sp.Poly(sp.expand(expr), lv, lw)
    val = 0
    for (a,b), coeff in poly.terms():
        val += coeff*tab.get((a,b), 0)
    return sp.expand(val)

def selfdual_coeffs(N):
    """c_n of V(eta) (self-dual), n=1..N, exact. Also returns e_n coefficient functions."""
    g = G_coeffs(N)
    A = lv + lw
    h = {}
    for n in range(1, N+1):
        h[n] = sp.expand(sum((-1)**(n-j)*g[j] for j in range(0, n+1)) + (-1)**n*A + sp.Rational((-1)**n, n))
    Hser = sum(h[n]*eta**n for n in range(1, N+1))
    e = series_exp_var(Hser, eta, N)
    return h, e

def general_coeffs(N, alpha, beta_, gamma, P):
    """coefficient functions e_n (n=0..N) of I(v,eta) = P(eta)(v(1-v))^{alpha(eta)}(v^{beta}+(1-v)^{beta})^{gamma}
       alpha, beta_, gamma, P: sympy series polys in eta."""
    g = G_coeffs(N+2)
    A = lv + lw
    # ln I = ln P + alpha*A + gamma*(ln2 + sum g_j beta^j)
    al = [alpha.coeff(eta,k) for k in range(N+1)]
    be = [beta_.coeff(eta,k) for k in range(N+1)]
    ga = [gamma.coeff(eta,k) for k in range(N+1)]
    lnP = sp.log(P)
    lnPser = sp.series(lnP, eta, 0, N+1).removeO()
    lnI = [sp.expand(lnPser.coeff(eta,k)) for k in range(N+1)]
    for n in range(1, N+1):  # alpha*A contribution
        lnI[n] = lnI[n] + al[n]*A
    # gamma * G(beta)
    # beta^j series
    bpow = {0: [sp.Integer(1)]+[sp.Integer(0)]*N}
    prev = [sp.Integer(1)]+[sp.Integer(0)]*N
    for j in range(1, N+1):
        cur = [0]*(N+1)
        for a in range(N+1):
            if prev[a]==0: continue
            for b in range(1, N+1-a):
                cur[a+b] = sp.expand(cur[a+b] + prev[a]*be[b])
        bpow[j] = cur
        prev = cur
    Gser = [sp.Integer(0)]*(N+1)
    Gser[0] = L2
    for n in range(1, N+1):
        s = 0
        for j in range(1, n+1):
            s = s + g[j]*bpow[j][n]
        Gser[n] = sp.expand(s)
    # convolve gamma * Gser
    for n in range(0, N+1):
        s = 0
        for a in range(0, n+1):
            s = s + ga[a]*Gser[n-a]
        lnI[n] = sp.expand(lnI[n] + s)
    factor = sp.exp(lnI[0])
    Hser = sum(lnI[n]*eta**n for n in range(1, N+1))
    e = series_exp_var(Hser, eta, N)
    e = [sp.expand(factor*ei) for ei in e]
    return lnI, e
