
# step_07.py : circle end self-dual V(s), q = 2/(1-2s), V = pi + a2 s^2 + a4 s^4 + ... EXACT
import sympy as sp
import av9core as C
import time
t0=time.time()
eta, t, x, y = C.eta, C.t, C.x, C.y
lv, lw = C.lv, C.lw
s, v = sp.symbols('s v')
pg = sp.polygamma
N = 8

def beta_half(j, maxdeg):
    """series of B(x+j+1/2, y+1/2)/B(j+1/2,1/2) to total degree maxdeg"""
    terms = 0
    for k in range(1, maxdeg+1):
        terms += (pg(k-1, j+sp.Rational(1,2))*x**k + pg(k-1, sp.Rational(1,2))*y**k
                  - pg(k-1, j+1)*(x+y)**k)/sp.factorial(k)
    terms = sp.expand(terms.subs({x: t*x, y: t*y}))
    e = C.series_exp_var(terms, t, maxdeg)
    return sum(e[k]*t**k for k in range(maxdeg+1)).subs(t,1)

Bj = {}
def get_Bj(j):
    if j not in Bj:
        Bj[j] = sp.expand(beta_half(j, N))
    return Bj[j]

def int_half(j, a, b):
    """int_0^1 v^j (v(1-v))^{-1/2} lv^a lw^b dv, exact"""
    B = get_Bj(j)
    coeff = sp.expand(B).coeff(x, a).coeff(y, b)
    return sp.expand(sp.pi*sp.binomial(2*j, j)/4**j * sp.factorial(a)*sp.factorial(b) * coeff)

def integrate_vpoly(expr):
    """expr: polynomial in v, lv, lw (times (v(1-v))^{-1/2} understood)."""
    expr = sp.expand(expr)
    total = 0
    # separate v-powers: write as polynomial in v with coeffs in lv,lw
    polyv = sp.Poly(expr, v)
    for (j,), cf in polyv.terms():
        pl = sp.Poly(sp.expand(cf), lv, lw)
        for (a,b), c2 in pl.terms():
            total += c2*int_half(j, a, b)
    return sp.expand(total)

# build coefficient functions in s
S1 = sum((4*s/(1-2*s))**k * (v*lv**k + (1-v)*lw**k)/sp.factorial(k) for k in range(1, N+1))
S1 = sp.series(S1, s, 0, N+1).removeO()  # S - 1
lnS = C.series_log1p_var(S1, s, N)
lnS_ser = sum(lnS[k]*s**k for k in range(1, N+1))
alpha = -sp.Rational(1,2) - s
gamma = sp.Rational(1,2) - s
lnI = sp.log(1-2*s) + alpha*(lv+lw) + gamma*lnS_ser
lnI = sp.series(sp.expand(lnI), s, 0, N+1).removeO()
e = C.series_exp_var(lnI - lnI.coeff(s,0), s, N)
factor = sp.exp(lnI.coeff(s,0))   # = exp(0) = 1? alpha(0)L etc: lnI(0) = 0 + (-1/2)L + (1/2)lnS(0); lnS(0)=0
print("lnI(0) =", sp.expand(lnI.coeff(s,0)))
# NOTE lnI(0) contains -L/2 (the base measure) -> keep it in e_0, do NOT exponentiate into factor
# We instead treat base (v(1-v))^{-1/2} separately: e coefficients multiply it.
a = {}
for n in range(0, N+1):
    en = e[n]  # coefficient of s^n, with base (v(1-v))^{-1/2} stripped
    val = integrate_vpoly(en)
    a[n] = sp.simplify(val)
    print(f"a_{n} done t={time.time()-t0:.1f}  a_{n} =", a[n], flush=True)
for n in range(0, N+1):
    print(f"a_{n}/pi =", sp.simplify(a[n]/sp.pi), "  numeric:", sp.N(a[n], 20))
