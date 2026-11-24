
# step_07b.py : post-process a_n with expand_func; verify odd vanish; identify even coefficients
import sympy as sp
import av9core as C
import time, pickle
t0=time.time()
lv, lw = C.lv, C.lw
s, v = sp.symbols('s v')
pg = sp.polygamma
N = 8

def pgeval(m, z):
    """polygamma(m, z) for z = half-integer or integer, fully evaluated via recurrence from 1/2 or 1"""
    z = sp.Rational(z)
    if sp.denom(z) == 2:
        val = sp.expand_func(pg(m, sp.Rational(1,2)))
        j = int(z - sp.Rational(1,2))
        for i in range(j):
            val += (-1)**m * sp.factorial(m) / (sp.Rational(1,2)+i)**(m+1)
        return sp.simplify(val)
    else:
        return sp.expand_func(pg(m, z))

def beta_half(j, maxdeg):
    terms = 0
    for k in range(1, maxdeg+1):
        terms += (pgeval(k-1, j+sp.Rational(1,2))*C.x**k + pgeval(k-1, sp.Rational(1,2))*C.y**k
                  - pgeval(k-1, j+1)*(C.x+C.y)**k)/sp.factorial(k)
    terms = sp.expand(terms.subs({C.x: C.t*C.x, C.y: C.t*C.y}))
    e = C.series_exp_var(terms, C.t, maxdeg)
    return sum(e[k]*C.t**k for k in range(maxdeg+1)).subs(C.t,1)

Bj = {}
def get_Bj(j):
    if j not in Bj: Bj[j] = sp.expand(beta_half(j, N))
    return Bj[j]

def int_half(j, a, b):
    B = get_Bj(j)
    coeff = sp.expand(B).coeff(C.x, a).coeff(C.y, b)
    return sp.expand(sp.pi*sp.binomial(2*j, j)/4**j * sp.factorial(a)*sp.factorial(b) * coeff)

def integrate_vpoly(expr):
    expr = sp.expand(expr)
    total = 0
    polyv = sp.Poly(expr, v)
    for (j,), cf in polyv.terms():
        pl = sp.Poly(sp.expand(cf), lv, lw)
        for (a,b), c2 in pl.terms():
            total += c2*int_half(j, a, b)
    return sp.expand(total)

S1 = sum((4*s/(1-2*s))**k * (v*lv**k + (1-v)*lw**k)/sp.factorial(k) for k in range(1, N+1))
S1 = sp.series(S1, s, 0, N+1).removeO()
lnS = C.series_log1p_var(S1, s, N)
lnS_ser = sum(lnS[k]*s**k for k in range(1, N+1))
alpha = -sp.Rational(1,2) - s
gamma = sp.Rational(1,2) - s
lnI = sp.log(1-2*s) + alpha*(lv+lw) + gamma*lnS_ser
lnI = sp.series(sp.expand(lnI), s, 0, N+1).removeO()
e = C.series_exp_var(lnI - lnI.coeff(s,0), s, N)

a = {}
for n in range(0, N+1):
    a[n] = sp.simplify(integrate_vpoly(e[n]))
    print(f"a_{n} =", a[n], f"   (t={time.time()-t0:.0f}s)", flush=True)
    print("      numeric:", sp.N(a[n], 25))
pickle.dump({n: str(a[n]) for n in a}, open('scratchpad/circle_selfdual.pkl','wb'))
