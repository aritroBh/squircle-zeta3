
step09 = r'''
# step_09.py : bivariate circle end V(K_{p1},K_{r2}), p1=2/(1+2 s1), r2=2/(1+2 s2), total order <=5
# ln I = ln(1+2 s2) + (-1/2-s2)L + (1/2-s1) ln S ; S = v e^{w lv}+(1-v) e^{w lw}, w = 2(s1-s2)/(1-2 s1)
import sympy as sp
import av9core as C
import time, pickle
t0=time.time()
lv, lw, t = C.lv, C.lw, C.t
s1, s2, v = sp.symbols('s1 s2 v')
pg = sp.polygamma
N = 5

def pgeval(m, z):
    z = sp.Rational(z)
    if sp.denom(z) == 2:
        val = sp.expand_func(pg(m, sp.Rational(1,2)))
        for i in range(int(z - sp.Rational(1,2))):
            val += (-1)**m * sp.factorial(m) / (sp.Rational(1,2)+i)**(m+1)
        return sp.simplify(val)
    return sp.expand_func(pg(m, z))

Bj = {}
def get_Bj(j):
    if j not in Bj:
        terms = 0
        for k in range(1, N+1):
            terms += (pgeval(k-1, j+sp.Rational(1,2))*C.x**k + pgeval(k-1, sp.Rational(1,2))*C.y**k
                      - pgeval(k-1, j+1)*(C.x+C.y)**k)/sp.factorial(k)
        terms = sp.expand(terms.subs({C.x: C.t*C.x, C.y: C.t*C.y}))
        e = C.series_exp_var(terms, C.t, N)
        Bj[j] = sp.expand(sum(e[k]*C.t**k for k in range(N+1)).subs(C.t,1))
    return Bj[j]

def int_half(j, a, b):
    B = get_Bj(j)
    coeff = sp.expand(B).coeff(C.x, a).coeff(C.y, b)
    return sp.expand(sp.pi*sp.binomial(2*j, j)/4**j * sp.factorial(a)*sp.factorial(b) * coeff)

def integrate_vpoly(expr):
    expr = sp.expand(expr)
    total = 0
    for (j,), cf in sp.Poly(expr, v).terms():
        for (a,b), c2 in sp.Poly(sp.expand(cf), lv, lw).terms():
            total += c2*int_half(j, a, b)
    return sp.expand(total)

# series in t (scaling s1,s2)
w = 2*(t*s1-t*s2)*sp.series(1/(1-2*t*s1), t, 0, N+1).removeO()
w = sp.expand(w)
S1 = 0
for k in range(1, N+1):
    S1 += (v*lv**k + (1-v)*lw**k)/sp.factorial(k) * w**k
    S1 = sp.expand(S1)
    S1 = sum(S1.coeff(t,m)*t**m for m in range(1, N+1))
lnS = C.series_log1p_var(S1, t, N)
bracket = sp.log(1+2*t*s2) - t*s2*(lv+lw) + (sp.Rational(1,2)-t*s1)*sum(lnS[k]*t**k for k in range(1,N+1))
bracket = sp.series(sp.expand(bracket), t, 0, N+1).removeO()
e = C.series_exp_var(bracket - bracket.coeff(t,0), t, N)
b = {}
for n in range(1, N+1):
    en = sp.expand(e[n])
    for (j1,j2), cf in sp.Poly(en, s1, s2).terms():
        if j1+j2 != n: continue
        b[(j1,j2)] = sp.simplify(integrate_vpoly(cf))
    print(f"order {n} done t={time.time()-t0:.0f}s", flush=True)
for k in sorted(b):
    print(f"b_{k} =", b[k])
    print("    numeric:", sp.N(b[k], 20))
pickle.dump({str(k): str(b[k]) for k in b}, open('scratchpad/circle_bivariate.pkl','wb'))
'''
open('step_09.py','w').write(step09)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','step_09.py'],capture_output=True,text=True,timeout=1700)
print(r.stdout); print(r.stderr[-800:])
