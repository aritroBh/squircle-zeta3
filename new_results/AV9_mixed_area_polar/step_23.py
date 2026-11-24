
step09b = r'''
# step_09b.py : bivariate circle end, optimized
import sympy as sp
import av9core as C
import time, pickle
T=time.time()
def log(*a): print(f"[{time.time()-T:6.1f}]", *a, flush=True)
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
        return sp.expand(val)
    return sp.expand_func(pg(m, z))

Bjd = {}
def get_Bjd(j):
    if j not in Bjd:
        terms = 0
        for k in range(1, N+1):
            terms += (pgeval(k-1, j+sp.Rational(1,2))*C.x**k + pgeval(k-1, sp.Rational(1,2))*C.y**k
                      - pgeval(k-1, j+1)*(C.x+C.y)**k)/sp.factorial(k)
        terms = sp.expand(terms.subs({C.x: C.t*C.x, C.y: C.t*C.y}))
        e = C.series_exp_var(terms, C.t, N)
        B = sp.expand(sum(e[k]*C.t**k for k in range(N+1)).subs(C.t,1))
        Bjd[j] = {mn: cf for mn, cf in sp.Poly(B, C.x, C.y).terms()}
        log("Bj", j, "built")
    return Bjd[j]

def int_half(j, a, b):
    coeff = get_Bjd(j).get((a,b), 0)
    return sp.pi*sp.binomial(2*j, j)/4**j * sp.factorial(a)*sp.factorial(b) * coeff

def integrate_vpoly(expr):
    expr = sp.expand(expr)
    total = 0
    for (j,), cf in sp.Poly(expr, v).terms():
        for (a,b), c2 in sp.Poly(sp.expand(cf), lv, lw).terms():
            if c2 == 0: continue
            total += c2*int_half(j, a, b)
    return sp.expand(total)

# series part (from dbg2, known ~50s)
w = sp.expand(2*(t*s1-t*s2)*sum((2*t*s1)**m for m in range(N+1)))
S1 = 0
for k in range(1, N+1):
    S1 += (v*lv**k + (1-v)*lw**k)/sp.factorial(k) * w**k
    S1 = sp.expand(S1)
    S1 = sum(S1.coeff(t,m)*t**m for m in range(1, N+1))
lnS = C.series_log1p_var(S1, t, N); log("lnS")
logp = sum(sp.Rational((-1)**(m+1),m)*(2*t*s2)**m for m in range(1,N+1))
bracket = logp - t*s2*(lv+lw) + (sp.Rational(1,2)-t*s1)*sum(lnS[k]*t**k for k in range(1,N+1))
bracket = sp.expand(bracket)
bracket = sum(bracket.coeff(t,m)*t**m for m in range(1,N+1))
e = C.series_exp_var(bracket, t, N); log("exp")

b = {}
for n in range(1, N+1):
    valn = integrate_vpoly(e[n])
    for (j1,j2), cf in sp.Poly(valn, s1, s2).terms():
        b[(j1,j2)] = sp.expand(cf)
    log(f"order {n} integrated")
for k in sorted(b):
    print(f"b_{k} =", b[k])
    print("    numeric:", sp.N(b[k], 20))
pickle.dump({str(k): str(b[k]) for k in b}, open('scratchpad/circle_bivariate.pkl','wb'))
'''
open('step_09b.py','w').write(step09b)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','step_09b.py'],capture_output=True,text=True,timeout=1500)
print(r.stdout); print(r.stderr[-800:])
