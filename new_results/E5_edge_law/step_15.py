code = r'''# step_02.py: FULL N=6 edge d-series with lumped non-local symbols tau1..5, gam1..5 (as E1 s01).
# Documents E5(t) structure; verifies the cross coefficients g_ab (a+b=5) are tau/gam-free,
# hence equal the reduced ones from step_01. Also extracts g50, g05 (pure) for completeness.
import sympy as sp, pickle, time

t, al, be, L, h, s0 = sp.symbols('t alpha beta L h s0', real=True)
l0, l1 = sp.symbols('l0 l1')
tau = sp.symbols('tau1:6'); gam = sp.symbols('gam1:6')
N = 6

def mul(a, b):
    return [sp.expand(sum(a[i]*b[k-i] for i in range(k+1))) for k in range(N)]
def inv(a):
    r = [sp.Integer(1)] + [0]*(N-1)
    for k in range(1, N):
        r[k] = sp.expand(-sum(a[i]*r[k-i] for i in range(1, k+1)))
    return r
def log1(a):
    x = [0] + a[1:]
    res = [0]*N; pw = [sp.Integer(1)] + [0]*(N-1)
    for j in range(1, N):
        pw = mul(pw, x)
        res = [sp.expand(res[k] + sp.Rational((-1)**(j+1), j)*pw[k]) for k in range(N)]
    return res
def exp0(a):
    res = [sp.Integer(1)] + [0]*(N-1); pw = [sp.Integer(1)] + [0]*(N-1)
    for j in range(1, N):
        pw = mul(pw, a)
        res = [sp.expand(res[k] + pw[j-1]/sp.factorial(j)) if False else sp.expand(res[k] + pw[k]/sp.factorial(j)) for k in range(N)]
    return res

t0 = time.time()
sig = s0 + L*t
T = [sp.Integer(1)] + [al*t*l0**k + be*(1-t)*l1**k + tau[k-1] for k in range(1, N)]
M = [sp.Integer(0)] + [al*l0**k - be*l1**k + gam[k-1] for k in range(1, N)]
F = [sp.expand(T[k]/sp.factorial(k)) for k in range(N)]
Num = [sp.expand(M[k]/sp.factorial(k)) for k in range(N)]
Phi = mul(Num, inv(F))
oneoverp = [sp.Integer((-1)**k) for k in range(N)]
Fpow = exp0([sp.expand(-c) for c in mul(oneoverp, log1(F))])
Phi2 = mul(Phi, Phi)
X = [sp.expand(-2*sig*Phi[k]/L + (sig**2 + h**2)*Phi2[k]/L**2) for k in range(N)]
sq = [sp.Integer(1)] + [0]*(N-1); pw = [sp.Integer(1)] + [0]*(N-1)
for j in range(1, N):
    pw = mul(pw, X)
    sq = [sp.expand(sq[k] + sp.binomial(sp.Rational(1, 2), j)*pw[k]) for k in range(N)]
S = [sp.expand(L*c) for c in sq]
E = mul(Fpow, S)
E = [sp.expand(c) for c in E]
print("build time: %.1f s" % (time.time()-t0))
E5 = E[5]
print("E5 terms:", len(sp.Add.make_args(E5)))
P5 = sp.Poly(E5, l0, l1, *tau, *gam)
print("total degree in all log symbols:", P5.total_degree())
ok = True
for a in range(1, 5):
    b = 5 - a
    gab = sp.expand(E5.coeff(l0, a).coeff(l1, b))
    dep = gab.free_symbols & set(list(tau) + list(gam))
    print(f"g{a}{b}: tau/gam dependence:", dep if dep else "none")
    ok = ok and not dep
# compare with reduced pickle
R = pickle.load(open('steps/step_01_reduced.pkl', 'rb'))
sy = {'t': t, 'alpha': al, 'beta': be, 'L': L, 'h': h, 's0': s0}
for a in range(1, 5):
    b = 5 - a
    red = sp.sympify(R['g'][f'({a}, {b})'], locals=sy)
    full = sp.expand(E5.coeff(l0, a).coeff(l1, b))
    print(f"g{a}{b} reduced==full:", sp.simplify(red - full) == 0)
# pure coefficients g50, g05 (for the record)
g50 = sp.expand(E5.coeff(l0, 5).coeff(l1, 0)); g05 = sp.expand(E5.coeff(l0, 0).coeff(l1, 5))
print("g50 tau/gam dep:", g50.free_symbols & set(list(tau)+list(gam)), " terms:", len(sp.Add.make_args(g50)))
# count terms of E5 by log-degree
from collections import Counter
deg = Counter()
for mon, c in P5.terms():
    deg[sum(mon)] += 1
print("E5 term count by total log-degree:", dict(sorted(deg.items())))
# how many distinct tau/gam monomials appear
tg = Counter()
for mon, c in P5.terms():
    key = tuple(m for m, v in zip(list(tau)+list(gam), mon[2:]) if v)
    tg[key] += 1
print("number of distinct (tau,gam) exponent patterns in E5:", len(tg))
pickle.dump(dict(E5=str(E5)), open('steps/step_02_full.pkl', 'wb'))
with open('steps/step_02_out.txt', 'w') as f:
    f.write(f"E5 total terms {len(sp.Add.make_args(E5))}; by log-degree {dict(sorted(deg.items()))}\n")
    f.write(f"g50 = {g50}\ng05 = {g05}\n")
print("saved steps/step_02_full.pkl / step_02_out.txt")
'''
open('steps/step_02.py','w').write(code)
import subprocess, time
t0=time.time()
r=subprocess.run(['<project>/.venv/bin/python','steps/step_02.py'],capture_output=True,text=True,timeout=3000)
print("wall: %.1f s"%(time.time()-t0))
print(r.stdout); print("ERR:", r.stderr[-2000:])
