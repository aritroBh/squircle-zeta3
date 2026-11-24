step24 = r'''"""step_24.py -- identify J(r) = slope * 2c^2/(1-c^2), weight-3 basis, many angles + holdouts."""
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 80

def Jval(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def lmt(w): return mp.log((1-w)*(w+r)/(2*(1-r)*w))
    def f(w):
        return (1/(w+1) - 1/(w+r))*lt(w)**2 + (1/w - 2/(w+1))*lt(w)*lmt(w)
    return mp.quad(f, [r, mp.sqrt(r), 1])

angles = [mp.mpf(1)/2, mp.mpf(3)/5, mp.mpf(1)/3, mp.mpf(4)/5, mp.mpf(2)/7, mp.mpf(5)/6,
          mp.mpf(1)/5, mp.mpf(3)/7, mp.mpf(5)/8, mp.mpf(2)/9, mp.mpf(5)/7, mp.mpf(7)/9, mp.mpf(1)/6, mp.mpf(4)/9]
holdout = [mp.mpf(1)/7, mp.mpf(2)/3, 1/mp.sqrt(2), mp.mpf(3)/8]
print("evaluating J...")
Js = [Jval(c) for c in angles]
Jh = [Jval(c) for c in holdout]

def basis(c):
    r = (1-c)/(1+c)
    L = mp.log; A3 = lambda x: mp.polylog(3, x); A2 = lambda x: mp.polylog(2, x)
    lr, l1r, lmr, l2 = L(r), L(1+r), L(1-r), L(2)
    return [('zeta3', mp.zeta(3)),
            ('pi^2 lr', mp.pi**2*lr), ('pi^2 l1r', mp.pi**2*l1r), ('pi^2 lmr', mp.pi**2*lmr), ('pi^2 l2', mp.pi**2*l2),
            ('lr^3', lr**3), ('l1r^3', l1r**3), ('lmr^3', lmr**3), ('l2^3', l2**3),
            ('lr^2 l1r', lr**2*l1r), ('lr l1r^2', lr*l1r**2), ('lr^2 lmr', lr**2*lmr), ('l1r^2 lmr', l1r**2*lmr),
            ('lr l2^2', lr*l2**2), ('lr^2 l2', lr**2*l2), ('l1r l2^2', l1r*l2**2), ('lmr l2^2', lmr*l2**2),
            ('Li3(r)', A3(r)), ('Li3(-r)', A3(-r)), ('Li3(c)', A3(c)), ('Li3(-c)', A3(-c)),
            ('Li3((1-c)/2)', A3((1-c)/2)), ('Li3((1+c)/2)', A3((1+c)/2)), ('Li3(1-r)', A3(1-r)), ('Li3(1-c)', A3(1-c)),
            ('Li2(r) lr', A2(r)*lr), ('Li2(r) l1r', A2(r)*l1r), ('Li2(r) lmr', A2(r)*lmr), ('Li2(r) l2', A2(r)*l2),
            ('Li2(-r) lr', A2(-r)*lr), ('Li2(c) l2', A2(c)*l2), ('Li2(-c) l2', A2(-c)*l2), ('Li2(c) lc', A2(c)*L(c))]

C = [[x for _, x in basis(c)] for c in angles]
names = [n for n, _ in basis(angles[0])]
M, N = len(angles), len(names)
def n2(v): return mp.sqrt(sum(x*x for x in v))
A = [mp.matrix([C[i][j] for i in range(M)]) for j in range(N)]
for j in range(N): A[j] = A[j]/n2(A[j])
sel = []
for j in range(N):
    v = A[j].copy()
    for s in sel: v = v - (v.T*A[s])[0]*A[s]
    if n2(v) > mp.mpf('1e-55'):
        for s in sel: v = v - (v.T*A[s])[0]*A[s]
        v = v/n2(v); A[j] = v; sel.append(j)
print("indep subset:", len(sel), "of", N)
Asel = mp.matrix(M, len(sel))
for k, j in enumerate(sel):
    for i in range(M): Asel[i, k] = C[i][j]
Tv = mp.matrix(Js)
q = mp.lu_solve(Asel.T*Asel, Asel.T*Tv)
resid = Asel*q - Tv
print("LS residual:", mp.nstr(max(abs(x) for x in resid), 3))
coef = {names[j]: q[k] for k, j in enumerate(sel) if abs(q[k]) > mp.mpf('1e-55')}
for n in coef: print(f"   {n}: {mp.nstr(coef[n], 35)}")
for c, Jtrue in zip(holdout, Jh):
    bmap = {n: x for n, x in basis(c)}
    pred = sum(coef[n]*bmap[n] for n in coef)
    print(f"   holdout c={mp.nstr(c,9)}: diff={mp.nstr(abs(Jtrue-pred),3)}")
print("STEP24 DONE")
'''
open('step_24.py','w').write(step24)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'step_24.py'], capture_output=True, text=True, timeout=1700)
print(r.stdout)
print(r.stderr[-1500:] if r.returncode else "RC=0")
