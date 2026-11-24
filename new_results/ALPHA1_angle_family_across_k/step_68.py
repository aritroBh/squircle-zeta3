step27 = r'''"""step_27.py -- identify J(r) = slope*2c^2/(1-c^2) in a pure log-zeta basis (LLL), then verify."""
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 110
from lllkernel import int_kernel

def Jval(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def lmt(w): return mp.log((1-w)*(w+r)/(2*(1-r)*w))
    def f(w):
        return (1/(w+1) - 1/(w+r))*lt(w)**2 + (1/w - 2/(w+1))*lt(w)*lmt(w)
    return mp.quad(f, [r, mp.sqrt(r), 1])

def basisJ(c):
    r = (1-c)/(1+c)
    L = mp.log
    lc, l1c, lmc, l2 = L(c), L(1+c), L(1-c), L(2)
    X = {'c': lc, '1+c': l1c, '1-c': lmc, '2': l2}
    out = [('zeta3', mp.zeta(3))]
    for nx, x in X.items():
        out.append((f'pi^2 l{nx}', mp.pi**2*x))
        out.append((f'l{nx}^3', x**3))
    keys = list(X)
    for i in range(len(keys)):
        for j in range(len(keys)):
            if i != j:
                out.append((f'l{keys[i]}^2 l{keys[j]}', X[keys[i]]**2*X[keys[j]]))
    out.append(('lc l1c l2', lc*l1c*l2))
    out.append(('lc lmc l2', lc*lmc*l2))
    out.append(('l1c lmc l2', l1c*lmc*l2))
    out.append(('lc l1c lmc', lc*l1c*lmc))
    return out

angles = [mp.mpf(1)/2, mp.mpf(3)/5, mp.mpf(1)/3, mp.mpf(4)/5, mp.mpf(2)/7, mp.mpf(5)/6, mp.mpf(1)/5,
          mp.mpf(3)/7, mp.mpf(5)/8, mp.mpf(2)/9, mp.mpf(5)/7, mp.mpf(1)/6, mp.mpf(3)/8, mp.mpf(4)/9]
print("evaluating J...")
Js = [Jval(c) for c in angles]
B = [[x for _, x in basisJ(c)] for c in angles]
names = ['J'] + [n for n, _ in basisJ(angles[0])]
A = [[Js[i]] + B[i] for i in range(len(angles))]
for q, tail in int_kernel(A, Kexp=60)[:6]:
    if any(q) and tail < 10**20:
        print("  q:", [(cf, names[i]) for i, cf in enumerate(q) if cf], " tail:", tail)
print("STEP27 DONE")
'''
open('step_27.py','w').write(step27)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'step_27.py'], capture_output=True, text=True, timeout=1700)
print(r.stdout)
print(r.stderr[-1500:] if r.returncode else "RC=0")
