step29 = r'''"""step_29.py -- slope identification: LLL with M=30 angles > N columns, J(r) (not slope), 90 dps."""
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 90
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
    L = mp.log; A3 = lambda x: mp.polylog(3, x); A2 = lambda x: mp.polylog(2, x)
    lc, l1c, lmc, l2 = L(c), L(1+c), L(1-c), L(2)
    Lc = logc = lc
    out = [('zeta3', mp.zeta(3))]
    for nm, a in [('c', c), ('-c', -c), ('1-c', 1-c), ('1/(1+c)', 1/(1+c)), ('(1+c)/2', (1+c)/2),
                  ('(1-c)/(1+c)', (1-c)/(1+c)), ('-(1-c)/(1+c)', -(1-c)/(1+c)), ('2c/(1+c)', 2*c/(1+c)), ('1/2', mp.mpf(1)/2)]:
        out.append((f'Li3({nm})', A3(a)))
    for nm, x in [('2', l2), ('c', lc), ('1+c', l1c), ('1-c', lmc)]:
        out.append((f'pi^2 l{nm}', mp.pi**2*x))
        out.append((f'l{nm}^3', x**3))
    for n1, x1 in [('2', l2), ('c', lc), ('1+c', l1c), ('1-c', lmc)]:
        for n2, x2 in [('2', l2), ('c', lc), ('1+c', l1c), ('1-c', lmc)]:
            if n1 != n2:
                out.append((f'l{n1}^2 l{n2}', x1**2*x2))
    out.append(('Li2(c) lc', A2(c)*lc)); out.append(('Li2(c) l2', A2(c)*l2))
    out.append(('Li2(-c) lc', A2(-c)*lc)); out.append(('Li2(-c) l2', A2(-c)*l2))
    out.append(('lc l1c l2', lc*l1c*l2)); out.append(('lc lmc l2', lc*lmc*l2))
    return out

angles = [mp.mpf(n)/d for n, d in [(1,2),(3,5),(1,3),(4,5),(2,7),(5,6),(1,5),(3,7),(5,8),(2,9),(5,7),(1,6),
          (3,8),(4,9),(1,7),(2,3),(1,4),(3,10),(5,9),(4,7),(5,11),(6,7),(2,5),(7,9),(1,8),(3,11),(5,12),(7,10),(2,11),(1,9)]]
print("evaluating J at", len(angles), "angles...")
Js = [Jval(c) for c in angles]
B = [[x for _, x in basisJ(c)] for c in angles]
names = ['J'] + [n for n, _ in basisJ(angles[0])]
A = [[Js[i]] + B[i] for i in range(len(angles))]
print("N columns:", len(names), "; M angles:", len(angles))
for q, tail in int_kernel(A, Kexp=50)[:8]:
    if any(q):
        short = all(abs(cf) < 10**12 for cf in q)
        print("  q:", [(cf, names[i]) for i, cf in enumerate(q) if cf][:12], "..." if len([cf for cf in q if cf])>12 else "", " tail:", tail, " SHORT" if short else "")
print("STEP29 DONE")
'''
open('step_29.py','w').write(step29)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'step_29.py'], capture_output=True, text=True, timeout=1700)
print(r.stdout)
print(r.stderr[-1500:] if r.returncode else "RC=0")
