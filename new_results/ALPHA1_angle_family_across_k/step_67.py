step26 = r'''"""step_26.py -- weight-three slope: LLL kernel, broad basis, 14 angles, 110 dps."""
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
def slope(c): return (1-c**2)/(2*c**2)*Jval(c)

def basis3(c):
    L = mp.log; A3 = lambda x: mp.polylog(3, x); A2 = lambda x: mp.polylog(2, x)
    lc, l1c, lmc, l2 = L(c), L(1+c), L(1-c), L(2)
    return [('zeta3', mp.zeta(3)),
            ('Li3(c)', A3(c)), ('Li3(-c)', A3(-c)), ('Li3(1-c)', A3(1-c)), ('Li3(1/(1+c))', A3(1/(1+c))),
            ('Li3((1+c)/2)', A3((1+c)/2)), ('Li3((1-c)/(1+c))', A3((1-c)/(1+c))), ('Li3(-(1-c)/(1+c))', A3(-(1-c)/(1+c))),
            ('Li3(2c/(1+c))', A3(2*c/(1+c))), ('Li3(1/2)', A3(mp.mpf(1)/2)),
            ('pi^2 l2', mp.pi**2*l2), ('pi^2 lc', mp.pi**2*lc), ('pi^2 l1c', mp.pi**2*l1c), ('pi^2 lmc', mp.pi**2*lmc),
            ('l2^3', l2**3), ('lc^3', lc**3), ('l1c^3', l1c**3), ('lmc^3', lmc**3),
            ('l2^2 lc', l2**2*lc), ('l2 lc^2', l2*lc**2), ('lc^2 l1c', lc**2*l1c), ('lc l1c^2', lc*l1c**2),
            ('l1c^2 lmc', l1c**2*lmc), ('lc^2 lmc', lc**2*lmc), ('l1c lmc^2', l1c*lmc**2), ('l2 lc l1c', l2*lc*l1c),
            ('l2 l1c^2', l2*l1c**2), ('l2^2 l1c', l2**2*l1c), ('l2 lmc^2', l2*lmc**2), ('l2^2 lmc', l2**2*lmc)]

angles = [mp.mpf(1)/2, mp.mpf(3)/5, mp.mpf(1)/3, mp.mpf(4)/5, mp.mpf(2)/7, mp.mpf(5)/6, mp.mpf(1)/5,
          mp.mpf(3)/7, mp.mpf(5)/8, mp.mpf(2)/9, mp.mpf(5)/7, mp.mpf(1)/6, mp.mpf(3)/8, mp.mpf(4)/9]
print("evaluating slope at", len(angles), "angles...")
SL = [slope(c) for c in angles]
B = [[x for _, x in basis3(c)] for c in angles]
names = ['slope'] + [n for n, _ in basis3(angles[0])]
A = [[SL[i]] + B[i] for i in range(len(angles))]
print("kernel:")
found = []
for q, tail in int_kernel(A, Kexp=60)[:6]:
    if any(q):
        print("  q:", [(cf, names[i]) for i, cf in enumerate(q) if cf], " tail:", tail)
        found.append(q)
print("STEP26 DONE")
'''
open('step_26.py','w').write(step26)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'step_26.py'], capture_output=True, text=True, timeout=1700)
print(r.stdout)
print(r.stderr[-1500:] if r.returncode else "RC=0")
