import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
step26 = r'''
import json
from mpmath import mp, mpf, zeta, fabs, log, matrix, lu_solve, nstr
mp.dps = 80
tab = json.load(open('b3_hifi.json'))
ps = [24, 48, 64, 256, 1024, 4096, 16384]
z3 = zeta(3)
data = []
for p in ps:
    v = mpf(tab[str(p)]['B3'])
    data.append((p, v - 2*z3))
    print("p=%6d  B3-2zeta3 = %s" % (p, nstr(v-2*z3, 18)))
# model: f = (A2 + B2 ln p) e^2 + (A3 + B3 ln p) e^3 + (A4 + B4 ln p) e^4 + (A5 + B5 ln p) e^5
unk = ['A2','B2','A3','B3','A4','B4','A5','B5']
A = []; b = []
for (p, f) in data:
    e = mpf(1)/p; L = log(p)
    A.append([e**2, L*e**2, e**3, L*e**3, e**4, L*e**4, e**5, L*e**5])
    b.append(f)
Am = matrix(A); bm = matrix(b)
sol = lu_solve(Am.T*Am, Am.T*bm)
res = Am*sol-bm
print("residual max:", nstr(max(fabs(x) for x in res),4))
for name, s in zip(unk, sol):
    print(name, "=", nstr(s, 40))
print()
print("B2 vs -2 zeta(2):", nstr(sol[1] + 2*zeta(2), 10))
print("A2 candidates: 2*zeta(3)+? ; zeta(2)^2/2? =", nstr(zeta(2)**2/2, 15), " ; 3*zeta(2)?=", nstr(3*zeta(2),12))
print("A2 =", nstr(sol[0], 30))
# try PSLQ on A2 and B2
from mpmath import pslq
for nm, val in [('A2', sol[0]), ('B2', sol[1]), ('A3', sol[2]), ('B3', sol[3])]:
    basis = [mpf(1), zeta(2), zeta(3), zeta(2)**2, zeta(2)*zeta(3), zeta(3)**2, zeta(4), zeta(5),
             log(2), zeta(2)*log(2), zeta(3)*log(2), log(2)**2]
    rel = pslq([val]+basis, tol=mpf('1e-20'), maxcoeff=10**8, maxsteps=1000)
    print(nm, "pslq:", rel)
'''
open('step_26.py','w').write(step26)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_26.py'],
                   capture_output=True, text=True, timeout=900)
print(r.stdout); print("STDERR:", r.stderr[-800:] if r.returncode!=0 else "(ok)")
