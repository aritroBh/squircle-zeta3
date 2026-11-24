import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
step27 = r'''
import json
from mpmath import mp, mpf, zeta, fabs, log, matrix, lu_solve, nstr, pslq
mp.dps = 80
tab = json.load(open('b3_hifi.json'))
z3 = zeta(3)
g = {}
for p in [24, 48, 64, 256, 1024, 4096, 16384]:
    v = mpf(tab[str(p)]['B3'])
    g[p] = (2*z3 - v)*p**2    # positive
def solve4(pvals, unknowns):
    # model g = A2 + B2 ln p + (A3+B3 ln p)/p  [4 unk] or + (A4+B4 ln p)/p^2 [6 unk]
    A = []; b = []
    for p in pvals:
        L = log(p); e = mpf(1)/p
        row = [mpf(1), L]
        if len(unknowns) >= 4: row += [e, L*e]
        if len(unknowns) >= 6: row += [e**2, L*e**2]
        A.append(row); b.append(g[p])
    Am = matrix(A); bm = matrix(b)
    return lu_solve(Am, bm)
sol = solve4([256, 1024, 4096, 16384], 4)
A2, B2, A3, B3 = sol
print("4-point exact solve (p=256..16384):")
print("A2 =", nstr(A2, 40))
print("B2 =", nstr(B2, 40))
print("A3 =", nstr(A3, 40))
print("B3 =", nstr(B3, 40))
def model(p, sol):
    L = log(p); e = mpf(1)/p
    s = sol
    val = s[0] + s[1]*L
    if len(s) >= 4: val += (s[2] + s[3]*L)*e
    if len(s) >= 6: val += (s[4] + s[5]*L)*e**2
    return val
for p in [24, 48, 64, 256]:
    print("predict check p=%5d: model-g = %s" % (p, nstr(model(p, sol)-g[p], 8)))
print()
print("B2 vs -?  sign: g has +B2 ln p: B2 vs 2 zeta(2):", nstr(B2 - 2*zeta(2), 12))
print("A2 candidates:")
basis = [mpf(1), zeta(2), zeta(3), zeta(2)**2, zeta(2)*zeta(3), zeta(3)**2, zeta(4), zeta(5), log(2), zeta(2)*log(2), zeta(3)*log(2), log(2)**2, log(2)*zeta(4)]
for nm, val in [('A2', A2), ('B2', B2), ('A3', A3), ('B3', B3)]:
    rel = pslq([val]+basis, tol=mpf('1e-24'), maxcoeff=10**9, maxsteps=1000)
    print(" ", nm, "=", nstr(val, 30), " pslq:", rel)
# 6-unknown solve on 6 largest
sol6 = solve4([64, 256, 1024, 4096, 16384, 24], 6)
print("6-point solve:")
for nm, s in zip(['A2','B2','A3','B3','A4','B4'], sol6):
    print(nm, "=", nstr(s, 30))
'''
open('step_27.py','w').write(step27)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_27.py'],
                   capture_output=True, text=True, timeout=900)
print(r.stdout); print("STDERR:", r.stderr[-800:] if r.returncode!=0 else "(ok)")
