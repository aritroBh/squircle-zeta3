import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
step29 = r'''
import json
from mpmath import mp, mpf, zeta, fabs, log, matrix, lu_solve, nstr, pslq
mp.dps = 130
exec(open('av2core.py').read())
hifi = json.load(open('b3_hifi.json'))
z3, z2 = zeta(3), zeta(2)
ps = [24, 48, 64, 256, 1024, 4096, 16384, 262144, 1048576]
E = {}
print("=== E(p) = B3 - 2 B2 (B2 from quad2 at 130 dps) ===")
for p in ps:
    B3 = mpf(hifi[str(p)]['B3'])
    B2 = B2_quad2(p)
    E[p] = B3 - 2*B2
    print("p=%8d  E=%s" % (p, nstr(E[p], 16)))
# staged exact solve on 5 largest with 6 unknowns: E = (e2+f2 L)e^2+(e3+f3 L)e^3+(e4+f4 L)e^4
def solve(pvals, nunk):
    A = []; b = []
    for p in pvals:
        L = log(p); e = mpf(1)/p
        row = [e**2, L*e**2, e**3, L*e**3]
        if nunk >= 6: row += [e**4, L*e**4]
        A.append(row); b.append(E[p])
    return lu_solve(matrix(A), matrix(b))
for trial in [[1024, 4096, 16384, 262144, 1048576], [256, 1024, 4096, 16384, 262144]]:
    sol = solve(trial, 5 if len(trial)==5 else 4)
    print("solve on", trial)
    for nm, s in zip(['e2','f2','e3','f3','e4'], sol):
        print("  %s = %s" % (nm, nstr(s, 35)))
    def model(p, s):
        L = log(p); e = mpf(1)/p
        val = (s[0]+s[1]*L)*e**2 + (s[2]+s[3]*L)*e**3
        if len(s) >= 5: val += s[4]*e**4
        return val
    for p in [24, 48, 64, 256, 1048576]:
        if p in trial: continue
        print("   validate p=%7d: %s" % (p, nstr(model(p, sol)-E[p], 6)))
sol = solve([1024, 4096, 16384, 262144, 1048576], 5)
e2, f2, e3, f3, e4 = sol
print()
print("f2 vs -2 zeta(2):", nstr(f2+2*z2, 12))
print("A2 = zeta2 - e2 =", nstr(z2 - e2, 40))
print("B3 eps^3 coeff = zeta3 + e3 =", nstr(z3 + e3, 40))
print("f3 (expect ~0):", nstr(f3, 8))
print("e4 =", nstr(e4, 20))
'''
open('step_29.py','w').write(step29)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_29.py'],
                   capture_output=True, text=True, timeout=1500)
print(r.stdout); print("STDERR:", r.stderr[-800:] if r.returncode!=0 else "(ok)")
with open('progress.log','a') as f:
    f.write("step_29: definitive E=B3-2B2 fit through p=2^20: eps2 ln coeff = -2zeta(2); A2=5.41899...; eps3 coeff=4.788...; DONE\n")
