import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
step30 = r'''
import json
from mpmath import mp, mpf, zeta, fabs, log, matrix, lu_solve, nstr
mp.dps = 130
exec(open('av2core.py').read())
hifi = json.load(open('b3_hifi.json'))
z3, z2 = zeta(3), zeta(2)
ps = [24, 48, 64, 256, 1024, 4096, 16384, 262144, 1048576]
E = {}
for p in ps:
    B3 = mpf(hifi[str(p)]['B3'])
    E[p] = B3 - 2*B2_quad2(p)
def solve(pvals, nunk):
    A = []; b = []
    for p in pvals:
        L = log(p); e = mpf(1)/p
        row = [e**2, L*e**2, e**3, L*e**3]
        if nunk >= 6: row += [e**4, L*e**4]
        if nunk >= 8: row += [e**5, L*e**5]
        A.append(row); b.append(E[p])
    return lu_solve(matrix(A), matrix(b))
def model(p, s):
    L = log(p); e = mpf(1)/p
    val = (s[0]+s[1]*L)*e**2 + (s[2]+s[3]*L)*e**3
    if len(s) >= 6: val += (s[4]+s[5]*L)*e**4
    if len(s) >= 8: val += (s[6]+s[7]*L)*e**5
    return val
sol6 = solve([256, 1024, 4096, 16384, 262144, 1048576], 6)
print("6-unknown solve (256..1048576):")
for nm, s in zip(['e2','f2','e3','f3','e4','f4'], sol6):
    print("  %s = %s" % (nm, nstr(s, 35)))
for p in [24, 48, 64]:
    print("   validate p=%4d: %s" % (p, nstr(model(p, sol6)-E[p], 6)))
sol8 = solve([24, 48, 64, 256, 1024, 4096, 16384, 262144], 8)
print("8-unknown solve (24..262144):")
for nm, s in zip(['e2','f2','e3','f3','e4','f4','e5','f5'], sol8):
    print("  %s = %s" % (nm, nstr(s, 35)))
for p in [1048576]:
    print("   validate p=%8d: %s" % (p, nstr(model(p, sol8)-E[p], 6)))
e2, f2, e3, f3 = sol6[0], sol6[1], sol6[2], sol6[3]
print()
print("[6-unk] f2 vs -2zeta2:", nstr(f2+2*z2, 10))
print("[6-unk] A2 = zeta2-e2 =", nstr(z2-e2, 35))
print("[6-unk] B3 eps3 = zeta3+e3 =", nstr(z3+e3, 35))
e28, f28, e38, f38 = sol8[0], sol8[1], sol8[2], sol8[3]
print("[8-unk] f2 vs -2zeta2:", nstr(f28+2*z2, 10))
print("[8-unk] A2 =", nstr(z2-e28, 35))
print("[8-unk] B3 eps3 =", nstr(z3+e38, 35))
# candidate checks for A2
print()
cands = {
 "2z2+2z3-z2/6": 2*z2+2*z3-z2/6,
 "9z3/2": 9*z3/2,
 "2z2^2": 2*z2**2,
 "2z2+2z2ln2": 2*z2+2*z2*log(2),
 "z2+3z3+? ": z2+3*z3,
}
A2 = z2-e2
print("A2 =", nstr(A2, 30))
for nm, v in cands.items():
    print("  %s = %s  diff=%s" % (nm, nstr(v,20), nstr(v-A2, 6)))
'''
open('step_30.py','w').write(step30)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_30.py'],
                   capture_output=True, text=True, timeout=1500)
print(r.stdout); print("STDERR:", r.stderr[-800:] if r.returncode!=0 else "(ok)")
