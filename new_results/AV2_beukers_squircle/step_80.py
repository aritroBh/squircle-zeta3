import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
step42 = r'''
import json
from mpmath import mp, mpf, zeta, fabs, nstr
mp.dps = 100
exec(open('av2core.py').read())
from math import comb, gcd
from functools import reduce
def lcm(a, b): return a*b//gcd(a, b)
def L(n): return reduce(lcm, range(1, n+1), 1)
def apery(n): return sum(comb(n,k)**2 * comb(n+k,k)**2 for k in range(n+1))

results = {}
for p in [1, 2, 4, 8, 32]:
    Mdict = {}
    for j in range(13):
        for k in range(j, 13):
            Mdict[(j,k)] = M_moment(j, k, p)
    B2p, _ = B2_series(p)
    print("=== p = %s, B2(p) = %s ===" % (p, nstr(B2p, 30)))
    rows = []
    for n in range(0, 13):
        Ln = L_legendre(n, p, Mdict)
        an = apery(n)
        r_n = Ln - 2*an*B2p
        # rational test: d_n^3 r_n near integer?
        dn3 = L(n)**3
        test = dn3*r_n
        nearest = round(test)
        err = fabs(test-nearest)
        rate = fabs(Ln)**(1.0/n) if n else mpf(0)
        print("n=%2d L_n=%s |L|^(1/n)=%s" % (n, nstr(Ln, 15), nstr(rate, 8)))
        print("     r_n=L_n-2*Ap_n*B2(p)=%s ; d_n^3 r_n - nearest int = %s" % (nstr(r_n, 15), nstr(err, 4)))
        rows.append(dict(n=n, Ln=mp.nstr(Ln, 50), rn=mp.nstr(r_n, 50), rate=mp.nstr(rate, 15)))
    results[str(p)] = rows
    json.dump(results, open('legendre_finite.json','w'), indent=1)
print("saved")
'''
open('step_42.py','w').write(step42)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_42.py'],
                   capture_output=True, text=True, timeout=1750)
print(r.stdout); print("STDERR:", r.stderr[-500:] if r.returncode!=0 else "(ok)")
with open('progress.log','a') as f:
    f.write("step_42: Legendre L_n(p) n<=12 at p=1,2,4,8,32; ansatz tests + decay rates. DONE\n")
