import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
dbg = r'''
import sys
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
mp.mp.dps = 1000
vals = [mp.zeta(5), mp.pi**2*mp.log(2)**3, mp.log(mp.mpf('1.5'))**5, mp.pi**4*mp.log(2)]
a_true = [1048576, -129600, 7776, -349920]
combo = sum(mp.mpf(ai)*v for ai,v in zip(a_true,vals))
vl = vals + [combo]
rows, dt = W.lll_rows(vl, 950)
import math
def norm(r): return math.sqrt(sum(x*x for x in r[:-1]) + (r[-1]/1e300)**2)
for r in sorted(rows, key=lambda r: sum(x*x for x in r[:-1]))[:4]:
    a=r[:-1]
    print("row a:", a, " last:", r[-1])
    print("  matches planted (up to sign):", a==a_true+[-1] or a==[-x for x in a_true]+[1])
'''
open('dbg7.py','w').write(dbg)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','dbg7.py'], capture_output=True, text=True, cwd='.')
print(r.stdout); print("ERR:", r.stderr[-800:])
