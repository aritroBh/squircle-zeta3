dbg3 = r'''
import sys
sys.path.insert(0, '.')
exec(open('step_30.py').read().split('# numeric evaluation')[0])
import mpmath as mp
mp.mp.dps = 50
print("Jc terms:")
for (w, mono), v in sorted(Jc.items()):
    print("   ", w, mono, v)
print("trailing-zero violation:", [w for (w, mono) in Jc if w and w[-1] == 0])
from hpl_tools import hpl_numeric_all
allw = sorted({w for (w, mono) in Jc}, key=len)
NT = int(mp.mp.dps*3.4/0.1505) + 40
CONST = {(): mp.mpf(1), ('z2',): mp.zeta(2), ('z3',): mp.zeta(3), ('z4',): mp.zeta(4)}
for c in (mp.mpf(1)/2, mp.mpf(3)/5):
    vals = hpl_numeric_all(allw, c, NT)
    tot = mp.fsum(mp.mpf(v.numerator)/v.denominator*CONST[mono]*vals[w] for (w, mono), v in Jc.items())
    print(f"c={mp.nstr(c,6)}: Jc(HPL)={mp.nstr(tot,30)}")
'''
open('dbg3.py','w').write(dbg3)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'dbg3.py'], capture_output=True, text=True, timeout=1500)
print(r.stdout)
print(r.stderr[-1500:] if r.returncode else "RC=0")
