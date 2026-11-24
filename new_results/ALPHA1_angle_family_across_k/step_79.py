dbg4 = r'''
import sys
sys.path.insert(0, '.')
exec(open('step_30.py').read().split('# numeric evaluation')[0])
import mpmath as mp
mp.mp.dps = 60
from hpl_tools import hpl_numeric_all
allw = sorted({w for (w, mono) in Jc}, key=len)
NT = int(mp.mp.dps*3.4/0.1505) + 40
CONST = {(): mp.mpf(1), ('z2',): mp.zeta(2), ('z3',): mp.zeta(3), ('z4',): mp.zeta(4)}
def Jhpl(c):
    vals = hpl_numeric_all(allw, c, NT)
    return mp.fsum(mp.mpf(v.numerator)/v.denominator*CONST[mono]*vals[w] for (w, mono), v in Jc.items())
lam = mp.log(2)
# squircle check
Js = Jhpl(1/mp.sqrt(2))
print("J(squircle) =", mp.nstr(Js, 40))
print("expected    =", mp.nstr(-(mp.zeta(3)/4 + mp.pi**2*lam/12 + lam**3/6), 40))
# classical candidate: J = z2 log(1-s2) - log(1-s2) Li2(s2) + 2 z3 - 2 Li3(1-s2) + 2 log(1-s2) Li2(1-s2) + log^2(1-s2) log(s2)
def Jclass(c):
    s2 = c*c
    L = mp.log(1-s2)
    return (mp.zeta(2)*L - L*mp.polylog(2, s2) + 2*mp.zeta(3) - 2*mp.polylog(3, 1-s2)
            + 2*L*mp.polylog(2, 1-s2) + L**2*mp.log(s2))
for c in (mp.mpf(1)/2, mp.mpf(3)/5, 1/mp.sqrt(2)):
    print(f"c={mp.nstr(c,8)}: HPL={mp.nstr(Jhpl(c),30)}")
    print(f"          class={mp.nstr(Jclass(c),30)}  diff={mp.nstr(abs(Jhpl(c)-Jclass(c)),3)}")
'''
open('dbg4.py','w').write(dbg4)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'dbg4.py'], capture_output=True, text=True, timeout=1500)
print(r.stdout)
print(r.stderr[-1500:] if r.returncode else "RC=0")
