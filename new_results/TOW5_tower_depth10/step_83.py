import subprocess, time, os
work = '<project>/research/08_theory/TOW5_tower_depth10'
PY = '<project>/.venv/bin/python'
script = '''
import sys
sys.path.insert(0, WORK)
import sympy as sp
import tow5lib as T5
d, m = T5.d, T5.m
i = sp.I
for dl in [2, 3, 5]:
    U = sp.together(sum((-1)**s * T5.forced(dl, s) for s in range(0, dl//2+1)))
    a = sp.together(T5.arho(dl))
    nn, dd = sp.fraction(sp.together(U - a))
    print(f"delta={dl}: nn = {sp.expand(nn)}, dd = {sp.factor(dd)}")
    print(f"   nn/4 = {sp.expand(nn/4) if dl==2 else ''}")
'''
script = script.replace('WORK', repr(work))
open(work+'/taskE_dbg.py','w').write(script)
r = subprocess.run(['nice','-n','10',PY,'taskE_dbg.py'], cwd=work, capture_output=True, text=True, timeout=1700)
print(r.stdout)
print("ERR", r.stderr[-300:] if r.stderr else "")
