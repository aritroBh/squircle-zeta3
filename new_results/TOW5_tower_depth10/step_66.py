import subprocess, time, os
work = '<project>/research/08_theory/TOW5_tower_depth10'
PY = '<project>/.venv/bin/python'
code = '''
import sys
sys.path.insert(0, %r)
import sympy as sp
import tow5lib as T5
d, m = T5.d, T5.m
for delta in [9, 10]:
    a = T5.arho(delta)
    nn, dd = sp.fraction(sp.together(a))
    print(f"a_{delta} = N_{delta}/({sp.factor(dd)})")
    print(f"  N_{delta} = {sp.expand(nn)}")
print()
for dl in [9, 10]:
    for s in range(1, dl//2+1):
        F = T5.forced(dl, s)
        nn, dd = sp.fraction(sp.together(F))
        print(f"F({dl},{s}) = ({sp.expand(nn)}) / ({sp.factor(dd)})")
''' % (work,)
r = subprocess.run(['nice','-n','10',PY,'-c',code], capture_output=True, text=True, timeout=1700)
print(r.stdout)
print("ERR", r.stderr[-500:] if r.stderr else "")
