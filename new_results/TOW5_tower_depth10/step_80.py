import subprocess, time, os
work = '<project>/research/08_theory/TOW5_tower_depth10'
PY = '<project>/.venv/bin/python'
code = '''
import sys, time
sys.path.insert(0, %r)
import sympy as sp
import tow5lib as T5
d, m = T5.d, T5.m
for dl in range(2, 11):
    U = sp.together(sum((-1)**s * T5.forced(dl, s) for s in range(0, dl//2+1)))
    a = sp.together(T5.arho(dl))
    defect = sp.together(U - a)
    nn, dd = sp.fraction(defect)
    nne = sp.expand(nn)
    # factor the numerator
    fac = sp.factor(nne)
    print(f"delta={dl}: U-a numerator = {fac}")
''' % (work,)
r = subprocess.run(['nice','-n','10',PY,'-c',code], capture_output=True, text=True, timeout=1700)
print(r.stdout)
print("ERR", r.stderr[-500:] if r.stderr else "")
