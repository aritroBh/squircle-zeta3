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
allok = True
for dl in range(2, 11):
    U = sp.together(sum((-1)**s * T5.forced(dl, s) for s in range(0, dl//2+1)))
    a = sp.together(T5.arho(dl))
    nn, dd = sp.fraction(sp.together(U - a))
    C = 2*dl*(dl-1)
    if dl % 2 == 0:
        pre = (-1)**(dl//2+1) * C
    else:
        pre = (-1)**((dl-1)//2) * C * i * (d - (dl+3)//2)
    q, r = sp.div(sp.Poly(nn, d, m), sp.Poly(sp.expand(pre), d, m))
    ok = (sp.expand(r.as_expr()) == 0)
    allok &= ok
    # residual polynomial degree
    deg = sp.degree(sp.expand(nn/pre), d) if ok else -1
    print(f"delta={dl}: nn = {pre} * Q, exact division remainder 0: {ok}; deg Q in d: {deg}", flush=True)
print("ALL exact:", allok)
'''
script = script.replace('WORK', repr(work))
open(work+'/taskE_defect.py','w').write(script)
r = subprocess.run(['nice','-n','10',PY,'taskE_defect.py'], cwd=work, capture_output=True, text=True, timeout=1700)
print(r.stdout)
print("ERR", r.stderr[-300:] if r.stderr else "")
