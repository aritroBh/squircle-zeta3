import subprocess, time, os
work = '<project>/research/08_theory/TOW5_tower_depth10'
PY = '<project>/.venv/bin/python'
code = '''
import sys, time
sys.path.insert(0, %r)
import sympy as sp
import tow5lib as T5
d, m = T5.d, T5.m
i = sp.I
rows = []
for dl in range(2, 11):
    U = sp.together(sum((-1)**s * T5.forced(dl, s) for s in range(0, dl//2+1)))
    a = sp.together(T5.arho(dl))
    nn, dd = sp.fraction(sp.together(U - a))
    # hypothesized prefactor
    C = 2*dl*(dl-1)
    if dl % 2 == 0:
        pre = (-1)**(dl//2+1) * C
        Q = sp.together(nn / pre)
    else:
        pre = (-1)**((dl-1)//2) * C * i * (d - (dl+3)//2)
        Q = sp.together(nn / pre)
    qn, qd = sp.fraction(Q)
    ispoly = (sp.expand(qd) == sp.expand(dd)) or (qd == 1)
    # actually Q should be nn/prefactor; check qd divides dd (defect denominator) -- Q is polynomial iff qd==1 after removing dd... 
    # Simpler: check nn/pre is a polynomial
    Ptest = sp.expand(nn / pre)
    ispoly = (sp.denom(sp.together(nn/pre)) == 1)
    rows.append((dl, sp.simplify(nn/pre) is not None, ispoly))
    print(f"delta={dl}: prefactor {sp.nsimplify(C)} * {'(-1)^(dl/2+1)' if dl%2==0 else '(-1)^((dl-1)/2) i (d-(dl+3)/2)'}; defect/prefactor is polynomial: {ispoly}", flush=True)
''' % (work,)
r = subprocess.run(['nice','-n','10',PY,'-c',code], capture_output=True, text=True, timeout=1700)
print(r.stdout)
print("ERR", r.stderr[-300:] if r.stderr else "")
