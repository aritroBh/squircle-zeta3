import subprocess, time, os
work = '<project>/research/08_theory/TOW5_tower_depth10'
PY = '<project>/.venv/bin/python'
code = '''
import sys, time
sys.path.insert(0, %r)
import sympy as sp
import tow5lib as T5
d, m = T5.d, T5.m
for dl in [9, 10]:
    U = sp.together(sum((-1)**s * T5.forced(dp_, s) for dp_, s in [(dl, s) for s in range(0, dl//2+1)]))
    a = sp.together(T5.arho(dl))
    un, ud = sp.fraction(U)
    an, ad = sp.fraction(a)
    print(f"delta={dl}: den U == den a : {sp.expand(ud-ad)==0}; ", end="")
    # compare leading numerator terms in d
    Pu = sp.Poly(sp.expand(un), d); Pa = sp.Poly(sp.expand(an), d)
    du = sp.degree(sp.expand(un), d)
    same_lead = all(sp.expand(Pu.nth(dd_) - Pa.nth(dd_))==0 for dd_ in range(du, du-3, -1))
    print(f"numerator degree {du}; top-3 d-coeffs of U == a : {same_lead}")
    # first differing power
    for dd_ in range(du, -1, -1):
        if sp.expand(Pu.nth(dd_) - Pa.nth(dd_)) != 0:
            print(f"   first d-power where U != a: d^{dd_}")
            break
''' % (work,)
r = subprocess.run(['nice','-n','10',PY,'-c',code], capture_output=True, text=True, timeout=1700)
print(r.stdout)
print("ERR", r.stderr[-300:] if r.stderr else "")
