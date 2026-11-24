work = '<project>/research/08_theory/TOW5_tower_depth10'
# Update tow5lib.py: R5 from N5_proved.pkl (clean form N_5/(3 Pi_5))
src = open(work+'/tow5lib.py').read()
old = '''# ---- R_5 (PROVED, prove_S5.py) ----
_raw5 = pickle.load(open(os.path.join(_HERE, 'R5_raw.pkl'), 'rb'))
_num5 = sp.sympify(_raw5['num5'], locals={'d': d, 'm': m, 'I': sp.I})
_den5 = sp.sympify(_raw5['den5'], locals={'d': d, 'm': m, 'I': sp.I})
R5 = sp.cancel(_num5 / _den5)'''
new = '''# ---- R_5 (PROVED, prove_S5.py): R_5 = N_5/(3 Pi_5) ----
_p5 = pickle.load(open(os.path.join(_HERE, 'N5_proved.pkl'), 'rb'))
_N5 = sp.sympify(_p5['N5'], locals={'d': d, 'm': m, 'I': sp.I})
def _Pi5(dd=None):
    dd = dd if dd is not None else d
    out = sp.Integer(1)
    for j in range(1, 6): out *= (2*dd - 2*j - 1)
    return out
R5 = sp.cancel(_N5 / (3*_Pi5()))'''
assert old in src
src = src.replace(old, new)
open(work+'/tow5lib.py','w').write(src)
print("tow5lib updated")

# Re-run rho(9), rho(10) check to confirm nothing changed
import subprocess, time
PY = '<project>/.venv/bin/python'
code = '''
import sys, pickle, time
from fractions import Fraction as Fr
sys.path.insert(0, %r)
import sympy as sp
import tow5lib as T5
d_, m_ = T5.d, T5.m
i = sp.I
def qsi_div(x, y):
    def gm(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
    def gd(a,b):
        den = b[0]**2+b[1]**2; return ((a[0]*b[0]+a[1]*b[1])/den, (a[1]*b[0]-a[0]*b[1])/den)
    A=(x[0],x[2]); B=(x[1],x[3]); C=(y[0],y[2]); D=(y[1],y[3])
    C2 = (gm(C,C)[0]-2*gm(D,D)[0], gm(C,C)[1]-2*gm(D,D)[1])
    nA = (gm(A,C)[0]-2*gm(B,D)[0], gm(A,C)[1]-2*gm(B,D)[1])
    nB = (gm(B,C)[0]-gm(A,D)[0], gm(B,C)[1]-gm(A,D)[1])
    return (gd(nA, C2), gd(nB, C2))
RHO = pickle.load(open(%r + "/data_rho_D10_d7_40.pkl","rb"))
t0=time.time()
for delta in [9, 10]:
    aa = T5.arho(delta)
    bad = 0; n = 0
    for dd in RHO:
        if 2*dd-4 < delta: continue
        for key in RHO[dd]:
            p,q = key
            rv = qsi_div(RHO[dd][key][delta], RHO[dd][key][0])
            cv = sp.expand(aa.subs({d_: dd, m_: p}))
            exp = (Fr(sp.re(cv)), Fr(sp.im(cv)))
            n += 1
            if rv[1] != (Fr(0),Fr(0)) or rv[0] != exp: bad += 1
    print(f"rho({delta})/rho(0) == a_{delta} (tow5lib): {n} cells, {bad} bad ({round(time.time()-t0,1)}s)", flush=True)
''' % (work, work)
t0=time.time()
r = subprocess.run(['nice','-n','10',PY,'-c',code], capture_output=True, text=True, timeout=1700)
print(r.stdout, "ERR", r.stderr[-500:] if r.stderr else "", "elapsed", round(time.time()-t0,1))
