work = '<project>/research/08_theory/TOW3R_review/'
src = '''import sys, os, time
work = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work+'/tow3_copy'); os.chdir(work+'/tow3_copy')
import sympy as sp
import tow2lib as T
d, m = T.d, T.m
import pickle
L = pickle.load(open('L_parsed.pkl','rb'))
# correct equality test: difference as element of Q(i)(d,m)[sqrt2]; evaluate sqrt2 to 60 digits
s2 = sp.sqrt(2)
for r in range(1, 7):
    t0 = time.time()
    diff = sp.expand(L[r] - T.L(r))
    # reduce: write diff = A + B*sqrt(2) by collecting powers of sqrt(2)
    pol = sp.Poly(diff, s2)
    terms = pol.terms()
    ok = all(sp.expand(c)==0 for ((e,),c) in terms if e%2==0) # even powers give rational*2^k
    parts = {}
    for ((e,), c) in terms:
        parts.setdefault(e%2, sp.Integer(0))
        parts[e%2] += c * 2**(e//2)
    ok = (sp.expand(parts.get(0,0))==0) and (sp.expand(parts.get(1,0))==0)
    print(f'r={r}: L_parsed == tow2lib.L(r) (Q(i)(d,m)[sqrt2] reduction): {ok}  ({time.time()-t0:.1f}s)', flush=True)
    if not ok:
        print('  residual even part:', str(parts.get(0,0))[:200]); print('  residual sqrt2 coeff:', str(parts.get(1,0))[:200])
# high-precision numeric spot confirmation r=1..8 at random (d,m)
import pickle, random
from mpmath import mp, mpf
mp.dps = 60
Ls = pickle.load(open('L_simplified.pkl','rb'))
random.seed(3)
for r in range(1,9):
    Lr = sp.sympify(Ls[r], locals={'d': d, 'm': m, 'I': sp.I})
    ok = True
    for _ in range(5):
        dv, mv = random.randint(8,60), random.randint(1,7)
        v1 = complex(sp.N(Lr.subs({d:dv,m:mv}), 50))
        v2 = complex(sp.N(T.L(r).subs({d:dv,m:mv}), 50))
        ok &= abs(v1-v2) < 1e-40
    print(f'r={r}: L_simplified == tow2lib.L(r) numeric 50-digit at 5 pts: {ok}', flush=True)
'''
open(work+'step_03.py','w').write(src)
import subprocess, datetime
PY='<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,work+'step_03.py'],capture_output=True,text=True,timeout=1700)
print(r.stdout, r.stderr[-1500:])
open(work+'step_03.out','w').write(r.stdout+r.stderr)
with open(work+'progress.log','a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} step_03 redo L equality with sqrt2-aware reduction\n")
