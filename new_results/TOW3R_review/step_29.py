work = '<project>/research/08_theory/TOW3R_review/'
src = '''import sys, os, time
work = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work+'/tow3_copy'); os.chdir(work+'/tow3_copy')
import sympy as sp
import tow2lib as T
d, m = T.d, T.m
import pickle
L = pickle.load(open('L_parsed.pkl','rb'))
print('L_parsed keys:', sorted(L.keys()))
for r in range(1, 7):
    t0 = time.time()
    diff = sp.expand(L[r] - T.L(r))
    print(f'r={r}: L_parsed == tow2lib.L(r) symbolically: {diff==0}  ({time.time()-t0:.1f}s)', flush=True)
# also L_simplified for r=7,8 (not needed at depth<=6; random points only, 20 pts)
Ls = pickle.load(open('L_simplified.pkl','rb'))
import random
random.seed(7)
for r in [7,8]:
    Lr = sp.sympify(Ls[r], locals={'d': d, 'm': m, 'I': sp.I})
    ok = True
    for _ in range(20):
        dv, mv = random.randint(8,60), random.randint(1,7)
        ok &= bool(sp.expand(Lr.subs({d:dv,m:mv}) - T.L(r).subs({d:dv,m:mv}))==0)
    print(f'r={r}: L_simplified == tow2lib.L(r) at 20 random points: {ok}', flush=True)
'''
open(work+'step_02.py','w').write(src)
import subprocess, datetime
PY='<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,work+'step_02.py'],capture_output=True,text=True,timeout=1700)
print(r.stdout, r.stderr[-2000:])
with open(work+'progress.log','a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} step_02 L_parsed==tow2lib.L(r) symbolic r=1..6 -> see step_02.out\n")
open(work+'step_02.out','w').write(r.stdout+r.stderr)
