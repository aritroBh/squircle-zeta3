import sys, os, subprocess, datetime, time, shutil
work = '<project>/research/08_theory/TOW3R_review/'
base = '<project>/research/08_theory/'
# copy TOW2 base scripts used at depths<=4
for f in ['prove_adj34.py','step_48.py','step_49.py','step_52.py']:
    if os.path.exists(base+'TOW3_tower_depth6/'+f):
        shutil.copy(base+'TOW3_tower_depth6/'+f, work+'tow3_copy/'+f)
    else:
        print('missing', f)
PY='<project>/.venv/bin/python'
# quick: conj(tau) rules + check_tau
src = '''import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/tow3_copy'); os.chdir(os.path.dirname(os.path.abspath(__file__))+'/tow3_copy')
import sympy as sp, tow2lib as T
print('check_tau(n<=7):', T.check_tau(7))
ok = True
for n in range(0, 8):
    c1 = sp.simplify(sp.conjugate(T.tau1(n)) - (-1)**n * T.tau2(n))
    c2 = sp.simplify(sp.conjugate(T.tau2(n)) - (-1)**n * T.tau1(n))
    ok &= (c1==0) and (c2==0)
    if c1!=0 or c2!=0: print('FAIL', n, c1, c2)
print('conj(tau1_n)=(-1)^n tau2_n and conj(tau2_n)=(-1)^n tau1_n for n=0..7:', ok)
'''
open(work+'step_15.py','w').write(src)
r = subprocess.run(['nice','-n','10',PY,work+'step_15.py'],capture_output=True,text=True,timeout=600)
print(r.stdout, r.stderr[-500:])
open(work+'step_15.out','w').write(r.stdout+r.stderr)
with open(work+'progress.log','a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} step_15 tau conjugation rules n<=7\n")
