forced910 = r'''"""forced910.py -- FORCED ascending recursion at depths 9, 10, symbolic in (d,m), fit-free.
Builds on forced_F78.pkl (delta <= 8, PROVED closing). F(delta,0) := a_delta (tow5lib, PROVED),
F(delta,s+1) := Lambda(delta,s) - F(delta,s). Closing: F(delta, floor(delta/2)+1) == 0.
Usage: python forced910.py 9   |   python forced910.py 10
Saves the full family (delta <= 10, cancel'd) to forced_F910.pkl incrementally.
"""
import sys, os, time, pickle
sys.dont_write_bytecode = True
work = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work)
import sympy as sp
import tow5lib as T5
import tow2lib as T2
d, m = T2.d, T2.m
I = sp.I
tau1 = T5.tau1f; tau2 = T5.tau2f

which = int(sys.argv[1])

# load full existing family: delta<=8 from forced_F78.pkl, plus any already-built 9/10
Fz = {}
raw78 = pickle.load(open(os.path.join(work, 'forced_F78.pkl'), 'rb'))
for k, v in raw78.items():
    Fz[k] = sp.sympify(v, locals={'d': d, 'm': m, 'I': sp.I})
outpath = os.path.join(work, 'forced_F910.pkl')
if os.path.exists(outpath):
    raw9 = pickle.load(open(outpath, 'rb'))
    for k, v in raw9.items():
        Fz[k] = sp.sympify(v, locals={'d': d, 'm': m, 'I': sp.I})

def F(dp, s):
    if s > dp//2: return sp.Integer(0)
    return Fz[(dp, s)]

def Lam(delta, s):
    tot = sp.Integer(0)
    for dp in range(0, delta+1):
        fA = F(dp, s).subs({d: d+1, m: m+1}, simultaneous=True)
        fB = F(dp, s).subs({d: d+1, m: m}, simultaneous=True)
        tot += I*(2*d-1)/(2*d-dp-1) * (tau1(delta-dp)*fA - tau2(delta-dp)*fB)
    return tot

def save():
    pickle.dump({k: str(v) for k, v in Fz.items()}, open(outpath, 'wb'))

t0 = time.time()
for delta in range(9, which+1):
    if (delta, 0) in Fz and (delta, delta//2) in Fz and Fz.get(('closed', delta), False):
        print(f'delta={delta} already done, skipping', flush=True); continue
    Fz[(delta, 0)] = sp.cancel(T5.arho(delta))
    for s in range(0, delta//2+1):
        ts = time.time()
        nxt = sp.together(Lam(delta, s) - F(delta, s))
        num = sp.expand(sp.numer(nxt))
        is0 = (num == 0)
        if s+1 <= delta//2:
            Fz[(delta, s+1)] = sp.cancel(nxt)
            print(f'  F({delta},{s+1}) built ({round(time.time()-ts,1)}s, total {round(time.time()-t0,1)}s)', flush=True)
        else:
            Fz[('closed', delta)] = is0
            print(f'  CLOSING delta={delta}: F({delta},{s+1}) numerator == 0 : {is0} ({round(time.time()-ts,1)}s, total {round(time.time()-t0,1)}s)', flush=True)
        save()
    print(f'delta={delta} done ({round(time.time()-t0,1)}s)', flush=True)
closed = {dl: Fz.get(('closed', dl)) for dl in range(9, which+1)}
print(f'VERDICT non-adjacent depth {which}: closing = {closed}', flush=True)
'''
work = '<project>/research/08_theory/TOW5_tower_depth10'
with open(work+'/forced910.py','w') as f:
    f.write(forced910)
# Also need forced_F78.pkl in the workdir (copy from TOW4)
import shutil
shutil.copy('<project>/research/08_theory/TOW4_tower_depth8/forced_F78.pkl', work+'/forced_F78.pkl')
print("forced910.py written, forced_F78.pkl copied")

import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run(['nice','-n','10',PY,'forced910.py','9'], cwd=work, capture_output=True, text=True, timeout=1700)
print(r.stdout)
print("ERR", r.stderr[-800:] if r.stderr else "")
print("elapsed", round(time.time()-t0,1), "rc", r.returncode)
with open(work+'/progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: forced910.py 9 rc={r.returncode} ({round(time.time()-t0,1)}s)\n")
