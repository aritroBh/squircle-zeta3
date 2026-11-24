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
def U(dp): return sum((-1)**s * T5.forced(dp, s) for s in range(0, dp//2+1))
def feedU(delta):
    tot = sp.Integer(0)
    for dp in range(0, delta+1):
        uA = U(dp).subs({d: d+1, m: m+1}, simultaneous=True)
        uB = U(dp).subs({d: d+1, m: m}, simultaneous=True)
        tot += i*(2*d-1)/(2*d-dp-1) * (T5.tau1f(delta-dp)*uA - T5.tau2f(delta-dp)*uB)
    return tot
for dl in [9, 10]:
    t0=time.time()
    rem = sp.cancel(sp.together(feedU(dl) - T5.arho(dl)))
    print(f"MASTER x=-1 at delta={dl}: feed(U(-1)) - a_delta == 0 : {rem==0}  ({round(time.time()-t0,1)}s)", flush=True)
''' % (work,)
t0=time.time()
r = subprocess.run(['nice','-n','10',PY,'-c',code], capture_output=True, text=True, timeout=1700)
print(r.stdout)
print("ERR", r.stderr[-500:] if r.stderr else "", "elapsed", round(time.time()-t0,1))
with open(work+'/progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: MASTER x=-1 at delta=9,10 rc={r.returncode} ({round(time.time()-t0,1)}s)\n")
