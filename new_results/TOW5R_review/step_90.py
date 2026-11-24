code = '''"""step_19b (TOW5R): MASTER identity feed(U(-1)) - a_delta == 0 at delta = 9, 10 (own rebuild).
U(dp) = sum_s (-1)^s F(dp,s); feed = i(2d-1) sum_dp [tau1(delta-dp) U(dp)(d+1,m+1)
        - tau2(delta-dp) U(dp)(d+1,m)]/(2d-dp-1). Exact cancel'd difference must be 0.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/tow5_copy')
import sympy as sp
import tow2lib as T2, tow5lib as T5
d, m = T2.d, T2.m
i = sp.I
t0 = time.time()
def U(dp): return sum((-1)**s * T5.forced(dp, s) for s in range(0, dp//2+1))
def feedU(delta):
    tot = sp.Integer(0)
    for dp in range(0, delta+1):
        uA = U(dp).subs({d: d+1, m: m+1}, simultaneous=True)
        uB = U(dp).subs({d: d+1, m: m}, simultaneous=True)
        tot += i*(2*d-1)/(2*d-dp-1) * (T5.tau1f(delta-dp)*uA - T5.tau2f(delta-dp)*uB)
    return tot
ok = True
for dl in [9, 10]:
    t1 = time.time()
    r = sp.cancel(sp.together(feedU(dl) - T5.arho(dl)))
    ok &= (r == 0)
    print(f'MASTER delta={dl}: feed(U(-1)) - a_delta == 0 : {r == 0} ({round(time.time()-t1,1)}s)', flush=True)
print('STEP_19B VERDICT:', 'PASS' if ok else 'FAIL', f'({round(time.time()-t0,1)}s)', flush=True)
assert ok
'''
wd = '<project>/research/08_theory/TOW5R_review'
open(wd+'/step_19b.py','w').write(code)
import subprocess, datetime
PY = '<project>/.venv/bin/python'
t0 = datetime.datetime.now()
r = subprocess.run(['nice','-n','10',PY,'step_19b.py'], cwd=wd, capture_output=True, text=True, timeout=1700)
print(r.stdout)
print('STDERR tail:', r.stderr[-300:], 'rc:', r.returncode, 'elapsed:', datetime.datetime.now()-t0)
with open(wd+'/progress.log','a') as f:
    f.write(f"{datetime.datetime.now():%H:%M:%S} step_19b: MASTER x=-1 identity dl=9,10 rc={r.returncode} ({(datetime.datetime.now()-t0).seconds}s)\n")
