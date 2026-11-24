
import os
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
step = r'''"""WT6S step E (part 1): N_{(i,j),r} to ~150 digits, two independent quadrature methods.

N_{(i,j),r} = int_b^a l1^i l2^j dv/(v-r), l1=log t, l2=log(1-t), t=(v^2+2v-1)/(4v).
Method RAW: tanh-sinh over [b,1/2,1,2,a] (M5a setup).
Method SUB: split [b,1]+[1,a], power substitution v=b+(1-b)u^p / v=a-(a-1)u^p (p=12)
  removing the endpoint log singularities, tanh-sinh on [0,1].
Cross-check agreement; store n6_values.json."""
import json, time
import mpmath as mp

mp.mp.dps = 200
r2 = mp.sqrt(2); aN = 1 + r2; bN = r2 - 1
ONE = mp.mpf(1)
POLE = {"0": mp.mpf(0), "a": aN, "-a": -aN, "-b": -bN}
KEYS = [(5, 0), (4, 1), (3, 2)]
PRS = ["0", "a", "-a", "-b"]
P = 12

def l1(v):
    return mp.log((v - bN) * (v + aN) / (4 * v))
def l2(v):
    return mp.log((aN - v) * (v + bN) / (4 * v))

def N_raw(i, j, rname):
    r = POLE[rname]
    def f(v):
        return l1(v)**i * l2(v)**j / (v - r)
    val, err = mp.quad(f, [bN, mp.mpf(1)/2, 1, 2, aN], error=True, maxdegree=14)
    return val, err

def N_sub(i, j, rname):
    r = POLE[rname]
    # piece 1: v = b + (1-b) u^p, u in [0,1]
    def f1(u):
        v = bN + (1 - bN) * u**P
        dv = (1 - bN) * P * u**(P - 1)
        return l1(v)**i * l2(v)**j / (v - r) * dv
    # piece 2: v = a - (a-1) u^p, integrate u 0->1 (dv negative, take -)
    def f2(u):
        v = aN - (aN - 1) * u**P
        dv = -(aN - 1) * P * u**(P - 1)
        return -l1(v)**i * l2(v)**j / (v - r) * dv
    v1, e1 = mp.quad(f1, [0, 1], error=True, maxdegree=14)
    v2, e2 = mp.quad(f2, [0, 1], error=True, maxdegree=14)
    return v1 + v2, max(e1, e2)

def main():
    t0 = time.time()
    vals = {}
    for (i, j) in KEYS:
        for rname in PRS:
            vr, er = N_raw(i, j, rname)
            vs, es = N_sub(i, j, rname)
            agree = abs(vr - vs)
            vals[f"{i},{j}|{rname}"] = mp.nstr(vs, 160)
            print(f"N_{{({i},{j}),{rname}}}: sub-raw |d|={mp.nstr(agree,3)} "
                  f"err_r={mp.nstr(er,2)} err_s={mp.nstr(es,2)} ({time.time()-t0:.0f}s)", flush=True)
    json.dump({'dps': 160, 'values': vals}, open('n6_values.json', 'w'), indent=1)
    with open('progress.log', 'a') as f:
        f.write('step_E1: 12 N_{(i,j),r} at 160 digits, raw vs sub cross-checked\n')
    print("saved n6_values.json")

if __name__ == '__main__':
    main()
'''
open(os.path.join(WD, 'step_07.py'), 'w').write(step)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice', '-n', '10', PY, '-B', 'step_07.py'], capture_output=True, text=True, cwd=WD, timeout=1500)
print(r.stdout)
print("STDERR:", r.stderr[-2000:])
