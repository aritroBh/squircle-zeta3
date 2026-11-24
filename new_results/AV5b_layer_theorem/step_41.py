import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
step = r'''
import mpmath as mp
mp.mp.dps = 60
import json, sys

def S(phi, m, p):
    return mp.fsum(abs(mp.cos(phi - i*mp.pi/m))**p for i in range(1, m+1))
def rho(phi, m, p):
    return S(phi, m, p)**(-mp.mpf(1)/p)
def speed(phi, m, p):
    r = rho(phi, m, p); dr = mp.diff(lambda ph: rho(ph, m, p), phi)
    return mp.sqrt(r*r + dr*dr)
def P_of_p(m, p):
    ks = [(mp.pi/2 + mp.pi*i/m) % (2*mp.pi) for i in range(2*m)]
    ms = [(mp.pi*i/m) % (2*mp.pi) for i in range(2*m)]
    pts = sorted(set(ks+ms))
    pts = [x for x in pts if mp.mpf('1e-12') < x < 2*mp.pi-mp.mpf('1e-12')]
    b = [mp.mpf(0)] + pts + [2*mp.pi]
    tot = mp.mpf(0)
    for a, bb in zip(b[:-1], b[1:]):
        mid = (a+bb)/2
        tot += mp.quad(lambda ph: speed(ph, m, p), [a, mid, bb])
    return tot

if __name__ == "__main__":
    mm = int(sys.argv[1])
    P1 = P_of_p(mm, mp.mpf(1))
    out = []
    for dstr in sys.argv[2].split(","):
        d = mp.mpf(dstr)
        v = P_of_p(mm, 1+d) - P1
        out.append((dstr, v))
        print(f"m={mm} d={dstr}: dP={mp.nstr(v,35)}", flush=True)
    json.dump([[d, str(v)] for d, v in out], open(f"ckpt_true_m{mm}.json","w"))
'''
with open("step_30.py","w") as f: f.write(step)
import subprocess, time
PY = "<project>/.venv/bin/python"
t0 = time.time()
r = subprocess.run([PY,"step_30.py","3","0.048,0.024,0.012,0.006,0.003,0.0015,0.00075,0.000375"],capture_output=True,text=True,timeout=1400)
print(r.stdout, r.stderr[-300:])
print("elapsed", time.time()-t0)
