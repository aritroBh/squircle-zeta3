import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
step = r'''
import mpmath as mp
mp.mp.dps = 55
import json, sys

exec(open("av5b_locmodel.py").read().split("def extract")[0])

def dP_loc(mm, d):
    Zo = Loc(mm)
    p = 1+d
    W = mp.pi/(2*mm)
    tau = d/(1+d)
    f = lambda u: Zo.h(u, p) - Zo.tent(u)
    tot = mp.mpf(0)
    for a, b in [(-W, -20*tau), (-20*tau, mp.mpf(0)), (mp.mpf(0), 20*tau), (20*tau, W)]:
        tot += mp.quad(f, [a, b])
    return 2*mm*tot

if __name__ == "__main__":
    mm = int(sys.argv[1])
    out = []
    for dstr in sys.argv[2].split(","):
        d = mp.mpf(dstr)
        v = dP_loc(mm, d)
        out.append((dstr, v))
        print(f"m={mm}loc d={dstr}: dP={mp.nstr(v,35)}", flush=True)
    json.dump([[d, str(v)] for d, v in out], open(f"ckpt_dPloc_m{mm}.json","w"))
'''
with open("av5b_dploc.py","w") as f: f.write(step)
import subprocess, time
PY = "<project>/.venv/bin/python"
t0=time.time()
r = subprocess.run([PY,"av5b_dploc.py","3","0.048,0.024,0.012,0.006,0.003,0.0015,0.00075"],capture_output=True,text=True,timeout=1400)
print(r.stdout, r.stderr[-300:])
print("elapsed",time.time()-t0)
