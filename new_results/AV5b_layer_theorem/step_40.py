import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
step = r'''
import mpmath as mp
import json
mp.mp.dps = 80

def load(mm, full=0):
    data = json.load(open(f"ckpt_arc_m{mm}_f{full}.json"))
    pts = [(mp.mpf(d), mp.mpf(v)) for d, v in data]
    pts.sort(key=lambda x: -x[0])
    return pts[:10]

# m=3: subtract A2 d^2 (A2 from before), look at (DP - A2 d^2)/d^3 at each d
pts = load(3)
A2 = mp.mpf("4.26972531291317353796562809865")
print("m=3: (DP - A2 d^2)/d^3 at each d:")
for d, v in pts:
    print(f"  d={mp.nstr(d,8)}: {mp.nstr((v - A2*d**2)/d**3, 20)}")
# and (DP)/d^2 raw:
print("m=3: DP/d^2 at each d:")
for d, v in pts:
    print(f"  d={mp.nstr(d,8)}: {mp.nstr(v/d**2, 20)}")
'''
with open("step_29.py","w") as f: f.write(step)
exec(step)
