"""step_06.py -- MODEL (a) h=1/64 for convergence trend."""
import os, json, time
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 80
from k4c_driver import chamber_S
s2 = mp.sqrt(2)
NA = [[1/s2,1/s2,0],[1/s2,-1/s2,0],[1/s2,0,1/s2],[1/s2,0,-1/s2],[0,1/s2,1/s2],[0,1/s2,-1/s2]]
t0=time.time()
S, st = chamber_S(NA, 4, mpf(1)/64, mpf(1)/64, mpf('4.5'), mpf('4.5'), HI2=mpf(1)/64)
el = time.time()-t0
print(f"h64 ({el:.0f}s) stats={st}", flush=True)
for k in range(5):
    print(f"   [d^{k}]S = {mp.nstr(S[k],55)}", flush=True)
json.dump({"Sdk":[mp.nstr(v,70) for v in S],"stats":st}, open("modelA_h64.json","w"), indent=1)
J48 = json.load(open("modelA_h48.json")); S48=[mpf(s) for s in J48["Sdk"]]
J32 = json.load(open("modelA_h32.json")); S32=[mpf(s) for s in J32["Sdk"]]
print("convergence trend (diffs):")
for k in range(5):
    d1 = abs(S48[k]-S32[k]); d2 = abs(S[k]-S48[k])
    r = d1/d2 if d2>0 else mpf('inf')
    print(f"[d^{k}] |h48-h32|={mp.nstr(d1,3)}  |h64-h48|={mp.nstr(d2,3)}  ratio={mp.nstr(r,3)}")
