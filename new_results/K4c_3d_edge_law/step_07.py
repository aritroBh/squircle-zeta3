"""step_07.py -- MODEL (a) h=1/96 final."""
import os, json, time
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 80
from k4c_driver import chamber_S
s2 = mp.sqrt(2)
NA = [[1/s2,1/s2,0],[1/s2,-1/s2,0],[1/s2,0,1/s2],[1/s2,0,-1/s2],[0,1/s2,1/s2],[0,1/s2,-1/s2]]
t0=time.time()
S, st = chamber_S(NA, 4, mpf(1)/96, mpf(1)/96, mpf('4.5'), mpf('4.5'), HI2=mpf(1)/96)
el = time.time()-t0
print(f"h96 ({el:.0f}s) stats={st}", flush=True)
for k in range(5):
    print(f"   [d^{k}]S = {mp.nstr(S[k],58)}", flush=True)
json.dump({"Sdk":[mp.nstr(v,70) for v in S],"stats":st}, open("modelA_h96.json","w"), indent=1)
J64 = json.load(open("modelA_h64.json")); S64=[mpf(s) for s in J64["Sdk"]]
J48 = json.load(open("modelA_h48.json")); S48=[mpf(s) for s in J48["Sdk"]]
print("trend:")
for k in range(5):
    d1 = abs(S48[k]-S64[k]); d2 = abs(S64[k]-S[k])
    r = d1/d2 if d2>0 else mpf('inf')
    est = d2/(r-1) if r>1 else d2
    print(f"[d^{k}] |h64-h48|={mp.nstr(d1,3)} |h96-h64|={mp.nstr(d2,3)} ratio={mp.nstr(r,3)} -> h96 err est {mp.nstr(est,3)}")
