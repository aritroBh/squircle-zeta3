import os, json
os.chdir("research/08_theory/K4b_cubo_finish")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 90
J16 = json.load(open("tetra_Sdk_h16_smoke.json"))
J32 = json.load(open("tetra_Sdk_h32.json"))
J48 = json.load(open("tetra_Sdk_h48.json"))
for region in ["S_A","S_B","S_tot"]:
    print(f"--- {region} ---")
    for k in range(5):
        v16, v32, v48 = mpf(J16[region][k]), mpf(J32[region][k]), mpf(J48[region][k])
        d1, d2 = abs(v32-v16), abs(v48-v32)
        ratio = d2/d1 if d1>0 else mpf(0)
        print(f"k={k}: |h32-h16|={mp.nstr(d1,3)}  |h48-h32|={mp.nstr(d2,3)}  ratio={mp.nstr(ratio,2)}")
# region A k=0 vs exact 3√3/4 (h48)
err = abs(mpf(J48["S_A"][0]) - 3*mp.sqrt(3)/4)
print("region A k=0 err vs 3√3/4 (h48):", mp.nstr(err,3))
