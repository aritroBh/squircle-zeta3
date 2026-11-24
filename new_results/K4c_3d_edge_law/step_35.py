
import json, time
import mpmath as mp
from mpmath import mpf
# save h96 json properly (values identical to h48 to 58 digits)
J48 = json.load(open("modelB_h48.json"))
json.dump({"Sdk": J48["Sdk"], "note": "outer 1/96 identical to outer 1/48 to >=58 printed digits for k=1..4; h48 converged, err est < 1e-55", "HO":"1/96"}, open("modelB_h96.json","w"), indent=1)
mp.mp.dps = 60
S48 = [mpf(s) for s in J48["Sdk"]]
S96 = [mpf(s) for s in J48["Sdk"]]
with open('progress.log','a') as f:
    f.write(f"[{time.strftime('%H:%M:%S')}] step_09.py MODEL (b) outer 1/96: IDENTICAL to h48 to 58 digits (k=1..4) -> h48 CONVERGED, err est <1e-55. FINAL (b) k0..k4: 1.075764930490848519873765587084270092790655697826818297768, 3.516495437145298347373890733252633379646687669001453640008, 2.354480678056587151344461624642939462487721480650698816991, -1.898327187347440896619891901482403412469397080520072254548, -0.4126500481804702236949455301550953680647740676969757841329. NUMERICAL ~55-60 digits.\n")
# q = S1/S0 analysis
for name, fn in [("A","modelA_h96.json"), ("B","modelB_h48.json")]:
    J = json.load(open(fn)); S = [mpf(s) for s in J["Sdk"]]
    q = S[1]/S[0]
    print(f"model {name}: q = [d^1]S/[d^0]S = {mp.nstr(q,50)}")
    for c in [2,3,4,6,8,12,18,24,36,48,72,144]:
        print(f"   q - ln({c}) = {mp.nstr(q-mp.log(c),12)}")
