import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
import mpmath as mp, json
mp.mp.dps = 50
prof = [(mp.mpf(s), mp.mpf(a), mp.mpf(b), mp.mpf(c)) for s,a,b,c in json.load(open("ckpt_prof_m3.json"))]
alpha = 1/mp.sqrt(3); ln2 = mp.log(2)
print("s: D1_extracted vs ln2 + (1/2)log(1+e^{-4 alpha |s|})")
for s, a, b, c in prof:
    if s in [mp.mpf(x)/2 for x in range(0, 13)]:
        pred = ln2 + mp.mpf('0.5')*mp.log(1+mp.exp(-4*alpha*abs(s)))
        print(f"  s={mp.nstr(s,5)}: ext={mp.nstr(a,18)} pred={mp.nstr(pred,18)} diff={mp.nstr(a-pred,6)}")
# also check alternative: ln2 + log(1+e^{-2 alpha |s|}) - log c with c=1/2:
print("alternative A=1, rate 2alpha:")
for s, a, b, c in prof:
    if s in [mp.mpf(x)/2 for x in range(0, 13)]:
        pred = ln2 + mp.log(1+mp.exp(-2*alpha*abs(s))) 
        print(f"  s={mp.nstr(s,5)}: ext={mp.nstr(a,18)} pred={mp.nstr(pred,18)} diff={mp.nstr(a-pred,6)}")
