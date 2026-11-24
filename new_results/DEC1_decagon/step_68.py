import os, pickle, mpmath as mp, json
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 100
gold = pickle.load(open('P3_660.pkl','rb'))['P3']
v=json.load(open('research/08_theory/D15b_polygon_law/vals_m2_8_200.json'))
t5 = mp.mpf(v['5'][3]); t2 = mp.mpf(v['2'][3])
print("|gold - D15b 200d value| =", mp.nstr(abs(gold-t5),8))
sq = mp.sqrt(2)
P3_m2 = 6*sq*(mp.zeta(3)-mp.pi**2/3)
print("m=2: 6 sqrt2 (z3 - pi^2/3) =", mp.nstr(P3_m2,40))
print("m=2 json value          =", mp.nstr(t2,40))
print("diff:", mp.nstr(abs(P3_m2-t2),6))
with open('progress.log','a') as f:
    f.write(f"step: gold vs D15b 200d: {mp.nstr(abs(gold-t5),4)}; m=2 closed form check diff {mp.nstr(abs(P3_m2-t2),4)}\n")