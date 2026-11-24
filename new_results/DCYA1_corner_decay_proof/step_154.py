import json, sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
res=json.load(open('ei_sweep.json'))
# profile: bin by N ranges and show max in each
import math
res.sort()
# sample every ~Nth cell
for target in [100,200,400,800,1600,3200,6400,10000,15000,20000,30000,50000,100000,200000,400000,700000,1000000]:
    # find cell containing target
    cand=[r for r in res if r[0]<=target<=r[1]]
    if cand:
        print(f"N~{target}: E_I/sc={cand[0][2]:.2f} (cell {cand[0][0]},{cand[0][1]})")
# also thin cell at 20000
from arb_ei import E_I_sc, iv
print("thin [20000,20050]:", E_I_sc(iv("20000","20050")))
print("thin [1600,1604]:", E_I_sc(iv("1600","1604")))
