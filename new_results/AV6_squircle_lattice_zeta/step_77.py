code = r'''
from squircle import setdps
from Uq_direct import Uq_direct
from mpmath import mp, nstr, power
import time, json
setdps(80)
s = 4
t0=time.time()
data = {}
for p_ in [32, 48, 64, 96, 128, 192, 256]:
    K = 6*p_
    row = []
    for q in [1,2,3]:
        ta=time.time()
        u = Uq_direct(p_, s, q, K=K, R=10)
        row.append(str(u))
    data[p_] = row
    print(f"p={p_}: {[nstr(mpf(x),12) for x in row]}  ({time.time()-t0:.0f}s)")
json.dump({str(k):v for k,v in data.items()}, open("Uq_data_s4_v2.json","w"))
print("total", time.time()-t0)
'''
exec(code)
with open("progress.log","a") as f: f.write("step_44: I0q sign bug fixed; Uq data recomputed\n")
