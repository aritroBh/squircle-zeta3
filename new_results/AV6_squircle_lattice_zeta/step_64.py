code = r'''
from Uq_direct import *
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, mpf, zeta, gamma
import time, json
setdps(100)
s = 4
t0=time.time()
data = {int(k): v for k, v in json.load(open("Uq_data_s4.json")).items()}
for p_ in [160, 192, 256]:
    K = 6*p_
    row = []
    for q in [1,2,3]:
        ta=time.time()
        u = Uq_direct(p_, s, q, K=K, R=10)
        row.append(str(u))
        print(f"p={p_} q={q}: {nstr(u, 20)}  ({time.time()-ta:.0f}s)")
    data[p_] = row
json.dump({str(k):v for k,v in data.items()}, open("Uq_data_s4.json","w"))
print("total", time.time()-t0)
'''
exec(code)
