code = r'''
from Uq_direct import *
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, mpf
import time, json
setdps(80)
s = 4
t0=time.time()
data = {int(k): v for k, v in json.load(open("Uq_data_s4.json")).items()}
for p_ in [384, 512]:
    K = 6*p_
    ta=time.time()
    u = Uq_direct(p_, s, 1, K=K, R=10)
    data[p_] = [str(u)]
    print(f"p={p_}: U1={nstr(u, 20)}  D={nstr(u*power(p_,s-1), 20)}  ({time.time()-ta:.0f}s)")
json.dump({str(k):v for k,v in data.items()}, open("Uq_data_s4.json","w"))
print("total", time.time()-t0)
'''
exec(code)
