import time, sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=35
import corner_series as cs
t=time.time()
cc=cs.corner_coeffs(mp.pi/2, 100)
print("corner_coeffs N=100 time:",round(time.time()-t,1),"s")
for k in [60,80,100]:
    print(f"M_{k}(pi/2) = {mp.nstr(cc[k],20)}")
# known M_100 from K2 = -1.772316807359054100088703e-10
print("K2 says M_100 = -1.772316807359054100088703e-10")
