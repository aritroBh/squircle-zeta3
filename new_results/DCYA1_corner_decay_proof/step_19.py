import time, sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps = 40
import corner_series as cs

# Test: M_0 at pi/2 (known: M_0(pi/2) = -(1+c)log(1+c)-(1-c)log(1-c), c=cos(pi/4)=sqrt2/2)
t=time.time()
a = mp.pi/2
coeffs = cs.corner_coeffs(a, 3)
print("time for N=3:", round(time.time()-t,2),"s")
for k,v in enumerate(coeffs):
    print(f"M_{k}(pi/2) = {mp.nstr(v,30)}")
# Check M_0 formula
c = mp.cos(mp.pi/4)
M0_pred = -(1+c)*mp.log(1+c)-(1-c)*mp.log(1-c)
print("M_0 formula:", mp.nstr(M0_pred,30))
print("diff:", mp.nstr(coeffs[0]-M0_pred,5))
