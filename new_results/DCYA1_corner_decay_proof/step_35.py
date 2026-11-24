import time, sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
import corner_series as cs
# Test M_60(pi/2) at increasing dps to find where corner_coeffs stabilizes
# True M_60 ~ -9.36e-9 (from FFT). Let me check corner_coeffs.
for dps in [35,50,70,100]:
    mp.mp.dps=dps
    t=time.time()
    cc=cs.corner_coeffs(mp.pi/2, 60)
    print(f"dps={dps}: M_60={mp.nstr(cc[60],18)} M_40={mp.nstr(cc[40],14)} t={round(time.time()-t,1)}s")
