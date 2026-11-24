import time, sys
sys.path.insert(0, '.')
import mpmath as mp
import cor7_lib as L

c3 = mp.mpf(1)/3
for dps in (400, 700):
    mp.mp.dps = dps
    t0 = time.time()
    M3q = L.corner_coeffs_hp(2*mp.acos(c3), 3)[3]
    M3ref, _, _, _ = L.M3_weights(c3)
    print(f"dps {dps}: quad {time.time()-t0:.0f}s  |quad-HPL| = {mp.nstr(abs(M3q-M3ref),3)}", flush=True)
