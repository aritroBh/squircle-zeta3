# ======== STEP 1 indep. check: Cauchy formula for a_n from P_direct quadrature (complex p) ========
from mpmath import mp, mpf, mpc, quad, sqrt, exp, log, pi, cos, sin
import json, time
mp.dps = 70
def P_complex(p):
    e = 1/p
    f = lambda t: sqrt((e*t**(e-1))**2 + (e*(1-t)**(e-1))**2)
    return 8*quad(f, [0, mpf(1)/16, mpf(1)/4, mpf(1)/2])
t0 = time.time()
print("timing one complex P:", P_complex(1+mpc('0.09')), "  sec:", round(time.time()-t0, 2))
