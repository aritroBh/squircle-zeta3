import json, time
from mpmath import mp, mpf, log, pi, sqrt, exp, factorial, sin, cos
t0=time.time()
mp.dps = 120
z0 = (-log(2)+pi*1j)/2
R = abs(z0); th = mp.arg(z0)
j120 = json.load(open("../SRF3_octahedron_radius_proof/SRF3_hi_coeffs_K120.json"))
A_tau = [mpf(x) for x in j120['eta']]
A_d   = [mpf(x) for x in j120['d']]
b_tau = [A_tau[k]/factorial(k) for k in range(121)]
b_d   = [A_d[k]/factorial(k) for k in range(121)]

# ---------- (A) free 5-parameter fit: b_k = C k^beta Rf^{-k} cos(k*thf + phi), k in [70,120]
import itertools
def resid2(params, b):
    C, beta, Rf, thf, phi = params
    s = mpf(0)
    for k in range(70,121):
        m = C*k**beta*Rf**(-k)*cos(k*thf+phi)
        s += ((m-b[k])/b[k])**2
    return float(s)
# grid + refine
best=None
for beta0 in [-1.55,-1.5,-1.45]:
    for Rf0 in [1.58,1.61,1.64]:
        for thf0 in [1.76,1.79,1.82]:
            for phi0 in [-2.4,-2.25,-2.1]:
                v = resid2([1.4,beta0,Rf0,thf0,phi0], b_tau)
                if best is None or v<best[0]: best=(v,[1.4,beta0,Rf0,thf0,phi0])
print("grid best:", best[0], best[1])
# simplex refine (Nelder-Mead by hand: use mp findroot on gradient? use scipy)
import numpy as np
from scipy.optimize import minimize
x0 = np.array(best[1], dtype=float)
def obj(x): return resid2(list(x), b_tau)
r = minimize(obj, x0, method='Nelder-Mead', options={'xatol':1e-12,'fatol':1e-30,'maxiter':20000})
C,beta,Rf,thf,phi = r.x
print(f"FREE FIT tau: C={C:.10f} beta={beta:.8f} R={Rf:.10f} theta={thf:.10f} phi={phi:.10f}  resid={r.fun:.3e}")
print(f"theory:       C={float(sqrt(3*pi*R/8)):.10f} beta=-1.5 R={float(R):.10f} theta={float(th):.10f} phi={float(pi-th/2):.10f}")
print(f"A=1/R fit: {1/Rf:.10f} vs theory {float(1/R):.10f}  [{time.time()-t0:.0f}s]")
