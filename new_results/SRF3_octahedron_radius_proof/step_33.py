
import os, json, math
import numpy as np
import mpmath as mp
mp.mp.dps = 60
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
h = json.load(open(wd+"/SRF3_hi_coeffs.json"))
print("K =", h['K'])
tau = [mp.mpf(x) for x in h['eta']]
K = h['K']
z0 = (-mp.log(2)+mp.pi*1j)/2
R = abs(z0); th = float(mp.arg(z0))
# normalized c_k = (a_k/k!) R^k, in float for LSQ (values are small enough)
c = [float(tau[k]/mp.factorial(k)*R**k) for k in range(K+1)]

def fit_model(ks, ncorr, beta=-1.5):
    M = np.zeros((len(ks), 2*ncorr)); y = np.array([c[k] for k in ks])
    for i,k in enumerate(ks):
        for j in range(ncorr):
            M[i,2*j]   = k**(beta-j)*math.cos(th*k)
            M[i,2*j+1] = k**(beta-j)*math.sin(th*k)
    x, res, *_ = np.linalg.lstsq(M, y, rcond=None)
    r = y - M@x
    C = math.hypot(x[0],x[1]); phi = math.atan2(x[1],x[0])
    return C, phi, float(np.sqrt(np.mean(r*r))), float(np.max(np.abs(r)))

print("=== K=120 fits: C, phi with ncorr corrections ===")
for k0 in [40,60,80,100]:
    for ncorr in [2,4,6]:
        C,phi,rms,mx = fit_model(list(range(k0,K+1)), ncorr)
        print(f"k0={k0} ncorr={ncorr}: C={C:.10f} phi={phi:.10f} rms={rms:.1e}")
    print()

# beta determination with K=120: scan beta, best rms
print("=== beta scan (window 60..120, ncorr=3) ===")
for beta in [-1.7,-1.6,-1.55,-1.52,-1.51,-1.5,-1.49,-1.48,-1.45,-1.4,-1.3]:
    C,phi,rms,mx = fit_model(list(range(60,K+1)), 3, beta=beta)
    print(f"beta={beta}: rms={rms:.3e} C={C:.8f}")
