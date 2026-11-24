
import os, json, math
import numpy as np
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
base = "<project>"
h1 = json.load(open(os.path.join(base,"research/10_frontier/H1_octahedron_end/hi_coeffs_n3_K80.json")))
tau = [mp.mpf(x) for x in h1['eta']] if False else None
import mpmath as mp
mp.mp.dps = 80
tau = [mp.mpf(x) for x in h1['eta']]
K = 80
z0 = complex(-math.log(2)/2, math.pi/2)
R = abs(z0); th = math.atan2(z0.imag, z0.real)
b = [float(tau[k]/mp.factorial(k)*R**k) for k in range(K+1)]  # c_k

def fit_model(ks, ncorr):
    # c_k = sum_{j=0}^{ncorr-1} k^{-3/2-j} (A_j cos(th k) + B_j sin(th k))
    M = np.zeros((len(ks), 2*ncorr)); y = np.array([b[k] for k in ks])
    for i,k in enumerate(ks):
        for j in range(ncorr):
            M[i,2*j]   = k**(-1.5-j)*math.cos(th*k)
            M[i,2*j+1] = k**(-1.5-j)*math.sin(th*k)
    x, res, *_ = np.linalg.lstsq(M, y, rcond=None)
    r = y - M@x
    C = math.hypot(x[0],x[1]); phi = math.atan2(x[1],x[0])
    return C, phi, float(np.sqrt(np.mean(r*r))), float(np.max(np.abs(r)))

print("C from fits with ncorr correction terms, varying window:")
for k0 in [20,30,40,50,60]:
    row=[]
    for ncorr in [1,2,3,4]:
        C,phi,rms,mx = fit_model(list(range(k0,K+1)), ncorr)
        row.append((C,phi,rms))
    print(f"k0={k0}: " + "  ".join(f"C={C:.6f}(rms{rms:.1e})" for C,phi,rms in row))

# predicted phase phi from theory: c_k ~ Re[D' e^{-ik th}] form; print fitted phi
C,phi,rms,mx = fit_model(list(range(40,K+1)), 3)
print("\nbest fit k0=40 ncorr=3: C=%.8f phi=%.8f rms=%.2e"%(C,phi,rms))
# theoretical candidates
pred1 = 2*math.sqrt(6*R/math.pi)         # 2*sqrt(6R/pi) = 3.5056 (earlier pred)
pred2 = pred1*abs(complex(-0.5*math.log(2),math.pi/2))**0/1
print("pred |D_full|/sqrt(pi) style C_pred = 2*sqrt(6R/pi) =", pred1)
print("C/pred =", C/pred1)
