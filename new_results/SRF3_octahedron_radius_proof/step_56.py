
import os, json, math
import numpy as np
import mpmath as mp
mp.mp.dps = 60
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
h = json.load(open(wd+"/SRF3_hi_coeffs_K120.json"))
tau = [mp.mpf(x) for x in h['eta']]
K = h['K']
z0 = (-mp.log(2)+mp.pi*1j)/2
R = abs(z0); th = float(mp.arg(z0))
b = [float(tau[k]/mp.factorial(k)*R**k) for k in range(K+1)]

def rms_beta(beta, ks, ncorr=3):
    M = np.zeros((len(ks),2*ncorr)); y = np.array([b[k] for k in ks])
    for i,k in enumerate(ks):
        for j in range(ncorr):
            M[i,2*j]   = k**(beta-j)*math.cos(th*k)
            M[i,2*j+1] = k**(beta-j)*math.sin(th*k)
    x,*_ = np.linalg.lstsq(M,y,rcond=None)
    return float(np.sqrt(np.mean((y-M@x)**2)))

def beta_est(k0):
    ks = list(range(k0,K+1))
    betas = np.arange(-1.505,-1.494,0.0002)
    rms = [rms_beta(bt,ks) for bt in betas]
    i = int(np.argmin(rms))
    # parabolic interpolation
    y1,y2,y3 = rms[i-1],rms[i],rms[i+1]
    x1,x2,x3 = betas[i-1],betas[i],betas[i+1]
    denom = (y1-y3)
    delta = 0.5*(y1-y3)/(y1-2*y2+y3)*(x3-x1)/2
    return betas[i]+delta, rms[i]

print("beta estimates (window, ncorr=3):")
for k0 in [40,50,60,70,80,90]:
    be, r = beta_est(k0)
    print(f"  k0={k0}: beta={be:.5f}  (rms={r:.1e})")
