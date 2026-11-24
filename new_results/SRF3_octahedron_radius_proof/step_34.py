
import os, json, math
import numpy as np
import mpmath as mp
mp.mp.dps = 50
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
h = json.load(open(wd+"/SRF3_hi_coeffs.json"))
tau = [mp.mpf(x) for x in h['eta']]; dco = [mp.mpf(x) for x in h['d']]
K = h['K']
z0 = (-mp.log(2)+mp.pi*1j)/2
R = abs(z0); th = float(mp.arg(z0))

# --- joint beta fit: c_k = k^beta [A0 cos + B0 sin + k^-1(A1 cos+B1 sin) + k^-2(...)]
def fit_beta(ks, ncorr):
    best=None
    for beta in np.arange(-1.53,-1.46,0.0005):
        M = np.zeros((len(ks), 2*ncorr)); y = np.array([float(tau[k]/mp.factorial(k)*R**k) for k in ks])
        for i,k in enumerate(ks):
            for j in range(ncorr):
                M[i,2*j]   = k**(beta-j)*math.cos(th*k)
                M[i,2*j+1] = k**(beta-j)*math.sin(th*k)
        x, *_ = np.linalg.lstsq(M, y, rcond=None)
        rms = float(np.sqrt(np.mean((y-M@x)**2)))
        if best is None or rms<best[0]: best=(rms,beta,x)
    return best

for k0 in [50,70,90]:
    rms,beta,x = fit_beta(list(range(k0,K+1)),3)
    print(f"k0={k0}: beta={beta:.4f} rms={rms:.2e}")

# --- mpmath high-precision LSQ for C on tau series (beta,theta fixed exact)
def mlsq(ks, ncorr, aco):
    y = [aco[k]/mp.factorial(k)*R**k for k in ks]
    nrow=2*ncorr
    S=[[mp.mpf(0)]*nrow for _ in range(nrow)]; T=[mp.mpf(0)]*nrow
    for i,k in enumerate(ks):
        bas=[]
        for j in range(ncorr):
            bas += [k**(-mp.mpf(1.5)-j)*mp.cos(th*k), k**(-mp.mpf(1.5)-j)*mp.sin(th*k)]
        for a in range(nrow):
            T[a]+=bas[a]*y[i]
            for b in range(nrow): S[a][b]+=bas[a]*bas[b]
    x = mp.lu_solve(mp.matrix(S), mp.matrix(T))
    return x

C_pred = mp.sqrt(3*mp.pi*R/8)
phi_pred = th/2-mp.pi
print("\nC_pred =", mp.nstr(C_pred,20))
print("phi_pred = theta/2 - pi =", mp.nstr(phi_pred,15))
for nm, aco in [("tau",tau),("d",dco)]:
    for k0,nc in [(60,4),(80,5),(90,5)]:
        x = mlsq(list(range(k0,K+1)), nc, aco)
        C = mp.sqrt(x[0]**2+x[1]**2); phi = mp.atan2(x[1],x[0])
        print(f"{nm}: k0={k0} nc={nc}: C={mp.nstr(C,16)} phi={mp.nstr(phi,14)} C/Cpred={mp.nstr(C/C_pred,8)}")
