
import os, json, math
import mpmath as mp
mp.mp.dps = 50
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
h = json.load(open(wd+"/SRF3_hi_coeffs.json"))
tau = [mp.mpf(x) for x in h['eta']]; dco = [mp.mpf(x) for x in h['d']]
K = h['K']
z0 = (-mp.log(2)+mp.pi*1j)/2
R = abs(z0); th = mp.arg(z0)
C_pred = mp.sqrt(3*mp.pi*R/8); phi = th/2-mp.pi
thf = float(th); phif = float(phi)

# normalized d-coefficients
print("k :  dco[k]/(k! R^-k k^-3/2)  vs  sqrt(2)C cos(k th + phi + pi/2)")
for k in range(100, 121, 4):
    norm = dco[k]/(mp.factorial(k)*R**(-k)*mp.mpf(k)**(-1.5))
    pred = mp.sqrt(2)*C_pred*mp.cos(k*th+phi+mp.pi/2)
    pred2 = C_pred*mp.re(mp.e**(1j*(k*thf+phif))*mp.e**(-mp.conj(z0)))
    print(f"k={k}: norm={mp.nstr(norm,8)}  sqrt2pred={mp.nstr(pred,8)}  ez_pred={mp.nstr(pred2,8)}  ratio={mp.nstr(norm/pred,5)}")
# Also: maybe the fitted phi_d was polluted; do a fresh careful LSQ on d over k=100..120
import numpy as np
def mlsq(ks, ncorr, aco):
    y = [aco[k]/mp.factorial(k)*R**k for k in ks]
    nrow=2*ncorr
    S=[[mp.mpf(0)]*nrow for _ in range(nrow)]; T=[mp.mpf(0)]*nrow
    for i,k in enumerate(ks):
        bas=[]
        for j in range(ncorr):
            bas += [k**(-mp.mpf(1.5)-j)*mp.cos(thf*k), k**(-mp.mpf(1.5)-j)*mp.sin(thf*k)]
        for a in range(nrow):
            T[a]+=bas[a]*y[i]
            for b in range(nrow): S[a][b]+=bas[a]*bas[b]
    return mp.lu_solve(mp.matrix(S), mp.matrix(T))
x = mlsq(list(range(100,121)), 3, dco)
Cd = mp.sqrt(x[0]**2+x[1]**2); phid = mp.atan2(x[1],x[0])
print("\nd-series k0=100 nc=3: C_d=", mp.nstr(Cd,14), "phi_d=", mp.nstr(phid,12), "C_d/C=", mp.nstr(Cd/C_pred,8))
resid = [dco[k]/mp.factorial(k)*R**k - (x[0]*k**(-1.5)*mp.cos(thf*k)+x[1]*k**(-1.5)*mp.sin(thf*k)+x[2]*k**(-2.5)*mp.cos(thf*k)+x[3]*k**(-2.5)*mp.sin(thf*k)+x[4]*k**(-3.5)*mp.cos(thf*k)+x[5]*k**(-3.5)*mp.sin(thf*k)) for k in range(100,121)]
print("max resid:", mp.nstr(max(abs(r) for r in resid),4))
