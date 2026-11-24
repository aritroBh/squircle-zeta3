
import os, json
import mpmath as mp
mp.mp.dps = 60
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
h = json.load(open(wd+"/SRF3_hi_coeffs.json"))
tau = [mp.mpf(x) for x in h['eta']]; dco = [mp.mpf(x) for x in h['d']]
K = h['K']
z0 = (-mp.log(2)+mp.pi*1j)/2
R = abs(z0); th = float(mp.arg(z0))

# fit full tau model: a_k = k! R^-k k^-3/2 sum_{r=0}^{nc-1} k^-r (A_r cos + B_r sin)
nc = 6
ks = list(range(60, K+1))
nrow = 2*nc
S=[[mp.mpf(0)]*nrow for _ in range(nrow)]; T=[mp.mpf(0)]*nrow
for k in ks:
    y = tau[k]/mp.factorial(k)*R**k
    bas=[]
    for r in range(nc):
        bas += [k**(-1.5-r)*mp.cos(th*k), k**(-1.5-r)*mp.sin(th*k)]
    for a in range(nrow):
        T[a]+=bas[a]*y
        for b in range(nrow): S[a][b]+=bas[a]*bas[b]
x = mp.lu_solve(mp.matrix(S), mp.matrix(T))
C_tau = mp.sqrt(x[0]**2+x[1]**2); phi_tau = mp.atan2(-x[1],x[0])  # model: x0 cos + x1 sin = C cos(kth - ph)... careful
# c_k = x0 cos(kth)+x1 sin(kth) = C cos(kth - phi) with phi=atan2(x1,x0)
phi_std = mp.atan2(x[1],x[0])
print("tau model: C=",mp.nstr(C_tau,14)," phi(atan2)=",mp.nstr(phi_std,14), " resid max:",
      mp.nstr(max(abs(tau[k]/mp.factorial(k)*R**k - sum((x[2*r]*mp.cos(th*k)+x[2*r+1]*mp.sin(th*k))*k**(-1.5-r) for r in range(nc))) for k in ks),3))

def model_tau(j):
    return mp.factorial(j)*R**(-j)*sum((x[2*r]*mp.cos(th*j)+x[2*r+1]*mp.sin(th*j))*mp.mpf(j)**(-1.5-r) for r in range(nc))

# exact mixing with the model -> model_d, compare with actual d coefficients
errs=[]
for k in range(70, K+1):
    md = sum((-1)**m*mp.binomial(k-1,m)*model_tau(k-m) for m in range(0,k-40))
    errs.append(abs(md-dco[k])/abs(dco[k]))
print("model-mixing vs actual d: worst rel err (k=70..120):", mp.nstr(max(errs),3))

# now extract C_d, phi_d from the model_d over k=100..120 via projection
ksd = list(range(100,121))
y = []
for k in ksd:
    md = sum((-1)**m*mp.binomial(k-1,m)*model_tau(k-m) for m in range(0,k-40))
    y.append(md/mp.factorial(k)*R**k*mp.mpf(k)**1.5)
X0 = sum(yi*mp.cos(th*k) for yi,k in zip(y,ksd))/len(ksd)*2
X1 = sum(yi*mp.sin(th*k) for yi,k in zip(y,ksd))/len(ksd)*2
Cd = mp.sqrt(X0**2+X1**2); phid = mp.atan2(X1,X0)
print("C_d from model:", mp.nstr(Cd,14), "phi_d =", mp.nstr(phid,14))
print("sqrt2*C_tau =", mp.nstr(mp.sqrt(2)*C_tau,14))
print("phi+pi/2 =", mp.nstr(phi_std+mp.pi/2,14), " phi+3pi/2-2pi =", mp.nstr(phi_std+3*mp.pi/2-2*mp.pi,14))
