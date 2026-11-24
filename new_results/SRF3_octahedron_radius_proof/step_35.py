
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

# verify exact mixing c_k(d) = sum_j a_j (-1)^{k-j} C(k-1,j-1)
for k in [10, 60, 100]:
    ck = sum(tau[j]*(-1)**(k-j)*mp.binomial(k-1,j-1) for j in range(1,k+1))
    print(f"k={k}: mixing vs actual d-coef: {mp.nstr(abs(ck-dco[k])/abs(dco[k]),3)}")

# asymptotic prediction: c_k(d) ~ C_tau k! R^-k k^-3/2 Re[e^{i(k th+phi)}(1-z0bar)]
zb = mp.conj(z0)
for k in [80,100,110]:
    pred = C_pred*mp.factorial(k)*R**(-k)*k**(-1.5)*mp.re(mp.e**(1j*(k*float(th)+float(phi)))*(1-zb))
    print(f"k={k}: pred={mp.nstr(pred,10)} actual={mp.nstr(dco[k],10)} ratio={mp.nstr(dco[k]/pred,8)}")
# fitted said C_d = sqrt(2)*C_tau, phi_d = phi+3pi/2: test
for k in [80,100,110]:
    pred2 = mp.sqrt(2)*C_pred*mp.factorial(k)*R**(-k)*k**(-1.5)*mp.cos(k*th+phi+3*mp.pi/2)
    print(f"k={k}: sqrt2-pred={mp.nstr(pred2,10)} ratio={mp.nstr(dco[k]/pred2,8)}")
