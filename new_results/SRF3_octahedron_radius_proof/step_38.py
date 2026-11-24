
import os, json
import mpmath as mp
mp.mp.dps = 50
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
h = json.load(open(wd+"/SRF3_hi_coeffs.json"))
tau = [mp.mpf(x) for x in h['eta']]; dco = [mp.mpf(x) for x in h['d']]
K = h['K']
z0 = (-mp.log(2)+mp.pi*1j)/2
R = abs(z0); th = mp.arg(z0)
C = mp.sqrt(3*mp.pi*R/8); phi = th/2-mp.pi

def a_asymp(j):
    return C*mp.factorial(j)*R**(-j)*mp.mpf(j)**(-1.5)*mp.cos(j*th+phi)

k = 100
# term-by-term: actual tau vs asymp
print(" m : term_actual   term_asymp   (signed, normalized by k!R^-k k^-3/2)")
nrm = mp.factorial(k)*R**(-k)*mp.mpf(k)**(-1.5)
for m in range(0, 10):
    tA = (-1)**m*mp.binomial(k-1,m)*tau[k-m]/nrm
    tS = (-1)**m*mp.binomial(k-1,m)*a_asymp(k-m)/nrm
    print(f"m={m}: {mp.nstr(tA,8)}   {mp.nstr(tS,8)}")
sumA = sum((-1)**m*mp.binomial(k-1,m)*tau[k-m] for m in range(0,30))
print("sum actual (m<30)/dco[k]:", mp.nstr(sumA/dco[k],8))
sumS = sum((-1)**m*mp.binomial(k-1,m)*a_asymp(k-m) for m in range(0,30))
print("sum asymp/dco[k]:", mp.nstr(sumS/dco[k],8))
print("sum actual vs dco residual:", mp.nstr(abs(sumA-dco[k])/abs(dco[k]),3))
# the law:
law = mp.sqrt(2)*C*nrm*mp.cos(k*th+phi+mp.pi/2)
print("law =", mp.nstr(law,10), " dco =", mp.nstr(dco[k],10), " sumS =", mp.nstr(sumS,10))
