
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
thf=float(th); phif=float(phi)

# asymptotic model of a_j (tau) used inside the exact mixing
def a_asymp(j):
    return C*mp.factorial(j)*R**(-j)*mp.mpf(j)**(-1.5)*mp.cos(j*th+phi)

print("ratio actual d-coef / [sum_m (-1)^m C(k-1,m) a_asymp(k-m)]  and / (sqrt2 law):")
for k in [60,80,100,120]:
    s = mp.mpf(0)
    for m in range(0, k-5):  # j=k-m >= 5 (asymp ok there); remaining j tiny
        s += (-1)**m*mp.binomial(k-1,m)*a_asymp(k-m)
    law = mp.sqrt(2)*C*mp.factorial(k)*R**(-k)*mp.mpf(k)**(-1.5)*mp.cos(k*th+phi+mp.pi/2)
    lawm = -law
    print(f"k={k}: actual/sum={mp.nstr(dco[k]/s,7)}  actual/law={mp.nstr(dco[k]/law,7)}")
# check convergence of the sum to the sqrt2 law
for k in [60,120]:
    s = mp.mpf(0)
    for m in range(0, k-5):
        s += (-1)**m*mp.binomial(k-1,m)*a_asymp(k-m)
    law = mp.sqrt(2)*C*mp.factorial(k)*R**(-k)*mp.mpf(k)**(-1.5)*mp.cos(k*th+phi+mp.pi/2)
    print(f"k={k}: sum/law={mp.nstr(s/law,7)}")
