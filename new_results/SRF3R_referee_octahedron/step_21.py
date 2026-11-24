import json
from mpmath import mp, mpf, mpc, log, pi, sqrt, factorial, re as mpre, im as mpim
mp.dps = 80
z0 = (-log(2)+pi*1j)/2
R = abs(z0); th = mp.arg(z0)

srf3 = json.load(open("../SRF3_octahedron_radius_proof/SRF3_hi_coeffs_K120.json"))
srf1 = json.load(open("../SRF1_radii_3d/SRF1_n3_K60.json"))

# parse
tau3 = [mpf(x) for x in srf3['eta']]
d3   = [mpf(x) for x in srf3['d']]
tau1 = [mpf(x) for x in srf1['tau']]
d1   = [mpf(x) for x in srf1['d']]

# do the files agree?
mx_tau = max(abs(tau3[k]-tau1[k]) for k in range(61))
mx_d   = max(abs(d3[k]-d1[k]) for k in range(61))
print("max |SRF3 tau - SRF1 tau| (k<=60):", mp.nstr(mx_tau,5))
print("max |SRF3 d   - SRF1 d  | (k<=60):", mp.nstr(mx_d,5))

# ratio test for d-series: a_k(d)/k! / (k^-3/2 * Im[z0^{1/2-k}])  with principal branch
print("\n k | ratio_d (z0) | ratio_d (zbar0) | ratio_tau (claim: -sqrt(3pi/8)*Re[z0^{1/2-k}])")
for k in range(52, 61):
    bk = d1[k]/factorial(k)
    base = mpf(k)**(-mpf('1.5'))
    r1 = bk/(base*mpim(z0**(mpf('0.5')-k)))
    r2 = bk/(base*mpim(mp.conj(z0)**(mpf('0.5')-k)))
    ak = tau1[k]/factorial(k)
    r3 = ak/(base*mpre(z0**(mpf('0.5')-k)))
    print(k, mp.nstr(r1,12), "|", mp.nstr(r2,12), "|", mp.nstr(r3,12))
print("\n+sqrt(3pi/4) =", mp.nstr(sqrt(3*pi/4),15), " -sqrt(3pi/4) =", mp.nstr(-sqrt(3*pi/4),15))
print("-sqrt(3pi/8) =", mp.nstr(-sqrt(3*pi/8),15))
