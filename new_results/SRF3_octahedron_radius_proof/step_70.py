
import os, json
import mpmath as mp
mp.mp.dps = 60
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
h = json.load(open(wd+"/SRF3_hi_coeffs_K120.json"))
dco = [mp.mpf(x) for x in h['d']]
K = h['K']
z0 = (-mp.log(2)+mp.pi*1j)/2
C2 = mp.sqrt(3*mp.pi/4)
print("target: ratios oscillate about -sqrt(3 pi/4) =", mp.nstr(-C2,10))
print("k : b_k(d) / (k^-3/2 Im[z0^{1/2-k}])")
rats = []
for k in list(range(52,61))+[80,100,110,120]:
    r = dco[k]/(mp.factorial(k)*mp.mpf(k)**(-1.5)*mp.im(z0**(mp.mpf('0.5')-k)))
    rats.append(float(r))
    print(f"  k={k}: {mp.nstr(r,10)}")
import statistics
print("mean of k=52..60:", round(statistics.mean(rats[:9]),4), " (theory -1.53499)")
# and the z0bar form (should be +):
for k in [56,58,60]:
    r = dco[k]/(mp.factorial(k)*mp.mpf(k)**(-1.5)*mp.im(mp.conj(z0)**(mp.mpf('0.5')-k)))
    print(f"  z0bar form k={k}: {mp.nstr(r,10)} (theory +1.53499)")
