
import os, json
import mpmath as mp
mp.mp.dps = 60
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
h = json.load(open(wd+"/SRF3_hi_coeffs_K120.json"))
tau = [mp.mpf(x) for x in h['eta']]; dco = [mp.mpf(x) for x in h['d']]
K = h['K']
z0 = (-mp.log(2)+mp.pi*1j)/2
R = abs(z0); th = mp.arg(z0); thf=float(th)
C1 = mp.sqrt(3*mp.pi/8); C2 = mp.sqrt(3*mp.pi/4)
def law_tau(k): return -C1*mp.factorial(k)*mp.mpf(k)**(-mp.mpf('1.5'))*mp.re(z0**(mp.mpf('0.5')-k))
def law_d(k):   return +C2*mp.factorial(k)*mp.mpf(k)**(-mp.mpf('1.5'))*mp.im(z0**(mp.mpf('0.5')-k))
# fit correction coefficients (1/k term) to show convergence: ratio / (1 + e1/k)
print("k | tau law ratio | tau ratio/(1+e1/k) | d law ratio | d ratio/(1+e1d/k)")
# estimate e1 from k=100,110
for kfix in [100]:
    pass
rows=[]
for k in [60,80,100,110,120]:
    rt = tau[k]/law_tau(k); rd = dco[k]/law_d(k)
    rows.append((k, rt, rd))
# estimate e1 from last two
k1,k2 = 110,120
e1t = (rows[4][1]-rows[3][1])/(1/k2-1/k1)*0  # just present ratios with a Richardson extrapolated value
# Richardson extrapolation of the ratios to k=inf (model r_k = 1 + a/k + b/k^2)
import numpy as np
def richardson(vals, ks):
    A = np.array([[1,1/k,1/k**2] for k in ks], float); b = np.array(vals,float)
    x,*_ = np.linalg.lstsq(A,b,rcond=None); return x[0]
rt_ex = richardson([float(r[1]) for r in rows],[r[0] for r in rows])
rd_ex = richardson([float(r[2]) for r in rows],[r[0] for r in rows])
for k, rt, rd in rows:
    print(f"k={k}: {mp.nstr(rt,12)}  {mp.nstr(rd,12)}")
print("Richardson extrapolated ratios to k=inf: tau:", round(rt_ex,8), " d:", round(rd_ex,8))
print("\nC_tau exact: 1.3766110843131085158 ; C_d = sqrt(3 pi R/4):", mp.nstr(mp.sqrt(3*mp.pi*R/4),18))
print("A = 1/R =", mp.nstr(1/R,18), " theta =", mp.nstr(th,16))
