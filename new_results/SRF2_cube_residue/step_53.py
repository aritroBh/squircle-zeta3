import json, mpmath as mp
mp.mp.dps = 60
d = json.load(open('<project>/research/08_theory/SRF2_cube_residue/EV_coeffs_K40.json'))
V = [mp.mpf(x) for x in d['V']]
R = -1-2*mp.pi**mp.mpf('1.5')/mp.gamma(mp.mpf('0.25'))**2
rho2 = mp.mpf('-5.3806643388446858024')

print("k      delta_k^(2) = delta_k - rho2(-1)^k/2^{k+1}      sign   log10|.|")
for k in range(5, 41):
    dk = V[k] - R*(-1)**k - rho2*(-1)**k/2**(k+1)
    if dk != 0:
        print(f"{k:3d}  {mp.nstr(dk, 14)}   {'+' if dk>0 else '-'}   {mp.nstr(mp.log10(abs(dk)), 6)}")
