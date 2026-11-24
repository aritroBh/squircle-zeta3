import json, mpmath as mp
mp.mp.dps = 60
d = json.load(open('<project>/research/08_theory/SRF2_cube_residue/EV_coeffs_K40.json'))
V = [mp.mpf(x) for x in d['V']]
R = -1-2*mp.pi**mp.mpf('1.5')/mp.gamma(mp.mpf('0.25'))**2
I2 = mp.quad(lambda a: (1+a)*a**3/(mp.sqrt(1+a**6)+1), [0, mp.mpf('0.3'), 1])
rho2 = 3*(-2+I2)
print("k   (-1)^k V_k - R          (-1)^k delta*2^{k+1}   rho2=-5.380664")
for k in range(3, 41):
    dk = (-1)**k*V[k] - R
    print(f"{k:3d}  {mp.nstr(dk, 14)}   {mp.nstr(dk*2**(k+1), 12)}")
