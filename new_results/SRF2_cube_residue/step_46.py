import sys, json
sys.path.insert(0, '<project>/research/08_theory/SRF2_cube_residue')
import mpmath as mp
mp.mp.dps = 60
d = json.load(open("<project>/research/08_theory/SRF1_radii_3d/EV_coeffs_K30.json"))
V = [mp.mpf(x) for x in d["V"]]
R = -1-2*mp.pi**mp.mpf('1.5')/mp.gamma(mp.mpf('0.25'))**2

# I2' = int_0^1 a^{-3}[(1+a)sqrt(1+a^6)-1-a] da  (stable form)
I2 = mp.quad(lambda a: (1+a)*a**3/(mp.sqrt(1+a**6)+1), [0, mp.mpf('0.3'), 1])
rho_m2 = 3*(-2 + I2)
print("I2' =", mp.nstr(I2, 25))
print("rho_{-2} = 3(-2+I2') =", mp.nstr(rho_m2, 25))

# delta_k from data, fit rho_{-2}: delta_k ~ rho_{-2} (-1)^k / 2^{k+1}
print("\n k      (-1)^k delta_k      rho fit (k)     ratio")
prev = None
for k in range(3, 31):
    dk = V[k] - R*(-1)**k
    sk = (-1)**k*dk
    rho_fit = sk*2**(k+1)
    print(f"{k:3d}  {mp.nstr(sk, 12)}   {mp.nstr(rho_fit, 10)}" + (f"   {mp.nstr(prev/sk,6)}" if prev else ""))
    prev = sk
