import json, mpmath as mp
mp.mp.dps = 60
d = json.load(open('<project>/research/08_theory/SRF2_cube_residue/EV_coeffs_K40.json'))
V = [mp.mpf(x) for x in d['V']]
R = -1-2*mp.pi**mp.mpf('1.5')/mp.gamma(mp.mpf('0.25'))**2
rho2 = mp.mpf('-5.3806643388446858024')

# rho_{-3} = binom(5,3) K^{(3)}(-3^+) = 10(-8/3 + I3), I3 = int_0^1 (1+a)^2 a^4/(sqrt(1+a^8)+1) da
I3 = mp.quad(lambda a: (1+a)**2*a**4/(mp.sqrt(1+a**8)+1), [0, mp.mpf('0.3'), 1])
rho3 = 10*(-mp.mpf(8)/3 + I3)
print("I3 =", mp.nstr(I3, 20), "  rho_{-3} =", mp.nstr(rho3, 20))

# fit rho3 from data: delta_k^{(2)} ~ (-1)^k rho3/3^{k+1} + B_k ; use k=8..16, model B_k = b0 + b1*(-1)^k/4^{k}... 
# simple: use pairs (k, k+1) to eliminate smooth B_k: d_k + d_{k+1} ~ (-1)^k rho3(1/3^{k+1} - 1/3^{k+2}) = (-1)^k rho3 * 2/3^{k+2}
for k in range(8, 17):
    dk  = V[k]   - R*(-1)**k   - rho2*(-1)**k/2**(k+1)
    dk1 = V[k+1] - R*(-1)**(k+1) - rho2*(-1)**(k+1)/2**(k+2)
    s = dk + dk1
    est = (-1)**k*s*3**(k+2)/2
    print(f"k={k:2d}: sum={mp.nstr(s,12)}  rho3 est={mp.nstr(est,10)}")
