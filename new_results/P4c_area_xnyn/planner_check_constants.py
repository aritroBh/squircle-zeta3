# Planner check of P4c: planner's own independent 25-digit n*E(n) values (from the P4b brief) vs A(n) = 1 + sqrt2/n + R/n^2 + R3/n^3 + R4/n^4
import mpmath as mp
mp.mp.dps = 30
lam, z2, z3 = mp.log(2), mp.pi**2/6, mp.zeta(3)
beta = lam/2 - 1
gam = -(1+mp.sqrt(2)) + lam/2 - mp.pi**2/12 - lam**2/8
b3 = (1-lam)/2; c3 = lam**2/4 - lam/2 + 1 + z2
d3 = -lam**3/24 + lam**2/8 - lam/2 - z2*lam/2 + z2/2 - z3/2 + 1 + 3*mp.sqrt(2)/2
nE = {10**3: mp.mpf('1.3830271752121981069910'), 10**4: mp.mpf('1.4090781805833874542671'), 10**5: mp.mpf('1.4134461509073296926161')}
print("beta", mp.nstr(beta, 20), "gamma", mp.nstr(gam, 20), "d3", mp.nstr(d3, 20))
for n, v in nE.items():
    L = mp.log(n)
    R = -L**2/2 + beta*L + gam
    R3 = L**3/3 + b3*L**2 + c3*L + d3
    data = n*(v - mp.sqrt(2))            # = R + R3/n + R4/n^2 + ...
    r1 = data - R; r2 = n*(data - R) - R3
    print(n, " n(nE-sqrt2)-R =", mp.nstr(r1, 8), "  (R3/n =", mp.nstr(R3/n, 8), ")   next level n*r1 - R3 =", mp.nstr(r2, 6), " vs -(1/8)log^4 n / n =", mp.nstr(-L**4/8/n, 6))
