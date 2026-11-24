# Planner check of D16: P/sqrt2 in tau = 1-1/p (pure weight) vs my own polar perimeter of the squircle.
import mpmath as mp
exec(open("d15_check/check_d15.py").read().split("for m in (2, 3):")[0])
mp.mp.dps = 45
z2, z3, z4, z5, z7, pi = [mp.zeta(k) for k in (2, 3, 4, 5, 7)] + [mp.pi]
c = [4, 0, z2, z3, -mp.mpf(33)/16*z4, -(mp.mpf(3)/4*z5 + mp.mpf(7)/4*z2*z3),
     mp.mpf(577)/34560*pi**6 - mp.mpf(7)/8*z3**2,
     mp.mpf(129)/16*z7 + mp.mpf(47)/32*pi**2*z5 + mp.mpf(119)/640*pi**4*z3,
     -mp.mpf(2891923)/116121600*pi**8 + mp.mpf(169)/192*pi**2*z3**2 + mp.mpf(141)/16*z3*z5]
for p in (mp.mpf('1.04'), mp.mpf('1.02'), mp.mpf('1.01')):
    t = 1 - 1/p
    P = perim(2, p)/mp.sqrt(2)
    for K in (5, 8):
        S = sum(c[k]*t**k for k in range(K+1))
        print("p", p, "K", K, "  (P/sqrt2 - series)/tau^(K+1) =", mp.nstr((P - S)/t**(K+1), 8))
