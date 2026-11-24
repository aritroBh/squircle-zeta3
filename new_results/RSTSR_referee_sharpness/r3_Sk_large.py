"""Part B: S_k (q=2 term) by own Arb integral for k up to 10^6 (rigorous enclosures).
rho_k = Z_k / Zpred, Zpred = -(pi^{3/2}/(2N)) e^{-th} e^{i(th - 3pi/4)},  th = sqrt(2 pi N).
Checks (i): rho_k -> 1; N e^{th} S_k/sqrt(pi) - cos(th - pi/4) = O(N^{-1/2});
and the varpi2 claim: (rho_k - 1) sqrt N -> varpi2 (Richardson fit in powers of N^{-1/2})."""
import sys
from rcommon import *
setup(256)
pi, om, a = consts()
v2 = acb(arb(2).sqrt()*(12+33*pi+2*pi**2), arb(2).sqrt()*(12-33*pi+2*pi**2))/(48*pi.sqrt())
print("varpi2 (RSTS closed form) =", v2.str(20))
ks = [100, 400, 1000, 2500, 10**4, 4*10**4, 16*10**4, 64*10**4, 256*10**4]
rows = []
print("k | S_k | N e^th S_k/sqrt(pi) - cos(th-pi/4) (times sqrtN) | rho_k | (rho-1)sqrtN | (rho-1-v2/sqrtN)N")
for k in ks:
    N = k+1; th = (2*pi*N).sqrt()
    zp = -(pi**1.5/(2*N))*(-th).exp()*(acb(0, 1)*(th - 3*pi/4)).exp()
    lam = float((N/pi).sqrt()); wid = 2*lam/float((pi*N).sqrt())**0.5   # saddle width in r
    Rk = int((2*float(th) + 120)/float(a)) + 10
    Z = Z2(k, atol=abs(zp)*arb("1e-30"), nsplit=max(32, int(2*(Rk-24)/wid)))
    S = 2/pi*Z.imag
    dev = (N*th.exp()*S/pi.sqrt() - (th-pi/4).cos())*arb(N).sqrt()
    rho = Z/zp
    c1 = (rho-1)*arb(N).sqrt(); c2 = (rho-1-v2/arb(N).sqrt())*N
    rows.append((N, rho))
    print(k, S.str(12), dev.str(8), rho.str(14), c1.str(12), c2.str(8), flush=True)
# Richardson: (rho-1) sqrtN = v + c/sqrtN + d/N + e/N^{3/2} fitted on the last 4 points (exact solve)
import mpmath as mp
mp.mp.dps = 40
pts = rows[-4:]
A = mp.matrix([[1, 1/mp.sqrt(n), mp.mpf(1)/n, 1/mp.mpf(n)**1.5] for n, _ in pts])
for part in ("real", "imag"):
    b = mp.matrix([mp.mpf(getattr((r-1)*arb(n).sqrt(), part).mid().str(30)) for n, r in pts])
    sol = mp.lu_solve(A, b)
    print(f"Richardson {part}: varpi2 estimate {mp.nstr(sol[0], 14)}   next coeff {mp.nstr(sol[1], 8)}")
print("RSTS varpi2:", v2.real.str(14), v2.imag.str(14))
