"""Numerical confirmation.
Part A (k<=400): exact M_k (80-digit balls, D13f) minus M_k^main (quadrature) = rest; compare with S_k and the
  leading law L2 = (sqrt(pi)/N) e^{-th2} cos(th2 - pi/4), and with the varpi2-corrected law.
Part B (large k): S_k by quadrature vs Zpred; checks rho_k -> 1 and (rho_k - 1 - varpi2/sqrt N)*N bounded."""
import mpmath as mp
from common import *
mp.mp.dps = 60
sp_ = mp.sqrt(pi)
varpi2 = mp.sqrt(2)*((12+33*pi+2*pi**2) + 1j*(12-33*pi+2*pi**2))/(48*sp_)
print("sqrt(pi) =", mp.nstr(sp_, 40)); print("varpi2   =", mp.nstr(varpi2, 30))
M = load_M()
print("\nPart A.  nr = N e^{th2} rest / sqrt(pi);  U = rest - S_k;  u = |U| N e^{sqrt(3 pi N)}/sqrt(pi)")
print(f"{'k':>4} {'rest':>22} {'rest/S_k-1':>12} {'nr':>10} {'cos(th2-pi/4)':>13} {'nr(varpi2)':>10} {'u':>8}")
for k in range(100, 401, 25):
    N = k+1; th = mp.sqrt(2*pi*N)
    rest = M[k]-Mmain(k); s = S2(k); U = rest-s
    nr = N*mp.exp(th)*rest/sp_
    pred2 = mp.re(mp.expj(th-pi/4)*(1+varpi2/mp.sqrt(N)))
    u = abs(U)*N*mp.exp(mp.sqrt(3*pi*N))/sp_
    print(f"{k:4d} {mp.nstr(rest,15):>22} {mp.nstr(rest/s-1,4):>12} {mp.nstr(nr,7):>10} {mp.nstr(mp.cos(th-pi/4),7):>13} {mp.nstr(pred2,7):>10} {mp.nstr(u,4):>8}")
mp.mp.dps = 40
print("\nPart B.  rho_k = Z_k / Zpred,  Zpred = -(pi^{3/2}/(2N)) e^{-th2} e^{i(th2-3pi/4)}")
print(f"{'k':>8} {'rho_k':>40} {'(rho-1)sqrtN':>34} {'(rho-1-varpi2/sqrtN) N':>34}")
for k in [100, 400, 1600, 6400, 25600, 102400, 409600, 1638400, 6553600]:
    N = k+1; th = mp.sqrt(2*pi*N)
    zp = -(pi**1.5/(2*N))*mp.exp(-th)*mp.expj(th-3*pi/4)
    r = Z2(k)/zp
    print(f"{k:8d} {mp.nstr(r,14):>40} {mp.nstr((r-1)*mp.sqrt(N),10):>34} {mp.nstr((r-1-varpi2/mp.sqrt(N))*N,8):>34}")
