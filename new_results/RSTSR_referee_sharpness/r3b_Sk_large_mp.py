"""Part B (non-rigorous, independent of RSTS code): S_k for k up to 6.5e6 by mpmath quadrature in log form.
Z_k = int_0^inf H(s) e^{i pi s/2} A2_k(s) ds on s = r w, w = e^{i pi/4};
log integrand = loggamma(s/2) - loggamma(s/2+3/2) + log(sqrt(pi)/4) + i pi s/2 + log(-2(s+1)) + (k-1)log s - (k+1)log(s+2).
Breakpoints: saddle r* = 2 lam, width wid = 2 lam / sqrt(Lam); r in [r*/8, r* + 40 wid] (outside: < e^{-40}-relative,
by Re psi(r) - Re psi(1) = -(r-1)^2/(sqrt2 r) and region (c) of RSTS)."""
import mpmath as mp
mp.mp.dps = 40
pi = mp.pi; w = mp.expjpi(mp.mpf(1)/4)
v2 = mp.sqrt(2)*((12+33*pi+2*pi**2) + 1j*(12-33*pi+2*pi**2))/(48*mp.sqrt(pi))
def Z2(k):
    N = k+1; lam = mp.sqrt(N/pi); Lam = mp.sqrt(pi*N); rs = 2*lam; wid = rs/mp.sqrt(Lam)
    th = mp.sqrt(2*pi*N)
    def g(r):
        s = r*w
        L = (mp.loggamma(s/2) - mp.loggamma(s/2 + mp.mpf(3)/2) + mp.log(mp.sqrt(pi)/4) + 1j*pi*s/2
             + mp.log(-2*(s+1)) + (k-1)*mp.log(s) - (k+1)*mp.log(s+2))
        return mp.exp(L + th)*w                     # scaled by e^{th}
    lo = max(rs/8, rs - 40*wid); hi = rs + 60*wid
    pts = [lo + (hi-lo)*j/80 for j in range(81)]
    return mp.quad(g, pts)*mp.exp(-th)
print("k | (rho-1)sqrtN | (rho-1-v2/sqrtN) N | N e^th S/sqrt(pi) - cos(th-pi/4) times sqrtN")
rows = []
for k in [100, 400, 1600, 6400, 25600, 102400, 409600, 1638400, 6553600]:
    N = k+1; th = mp.sqrt(2*pi*N)
    zp = -(pi**1.5/(2*N))*mp.exp(-th)*mp.expj(th-3*pi/4)
    Z = Z2(k); rho = Z/zp; S = 2/pi*mp.im(Z)
    rows.append((N, rho))
    print(k, mp.nstr((rho-1)*mp.sqrt(N), 12), mp.nstr((rho-1-v2/mp.sqrt(N))*N, 8),
          mp.nstr((N*mp.exp(th)*S/mp.sqrt(pi) - mp.cos(th-pi/4))*mp.sqrt(N), 8), flush=True)
pts = rows[-4:]
A = mp.matrix([[1, 1/mp.sqrt(n), mp.mpf(1)/n, 1/mp.mpf(n)**1.5] for n, _ in pts])
for part in (mp.re, mp.im):
    sol = mp.lu_solve(A, mp.matrix([part((r-1)*mp.sqrt(n)) for n, r in pts]))
    print("Richardson varpi2 estimate:", mp.nstr(sol[0], 14), " next coeff:", mp.nstr(sol[1], 8))
print("RSTS varpi2:", mp.nstr(v2, 16))
