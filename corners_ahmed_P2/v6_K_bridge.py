"""V6. Bridge between Ahmed's power-corner constant C(theta) and the squircle corner constant K(alpha).

Both are moments of the same density h_c (c = cos theta, alpha = pi - theta = interior angle):
  C(theta)      = int_0^oo [h_c(s) - (1-c) 1_{s>1}] ds/s                 (regularised log-moment)
  -2 K(pi-th)   = int_0^oo h_c(s) ds/(s(1+s))                            (logistic moment)
with K(alpha) = -(1+k)log(1+k) - (1-k)log(1-k), k = cos(alpha/2) = sin(theta/2).

Mellin-Barnes (Parseval with the Mellin transform of h_c, V5a):
  -2K = (1/2 pi i) int_{(sigma)} pi^2 Lambda_w / sin^2(pi w) dw,  -1<sigma<0,
so that, shifting the contour past the double poles w = 0..M:
  -2K + sum_{m=0}^{M} Lambda'_m = I_M := (1/2 pi i) int_{(M+1/2)} pi^2 Lambda_w/sin^2(pi w) dw   (exact for each M)
and, letting M -> oo,   2 K(pi - theta) = sum_{m>=0} Lambda'_m(theta),   while   C(theta) = -Lambda'_0(theta).
Lambda'_m = d/dnu Lambda_nu at nu = m is computed by an exact series (no numerical differentiation).
"""
import mpmath as mp
from common import *

setdps(50)


def Lam(nu, th):
    z = mp.sin(th/2)**2
    return 2*z*mp.hyp2f1(-nu, 1 + nu, 2, z)


def dLam_int(m, z, tol):
    """d/dnu [2z 2F1(-nu,1+nu;2;z)] at nu = m (integer >= 0), exact series."""
    m = int(m)
    s1 = mp.mpf(0)
    for k in range(1, m + 1):
        pk = (-1)**k*mp.factorial(m + k)/mp.factorial(m - k)    # (-m)_k (m+1)_k
        H = mp.fsum(mp.mpf(1)/(m - j) + mp.mpf(1)/(m + 1 + j) for j in range(k))
        s1 += pk*H*z**k/(mp.factorial(k + 1)*mp.factorial(k))
    s2 = mp.mpf(0)
    k = m + 1
    t = mp.factorial(m + k)*z**k/(mp.factorial(k + 1)*mp.factorial(k))   # (k-1-m)! = 0! = 1
    while True:
        s2 += t
        # ratio t_{k+1}/t_k = (k-m)(m+k+1) z / ((k+2)(k+1))
        t = t*(k - m)*(m + k + 1)*z/((k + 2)*(k + 1))
        k += 1
        if abs(t) < tol*abs(s2) and k > m + 5:
            break
    return 2*z*(s1 - (-1)**m*s2)


def dLam_fast(m, c, P, H):
    """Lambda'_m via Lambda_nu = (P_{nu-1} - P_{nu+1})/(2nu+1) and the degree derivative at integers
       d/dnu P_nu |_{nu=n} = P_n(c)[log((1+c)/2) + 2(H_{2n} - H_n)] + 2 sum_{k<n} (-1)^{n+k}(2k+1)P_k(c)/((n-k)(n+k+1))
       (Szmytkowski 2006; re-verified against numerical differentiation in V6e).  P, H: precomputed tables."""
    def dP(n):
        if n < 0:            # P_mu = P_{-mu-1}  =>  dP at mu = n equals -dP at -n-1
            return -dP(-n - 1)
        s = P[n]*(mp.log((1 + c)/2) + 2*(H[2*n] - H[n]))
        s += 2*mp.fsum((-1)**(n + k)*(2*k + 1)*P[k]/((n - k)*(n + k + 1)) for k in range(n))
        return s
    Lm = (P[m - 1] - P[m + 1])/(2*m + 1) if m >= 1 else 1 - c
    return (dP(m - 1) - dP(m + 1))/(2*m + 1) - 2*Lm/(2*m + 1)


def tables(c, N):
    P = [mp.mpf(1), c]
    for k in range(1, N + 1):
        P.append(((2*k + 1)*c*P[k] - k*P[k - 1])/(k + 1))
    H = [mp.mpf(0)]
    for k in range(1, 2*N + 4):
        H.append(H[-1] + mp.mpf(1)/k)
    return P, H


def Kalpha(th):
    k = mp.sin(th/2)
    return -(1 + k)*mp.log(1 + k) - (1 - k)*mp.log(1 - k)


angles = [('pi/3', mp.pi/3, 6), ('pi/2', mp.pi/2, 4), ('2pi/3', 2*mp.pi/3, 3), ('3pi/4', 3*mp.pi/4, 8)]
tol = mp.mpf(10)**-55

print("V6a: the two corner constants as moments of the same density h_c (quadrature vs closed forms)")
for name, th, q in angles:
    c = mp.cos(th)
    z = mp.sin(th/2)**2
    C = 2*z + 2*(1 - z)*mp.log(1 - z)
    Creg = mp.quad(lambda s: hc(s, c)/s, [0, 1]) + mp.quad(lambda s: (hc(s, c) - (1 - c))/s, [1, 10, mp.inf])
    Km = -mp.quad(lambda s: hc(s, c)/(s*(1 + s)), [0, 1, 10, mp.inf])/2
    print(f"  theta={name:>6}: C={mp.nstr(C,40)} |C-reg.int|={mp.nstr(abs(C-Creg),3)} ;"
          f" K(pi-theta)={mp.nstr(Kalpha(th),40)} |K-int|={mp.nstr(abs(Kalpha(th)-Km),3)}")

print("\nV6b: C(theta) = -Lambda'_0(theta)  (exact series vs closed form, and vs mp.diff)")
for name, th, q in angles:
    z = mp.sin(th/2)**2
    C = 2*z + 2*(1 - z)*mp.log(1 - z)
    d0 = dLam_int(0, z, tol)
    dd = mp.diff(lambda nu: Lam(nu, th), 0)
    print(f"  theta={name:>6}: -Lambda'_0={mp.nstr(-d0,40)}  |+C|={mp.nstr(abs(d0+C),3)}  |series - mp.diff|={mp.nstr(abs(d0-dd),3)}")

print("\nV6e: fast degree-derivative formula vs exact hypergeometric series for Lambda'_m (m <= 15)")
for name, th, q in angles:
    z = mp.sin(th/2)**2
    c = mp.cos(th)
    P, H = tables(c, 40)
    worst = max(abs(dLam_fast(m, c, P, H) - dLam_int(m, z, tol)) for m in range(0, 16))
    print(f"  theta={name:>6}: max_m<=15 |fast - exact series| = {mp.nstr(worst,3)}")

print("\nV6d: infinite residue sum  sum_{m>=0} Lambda'_m = 2K(pi-theta).")
print("     Partial sums S_M along complete periods M = qj-1 obey S_M - S ~ M^(-3/2) sum_i c_i M^(-i);")
print("     S is extrapolated by solving for (S, c_0..c_{p-1}) from p+1 values of M <= 6000 (dps 60).")
mp.mp.dps = 60
for name, th, q in angles:
    c = mp.cos(th)
    Mmax = 6000
    P, H = tables(c, Mmax + 2)
    S = []
    acc = mp.mpf(0)
    for m in range(Mmax + 1):
        acc += dLam_fast(m, c, P, H)
        S.append(acc)
    K2 = 2*Kalpha(th)
    for M in [100, 1000, 6000]:
        print(f"  theta={name:>6}: partial sum to M={M}: error {mp.nstr(S[M]-K2,3)}   error*M^1.5 = {mp.nstr((S[M]-K2)*mp.mpf(M)**1.5,4)}")
    for pp in [6, 9, 11]:
        Ms = [q*j - 1 for j in range(Mmax//q - pp*10, Mmax//q + 1, 10)][-(pp + 1):]
        mp.mp.dps = 90      # linear solve in higher precision (basis scaled by Mmax for conditioning)
        A = mp.matrix([[1] + [(mp.mpf(Mmax)/M)**(mp.mpf(3)/2 + i) for i in range(pp)] for M in Ms])
        x = mp.qr_solve(A, mp.matrix([S[M] for M in Ms]))[0]
        mp.mp.dps = 60
        print(f"  theta={name:>6}: extrapolated sum (p={pp:2d}) = {mp.nstr(x[0],45)}   |sum - 2K| = {mp.nstr(abs(x[0]-K2),3)}")
    print(f"  theta={name:>6}: 2K(pi-theta) = {mp.nstr(K2,45)};  C + 2K = sum_(m>=1) Lambda'_m = {mp.nstr(K2 + (2*mp.sin(th/2)**2 + 2*mp.cos(th/2)**2*mp.log(mp.cos(th/2)**2)),30)}")
mp.mp.dps = 50

print("\nV6c: exact finite-M identity  -2K + sum_{m<=M} Lambda'_m = I_M  (I_M by quadrature on Re w = M+1/2)")
for name, th, q in angles:
    z = mp.sin(th/2)**2
    K2 = 2*Kalpha(th)
    for M in [-1, 0, 3, 10]:
        sig = mp.mpf(M) + mp.mpf(1)/2
        f = lambda t: mp.pi**2*Lam(mp.mpc(sig, t), th)/mp.sin(mp.pi*mp.mpc(sig, t))**2
        T = 120/(2*mp.pi - th)
        I = 2*mp.re(mp.quad(f, [0] + [T*k/12 for k in range(1, 13)]))/(2*mp.pi)
        lhs = -K2 + mp.fsum(dLam_int(m, z, tol) for m in range(0, M + 1))
        print(f"  theta={name:>6} M={M:>3}: I_M={mp.nstr(I,30):>36}  |lhs - I_M|={mp.nstr(abs(lhs-I),3)}")

