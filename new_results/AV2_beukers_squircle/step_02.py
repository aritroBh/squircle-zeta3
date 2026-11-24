
import json
from mpmath import mp, mpf, zeta, pi, log, polylog, quad, fabs, euler
mp.dps = 100
exec(open('av2core.py').read())

print("=== 1. phi(w) small-w structure: is phi_2 really 0? ===")
for wtest in ['0.5', '0.2', '0.1', '0.01']:
    w = mpf(wtest)
    print("phi(%s) = %s" % (wtest, mp.nstr(phi(w), 30)))

# least-squares fit phi(w) - 1 = p3 w^3 + p4 w^4 + p5 w^5 + p6 w^6 + p7 w^7 on tiny w
import mpmath
def fit_phi():
    ws = [mpf(10)**(-4) * (1 + k*mpf('0.23')) for k in range(8)]
    # model: phi(w)-1 = c3 w^3 + c4 w^4 + c5 w^5 + c6 w^6
    A = []; b = []
    for w in ws:
        A.append([w**3, w**4, w**5, w**6])
        b.append(phi(w) - 1)
    Am = mpmath.matrix(A); bm = mpmath.matrix(b)
    sol, res = mpmath.lu_solve(Am.T*Am, Am.T*bm), None
    return [sol[i] for i in range(4)]
c3, c4, c5, c6 = fit_phi()
z2, z3, z4, z5, z6 = [zeta(k) for k in range(2, 7)]
print("numeric phi_3 =", mp.nstr(c3, 25), " predicted -zeta(3) =", mp.nstr(-z3, 25))
print("numeric phi_4 =", mp.nstr(c4, 25), " predicted 9/4 zeta(4) =", mp.nstr(mpf(9)/4*z4, 25))
print("numeric phi_5 =", mp.nstr(c5, 25), " predicted 3 z2 z3 - 9 z5 =", mp.nstr(3*z2*z3-9*z5, 25))
print("numeric phi_6 =", mp.nstr(c6, 25), " predicted 79/35 z2^3 + 2 z3^2 =", mp.nstr(mpf(79)/35*z2**3+2*z3**2, 25))

print("=== 2. B2(p) two methods, quick checks ===")
for p in [1, 2, 4, 16]:
    s, nn = B2_series(p)
    q = B2_quad(p)
    print("p=%3d  series(%d terms)=%s" % (p, nn, mp.nstr(s, 40)))
    print("       quad           =%s   |diff|=%s" % (mp.nstr(q, 40), mp.nstr(fabs(s-q), 5)))

print("=== 3. large p: B2 vs zeta(3) ===")
for p in [256, 4096]:
    q = B2_quad(p)
    print("p=%5d  B2=%s" % (p, mp.nstr(q, 50)))
    print("        B2-zeta(3)=%s" % mp.nstr(q - zeta(3), 10))

print("=== 4. C2 = int_0^inf (phi(w)-1) w^{-3} dw (predicted B2 = zeta3 + C2 eps^2 + zeta3/2 eps^3 + ...) ===")
def C2integrand(w):
    if w == 0:
        return -zeta(3)
    return (phi(w) - 1)/w**3
C2 = quad(C2integrand, [0, 0.25, 1, 4, 40, 200, mp.inf])
print("C2 =", mp.nstr(C2, 50))
# compare with B2(p) - zeta(3) at large p
for p in [4096, 65536]:
    q = B2_quad(p)
    eps = mpf(1)/p
    resid = q - zeta(3)
    print("p=%d: resid=%s ; resid/eps^2=%s" % (p, mp.nstr(resid, 12), mp.nstr(resid/eps**2, 12)))
