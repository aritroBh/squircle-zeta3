"""V3. Taylor coefficients a_j(theta) of log R_theta(eps), closed forms with polylogarithms,
and the all-order multiple-polylogarithm expansion.

Reference values: Taylor coefficients of R_theta(eps)/(2z) computed by the trapezoidal rule on the
circle |eps| = 1/2 (R_theta is analytic in |eps|<1), using mpmath's hyp2f1 at complex eps; then the
power-series logarithm.  This is independent of the closed forms being tested.

Closed forms (z = sin^2(th/2), w = cot^2(th/2) = (1-z)/z, l = log(1-z) = 2 log cos(th/2)):
  a1 = -1 - w l
  a2 = zeta(2) + 1/2 - Li2(z) - (w l)^2/2
  a3 = 2 w S12(z) - (1 + (w l)^3)/3
  a4 = zeta(4)/2 - 2 S22(z) + 2 w^2 l S12(z) + (1 - (w l)^4)/4
with Nielsen polylogarithms S12(z) = 1/2 int_0^z log^2(1-t) dt/t,  S22(z) = int_0^z S12(t) dt/t.

All-order: (1+eps) 2F1(-eps,1+eps;2;z) = sum_r (-eps^2)^r [ Li_{2^r}(z) + eps w Li_{1,2^r}(z) ],
  Li_{2^r}(z) = sum_{k1>...>kr>=1} z^k1/(k1^2...kr^2),  Li_{1,2^r}(z) = sum_{k>k1>...>kr} z^k/(k k1^2...kr^2).
"""
import mpmath as mp
from common import *

setdps(70)
angles = [('pi/6', mp.pi/6), ('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3),
          ('3pi/4', 3*mp.pi/4), ('5pi/6', 5*mp.pi/6), ('1', mp.mpf(1)), ('2.9', mp.mpf('2.9'))]
JMAX = 6


def taylor_R_over_2z(th, N=512, rho=mp.mpf(1)/2):
    z = mp.sin(th/2)**2
    vals = []
    for k in range(N):
        e = rho*mp.expjpi(2*mp.mpf(k)/N)
        vals.append((mp.pi*e/mp.sin(mp.pi*e))*mp.hyp2f1(-e, 1 + e, 2, z))
    coeffs = []
    for j in range(JMAX + 1):
        s = mp.fsum(vals[k]*mp.expjpi(-2*mp.mpf(j*k)/N) for k in range(N))
        coeffs.append(mp.re(s/N/rho**j))
    return coeffs


def series_log(cf):
    # log of a power series with cf[0] = 1 : b' = c'/c
    b = [mp.mpf(0)]*len(cf)
    for j in range(1, len(cf)):
        b[j] = cf[j] - mp.fsum(k*b[k]*cf[j - k] for k in range(1, j))/j
    return b


from closed_forms import S12, S12_classical, S22, closed_a


def mpl_tables(z, R, tol):
    """Li_{2^r}(z) and Li_{1,2^r}(z) for r = 0..R by nested sums (e_r = elementary symmetric in 1/j^2)."""
    e = [mp.mpf(1)] + [mp.mpf(0)]*R           # e_r(k-1)
    L2 = [mp.mpf(0)]*(R + 2)                  # Li_{2^(r+1)}
    L1 = [mp.mpf(0)]*(R + 1)                  # Li_{1,2^r}
    zk = mp.mpf(1)
    k = 1
    while True:
        zk *= z
        for r in range(R + 1):
            L1[r] += zk*e[r]/k
            L2[r + 1] += zk*e[r]/k**2
        # update e_r(k) = e_r(k-1) + e_{r-1}(k-1)/k^2   (descending r)
        for r in range(R, 0, -1):
            e[r] += e[r - 1]/k**2
        if zk < tol:
            break
        k += 1
    L2[0] = mp.mpf(1)
    return L2[:R + 1], L1


print("V3a: Taylor coefficients a_j of log R_theta(eps) (Cauchy/trapezoid, N=512) vs closed forms (dps=70)")
worst = [0]*5
store = {}
for name, th in angles:
    cf = taylor_R_over_2z(th)
    b = series_log(cf)
    b[1] -= 0  # log R = log(2z) + sum_{j>=1} b_j eps^j
    ca = closed_a(th)
    store[name] = b
    print(f" theta={name}")
    for j in range(1, 5):
        d = abs(b[j] - ca[j - 1])
        worst[j] = max(worst[j], d)
        print(f"   a{j} = {mp.nstr(b[j],50):>58}   |closed - Cauchy| = {mp.nstr(d,3)}")
    print(f"   a5 = {mp.nstr(b[5],40)}   a6 = {mp.nstr(b[6],40)}   (numerical only)")
print("worst residuals a1..a4:", [mp.nstr(x, 3) for x in worst[1:]])

print("\nV3b: Nielsen S12: quadrature vs classical reduction (Li3(1-z), Li2(1-z), logs, zeta(3))")
for name, th in angles[:6]:
    z = mp.sin(th/2)**2
    print(f"   theta={name:>6}  S12(z)={mp.nstr(S12(z),45)}  |diff|={mp.nstr(abs(S12(z)-S12_classical(z)),3)}")

print("\nV3c: theta = pi/2 consistency with Ahmed's Lemma 4 (a1=log2-1, a2=pi^2/12+1/2, a3=zeta(3)/4-1/3)"
      "\n     and the Gamma-function value a4 = zeta(4)/4 + 1/4 (from Ahmed's R)")
ca = closed_a(mp.pi/2)
refs = [mp.log(2) - 1, mp.pi**2/12 + mp.mpf(1)/2, mp.zeta(3)/4 - mp.mpf(1)/3, mp.zeta(4)/4 + mp.mpf(1)/4]
for j in range(4):
    print(f"   a{j+1}: closed={mp.nstr(ca[j],45)}  Ahmed/Gamma={mp.nstr(refs[j],45)}  |diff|={mp.nstr(abs(ca[j]-refs[j]),3)}")
z = mp.mpf(1)/2
print("   implied identity S22(1/2) = zeta(4)/8 - log2 zeta(3)/8 + log^4 2/24 :  |diff| =",
      mp.nstr(abs(S22(z) - (mp.zeta(4)/8 - mp.log(2)*mp.zeta(3)/8 + mp.log(2)**4/24)), 3))

print("\nV3d: all-order MPL expansion of G(eps) = (1+eps) 2F1(-eps,1+eps;2;z), truncated at r = 30")
for name, th in angles:
    z = mp.sin(th/2)**2
    w = (1 - z)/z
    L2, L1 = mpl_tables(z, 30, mp.mpf(10)**-75)
    for e in [mp.mpf('0.3'), mp.mpf('0.9'), mp.mpf('2.5')]:
        G = (1 + e)*mp.hyp2f1(-e, 1 + e, 2, z)
        Gs = mp.fsum((-e*e)**r*(L2[r] + e*w*L1[r]) for r in range(31))
        print(f"   theta={name:>6} eps={mp.nstr(e,2):>4}  G={mp.nstr(G,40):>44}  |series - G| = {mp.nstr(abs(G-Gs),3)}")
    # Legendre function itself: P_eps(1-2z) = sum_r (-eps^2)^r [Li_{2^r}(z) - eps Li_{1,2^r}(z)]
    e = mp.mpf('0.7')
    P = mp.legendre(e, 1 - 2*z)
    Ps = mp.fsum((-e*e)**r*(L2[r] - e*L1[r]) for r in range(31))
    print(f"   theta={name:>6} eps=0.7   P_eps(cos th)={mp.nstr(P,40)}  |series - P| = {mp.nstr(abs(P-Ps),3)}")
    # MPL reductions at weight 3, 4
    print(f"      Li_(1,2)(z) + log(1-z)Li2(z) + 2 S12(z) = {mp.nstr(L1[1] + mp.log(1-z)*mp.polylog(2,z) + 2*S12_classical(z),3)};"
          f"   Li_(2,2)(z) - Li2(z)^2/2 + 2 S22(z) = {mp.nstr(L2[2] - mp.polylog(2,z)**2/2 + 2*S22(z),3)}")

print("\nV3e: golden corner  sin^2(th/2) = 1/phi  (cos th = 2 - sqrt5, th = 103.65 deg): all of a1, a2, a3 classical")
phi = (1 + mp.sqrt(5))/2
th = mp.acos(2 - mp.sqrt(5))
z = mp.sin(th/2)**2
print("   z - 1/phi =", mp.nstr(z - 1/phi, 3), "  theta (deg) =", mp.nstr(th*180/mp.pi, 12))
lp = mp.log(phi)
g_a1 = -1 + 2*lp/phi
g_a2 = mp.pi**2/15 + mp.mpf(1)/2 + lp**2 - 2*lp**2/phi**2
g_a3 = 2*mp.zeta(3)/(5*phi) - 4*lp**3/(3*phi) - mp.mpf(1)/3 + 8*lp**3/(3*phi**3)
cf = series_log(taylor_R_over_2z(th))
for j, g in enumerate([g_a1, g_a2, g_a3], start=1):
    print(f"   a{j}: golden closed form = {mp.nstr(g,50)}   |diff vs Cauchy| = {mp.nstr(abs(g-cf[j]),3)}")
print("   S12(1/phi) - (zeta(3)/5 - 2/3 log^3 phi) =", mp.nstr(S12(1/phi) - (mp.zeta(3)/5 - 2*lp**3/3), 3))
print("   Li2(1/phi) - (pi^2/10 - log^2 phi)        =", mp.nstr(mp.polylog(2, 1/phi) - (mp.pi**2/10 - lp**2), 3))

print("\nV3f: supplementary-angle reflection  a2(th) + a2(pi-th) = zeta(2) + 1 + log z log(1-z) - (w l)^2/2 - (l'/w)^2/2")
for name, th in angles[:6]:
    z = mp.sin(th/2)**2
    w = (1 - z)/z
    l, lp_ = mp.log(1 - z), mp.log(z)
    lhs = closed_a(th)[1] + closed_a(mp.pi - th)[1]
    rhs = mp.zeta(2) + 1 + lp_*l - (w*l)**2/2 - (lp_/w)**2/2
    print(f"   theta={name:>6}  lhs={mp.nstr(lhs,40)}  |lhs-rhs|={mp.nstr(abs(lhs-rhs),3)}")
