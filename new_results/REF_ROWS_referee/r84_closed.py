# R84: evaluate the 34-term table of OCT6 (typed by hand from OCT6_octagon_closed/result.md),
# compare with own P'''(1); check C4, zeta(3) bookkeeping, the A-identity, and the Li_{2,1} reductions.
import sys
from mpmath import mp, mpf, sqrt, log, pi, zeta, polylog, nstr, sin, cos, fsum
mp.dps = 45
sys.path.insert(0, '.')
from r84_octagon import coeffs
r2 = sqrt(2); u = sqrt(2+r2)
x1 = 1-r2; x2 = 1-r2/2; k = r2-1; d = 2-r2; h = r2/2; mh = -r2/2
l2 = log(2); lu = log(1+r2); z2 = pi**2/6; z3 = zeta(3)
L2 = lambda z: polylog(2, z); L3 = lambda z: polylog(3, z)
Q = lambda a, b: mpf(a) + mpf(b)*r2
F = mpf
terms = [
 (Q(5,-F(7)/2), l2), (Q(3,-2), lu),
 (Q(-F(3)/4,F(1)/2), L2(d)), (Q(F(3)/4,-F(1)/2), L2(k)), (Q(-F(7)/4,F(5)/4), L2(x1)), (Q(F(7)/4,-F(5)/4), L2(x2)),
 (Q(-F(77)/8,F(27)/4), l2*lu), (Q(-F(389)/32,F(275)/32), l2**2), (Q(-3,2), lu**2), (Q(-F(3)/4,F(1)/2), z2),
 (Q(F(3)/8,-F(1)/4), L3(d)), (Q(-F(3)/8,F(1)/4), L3(k)), (Q(-F(9)/2,3), L3(mh)), (Q(-F(155)/32,F(103)/32), L3(h)),
 (Q(F(33)/32,-F(21)/32), L3(x1)), (Q(-F(79)/32,F(59)/32), L3(x2)), (Q(-F(17)/16,F(49)/64), l2*L2(d)),
 (Q(F(11)/32,-F(15)/64), l2*L2(k)), (Q(F(41)/32,-F(7)/8), l2*L2(mh)), (Q(-F(41)/32,F(7)/8), l2*L2(h)),
 (Q(-F(5)/32,F(3)/32), l2*L2(x1)), (Q(-F(41)/16,F(57)/32), l2*L2(x2)), (Q(F(63)/32,-F(87)/64), l2*lu**2),
 (Q(F(101)/32,-F(71)/32), l2*z2), (Q(F(443)/64,-F(631)/128), l2**2*lu), (Q(F(8635)/768,-F(3053)/384), l2**3),
 (Q(F(3)/8,-F(1)/4), lu*L2(d)), (Q(-F(3)/8,F(1)/4), lu*L2(k)), (Q(-F(5)/32,F(1)/8), lu*L2(h)),
 (Q(F(63)/32,-F(45)/32), lu*L2(x1)), (Q(F(1)/16,-F(1)/32), lu*L2(x2)), (Q(F(17)/32,-F(3)/8), lu*z2),
 (Q(F(1)/2,-r2/3/r2), lu**3), (Q(F(29)/32,-F(39)/64), z3)]
# note lu^3 coefficient is 1/2 - sqrt2/3
terms[32] = (mpf(1)/2 - r2/3, lu**3)
print("number of terms:", len(terms))
S = fsum(c*v for c, v in terms)
closed = 48*u*S
own = 6*coeffs(4)[3]
print("closed form 48uS =", nstr(closed, 40))
print("own P'''(1)      =", nstr(own, 40))
print("difference:", nstr(abs(closed-own), 3))
# zeta(3) bookkeeping
zc = 48*Q(F(29)/32, -F(39)/64)
print("zeta3 coeff /u:", nstr(zc, 20), " printed 87/2-117sqrt2/4:", nstr(F(87)/2-117*r2/4, 20))
m = 4; s_ = sin(pi/(2*m)); c_ = cos(pi/(2*m))
Cm = 48*m*s_**4*c_*(3*s_**2-1)
print("C4 (thm:Cm) =", nstr(Cm, 25), " printed (54-39sqrt2)u =", nstr((54-39*r2)*u, 25))
print("nonlocal: total - C4/u =", nstr(zc-(54-39*r2), 20), " printed (39sqrt2-42)/4 =", nstr((39*r2-42)/4, 20))
# A-identity  A(x) = sum H_n x^n/n^2  for 0<x<1
def A(x):
    s = mpf(0); H = mpf(0); n = 1; t = x
    while True:
        H += mpf(1)/n; term = H*t/n**2; s += term
        if abs(term) < mpf(10)**(-mp.dps-5): break
        n += 1; t *= x
    return s
for x in [mpf('0.3'), mpf('0.5'), mpf('0.7'), x2]:
    lhs = A(x) - A(1-x)
    rhs = 2*L3(x)-2*L3(1-x)+log(1-x)*L2(1-x)-log(x)*L2(x)-log(x)**2*log(1-x)/2+log(1-x)**2*log(x)/2
    print("A-identity x=", nstr(x, 6), " diff", nstr(abs(lhs-rhs), 3))
# Li_{2,1}(x,y) = sum_{n>m>=1} x^n y^m/(n^2 m), printed definition; compare with PRV1 forms
def Li21(x, y, nmax=None):
    s = mpf(0); inner = mpf(0); xn = mpf(1); ym = mpf(1)
    n = 1
    while True:
        xn *= x
        term = xn*inner/n**2
        s += term
        if n > 20 and abs(term) < mpf(10)**(-mp.dps-5): break
        ym *= y; inner += ym/n
        n += 1
    return s
Lam = Li21(x1, mh); sw = Li21(x2, -r2)
Lam_cf = F(3)/4*z3+L3(k)+2*L3(x1)-L3(d)+l2/2*L2(k)+l2/4*(L2(x2)-L2(h))
sw_cf = L3(d)+2*L3(x2)-L3(k)-4*L3(h)-4*L3(mh)-l2/2*(L2(k)+L2(x2))+l2**2/8*(l2/2-lu)
print("Lambda", nstr(Lam, 30), " vs PRV1 form diff", nstr(abs(Lam-Lam_cf), 3))
print("swap  ", nstr(sw, 30), " vs PRV1 form diff", nstr(abs(sw-sw_cf), 3))
