"""V2. Exact length formula vs direct numerical arc length.

  L_theta(n) = 1 + cos th + R_theta(eps) n^(-eps) - sum_{m>=2} gamma_m / ( n^(m-1) ((m-1)(n-1) - 1) ),
  eps = 1/(n-1),  R_theta = 2 sin^2(th/2) (pi eps/sin pi eps) 2F1(-eps,1+eps;2;sin^2(th/2)),
  gamma_m = (-1)^m (P_{m-2}(cos th) - P_m(cos th))/(2m-1).

Direct value A: literal quadrature of int_0^1 |u + n x^(n-1) v| dx (x-variable, breakpoints).
Direct value B: the defect identity 2 - L = eps n^-eps int_0^n h_c(s) s^(eps-1) ds (for very large n).
Also: n in (1,2) (eps > 1) by analytic continuation; the removable singularity at n = 2; theta -> pi.
"""
import mpmath as mp
from common import *

setdps(50)
angles = [('pi/6', mp.pi/6), ('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3),
          ('3pi/4', 3*mp.pi/4), ('5pi/6', 5*mp.pi/6), ('2.9', mp.mpf('2.9'))]
ns = [mp.mpf('1.7'), mp.mpf('2.5'), mp.mpf(3), mp.mpf(5), mp.mpf('7.5'), mp.mpf(10), mp.mpf(100), mp.mpf(1000)]

print("V2a: exact formula vs literal x-quadrature of the arc length (dps=50)")
print(f"{'theta':>6} {'n':>6} {'L_theta(n) exact formula':>48} {'|exact - quad|':>15} {'sig.digits':>10}")
worst = 0
for name, th in angles:
    for n in ns:
        E = L_exact(n, th)
        Q = L_direct_x(n, th)
        d = abs(E - Q)
        worst = max(worst, d)
        print(f"{name:>6} {mp.nstr(n,5):>6} {mp.nstr(E,45):>48} {mp.nstr(d,3):>15} {mp.nstr(digits(E,Q),4):>10}")
print("worst |exact - quad| :", mp.nstr(worst, 3))

print("\nV2b: very large n: exact formula vs defect-identity quadrature (s = e^u), dps=50")
worst = 0
for name, th in angles[:6]:
    for n in [mp.mpf(10)**4, mp.mpf(10)**6, mp.mpf(10)**10]:
        E = L_exact(n, th)
        Q = L_direct_s(n, th)
        d = abs(E - Q)
        worst = max(worst, d)
        print(f"{name:>6} n=1e{int(mp.log10(n)):<3d} L={mp.nstr(E,45)}  |exact - quad|={mp.nstr(d,3)}")
print("worst:", mp.nstr(worst, 3))

print("\nV2c: removable singularity at n = 2 (eps = 1): formula at n = 2 +- delta averaged vs elementary L_theta(2)")
setdps(140)
for name, th in angles[:6]:
    c = mp.cos(th)
    # L(2) = int_0^1 sqrt(1 + 4 c x + 4 x^2) dx  (elementary)
    F = lambda x: ((2*x + c)*mp.sqrt(4*x*x + 4*c*x + 1)/4 + (1 - c*c)/4*mp.log(2*x + c + mp.sqrt(4*x*x + 4*c*x + 1)))
    L2 = F(1) - F(0)
    delta = mp.mpf(10)**-30
    Lf = (L_exact(2 + delta, th) + L_exact(2 - delta, th))/2
    print(f"  theta={name:>6}  L(2)={mp.nstr(L2,40)}  |formula - elementary| = {mp.nstr(abs(Lf-L2),3)}")
setdps(50)

print("\nV2d: hairpin limit theta -> pi: R_pi(eps) = 2/(1+eps), gamma_m = 0, L_pi(n) = 2 (1-1/n) n^(-1/(n-1))")
for n in [mp.mpf(3), mp.mpf(10), mp.mpf(1000)]:
    eps = 1/(n - 1)
    th = mp.pi - mp.mpf(10)**-40
    lhs = L_exact(n, th)
    rhs = 2*(1 - 1/n)*n**(-eps)
    print(f"  n={mp.nstr(n,5):>6}  formula(theta=pi-1e-40)={mp.nstr(lhs,40)}  closed={mp.nstr(rhs,40)}  diff={mp.nstr(abs(lhs-rhs),3)}")
for e in [mp.mpf('0.3'), mp.mpf('0.77')]:
    print(f"  eps={e}: R_theta at theta=pi-1e-40: {mp.nstr(R_hyp(e, mp.pi - mp.mpf(10)**-40),35)}  2/(1+eps)={mp.nstr(2/(1+e),35)}")

print("\nV2e: theta = pi/2 vs Ahmed's Theorem 2 (exact formula) at n = 3, 10, 100")
for n in [mp.mpf(3), mp.mpf(10), mp.mpf(100)]:
    eps = 1/(n - 1)
    ah = 1 + R_ahmed(eps)/n**eps - mp.nsum(lambda k: mp.binomial(mp.mpf(1)/2, k)/(n**(2*k - 1)*((2*k - 1)*(n - 1) - 1)), [1, mp.inf])
    E = L_exact(n, mp.pi/2)
    print(f"  n={mp.nstr(n,4):>4}  ours={mp.nstr(E,40)}  Ahmed={mp.nstr(ah,40)}  diff={mp.nstr(abs(E-ah),3)}")
n = mp.mpf(3)
print("  L(3) (Ahmed Example 1) = 1.547865... :", mp.nstr(L_exact(n, mp.pi/2), 20))
