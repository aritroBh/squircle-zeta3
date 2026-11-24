"""V7. (a) Mehler-Dirichlet (angular) forms: with K_th(psi) = sqrt(2(cos psi - cos th)),
        Lambda_nu(th)  = (2/pi) int_0^th cos((nu+1/2) psi) K_th(psi) dpsi
        C(th)          = (2/pi) int_0^th psi sin(psi/2) K_th(psi) dpsi                 (Ahmed's constant)
        -2 K(pi - th)  = (2/pi) int_0^th psi /(2 sin(psi/2)) K_th(psi) dpsi            (squircle corner constant)
        2 sin^2(th/2)  = (2/pi) int_0^th cos(psi/2) K_th(psi) dpsi                    (log coefficient)
    (b) J(psi) = (pi/2) int cosh(psi t)/cosh^2(pi t) dt = psi/(2 sin(psi/2))   (the kernel used in the proof)
    (c) scaled corner x u + a x^n v: L = a + c + R(eps)(an)^-eps + sum_m (-1)^m eps/(m-eps) Lambda_m (an)^-m
    (d) second-order form of Ahmed's Theorem C
    (e) table of constants at the four test angles, 40 digits
"""
import mpmath as mp
from common import *
from closed_forms import closed_a

setdps(50)
angles = [('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3), ('3pi/4', 3*mp.pi/4)]


def Lam(nu, th):
    z = mp.sin(th/2)**2
    return 2*z*mp.hyp2f1(-nu, 1 + nu, 2, z)


def Kal(th):
    k = mp.sin(th/2)
    return -(1 + k)*mp.log(1 + k) - (1 - k)*mp.log(1 - k)


def Cth(th):
    z = mp.sin(th/2)**2
    return 2*z + 2*(1 - z)*mp.log(1 - z)


print("V7a: angular (Mehler-Dirichlet) representations, residuals")
for name, th in angles + [('1', mp.mpf(1)), ('2.9', mp.mpf('2.9'))]:
    Kf = lambda p: mp.sqrt(2*(mp.cos(p) - mp.cos(th)))
    ang = lambda g: 2/mp.pi*mp.quad(lambda p: g(p)*Kf(p), [0, th/2, th])
    r1 = max(abs(ang(lambda p: mp.cos((nu + mp.mpf(1)/2)*p)) - Lam(nu, th)) for nu in [mp.mpf('0.3'), mp.mpf('2.5'), mp.mpf(3)])
    r2 = abs(ang(lambda p: p*mp.sin(p/2)) - Cth(th))
    r3 = abs(ang(lambda p: p/(2*mp.sin(p/2)) if p != 0 else mp.mpf(1)) + 2*Kal(th))
    r4 = abs(ang(lambda p: mp.cos(p/2)) - 2*mp.sin(th/2)**2)
    print(f"  theta={name:>6}: |Lambda|={mp.nstr(r1,3)}  |C|={mp.nstr(r2,3)}  |-2K|={mp.nstr(r3,3)}  |2sin^2|={mp.nstr(r4,3)}")

print("\nV7b: kernel J(psi) = (pi/2) int_R cosh(psi t)/cosh^2(pi t) dt vs psi/(2 sin(psi/2))")
for p in [mp.mpf('0.5'), mp.mpf(2), mp.mpf(3), mp.mpf('5.5')]:
    J = mp.pi*mp.quad(lambda t: mp.cosh(p*t)/mp.cosh(mp.pi*t)**2, [0, 5, 20, mp.inf])
    print(f"  psi={mp.nstr(p,3):>4}: J={mp.nstr(J,40)}  |J - psi/(2sin(psi/2))|={mp.nstr(abs(J - p/(2*mp.sin(p/2))),3)}")

print("\nV7c: scaled corner x u + a x^n v  (exact formula vs literal quadrature)")
for name, th in angles:
    c = mp.cos(th)
    for a, n in [(mp.mpf(3), mp.mpf(5)), (mp.mpf('0.5'), mp.mpf(7)), (mp.mpf('0.25'), mp.mpf('9.5'))]:
        eps = 1/(n - 1)
        an = a*n
        # tail sum_{m>=1} (-1)^m eps/(m-eps) Lambda_m (an)^-m, Lambda_m = (P_{m-1}-P_{m+1})/(2m+1), |Lambda_m| <= 2
        P = [mp.mpf(1), c]
        tail = mp.mpf(0)
        m = 1
        while True:
            while len(P) <= m + 1:
                k = len(P) - 1
                P.append(((2*k + 1)*c*P[k] - k*P[k - 1])/(k + 1))
            tail += (-1)**m*eps/(m - eps)*(P[m - 1] - P[m + 1])/(2*m + 1)*an**(-m)
            if m > 3 and 2*eps/(m - eps)*an**(-m)*an/(an - 1) < mp.mpf(10)**(-mp.mp.dps - 3):
                break
            m += 1
        E = a + c + R_hyp(eps, th)*an**(-eps) + tail
        f = lambda x: mp.sqrt(1 + 2*c*a*n*x**(n - 1) + (a*n*x**(n - 1))**2)
        pts = [0] + [((mp.mpf(2)**j)/(a*n))**eps for j in range(-30, int(mp.log(a*n, 2)) + 1)] + [1]
        pts = sorted(set([p for p in pts if 0 <= p <= 1]))
        Q = mp.quad(f, pts)
        print(f"  theta={name:>6} a={mp.nstr(a,3):>4} n={mp.nstr(n,3):>4}: L={mp.nstr(E,40)}  |exact - quad|={mp.nstr(abs(E-Q),3)}")

print("\nV7d: Theorem C to second order: 2 - L = (2z log n + C)/(n-1) - (z (log n - a1)^2 + Q2)/(n-1)^2 + O(log^3 n/n^3)")
for name, th in angles:
    z = mp.sin(th/2)**2
    A = closed_a(th)
    Q2 = 2*z*(A[1] - (1 - z))
    for n in [mp.mpf(10)**4, mp.mpf(10)**6, mp.mpf(10)**8]:
        eps = 1/(n - 1)
        lg = mp.log(n)
        D = 2 - L_exact(n, th)
        approx = (2*z*lg + Cth(th))*eps - (z*(lg - A[0])**2 + Q2)*eps**2
        print(f"  theta={name:>6} n=1e{int(mp.log10(n))}: error={mp.nstr(D-approx,4):>11}   error*n^3/log^3 n={mp.nstr((D-approx)*n**3/lg**3,5)}")

print("\nV7e: table of constants (40 digits)")
phi = (1 + mp.sqrt(5))/2
for name, th in angles + [('golden', mp.acos(2 - mp.sqrt(5)))]:
    z = mp.sin(th/2)**2
    A = closed_a(th)
    Q2 = 2*z*(A[1] - (1 - z))
    Q3c = 2*z*(A[1]*(2 + A[0]) + A[2] - (3 - z - 2*z**2)/2)
    print(f"  theta={name}")
    print(f"    2sin^2(th/2) = {mp.nstr(2*z,40)}")
    print(f"    C(theta)     = {mp.nstr(Cth(th),40)}")
    print(f"    kappa=e^a1   = {mp.nstr(mp.exp(A[0]),40)}")
    for j in range(4):
        print(f"    a{j+1}           = {mp.nstr(A[j],40)}")
    print(f"    Q2           = {mp.nstr(Q2,40)}")
    print(f"    Q3 = {mp.nstr(Q3c,40)} - ({mp.nstr(2*z*A[1],40)}) log n")
    print(f"    K(pi-theta)  = {mp.nstr(Kal(th),40)}")
