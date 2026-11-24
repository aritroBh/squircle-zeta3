"""hexagon_p1_check.py -- Theorem 12 and Observation 2 of "Investigating zeta(3)".

Smoothed regular 2m-gons K_p = {sum_{i=1}^m |<n_i, x>|^p <= 1}, n_i = (cos(i pi/m), sin(i pi/m));
m = 2 is the squircle. The perimeter is computed in polar coordinates (with the kinks split), the
value P(1) = 4m(sec(pi/2m) - cos(pi/2m)) is checked, and for the hexagon (m = 3) the derivatives
P'(1), P''(1), P'''(1) are read off an exact degree-14 interpolation at p = 1 + k h, |k| <= 7, and
compared with
    P'(1)   = 2 sqrt3 log 2                                                     (Theorem 12, proved)
    P''(1)  = (sqrt3/8) pi^2 + 2 sqrt3 log^2 2 - 4 sqrt3 log 2                   (Observation 2, PSLQ)
    P'''(1) = (sqrt3/8)(16 log^3 2 + 3 pi^2 log 2 - 96 log^2 2 + 96 log 2 - 9 zeta(3) - 6 pi^2)
The finite-difference check gives about 19, 17 and 13 digits. The closed forms for P'' and P''' were
found by PSLQ at 170 digits (largest coefficient 96, residual 2e-139); they are numerical, not proved.
Run: python hexagon_p1_check.py   (a few seconds)
"""
import mpmath as mp
mp.mp.dps = 50
def perim(m, p):
    ang = [i*mp.pi/m for i in range(1, m+1)]
    def f(th):
        c = [mp.cos(th - a) for a in ang]; s = [mp.sin(th - a) for a in ang]
        S = mp.fsum(abs(x)**p for x in c)
        Sp = mp.fsum(-p*abs(x)**(p-1)*mp.sign(x)*y for x, y in zip(c, s))
        rho = S**(-1/p); drho = -(1/p)*S**(-1/p-1)*Sp
        return mp.sqrt(rho**2 + drho**2)
    kinks = sorted(set([mp.mpf(0), mp.pi] + [ (a + mp.pi/2) % mp.pi for a in ang]))
    return 2*mp.quad(f, kinks)
for m in (2, 3):
    print("m", m, "P(1)", mp.nstr(perim(m, 1), 20), "closed form", mp.nstr(4*m*(mp.sec(mp.pi/(2*m)) - mp.cos(mp.pi/(2*m))), 20))
m = 3; h = mp.mpf('0.004'); K = 7
xs = [k*h for k in range(-K, K+1)]
ys = [perim(m, 1 + x) for x in xs]
# least-squares-free: exact polynomial interpolation of degree 2K, read off Taylor coefficients
coeffs = mp.lu_solve(mp.matrix([[x**j for j in range(2*K+1)] for x in xs]), mp.matrix(ys))
d1, d2, d3 = coeffs[1], 2*coeffs[2], 6*coeffs[3]
r3, L = mp.sqrt(3), mp.log(2)
cl1 = 2*r3*L
cl2 = r3/8*mp.pi**2 + 2*r3*L**2 - 4*r3*L
cl3 = r3/8*(16*L**3 + 3*mp.pi**2*L - 96*L**2 + 96*L - 9*mp.zeta(3) - 6*mp.pi**2)
for nm, a, b in (("P'", d1, cl1), ("P''", d2, cl2), ("P'''", d3, cl3)):
    print(nm, "fit", mp.nstr(a, 16), "closed form", mp.nstr(b, 16), "rel diff", mp.nstr(abs(a-b)/abs(b), 3))
print("P(1.1) own", mp.nstr(perim(3, mp.mpf('1.1')), 30), "earlier 50-digit value 3.69840169493827649453206212737602")
