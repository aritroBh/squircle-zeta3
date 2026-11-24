# R101: even path |x|^s and its odd partner.
from mpmath import mp, mpf, quad, cos, sin, pi, gamma, zeta, log, nsum, inf, nstr, cot, quadosc
mp.dps = 30
def a(n, s):  # (2/pi) int_0^pi x^s cos(nx) dx
    return 2/pi*quad(lambda x: x**s*cos(n*x), [0]+[k*pi/n for k in range(1, n+1)])
def bt(n, s):
    return 2/pi*quad(lambda x: x**s*sin(n*x), [0]+[k*pi/n for k in range(1, n+1)])
# asymptotics: a_n(s) ~ -(2/pi)Gamma(s+1) sin(pi s/2) n^{-s-1} + (-1)^n sum_j c_j n^{-2j}
# c_1 from IBP: (2/pi) * s pi^{s-1} ; c_2 = -(2/pi) s(s-1)(s-2) pi^{s-3}
for s in [mpf('2.5'), mpf('1.3'), mpf(2), mpf('2.0001')]:
    for n in [200, 201]:
        corner = -2/pi*gamma(s+1)*sin(pi*s/2)*mpf(n)**(-s-1)
        c1 = 2/pi*s*pi**(s-1); c2 = -2/pi*s*(s-1)*(s-2)*pi**(s-3)
        c3 = 2/pi*s*(s-1)*(s-2)*(s-3)*(s-4)*pi**(s-5)
        bnd = (-1)**n*(c1/mpf(n)**2 + c2/mpf(n)**4 + c3/mpf(n)**6)
        val = a(n, s)
        print("s=%s n=%d  a_n=%s  (a_n - corner - 3 boundary terms)*n^8 = %s" % (nstr(s, 6), n, nstr(val, 15), nstr((val-corner-bnd)*mpf(n)**8, 6)))
# s = 2 exactly: a_n(2) = 4(-1)^n/n^2 (no n^-3)
print("a_7(2) - 4(-1)^7/49 =", nstr(a(7, 2)+mpf(4)/49, 3))
# odd partner at s=2: bt_n(2) = 2 pi (-1)^{n+1}/n - (4/pi)(1-(-1)^n)/n^3
for n in [3, 4, 9]:
    print("bt_%d(2) closed-form check:" % n, nstr(bt(n, 2) - (2*pi*(-1)**(n+1)/n - 4/pi*(1-(-1)**n)/mpf(n)**3), 3))
mp.dps = 40
W2 = 2*pi*log(2) - 7*zeta(3)/pi
s_series = nsum(lambda n: 2*pi*(-1)**(n+1)/n - 4/pi*(1-(-1)**n)/n**3, [1, inf])
s_int = quad(lambda t: t**2*cot(t/2), [0, pi])/pi
print("sum bt_n(2) (series) =", nstr(s_series, 35))
print("(1/pi) int t^2 cot(t/2) =", nstr(s_int, 35))
print("2 pi log2 - 7 zeta3/pi =", nstr(W2, 35))
# cross parity: a_n(1) = (2/pi)((-1)^n-1)/n^2, bt_n(1) = 2(-1)^{n+1}/n
# nsum's extrapolation fails on the even/odd pattern; sum over odd n = 2k+1 explicitly (even n give 0)
cp = nsum(lambda k: (2/pi)*((-1)**(2*k+1)-1)/(2*k+1)**2*2*(-1)**(2*k+2)/(2*k+1), [0, inf])
print("sum a_n(1) bt_n(1) =", nstr(cp, 35), "  -7zeta3/pi =", nstr(-7*zeta(3)/pi, 35))
mp.dps = 30
print("bt_3(5/2) =", nstr(bt(3, mpf(5)/2), 25))
