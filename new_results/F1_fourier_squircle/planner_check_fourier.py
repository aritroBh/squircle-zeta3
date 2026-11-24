# Planner check of F1: Fourier coefficients a_k(p) = (1/pi) * int_0^{2pi} h_p(t) cos(4kt) dt of the squircle support function
# h_p = (|cos|^q + |sin|^q)^(1/q), q = p/(p-1), vs (-1)^k (sqrt2/pi)[-4/(16k^2-1) + z2 tau^2 + z3 tau^3 - (33/16+7k^2) z4 tau^4], tau = 1-1/p.
import mpmath as mp
mp.mp.dps = 30
def ak(p, k):
    q = p/(p-1)
    f = lambda t: (mp.cos(t)**q + mp.sin(t)**q)**(1/q)*mp.cos(4*k*t)
    return 8/mp.pi*mp.quad(f, mp.linspace(0, mp.pi/4, 9) + [mp.pi/2])/1  # 8 * int_0^{pi/4}... see below
def ak2(p, k):
    q = p/(p-1)
    f = lambda t: (mp.cos(t)**q + mp.sin(t)**q)**(1/q)*mp.cos(4*k*t)
    return (8/mp.pi)*mp.quad(f, mp.linspace(0, mp.pi/4, 17))
def pred(tau, k):
    return (-1)**k*mp.sqrt(2)/mp.pi*(-4/mp.mpf(16*k*k-1) + mp.zeta(2)*tau**2 + mp.zeta(3)*tau**3 - (mp.mpf(33)/16 + 7*k*k)*mp.zeta(4)*tau**4)
print("check p=2 a_0*pi = 2pi:", mp.nstr(ak2(2, 0)*mp.pi, 15))
for k in (0, 1, 2, 3):
    rows = []
    for tau in (mp.mpf('0.04'), mp.mpf('0.02'), mp.mpf('0.01')):
        p = 1/(1 - tau)
        a = ak2(p, k); r = a - pred(tau, k)
        rows.append((tau, r))
    print("k", k, " resid/tau^5:", [mp.nstr(r/t**5, 6) for t, r in rows], " resid/tau^4:", [mp.nstr(r/t**4, 4) for t, r in rows])
