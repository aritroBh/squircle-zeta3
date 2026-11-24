# numerical test of the first link: Theta_n(k^2) (from the stored D~ zeta-product data) against the actual Fourier coefficient of the support function
from common import *
import mpmath as mp
mp.mp.dps = 45
Dtarget = pickle.load(open(T+'Dtilde_relations.pkl', 'rb'))
ZB2 = {2: [(2,)], 3: [(3,)], 4: [(4,)]}; ZBall = {**ZB2, **ZB}
def theta(n, u):
    s = mp.mpf(0)
    for m in range(n//2):
        c = Dtarget[f"{n},{m}"]
        Dv = sum(mp.mpf(F(x).numerator)/F(x).denominator*mp.fprod([mp.zeta(k) for k in ns]) for x, ns in zip(c, ZBall[n]))
        s += Dv*(-4*u)**m/fact(2*m)
    return s
def resid(k, tau):
    tau = mp.mpf(tau)
    def f(s):
        th = mp.atan(mp.e**(-tau*s))
        return mp.cos(th)*mp.cos(4*k*th)*((1+mp.e**(-s))**tau-1)*tau*mp.e**(-tau*s)/(1+mp.e**(-2*tau*s))
    Fm = (-1)**k*4*mp.sqrt(2)*mp.quad(f, [0, 2, 8, 24, 60, 120, 240])*(-1)**k*(-1)**k   # (pi/sqrt2)(a_k-A_k)(-1)^k /(-1)^k ... see below
    return Fm
for k in range(0, 4):
    print("k =", k)
    for tau in (mp.mpf(1)/16, mp.mpf(1)/32, mp.mpf(1)/64):
        # (-1)^k (pi/sqrt2)(a_k - A_k) = B*(pi/sqrt2) = 4*sqrt(2)*int(...)   (K3: a_k-A_k=(-1)^k B, B=(8/pi) int_0^{pi/4}(h_p-h_1)cos(4k.theta)(-1)^k...)
        th = lambda s: mp.atan(mp.e**(-tau*s))
        f = lambda s: mp.cos(th(s))*mp.cos(4*k*th(s))*((1+mp.e**(-s))**tau-1)*tau*mp.e**(-tau*s)/(1+mp.e**(-2*tau*s))
        val = 4*mp.sqrt(2)*mp.quad(f, [0, 2, 8, 24, 60, 120, 240])      # = (pi/sqrt2) * (8/pi)/... : (pi/sqrt2)*(4/pi)*2*int
        val *= (-1)**k    # a_k-A_k = (4/pi)*2*int_0^{pi/4}(h_p-h_1)cos4k th ; times (-1)^k(pi/sqrt2)... (sign (-1)^k cancels only for even k): keep explicit
        # direct: (pi/sqrt2)(-1)^k (a_k-A_k) with a_k-A_k = (8/pi) int_0^{pi/4} (h_p-h_1) cos(4k th) d th
        direct = (mp.pi/mp.sqrt(2))*(-1)**k*(8/mp.pi)*mp.quad(f, [0, 2, 8, 24, 60, 120, 240])
        R8 = direct - sum(theta(n, k*k)*tau**n for n in range(2, 9))
        R9 = R8 - theta(9, k*k)*tau**9
        print(f"  tau={mp.nstr(tau,4)}:  R8/tau^9 = {mp.nstr(R8/tau**9, 12)}  (Theta_9 = {mp.nstr(theta(9,k*k),12)});  R9/tau^10 = {mp.nstr(R9/tau**10, 8)}")
