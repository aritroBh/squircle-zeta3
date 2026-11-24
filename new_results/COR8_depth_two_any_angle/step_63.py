
import mpmath as mp
mp.mp.dps = 50
# PV integral definition of H_{0,0,0,1}(1.5) = int_0^1.5 log^3 t / (6(1-t)) dt (Hadamard PV at t=1)
def pv_li4(x):
    # H_{0,0,0,1}(x) = Li4(x) for x<=1; for x>1 use PV: subtract the pole regularisation
    f = lambda t: mp.log(t)**3/(6*(1-t))
    eps = mp.mpf('1e-20')
    # PV: int_0^{1-eps} + int_{1+eps}^x  + limit term
    a = mp.quad(f, [0, 1-eps])
    b = mp.quad(f, [1+eps, x])
    # the singular part: log^3(1)/(...) = 0 so PV ~ a+b + O(eps^2 log eps)
    return a + b
print("PV Li4(1.5):", mp.nstr(pv_li4(1.5), 30))
print("Re polylog :", mp.nstr(mp.re(mp.polylog(4, 1.5)), 30))
print("Im polylog :", mp.nstr(mp.im(mp.polylog(4, 1.5)), 15))
# check with the analytic inversion relation
x = mp.mpf('1.5')
inv = -mp.polylog(4, 1/x) - mp.pi**4/45 - (mp.pi**2/12)*mp.log(x)**2 + mp.log(x)**4/24
print("inversion Re value:", mp.nstr(mp.re(inv), 30))
