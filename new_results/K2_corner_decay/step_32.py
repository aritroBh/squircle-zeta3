
import mpmath as mp
mp.mp.dps = 60
L = lambda x: mp.log1p(mp.exp(-x))
# R_k = ((-1)^k/k!) int_0^inf [(x+L)^k/(1+e^{-x}) - x^k] dx
for k in range(0, 9):
    if k==0:
        f = lambda x: 1/(1+mp.exp(-x)) - 1
        Ik = mp.quad(f, [0, mp.inf])
        Rk = Ik
    else:
        f = lambda x: (x+L(x))**k/(1+mp.exp(-x)) - x**k
        Ik = mp.quad(f, [0, 5, 20, mp.inf])
        Rk = ((-1)**k)*Ik/mp.factorial(k)
    print("k=%d  int=%s  R_k=%s" % (k, mp.nstr(Ik,15), mp.nstr(Rk,15)))
