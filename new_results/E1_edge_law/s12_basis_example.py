# s12: the resummation caveat, made concrete.  I = int_0^1 log^2 t * log(1+t) dt is a non-local term of the
# type that occurs on an edge with r = 2.  Canonical (zeta(3)-free) series: expand log(1+t) about t = 1,
#   log(1+t) = log 2 + sum_k (-1)^{k+1} (-(1-t)/2)^k / k ... i.e. log 2 + log(1 - (1-t)/2),
# each term int log^2 t (1-t)^k dt is RATIONAL.  Classical evaluation: 4 log2 - 6 + pi^2/6 + (3/2) zeta(3).
import mpmath as mp
mp.mp.dps = 50
def Jt(k):  # int_0^1 log^2 t (1-t)^k dt = 2 * sum_{i=0}^k C(k,i)(-1)^i/(i+1)^3  (rational)
    return 2*mp.fsum([mp.binomial(k, i)*(-1)**i/mp.mpf(i+1)**3 for i in range(k+1)])
ser = mp.log(2)*Jt(0) - mp.fsum([Jt(k)/(k*mp.mpf(2)**k) for k in range(1, 400)])
quad = mp.quad(lambda t: mp.log(t)**2*mp.log(1+t), [0, 1])
classical = 4*mp.log(2) - 6 + mp.pi**2/6 + mp.mpf(3)/2*mp.zeta(3)
print("canonical series (rational terms, log 2 only):", mp.nstr(ser, 40))
print("quadrature                                    :", mp.nstr(quad, 40))
print("4log2 - 6 + pi^2/6 + (3/2)zeta(3)             :", mp.nstr(classical, 40))
print("Jt(k) rational check, k=3:", mp.identify(Jt(3)))
