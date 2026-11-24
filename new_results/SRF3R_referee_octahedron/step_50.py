import time
from mpmath import mp, mpf, quad, sqrt, log10
mp.dps = 32
t0=time.time()
def S_graph(p):
    p = mpf(p)
    def integrand(x, y):
        z = (1 - x**p - y**p)**(1/p)
        zx = -x**(p-1) * z**(1-p); zy = -y**(p-1) * z**(1-p)
        return sqrt(1 + zx**2 + zy**2)
    def inner(x):
        ymax = (1 - x**p)**(1/p)
        return ymax * quad(lambda v: integrand(x, ymax*v), [0, 1])
    return 8 * quad(inner, [0, 1])
S1 = S_graph(1.01)
print("d=0.01 graph:", mp.nstr(S1,30), f" [{time.time()-t0:.0f}s]")
t1=time.time()
def T_simplex(tau):
    tau = mpf(tau)
    f = lambda t1, t2: ((t1*t2*(1-t1-t2))**(-tau)
                        * sqrt(t1**(2*tau) + t2**(2*tau) + (1-t1-t2)**(2*tau)))
    def inner(t1):
        return (1-t1)*quad(lambda v: f(t1, (1-t1)*v), [0,1])
    return quad(inner, [0,1])
p = mpf('1.01'); tau = 1-1/p
ST = 8*(1-tau)**2*T_simplex(tau)
print("via T:", mp.nstr(ST,30), f"[{time.time()-t1:.0f}s]")
print("digits:", float(-log10(abs(S1-ST)/abs(ST))))
