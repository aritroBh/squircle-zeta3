import time
from mpmath import mp, mpf, quad, sqrt, log10, sin, cos, pi
t0=time.time()
mp.dps = 42
def S_radial(p):
    p = mpf(p)
    def Nrho(th, ph):
        st, ct = sin(th), cos(th)
        a = (st*cos(ph))**p + (st*sin(ph))**p + ct**p
        rho = a**(-1/p)
        dNth = p*(st*cos(ph))**(p-1)*ct*cos(ph) + p*(st*sin(ph))**(p-1)*ct*sin(ph) - p*ct**(p-1)*st
        dNph = p*(st*cos(ph))**(p-1)*(-st*sin(ph)) + p*(st*sin(ph))**(p-1)*(st*cos(ph))
        drhoth = -(1/p)*a**(-1/p-1)*dNth
        drhoph = -(1/p)*a**(-1/p-1)*dNph
        grad2 = drhoth**2 + (drhoph/st)**2 if st != 0 else drhoth**2
        return rho*sqrt(rho**2 + grad2)*st
    return 16 * quad(lambda ph: quad(lambda th: Nrho(th, ph), [0, pi/2]), [0, pi/4])
Sr = S_radial(1.01)
print("radial d=0.01:", mp.nstr(Sr,40), f"[{time.time()-t0:.0f}s]")
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
print("via T:", mp.nstr(ST,40), f"[{time.time()-t1:.0f}s]")
print("radial-vs-T digits:", float(-log10(abs(Sr-ST)/abs(ST))))
