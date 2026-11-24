import time
from mpmath import mp, mpf, quad, sqrt, log10, sin, cos, pi
mp.dps = 55
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
Sg = S_graph(1.01)
ST = mpf('6.997417989197152537546960063980097361868')
print("graph dps55 d=0.01:", mp.nstr(Sg,45), f"[{time.time()-t0:.0f}s]")
print("digits vs T/H1:", float(-log10(abs(Sg-ST)/abs(ST))))

# known cases for the radial formula
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
        st_ = st if st != 0 else mpf(10)**(-mp.dps)
        grad2 = drhoth**2 + (drhoph/st_)**2
        return rho*sqrt(rho**2 + grad2)*st
    return 16 * quad(lambda ph: quad(lambda th: Nrho(th, ph), [0, pi/2]), [0, pi/4])
t1=time.time()
print("p=2 radial:", mp.nstr(S_radial(2),30), "vs 4pi:", mp.nstr(4*pi,30), f"[{time.time()-t1:.0f}s]")
t2=time.time()
print("p=1 radial:", mp.nstr(S_radial(1),30), "vs 4sqrt3:", mp.nstr(4*sqrt(3),30), f"[{time.time()-t2:.0f}s]")
