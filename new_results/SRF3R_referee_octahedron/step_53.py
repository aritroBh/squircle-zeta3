import time
from mpmath import mp, mpf, quad, sqrt, log10, sin, cos, pi
mp.dps = 60
t0=time.time()
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
Sr = S_radial(1.01)
ST = mpf('6.997417989197152537546960063980097361868')
print("radial dps60:", mp.nstr(Sr,40), f"[{time.time()-t0:.0f}s]")
print("digits vs T/H1:", float(-log10(abs(Sr-ST)/abs(ST))))
