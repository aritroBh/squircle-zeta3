from mpmath import mp, mpf, pi, sin, cos, log, quad, odefun, polylog
mp.dps = 60
import sys, time; sys.path.insert(0,'.')
from squircle import H1u, H2u, H3u, L1

t0=time.time()
def rhs(u, Y):
    return [H1u(u), H2u(u), H3u(u)]
f = odefun(rhs, 0, [0,0,0], tol=mpf('1e-50'), degree=25)
v = f(1.5)
print("solve time", time.time()-t0)
def Cl2(theta):
    return (polylog(2, mp.e**(1j*theta))).imag
def L1_closed(phi):
    from mpmath import tan
    return mpf('0.5')*(tan(phi)*log(sin(phi)) + (Cl2(pi-2*phi)+Cl2(2*phi))/2) - phi/2
print("L1(1.5) err:", mp.nstr(v[0]-L1(1.5),5))
t0=time.time(); v2 = f(mpf('1.57')); print("eval at 1.57 time", time.time()-t0)
print("L1(1.57) err:", mp.nstr(v2[0]-L1(mpf('1.57')),5))
print("L2(1.57) ~", mp.nstr(v2[1],20))
