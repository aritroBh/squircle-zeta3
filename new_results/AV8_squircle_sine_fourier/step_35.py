from mpmath import mp, mpf, pi, sin, cos, log, quad, odefun, polylog, tan
mp.dps = 60
import sys, time; sys.path.insert(0,'.')
from squircle import H1u, H2u, H3u

t0=time.time()
f = odefun(lambda u,Y: [H1u(u),H2u(u),H3u(u)], 0, [0,0,0], tol=mpf('1e-50'), degree=25)
print("solve setup", time.time()-t0)
def Cl2(theta): return (polylog(2, mp.e**(1j*theta))).imag
def L1_closed(phi):
    return mpf('0.5')*(tan(phi)*log(sin(phi)) + (Cl2(pi-2*phi)+Cl2(2*phi))/2) - phi/2
t0=time.time()
errs=[]
for phi in [mpf('0.3'), mpf('0.9'), mpf('1.4'), mpf('1.5'), mpf('1.57')]:
    v = f(phi)
    errs.append(abs(v[0]-L1_closed(phi)))
print("eval time", time.time()-t0, "L1 errs:", [mp.nstr(e,4) for e in errs])
# check L2, L3 against slow quad
from squircle import L2, L3
for phi in [mpf('0.9'), mpf('1.5')]:
    v=f(phi)
    print(phi, "L2 err:", mp.nstr(v[1]-L2(phi),4), "L3 err:", mp.nstr(v[2]-L3(phi),4))
