
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi
mp.dps = 60

def saddle_xs(N, Z):
    return (-Z + sqrt(Z**2 + 4*N*Z))/2

N = 1500; Z = 1j*mp.pi/2; m = 0
x1 = saddle_xs(N, Z)
def f(x):
    return exp(-x)*(1+Z/x)**(-N)*(x+Z)**(m+mpf(3)/2)/x

d2 = N*(1/x1**2 - 1/(x1+Z)**2); ad2 = abs(d2)
v = exp(-mpc(0,1)*mp.arg(d2)/2)
phis = x1 + N*log(1+Z/x1)

# scan Re phi along path piece1 (0->x1) and piece2 (x1->x1+Tv)
print('scan t along 0->x1: Re(phi-phis)')
for t in [0.001,0.01,0.05,0.1,0.2,0.4,0.6,0.8,0.9,1.0]:
    x = t*x1
    dphi = (x + N*log(1+Z/x)) - phis
    print(' t=%.3f  Re=%.3f  Im=%.3f' % (t, float(mp.re(dphi)), float(mp.im(dphi))))
T = 11.0/sqrt(ad2)
print('T*v magnitude:', float(abs(T*v)), ' |x1| =', float(abs(x1)), ' arg v =', float(mp.arg(v)))
print('scan along x1->x1+T*v')
for t in [0.0,0.2,0.5,1.0]:
    x = x1 + t*T*v
    dphi = (x + N*log(1+Z/x)) - phis
    print(' t=%.2f  Re=%.3f  Im=%.3f' % (t, float(mp.re(dphi)), float(mp.im(dphi))))

I1 = x1*mp.quadgl(lambda t: f(t*x1), [0, mpf('0.2'), mpf('0.5'), mpf('0.8'), 1])
I2 = (T*v)*mp.quadgl(lambda t: f(x1+t*T*v), [0, mpf('0.3'), mpf('0.7'), 1])
x2 = x1+T*v
I3 = mp.quadgl(lambda s: f(x2+s), [0, mpf('5'), mpf('20'), mp.inf])
model = exp(-phis)*sqrt(2*pi/d2)*(x1+Z)**(m+mpf(3)/2)/x1
print('I1 =', mp.nstr(I1, 12))
print('I2 =', mp.nstr(I2, 12))
print('I3 =', mp.nstr(I3, 12))
print('model =', mp.nstr(model, 12))
print('(I1+I2+I3)/model =', mp.nstr((I1+I2+I3)/model, 12))
print('|model| =', mp.nstr(abs(model), 6))
