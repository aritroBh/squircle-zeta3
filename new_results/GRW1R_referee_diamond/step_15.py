import mpmath as mp
mp.mp.dps = 60

def I_of_eta(eta):
    p = 1/(1+eta)
    # perimeter of |x|^p+|y|^p=1, A=B=1: P = 4 int_0^1 sqrt(1 + (dy/dx)^2) dx, y=(1-x^p)^(1/p)
    # use the t-parametrization: x=t^(1/p), y=(1-t)^(1/p), speed=(1+eta)*sqrt(t^(2eta)+(1-t)^(2eta))
    f = lambda t: mp.sqrt(t**(2*eta) + (1-t)**(2*eta))
    P = 4*(1+eta)*mp.quad(f, [0, 0.5, 1])
    return P/(4*mp.sqrt(2)*(1+eta))

# Taylor coefficients at 0 by high-order finite differences via mpmath diff? Use Cauchy: sample on small circle
import numpy as np
for eta in [mp.mpf('0.001')]:
    print("I(0.001) =", mp.nstr(I_of_eta(eta), 30))
# f_1 = I'(0): central difference with h=1e-4 (cancellation: need ~ 60 digits; values ~1)
for h in ['1e-2','1e-4','1e-6']:
    hh = mp.mpf(h)
    d1 = (I_of_eta(hh)-I_of_eta(-hh))/(2*hh)
    d2 = (I_of_eta(hh)-2*I_of_eta(0)+I_of_eta(-hh))/hh**2
    print(h, "I'(0)~", mp.nstr(d1,25), " I''(0)~", mp.nstr(d2,25))
print("1-2ln2 =", mp.nstr(1-2*mp.log(2),25))
print("f_2 candidates: 1.4112335... =? ")
# f_2 = (1/2)I''(0)
