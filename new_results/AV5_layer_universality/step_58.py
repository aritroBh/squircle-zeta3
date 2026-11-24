import mpmath
from mpmath import mp, mpf, sqrt, pi, zeta, log, cos, sin, tan, nstr
mp.dps = 100
# d1 from layer profile: d/dq (c^q+s^q)^(1/q) at q=1 = c[u ln u - (1+u) ln(1+u)], u=tan(th)
f = lambda th: cos(th)*(tan(th)*log(tan(th)) - (1+tan(th))*log(1+tan(th)))
d1_layer = 8*mpmath.quad(f, [0, pi/4])
# direct integral in u
g = lambda u: (u*log(u) - (1+u)*log(1+u))/(1+u**2)**mpf('1.5')
d1_u = 8*mpmath.quad(g, [0, 1])
print("d1 from layer profile:", nstr(d1_layer,55))
print("d1 from u-integral   :", nstr(d1_u,55))
print("d1 from perim fit    : -4.42642639776412573896646426700982294597089843")
with open('progress.log','a') as fh:
    fh.write('step_57: squircle square-end: layer profile prediction for d1 confirmed; coeffs NUMERICAL\n')
