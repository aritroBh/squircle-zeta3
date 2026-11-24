import mpmath as mp
mp.mp.dps = 50
rho = mp.pi/2
def E_of_s(s):
    # E(is) = -(1/s^2)*(sqrt((1-exp(-2 i s))/(2 i s)) - 1 + i s/2), principal sqrt
    if s == 0:
        return mp.mpf(5)/24
    w = 1j*s
    q = (1-mp.e**(-2*w))/(2*w)
    return -(mp.sqrt(q) - 1 + w/2)/s**2

# scan
import numpy as np
ss = [i*2.9845130/4000 for i in range(0,4001)]
vals = []
maxv = 0; smax=0
for s in ss:
    v = abs(E_of_s(mp.mpf(s)))**2
    vals.append((s,v))
    if v > maxv: maxv, smax = v, s
print("max |E|^2 on grid:", mp.nstr(mp.mpf(maxv),25), "at s=", smax)
print("(5/24)^2 =", mp.nstr((mp.mpf(5)/24)**2,25))
# check derivative behavior near 0: values
for s in [0, 1e-3, 1e-2, 0.05, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 2.9, 2.9845]:
    e = E_of_s(mp.mpf(s))
    print(f"s={s:8.4f}  |E|={mp.nstr(abs(e),20)}  Re={mp.nstr(mp.re(e),15)} Im={mp.nstr(mp.im(e),15)}")
with open('progress.log','a') as f:
    f.write("step_02: numeric scan of |E(is)| on slit; see stdout\n")
