import mpmath
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
mp.dps = 50
def r_oct(phi, p):
    cc, ss = cos(phi), sin(phi)
    return (cc**p + ss**p + mpf(2)**(-p/2)*((cc+ss)**p + (cc-ss)**p))**(-1/p)
for p in [mpf(100), mpf(1000)]:
    r0 = r_oct(0, p)
    print(f"p={p}: h(0)-1 = {nstr(r0-1,8)}  (pinned, ~2^(1-p/2)/p)")
# hexagon pinning: h(pi/4) at p: (1/sqrt2)(1+2^{1-p})^{-1/p}
p = mpf(1000)
print("hex h(pi/4)-1/sqrt2 =", nstr((1/sqrt(2))*((1+2**(1-p))**(-1/p)-1),8))
with open('progress.log','a') as fh:
    fh.write('step_117: hexagon + octagon: linear 1/p deficit, non-zeta constants, pinned edge normals\n')
