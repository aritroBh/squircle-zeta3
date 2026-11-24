
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 40
from fractions import Fraction as Fr
import reflib as R

zeta = mp.e**(mp.pi*mp.j/3)
# claimed pullback for map 0, sig=1: letter(4) - letter(2)
claimed = R.hpl_form_pullback(0, 1)
print("claimed letter vec:", claimed)
def letterval(j, u):
    return zeta**j / (1 - zeta**j * u)
u = mp.mpf(1)/3
cv = sum(c * letterval(w[0], u) for w, c in claimed.items())
print("claimed value at u=1/3:", mp.nstr(cv, 15))
# direct: g = f^{-1}, f(z) = (zeta^2 z + 1)/(zeta^4 z + 1)
f = lambda z: (zeta**2*z + 1)/(zeta**4*z + 1)
# inverse: u = f(z) -> z = g(u): solve u(zeta^4 z + 1) = zeta^2 z + 1: z(u zeta^4 - zeta^2) = 1 - u
g = lambda u: (1 - u)/(u*zeta**4 - zeta**2)
h = mp.mpf("1e-8")
dlog_gm1 = (mp.log(g(u+h)-1) - mp.log(g(u)-1))/h
print("-d/du log(g(u)-1) at u=1/3:", mp.nstr(-dlog_gm1, 15))
# also check the actual map used by step_10: mobius_from_3 + inv_map, numerically
# build M from images f(0)=1, f(1)=zeta^2, f(inf)=zeta^4 and invert
# g should satisfy g(f(z)) = z:
print("consistency g(f(1/2)) =", mp.nstr(g(f(mp.mpf(1)/2)), 10))
