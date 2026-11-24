
import sympy as sp

# ---------- Symbols ----------
sig = sp.Symbol('sig', positive=True)   # sigma = 1/sqrt(N)
Z   = sp.Symbol('Z')                    # z_s = s*i*alpha
g1, g2 = sp.symbols('g1 g2')            # local series: Btilde(z_s+zeta) = ... + gam*zeta^1/2*(1+g1*zeta+g2*zeta^2+...)
tau = sp.Symbol('tau')
m   = sp.Symbol('m', nonnegative=True, integer=True)

ORD = 6   # expand in sigma up to this order
# ---------- saddle x_s : x_s(x_s+Z) = N Z,  x_s = sqrt(NZ) y ----------
sqZ = sp.sqrt(Z)
y = (1 - sqZ/2*sig + Z/8*sig**2 - Z**2/128*sig**4 + Z**3/1024*sig**6)  # sqrt(1+Z sig^2/4)-sig sqrtZ/2
y = sp.series(y, sig, 0, ORD+1).removeO()
xs = sqZ * sig**-1 * y          # x_s
xs = sp.expand(xs)

# ---------- phi(x) = x + N log(1+Z/x), N = 1/sig^2 ; phi_s series ----------
logy = sp.series(sp.log(1 + (y-1)), sig, 0, ORD+1).removeO()   # log y
phi_s = sp.expand(sqZ*sig**-1*y - 2*sig**-2*logy)
phi_s = sp.series(phi_s, sig, 0, 3).removeO()   # through sig^2 (relative 1/N)
print("phi_s =", sp.expand(phi_s))
# check vs D13b: 2 sqrt(NZ) - Z/2 + Z^{3/2}/(12 sqrtN) - Z^2/(64N)
