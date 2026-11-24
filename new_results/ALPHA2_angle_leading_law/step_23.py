import sympy as sp
u = sp.symbols('u')
k = sp.symbols('kappa', positive=True)
S = sp.sqrt(u**2 + 4*k**2); D = sp.sqrt(1+4*k**2)
# The four pole measures mu = dz/(z-z0): derived formulas
mu_w1   = sp.Rational(1,2)*(1/S + 1/(u-1) + D/(S*(u-1)))   # pole w=1  (z=gamma)
mu_wm1  = sp.Rational(1,2)*(1/S + 1/(u+1) - D/(S*(u+1)))   # pole w=-1 (z=-gamma)
mu_wr   = sp.Rational(1,2)*(1/S + 1/(u-1) - D/(S*(u-1)))   # pole w=-r (z=-s)
mu_w0   = 1/S                                              # pole w=0
# verify against direct computation dz/(z-z0) with z=(u+S)/(2k)
z = (u+S)/(2*k)
gamma = (D+1)/(2*k); sroot = (D-1)/(2*k)
dz_du = sp.simplify(sp.diff(z, u))
def direct(z0):
    return sp.simplify(dz_du/(z - z0))
for name, z0, mu in [('w=1',gamma,mu_w1),('w=-1',-gamma,mu_wm1),('w=-r',-sroot,mu_wr),('w=0',0,mu_w0)]:
    diff = sp.simplify(direct(z0) - mu)
    print(name, "check:", diff)
