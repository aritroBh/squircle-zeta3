import sympy as sp
u, kappa = sp.symbols('u kappa', positive=True)
# Corner substitution: u = kappa*(z - 1/z), z = (u + S)/(2*kappa), S = sqrt(u^2+4k^2)
S = sp.sqrt(u**2 + 4*kappa**2)
z = (u + S)/(2*kappa)
D = sp.sqrt(1 + 4*kappa**2)   # = gamma*kappa*2 - 1 ... gamma = (D+1)/(2k), s = (D-1)/(2k)
gamma = (D+1)/(2*kappa); sroot = (D-1)/(2*kappa)
# measures: dz/(z-z0) = z/(S*(z-z0)) du ; use identities (z-gamma)(z+sroot) = (u-1)z/kappa etc.
# Check the product identities symbolically (as rational functions of z, eliminating S)
zsym = sp.symbols('z')
# u = kappa*(z-1/z) -> verify (z-gamma)(z+sroot) = (u-1)*z/kappa
lhs1 = (zsym-gamma)*(zsym+sroot) - ((kappa*(zsym-1/zsym)-1)*zsym/kappa)
print("(z-g)(z+s)-(u-1)z/k =", sp.simplify(lhs1.expand()))
lhs2 = (zsym+gamma)*(zsym-sroot) - ((kappa*(zsym-1/zsym)+1)*zsym/kappa)
print("(z+g)(z-s)-(u+1)z/k =", sp.simplify(lhs2.expand()))
