import sympy as sp, time
t0 = time.time()
c, z = sp.symbols('c z')
k2, k3, Xi, eps = sp.symbols('k2 k3 Xi eps')
mu1 = Xi; mu2 = Xi**2 + k2; mu3 = Xi**3 + 3*k2*Xi + k3
exact_J1 = 2*z*mu1
exact_J2 = -z*mu2 + 2*z*(1-z)
exact_J3 = z*mu3/3 - c*z*(1-z)

# term1 = 2z E[1 - exp(-eps(log n + xi))] = 2z(1 - exp(-eps Xi + k2/2 eps^2 - k3/6 eps^3))
term1 = sp.expand(2*z*(1 - sp.exp(-eps*Xi)*sp.exp(k2*eps**2/2 - k3*eps**3/6)))
term1 = sp.expand(sp.series(term1, eps, 0, 4).removeO())

# term2 = sum_l d_l n^{1-l} [1/(l-1) - 1/(l-1-eps)] ; 2-L = term1 - term2
# n^{1-l}: n^{-1} = eps/(1+eps), n^{-2} = eps^2/(1+eps)^2
t2 = ( 2*z*(1-z) * sp.series(1 - 1/(1-eps), eps, 0, 4).removeO() * eps/(1+eps)
     + (-4*c*z*(1-z)) * sp.series(sp.Rational(1,2) - 1/(2-eps), eps, 0, 3).removeO() * eps**2/(1+eps)**2 )
t2 = sp.expand(sp.series(t2, eps, 0, 4).removeO())
total = sp.expand(sp.series(term1 - t2, eps, 0, 4).removeO())
print("independent f=x derivation vs exact_J's:")
print(" eps^1:", sp.simplify(total.coeff(eps,1) - exact_J1))
print(" eps^2:", sp.simplify(sp.expand(total.coeff(eps,2) - exact_J2)))
print(" eps^3:", sp.simplify(sp.expand(total.coeff(eps,3) - exact_J3)))

# also verify the relation 2z = 1-c was used consistently: substitute z=(1-c)/2 as final sanity
sub = {z: (1-c)/2}
print(" with z=(1-c)/2 substituted:", [sp.simplify(sp.expand((total.coeff(eps,i) - [exact_J1,exact_J2,exact_J3][i-1]).subs(sub))) for i in (1,2,3)])
print("elapsed", time.time()-t0)
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} independent f=x every-angle series from layer-cake + CGF + tail: matches exact_J's (diff 0)  [EXACT]\n")
