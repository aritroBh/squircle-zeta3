import sympy as sp, pickle
# General edge, one nonlocal term. Parameters:
A0, A1, B, kap = sp.symbols('A0 A1 B kap')          # term values; constraint A0 = A1 - B*kap (checked per instance)
L2s, q0, V02 = sp.symbols('L2 q0 V02')               # yp.yp, y.yp at t=0, y.y at t=0
lA0, lA1, lB = sp.symbols('lA0 lA1 lB')              # logs of term values
t, X, Y, W = sp.symbols('t X Y W')
d = sp.symbols('d')
Ls = sp.sqrt(L2s)

a0 = A0*t; a1 = A1*(1-t); a2 = B*(1+kap*t)
LA0 = X + lA0; LA1 = Y + lA1; LA2 = W + lB
# F coefficients
f1 = a0*LA0 + a1*LA1 + a2*LA2
f2 = (a0*LA0**2 + a1*LA1**2 + a2*LA2**2)/2
f3 = (a0*LA0**3 + a1*LA1**3 + a2*LA2**3)/6
# S coefficients (a0' = A0, a1' = -A1, a2' = B*kap)
s1 = A0*LA0 - A1*LA1 + B*kap*LA2
s2 = A0*LA0**2 - A1*LA1**2 + B*kap*LA2**2
s3 = A0*LA0**3 - A1*LA1**3 + B*kap*LA2**3
phi1 = s1
phi2 = s2/2 - s1*f1
phi3 = s3/6 - s2*f1/2 + s1*(f1**2 - f2)
el1 = f1
el2 = f2 - f1**2/2
el3 = f3 - f1*f2 + f1**3/3
m1 = -2*el1
m2 = 2*el1 - 2*el2 + 2*el1**2
m3 = (-2*el1 + 2*el2 - 2*el3) + (-2*el1)*(2*el1 - 2*el2) + (-2*el1)**3/6
ypyp = L2s; yyp = q0 + L2s*t; yy = V02 + 2*q0*t + L2s*t**2
h1 = -2*yyp*phi1
h2 = -2*yyp*phi2 + yy*phi1**2
h3 = -2*yyp*phi3 + 2*yy*phi1*phi2
g1 = m1*L2s + h1
g2 = m2*L2s + m1*h1 + h2
g3 = m3*L2s + m2*h1 + m1*h2 + h3
E3 = sp.expand(g3/(2*Ls) - g1*g2/(8*Ls**3) + g1**3/(16*Ls**5))
print("E3 computed; expanding...")
E3e = sp.expand(E3)
P = sp.Poly(E3e, t, X, Y, W)
mons = P.terms()
print("total monomials:", len(mons))
nl = [((n,a,b,c), co) for (n,a,b,c), co in mons if c >= 1]
loc = [((n,a,b,c), co) for (n,a,b,c), co in mons if c == 0]
print("nonlocal monomials:", len(nl), " local:", len(loc))
pickle.dump({'nl': nl, 'loc': loc}, open('E3_general.pkl','wb'))
# check t-degree cap
print("max t-degree nonlocal:", max(m[0][0] for m in nl))
