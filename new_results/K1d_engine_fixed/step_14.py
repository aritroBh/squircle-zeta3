step_40 = r'''
import sympy as sp, pickle, time
t0=time.time()
# Rebuild of K1_nonlocal_law/step_33.py with the ONE-LINE FIX: 1/8 -> 1/4 in final sqrt assembly
A0, A1, B, kap = sp.symbols('A0 A1 B kap')
L2s, q0, V02 = sp.symbols('L2 q0 V02')
lA0, lA1, lB = sp.symbols('lA0 lA1 lB')
t, X, Y, W = sp.symbols('t X Y W')
Ls = sp.sqrt(L2s)

a0 = A0*t; a1 = A1*(1-t); a2 = B*(1+kap*t)
LA0 = X + lA0; LA1 = Y + lA1; LA2 = W + lB
f1 = a0*LA0 + a1*LA1 + a2*LA2
f2 = (a0*LA0**2 + a1*LA1**2 + a2*LA2**2)/2
f3 = (a0*LA0**3 + a1*LA1**3 + a2*LA2**3)/6
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
E3 = sp.expand(g3/(2*Ls) - g1*g2/(4*Ls**3) + g1**3/(16*Ls**5))   # FIXED: 1/4
print("E3 computed (fix applied); expanding...", "t=%.1fs"%(time.time()-t0))
E3e = sp.expand(E3)
P = sp.Poly(E3e, t, X, Y, W)
mons = P.terms()
print("total monomials:", len(mons))
nl = [((n,a,b,c), co) for (n,a,b,c), co in mons if c >= 1]
loc = [((n,a,b,c), co) for (n,a,b,c), co in mons if c == 0]
print("nonlocal monomials:", len(nl), " local:", len(loc))
print("max t-degree nonlocal:", max(m[0][0] for m in nl))
pickle.dump({'nl': nl, 'loc': loc}, open('E3_fixed.pkl','wb'))
print("saved E3_fixed.pkl   (%.1f s)"%(time.time()-t0))
'''
open('step_40.py','w').write(step_40)
print("written step_40.py")
