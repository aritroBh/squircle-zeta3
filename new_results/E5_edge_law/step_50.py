code = r'''# step_08.py: EXACT P^(k)(1), k<=5, for the rhombus |x|^p + |y/2|^p <= 1 (vertices +-(1,0), +-(0,2);
# power convention; single edge type, geometry sqrt5). Checks vs the edge laws; also prepares the
# PSLQ target P^(5)(1)/sqrt5 (a rational combination of {1,z2,z3,z2^2,z2z3,z5}).
import sympy as sp, pickle, sys
sys.path.insert(0, 'steps')
from step_03 import J, z2, z3, z5

t, al, be, L, h, s0 = sp.symbols('t alpha beta L h s0', real=True)
l0, l1, X, Y = sp.symbols('l0 l1 X Y')
R = pickle.load(open('steps/step_01_reduced.pkl', 'rb'))
loc = {'t': t, 'alpha': al, 'beta': be, 'L': L, 'h': h, 's0': s0, 'l0': l0, 'l1': l1}
E = [sp.sympify(s, locals=loc) for s in R['E']]

def edge_exact_int(aval):
    out = []
    for k in range(6):
        ex = sp.expand(E[k].subs({al: aval['alpha'], be: aval['beta'], L: aval['L'],
                                  h: aval['h'], s0: aval['s0'], l0: X, l1: Y}))
        P = sp.Poly(ex, t, X, Y)
        tot = 0
        for (n, a, b), c in P.terms():
            tot += c*J(n, 0, a, b)
        out.append(sp.expand(tot))
    return out

r5 = sp.sqrt(5)
eA = dict(alpha=1, beta=1, L=r5, h=2/r5, s0=-1/r5)   # edge (1,0)->(0,2)
IA = edge_exact_int(eA)
Pk = [sp.expand(sp.factorial(k)*4*IA[k]) for k in range(6)]
for k in range(6):
    print("P^(%d)(1) =" % k, sp.collect(Pk[k], [z5, z2*z3, z3, z2]))
z5c  = sp.simplify(sp.expand(Pk[5]).coeff(z5))
z23c = sp.simplify(sp.expand(Pk[5]).coeff(z2*z3))

u, v = sp.symbols('u v')
w1 = al*u + be*v; e = al*be*u*v; q2m = al**2*u**2 + be**2*v**2
q2 = al**2 + be**2; p = al*be; s1 = al + be
A5  = 20*w1**3 - 20*e*w1 - 32*q2m - 48*e + 15*w1 - 2
B5  = -15*q2*w1 + 15*p*w1 + 5*p*s1 + 8*q2 - 12*p
A23 = -60*e*w1 + 12*q2m + 72*e - 12*w1 + 2
B23 = -45*p*w1 + 15*p*s1 - 3*q2 + 18*p
Z5law  = 360*al*be*h**2/L*(A5 + h**2/L**2*B5)
Z23law = 120*al*be*h**2/L*(A23 + h**2/L**2*B23)
u0 = -eA['s0']/eA['L']; v0 = 1 - u0
sub = {al: 1, be: 1, L: eA['L'], h: eA['h'], u: u0, v: v0}
Z5sum  = sp.simplify(4*Z5law.subs(sub)); Z23sum = sp.simplify(4*Z23law.subs(sub))
print("z5  coeff:", z5c, " 4 Z_e^(5):", Z5sum, " diff:", sp.simplify(z5c - Z5sum))
print("z2z3 coeff:", z23c, " 4 Z_e^(2,3):", Z23sum, " diff:", sp.simplify(z23c - Z23sum))
# PSLQ target: P5/sqrt(5) as rational vector on {1,z2,z3,z2^2,z2z3,z5}
P5s = sp.expand(Pk[5]/r5)
vec = [sp.simplify(P5s.coeff(b_)) for b_ in [1, z2, z3, z2**2, z2*z3, z5]]
print("P^(5)(1)/sqrt5 coefficients on [1,z2,z3,z2^2,z2z3,z5]:", vec)
pickle.dump(dict(Pk=[str(x) for x in Pk], vec=[str(x) for x in vec]), open('steps/step_08_exact.pkl', 'wb'))
with open('steps/step_08_out.txt', 'w') as f:
    for k in range(6):
        f.write("P^(%d)(1) = %s\n" % (k, sp.collect(Pk[k], [z5, z2*z3, z3, z2])))
    f.write("vec P5/sqrt5 on [1,z2,z3,z2^2,z2z3,z5] = %s\n" % vec)
print("saved steps/step_08_exact.pkl / step_08_out.txt")
'''
open('steps/step_08.py','w').write(code)
import subprocess, time
t0=time.time()
r=subprocess.run(['<project>/.venv/bin/python','steps/step_08.py'],capture_output=True,text=True,timeout=3000)
print("wall: %.1f s"%(time.time()-t0))
print(r.stdout); print("ERR:", r.stderr[-2000:])
