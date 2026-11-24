# step_06.py: EXACT P^(k)(1), k<=5, for the sheared parallelogram with vertices +-(2,0), +-(1,2)
# (power convention K_p = { |<v0,x>|^p + |<v1,x>|^p <= 1 }, v0=(0,1/2), v1=(-1/2,1/4)).
# Parallelograms are fully local (m=2: no non-local normals, tau=gam=0; mu0=mu1=0 in power
# convention), so the exact local integral of the reduced E-series IS the exact P^(k)(1).
# Checks: z5 coefficient == sum_e Z_e^(5); z2z3 == sum_e Z_e^(2,3); z3 of P3(1) == 6S (E1).
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

r5, r13 = sp.sqrt(5), sp.sqrt(13)
eA = dict(alpha=1, beta=1, L=r5, h=4/r5, s0=-2/r5)
eB = dict(alpha=1, beta=1, L=r13, h=4/r13, s0=-7/r13)
IA = edge_exact_int(eA); IB = edge_exact_int(eB)
Pk = [sp.expand(sp.factorial(k)*(2*IA[k] + 2*IB[k])) for k in range(6)]
for k in range(6):
    print("P^(%d)(1) =" % k, sp.collect(Pk[k], [z5, z2*z3, z3, z2]))
z5c  = sp.simplify(sp.expand(Pk[5]).coeff(z5))
z23c = sp.simplify(sp.expand(Pk[5]).coeff(z2*z3))
z3c  = sp.simplify(sp.expand(Pk[3]).coeff(z3))
z2of3 = sp.simplify(sp.expand(Pk[3]).coeff(z2))

u, v = sp.symbols('u v')
w1 = al*u + be*v; e = al*be*u*v; q2m = al**2*u**2 + be**2*v**2
q2 = al**2 + be**2; p = al*be; s1 = al + be
A5  = 20*w1**3 - 20*e*w1 - 32*q2m - 48*e + 15*w1 - 2
B5  = -15*q2*w1 + 15*p*w1 + 5*p*s1 + 8*q2 - 12*p
A23 = -60*e*w1 + 12*q2m + 72*e - 12*w1 + 2
B23 = -45*p*w1 + 15*p*s1 - 3*q2 + 18*p
Z5law  = 360*al*be*h**2/L*(A5 + h**2/L**2*B5)
Z23law = 120*al*be*h**2/L*(A23 + h**2/L**2*B23)
def laws(aval):
    u0 = -aval['s0']/aval['L']; v0 = 1 - u0
    sub = {al: aval['alpha'], be: aval['beta'], L: aval['L'], h: aval['h'], u: u0, v: v0}
    return sp.simplify(Z5law.subs(sub)), sp.simplify(Z23law.subs(sub))
Z5sum  = sp.simplify(2*laws(eA)[0] + 2*laws(eB)[0])
Z23sum = sp.simplify(2*laws(eA)[1] + 2*laws(eB)[1])
print("z5  coeff of P^(5)(1):", z5c, " sum_e Z_e^(5):", Z5sum, " diff:", sp.simplify(z5c - Z5sum))
print("z2z3 coeff of P^(5)(1):", z23c, " sum_e Z_e^(2,3):", Z23sum, " diff:", sp.simplify(z23c - Z23sum))
S = sp.simplify(2*(eA['h']**2/eA['L']) + 2*(eB['h']**2/eB['L']))
print("E1 check: z3 coeff P3(1) =", z3c, " 6S =", sp.simplify(6*S), " diff:", sp.simplify(z3c - 6*S))
print("E1 check: z2 coeff P3(1) =", z2of3, " +12 z2 S gives:", sp.simplify(z2of3 + 12*S))
pickle.dump(dict(Pk=[str(x) for x in Pk], z5c=str(z5c), z23c=str(z23c),
                 Z5sum=str(Z5sum), Z23sum=str(Z23sum)), open('steps/step_06_exact.pkl', 'wb'))
with open('steps/step_06_out.txt', 'w') as f:
    for k in range(6):
        f.write("P^(%d)(1) = %s\n" % (k, sp.collect(Pk[k], [z5, z2*z3, z3, z2])))
    f.write("z5 coeff: %s\nz2z3 coeff: %s\nZ5sum: %s\nZ23sum: %s\n" % (z5c, z23c, Z5sum, Z23sum))
print("saved steps/step_06_exact.pkl / step_06_out.txt")
