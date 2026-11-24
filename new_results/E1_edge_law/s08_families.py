# s08: closed forms for the two fully-local families (power convention, sum_j |<v_j,x>|^p <= 1):
#  (a) every parallelogram (m = 2): on every edge alpha = beta = 1, mu = 0, no non-local normals;
#  (b) every affine-regular hexagon (m = 3): alpha = beta = 1/2 on every edge (affine invariant), the third
#      normal's term is the constant 1/2, mu0 = mu1 = mu2 = -log 2, gam = 0.
# Per-edge integrals int_0^1 E_k dt are exact functions of (L, h, u0), u0 = -s0/L.
import sympy as sp, pickle
from s03_beta_table import J, z2, z3
D = pickle.load(open('s01_E.pkl', 'rb')); E = D['E']
t, al, be, L, h, s0 = sp.symbols('t alpha beta L h s0', real=True)
l0, l1 = sp.symbols('l0 l1'); tau = sp.symbols('tau1:4'); gam = sp.symbols('gam1:4')
X, Y, u0, lg2 = sp.symbols('X Y u0 log2')
def edge_int(k, A, mu, a2):
    sub = {al: A, be: A, l0: X + mu, l1: Y + mu, s0: -u0*L}
    for kk in range(1, 4):
        sub[tau[kk-1]] = a2*mu**kk; sub[gam[kk-1]] = 0
    P = sp.Poly(sp.expand(E[k].xreplace(sub)), t, X, Y)
    return sp.factor(sp.expand(sum(c*J(n, 0, a, b) for (n, a, b), c in P.terms())))
res = {}
for name, A, mu, a2 in (('parallelogram', 1, 0, 0), ('affine_regular_hexagon', sp.Rational(1, 2), -lg2, sp.Rational(1, 2))):
    print("==", name)
    res[name] = []
    for k in range(4):
        I = edge_int(k, A, mu, a2)
        res[name].append(I)
        print(f"  int_0^1 E_{k} dt =", sp.collect(sp.expand(I), [z3, z2, lg2]))
pickle.dump(res, open('s08_families.pkl', 'wb'))
# rhombus |x/a|^p+|y/b|^p: 4 edges, L = r, h = ab/r, u0 = a^2/r^2 (V0 = (a,0) -> V1 = (0,b))
a, b = sp.symbols('a b', positive=True); r = sp.sqrt(a**2 + b**2)
P3 = sp.simplify(4*6*res['parallelogram'][3].subs({L: r, h: a*b/r, u0: a**2/r**2}))
print("rhombus P'''(1) =", sp.factor(P3), "; R2: 24a^2b^2/r^3 (z3 - pi^2/3) -> diff",
      sp.simplify(P3.subs(z2, sp.pi**2/6) - 24*a**2*b**2/r**3*(z3 - sp.pi**2/3)))
P2 = sp.simplify(4*2*res['parallelogram'][2].subs({L: r, h: a*b/r, u0: a**2/r**2}))
print("rhombus P''(1) =", sp.factor(P2.subs(z2, sp.pi**2/6)), "; R2: 4 pi^2 a^2 b^2/(3 r^3)")
with open('s08_families.txt', 'w') as f:
    for k, v in res.items():
        f.write(k + '\n' + '\n'.join(f'  int E_{i} dt = {x}' for i, x in enumerate(v)) + '\n')
