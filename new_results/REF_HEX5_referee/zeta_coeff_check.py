"""thm:hexagon5: in the 307-monomial form, the standalone zeta(5) and zeta(2)zeta(3) coefficients
(sum over x-powers of the monomials with no HPL factor) equal the printed Z_e^{(5)}, Z_e^{(2,3)}
(checked symbolically after A0 = 1 - B(1+kap), with the referee's derived law).  Also lists every
monomial that contains zeta(5) or zeta(2)zeta(3) or zeta(3) times an HPL."""
import pickle, sympy as sp, random
from local_law import printed, A0 as sA0, A1 as sA1, L as sL, L0 as sL0, h as sh
TAB = '<project>/research/08_theory/K1g_hexagon_order5_assembly/coef_acc_raw.pkl'
coef = pickle.load(open(TAB, 'rb'))
A0, A1, B, kap, L2, q0, V02, lA0, lA1, lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
z5keys = [k for k in coef if k[0] == () and k[2] == (0, 0, 1)]
z23keys = [k for k in coef if k[0] == () and k[2] == (1, 1, 0)]
print('standalone zeta5 monomials (x-powers):', [k[1] for k in z5keys])
print('standalone zeta2zeta3 monomials (x-powers):', [k[1] for k in z23keys])
other = sorted(set((k[0] != (), k[2]) for k in coef))
print('zeta patterns present (has HPL factor, (z2,z3,z5) powers):', other)
P5, P23 = printed()
x = -kap
random.seed(1)
for trial in range(4):
    # random valid rational edge data
    Bv = sp.Rational(random.randint(1, 9), 10); kv = sp.Rational(random.randint(-9, 30), 10)
    A0v = 1 - Bv * (1 + kv); A1v = 1 - Bv
    L2v = sp.Rational(random.randint(5, 40), 7); q0v = sp.Rational(random.randint(-20, 20), 9)
    V02v = q0v ** 2 / L2v + sp.Rational(random.randint(1, 30), 11)
    sub = {A0: A0v, A1: A1v, B: Bv, kap: kv, L2: L2v, q0: q0v, V02: V02v,
           lA0: sp.Symbol('u0'), lA1: sp.Symbol('u1'), lB: sp.Symbol('u2')}
    c5 = sp.simplify(120 * sum(coef[k].subs(sub) * (-kv) ** k[1] for k in z5keys))
    c23 = sp.simplify(120 * sum(coef[k].subs(sub) * (-kv) ** k[1] for k in z23keys))
    Lv = sp.sqrt(L2v)
    lsub = {sA0: A0v, sA1: A1v, sL: Lv, sL0: -q0v / Lv, sh: sp.sqrt(V02v - q0v ** 2 / L2v)}
    print('trial %d: B=%s kap=%s  [z5] table - printed = %s ; [z2z3] table - printed = %s' % (
        trial, Bv, kv, sp.simplify(c5 - P5.subs(lsub)), sp.simplify(c23 - P23.subs(lsub))))
