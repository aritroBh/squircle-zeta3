"""Is L^5 * coefficient a polynomial?  Expand L2^(5/2)*coef and look for negative powers of L2;
also substitute L2 -> s^2 and check the order of the pole at s = 0 numerically (random other data)."""
import pickle, sympy as sp, random, time
TAB = '<project>/research/08_theory/K1g_hexagon_order5_assembly/coef_acc_raw.pkl'
coef = pickle.load(open(TAB, 'rb'))
A0, A1, B, kap, L2, q0, V02, lA0, lA1, lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
s = sp.Symbol('s', positive=True)
random.seed(5)
worst = {}
t0 = time.time()
for i, (k, e) in enumerate(coef.items()):
    sub = {A0: sp.Rational(random.randint(1, 50), 17), A1: sp.Rational(random.randint(1, 50), 19), B: sp.Rational(random.randint(1, 50), 23),
           kap: sp.Rational(random.randint(-50, 50), 29), q0: sp.Rational(random.randint(-50, 50), 31), V02: sp.Rational(random.randint(1, 50), 37),
           lA0: sp.Rational(random.randint(-50, 50), 41), lA1: sp.Rational(random.randint(-50, 50), 43), lB: sp.Rational(random.randint(-50, 50), 47), L2: s ** 2}
    f = sp.expand(e.subs(sub))
    # lowest power of s
    pw = min(sp.Poly(sp.expand(f * s ** 20), s).monoms())[0] - 20 if f != 0 else None
    worst[k] = pw
print('time %.0f s' % (time.time() - t0))
import collections
print('lowest power of L in each coefficient (random data):', sorted(collections.Counter(worst.values()).items(), key=lambda x: (x[0] is None, x[0])))
