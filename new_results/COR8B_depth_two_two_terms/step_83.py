
import pickle
import sympy as sp
c = sp.symbols('c')
vals = pickle.load(open('sunit_args_fixed.pkl','rb'))
safe, unsafe = [], []
for arg_s in vals:
    g = sp.sympify(arg_s)
    pts = [sp.N(g.subs(c, sp.Rational(k)), 20) for k in ('1/5','1/2','4/5')]
    ok = all(abs(float(p)) <= 1.0 for p in pts)
    # check boundaries too
    try:
        lim0 = sp.N(g.subs(c, sp.Rational('1/100')), 20)
        lim1 = sp.N(g.subs(c, sp.Rational('99/100')), 20)
        ok = ok and abs(float(lim0)) <= 1.0 and abs(float(lim1)) <= 1.0
    except Exception:
        ok = False
    (safe if ok else unsafe).append((arg_s, [str(p) for p in pts]))
print("SAFE args:", len(safe))
for a,p in safe: print("   ", a, p)
print("UNSAFE:", len(unsafe))
for a,p in unsafe: print("   ", a, p)
