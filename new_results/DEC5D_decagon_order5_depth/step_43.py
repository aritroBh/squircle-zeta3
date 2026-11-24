
import pickle, time
import sympy as sp
t = sp.Symbol('t')
el0, el1, W1, W4, L2, LF = sp.symbols('el0 el1 W1 W4 L2 LF')
E5d = sp.sympify(pickle.load(open("E5_dec_subs.pkl","rb")),
                 locals={'t':t,'el0':el0,'el1':el1,'W1':W1,'W4':W4,'L2':L2,'LF':LF})
t0 = time.time()
E5e = sp.expand(E5d)
print("expand %.1fs, terms %d" % (time.time()-t0, len(sp.Add.make_args(E5e))))
pickle.dump(str(E5e), open("E5_dec_expanded.pkl","wb"))
logs = [el0, el1, W1, W4, L2, LF]
t0 = time.time()
P5 = sp.Poly(E5e, *logs)
mons = P5.terms()
print("Poly terms %.1fs, distinct monomials: %d" % (time.time()-t0, len(mons)))
from collections import Counter
degcount = Counter()
for mon, c in mons:
    degcount[sum(mon)] += 1
print("by total log-degree:", dict(sorted(degcount.items())))
pickle.dump([(mon, str(c)) for mon, c in mons], open("E5_dec_poly.pkl","wb"))
print("saved E5_dec_poly.pkl")
