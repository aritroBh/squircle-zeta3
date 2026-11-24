import sys, os, time, pickle
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
os.nice(10)
import sympy as sp
import tow2lib as T

Ls = pickle.load(open('L_simplified.pkl','rb'))
d, m = T.d, T.m
L = {}
for r in range(1, 9):
    L[r] = sp.sympify(Ls[r], locals={'d': d, 'm': m, 'I': sp.I})
print('L parsed')
# sanity check L[5], L[6] against tow2lib's L(r) at random rational points
import random
random.seed(1)
for r in [5, 6]:
    for _ in range(3):
        dv, mv = random.randint(6, 30), random.randint(1, 5)
        v1 = sp.Rational(L[r].subs({d: dv, m: mv}))
        v2 = sp.Rational(sp.expand(T.L(r)).subs({d: dv, m: mv}))
        assert v1 == v2, (r, dv, mv, v1-v2)
print('L[5], L[6] match tow2lib at random points: OK')
pickle.dump(L, open('L_parsed.pkl','wb'))
