import sys, os, time, pickle
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
os.nice(10)
import sympy as sp
import tow2lib as T
Ls = pickle.load(open('L_simplified.pkl','rb'))
d, m = T.d, T.m
L = {r: sp.sympify(Ls[r], locals={'d': d, 'm': m, 'I': sp.I}) for r in range(1, 9)}
import random
random.seed(1)
for r in [5, 6]:
    for _ in range(3):
        dv, mv = random.randint(6, 30), random.randint(1, 5)
        diff = sp.expand(L[r].subs({d: dv, m: mv}) - T.L(r).subs({d: dv, m: mv}))
        assert diff == 0, (r, dv, mv, diff)
print('L[5], L[6] match tow2lib at random points: OK')
pickle.dump(L, open('L_parsed.pkl','wb'))
# Now Bell-combine for a_5, a_6 using parsed L (cheap)
from collections import Counter
def parts(rem, minpart, cur):
    if rem == 0:
        yield tuple(sorted(cur)); return
    for p in range(minpart, rem+1):
        yield from parts(rem-p, p, cur+[p])
def bell(delta, L):
    tot = sp.Integer(0)
    for part in parts(delta, 1, []):
        C = Counter(part)
        coef = sp.Integer(1)
        for sz, mult in C.items():
            coef *= (L[sz]**mult) / sp.factorial(mult)
        tot += coef
    return tot
t0 = time.time()
a5 = bell(5, L)
print('a5 combined', round(time.time()-t0,1), flush=True)
a5s = sp.cancel(a5)
print('a5 cancelled', round(time.time()-t0,1), flush=True)
print('a5 =', a5s)
