import sys, os, time, pickle
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
os.nice(10)
import sympy as sp
import tow2lib as T
L = pickle.load(open('L_parsed.pkl','rb'))
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
a5 = sp.cancel(bell(5, L))
print('a5:', round(time.time()-t0,1), 's', flush=True)
n5, d5 = sp.fraction(a5)
print('num terms:', len(sp.Poly(sp.expand(n5), T.d, T.m).terms()), ' den =', sp.factor(d5))
print(a5)
