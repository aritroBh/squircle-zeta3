"""a02_involution_obstruction.py -- pull the corner weight-4 obstruction back to tau = tan(alpha/4)
and test its behaviour under the complementary-angle involution tau -> (1-tau)/(1+tau) (c <-> s),
which fixes alpha = pi/2 (tau = sqrt2-1).

tau-alphabet (function symbol, constants dropped): t0=tau, t1=1-tau, t2=1+tau, t3=1+tau^2.
c = (1-tau^2)/(1+tau^2):  [c] = t1+t2-t3,  [1-c] = 2 t0 - t3,  [1+c] = -t3
s = 2tau/(1+tau^2):       [s] = t0-t3,     [1-s] = 2 t1 - t3,  [1+s] = 2 t2 - t3
Exact rational arithmetic.  Also checks that Li4/product symbols are killed (controls).
"""
import json
from fractions import Fraction as Fr
from collections import defaultdict
from itertools import product

G1 = "research/10_frontier/G1_zeta3_any_angle/"
fit = json.load(open(G1 + "hpl_fit_k3_(0,0,0).json"))

def parse(name):
    pre, w = name.split(' H')
    return int(pre.replace('c^', '')), tuple(int(x) for x in w.strip('[]').split(','))

F = defaultdict(lambda: defaultdict(Fr))
for name, v in fit:
    if name in ('1', 'c'): continue
    e, w = parse(name)
    if len(w) == 4: F[e][w] += Fr(v)

SGN = {0: 1, 1: -1, -1: 1}
def csymbol(lin):            # symbol in c-letters 0=[c],1=[1-c],2=[1+c]
    out = defaultdict(Fr)
    IDX = {0: 0, 1: 1, -1: 2}
    for w, v in lin.items():
        sg = 1
        for a in w: sg *= SGN[a]
        out[tuple(IDX[a] for a in reversed(w))] += v*sg
    return out

def pull(S, letters):        # letters: list of dict tau-index -> coeff for c-letters 0,1,2
    out = defaultdict(Fr)
    for word, v in S.items():
        for ch in product(*(letters[i].items() for i in word)):
            co = v
            for _, d in ch: co *= d
            out[tuple(i for i, _ in ch)] += co
    return {k: x for k, x in out.items() if x}

def wedge(a, b):
    if a == b: return None, 0
    return ((a, b), 1) if a < b else ((b, a), -1)

def delta22(S):
    out = defaultdict(Fr)
    for (a, b, c, d), v in S.items():
        p, s1 = wedge(a, b); q, s2 = wedge(c, d)
        if s1 and s2:
            key, s3 = wedge(p, q)
            if s3: out[key] += v*s1*s2*s3
    return {k: x for k, x in out.items() if x}

def add(*ds, coeffs=None):
    out = defaultdict(Fr)
    for i, d in enumerate(ds):
        cf = 1 if coeffs is None else coeffs[i]
        for k, v in d.items(): out[k] += cf*v
    return {k: v for k, v in out.items() if v}

CL = [{1: 1, 2: 1, 3: -1}, {0: 2, 3: -1}, {3: -1}]      # [c],[1-c],[1+c]
SL = [{0: 1, 3: -1}, {1: 2, 3: -1}, {2: 2, 3: -1}]      # [s],[1-s],[1+s]

for e in sorted(F):
    S = csymbol(F[e])
    obc = delta22(pull(S, CL)); obs = delta22(pull(S, SL))
    print(f"e={e}: ob(F_e o c) = {obc}")
    print(f"      ob(F_e o s) = {obs}")
    print(f"      sum  = {add(obc, obs)}")
    print(f"      diff = {add(obc, obs, coeffs=[1, -1])}")

# controls: Li4 of tau-arguments and products are killed
def li4sym(one_minus, arg):
    out = defaultdict(Fr)
    for i, a in one_minus.items():
        for j, b in arg.items():
            for k, cc in arg.items():
                for l, d in arg.items():
                    out[(i, j, k, l)] += -a*b*cc*d
    return out
tests = {'tau': ({1: 1}, {0: 1}), '-tau': ({2: 1}, {0: 1}), '-tau^2': ({3: 1}, {0: 2}),
         'c': ({0: 2, 3: -1}, {1: 1, 2: 1, 3: -1}), '(1-tau)/(1+tau)': ({0: 1, 2: -1}, {1: 1, 2: -1})}
print("controls (must be {}):", [delta22(li4sym(*v)) for v in tests.values()])
