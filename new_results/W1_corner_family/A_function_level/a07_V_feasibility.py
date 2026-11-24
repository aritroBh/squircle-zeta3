"""a07_V_feasibility.py -- symbol-level feasibility of the "vanishing pair" rescue for the sqrt2-part V.

V(c) (exact HPL form, a04) has Lambda^2B_2 obstruction (1/2){c^2}^{-L}, L = (1-c)/(1+c).
N(x,y) := G(1/x,0,0,1/y;1)  has obstruction {x}^{y}  (a06).  Since {1-c^2} = -{c^2} and {(1+c)/2} = {-L},
   Theta_V(c) := (1/4)[ N(c^2,(1+c)/2) - N(1-c^2,(1+c)/2) ]
has obstruction (1/2){c^2}^{-L} and vanishes identically at c = 1/sqrt2 (where c^2 = 1-c^2).
This script (exact, symbol level):
  1. checks delta22(V - Theta_V) = 0;
  2. projects S(V - Theta_V) to the Lie coalgebra (DGR rho-operator, kills products) and solves
     rho S = sum_j beta_j rho S(Li4(f_j)) over admissible arguments f_j (f, 1-f in the alphabet);
  3. reports the point values f_j(1/sqrt2) of the arguments that are needed.
"""
import itertools, pickle, sys
from fractions import Fraction as Fr
from collections import defaultdict
import sympy as sp
from sym2 import SymbolRing, delta22, add

c = sp.Symbol('c')
R = SymbolRing([c])
for f in (c, 1 - c, 1 + c, 1 + 2*c, 1 - 2*c):
    R.letter_index(f)
print('alphabet:', R.alphabet)

# ---- V(c) symbol from the exact HPL form
Vd = pickle.load(open('a04_rfamily.pkl', 'rb'))['V']
IDX = {0: R.letter_index(c), 1: R.letter_index(1 - c), -1: R.letter_index(1 + c)}
SGN = {0: 1, 1: -1, -1: 1}
SV = defaultdict(Fr)
for (w, mono), v in Vd.items():
    if mono != () or len(w) != 4:
        continue
    sg = 1
    for a in w:
        sg *= SGN[a]
    SV[tuple(IDX[a] for a in reversed(w))] += v*sg
SV = {k: v for k, v in SV.items() if v}

def Nsym(X, Y):
    return R.gpl_symbol((1/X, 0, 0, 1/Y))
ST = add(Nsym(c**2, (1 + c)/2), Nsym(1 - c**2, (1 + c)/2), coeffs=[Fr(1, 4), Fr(-1, 4)])
print('delta22(V)       =', delta22(SV))
print('delta22(Theta_V) =', delta22(ST))
SQ = add(SV, ST, coeffs=[1, -1])
print('delta22(V - Theta_V) =', delta22(SQ), '  (alphabet size now %d)' % len(R.alphabet))
assert not delta22(SQ)

# ---- rho projection (Duhr-Gangl-Rhodes) : kills shuffle products
def rho(S):
    memo = {}
    def r(word):
        if word in memo:
            return memo[word]
        n = len(word)
        if n == 1:
            res = {word: Fr(1)}
        else:
            res = defaultdict(Fr)
            for w, v in r(word[:-1]).items():
                res[w + (word[-1],)] += Fr(n - 1, n)*v
            for w, v in r(word[1:]).items():
                res[w + (word[0],)] -= Fr(n - 1, n)*v
            res = {k: v for k, v in res.items() if v}
        memo[word] = res
        return res
    out = defaultdict(Fr)
    for w, v in S.items():
        for u, x in r(w).items():
            out[u] += v*x
    return {k: v for k, v in out.items() if v}

def li4_symbol(f):
    A = R.logs(1 - f); B = R.logs(f)
    out = defaultdict(Fr)
    for i, x in A.items():
        for j, y in B.items():
            for k, z in B.items():
                for l, t in B.items():
                    out[(i, j, k, l)] += -x*y*z*t
    return {k: v for k, v in out.items() if v}

# ---- admissible arguments in the alphabet {c,1-c,1+c,1+2c,1-2c} (+ +-2^k)  [fast exact enumerator]
from admissible import enumerate_admissible
E = int(sys.argv[1]) if len(sys.argv) > 1 else 2
roots = [Fr(0), Fr(1), Fr(-1), Fr(-1, 2), Fr(1, 2)]
adm = []
for s, k, ex in enumerate_admissible(roots, E=E, K=3):
    f = sp.Integer(s)*sp.Integer(2)**k
    for r, e in zip(roots, ex):
        f *= (c - sp.Rational(r.numerator, r.denominator))**e
    adm.append(sp.cancel(f))
print('admissible arguments found:', len(adm))
pt = 1/sp.sqrt(2)
NICE = {sp.nsimplify(v) for v in (pt, -pt, sp.Rational(1, 2), 1 - pt, 2 - sp.sqrt(2), sp.sqrt(2) - 1,
        -(sp.sqrt(2) - 1), -(sp.sqrt(2) - 1)**2, (sp.sqrt(2) - 1)**2, (sp.sqrt(2) - 1)**4, -1)}
def ptval(f):
    return sp.nsimplify(sp.radsimp(f.subs(c, pt)))
cols = []
for f in adm:
    v = ptval(f)
    rs = rho(li4_symbol(f))
    if rs:
        cols.append((f, v, rs))
target = rho(SQ)
keys = sorted(set(target) | set().union(*[set(r) for _, _, r in cols]))
M = sp.Matrix([[r.get(k, 0) for _, _, r in cols] for k in keys])
bvec = sp.Matrix([target.get(k, 0) for k in keys])
print('rho-space dimension touched:', len(keys), ' columns:', len(cols), ' rank:', M.rank())
aug = M.row_join(bvec)
print('rank with target:', aug.rank(), '(solvable iff equal)')
# solve preferring nice columns: first try nice-only
nice_cols = [i for i, (f, v, r) in enumerate(cols) if v in NICE or (v != 0 and 1/v in NICE and v < 0)]
Mn = M[:, nice_cols]
print('nice columns:', len(nice_cols), ' rank', Mn.rank(), ' rank with target', Mn.row_join(bvec).rank())
pickle.dump({'adm': [str(f) for f in adm]}, open('a07_adm.pkl', 'wb'))
for i in nice_cols:
    print('   nice:', cols[i][0], '->', cols[i][1])
