import pickle, sympy as sp, itertools
from fractions import Fraction as Fr
import sys; sys.path.insert(0,".")
from oct8 import vec_of_rational, vadd, vmul

sols = pickle.load(open("sols235.pkl","rb"))
letters = [2,3,5]
bidx = {t:i for i,t in enumerate(itertools.product(letters, repeat=2))}

def coord2(u, v):
    w = sp.zeros(1, 9)
    for a,ca in u.items():
        for b,cb in v.items():
            w[bidx[(a,b)]] += sp.Rational(ca*cb)
    return w

# candidate weight-2 symbol tensors (1-z) (x) z with z,1-z in V (up to torsion)
cands = list(sols) + [Fr(2),Fr(3),Fr(5),Fr(1,2),Fr(1,3),Fr(1,5),Fr(-1,2),Fr(-1,3),Fr(-1,5),Fr(-2),Fr(-4)]
W2 = []
for z in cands:
    if z == 0 or z == 1: continue
    W2.append(coord2(vec_of_rational(1-z), vec_of_rational(z)))
# symmetric products a(x)b + b(x)a (log products)
for a in letters:
    for b in letters:
        W2.append(coord2({a:Fr(1)},{b:Fr(1)}) + coord2({b:Fr(1)},{a:Fr(1)}))
M = sp.Matrix.hstack(*[w.T for w in W2])
print("rank W2 (available weight-2 symbols) in V(x)V (dim 9):", M.rank())

# also rank of {(1-z)(x)z} alone
M1 = sp.Matrix.hstack(*[coord2(vec_of_rational(1-z), vec_of_rational(z)).T for z in cands if z!=0 and z!=1])
print("rank of {(1-z)(x)z} alone:", M1.rank())

# The correct Decomp = shuffle(V (x) W2): generators a sh w2 for a in letters, w2 basis of W2
# find a basis of W2
piv = M.rref()[1]
print("W2 basis pivot cols:", piv)
