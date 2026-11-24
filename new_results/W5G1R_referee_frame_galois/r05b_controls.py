"""r05b: sensitivity controls for r05(c): X.pi and (e sh Z).pi are nonzero mod p, so omitting either breaks r05(c)."""
from rcommon import *
import pickle, random
L = 'r05b.log'
W5 = words(5); IDX = {w: i for i, w in enumerate(W5)}
X = pickle.load(open(HERE + '/r05_X.pkl', 'rb'))['X']; Z = pickle.load(open(HERE + '/r04_Z.pkl', 'rb'))
p = 1000000007; rnd = random.Random(12345 + p % 1000); piv = [rnd.randrange(p) for _ in W5]
def m(fr): fr = F(fr); return fr.numerator % p * pow(fr.denominator % p, -1, p) % p
def proj(d): return sum(m(c) * piv[IDX[w]] for w, c in d.items()) % p
for c in ('rat', 'sq'):
    log(f'{c}: X.pi = {proj(X[c])}, (e sh Z).pi = {proj(shv({(2,): 1, (6,): -1}, Z[c]))} (mod {p}; same pi as r05)', L)
