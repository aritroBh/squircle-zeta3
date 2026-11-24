"""a06_find_N.py -- find a 2-variable GPL N(a,b) (letters {0,a,b}, endpoint 1) whose Lambda^2B_2 obstruction
is a single wedge {X}^{Y} of B_2-elements X ~ a, Y ~ b, i.e. a weight-4 function realising {x}^{y}."""
import itertools, sympy as sp
from sym2 import SymbolRing, delta22, b2, wedge2, add
a, b = sp.symbols('a b')
R = SymbolRing([a, b])
for f in (a, b, 1 - a, 1 - b, a - b):
    R.letter_index(f)
print('alphabet:', R.alphabet)
targets = {'{1/a}^{1/b}': wedge2(b2(R, 1/a), b2(R, 1/b)),
           '{a/b}^{1/b}': wedge2(b2(R, a/b), b2(R, 1/b)),
           '{1/a}^{a/b}': wedge2(b2(R, 1/a), b2(R, a/b))}
for k, v in targets.items():
    print('target', k, v)
res = {}
for w in itertools.product([0, a, b], repeat=4):
    if w[-1] == 0 or w[0] == 0 and all(x == 0 for x in w):
        continue
    d = delta22(R.gpl_symbol(w))
    if d:
        res[w] = d
print(len(res), 'words with nonzero delta22')
for w in [(0, 0, a, b), (0, a, 0, b), (a, 0, 0, b), (0, 0, b, a), (0, b, 0, a)]:
    print(w, res.get(w))
