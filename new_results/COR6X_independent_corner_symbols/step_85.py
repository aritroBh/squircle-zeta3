t = r'''
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
import sympy as sp, itertools

T = {(2,2,2,3): sp.Rational(5,2), (2,2,3,2): -1, (2,3,2,2): -2, (2,3,2,3): sp.Rational(1,4),
     (3,2,2,2): 4, (3,2,2,3): sp.Rational(-9,4), (3,2,3,2): sp.Rational(-1,2)}
T = {tuple(map(sp.Integer, k)): v for k, v in T.items()}

primes = [sp.Integer(2), sp.Integer(3)]
basis = list(itertools.product(primes, repeat=4))
idx = {m: i for i, m in enumerate(basis)}
rows = []
for m1 in itertools.product(primes, repeat=1):
    for m2 in itertools.product(primes, repeat=3):
        r = [0]*16
        for w in shuffle_words(m1, m2): r[idx[w]] += 1
        rows.append(r)
for m1 in itertools.product(primes, repeat=2):
    for m2 in itertools.product(primes, repeat=2):
        r = [0]*16
        for w in shuffle_words(m1, m2): r[idx[w]] += 1
        rows.append(r)
for m1 in itertools.product(primes, repeat=3):
    for m2 in itertools.product(primes, repeat=1):
        r = [0]*16
        for w in shuffle_words(m1, m2): r[idx[w]] += 1
        rows.append(r)
P = sp.Matrix(rows)
print("shuffle subspace dim:", P.rank())
dv = sp.Matrix([T.get(m, 0) for m in basis])
print("T in shuffle subspace:", P.T.rank() == P.T.row_join(dv).rank())
# also print the coordinate form if decomposable
import sympy as sp
if P.T.rank() == P.T.row_join(dv).rank():
    sol = sp.linsolve((P.T, dv))
    print("decomposition found:", sol)
'''
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=300)
print(r.stdout); print("ERR:", (r.stderr or "")[-800:])
