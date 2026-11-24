"""classical_nullspace.py -- find the functional relations (null vectors) among the 45 classical basis
functions of classical_dict.py, so that a dependent function can be dropped."""
from fractions import Fraction as Fr
import mpmath as mp
import classical_dict as cd
mp.mp.dps = 130
B = cd.basis_funcs()
xs = [mp.mpf(3)/100 + mp.mpf(i)*mp.mpf(42)/100/69 for i in range(70)]
A = mp.matrix(len(xs), len(B))
for i, x in enumerate(xs):
    for j, (_, fn) in enumerate(B):
        A[i, j] = fn(x)
U, S, V = mp.svd_r(A)
S = list(S)
small = [k for k, s in enumerate(S) if s < mp.mpf(10)**-60]
print('singular values below 1e-60:', [mp.nstr(S[k], 3) for k in small])
for k in small:
    v = [V[k, j] for j in range(len(B))]
    piv = max(range(len(B)), key=lambda j: abs(v[j]))
    v = [x/v[piv] for x in v]
    terms = [(B[j][0], Fr(str(mp.nstr(v[j], 50))).limit_denominator(1000)) for j in range(len(B)) if abs(v[j]) > mp.mpf(10)**-50]
    print('RELATION:', ' + '.join(f'({c})*{n}' for n, c in terms))
