import mpmath as mp
mp.mp.dps = 200
from flint import fmpz, fmpz_mat
import sympy as sp

z3 = mp.zeta(3)
v3m = -mp.mpf('0.75')*z3
v12 = z3/8
vals = [z3, v3m, v12]
C = 10**150
rows = []
for i,v in enumerate(vals):
    r = [0]*5
    r[i] = 1
    r[3] = int(mp.nint(C*v.real))
    r[4] = int(mp.nint(C*v.imag))
    rows.append(r)
Msym = sp.Matrix(rows)
print('sympy lll...')
R = Msym.lll(delta=sp.Rational(99,100))
for r in R.tolist():
    res = sum(c*v for c,v in zip(r[:3], vals))
    print(r[:3], 'resid:', mp.nstr(abs(res),3))

print('flint variants:')
M = fmpz_mat(3, 5)
for i, r in enumerate(rows):
    for c in range(5):
        M[i,c] = fmpz(r[c])
import inspect
try:
    R2 = M.lll(delta=0.99, eta=0.51, algorithm='proved')
    for r in range(3):
        print([int(R2[r,c]) for c in range(3)])
except Exception as e:
    print('proved failed:', e)
