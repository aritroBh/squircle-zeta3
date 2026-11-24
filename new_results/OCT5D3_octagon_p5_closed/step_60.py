import mpmath as mp
from flint import fmpz_mat
mp.mp.dps = 100
r2 = mp.sqrt(2); lam = mp.log(2); L = mp.log(1+r2)
T = (4800-3360*r2)*lam + (2880-1920*r2)*L
vals = [T, lam, r2*lam, L, r2*L]
n = len(vals)
C = mp.mpf(10)**90
rows = []
for i in range(n):
    row = [0]*(n+1); row[i] = 1
    row[n] = int(mp.nint(C*vals[i]))
    rows.append(row)
A = fmpz_mat(rows)
M = A.lll()
print("type:", type(M))
for i in range(n):
    print([int(M[i,j]) for j in range(n)])
