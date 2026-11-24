import mpmath as mp
from flint import fmpz_mat
mp.mp.dps = 800
r2 = mp.sqrt(2); lam = mp.log(2); L = mp.log(1+r2)
T = (4800-3360*r2)*lam + (2880-1920*r2)*L
vals = [T, lam, r2*lam, L, r2*L]
n = len(vals)
C = mp.mpf(10)**736
rows = []
for i in range(n):
    row = [0]*(n+1); row[i] = 1
    row[n] = int(mp.nint(C*vals[i]))
    rows.append(row)
M = fmpz_mat(rows).lll()
print("first row:", [int(M[0,j]) for j in range(n)])
# norm of first row
import math
print("log10 norm row0:", math.log10(sum(int(M[0,j])**2 for j in range(n)))/2 if any(M[0,j] for j in range(n)) else 0)
