import mpmath as mp
from flint import fmpz_mat
mp.mp.dps = 100
# small sanity test: z3, z3*2, z3l2 with a planted relation 2*(z3*2) - 4*z3 = 0
z3 = mp.zeta(3); l2 = mp.log(2); r5 = mp.sqrt(5)
vals = [z3, r5*z3, 2*z3, r5*2*z3, z3*l2, r5*z3*l2]
names = ['z3','2z3','z3l2']
n = 3; C = mp.mpf(10)**90
rows = []
for i in range(n):
    row = [0]*(n+2); row[i] = 1
    row[n] = int(mp.nint(C*vals[2*i])); row[n+1] = int(mp.nint(C*vals[2*i+1]))
    rows.append(row)
M = fmpz_mat(rows)
print("shape:", M.nrows(), M.ncols())
Mr = M.lll()
for i in range(Mr.nrows()):
    v = [int(Mr[i,j]) for j in range(n)]
    tail = [int(Mr[i,j]) for j in range(n, n+2)]
    res1 = sum(mp.mpf(v[j])*vals[2*j] for j in range(n))
    print("row", i, "coefs", v, "tail", tail, "resid", mp.nstr(abs(res1),3))
