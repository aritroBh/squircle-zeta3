import mpmath as mp, pickle, time
import sympy as sp
mp.mp.dps = 150
r3 = mp.sqrt(3)
X = {1: 1-2*r3/3, 2: 1-r3, 5: 1-r3/3, 6: 1-r3/2}
IDX=[1,2,5,6]
l2,l3,lu,z2,z3 = mp.log(2),mp.log(3),mp.log(2-r3),mp.pi**2/6,mp.zeta(3)
vals = [z3, z2*l2, l2**3, l3**3]
N = len(vals)
C = 10**120
rows = []
for i in range(N):
    rr_ = [0]*N + [int(mp.nint(C*vals[i]))]; rr_[i] = 1
    rows.append(rr_)
M = sp.Matrix(rows)
print("input matrix:")
print(M)
BM = M.lll(delta=sp.Rational(3,4))
print("LLL output:")
print(BM)
for i in range(BM.rows):
    vc = [int(BM[i,j]) for j in range(N)]
    resid = sum(mp.mpf(vc[q])*vals[q] for q in range(N))
    print(vc, "resid:", mp.nstr(resid,4), "lastcol:", BM[i,N])
