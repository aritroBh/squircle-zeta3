
import pickle, time
import mpmath as mp
mp.mp.dps = 90
DS = {eval(k): mp.mpf(v) for k,v in pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/series_D_v3.pkl','rb'))['DS'].items()}
B = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/invisible_series.pkl','rb'))
basis = [(nm, {eval(k): mp.mpf(v) for k, v in S_.items()}) for nm, S_ in B]
# allow e in {-1, 0}: unknowns = 2 * 149 = 298
unk = []
for nm, S_ in basis:
    unk.append((nm, -1, S_))
    unk.append((nm, 0, S_))
print(len(unk), "unknowns")
# grid: all (k,l) in D
grid = sorted(set(DS) | {k for _, e, S_ in unk for k in S_})
grid = [(k, l) for (k, l) in grid if k <= 120]
print(len(grid), "equations")
# build matrix
t0 = time.time()
rows = []
for (k, l) in grid:
    row = []
    for nm, e, S_ in unk:
        row.append(S_.get((k - e, l), mp.mpf(0)))
    rows.append(row)
b = [DS.get((k, l), mp.mpf(0)) for (k, l) in grid]
print("matrix built", time.time()-t0)
# least squares via normal equations at grid points with k <= 120
import numpy as np
# use mpmath for precision: normal equations
N = len(unk)
AtA = mp.matrix(N, N)
Atb = mp.matrix(N, 1)
for a in range(N):
    col_a = [rows[i][a] for i in range(len(grid))]
    Atb[a, 0] = sum(col_a[i]*b[i] for i in range(len(grid)))
    for bb in range(a, N):
        col_b = [rows[i][bb] for i in range(len(grid))]
        s = sum(col_a[i]*col_b[i] for i in range(len(grid)))
        AtA[a, bb] = s
        AtA[bb, a] = s
print("normal eqs", time.time()-t0)
q = mp.lu_solve(AtA, Atb)
res = []
for i in range(len(grid)):
    s = b[i] - sum(q[a,0]*rows[i][a] for a in range(N))
    res.append(abs(s))
print("max residual:", mp.nstr(max(res), 3), "time", time.time()-t0)
big = [(unk[a][0], unk[a][1], q[a,0]) for a in range(N) if abs(q[a,0]) > mp.mpf('1e-30')]
print(len(big), "nonzero")
for nm, e, v in big:
    print(f"  {mp.nstr(v, 20)}  c^{e}  {nm}")
pickle.dump({'q': [str(q[a,0]) for a in range(N)]}, open('fit_series_q.pkl','wb'))
