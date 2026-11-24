
import pickle, time
import mpmath as mp
mp.mp.dps = 90
DS = {eval(k): mp.mpf(v) for k,v in pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/series_D_v3.pkl','rb'))['DS'].items()}
B = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/invisible_series.pkl','rb'))['basis']
basis = [(nm, {eval(k): mp.mpf(v) for k, v in S_.items()}) for nm, S_ in B]
unk = []
for nm, S_ in basis:
    unk.append((nm, -1, S_))
    unk.append((nm, 0, S_))
grid = sorted(set(DS) | {k for _, e, S_ in unk for k in S_})
grid = [(k, l) for (k, l) in grid if k <= 120]
print(len(unk), "unknowns,", len(grid), "equations")
t0 = time.time()
N = len(unk)
G = len(grid)
# build column arrays
cols = []
for nm, e, S_ in unk:
    cols.append([S_.get((k - e, l), mp.mpf(0)) for (k, l) in grid])
b = [DS.get((k, l), mp.mpf(0)) for (k, l) in grid]
AtA = mp.matrix(N, N)
Atb = mp.matrix(N, 1)
for a in range(N):
    ca = cols[a]
    Atb[a, 0] = sum(ca[i]*b[i] for i in range(G))
    for bb in range(a, N):
        cb = cols[bb]
        s = sum(ca[i]*cb[i] for i in range(G))
        AtA[a, bb] = s
        AtA[bb, a] = s
print("normal eqs built", time.time()-t0)
q = mp.lu_solve(AtA, Atb)
maxres = max(abs(b[i] - sum(q[a,0]*cols[a][i] for a in range(N))) for i in range(G))
print("max residual:", mp.nstr(maxres, 3), "time", time.time()-t0)
big = [(unk[a][0], unk[a][1], q[a,0]) for a in range(N) if abs(q[a,0]) > mp.mpf('1e-25')]
print(len(big), "nonzero coefficients")
for nm, e, v in big:
    print(f"  {mp.nstr(v, 20)}  c^{e}  {nm}")
pickle.dump({'q': [str(q[a,0]) for a in range(N)], 'unk': [(u[0], u[1]) for u in unk]}, open('fit_series_q.pkl','wb'))
