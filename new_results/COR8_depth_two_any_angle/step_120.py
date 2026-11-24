
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
N = len(unk); G = len(grid)
cols = []
for nm, e, S_ in unk:
    cols.append([S_.get((k - e, l), mp.mpf(0)) for (k, l) in grid])
colnorm = [mp.sqrt(sum(x**2 for x in col)) for col in cols]
b = [DS.get((k, l), mp.mpf(0)) for (k, l) in grid]
def solve_ls(S):
    if not S:
        return [], mp.sqrt(sum(x**2 for x in b))
    A = mp.matrix(G, len(S))
    for a, j in enumerate(S):
        for i in range(G):
            A[i, a] = cols[j][i]
    bb = mp.matrix(G, 1)
    for i in range(G): bb[i, 0] = b[i]
    q = list(mp.lu_solve(A.T * A, A.T * bb))
    q = [x[0] if hasattr(x, '__getitem__') else x for x in q]
    res = [b[i] - sum(q[a]*cols[S[a]][i] for a in range(len(S))) for i in range(G)]
    return q, mp.sqrt(sum(x**2 for x in res))
S = []
remaining = set(j for j in range(N) if colnorm[j] > mp.mpf('1e-60'))
q = []
t0 = time.time()
for step in range(60):
    q, nrm = solve_ls(S)
    cur = [b[i] - sum(q[a]*cols[S[a]][i] for a in range(len(S))) for i in range(G)] if S else b[:]
    curn = mp.sqrt(sum(x**2 for x in cur))
    best = None
    for j in remaining:
        num = abs(sum(cols[j][i]*cur[i] for i in range(G)))
        sc = num/(colnorm[j]*curn)
        if best is None or sc > best[1]:
            best = (j, sc)
    S.append(best[0]); remaining.discard(best[0])
    q, nrm = solve_ls(S)
    print(f"step {step}: + {unk[best[0]][0]} c^{unk[best[0]][1]}, |res| = {mp.nstr(nrm, 3)}", flush=True)
    if nrm < mp.mpf('1e-50'):
        break
pickle.dump({'S': S, 'q': [str(x) for x in q], 'unk': [(u[0], u[1]) for u in unk]}, open('fit_series_greedy.pkl','wb'))
print("time", time.time()-t0)
