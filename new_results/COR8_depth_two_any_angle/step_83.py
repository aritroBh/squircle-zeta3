
import pickle
import mpmath as mp
mp.mp.dps = 50
R = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/fit_resid.pkl','rb'))
B = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/fit_basis.pkl','rb'))
rv = [mp.mpf(x) for x in R['rv']]
M = [[mp.mpf(x) for x in col] for col in B['B']]
names = B['names']
npts = len(rv)
def solve_ls(S):
    if not S:
        return [], mp.sqrt(sum(r**2 for r in rv))
    AtA = mp.matrix(len(S), len(S))
    Atb = mp.matrix(len(S), 1)
    for a in range(len(S)):
        Atb[a,0] = sum(M[S[a]][i]*rv[i] for i in range(npts))
        for bb in range(len(S)):
            AtA[a,bb] = sum(M[S[a]][i]*M[S[bb]][i] for i in range(npts))
    q = list(mp.lu_solve(AtA, Atb))
    res = [rv[i] - sum(q[a]*M[S[a]][i] for a in range(len(S))) for i in range(npts)]
    return q, mp.sqrt(sum(r**2 for r in res))
S = []
remaining = set(range(len(names)))
q = []
for step in range(30):
    q, nrm = solve_ls(S)
    cur_res = [rv[i] - sum(q[a]*M[S[a]][i] for a in range(len(S))) for i in range(npts)] if S else rv[:]
    best = None
    for j in remaining:
        num = abs(sum(M[j][i]*cur_res[i] for i in range(npts)))
        den = mp.sqrt(sum(M[j][i]**2 for i in range(npts)) * (sum(r**2 for r in cur_res) or 1))
        score = num/den
        if best is None or score > best[1]:
            best = (j, score)
    S.append(best[0])
    remaining.discard(best[0])
    q, nrm = solve_ls(S)
    print(f"step {step}: + {names[best[0]]}, |res| = {mp.nstr(nrm,3)}", flush=True)
    if nrm < mp.mpf('1e-38'):
        break
pickle.dump({'S': S, 'q': [str(x) for x in q]}, open('fit_greedy.pkl','wb'))
