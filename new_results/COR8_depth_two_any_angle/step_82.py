
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
    # least squares q minimizing |sum q_j M_j - rv| over j in S
    A = mp.matrix([[M[j][i] for j in S] for i in range(npts)])
    b = mp.matrix(rv)
    # normal equations
    AtA = mp.matrix(len(S), len(S))
    for a in range(len(S)):
        for bb in range(len(S)):
            AtA[a,bb] = sum(A[i,a]*A[i,bb] for i in range(npts))
    Atb = mp.matrix([sum(A[i,a]*b[i] for i in range(npts)) for a in range(len(S))])
    q = mp.lu_solve(AtA, Atb)
    res = [b[i] - sum(q[a]*A[i,a] for a in range(len(S))) for i in range(npts)]
    nrm = mp.sqrt(sum(r**2 for r in res))
    return q, nrm
S = []
remaining = set(range(len(names)))
q = []
for step in range(25):
    best = None
    for j in remaining:
        # cheap: correlation of M_j with current residual
        _, nrm = solve_ls(S)
        cur_res = [rv[i] - sum(q[a]*M[S[a]][i] for a in range(len(S))) for i in range(npts)] if S else rv
        num = abs(sum(M[j][i]*cur_res[i] for i in range(npts)))
        den = mp.sqrt(sum(M[j][i]**2 for i in range(npts)) * sum(r**2 for r in cur_res))
        score = num/den if den else 0
        if best is None or score > best[1]:
            best = (j, score)
    S.append(best[0])
    remaining.discard(best[0])
    q, nrm = solve_ls(S)
    print(f"step {step}: + {names[best[0]]}, residual norm {mp.nstr(nrm,3)}")
    if nrm < mp.mpf('1e-35'):
        break
pickle.dump({'S': S, 'q': [str(x) for x in q]}, open('fit_greedy.pkl','wb'))
