
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
    AtA = mp.matrix(len(S), len(S)); Atb = mp.matrix(len(S), 1)
    for a in range(len(S)):
        Atb[a,0] = sum(M[S[a]][i]*rv[i] for i in range(npts))
        for bb in range(len(S)):
            AtA[a,bb] = sum(M[S[a]][i]*M[S[bb]][i] for i in range(npts))
    try:
        q = list(mp.lu_solve(AtA, Atb))
    except Exception:
        return None, mp.mpf('inf')
    res = [rv[i] - sum(q[a]*M[S[a]][i] for a in range(len(S))) for i in range(npts)]
    return q, mp.sqrt(sum(r**2 for r in res))
idx = {nm: j for j, nm in enumerate(names)}
def test(label, pred):
    S = [j for nm, j in idx.items() if pred(nm)]
    q, nrm = solve_ls(S)
    print(f"{label}: {len(S)} funcs, |res| = {mp.nstr(nrm, 3)}")
    return q, S
# planner's basis: constants, z3*log, z2*Li2, z2*loglog, z2l2*log, l2^3*log
test("planner", lambda nm: nm in ('z4','z3*l2','z2*l2^2','l2^4') or nm.startswith('z3*log') or nm.startswith('z2*Li2') or nm.startswith('z2*log') or nm.startswith('z2*l2*log') or nm.startswith('l2^3*log'))
# + log2*Li3
test("planner+l2Li3", lambda nm: nm in ('z4','z3*l2','z2*l2^2','l2^4') or nm.startswith('z3*log') or nm.startswith('z2*Li2') or nm.startswith('z2*log') or nm.startswith('z2*l2*log') or nm.startswith('l2^3*log') or nm.startswith('l2*Li3'))
# everything except l2*Li2*log and l2*log^3
test("no-l2Li2log", lambda nm: not nm.startswith('l2*Li2') and not nm.startswith('l2*log'))
