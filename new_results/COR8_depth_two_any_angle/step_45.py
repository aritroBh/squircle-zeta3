
import pickle
from fractions import Fraction as Fr
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve81.pkl','rb'))
li4v, prodv, Rv = D['li4v'], D['prodv'], D['Rv']
li4n = list(li4v.keys()); prodn = list(prodv.keys())
alln = [('Li4', k) for k in li4n] + [('prod', k) for k in prodn]
allv = [li4v[k] for k in li4n] + [prodv[k] for k in prodn]
print(len(alln), "columns")
# exact rref of the 81 x 1149 matrix (columns ordered as given), for all three RHS at once
M = [[allv[j][i] for j in range(len(alln))] + [Rv[-1][i], Rv[0][i], Rv[1][i]] for i in range(81)]
# forward elimination to row echelon (not reduced), tracking pivots
import time
t0 = time.time()
basis = []  # (row, pivot)
for r in M:
    for b, pw in basis:
        if r[pw]:
            f = r[pw]
            r = [x - f*y for x, y in zip(r, b)]
    nz = next((j for j in range(len(alln)) if r[j]), None)
    if nz is None:
        continue
    f = r[nz]
    basis.append(([x/f for x in r], nz))
pivots = [pw for _, pw in basis]
print("rank:", len(basis), "time", time.time()-t0)
# back-substitute: for each RHS, particular solution with free=0, and check consistency
sols = {}
for k, E in enumerate((-1,0,1)):
    rhs_idx = len(alln) + k
    # check consistency: any row in M (after elim) with all-zero in cols but nonzero in rhs?
    sol = [Fr(0)]*len(alln)
    for b, pw in reversed(list(zip(basis, pivots))):
        val = b[rhs_idx] - sum(b[j]*sol[j] for j in range(pw+1, len(alln)) if sol[j])
        sol[pw] = val
    sols[E] = sol
    used = [(alln[j], sol[j]) for j in range(len(alln)) if sol[j]]
    nli4 = sum(1 for x in used if x[0][0]=='Li4')
    print(f"E={E}: {len(used)} terms used ({nli4} Li4, {len(used)-nli4} products)")
    print("   Li4:", [(x[0][1], str(x[1])) for x in used if x[0][0]=='Li4'])
pickle.dump({'sols':sols,'alln':alln}, open('full_sols_raw.pkl','wb'))
