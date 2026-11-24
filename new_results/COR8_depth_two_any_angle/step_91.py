
import pickle, time
from fractions import Fraction as Fr
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve81_v3.pkl','rb'))
li4v, prodv, Rv = D['li4v'], D['prodv'], D['Rv']
li4n = list(li4v.keys()); prodn = list(prodv.keys())
alln = [('Li4', k) for k in li4n] + [('prod', k) for k in prodn]
allv = [li4v[k] for k in li4n] + [prodv[k] for k in prodn]
M = [[allv[j][i] for j in range(len(alln))] + [Rv[-1][i], Rv[0][i], Rv[1][i]] for i in range(81)]
t0 = time.time()
basis = []
for r in M:
    for b, pw in basis:
        if r[pw]:
            f = r[pw]; r = [x - f*y for x, y in zip(r, b)]
    nz = next((j for j in range(len(alln)) if r[j]), None)
    if nz is None: continue
    f = r[nz]; basis.append(([x/f for x in r], nz))
sols = {}
for k, E in enumerate((-1,0,1)):
    rhs = len(alln) + k
    sol = [Fr(0)]*len(alln)
    for b, pw in reversed(basis):
        sol[pw] = b[rhs] - sum(b[j]*sol[j] for j in range(pw+1, len(alln)) if sol[j])
    sols[E] = sol
for E in (-1,0,1):
    acc = [Fr(0)]*81
    for j in range(len(alln)):
        if sols[E][j]:
            acc = [x + sols[E][j]*y for x, y in zip(acc, allv[j])]
    print(f"E={E} residual:", max(abs(x-y) for x,y in zip(Rv[E], acc)), end='')
    used = [(alln[j], sols[E][j]) for j in range(len(alln)) if sols[E][j]]
    print(f"  {len(used)} terms ({sum(1 for x in used if x[0][0]=='Li4')} Li4)")
pickle.dump({'sols':sols,'alln':alln}, open('full_sols_v3.pkl','wb'))
for j, (kind, nm) in enumerate(alln):
    if any(sols[E][j] for E in (-1,0,1)) and kind=='Li4':
        print(f"  a[{nm}](c) =", {E: str(sols[E][j]) for E in (-1,0,1) if sols[E][j]})
