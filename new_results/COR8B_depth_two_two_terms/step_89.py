
import pickle
from fractions import Fraction as Fr
dataH = pickle.load(open('hpl_delta_all.pkl','rb'))
ALLCOMP = sorted({(i,j) for i in range(6) for j in range(6) if i<j})
cidx = {k:i for i,k in enumerate(ALLCOMP)}
def vec_of(dd):
    v = [Fr(0)]*len(ALLCOMP)
    for k,val in dd.items():
        kk=(min(k),max(k)); sgn=1 if k==kk else -1
        v[cidx[kk]] += sgn*val
    return tuple(v)
# -R components (from DT_A_residual_2kept.json): (mp^cp): -(1/c-1), (mp^cq): -(1/c+1), (mq^cp): -(-1/c+1), (mq^cq): -(-1/c-1)
# i.e. -R = (mp^cp)(1-1/c) + (mp^cq)(-1-1/c) + (mq^cp)(1-1/c) + (mq^cq)(1+1/c)
R0 = [Fr(0)]*len(ALLCOMP); R1 = [Fr(0)]*len(ALLCOMP)
# R = R_{-1} * (1/c) + R_0 * 1 ; -R = -(R_{-1}/c + R_0)
# From earlier: resid d22(DT_A)-d22(Omega): (mp^cp){-1:1,0:-1}, (mp^cq){-1:1,0:1}, (mq^cp){-1:-1,0:1}, (mq^cq){-1:-1,0:-1}
Rm1 = [Fr(0)]*len(ALLCOMP); R0 = [Fr(0)]*len(ALLCOMP)
Rm1[cidx[(1,3)]] = 1; R0[cidx[(1,3)]] = -1
Rm1[cidx[(1,4)]] = 1; R0[cidx[(1,4)]] = 1
Rm1[cidx[(2,3)]] = -1; R0[cidx[(2,3)]] = 1
Rm1[cidx[(2,4)]] = -1; R0[cidx[(2,4)]] = -1
targ0 = tuple(-x for x in R0); targ1 = tuple(-x for x in Rm1)
vecs, labels = [], []
for (g,w,dd) in dataH:
    vecs.append(vec_of(dd)); labels.append((g,w))
def rank(rows):
    rows = [list(r) for r in rows if any(r)]
    if not rows: return 0
    m,n=len(rows),len(rows[0]); r=0
    for col in range(n):
        piv=None
        for i in range(r,m):
            if rows[i][col]: piv=i;break
        if piv is None: continue
        rows[r],rows[piv]=rows[piv],rows[r]
        pv=rows[r][col]
        for i in range(m):
            if i!=r and rows[i][col]:
                f=rows[i][col]/pv
                rows[i]=[a-f*b for a,b in zip(rows[i],rows[r])]
        r+=1
        if r==m: break
    return r
# find solutions for targ0 and targ1 with shared support (small): try singles, pairs, triples
import itertools, time
t0=time.time()
def solve_with_basis(idxset, targ):
    # exact solve via Fraction GE
    rows = [vecs[i] for i in idxset]
    # transpose: solve M x = targ
    import sympy as sp
    M = sp.Matrix([[vecs[j][i] for j in idxset] for i in range(len(ALLCOMP))])
    tv = sp.Matrix(list(targ))
    try:
        sol = M.LUsolve(tv)
        return [Fr(sol[i]) for i in range(len(idxset))]
    except Exception:
        return None
found = None
for k in range(1,5):
    cnt=0
    for subset in itertools.combinations(range(len(vecs)), k):
        # quick rank screen
        sub = [vecs[i] for i in subset]
        if rank(sub) < k: continue
        if rank(sub+[targ0]) != rank(sub): continue
        if rank(sub+[targ1]) != rank(sub): continue
        s0 = solve_with_basis(subset, targ0)
        s1 = solve_with_basis(subset, targ1)
        if s0 and s1:
            found = (subset, s0, s1); break
    print("size",k,"done, found:", found is not None, "time", time.time()-t0)
    if found: break
if found:
    subset, s0, s1 = found
    print("\nX = combination of HPLs with d22(X) = -R:")
    for i, a0, a1 in zip(subset, s0, s1):
        g, w = labels[i]
        print(f"   H_{{{w}}}({g})   coeff: (c-part {a1}) + (const {a0})")
pickle.dump(found, open('X_hpl_solution.pkl','wb'))
