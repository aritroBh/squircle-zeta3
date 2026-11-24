
import pickle, itertools
from fractions import Fraction as Fr
import sys; sys.path.insert(0,'.')
from cor8b_sym import WNAME
import sympy as sp

data = pickle.load(open('pairs_delta_all.pkl','rb'))
COMP = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(1,4),(1,5),(2,4),(2,5),(3,4),(3,5),(4,5),(0,4),(0,5)]
COMP = sorted({(i,j) for i in range(6) for j in range(6) if i<j})
cidx = {k:i for i,k in enumerate(COMP)}

terms = []   # (kind, x0, y0, vector over COMP as Fractions)
for (x0,y0,d22,d31,d13) in data:
    for kind, dd in [('Li22',d22),('Li31',d31),('Li13',d13)]:
        v = [Fr(0)]*len(COMP)
        for k, val in dd.items():
            kk = (min(k),max(k))
            sgn = 1 if k==kk else -1
            v[cidx[kk]] += sgn*val
        terms.append((kind,x0,y0,tuple(v)))
print("total candidate terms:", len(terms))

T0 = [Fr(0)]*len(COMP); T1 = [Fr(0)]*len(COMP)
# target: (1-1/c)(cp^pq) + (-1-1/c)(cq^pq); cp^pq=(3,5), cq^pq=(4,5)
T0[cidx[(3,5)]] = 1; T0[cidx[(4,5)]] = -1
T1[cidx[(3,5)]] = -1; T1[cidx[(4,5)]] = -1

M = sp.Matrix([[terms[j][3][i] for j in range(len(terms))] for i in range(len(COMP))])
print("rank of all term vectors:", M.rank())
t0 = sp.Matrix(T0); t1 = sp.Matrix(T1)
aug0 = M.row_join(t0); aug1 = M.row_join(t1)
print("T0 in span?", aug0.rank()==M.rank(), " T1 in span?", aug1.rank()==M.rank())

# minimal support covering BOTH T0 and T1: try subsets of size 1,2,3,4
import itertools, time
t=time.time()
found = None
for k in range(1,5):
    sols = []
    for subset in itertools.combinations(range(len(terms)), k):
        Ms = sp.Matrix([[terms[j][3][i] for j in subset] for i in range(len(COMP))])
        if Ms.rank() < k: continue
        if Ms.row_join(t0).rank()!=Ms.rank(): continue
        if Ms.row_join(t1).rank()!=Ms.rank(): continue
        sols.append(subset)
        if len(sols)>=6: break
    print("size",k,"solutions found (capped at 6):", len(sols), "time", time.time()-t)
    if sols:
        found = sols
        break
if found:
    for subset in found:
        print("\n=== solution subset:", [(terms[j][0],terms[j][1],terms[j][2]) for j in subset])
        Ms = sp.Matrix([[terms[j][3][i] for j in subset] for i in range(len(COMP))])
        for name, tt in [('T0(coeff of c^0)',t0),('T1(coeff of c^-1)',t1)]:
            sol = Ms.LUsolve(tt) if Ms.rows>=Ms.cols else None
            try:
                sol = sp.solve_linear_system(Ms.row_join(tt), [sp.symbols(f'r{i}') for i in range(len(subset))])
            except Exception as e:
                sol = None
            print("  ", name, sol)
