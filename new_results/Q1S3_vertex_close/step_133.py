import pickle, itertools, time
from fractions import Fraction as Fr
import sympy as sp
from q1s_symbol import tadd, t_conj
B4=('a','p2','q1','q2')
BIDX=list(itertools.product(B4,repeat=3))
def coord(T):
    v=sp.zeros(64,1)
    for k,c in T.items(): v[BIDX.index(k)]=sp.Rational(c)
    return v
Kint=pickle.load(open('q3_kernel.pkl','rb'))['Kint']
names=pickle.load(open('q3_names.pkl','rb'))['names']
t0=time.time()
# particular solution: from q3_solve_min (55-col solve), embedded
solmin=pickle.load(open('q3_solve_min.pkl','rb'))['V3p']
sel=pickle.load(open('q3_sel.pkl','rb'))['sel']
# map sel entries to column indices of the 336-list
idx_of={}
for ci,(k,v) in enumerate(names):
    idx_of.setdefault((k,v), ci)
x0=[Fr(0)]*336
for (k,v),c in solmin:
    x0[idx_of[(k,v)]]=Fr(c)
D=2
x0i=[int(sp.nsimplify(D*c)) for c in x0]
print("particular sol x0' integral, max:", max(abs(e) for e in x0i))
# flint LLL on kernel basis
from flint import fmpz_mat
Kb=fmpz_mat([[int(e) for e in row] for row in Kint])
print("Kb:", Kb.nrows(), Kb.ncols(), " %.0fs"%(time.time()-t0))
t1=time.time()
Kr=fmpz_mat.lll(Kb)
print("LLL done %.0fs"%(time.time()-t1))
Kred=[[Kr[i,j] for j in range(336)] for i in range(281)]
pickle.dump({'Kred':Kred}, open('q3_Kred.pkl','wb'))
print("max entry after LLL:", max(max(abs(e) for e in row) for row in Kred))
