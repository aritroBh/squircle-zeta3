import pickle, time
from fractions import Fraction as Fr
import sympy as sp
from q1s_symbol import tadd, t_conj
from q1s_exact import fconj, feq
B4=('a','p2','q1','q2')
sols=pickle.load(open('sunit_li3.pkl','rb'))
pairs=[]; used=set()
for i,s in enumerate(sols):
    if i in used: continue
    j=next((jj for jj,s2 in enumerate(sols) if jj!=i and feq(s2['z'],fconj(s['z']))), None)
    used.add(i)
    if j is not None: used.add(j)
    pairs.append((i,j))
def height(i):
    return sum(abs(c) for c in list(sols[i]['vec_z'].values())+list(sols[i]['vec_1z'].values()))
Gp=sorted([(height(i),(i,j)) for i,j in pairs], key=lambda x:x[0])
names=[('li3',pr) for _,pr in Gp]+[('dec',(u,v,w)) for u in B4 for v in B4 for w in B4]
print("names:", len(names))
pickle.dump({'names':names}, open('q3_names.pkl','wb'))
Kint=pickle.load(open('q3_kernel.pkl','rb'))['Kint']
solmin=pickle.load(open('q3_solve_min.pkl','rb'))['V3p']
idx_of={(k,v):ci for ci,(k,v) in enumerate(names)}
x0=[Fr(0)]*336
missing=[]
for (k,v),c in solmin:
    if (k,v) in idx_of: x0[idx_of[(k,v)]]=Fr(c)
    else: missing.append((k,v))
print("missing:", missing)
x0i=[int(sp.nsimplify(2*c)) for c in x0]
print("max x0':", max(abs(e) for e in x0i))
from flint import fmpz_mat
t0=time.time()
Kb=fmpz_mat([[int(e) for e in row] for row in Kint])
Kr=fmpz_mat.lll(Kb)
Kred=[[Kr[i,j] for j in range(336)] for i in range(281)]
pickle.dump({'Kred':Kred,'x0i':x0i}, open('q3_Kred.pkl','wb'))
print("LLL done %.0fs, max entry:"%(time.time()-t0), max(max(abs(e) for e in row) for row in Kred))
