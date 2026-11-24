import pickle, itertools, time
from fractions import Fraction as Fr
import sympy as sp
B4=('a','p2','q1','q2')
Kint=pickle.load(open('q3_kernel.pkl','rb'))['Kint']
names=pickle.load(open('q3_names.pkl','rb'))['names']
names=[(k,v) if k=='li3' else ('dec',dv) for (k,v),dv in
       zip(names, [(u,v,w) for u in B4 for v in B4 for w in B4])]
pickle.dump({'names':names}, open('q3_names.pkl','wb'))
solmin=pickle.load(open('q3_solve_min.pkl','rb'))['V3p']
idx_of={}
for ci,(k,v) in enumerate(names):
    idx_of[(k,v)]=ci
x0=[Fr(0)]*336
missing=[]
for (k,v),c in solmin:
    if (k,v) in idx_of: x0[idx_of[(k,v)]]=Fr(c)
    else: missing.append((k,v))
print("missing:", missing)
D=2
x0i=[int(sp.nsimplify(D*c)) for c in x0]
print("x0' integral? ", all(sp.nsimplify(D*c).is_integer for c in x0), " max:", max(abs(e) for e in x0i))
from flint import fmpz_mat
t0=time.time()
Kb=fmpz_mat([[int(e) for e in row] for row in Kint])
Kr=fmpz_mat.lll(Kb)
Kred=[[Kr[i,j] for j in range(336)] for i in range(281)]
pickle.dump({'Kred':Kred,'x0i':x0i}, open('q3_Kred.pkl','wb'))
print("LLL done %.0fs, max entry:"%(time.time()-t0), max(max(abs(e) for e in row) for row in Kred))
