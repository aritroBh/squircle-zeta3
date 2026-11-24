import pickle, time
from fractions import Fraction as Fr
import sympy as sp
B4=('a','p2','q1','q2')
names=pickle.load(open('q3_names.pkl','rb'))['names']
# check current state
print(names[272], type(names[272][1]))
decs=[(u,v,w) for u in B4 for v in B4 for w in B4]
for i in range(64):
    names[272+i]=('dec',decs[i])
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
