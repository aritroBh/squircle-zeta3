
import pickle, time
from fractions import Fraction as F
import mzv_exact as M, engine as E
def table(W):
    if W not in E._TABLES:
        E._TABLES[W]=M.DSTableFlint(W) if W>=11 else M.DSTable(W)
    return E._TABLES[W]
E.table=table
result = pickle.load(open('partB_planecoords.pkl','rb'))

def peval2(poly, x, y):
    return sum(c*x**e1*y**e2 for (e1,e2),c in poly.items())

t0=time.time(); out=[]
mism=0
for plane in ['Pi1','Pi2']:
    for s in range(1,21):
        a=(s,s,1,1) if plane=='Pi1' else (s,1,s,1)
        r=E.reduced_coeffs(a,12)[12]
        for b in range(12):
            mine=peval2(result[plane][12][b], F(s), F(1))
            if mine!=r[b]:
                mism+=1
                print("MISMATCH", plane, s, b, r[b], mine)
        out.append((plane,s))
    print(f"{plane}: 20 points compared, t={time.time()-t0:.1f}s")
print("W=12 mismatches:", mism)
