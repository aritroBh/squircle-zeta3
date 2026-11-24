
import pickle, time
from fractions import Fraction as Fr
from collections import defaultdict
data = pickle.load(open('pairs_delta_all_v2.pkl','rb'))
USE = [(1,3),(1,4),(2,3),(2,4),(3,4),(3,5),(4,5)]  # mp^cp,mp^cq,mq^cp,mq^cq,cp^cq,cp^pq,cq^pq
uidx = {k:i for i,k in enumerate(USE)}
terms = []
for (x0,y0,d22,d31,d13) in data:
    for kind, dd in [('Li22',d22),('Li31',d31),('Li13',d13)]:
        v = [Fr(0)]*7
        for k,val in dd.items():
            kk=(min(k),max(k)); sgn=1 if k==kk else -1
            if kk in uidx: v[uidx[kk]] += sgn*val
        terms.append((kind,x0,y0,tuple(v)))
byvec = defaultdict(list)
for t in terms: byvec[t[3]].append((t[0],t[1],t[2]))
vecs = list(byvec.keys())
print("terms:", len(terms), "distinct:", len(vecs))
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
T0 = tuple(Fr(1) if u==(3,5) else Fr(-1) if u==(4,5) else Fr(0) for u in USE)
T1 = tuple(Fr(-1) if u in ((3,5),(4,5)) else Fr(0) for u in USE)
R = rank(vecs)
print("span rank:", R, " T0 in span:", rank(vecs+[T0])==R, " T1 in span:", rank(vecs+[T1])==R)
print("\ndistinct vectors and example terms:")
for v in vecs:
    ex = byvec[v][0]
    print("  ", tuple(int(x) for x in v), "  e.g.", ex, f"(+{len(byvec[v])-1} more)")
