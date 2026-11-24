
import pickle, time
from fractions import Fraction as Fr
data = pickle.load(open('pairs_delta_all.pkl','rb'))
COMP = sorted({(i,j) for i in range(6) for j in range(6) if i<j})
cidx = {k:i for i,k in enumerate(COMP)}
USE = [(1,3),(1,4),(2,3),(2,4),(3,4),(3,5),(4,5)]  # the 7 components that appear
uidx = {k:i for i,k in enumerate(USE)}
terms = []
for (x0,y0,d22,d31,d13) in data:
    for kind, dd in [('Li22',d22),('Li31',d31),('Li13',d13)]:
        v = [Fr(0)]*7
        for k,val in dd.items():
            kk=(min(k),max(k)); sgn=1 if k==kk else -1
            if kk in uidx: v[uidx[kk]] += sgn*val
        terms.append((kind,x0,y0,tuple(v)))
print("terms:", len(terms))
# unique vectors with their term labels (dedupe by vector)
from collections import defaultdict
byvec = defaultdict(list)
for t in terms: byvec[t[3]].append(t[:3])
vecs = list(byvec.keys())
print("distinct vectors:", len(vecs))
T0 = tuple(Fr(1) if u==(3,5) else Fr(-1) if u==(4,5) else Fr(0) for u in USE)
T1 = tuple(Fr(-1) if u in ((3,5),(4,5)) else Fr(0) for u in USE)

def rank(rows):
    # exact rank via Fraction GE
    rows = [list(r) for r in rows if any(r)]
    if not rows: return 0
    m, n = len(rows), len(rows[0]); r = 0
    for col in range(n):
        piv = None
        for i in range(r, m):
            if rows[i][col]: piv = i; break
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        pv = rows[r][col]
        for i in range(m):
            if i != r and rows[i][col]:
                f = rows[i][col]/pv
                rows[i] = [a - f*b for a,b in zip(rows[i], rows[r])]
        r += 1
        if r == m: break
    return r

t=time.time()
R = rank(vecs)
print("rank of distinct vectors:", R, "time", time.time()-t)
print("T0 in span:", rank(vecs+[T0])==R, " T1 in span:", rank(vecs+[T1])==R)

# terms lying in plane P=span(T0,T1)
inP = [i for i,v in enumerate(vecs) if rank([T0,T1,v])==2]
print("vectors lying in target plane P:", len(inP))
for i in inP[:10]: print("   ", vecs[i], byvec[vecs[i]][:3])
