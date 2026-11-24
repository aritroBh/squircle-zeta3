
import pickle
from fractions import Fraction as Fr
from collections import defaultdict
data = pickle.load(open('pairs_delta_all_v2.pkl','rb'))
USE = [(1,3),(1,4),(2,3),(2,4),(3,4),(3,5),(4,5)]
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
vecs = [v for v in byvec if any(v)]

def rref(rows, ncols):
    rows = [list(r) for r in rows]
    m = len(rows); r = 0; pivs = []
    for col in range(ncols):
        piv = None
        for i in range(r, m):
            if rows[i][col]: piv = i; break
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        pv = rows[r][col]
        rows[r] = [a/pv for a in rows[r]]
        for i in range(m):
            if i != r and rows[i][col]:
                f = rows[i][col]
                rows[i] = [a - f*b for a,b in zip(rows[i], rows[r])]
        pivs.append(col); r += 1
        if r == m: break
    return rows, pivs

# span basis
rows, pivs = rref([list(v) for v in vecs], 7)
basis = [tuple(rows[i]) for i in range(len(pivs))]
print("span basis (RREF):", basis)
# annihilator: solve A f = 0 where A rows = basis
A = [list(b) for b in basis]
# f in Q^7, A f = 0 -> compute null space
rr, pv = rref([list(b) for b in basis], 7)
free = [j for j in range(7) if j not in pv]
nulls = []
for j in free:
    f = [Fr(0)]*7; f[j] = Fr(1)
    for i,pc in enumerate(pv):
        f[pc] = -rr[i][j]
    nulls.append(tuple(f))
print("\nannihilator basis (4 functionals):")
for f in nulls: print("  ", tuple(int(x) if x==int(x) else x for x in f))
T0 = tuple(Fr(1) if u==(3,5) else Fr(-1) if u==(4,5) else Fr(0) for u in USE)
T1 = tuple(Fr(-1) if u in ((3,5),(4,5)) else Fr(0) for u in USE)
print("\nobstruction values:")
for f in nulls:
    print("  f(T0) =", sum(a*b for a,b in zip(f,T0)), " f(T1) =", sum(a*b for a,b in zip(f,T1)))
