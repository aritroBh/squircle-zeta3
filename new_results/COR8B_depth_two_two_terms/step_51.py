
import sys, pickle, itertools, time
sys.path.insert(0,'.')
from fractions import Fraction as Fr
br = pickle.load(open('basis_red81.pkl','rb'))
basis_red = [(n, r) for n,r in br]
piv_rows = {}
for name, row in basis_red:
    p = min(row)
    piv_rows[p] = (row, name)
T = pickle.load(open('targets81.pkl','rb'))
def reduce_track(tv):
    v = {i: x for i,x in enumerate(tv) if x}
    coef = {}
    while v:
        i = min(v)
        if i not in piv_rows: return v, coef
        row, name = piv_rows[i]
        f = v[i]/row[i]
        coef[name] = coef.get(name, Fr(0)) + f
        for j,x in row.items():
            if j==i: continue
            nv = v.get(j, Fr(0)) - f*x
            if nv: v[j] = nv
            elif j in v: del v[j]
        del v[i]
    return v, coef
for e in (-1,0,1):
    resid, coef = reduce_track(list(T[e]))
    print("e=",e,"residual coords:", len(resid))
    if resid:
        for i,x in sorted(resid.items())[:20]:
            print(f"    idx {i} (letters {i//27},{(i//9)%3},{(i//3)%3},{i%3}): {x}")
# also print the e=1 solution
resid1, coef1 = reduce_track(list(T[1]))
print("\ne=1 solution (39 terms):")
for n,c0 in sorted(coef1.items(), key=lambda kv: str(kv[0])):
    print("   ", c0, n)
