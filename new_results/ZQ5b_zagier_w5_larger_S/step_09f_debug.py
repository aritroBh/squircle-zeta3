
import sys, json
from fractions import Fraction as Ff
import flint
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
J = json.load(open('out/d13/candidates.json'))
cands=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),13), vx=c['vx'], v1m=c['v1m']) for c in J['cands']]
names=J['names']; d=len(names); N=len(cands)
sym_dicts=[W.dense_to_dict(W.symbol_liw_dense(c['vx'],c['v1m'],d,5)) for c in cands]
orbits = W.sym_orbits(d,5)
support=set()
for s in sym_dicts: support |= set(s.keys())
cols=[]
for skey, perms in orbits.items():
    sup=[t for t in perms if W.tidx(t,d) in support]
    if not sup: continue
    ref=skey
    if W.tidx(ref,d) in support:
        for t in sup:
            if t!=ref: cols.append((W.tidx(t,d), W.tidx(ref,d)))
        if len(sup)<len(perms): cols.append((W.tidx(ref,d), None))
    else:
        for t in sup: cols.append((W.tidx(t,d), None))
M=[[s.get(wi,0)-(s.get(ri,0) if ri is not None else 0) for wi,ri in cols] for s in sym_dicts]
p = (1<<61)-1
Mm = flint.nmod_mat([[x % p for x in row] for row in M], p)
R, rank_nmod = Mm.rref()
piv=[]
for i in range(R.nrows()):
    for j in range(R.ncols()):
        if R[i,j]!=0: piv.append(j); break
Ms=[[row[j] for j in piv] for row in M]
MT=[list(col) for col in zip(*Ms)]
NS, nulldim = flint.fmpz_mat(MT).nullspace()
print("nullity:", nulldim, "NS dims:", NS.nrows(), NS.ncols())
# count nonzero columns
nz_first = sum(1 for j in range(nulldim) if any(int(NS[i,j]) for i in range(NS.nrows())))
nz_last = sum(1 for j in range(NS.ncols()-nulldim, NS.ncols()) if any(int(NS[i,j]) for i in range(NS.nrows())))
print("nonzero first-263:", nz_first, " nonzero last-263:", nz_last)
# rank of first/last extraction
def rank_of(jrange):
    rows=[[int(NS[i,j]) for i in range(NS.nrows())] for j in jrange if any(int(NS[i,j]) for i in range(NS.nrows()))]
    if not rows: return 0, []
    return flint.fmpz_mat(rows).rank(), rows
r1,rows1 = rank_of(range(nulldim))
r2,rows2 = rank_of(range(NS.ncols()-nulldim, NS.ncols()))
print("rank first:", r1, "rows:", len(rows1), " rank last:", r2, "rows:", len(rows2))
