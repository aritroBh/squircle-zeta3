
import sys, json, time
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
print("rank", len(piv), flush=True)
Ms=[[row[j] for j in piv] for row in M]
MT=[list(col) for col in zip(*Ms)]   # 224 x 487
t0=time.time()
NS, nulldim = flint.fmpz_mat(MT).nullspace()
print(f"nullspace: {NS.nrows()}x{NS.ncols()} ({time.time()-t0:.0f}s)", flush=True)
# NS: right nullspace basis? rows are basis vectors of {x: MT x = 0}? check
n0 = NS.nrows(); n1 = NS.ncols(); print('nulldim', nulldim)
print("dims", n0, n1)
