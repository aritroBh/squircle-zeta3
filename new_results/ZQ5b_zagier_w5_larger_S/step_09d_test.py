
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
rank=len(piv)
Ms=[[row[j] for j in piv] for row in M]
MT=[list(col) for col in zip(*Ms)]
NS, nulldim = flint.fmpz_mat(MT).nullspace()
# try LAST nulldim columns as basis
Kcols=[]
for j in range(NS.ncols()-nulldim, NS.ncols()):
    v=[int(NS[i,j]) for i in range(NS.nrows())]
    Kcols.append(v)
# check they annihilate MT
ok=0
for v in Kcols:
    good=all(sum(MT[r][i]*v[i] for i in range(N))==0 for r in range(len(MT)))
    if good: ok+=1
print(f"last-cols convention: {ok}/{len(Kcols)} annihilate MT")
# also first nulldim cols if needed
Kcols2=[]
for j in range(nulldim):
    v=[int(NS[i,j]) for i in range(NS.nrows())]
    Kcols2.append(v)
ok2=sum(1 for v in Kcols2 if all(sum(MT[r][i]*v[i] for i in range(N))==0 for r in range(len(MT))))
print(f"first-cols convention: {ok2}/{len(Kcols2)} annihilate MT")
# pick the working one
K = Kcols if ok==len(Kcols) else Kcols2
# full check on all 2999 cols + SNF
bad=0
for v in K:
    for j in range(len(cols)):
        if sum(v[i]*M[i][j] for i in range(N))!=0: bad+=1; break
print("rows failing FULL M check:", bad, " of", len(K))
inv = W.snf_invariants(K)
print("SNF invariants:", sorted(set(inv)))
print("rank check: N - rank =", N-rank, " nulldim =", nulldim)
