
import sys, json, time
from fractions import Fraction as Ff
import flint
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
J = json.load(open('out/d13/candidates.json'))
cands=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),13), vx=c['vx'], v1m=c['v1m']) for c in J['cands']]
names=J['names']; d=len(names); N=len(cands)
t0=time.time()
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
K=[[int(NS[i,j]) for i in range(NS.nrows())] for j in range(nulldim)]
K=[W._primitive(r) for r in K if any(r)]
print(f"rank {rank} kernel dim {len(K)} (unsaturated) ({time.time()-t0:.0f}s)", flush=True)
# verify
bad=0
for v in K:
    for j in range(len(cols)):
        if sum(v[i]*M[i][j] for i in range(N))!=0: bad+=1; break
print("rows failing FULL check:", bad, flush=True)
# saturation index via SNF (just report)
t1=time.time()
Kred = flint.fmpz_mat(K).lll(gram='exact')
K=[[int(Kred[i,j]) for j in range(Kred.ncols())] for i in range(Kred.nrows())]
K=[r for r in K if any(r)]
print(f"LLL-reduced kernel {len(K)} rows ({time.time()-t1:.0f}s)", flush=True)
json.dump(dict(K=K, names=names, d=d, rank=rank, saturated=False), open('out/d13/kernel5.json','w'))
print("saved kernel (no comps yet)", flush=True)
# companions with progress
t1=time.time()
comps=[]
for ri,row in enumerate(K):
    stot=W.combo_symbol(row, sym_dicts, N)
    comp=W.companion(stot, d, 5, orbits)
    comps.append({str(ms): str(cf) for ms,cf in comp.items()})
    if ri%40==0: print(f"  comp {ri}/{len(K)} ({time.time()-t1:.0f}s)", flush=True)
print(f"companions done ({time.time()-t1:.0f}s)", flush=True)
json.dump(dict(K=K, comps=comps, names=names, d=d, rank=rank, saturated=False), open('out/d13/kernel5.json','w'))
print("saved out/d13/kernel5.json with companions", flush=True)
