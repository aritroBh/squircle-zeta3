import os, json
from fractions import Fraction as Ff
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
step = r'''
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
print(f"functional columns: {len(cols)} ({time.time()-t0:.0f}s)", flush=True)
M=[]
for s in sym_dicts:
    M.append([s.get(wi,0)-(s.get(ri,0) if ri is not None else 0) for wi,ri in cols])
print(f"M built {len(M)}x{len(cols)} ({time.time()-t0:.0f}s)", flush=True)
# mod-p RREF to find column basis
p = (1<<61)-1
Mm = flint.nmod_mat([[x % p for x in row] for row in M], p)
t1=time.time()
R = Mm.rref()
print(f"nmod rref ({time.time()-t1:.0f}s)", flush=True)
piv=[]
for i in range(R.nrows()):
    for j in range(R.ncols()):
        if R[i,j]!=0:
            piv.append(j); break
rank=len(piv)
print(f"rank(M) mod p = {rank}, pivot columns = {len(piv)} ({time.time()-t0:.0f}s)", flush=True)
Ms=[[row[j] for j in piv] for row in M]
t1=time.time()
K = W.int_left_kernel_hnf(Ms)
print(f"integer kernel dim={len(K)} (HNF {time.time()-t1:.0f}s)", flush=True)
# verify kernel rows annihilate ALL of M
bad=0
for c in K:
    for j in range(len(cols)):
        if sum(c[i]*M[i][j] for i in range(N))!=0: bad+=1; break
print(f"kernel rows failing full check: {bad}", flush=True)
# rank over Q sanity: kernel dim should be N - rank
print("N - rank =", N-rank, " kernel dim =", len(K), " match:", len(K)==N-rank)
inv = W.snf_invariants(K)
print("SNF invariants (all 1 => saturated):", sorted(set(inv)), flush=True)
json.dump(dict(K=K, names=names, d=d, rank=rank), open('out/d13/kernel5_raw.json','w'))
print("saved out/d13/kernel5_raw.json", flush=True)
'''
open('step_09b_d13_kernel.py','w').write(step)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_09b_d13_kernel.py'], capture_output=True, text=True, cwd='.', timeout=1700)
print(r.stdout); print("ERR:", r.stderr[-1500:])
