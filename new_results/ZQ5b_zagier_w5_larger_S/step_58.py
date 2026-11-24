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
print(f"D=13: N={N} d={d}, building weight-5 symbols...", flush=True)
t0=time.time()
sym_dicts=[W.dense_to_dict(W.symbol_liw_dense(c['vx'],c['v1m'],d,5)) for c in cands]
print(f"symbols built ({time.time()-t0:.0f}s), avg support {sum(len(s) for s in sym_dicts)//N}", flush=True)
orbits = W.sym_orbits(d,5)
support=set()
for s in sym_dicts: support |= set(s.keys())
# build functional matrix
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
print(f"quotient functional columns: {len(cols)} (d^5 - C(d+4,5) = {d**5-126})", flush=True)
M=[]
for s in sym_dicts:
    M.append([s.get(wi,0)-(s.get(ri,0) if ri is not None else 0) for wi,ri in cols])
print(f"M built ({time.time()-t0:.0f}s), calling HNF left-kernel...", flush=True)
t1=time.time()
K = W.int_left_kernel_hnf(M)
print(f"kernel dim={len(K)} (HNF {time.time()-t1:.0f}s)", flush=True)
# verify kernel
for c in K:
    for j in range(len(cols)):
        assert sum(c[i]*M[i][j] for i in range(N))==0
print("kernel verified (all rows annihilate quotient)", flush=True)
inv = W.snf_invariants(K)
print("SNF invariants (all 1 => saturated):", set(inv), flush=True)
# LLL-reduce kernel
t1=time.time()
Kred = flint.fmpz_mat(K).lll(gram='exact')
K=[[int(Kred[i,j]) for j in range(Kred.ncols())] for i in range(Kred.nrows())]
K=[r for r in K if any(r)]
print(f"LLL-reduced kernel: {len(K)} rows ({time.time()-t1:.0f}s)", flush=True)
# companions
t1=time.time()
comps=[]
for row in K:
    stot=W.combo_symbol(row, sym_dicts, N)
    comp=W.companion(stot, d, 5, orbits)  # ms(tuple) -> Fraction
    comps.append({str(ms): str(cf) for ms,cf in comp.items()})
print(f"companions computed ({time.time()-t1:.0f}s)", flush=True)
json.dump(dict(K=K, comps=comps, names=names, d=d), open('out/d13/kernel5.json','w'))
print("saved out/d13/kernel5.json", flush=True)
'''
open('step_09_d13_kernel.py','w').write(step)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_09_d13_kernel.py'], capture_output=True, text=True, cwd='.', timeout=1700)
print(r.stdout); print("ERR:", r.stderr[-1500:])
