import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
# add fast saturation + fast kernel to w5blib
s = open('w5blib.py').read()
fast = r'''

def saturate_rows_fast(B):
    """Saturated basis of span_Q(B) cap Z^n via HNF + flint inverse."""
    if not B: return B
    B=[_primitive(list(map(int,r))) for r in B]
    k=len(B); n=len(B[0])
    BT=[list(col) for col in zip(*B)]          # n x k
    H=flint.fmpz_mat(BT).hnf()
    H0=flint.fmpz_mat([[int(H[i,j]) for j in range(k)] for i in range(k)])
    for i in range(k,n):
        for j in range(k):
            assert int(H[i,j])==0, "hnf bottom not zero"
    H0i, den = H0.inv()                        # H0^{-1} = H0i/den
    Bm = flint.fmpz_mat(B)
    num = H0i.transpose() * Bm                 # k x n
    assert num.nrows()==k and num.ncols()==n
    rows=[]
    for i in range(k):
        row=[]
        for j in range(n):
            x = int(num[i,j])
            assert x % den == 0, "saturation not integral?!"
            row.append(x//den)
        rows.append(_primitive(row))
    return rows

def symbol_kernel_fast(M, N, ncols_quotient):
    """M: N x q integer matrix (quotient functionals applied to symbol rows).
    Returns saturated integer kernel basis of {c: c M = 0}."""
    MT=[list(col) for col in zip(*M)]
    NS, nulldim = flint.fmpz_mat(MT).nullspace()
    K=[]
    for j in range(NS.ncols()-nulldim, NS.ncols()):
        K.append([int(NS[i,j]) for i in range(NS.nrows())])
    K = [_primitive(r) for r in K]
    K = saturate_rows_fast(K)
    # verify
    for c in K:
        for j in range(ncols_quotient):
            assert sum(c[i]*M[i][j] for i in range(N))==0
    return K
'''
s += fast
open('w5blib.py','w').write(s)
print("added fast saturation + kernel")
# now finish D=13 kernel with saturation + SNF + LLL + companions
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
M=[[s.get(wi,0)-(s.get(ri,0) if ri is not None else 0) for wi,ri in cols] for s in sym_dicts]
p = (1<<61)-1
Mm = flint.nmod_mat([[x % p for x in row] for row in M], p)
R, rank_nmod = Mm.rref()
piv=[]
for i in range(R.nrows()):
    for j in range(R.ncols()):
        if R[i,j]!=0: piv.append(j); break
rank=len(piv)
print(f"rank {rank} ({time.time()-t0:.0f}s)", flush=True)
Ms=[[row[j] for j in piv] for row in M]
K = W.symbol_kernel_fast(Ms, N, len(piv))
print(f"saturated kernel dim={len(K)} ({time.time()-t0:.0f}s)", flush=True)
inv = W.snf_invariants(K)
print("SNF invariants:", sorted(set(inv)), flush=True)
# full check
bad=0
for v in K:
    for j in range(len(cols)):
        if sum(v[i]*M[i][j] for i in range(N))!=0: bad+=1; break
print("rows failing FULL check:", bad, flush=True)
# LLL reduce
t1=time.time()
Kred = flint.fmpz_mat(K).lll(gram='exact')
K=[[int(Kred[i,j]) for j in range(Kred.ncols())] for i in range(Kred.nrows())]
K=[r for r in K if any(r)]
print(f"LLL-reduced kernel ({time.time()-t1:.0f}s)", flush=True)
# companions
t1=time.time()
comps=[]
for row in K:
    stot=W.combo_symbol(row, sym_dicts, N)
    comp=W.companion(stot, d, 5, orbits)
    comps.append({str(ms): str(cf) for ms,cf in comp.items()})
print(f"companions ({time.time()-t1:.0f}s)", flush=True)
json.dump(dict(K=K, comps=comps, names=names, d=d, rank=rank), open('out/d13/kernel5.json','w'))
print("saved out/d13/kernel5.json", flush=True)
'''
open('step_09e_d13_kernel.py','w').write(step)
import subprocess
proc = subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','-u','step_09e_d13_kernel.py'],
                        stdout=open('step09e.log','w'), stderr=subprocess.STDOUT, cwd='.')
print("launched", proc.pid)
import time
time.sleep(120)
print(open('step09e.log').read())
print("running:", proc.poll() is None)
