import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
# Prepare a generic weight-5 pipeline module reused for D=29, D=21 (and rerun of D=13 if needed)
gen = r'''
# w5pipe.py -- generic weight-5 pipeline helpers for ZQ5b (enum already in w5blib)
import sys, json, time, pickle, itertools
from fractions import Fraction as Ff
import mpmath as mp
import flint
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD

def load_cands(dirn, D):
    J = json.load(open(f'{dirn}/candidates.json'))
    cands=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),D), vx=c['vx'], v1m=c['v1m'], h=c.get('h',0)) for c in J['cands']]
    return J, cands

def quotient_cols(sym_dicts, d):
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
    return cols, orbits

def kernel_and_comps(cands, d, tag, log=print):
    N=len(cands)
    t0=time.time()
    sym_dicts=[W.dense_to_dict(W.symbol_liw_dense(c['vx'],c['v1m'],d,5)) for c in cands]
    cols, orbits = quotient_cols(sym_dicts, d)
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
    log(f"[{tag}] rank {rank} kernel dim {len(K)} ({time.time()-t0:.0f}s)")
    bad=0
    for v in K:
        for j in range(len(cols)):
            if sum(v[i]*M[i][j] for i in range(N))!=0: bad+=1; break
    log(f"[{tag}] rows failing FULL check: {bad}")
    t1=time.time()
    Kred = flint.fmpz_mat(K).lll(gram='exact')
    K=[[int(Kred[i,j]) for j in range(Kred.ncols())] for i in range(Kred.nrows())]
    K=[r for r in K if any(r)]
    log(f"[{tag}] LLL-reduced kernel {len(K)} rows ({time.time()-t1:.0f}s)")
    t1=time.time()
    comps=[]
    for ri,row in enumerate(K):
        stot=W.combo_symbol(row, sym_dicts, N)
        comp=W.companion(stot, d, 5, orbits)
        comps.append({str(ms): str(cf) for ms,cf in comp.items()})
        if ri%60==0: log(f"[{tag}] comp {ri}/{len(K)} ({time.time()-t1:.0f}s)")
    log(f"[{tag}] companions done ({time.time()-t1:.0f}s)")
    return K, comps, sym_dicts, rank
'''
open('w5pipe.py','w').write(gen)
print("w5pipe.py written")
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','-c','import w5pipe; print("ok")'],capture_output=True,text=True,cwd='.')
print(r.stdout, r.stderr)
# poll d13
import time
print(open('step11b.log').read())
