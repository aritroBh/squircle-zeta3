import os
cwd = os.getcwd()
kscript = r'''import sys, json, time
from fractions import Fraction as Ff
import flint
sys.path.insert(0,'.')
import w5blib as W, w5pipe
from w5blib import QD
J, allc = w5pipe.load_cands('out/d29/candidates_S7_B5.json'.replace('candidates_S7_B5.json','candidates_S7_B5.json'), 29) if False else (None, None)
J = json.load(open('out/d29/candidates_S7_B5.json'))
allc=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),29), vx=c['vx'], v1m=c['v1m'], h=c['h']) for c in J['cands']]
names=J['names']; d=len(names)
print(f"d29-S7: total {len(allc)}, d={d}", flush=True)
def kdim(cands):
    N=len(cands)
    sym_dicts=[W.dense_to_dict(W.symbol_liw_dense(c['vx'],c['v1m'],d,5)) for c in cands]
    cols, orbits = w5pipe.quotient_cols(sym_dicts, d)
    M=[[s.get(wi,0)-(s.get(ri,0) if ri is not None else 0) for wi,ri in cols] for s in sym_dicts]
    p = (1<<61)-1
    Mm = flint.nmod_mat([[x % p for x in row] for row in M], p)
    R, rank_nmod = Mm.rref()
    rank=0
    for i in range(R.nrows()):
        if any(R[i,j] for j in range(R.ncols())): rank+=1
    return N - rank, rank
for hmax in (1,2,3,4,5):
    cands=[c for c in allc if c['h']<=hmax]
    t0=time.time()
    kd, rank = kdim(cands)
    print(f"  h<={hmax}: N={len(cands)}, kernel dim={kd}, rank={rank} ({time.time()-t0:.0f}s)", flush=True)
'''
open('step_11_d29_kdims.py','w').write(kscript)
import subprocess, time
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_11_d29_kdims.py'],
                   capture_output=True, text=True, cwd=cwd, timeout=2400)
print("elapsed", time.time()-t0)
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "")
