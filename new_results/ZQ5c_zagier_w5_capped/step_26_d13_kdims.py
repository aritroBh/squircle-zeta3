import sys, json, time
from fractions import Fraction as Ff
import flint
sys.path.insert(0,'.')
import w5blib as W, w5pipe
from w5blib import QD
J = json.load(open('out/d13/candidates.json'))
allc=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),13), vx=c['vx'], v1m=c['v1m'], h=c['h']) for c in J['cands']]
names=J['names']; d=len(names)
def keyf(c): return (c['h'], abs(float(c['x'].emb(40))))
sortedc = sorted(allc, key=keyf)
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
for cap in (300, 350, 400, 450, 487):
    cands = sortedc[:cap]
    t0=time.time()
    kd, rank = kdim(cands)
    print(f"  cap {cap}: kernel dim={kd}, rank={rank} ({time.time()-t0:.0f}s)", flush=True)
