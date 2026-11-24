# codimension of my DS+DUAL relation space at weights 2..7 vs Fibonacci/DG dims 2,3,5,8,13,21 ; also negative control for numerics
from mzv import *
import itertools
P = 2**61-1
def allconv(w): return [x for x in itertools.product((0,1,-1), repeat=w) if conv(x)]
def rels_at(W):
    cw = {w: allconv(w) for w in range(1, W)}
    out = []
    for a in range(1, W//2+1):
        for u in cw[a]:
            for v in cw[W-a]:
                out.append(DS(u, v))
    for w in allconv(W): out.append(DUAL(w))
    return out
def rank_mod(rows, cols):
    idx = {c: i for i, c in enumerate(cols)}
    piv = {}
    r = 0
    for row in rows:
        vec = {idx[k] % 10**9 if False else idx[k]: v % P for k, v in row.items() if v % P}
        while vec:
            lead = min(vec)
            if lead in piv:
                pr = piv[lead]; f = vec[lead]
                for k, v in pr.items():
                    nv = (vec.get(k, 0) - f*v) % P
                    if nv: vec[k] = nv
                    else: vec.pop(k, None)
            else:
                inv = pow(vec[lead], -1, P)
                piv[lead] = {k: v*inv % P for k, v in vec.items()}; r += 1; break
    return r
for W in range(2, 8):
    cols = allconv(W)
    R = rels_at(W)
    rk = rank_mod(R, cols)
    print(W, "columns", len(cols), "rank", rk, "codim", len(cols)-rk, flush=True)
