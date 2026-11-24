
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json, math
import mpmath as mp
import sympy as sp
import w5lib as W
import pipe

D = 8
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
K = [list(map(int, r)) for r in json.load(open(pipe.cand_dir(D)+"/kernel.json"))["K"]]
vals = pipe.load_vals(D, 1150)
N = len(cands)
li5v = [vals[f"li5:{i}"] for i in range(N)]
mp.mp.dps = 100
print("kernel rows (12x16):")
for j,row in enumerate(K):
    cv = sum(mp.mpf(c)*v for c,v in zip(row, li5v))
    nz = [(i,c) for i,c in enumerate(row) if c]
    print(f"K{j}: max|coef|={max(abs(x) for x in row)}  nonzero={nz[:6]}{'...' if len(nz)>6 else ''}  value={mp.nstr(cv,3)}")

for P in (520, 830):
    R = pipe.phase_lll(D, C, cands, K, vals, P=P, use_atoms=False)
    print(f"--- P={P}: rels={len(R['rels'])} zeros={len(R['zeros'])}")
    for mx,a,res in R['rels'][:3]:
        print(f"   REL a0={a[0]} max|a|={mx} resid={mp.nstr(res,2)} a[:8]={a[:8]}")
    for mx,a,res in R['zeros'][:5]:
        print(f"   ZERO max|a|={mx} resid={mp.nstr(res,2)} a[:8]={a[:8]}")
