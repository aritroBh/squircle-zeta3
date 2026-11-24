import os, sys, json; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log
mp.mp.dps = 300
consts = {k: mp.mpf(v) for k,v in json.load(open('purew3_300.json')).items()}
names = list(consts.keys()); xs = [consts[k] for k in names]
n = len(xs)
print("n constants:", n)

from flint import fmpz_mat, fmpz
def find_relations(vals, dps, scale_exp=None, thr_ratio=1e4):
    old=mp.mp.dps; mp.mp.dps=dps
    if scale_exp is None: scale_exp = dps-10
    B = mp.power(10, scale_exp)
    n=len(vals)
    rows=[]
    for i in range(n):
        row=[0]*n + [int(mp.nint(B*vals[i]))]
        row[i]=1
        rows.append(row)
    M=fmpz_mat(rows)
    R = M.lll()
    rels=[]
    for r in R.rows():
        v=list(r)
        last=abs(v[-1]); coef=v[:-1]
        maxc=max(abs(c) for c in coef)
        rels.append((last,maxc,coef))
    mp.mp.dps=old
    return rels

rels = find_relations(xs, 300, scale_exp=285)
# sort by last (residual), print smallest
rels.sort(key=lambda z:z[0])
print("smallest-residual vectors (last coord = |residual*B|, maxc):")
for last,maxc,coef in rels[:12]:
    print(f"   residB={last:>25}  maxc={maxc:<12} nz={sum(1 for c in coef if c!=0)}")
log("step: LLL on 33 pure weight-3 constants; printed smallest-residual vectors")
