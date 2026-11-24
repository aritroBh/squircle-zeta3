
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, rul1la

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')
log("master GF: D~(tau;m) = 2^(5/2) int (gd x)^{2m} H(x)[(1+e^{-x/tau})^tau -1]dx verified vs series (diffs = omitted n>=10 terms)")

t0 = time.time()
# Build relations at weight 6, families separately, get ranks (exact Fraction RREF)
W = 6
cols = sorted(rul1.conv_words_of_weight(W))
print("cols:", len(cols))
for fams in [('DS',), ('DS','DUAL'), ('DS','DUAL','RDS'), ('DS','DUAL','RDS','DIST')]:
    t = time.time()
    rels = rul1.relations_at_weight(W, families=fams)
    tbuild = time.time()-t
    rows = [rul1la.Row(r, {i: Fr(1)}) for i, (fam, tag, r) in enumerate(rels)]
    t = time.time()
    rref, piv = rul1la.eliminate(rows, cols)
    print(f"families={fams}: rels={len(rels)} build={tbuild:.1f}s rank={len(piv)} codim={len(cols)-len(piv)} rref={time.time()-t:.1f}s")
log("w=6 relation ranks computed (DS / +DUAL / +RDS / +DIST)")
