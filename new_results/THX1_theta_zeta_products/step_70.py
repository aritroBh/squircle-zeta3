
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, thx1_cert
import flint

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

rels9_data = pickle.load(open('rels9_dsdual.pkl','rb'))
rels = [{k: Fr(n, d) for k, (n, d) in r.items()} for fam, tag, r in rels9_data]
cols9 = pickle.load(open('cols9.pkl','rb'))
p0 = 4611686018427387847
t=time.time()
A, colpos = thx1_cert.build_nmod(rels, cols9, p0)
print("matrix built:", time.time()-t, A.nrows(), A.ncols(), flush=True)
t=time.time()
pivcols, _ = thx1_cert.rref_pivots(A)
print("rref(R): pivcols", len(pivcols), "in", time.time()-t, "s", flush=True)
pickle.dump({'pivcols': pivcols}, open('w9_pivcols.pkl','wb'))
log("w9 rref(R) done: rank %d" % len(pivcols))
t=time.time()
At = A.transpose()
pivrows, _ = thx1_cert.rref_pivots(At)
print("rref(R^T): pivrows", len(pivrows), "in", time.time()-t, "s", flush=True)
st = {'rank': len(pivcols), 'pivcols': pivcols, 'indep_rows': pivrows, 'colpos': colpos, 'prime': p0}
pickle.dump(st, open('struct_w9.pkl','wb'))
log("w9 structure complete: codim %d (DG 55)" % (len(cols9)-len(pivcols)))
