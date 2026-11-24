
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
print("matrix built:", time.time()-t, flush=True)
t=time.time()
At = A.transpose()
print("transpose:", time.time()-t, flush=True)
t=time.time()
pivrows, _ = thx1_cert.rref_pivots(At)
print("rref(R^T): indep rows", len(pivrows), "in %.0fs"%(time.time()-t), flush=True)
pickle.dump({'indep_rows': pivrows}, open('w9_indep_rows.pkl','wb'))
log("w9 rref(R^T) done: %d independent rows (%.0fs)" % (len(pivrows), time.time()-t))
