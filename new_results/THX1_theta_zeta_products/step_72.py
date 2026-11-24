
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
colpos = {c:i for i,c in enumerate(cols9)}
p0 = 4611686018427387847
B = pickle.load(open('w9_indep_rows.pkl','rb'))['indep_rows']
t=time.time()
AB = flint.nmod_mat(len(B), len(cols9), p0)
for i, bi in enumerate(B):
    for k, v in rels[bi].items():
        AB[i, colpos[k]] = (v.numerator % p0) * pow(v.denominator % p0, -1, p0) % p0
print("submatrix built:", time.time()-t, flush=True)
t=time.time()
pivcols, _ = thx1_cert.rref_pivots(AB)
print("rref(R_B): pivcols", len(pivcols), "in %.0fs"%(time.time()-t), flush=True)
st = {'rank': len(pivcols), 'pivcols': pivcols, 'indep_rows': B, 'colpos': colpos, 'prime': p0}
pickle.dump(st, open('struct_w9.pkl','wb'))
log("w9 structure complete: codim %d (DG 55)" % (len(cols9)-len(pivcols)))
