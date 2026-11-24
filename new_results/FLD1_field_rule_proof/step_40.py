import subprocess
WD='<project>/research/08_theory/FLD1_field_rule_proof'
dbg = r'''
import sys, pickle
sys.path.insert(0,'.')
from fractions import Fraction as Fr
import zsum
from hpl import coeffs_exact
store = pickle.load(open('hpl_constructed_K3.pkl','rb'))
def laurent(terms, lo, hi):
    # coefficients of c^j for j in [lo,hi]
    out={}
    for (e,w,q) in terms:
        ce = coeffs_exact(w, hi-e+2)
        for j in range(lo,hi+1):
            idx=j-e
            if 1<=idx<len(ce):
                out[j]=out.get(j,Fr(0))+q*ce[idx]
    return out
for (k,zk),terms in sorted(store.items()):
    L = laurent(terms,-8,0)
    nonpos = {j:v for j,v in sorted(L.items()) if v!=0}
    print(f'k={k} z={zk}: non-positive Laurent coeffs (should be all 0 for a clean power series):')
    print('     ', nonpos if nonpos else '{all zero}')
'''
open(WD+'/dbg8.py','w').write(dbg)
VENV='<project>/.venv/bin/python'
r=subprocess.run([VENV,'dbg8.py'],capture_output=True,text=True,cwd=WD)
print(r.stdout); print(r.stderr[-2000:] if r.returncode else '')