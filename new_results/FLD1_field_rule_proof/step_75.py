WD='<project>/research/08_theory/FLD1_field_rule_proof'
code = r'''"""step33_k6exact.py -- k=6 exact check: [c^N]F == mu_{6,n} (n<=60), non-positive & odd == 0."""
import pickle, json
from fractions import Fraction as Fr
import sys
sys.path.insert(0,'.')
from hpl import coeffs_exact
store = pickle.load(open('hpl_constructed_K6.pkl','rb'))
MU = json.load(open('mu_exact_K6_N60.json'))
allok=True
for zkey, terms in sorted(store.items()):
    zks=','.join(map(str,zkey))
    bad=0
    for N in range(-6,121):
        tot=Fr(0)
        for (e,w,q) in terms:
            if len(w)==0:
                if N==e: tot+=q
                continue
            idx=N-e
            if idx>=1: tot+=q*coeffs_exact(w,idx+2)[idx]
        target=Fr(0)
        if N>=2 and N%2==0: target=Fr(MU['6'][str(N//2)].get(zks,'0'))
        if tot!=target:
            bad+=1
            if bad<=2: print('   mismatch',zkey,N,tot,target)
    if bad: allok=False
    print(f'k=6 z={zkey}: [c^N] N=-6..120 == mu (non-pos & odd ==0): {bad==0}')
print('k=6 exact (n<=60):', allok)
'''
open(WD+'/step_33_k6exact.py','w').write(code)
import subprocess
VENV='<project>/.venv/bin/python'
r = subprocess.run(['nice',VENV,'step_33_k6exact.py'],capture_output=True,text=True,cwd=WD)
print(r.stdout); print(r.stderr[-1500:] if r.returncode else '')