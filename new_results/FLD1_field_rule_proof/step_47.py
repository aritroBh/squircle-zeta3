WD='<project>/research/08_theory/FLD1_field_rule_proof'
with open(WD+'/progress.log','a') as fh:
    fh.write("step26: constructed k=2,3 forms == fitted B2/B3 forms (Taylor agree to c^40)\n")
    fh.write("step27: VALIDATION k=4,5: exact mu n<=100 all zeta-monos OK; corner_series 4 angles >=40 digits (76-78 at 2pi3,pi2,3pi4; 42 at pi3, corner_series-limited)\n")
code = r'''"""step28_analyze.py -- weight bound, letters, denominators of the constructed HPL forms."""
import pickle
from fractions import Fraction as Fr
store = pickle.load(open('hpl_constructed_K5.pkl','rb'))
def wt(zkey): return sum((i+2)*e for i,e in enumerate(zkey))
for k in range(6):
    for (kk,zk),terms in sorted(store.items()):
        if kk!=k: continue
        zw = wt(zk)
        maxw = max((len(w) for (e,w,q) in terms), default=0)
        letters = set()
        denom_e = set()
        for (e,w,q) in terms:
            letters |= set(w)
            denom_e.add(e)
        viol = [(len(w),zw) for (e,w,q) in terms if zw+len(w) > k+1]
        print(f'k={k} z={zk} (wt {zw}): maxHPLwt={maxw}, wt(z)+maxHPLwt={zw+maxw} (bound {k+1}), '
              f'letters={sorted(letters)}, c-exponents={sorted(denom_e)}, weight-violations={len(viol)}')
# is -1 needed? count terms with -1
print()
for k in (2,3,4):
    n_with_minus1 = 0; tot=0
    for (kk,zk),terms in store.items():
        if kk!=k or zk!=(0,0,0,0): continue
        for (e,w,q) in terms:
            tot+=1
            if -1 in w: n_with_minus1+=1
    print(f'k={k} rational part: {n_with_minus1}/{tot} terms contain letter -1')
'''
open(WD+'/step_28_analyze.py','w').write(code)
import subprocess
VENV='<project>/.venv/bin/python'
r = subprocess.run(['nice',VENV,'step_28_analyze.py'],capture_output=True,text=True,cwd=WD)
print(r.stdout); print(r.stderr[-2000:] if r.returncode else '')