import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import symengine as S, engine as E, mzv_exact as M, pickle, random
forms=S.symbolic_forms(11)
print('forms computed, t=%.0fs'%(time.time()-t0))
# save
def ser(p): return {','.join(map(str,e)): str(c) for e,c in p.items()}
out={}
for W,fl in forms.items():
    for b,p in enumerate(fl):
        out[(W,b)]=ser(p)
pickle.dump(out, open('forms.pkl','wb'))
print('saved forms.pkl, t=%.0fs'%(time.time()-t0))
# cross-check vs single-direction at 6 random directions, all W 2..11
random.seed(2); ok=True
for _ in range(6):
    a=tuple(random.randint(1,5) for _ in range(4))
    r=E.reduced_coeffs(a,11)
    for W in range(2,12):
        T=E.table(W)
        for b in range(T.nullity):
            v=S.peval(forms[W][b], a)
            if v!=r[W][b]:
                ok=False; print('MISMATCH', a, W, b)
print('full cross-check:', 'PASS' if ok else 'FAIL', 't=%.0fs'%(time.time()-t0))
open('progress.log','a').write('step_24-25: symbolic engine Wmax=11 computed and saved forms.pkl; cross-check vs single-direction at 6 random directions all W: PASS (t=%.0fs)\n'%(time.time()-t0))
