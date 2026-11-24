import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import symengine as S, engine as E
# smoke: Wmax=6
forms=S.symbolic_forms(6)
# homogeneity check
for W,fl in forms.items():
    for b,p in enumerate(fl):
        for e in p: assert sum(e)==W, (W,b,e)
print('homogeneity OK')
# cross-check vs single-direction at several directions
import random
random.seed(1)
ok=True
for _ in range(5):
    a=tuple(random.randint(1,4) for _ in range(4))
    r=E.reduced_coeffs(a,6)
    for W in range(2,7):
        for b in range(len(r[W])):
            v=S.peval(forms[W][b], a)
            if v!=r[W][b]:
                ok=False; print('MISMATCH', a, W, b, v, r[W][b])
print('cross-check single-direction vs symbolic:', 'PASS' if ok else 'FAIL', 't=%.1fs'%(time.time()-t0))
