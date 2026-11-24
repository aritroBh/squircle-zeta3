import subprocess, os, time
script = '''
import os, sys, time
os.chdir('<project>/research/08_theory/CHN3_p5_locus')
sys.path.insert(0,'.')
exec(open('common.py').read())
import sympy as sp
from sympy.polys.subresultants_qq_zz import subresultants
t0=time.time()
C4,Q7 = get_C4_Q7()
C4_1 = sp.Poly(sp.expand(C4.subs(a3,1)), a5)
Q7_1 = sp.Poly(sp.expand(Q7.subs(a3,1)), a5)
print(f"setup {time.time()-t0:.1f}s; subresultants...", flush=True)
t0=time.time()
subs = subresultants(C4_1, Q7_1, a5)
print(f"subresultants done {time.time()-t0:.1f}s; degrees: {[sp.Poly(s,a5).degree() for s in subs]}", flush=True)
# degree-1 subresultant: S1 = A(a1,a2,a4) a5 + B(a1,a2,a4)
for s in subs:
    Ps = sp.Poly(s, a5)
    if Ps.degree()==1:
        A = sp.expand(Ps.nth(1)); B = sp.expand(Ps.nth(0))
        open('subres1_A.txt','w').write(str(A))
        open('subres1_B.txt','w').write(str(B))
        print("deg-1 subresultant saved; A terms:", len(sp.Poly(A,a1,a2,a4).terms()),
              " B terms:", len(sp.Poly(B,a1,a2,a4).terms()), flush=True)
    if Ps.degree()==0:
        R0 = sp.expand(Ps.nth(0))
        open('subres0.txt','w').write(str(R0))
        print("deg-0 subresultant (resultant) saved; terms:", len(sp.Poly(R0,a1,a2,a4).terms()), flush=True)
print("DONE", flush=True)
'''
open('run_subres.py','w').write(script)
PY = '<project>/.venv/bin/python'
logf = open('subres.log','w')
p = subprocess.Popen(['nice','-n','10',PY,'run_subres.py'], stdout=logf, stderr=subprocess.STDOUT,
                     cwd='<project>/research/08_theory/CHN3_p5_locus')
print("PID:", p.pid); time.sleep(25); print(open('subres.log').read())
