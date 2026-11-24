import os, json
from fractions import Fraction as Ff
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
step = r'''
import sys, json, time
from fractions import Fraction as Ff
import sympy as sp
sys.path.insert(0,'.')
S = json.load(open('d17_min_setup.json'))
N=S['N']; M=S['M']; L=1+N+M
c0=[Ff(x) for x in S['c_rel']]
Z=[[Ff(x) for x in row] for row in S['Z']]
nZ=len(Z)
def maxnum(c): return max(abs(x.numerator) for x in c)
def nli5(c): return sum(1 for i in range(N) if c[1+i])
def nterms(c): return sum(1 for x in c if x)
Zm = sp.Matrix(Z)          # 43 x L
used = [i for i in range(N) if c0[1+i]]
used_sorted = sorted(used, key=lambda i: -abs(c0[1+i].numerator))
# zero the top-43 largest-coef used li5 coords
I = sorted(1+used_sorted[:43])
ZI = Zm[:, I]              # 43 x 43
b = -sp.Matrix([c0[i] for i in I])
t0=time.time()
det = ZI.det()
print(f"det(ZI) = {det}  (nonzero: {det!=0})  ({time.time()-t0:.0f}s)")
if det:
    t = ZI.solve(b)        # 43-vector
    tt = [Ff(t[j]) for j in range(nZ)]
    c = list(c0)
    for j in range(nZ):
        if tt[j]:
            for k in range(L):
                if Z[j][k]: c[k]+= tt[j]*Z[j][k]
    print(f"after zeroing 43 largest: maxnum={maxnum(c)} li5={nli5(c)} terms={nterms(c)}")
    print("  zeros at chosen coords:", all(c[i]==0 for i in I))
    json.dump(dict(c=[str(x) for x in c]), open('d17_min_elim43.json','w'))
'''
open('step_04g_d17_elim.py','w').write(step)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_04g_d17_elim.py'], capture_output=True, text=True, cwd='.', timeout=1200)
print(r.stdout); print("ERR:", r.stderr[-1500:])
