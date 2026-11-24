
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
Zm = sp.Matrix(Z)
used = [i for i in range(N) if c0[1+i]]
used_sorted = sorted(used, key=lambda i: -abs(c0[1+i].numerator))
t0=time.time()
# greedy build independent column set among used li5 coords
Icols=[]   # atom coords (1+i)
Zcur = sp.Matrix.zeros(nZ,0)
for i in used_sorted:
    col = 1+i
    trial = Zcur.row_join(Zm[:, col])
    if trial.rank() > len(Icols):
        Icols.append(col); Zcur = trial
    if len(Icols) >= nZ: break
print(f"independent used-li5 columns: {len(Icols)} ({time.time()-t0:.0f}s)")
# solve t * Z_I = -c0_I
ZI = Zcur                      # nZ x r
b = -sp.Matrix([c0[i] for i in Icols])
# particular solution: ZI.T t = b  (r eqs, nZ unknowns)
part = ZI.T.solve(b) if len(Icols)==nZ else None
if part is None:
    # underdetermined: use gauss_jordan_solve
    tvec, params = ZI.T.gauss_jordan_solve(b)
    t0v = [Ff(tvec[j]) for j in range(nZ)]
    print("free params:", params.rows if params is not None else 0)
else:
    t0v=[Ff(part[j]) for j in range(nZ)]
c=list(c0)
for j in range(nZ):
    if t0v[j]:
        for k in range(L):
            if Z[j][k]: c[k]+=t0v[j]*Z[j][k]
print(f"after elim {len(Icols)}: maxnum={maxnum(c)} li5={nli5(c)} terms={nterms(c)} ({time.time()-t0:.0f}s)")
print("  chosen coords zeroed:", all(c[i]==0 for i in Icols))
json.dump(dict(c=[str(x) for x in c], Icols=Icols), open('d17_min_elim.json','w'))
# nullspace freedom: B rows (null of ZI.T)
ns = ZI.T.nullspace()
print("nullspace dim:", len(ns))
if ns:
    B = sp.Matrix.hstack(*ns).T   # (nZ-r) x nZ
    ZB = (B * Zm)                 # freedom rows in atom space
    ZB_rows = [[Ff(ZB[i,j]) for j in range(L)] for i in range(ZB.rows)]
    json.dump(dict(ZB=[[str(x) for x in row] for row in ZB_rows]), open('d17_min_elimZ.json','w'))
    print("saved freedom rows:", len(ZB_rows))
print(f"({time.time()-t0:.0f}s total)")
