step=r'''import sys, pickle
sys.path.insert(0,'.')
import sympy as sp
from fastred import *
from fastred import v
record = pickle.load(open('record_W2b.pkl','rb'))
S2 = sp.sqrt(2)

g20   = (1+S2)/4 - S2/(2*(v**2+1))
g11   = -S2*(v**2+2*v-1)/(4*(v**2+1))
g20_pf = pf_from_sympy(g20); g11_pf = pf_from_sympy(g11); ng20_pf = pf_neg(g20_pf)

bad=[]
for k in range(2,13):
    G = record[(k,(2,0))]['Gm']
    if pf_sub(G,g20_pf)[1] or pf_sub(G,g20_pf)[2] or not pf_sub(G,g20_pf)[3].is_zero() or pf_sub(G,g20_pf)[0]:
        bad.append(('Gm20',k))
print("Gm_20^(k) == g20 for k=2..12:", "ALL MATCH" if not bad else bad)
bad=[]
for k in range(2,13):
    G = record[(k,(1,1))]['Gm']
    tgt = ng20_pf if k==2 else g11_pf
    dff = pf_sub(G,tgt)
    nz = (dff[0] and any(not c.is_zero() for c in dff[0].values())) or any(not e.is_zero() for e in dff[1]) or any(not e.is_zero() for e in dff[2]) or not dff[3].is_zero()
    if nz: bad.append(('Gm11',k))
print("Gm_11^(2) == -g20, Gm_11^(k) == g11 (k=3..12):", "ALL MATCH" if not bad else bad)

# Final identity with symbolic eps (WT2 step_56/57 redone)
eps = sp.Symbol('eps')
a = 1+S2; b = S2-1; I = sp.I
R10   = (1-S2/2)/v - 2/(v+a)
Rw2   = 1/(2*v) + 1/(2*(v-a)) - 1/(v+a) - 1/(2*(v+b))
DL1s = 1/(v-b) + 1/(v+a) - 1/v
DL2s = 1/(v-a) + 1/(v+b) - 1/v
def F(R): return sp.together(R.subs(v, 1/v)/v**2)
Gm20 = eps**2/(1-eps) * g20
Gm11 = -eps**2 * g20 + eps**3/(1-eps) * g11
Wcas = -2*Gm20*DL1s - Gm11*DL2s - F(sp.together(Gm11*DL1s))
num = sp.expand(sp.together((1-eps)*Wcas))
c2 = sp.together(num.coeff(eps, 2)); c3 = sp.together(num.coeff(eps, 3))
rem = sp.together(num - c2*eps**2 - c3*eps**3)
print("(1-eps)Wcas has only eps^2, eps^3 terms:", rem == 0)
print("c2 == Rw2:", sp.simplify(sp.radsimp(sp.together(c2 - Rw2)))==0)
print("c3 == 0:", sp.simplify(sp.radsimp(sp.together(c3)))==0)
pf = sp.apart(sp.together(Wcas), v, extension=[S2, I])
ipoles = [t for t in sp.Add.make_args(sp.expand(pf)) if sp.denom(sp.together(t)).has(I)]
print("poles at +-i in Wcas (all eps):", ipoles if ipoles else "NONE")
print("=> W_10(eps) = eps R10 + eps^2/(1-eps) Rw2 ; [eps^k] SimplePart = Rw2 for k>=2")
'''
open('step_40.py','w').write(step)
import subprocess
r=subprocess.run(['<project>/.venv/bin/python','step_40.py'],capture_output=True,text=True)
print(r.stdout); print(r.stderr[-800:] if r.stderr else '')
