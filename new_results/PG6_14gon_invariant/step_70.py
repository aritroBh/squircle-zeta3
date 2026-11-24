import pickle, sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
import sympy as sp
g = S.init()
setup = pickle.load(open('q3setup.pkl','rb'))
LI3v=setup['LI3v']; Dg=setup['Dg']; basis3=setup['basis3']; CRs=setup['CRs']
bidx={t:i for i,t in enumerate(basis3)}
def t2v(T):
    v=[Fr(0)]*27
    for t,c in T.items():
        if t in bidx: v[bidx[t]]+=Fr(c)
    return v
Dgv=[t2v(T) for T in Dg]
Mcols = LI3v + Dgv
M = sp.Matrix(27, len(Mcols), lambda i,j: sp.Rational(Mcols[j][i]))
print("rank [li3|Decomp]:", M.rank())
# sigma 3 components
sigma = pickle.load(open('sigma14.pkl','rb'))
s_comp=[[Fr(0)]*27 for _ in range(3)]
for key,cFT in sigma.items():
    A,B,C,d = cFT
    if key in bidx:
        s_comp[0][bidx[key]] += Fr(A)/d
        s_comp[1][bidx[key]] += Fr(B)/d
        s_comp[2][bidx[key]] += Fr(C)/d
# solve each component: Mcomp = [LI3v | Dgv], target = s_comp[k]; keep only first 24 coeffs
sol_c = [[None]*24 for _ in range(3)]
rankM = M.rank()
ok=True
for k in range(3):
    tgt = sp.Matrix(27,1, lambda i,j: sp.Rational(s_comp[k][i]))
    if M.row_join(tgt).rank()!=rankM:
        print(f"component {k}: NOT in span (rank increases)"); ok=False; continue
    sol = M.gauss_jordan_solve(tgt)[0]
    for z in range(24):
        sol_c[k][z]=Fr(sp.Rational(sol[z]))
    print(f"component {k}: solved, residual 0")
print("all components in span:", ok)
pickle.dump(sol_c, open('sol_c.pkl','wb'))
open('progress.log','a').write('PG6 s_B13: sigma reduced to 24-CR li3 basis over F, all 3 components in span: %s\n'%ok)
