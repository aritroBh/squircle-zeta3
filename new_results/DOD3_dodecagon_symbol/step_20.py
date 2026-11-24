import sys, os, pickle, time
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
import sympy as sp
from dod3_symbol import *

sols = pickle.load(open('solsF_B8.pkl','rb'))
Q = Quotient()
t0 = time.time()
zs = sorted(sols.keys(), key=pfloat)
gens = []
for z in zs:
    gens.append(Q.coord(S_Li3(z)))
print("computed %d Li3 symbols, t=%.1fs"%(len(gens), time.time()-t0))
M = sp.Matrix.hstack(*[g.T for g in gens] + [d.T for d in Q.D])
r = M.rank()
print("rank [S(Li3(z)) : z 129 sols | Decomp] =", r, " (full = 27)")
pickle.dump({'zs':zs, 'gens':gens}, open('li3_gens.pkl','wb'))
