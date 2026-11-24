import sys, os, pickle
from fractions import Fraction as Fr
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
os.nice(10)
data = pickle.load(open('e61_fit_system.pkl','rb'))
rows, mons = data['rows'], data['mons']
# Build rational matrix: unknowns x = [re_0..re_20, im_0..im_20]
N = len(mons)
A = []; b = []
for row, (tre, tim) in rows:
    A.append([Fr(v) for v in row] + [Fr(0)]*N); b.append(tre)
    A.append([Fr(0)]*N + [Fr(v) for v in row]); b.append(tim)
import sympy as sp
M = sp.Matrix(A)
xs = M.gauss_jordan_solve(sp.Matrix(b))
print('solved; unique:', xs[1].shape[0] == 0 if hasattr(xs[1],'shape') else xs)
sol = xs[0]
coefs = {}
for idx,(i,j) in enumerate(mons):
    re = sp.Rational(sol[idx]); im = sp.Rational(sol[N+idx])
    if re or im: coefs[(i,j)] = (re, im)
print('nonzero coefs (d^i m^j):')
for k,v in coefs.items(): print(' ', k, v)
pickle.dump({'mons': mons, 'coefs': coefs}, open('e61_fit.pkl','wb'))
