import sys, os, pickle
from fractions import Fraction as Fr
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
os.nice(10)
import sympy as sp
samples = pickle.load(open('e61_samples.pkl','rb'))
def try_fit(deg, denfun, nrows=60):
    mons = [(i,j) for i in range(0,deg+1) for j in range(0,deg+1-i)]
    N = len(mons)
    A = []; b = []
    for dd, p, (re, im) in samples[:nrows]:
        D = denfun(dd)
        rowv = [dd**i * p**j for (i,j) in mons]
        A.append([Fr(v) for v in rowv] + [Fr(0)]*N); b.append(re*D)
        A.append([Fr(0)]*N + [Fr(v) for v in rowv]); b.append(im*D)
    M = sp.Matrix(A)
    try:
        xs = M.gauss_jordan_solve(sp.Matrix(b))
    except ValueError as e:
        return None, mons
    return xs[0], mons
pi2 = lambda dd: Fr((2*dd-3)*(2*dd-5))
pi3 = lambda dd: Fr((2*dd-3)*(2*dd-5)*(2*dd-7))
for deg, dname, df in [(4,'Pi2',pi2),(5,'Pi3',pi3),(6,'Pi3',pi3),(6,'Pi2*?',None)]:
    if df is None: continue
    sol, mons = try_fit(deg, df)
    print(f'deg {deg} / {dname}: {"SOLVED" if sol is not None else "no solution"}')
    if sol is not None:
        N = len(mons)
        for idx,(i,j) in enumerate(mons):
            re = sp.Rational(sol[idx]); im = sp.Rational(sol[N+idx])
            if re or im: print('   d^%d m^%d: %s + %s i' % (i,j,re,im))
        pickle.dump({'deg':deg,'mons':mons,'sol':[sp.Rational(sol[t]) for t in range(2*N)]}, open('e61_fit.pkl','wb'))
        break
