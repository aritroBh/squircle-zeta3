
import pickle, json, time
import mpmath as mp
mp.mp.dps = 55
from fractions import Fraction as Fr
import sympy as sp
c = sp.Symbol('c')
import itertools
C4 = "<project>/research/08_theory/COR4_M3alpha_classical/"
sa = json.load(open(C4+"sunit_args.json"))
args = [s for _,_,s in sa if s not in ('1','-1','2','1/2')]
argex = {s: sp.sympify(s) for s in args}
def lg(cc, i): return mp.log([cc, 1-cc, 1+cc][i])
def Li(n, ex, cc): return mp.re(mp.polylog(n, mp.mpf(sp.N(ex.subs(c,cc),40))))
L2 = mp.log(2)
basis = []
basis.append(('z4', lambda cc: mp.zeta(4)))
basis.append(('z3*l2', lambda cc: mp.zeta(3)*L2))
basis.append(('z2*l2^2', lambda cc: mp.zeta(2)*L2**2))
basis.append(('l2^4', lambda cc: L2**4))
for g, i in [('c',0),('1-c',1),('1+c',2)]:
    basis.append((f'z3*log({g})', lambda cc, i=i: mp.zeta(3)*lg(cc,i)))
    basis.append((f'z2*l2*log({g})', lambda cc, i=i: mp.zeta(2)*L2*lg(cc,i)))
    basis.append((f'l2^3*log({g})', lambda cc, i=i: L2**3*lg(cc,i)))
for s in args:
    ex = argex[s]
    basis.append((f'z2*Li2({s})', lambda cc, ex=ex: mp.zeta(2)*Li(2, ex, cc)))
    basis.append((f'l2^2*Li2({s})', lambda cc, ex=ex: L2**2*Li(2, ex, cc)))
    basis.append((f'l2*Li3({s})', lambda cc, ex=ex: L2*Li(3, ex, cc)))
    for h, j in [('c',0),('1-c',1),('1+c',2)]:
        basis.append((f'l2*Li2({s})*log({h})', lambda cc, ex=ex, j=j: L2*Li(2, ex, cc)*lg(cc,j)))
for (h1,i1),(h2,i2) in itertools.combinations_with_replacement([('c',0),('1-c',1),('1+c',2)],2):
    basis.append((f'z2*log({h1})*log({h2})', lambda cc, i1=i1, i2=i2: mp.zeta(2)*lg(cc,i1)*lg(cc,i2)))
    basis.append((f'l2^2*log({h1})*log({h2})', lambda cc, i1=i1, i2=i2: L2**2*lg(cc,i1)*lg(cc,i2)))
for trip in itertools.combinations_with_replacement([('c',0),('1-c',1),('1+c',2)],3):
    idx = [t[1] for t in trip]
    nm = '*'.join(f'log({t[0]})' for t in trip)
    basis.append((f'l2*{nm}', lambda cc, idx=idx: L2*mp.fprod([lg(cc,i) for i in idx])))
print(len(basis), "basis")
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/fit_resid.pkl','rb'))
pts = [mp.mpf(x) for x in D['pts']]
rv = [mp.mpf(x) for x in D['rv']]
t0 = time.time()
B = [[f(p_) for p_ in pts] for nm, f in basis]
print("basis evals", time.time()-t0)
pickle.dump({'B': [[str(x) for x in col] for col in B], 'names': [nm for nm,_ in basis]}, open('fit_basis.pkl','wb'))
