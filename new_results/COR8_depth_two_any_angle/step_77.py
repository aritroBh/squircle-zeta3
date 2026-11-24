
import pickle, json, sys, re, time
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
import mpmath as mp
mp.mp.dps = 60
from hpl_series import H_series, Li22
from fractions import Fraction as Fr
import sympy as sp
c = sp.Symbol('c')
C4 = "<project>/research/08_theory/COR4_M3alpha_classical/"
sa = json.load(open(C4+"sunit_args.json"))
args = [s for _,_,s in sa if s not in ('1','-1','2','1/2')]

t0=time.time()
# invisible basis definitions: list of (name, evaluator(cc))
L2 = mp.log(2)
def mk():
    basis = []
    basis.append(('z4', lambda cc: mp.zeta(4)))
    basis.append(('z3*l2', lambda cc: mp.zeta(3)*L2))
    basis.append(('z2*l2^2', lambda cc: mp.zeta(2)*L2**2))
    basis.append(('l2^4', lambda cc: L2**4))
    for g, gv in [('c', lambda x: x), ('1-c', lambda x: 1-x), ('1+c', lambda x: 1+x)]:
        basis.append((f'z3*log({g})', lambda cc, gv=gv: mp.zeta(3)*mp.log(gv(cc))))
        basis.append((f'z2*l2*log({g})', lambda cc, gv=gv: mp.zeta(2)*L2*mp.log(gv(cc))))
        basis.append((f'l2^3*log({g})', lambda cc, gv=gv: L2**3*mp.log(gv(cc))))
    import itertools
    logpairs = list(itertools.combinations_with_replacement([('c',0),('1-c',1),('1+c',2)],2))
    def lg(cc, i): return mp.log([cc, 1-cc, 1+cc][i])
    for s in args:
        ex = sp.sympify(s)
        basis.append((f'z2*Li2({s})', lambda cc, ex=ex: mp.zeta(2)*mp.re(mp.polylog(2, mp.mpf(sp.N(ex.subs(c,cc),40))))))
        basis.append((f'l2^2*Li2({s})', lambda cc, ex=ex: L2**2*mp.re(mp.polylog(2, mp.mpf(sp.N(ex.subs(c,cc),40))))))
        basis.append((f'l2*Li3({s})', lambda cc, ex=ex: L2*mp.re(mp.polylog(3, mp.mpf(sp.N(ex.subs(c,cc),40))))))
        for h, j in [('c',0),('1-c',1),('1+c',2)]:
            basis.append((f'l2*Li2({s})*log({h})', lambda cc, ex=ex, j=j: L2*mp.re(mp.polylog(2, mp.mpf(sp.N(ex.subs(c,cc),40))))*lg(cc,j)))
    for (h1,i1),(h2,i2) in logpairs:
        basis.append((f'z2*log({h1})*log({h2})', lambda cc, i1=i1, i2=i2: mp.zeta(2)*lg(cc,i1)*lg(cc,i2)))
        basis.append((f'l2^2*log({h1})*log({h2})', lambda cc, i1=i1, i2=i2: L2**2*lg(cc,i1)*lg(cc,i2)))
    logtrips = list(itertools.combinations_with_replacement([('c',0),('1-c',1),('1+c',2)],3))
    for trip in logtrips:
        nm = '*'.join(f'log({t[0]})' for t in trip)
        idx = [t[1] for t in trip]
        basis.append((f'l2*{nm}', lambda cc, idx=idx: L2*mp.prod([lg(cc,i) for i in idx])))
    return basis
basis = mk()
print(len(basis), "invisible basis functions, build time", time.time()-t0)
pickle.dump([b[0] for b in basis], open('invisible_names.pkl','wb'))
