import sys,time,pickle,os
sys.path.insert(0,'.')
import sympy as sp
from fractions import Fraction as Fr
import trunc_engine as te
from fastred import QSI, ZERO
W='../..'
d_,m_=sp.symbols('d m'); I=sp.I
A={dl:sp.sympify(pickle.load(open(f'{W}/a{dl}_closed.pkl','rb')),locals={'d':d_,'m':m_,'I':I}) for dl in (9,10)}
def qdiv(num,den):
    # (num/den) in Q(sqrt2,i): QSI p0..p3 ; use sympy for exactness via conj trick: convert to sympy
    return sp.nsimplify(0)
def to_sp(q):
    return q.p0+q.p1*sp.sqrt(2)+(q.p2+q.p3*sp.sqrt(2))*I
degs=[int(a) for a in sys.argv[1:]]
for dd in degs:
    t0=time.time()
    S2=te.stirling_table(dd+2)
    rho=te.rho_level(dd,10,S2)
    bad=0;n=0
    for key,ser in rho.items():
        p,q=key
        r0=to_sp(ser[0])
        for dl in (9,10):
            if 2*dd-4<dl: continue
            lhs=sp.radsimp(to_sp(ser[dl])/r0)
            lhs=sp.nsimplify(sp.simplify(lhs)) if False else sp.simplify(lhs)
            rhs=sp.simplify(A[dl].subs({d_:dd,m_:p}))
            n+=1
            if sp.simplify(lhs-rhs)!=0:
                bad+=1; print('BAD',dd,key,dl,lhs,rhs)
    print(f'd={dd}: {n} cells (all keys, delta=9,10) bad={bad} ({time.time()-t0:.1f}s)',flush=True)
