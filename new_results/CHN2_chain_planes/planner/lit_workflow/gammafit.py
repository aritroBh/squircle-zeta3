"""Fit log(lambda_pm) = sum_k (-1)^k zeta_k/k * sum_j m_j L_j(a)^k  with L_j = p*x+q*y (integer p,q) on Pi1, Pi2 and generic slice."""
import sys; sys.dont_write_bytecode=True
sys.path.insert(0,'.')
from fractions import Fraction as F
import sympy as sp
import loglam as LL
import phi
def Pk_table(pts,which):
    rows=[]
    for (x,y,a) in pts:
        X11,X21,lp,lm=LL.lam(a)
        X={'p':lp,'m':lm,'11':X11}[which]
        L=LL.logseries(X)
        P={}
        for k in range(2,9):
            v=L[k]
            assert all(t==0 for t in v[1:]), (a,k,v)
            P[k]=(-1)**k*k*v[0]
        rows.append(((x,y),P))
    return rows
def fit(rows,pmax=4,qmax=4,qmin=0,pmin=0):
    forms=[(p,q) for p in range(pmin,pmax+1) for q in range(qmin,qmax+1) if (p,q)!=(0,0)]
    A=[];b=[]
    for (x,y),P in rows:
        for k in range(2,9):
            A.append([sp.Integer(p*x+q*y)**k for (p,q) in forms]); b.append(sp.Rational(P[k].numerator,P[k].denominator))
    A=sp.Matrix(A); b=sp.Matrix(b)
    sol,params=A.gauss_jordan_solve(b)
    return forms,sol,params
def run(label,famfn,which):
    pts=[]
    for x in (1,2,3):
        for y in (1,2,3,4):
            if x==y and label.startswith('Pi'): pass
            pts.append((x,y,famfn(x,y)))
    rows=Pk_table(pts,which)
    try:
        forms,sol,params=fit(rows)
        print(label,which,'solution (nonzero multiplicities), free params:',list(params)) 
        print('   ',{f:s for f,s in zip(forms,sol) if s!=0})
    except Exception as e:
        print(label,which,'NO exact solution with p,q in 0..4:',type(e).__name__,str(e)[:80])
if __name__=='__main__':
    run('Pi1: a=(x,x,y,y)',lambda x,y:(x,x,y,y),'m')
    run('Pi1: a=(x,x,y,y)',lambda x,y:(x,x,y,y),'p')
    run('Pi2: a=(x,y,x,y)',lambda x,y:(x,y,x,y),'m')
    run('Pi2: a=(x,y,x,y)',lambda x,y:(x,y,x,y),'p')
