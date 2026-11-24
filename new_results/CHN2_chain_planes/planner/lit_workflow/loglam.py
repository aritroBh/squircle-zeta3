"""On Pi1/Pi2: lambda_pm = F^{(23)} +- F^{(32)} (s23=0 so F^{(2)}=1).  Is log(lambda_pm) a sum of single zetas only (Gamma-product)?"""
import sys; sys.dont_write_bytecode=True
sys.path.insert(0,'.')
from fractions import Fraction as F
import phi
E=phi.E; Mz=phi.M
Wmax=8
# product tables
mono={W:E.table(W).monomials for W in range(2,Wmax+1)}
_prod={}
def prod_table(a,b):
    if (a,b) in _prod: return _prod[(a,b)]
    T=E.table(a+b); tab=[[None]*len(mono[b]) for _ in mono[a]]
    for i,mi in enumerate(mono[a]):
        for j,mj in enumerate(mono[b]):
            tab[i][j]=T.reduce(Mz.wordcombo(tuple(mi)+tuple(mj)))
    _prod[(a,b)]=tab; return tab
def mul(X,Y):
    """X,Y: dict weight->vector (weights>=2). returns dict up to Wmax."""
    out={}
    for a,va in X.items():
        for b,vb in Y.items():
            if a+b>Wmax: continue
            tab=prod_table(a,b)
            r=out.setdefault(a+b,[F(0)]*len(mono[a+b]))
            for i,x in enumerate(va):
                if x==0: continue
                for j,y in enumerate(vb):
                    if y==0: continue
                    for k,c in enumerate(tab[i][j]): r[k]+=x*y*c
    return out
def addto(A,B,c):
    for w,v in B.items():
        r=A.setdefault(w,[F(0)]*len(v))
        for k in range(len(v)): r[k]+=c*v[k]
def logseries(X):
    L={}; P=dict(X)
    k=1
    while P:
        addto(L,P,F((-1)**(k+1),k))
        P=mul(P,X); k+=1
        P={w:v for w,v in P.items() if any(x!=0 for x in v)}
    return L
def lam(a):
    k=phi.kin_from_a(a); e0,e1=phi.e_mats(**k)
    out=phi.phi_first_column(e0,e1,Wmax)
    X11={W:out[W][0] for W in range(2,Wmax+1)}
    X21={W:out[W][1] for W in range(2,Wmax+1)}
    lp={W:[x+y for x,y in zip(X11[W],X21[W])] for W in X11}
    lm={W:[x-y for x,y in zip(X11[W],X21[W])] for W in X11}
    return X11,X21,lp,lm
def show(a,label):
    X11,X21,lp,lm=lam(a)
    print('==',label,'a=',a)
    for nm,X in (('log F11',X11),('log lam+ ',lp),('log lam- ',lm)):
        L=logseries(X)
        row=[]
        for W in range(2,Wmax+1):
            v=L.get(W,[F(0)]*len(mono[W]))
            # single-zeta content: index 0 of BASIS[W] is zeta_W for W in (2..4,6,7,8) ; weight 5: idx0=zeta5, idx1=z2z3 ; weight 6: idx0 = zeta6, idx1=z3^2
            nonsingle=[str(x) for x in v[1:]]
            row.append('W%d: z_W=%s other=%s'%(W,v[0],nonsingle))
        print(' ',nm); [print('     ',r) for r in row]
if __name__=='__main__':
    show((1,1,2,3),'generic positive (not on locus)')
    show((2,2,3,3),'Pi1')
    show((2,3,2,3),'Pi2')
