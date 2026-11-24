sym_src = r'''
from fractions import Fraction as Fr
import sys, math, pickle
import numpy as np, sympy as sp
sys.path.insert(0,'.')
import engine14 as E
AL=(Fr(0),Fr(1),Fr(0),1); UU=(Fr(-1),Fr(1),Fr(1),1)
def ft_norm(x):
    A,B,C,d=x; t=sp.symbols('t')
    f=sp.Rational(A)+sp.Rational(B)*t+sp.Rational(C)*t**2; g=t**3+t**2-2*t-1
    return Fr(sp.resultant(f,g,t))/(d**3)
LOG2=math.log(2); EMB=[2*math.cos(2*math.pi/7),2*math.cos(4*math.pi/7),2*math.cos(6*math.pi/7)]
UEMB=[e*e+e-1 for e in EMB]
def ft_embeds(x):
    A,B,C,d=x; A,B,C=float(A),float(B),float(C)
    return [(A+B*e+C*e*e)/d for e in EMB]
def factor_ft(x):
    if x[0]==0 and x[1]==0 and x[2]==0: return None
    N=abs(Fr(ft_norm(x))); num,den=N.numerator,N.denominator
    e0n=0
    while num%2==0 and num!=0: num//=2; e0n+=1
    e0d=0
    while den%2==0 and den!=0: den//=2; e0d+=1
    if num!=1 or den!=1: return ('NOT2',num,den)
    e3=e0n-e0d
    if e3%3!=0: return ('BAD3',e3)
    e0=e3//3
    zv=[zz/(2.0**e0) for zz in ft_embeds(x)]
    la=[math.log(abs(EMB[i])) for i in range(3)]; lu=[math.log(abs(UEMB[i])) for i in range(3)]
    lz=[math.log(abs(zz)) for zz in zv]
    M=np.array([[la[i],lu[i]] for i in range(3)]); b=np.array(lz)
    sol,_,_,_=np.linalg.lstsq(M,b,rcond=None); e1,e2=int(round(sol[0])),int(round(sol[1]))
    cand=E.fmul(E.fpow(AL,e1) if e1>=0 else E.fpow(E.finv(AL),-e1), E.fpow(UU,e2) if e2>=0 else E.fpow(E.finv(UU),-e2))
    cand=E.fs(cand,Fr(2)**e0)
    ratio=E.fdiv(x,cand)
    if E.feq0(E.fsub(ratio,E.ONE)): return (1,e0,e1,e2)
    if E.feq0(E.fadd(ratio,E.ONE)): return (-1,e0,e1,e2)
    return ('UNITFAIL',(e0,e1,e2))
def vec(x):
    """letter vector dict {0:e0,1:e1,2:e2} (letters 2,alpha,u), sign dropped."""
    f=factor_ft(x)
    assert f is not None and isinstance(f[0],int),(x,f)
    _,e0,e1,e2=f
    out={}
    if e0: out[0]=Fr(e0)
    if e1: out[1]=Fr(e1)
    if e2: out[2]=Fr(e2)
    return out
def vadd(*vs):
    r={}
    for v in vs:
        for k,c in v.items():
            r[k]=r.get(k,Fr(0))+c
            if r[k]==0: del r[k]
    return r
def vscale(v,q):
    q=Fr(q); return {k:c*q for k,c in v.items() if c*q}
def tensor(u,v,w,coeff=Fr(1)):
    S={}
    for a,ca in u.items():
        for b,cb in v.items():
            for c,cc in w.items():
                S[(a,b,c)]=S.get((a,b,c),Fr(0))+coeff*ca*cb*cc
    return {k:x for k,x in S.items() if x}
def sadd(*Ss):
    r={}
    for S in Ss:
        for k,c in S.items():
            r[k]=r.get(k,Fr(0))+c
            if r[k]==0: del r[k]
    return r
def sscale(S,q):
    q=Fr(q); return {k:c*q for k,c in S.items() if c*q}
def li3_symbol(z):
    """PGEN1R: S(Li3(z)) = z (x) z (x) (1-z)."""
    zv=vec(z); omz=vec(E.fsub(E.ONE,z))
    return tensor(zv,zv,omz)
'''
open('sym14.py','w').write(sym_src)
# quick test of li3_symbol
import sym14 as S
g=pickle.load(open('geom14.pkl','rb'))
for i in [1,2,5,6]:
    s=S.li3_symbol(g['X'][i])
    print("li3(x%d): %d terms"%(i,len(s)))
print("sym14.py base OK")
