
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

# ---------- Goncharov S2/S3 recursion over concrete points ----------
INF='INF'
def diff_vec(a,b):
    """letter vector of (a-b) (FT points; 0=ZERO,1=ONE,INF). None if a==b, 'INF' if inf involved."""
    if a==b: return None
    if a==INF or b==INF: return INF
    d=E.fsub(a,b)
    if E.feq0(E.fsub(d,E.ONE)) or E.feq0(E.fadd(d,E.ONE)): return {}  # +-1 torsion
    return vec(d)
def ratio_vec(nums, dens):
    inf=0; v={}
    for d in nums:
        if d==INF: inf+=1
        elif d is None: return ('ZERO',None)
        else: v=vadd(v,d)
    for d in dens:
        if d==INF: inf-=1
        elif d is None: return ('ZERO',None)
        else: v=vadd(v,vscale(d,-1))
    if inf!=0: return ('ZERO',None)
    if not v: return ('TORSION',None)
    return ('OK',v)
def tmul(T, vv):
    out={}
    for key,c in T.items():
        for l,e in vv.items():
            k2=key+(l,); out[k2]=out.get(k2,Fr(0))+c*e
            if out[k2]==0: del out[k2]
    return out
def s1_vec(x):
    """symbol of I(0;x;1)=log((1-x)/(-x)): [1-x]-[x]."""
    if E.feq0(x): return {'ZERO':-1}
    if E.feq0(E.fsub(x,E.ONE)): return {'ZERO':1}
    return vadd(vec(E.fsub(E.ONE,x)), vscale(vec(x),-1))
def S2p(x,y):
    """symbol of I(0;x,y;1) as 2-tensor dict."""
    out={}
    st,rv=ratio_vec([diff_vec(y,x)],[diff_vec(x,E.ZERO)])
    if st=='OK':
        for k,c in s1_vec(y).items():
            for l,e in rv.items(): out[(k,l)]=out.get((k,l),0)+c*e
    st,rv=ratio_vec([diff_vec(E.ONE,y)],[diff_vec(y,x)])
    if st=='OK':
        for k,c in s1_vec(x).items():
            for l,e in rv.items(): out[(k,l)]=out.get((k,l),0)+c*e
    return {k:Fr(v) for k,v in out.items() if v}
def S3(a,b,c):
    """symbol of I(0;a,b,c;1), distinct letters, convergent."""
    out={}
    def add(T2, nums, dens):
        st,rv=ratio_vec(nums,dens)
        if st=='OK':
            for (i,j),cf in T2.items():
                for l,e in rv.items():
                    out[(i,j,l)]=out.get((i,j,l),Fr(0))+cf*e
    add(S2p(b,c),[diff_vec(b,a)],[diff_vec(a,E.ZERO)])
    add(S2p(a,c),[diff_vec(c,b)],[diff_vec(b,a)])
    add(S2p(a,b),[diff_vec(E.ONE,c)],[diff_vec(c,b)])
    return {k:v for k,v in out.items() if v}
# ---------- atom symbols ----------
X={i:None for i in [1,2,5,6]}
def init():
    g=pickle.load(open('geom14.pkl','rb'))
    for i in [1,2,5,6]: X[i]=g['X'][i]
    return g
def sym_T3a(i,j,k):
    """T3a = -(1/x_i)[G(ri,rj,rk)+G(ri,rk,rj)]; return symbol with FT scalar -(1/x_i)."""
    ri,rj,rk = E.finv(X[i]),E.finv(X[j]),E.finv(X[k])
    S = sadd(S3(ri,rj,rk), S3(ri,rk,rj))
    return S, E.fneg(E.finv(X[i]))
def sym_T3p(i,j):
    """T3p = -(2/x_i) G(ri,rj,rj) = -(1/x_i) N(rj,ri); N = 2Li3(c)-2Li3(cX)+prod.
    c = x/(x-y), cX=(x-1)/(x-y), x=rj=1/x_j, y=ri=1/x_i."""
    x=E.finv(X[j]); y=E.finv(X[i])
    c=E.fdiv(x,E.fsub(x,y)); cX=E.fdiv(E.fsub(x,E.ONE),E.fsub(x,y))
    S = sscale(sadd(li3_symbol(c), sscale(li3_symbol(cX),-1)), 2)
    return S, E.fneg(E.finv(X[i]))
def sym_Ax(i):
    """A(x_i) = Li3(x_i) - Li3(1-x_i) + zeta3 + products."""
    return sadd(li3_symbol(X[i]), sscale(li3_symbol(E.fsub(E.ONE,X[i])),-1)), E.ONE
def sym_L3x(i):
    return li3_symbol(X[i]), E.ONE
def sym_Lam(i,j):
    """Lam = Li_{2,1}(x_i, x_j/x_i). S_Li21 (DOD3) converted by -tau (reverse tensor)."""
    x=X[i]; y=E.fdiv(X[j],X[i])
    L1=vec(E.fsub(E.fmul(x,y),E.ONE)); L2=vec(E.fmul(x,y)); L3=vec(E.fsub(y,E.ONE))
    L4=vec(E.fdiv(E.fmul(x,E.fsub(y,E.ONE)),E.fsub(x,E.ONE))); L5=vec(E.fsub(x,E.ONE)); L6=vec(x)
    Sdod = sadd(tensor(L1,L2,L3),tensor(L1,L4,L6),tensor(L3,L5,L6))
    # -tau: reverse tensor, negate
    out={}
    for (a,b,c),cf in Sdod.items():
        out[(c,b,a)]=out.get((c,b,a),Fr(0))-cf
        if out[(c,b,a)]==0: del out[(c,b,a)]
    return out, E.ONE
