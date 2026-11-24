import os, sys, subprocess
lib = r'''
# w5lib.py -- shared library for ZQ5: explicit weight-5 formulas for L(chi_D,5)
# Exact Q(sqrt D) arithmetic, S-unit factorisation, exceptional S-units,
# weight-5 symbols (exact), numerics (mpmath), LLL relation search (flint).
from fractions import Fraction
import itertools, math, time, os, json, pickle
import mpmath as mp
import sympy as sp
import flint

# ---------------- exact Q(sqrt D) ----------------
class QD:
    __slots__=('a','b','D')
    def __init__(s,a,b,D): s.a=Fraction(a); s.b=Fraction(b); s.D=int(D)
    def __add__(s,o):
        if isinstance(o,(int,Fraction)): o=QD(o,0,s.D)
        assert s.D==o.D; return QD(s.a+o.a,s.b+o.b,s.D)
    __radd__=__add__
    def __sub__(s,o):
        if isinstance(o,(int,Fraction)): o=QD(o,0,s.D)
        return QD(s.a-o.a,s.b-o.b,s.D)
    def __rsub__(s,o):
        if isinstance(o,(int,Fraction)): o=QD(o,0,s.D)
        return o-s
    def __neg__(s): return QD(-s.a,-s.b,s.D)
    def __mul__(s,o):
        if isinstance(o,(int,Fraction)): return QD(s.a*o,s.b*o,s.D)
        assert s.D==o.D
        return QD(s.a*o.a+s.b*o.b*s.D, s.a*o.b+s.b*o.a, s.D)
    __rmul__=__mul__
    def conj(s): return QD(s.a,-s.b,s.D)
    def norm(s): return s.a*s.a-s.D*s.b*s.b
    def inv(s):
        n=s.norm(); assert n!=0
        return QD(s.a/n,-s.b/n,s.D)
    def __truediv__(s,o):
        if isinstance(o,(int,Fraction)): return QD(s.a/o,s.b/o,s.D)
        return s*o.inv()
    def __pow__(s,n):
        n=int(n)
        if n==0: return QD(1,0,s.D)
        if n<0: return (s.inv())**(-n)
        r=QD(1,0,s.D); b=s
        while n:
            if n&1: r=r*b
            b=b*b; n>>=1
        return r
    def __eq__(s,o): return isinstance(o,QD) and s.D==o.D and s.a==o.a and s.b==o.b
    def __hash__(s): return hash((s.a,s.b,s.D))
    def emb(s,prec=80):
        old=mp.mp.dps; mp.mp.dps=max(prec,60)
        v=mp.mpf(s.a.numerator)/s.a.denominator + (mp.mpf(s.b.numerator)/s.b.denominator)*mp.sqrt(s.D)
        mp.mp.dps=old
        return +v
    def __repr__(s): return f"({s.a}{'+' if s.b>=0 else ''}{s.b}*w{s.D})"

def _v_p(fr,p):
    fr=Fraction(fr); v=0
    num,den=abs(fr.numerator),fr.denominator
    while num%p==0 and num>0: num//=p; v+=1
    while den%p==0: den//=p; v-=1
    return v

def _allowed(fr,primes):
    fr=Fraction(fr); rem=abs(fr.numerator)*fr.denominator
    remn=abs(fr.numerator); remd=fr.denominator
    for p in primes:
        while remn%p==0 and remn>0: remn//=p
        while remd%p==0: remd//=p
    return remn==1 and remd==1

# ---------------- per-D data ----------------
def letters_info(D):
    if D==8:
        w=QD(0,1,8)
        letters=[w,1+w]; names=["ls2","lE"]; disp=["log(sqrt2)","log(1+sqrt2)"]
        unit_idx=1; unit=1+w
        nonunit=[w]
        def pexps(N):
            if not _allowed(N,[2]): return None
            return [_v_p(N,2)]
    elif D==5:
        phi=QD(Fraction(1,2),Fraction(1,2),5); r5=QD(0,1,5)
        letters=[phi,r5,QD(2,0,5)]; names=["lphi","ls5","l2"]; disp=["log(phi)","log(sqrt5)","log(2)"]
        unit_idx=0; unit=phi; nonunit=[r5,QD(2,0,5)]
        def pexps(N):
            if not _allowed(N,[2,5]): return None
            v2=_v_p(N,2)
            if v2%2!=0: return None
            return [_v_p(N,5), v2//2]
    elif D==12:
        r3=QD(0,1,12); eps=2+r3; u=1+r3
        letters=[eps,r3,u]; names=["lE3","ls3","lu"]; disp=["log(2+sqrt3)","log(sqrt3)","log(1+sqrt3)"]
        unit_idx=0; unit=eps; nonunit=[r3,u]
        def pexps(N):
            if not _allowed(N,[2,3]): return None
            return [_v_p(N,3), _v_p(N,2)]
    elif D==13:
        r13=QD(0,1,13); eps=QD(Fraction(3,2),Fraction(1,2),13)
        letters=[eps,r13,QD(2,0,13)]; names=["lE13","ls13","l2"]; disp=["log((3+sqrt13)/2)","log(sqrt13)","log(2)"]
        unit_idx=0; unit=eps; nonunit=[r13,QD(2,0,13)]
        def pexps(N):
            if not _allowed(N,[2,13]): return None
            v2=_v_p(N,2)
            if v2%2!=0: return None
            return [_v_p(N,13), v2//2]
    else:
        raise ValueError(D)
    Sdesc={8:"S={2} (ramified; 2=(sqrt2)^2)",
           5:"S={5,2} (5=(sqrt5)^2 ramified, 2 inert)",
           12:"S={3,2} (3=(sqrt3)^2, 2=(1+sqrt3)^2 both ramified)",
           13:"S={13,2} (13=(sqrt13)^2 ramified, 2 inert)"}[D]
    return dict(D=D,letters=letters,names=names,disp=disp,unit_idx=unit_idx,
                unit=unit,nonunit=nonunit,pexps=pexps,Sdesc=Sdesc,
                unit_emb=unit.emb(150))

def factor_sunit(x,C,verb=False):
    """Return coords (aligned with C['letters']) of |x| if x is an S-unit, else None.
    (Sign is ignored: -1 is torsion.)"""
    if x==QD(0,0,C['D']): return None
    N=x.norm()
    e=C['pexps'](N)
    if e is None: return None
    y=x
    for g,ei in zip(C['nonunit'],e):
        y=y/(g**ei)
    if y.norm() not in (1,-1): return None
    old=mp.mp.dps; mp.mp.dps=150
    t=mp.log(abs(y.emb(150)))/mp.log(C['unit_emb'])
    k=int(mp.nint(t))
    ok = abs(t-k)<mp.mpf('1e-80')
    mp.mp.dps=old
    if not ok: return None
    eps=C['unit']
    if y==eps**k or y==-(eps**k):
        coords=list(e)
        coords.insert(C['unit_idx'],k)
        return coords
    return None

def enum_exceptional(C,B=8):
    D=C['D']; gens=C['letters']; d=len(gens)
    out={}
    zero=QD(0,0,D)
    for exps in itertools.product(range(-B,B+1),repeat=d):
        if all(e==0 for e in exps): continue
        xp=QD(1,0,D)
        for g,e in zip(gens,exps):
            if e: xp=xp*(g**e)
        for sgn in (1,-1):
            x=xp if sgn==1 else -xp
            ex=float(x.emb(60))
            if not (abs(ex)<1.0): continue
            y=1-x
            f=factor_sunit(y,C)
            if f is None: continue
            out[x]=dict(x=x,vx=list(exps),v1m=f,h=max(abs(e) for e in exps))
    res=sorted(out.values(), key=lambda r:(r['h'], r['vx']))
    return res

# ---------------- characters and L ----------------
def chi(D,a):
    if D==8: return {1:1,3:-1,5:-1,7:1}.get(a%8,0)
    if D==5: return {1:1,2:-1,3:-1,4:1}.get(a%5,0)
    if D==12: return {1:1,5:-1,7:-1,11:1}.get(a%12,0)
    if D==13:
        if a%13==0: return 0
        r=pow(a%13,6,13); return 1 if r==1 else -1
    raise ValueError(D)

def Lchi_hurwitz(D,P):
    old=mp.mp.dps; mp.mp.dps=P+25
    s=mp.mpf(0)
    for a in range(1,D):
        c=chi(D,a)
        if c: s+=c*mp.zeta(5,mp.mpf(a)/D)
    v=s/mp.power(D,5)
    mp.mp.dps=old
    return +v

def Lchi_nsum(D,P):
    old=mp.mp.dps; mp.mp.dps=P+25
    f=lambda n: chi(D,n)/mp.mpf(n)**5
    v=mp.nsum(f,[1,mp.inf],tol=mp.mpf(10)**(-P-10))
    mp.mp.dps=old
    return +v

def li5_series(x,P):
    """direct series sum x^k/k^5, |x|<1; guard 25 digits."""
    old=mp.mp.dps; mp.mp.dps=P+25
    s=mp.mpf(0); k=1; xk=x
    thresh=mp.mpf(10)**(-(P+15))
    ax=abs(x)
    while True:
        t=xk/mp.power(k,5)
        s+=t
        if abs(t)<thresh and k>30:
            # tail bound: |sum_{j>k}| <= |t| * something/(1-|x|)
            break
        xk*=x; k+=1
        if k>2000000: raise RuntimeError("series too slow")
    mp.mp.dps=old
    return +s

# ---------------- symbols ----------------
def tidx(t,d):
    i=0
    for x in t: i=i*d+x
    return i

def symbol_li5_dense(vx,v1m,d):
    """symbol of Li5(x): -(1-x) otimes x^4 ; dense list length d^5 (ints)."""
    vec=[0]*(d**5)
    for i0 in range(d):
        c0=-v1m[i0]
        if c0==0: continue
        base0=i0*(d**4)
        for t in itertools.product(range(d),repeat=4):
            cc=c0
            for j in t: cc*=vx[j]
            if cc: vec[base0+tidx(t,d)]+=cc
    return vec

def sym5_basis(d):
    """rows: symmetrised basis tensors of Sym^5, one per multiset."""
    rows=[]; keys=[]
    for ms in itertools.combinations_with_replacement(range(d),5):
        v=[0]*(d**5)
        for p in set(itertools.permutations(ms)):
            v[tidx(p,d)]+=1
        rows.append(v); keys.append(ms)
    return rows,keys

def prod_symbol_dense(ms,d):
    v=[0]*(d**5)
    for p in set(itertools.permutations(ms)):
        v[tidx(p,d)]+=1
    return v

def weight5_products(names):
    d=len(names); out=[]
    for ms in itertools.combinations_with_replacement(range(d),5):
        out.append(("log5",ms,"*".join(names[i] for i in ms)))
    for ms in itertools.combinations_with_replacement(range(d),3):
        out.append(("pi2log3",ms,"*".join([names[i] for i in ms]+["pi2"])))
    for i in range(d):
        out.append(("pi4log",(i,),"*".join([names[i],"pi2","pi2"])))
    return out

# ---------------- exact integer linear algebra ----------------
def _primitive(row):
    g=0
    for x in row: g=math.gcd(g,abs(int(x)))
    if g>1: row=[x//g for x in row]
    return row

def saturate_rows(B):
    """B: list of integer rows (k x n), rank k. Return saturated basis (span_Q(B) cap Z^n).
    Certificate: SNF invariants all 1 (checked by caller via snf_invariants)."""
    if not B: return B
    B=[_primitive(list(map(int,r))) for r in B]
    k=len(B); n=len(B[0])
    BT=[list(col) for col in zip(*B)]  # n x k
    H=flint.fmpz_mat(BT).hnf()
    H0=sp.Matrix(k,k,lambda i,j: sp.Integer(H[i,j]))
    for i in range(k,n):
        for j in range(k):
            assert H[i,j]==0, "hnf bottom not zero"
    Bm=sp.Matrix(B)
    Bsat=H0.T.inv()*Bm
    rows=[]
    for r in Bsat.tolist():
        rr=[sp.nsimplify(x) for x in r]
        assert all(x.q==1 for x in rr), "saturation not integral?!"
        rows.append([int(x) for x in rr])
    return rows

def snf_invariants(B):
    if not B: return []
    M=flint.fmpz_mat([list(map(int,r)) for r in B])
    S=M.snf()
    n=min(S.nrows(),S.ncols())
    return [int(S[i,i]) for i in range(n)]

def int_left_kernel(A):
    """A: list of integer rows (n x r). Return saturated Z-basis of {c in Z^n : c A = 0}."""
    Am=sp.Matrix(A)
    ns=Am.T.nullspace()
    B=[]
    for v in ns:
        den=1
        for x in v:
            x=sp.nsimplify(x); den=sp.ilcm(den,int(x.q))
        row=[int(sp.nsimplify(x)*den) for x in v]
        B.append(_primitive(row))
    if not B: return []
    return saturate_rows(B)

# ---------------- cobracket (3,2) ----------------
def wedge3(u,v,w,d):
    out={}
    for p in range(d):
        for q in range(p+1,d):
            for r in range(q+1,d):
                s=0
                for perm,sgn in [((p,q,r),1),((p,r,q),-1),((q,p,r),-1),((q,r,p),1),((r,p,q),1),((r,q,p),-1)]:
                    s+=sgn*u[perm[0]]*v[perm[1]]*w[perm[2]]
                if s: out[(p,q,r)]=out.get((p,q,r),0)+s
    return out

def wedge2(u,v,d):
    out={}
    for p in range(d):
        for q in range(p+1,d):
            s=u[p]*v[q]-u[q]*v[p]
            if s: out[(p,q)]=out.get((p,q),0)+s
    return out

def cobracket32_sparse(T,d,mode="full"):
    """T: dict (i0..i4)->coeff. Returns dict ((tri),(pair))->coeff in Lambda^3 (x) Lambda^2.
    mode 'full': sum over all C(5,3) index choices with shuffle signs; 'cut': only (0,1,2)|(3,4)."""
    from itertools import combinations
    out={}
    for tup,c in T.items():
        a=list(tup)
        if mode=="cut":
            choices=[((0,1,2),(3,4),1)]
        else:
            choices=[]
            for I in combinations(range(5),3):
                J=tuple(j for j in range(5) if j not in I)
                perm=list(I)+list(J)
                inv=sum(1 for x in range(5) for y in range(x+1,5) if perm[x]>perm[y])
                sgn=-1 if inv%2 else 1
                choices.append((tuple(I),J,sgn))
        for I,J,sgn in choices:
            u=[1 if a[I[0]]==i else 0 for i in range(d)]  # NOT correct; need letter vectors
        # placeholder (implemented properly in pipeline with letter-coordinate vectors)
    raise NotImplementedError
'''
open('w5lib.py','w').write(lib)
r = subprocess.run(['nice','-n','10',sys.executable,'-c',
    'import w5lib; print("import ok"); C=w5lib.letters_info(8); print(C["names"], C["Sdesc"]); '
    'print("hnf doc:", (w5lib.flint.fmpz_mat.hnf.__doc__ or "")[:200]); '
    'print("snf doc:", (w5lib.flint.fmpz_mat.snf.__doc__ or "")[:120])'],
    capture_output=True,text=True,timeout=300)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
with open('progress.log','a') as f: f.write("step_05_lib: w5lib.py written, import ok\n")
