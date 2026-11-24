
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
    def __rtruediv__(s,o):
        if isinstance(o,(int,Fraction)): return QD(o,0,s.D)*s.inv()
        return o*s.inv()
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
        return v
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
        w=QD(0,1,2)
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
        r3=QD(0,1,3); eps=2+r3; u=1+r3
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
    elif D==17:
        r=QD(0,1,17); eps=QD(4,1,17)
        pi1=QD(Fraction(3,2),Fraction(1,2),17); pi2_=QD(Fraction(3,2),Fraction(-1,2),17)
        letters=[eps,r,pi1,pi2_]; names=["lE17","ls17","lp","lp'"]; disp=["log(4+sqrt17)","log(sqrt17)","log((3+sqrt17)/2)","log((3-sqrt17)/2)"]
        unit_idx=0; unit=eps; nonunit=[r,pi1,pi2_]
        def pexps(N):
            return None  # unused; factor_override handles D=17
    elif D==21:
        r=QD(0,1,21); eps=QD(Fraction(5,2),Fraction(1,2),21)
        al=QD(Fraction(3,2),Fraction(1,2),21)
        letters=[eps,r,al,QD(2,0,21)]; names=["lE21","ls21","la","l2"]; disp=["log((5+sqrt21)/2)","log(sqrt21)","log((3+sqrt21)/2)","log(2)"]
        unit_idx=0; unit=eps; nonunit=[r,al,QD(2,0,21)]
        def pexps(N):
            if not _allowed(N,[2,3,7]): return None
            v2=_v_p(N,2)
            if v2%2!=0: return None
            v7=_v_p(N,7); v3=_v_p(N,3)
            return [v7, v3-v7, v2//2]
    elif D==29:
        r=QD(0,1,29); eps=QD(Fraction(5,2),Fraction(1,2),29)
        letters=[eps,r,QD(2,0,29)]; names=["lE29","ls29","l2"]; disp=["log((5+sqrt29)/2)","log(sqrt29)","log(2)"]
        unit_idx=0; unit=eps; nonunit=[r,QD(2,0,29)]
        def pexps(N):
            if not _allowed(N,[2,29]): return None
            v2=_v_p(N,2)
            if v2%2!=0: return None
            return [_v_p(N,29), v2//2]
    else:
        raise ValueError(D)
    Sdesc={8:"S={2} (ramified; 2=(sqrt2)^2)",
           5:"S={5,2} (5=(sqrt5)^2 ramified, 2 inert)",
           12:"S={3,2} (3=(sqrt3)^2, 2=(1+sqrt3)^2 both ramified)",
           13:"S={13,2} (13=(sqrt13)^2 ramified, 2 inert)",
           29:"S={29,2} (29=(sqrt29)^2 ramified, 2 inert)",
           17:"S={17,2} (17=(sqrt17)^2 ramified, 2 split: ((3+-sqrt17)/2))",
           21:"S={3,7,2} (3,7 ramified via sqrt21 and (3+sqrt21)/2, 2 inert)"}[D]
    return dict(D=D,letters=letters,names=names,disp=disp,unit_idx=unit_idx,
                unit=unit,nonunit=nonunit,pexps=pexps,Sdesc=Sdesc,
                unit_emb=unit.emb(150))

def _integral_half(x):
    # integrality in Q(sqrt d) with d = 1 mod 4: x = (a-b) + (2b) omega, omega=(1+sqrt d)/2
    return (x.a - x.b).denominator == 1 and (2*x.b).denominator == 1

def _val_at(x, pi):
    k=0; y=x
    while _integral_half(y/pi): y=y/pi; k+=1
    while _integral_half((1/y)/pi): y=y*pi; k-=1
    return k,y

def factor_sunit(x,C,verb=False):
    """Return coords (aligned with C['letters']) of |x| if x is an S-unit, else None.
    (Sign is ignored: -1 is torsion.)"""
    if x.a==0 and x.b==0: return None
    N=x.norm()
    if C['D']==17:
        if not _allowed(N,[2,17]): return None
        r=QD(0,1,17); pi1=QD(Fraction(3,2),Fraction(1,2),17); pi2_=QD(Fraction(3,2),Fraction(-1,2),17)
        two=QD(2,0,17)
        v17=_v_p(N,17); S=_v_p(N,2)
        y=x/(r**v17)
        m=0
        while not _integral_half(y*(two**m)): m+=1
        y=y*(two**m)
        c1=0
        while _integral_half(y/pi1): y=y/pi1; c1+=1
        e1=c1-m
        e2=S-c1+m
        y=y/(pi2_**(S+2*m-c1))
        if y.norm() not in (1,-1): return None
        old=mp.mp.dps; mp.mp.dps=150
        t=mp.log(abs(y.emb(150)))/mp.log(C['unit_emb'])
        k=int(mp.nint(t))
        ok = abs(t-k)<mp.mpf('1e-80')
        mp.mp.dps=old
        if not ok: return None
        eps=C['unit']
        if y==eps**k or y==-(eps**k):
            return [k,v17,e1,e2]
        return None
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
    rad=gens[0].D
    for exps in itertools.product(range(-B,B+1),repeat=d):
        if all(e==0 for e in exps): continue
        xp=QD(1,0,rad)
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
    if D in (13,29):
        if a%D==0: return 0
        r=pow(a%D,(D-1)//2,D); return 1 if r==1 else -1
    if D==17:
        if a%17==0: return 0
        r=pow(a%17,8,17); return 1 if r==1 else -1
    if D==21:
        if a%3==0 or a%7==0: return 0
        l3 = 1 if pow(a%3,1,3)==1 else -1
        l7 = 1 if pow(a%7,3,7)==1 else -1
        return l3*l7
    raise ValueError(D)

def Lchi_hurwitz(D,P):
    old=mp.mp.dps; mp.mp.dps=P+25
    s=mp.mpf(0)
    for a in range(1,D):
        c=chi(D,a)
        if c: s+=c*mp.zeta(5,mp.mpf(a)/D)
    v=s/mp.power(D,5)
    mp.mp.dps=old
    return v

def Lchi_nsum(D,P):
    old=mp.mp.dps; mp.mp.dps=P+25
    f=lambda n: chi(D,n)/mp.mpf(n)**5
    v=mp.nsum(f,[1,mp.inf],tol=mp.mpf(10)**(-P-10))
    mp.mp.dps=old
    return v

def li5_jonquiere(x,P):
    """Li5(x) for 0<x<1 near 1 via Jonquiere's expansion, t=-log(x), |t|<2pi."""
    old=mp.mp.dps; mp.mp.dps=P+30
    x=mp.mpf(x)
    t=-mp.log(x)
    H4=mp.mpf(25)/12
    s=t**4/24*(H4-mp.log(t))
    # sum over k>=0, k!=4: zeta(5-k)(-t)^k/k!   (fixed count; terms with zeta(even neg)=0 are skipped)
    for k in range(0, 320):
        if k==4: continue
        zk=mp.zeta(5-k)
        if zk:
            s+=zk*(-t)**k/mp.factorial(k)
    mp.mp.dps=old
    return s

def li5_smart(x,P):
    x=mp.mpf(x)
    if mp.mpf('0.9') < x < 1:
        return li5_jonquiere(x,P)
    return li5_series(x,P)

def li5_series(x,P):
    """direct series sum x^k/k^5, |x|<1; guard 25 digits."""
    old=mp.mp.dps; mp.mp.dps=P+25
    x=mp.mpf(x)
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
    return s

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
    # shuffle (symbol-of-product) convention: each distinct permutation p of ms
    # occurs prod_i mult_i! times (multiplicities in the shuffle of 1-letter words)
    from math import factorial
    cnt = {}
    for i in ms: cnt[i] = cnt.get(i, 0) + 1
    mult = 1
    for k in cnt.values(): mult *= factorial(k)
    v=[0]*(d**5)
    for p in set(itertools.permutations(ms)):
        v[tidx(p,d)] += mult
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

def int_left_kernel_hnf(A):
    """Saturated Z-basis of {c : c A = 0} via HNF with transformation (fast, flint)."""
    import flint as _fl
    M = _fl.fmpz_mat([list(map(int, r)) for r in A])
    H, U = M.hnf(transform=True)
    nr = M.nrows(); nc = M.ncols()
    rank = 0
    for i in range(nr):
        if any(H[i, j] for j in range(nc)): rank += 1
    K = []
    for i in range(rank, nr):
        K.append([int(U[i, j]) for j in range(nc if False else U.ncols())])
    # U is (nr x nr); rows rank..nr-1 satisfy u*M == 0
    return K

def int_left_kernel(A):
    """A: list of integer rows (n x r). Return saturated Z-basis of {c in Z^n : c A = 0}."""
    Am=sp.Matrix(A)
    ns=Am.T.nullspace()
    B=[]
    for v in ns:
        den=1
        for x in v:
            x=sp.Rational(x); den=sp.ilcm(den,int(x.q))
        row=[int(sp.Rational(x)*den) for x in v]
        B.append(_primitive(row))
    if not B: return []
    return saturate_rows(B)

# ---------------- cobracket (3,2) ----------------
def _perm_sign(seq):
    """sign of the permutation that sorts seq (seq must have distinct entries)."""
    inv = 0
    for i in range(len(seq)):
        for j in range(i+1, len(seq)):
            if seq[i] > seq[j]: inv += 1
    return -1 if inv % 2 else 1

def cobracket32(T, d, mode="full"):
    """T: dict (i0..i4) -> coeff (letters = formal basis of V, dim d).
    Returns dict ((triple),(pair)) -> coeff for delta_{3,2}: V^ot5 -> Lambda^3 V (x) Lambda^2 V.
    mode 'full': sum over all C(5,3)=10 (3+2)-shuffles with signs (Lie cobracket);
    mode 'cut' : plain deconcatenation (0,1,2)|(3,4) followed by antisymmetrisation."""
    out = {}
    for tup, c in T.items():
        if c == 0: continue
        if mode == "cut":
            splits = [((0,1,2),(3,4),1)]
        else:
            splits = []
            for I in itertools.combinations(range(5),3):
                J = tuple(j for j in range(5) if j not in I)
                sgn = _perm_sign(list(I)+list(J))
                splits.append((I,J,sgn))
        for I,J,sgn in splits:
            tri = [tup[I[0]], tup[I[1]], tup[I[2]]]
            if len(set(tri)) < 3: continue
            pair = [tup[J[0]], tup[J[1]]]
            if pair[0] == pair[1]: continue
            tsgn = _perm_sign(tri)
            psgn = 1 if pair[0] < pair[1] else -1
            key = (tuple(sorted(tri)), tuple(sorted(pair)))
            out[key] = out.get(key, 0) + c*sgn*tsgn*psgn
    return {k:v for k,v in out.items() if v}
