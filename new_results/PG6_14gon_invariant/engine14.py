
from fractions import Fraction as Fr
from math import lcm, comb
import pickle
# ---------- fast FT ----------
def ft(A=0,B=0,C=0,d=1): return (Fr(A),Fr(B),Fr(C),int(d))
ZERO=(Fr(0),)*3+(1,); ONE=(Fr(1),Fr(0),Fr(0),1); AL=(Fr(0),Fr(1),Fr(0),1)
def fadd(x,y):
    A1,B1,C1,d1=x; A2,B2,C2,d2=y; d=lcm(d1,d2)
    return (A1*(d//d1)+A2*(d//d2), B1*(d//d1)+B2*(d//d2), C1*(d//d1)+C2*(d//d2), d)
def fneg(x): A,B,C,d=x; return (-A,-B,-C,d)
def fsub(x,y): return fadd(x,fneg(y))
def fmul(x,y):
    A1,B1,C1,d1=x; A2,B2,C2,d2=y
    A=A1*A2+(B1*C2+C1*B2)-C1*C2
    B=A1*B2+B1*A2+2*(B1*C2+C1*B2)-C1*C2
    C=A1*C2+C1*A2+B1*B2-(B1*C2+C1*B2)+3*C1*C2
    return (A,B,C,d1*d2)
def fs(x,q):
    A,B,C,d=x; q=Fr(q); return (A*q.numerator,B*q.numerator,C*q.numerator,d*q.denominator)
def finv(x):
    A,B,C,d=x
    M=[[A,C,B-C],[B,A+2*C,2*B-C],[C,B-C,A-B+3*C]]
    import sympy as sp
    Ms=sp.Matrix([[sp.Rational(v) for v in r] for r in M]); col=Ms.inv()*sp.Matrix([sp.Rational(d),0,0])
    return (Fr(col[0]),Fr(col[1]),Fr(col[2]),1)
def fdiv(x,y): return fmul(x,finv(y))
def fpow(x,k):
    if k<0: return fpow(finv(x),-k)
    r=ONE
    for _ in range(k): r=fmul(r,x)
    return r
def feq0(x): return x[0]==0 and x[1]==0 and x[2]==0
import mpmath as mp
def fev(x,prec=60):
    old=mp.mp.dps; mp.mp.dps=max(prec,30)
    A,B,C,d=x; a=2*mp.cos(2*mp.pi/7)
    v=(mp.mpf(int(A.numerator))/A.denominator+mp.mpf(int(B.numerator))/B.denominator*a+mp.mpf(int(C.numerator))/C.denominator*a*a)/d
    mp.mp.dps=old; return v
# ---------- Expr: dict monomial->FT ; monomial=tuple((sym,exp),...) sorted ----------
def mkey(d): return tuple(sorted(d.items()))
def mdot(m1,m2):
    d=dict(m1)
    for s,e in m2:
        d[s]=d.get(s,0)+e
        if d[s]==0: del d[s]
    return mkey(d)
def e_const(c): return {():c} if not feq0(c) else {}
def e_sym(name): return {((name,1),):ONE}
def e_add(P,Q):
    R=dict(P)
    for m,c in Q.items():
        R[m]=fadd(R.get(m,ZERO),c)
        if feq0(R[m]): del R[m]
    return R
def e_addinto(R,m,c):
    R[m]=fadd(R.get(m,ZERO),c)
    if feq0(R[m]): del R[m]
def e_neg(P): return {m:fneg(c) for m,c in P.items()}
def e_sub(P,Q): return e_add(P,e_neg(Q))
def e_scale(P,c):  # c is FT
    if feq0(c): return {}
    return {m:fmul(co,c) for m,co in P.items()}
def e_scaleq(P,q): # rational
    if q==0: return {}
    return {m:fs(co,q) for m,co in P.items()}
def e_mul(P,Q):
    R={}
    for m1,c1 in P.items():
        for m2,c2 in Q.items():
            m=mdot(m1,m2); e_addinto(R,m,fmul(c1,c2))
    return R
def e_pow(P,k):
    R=e_const(ONE)
    for _ in range(k): R=e_mul(R,P)
    return R
# ---------- geometry data ----------
import pickle as pk
_g = pk.load(open('geom14.pkl','rb'))
X = {i:_g['X'][i] for i in [1,2,5,6]}
Km = {i:fneg(X[i]) for i in X}
invX = {i:finv(X[i]) for i in X}
inv1mX = {i:finv(fsub(ONE,X[i])) for i in X}
IDX=[1,2,5,6]
LAN={1:6,6:1,2:5,5:2}
# LOG1MX = log(1-x_i) as letter Expr : {1: la, 2: lu, 5: -lu, 6: -la}
LOG1MX = {1:e_sym('la'), 2:e_sym('lu'), 5:e_neg(e_sym('lu')), 6:e_neg(e_sym('la'))}
L1 = {i:e_neg(LOG1MX[i]) for i in IDX}
# atoms
L2x={i:e_sym('L2x%d'%i) for i in IDX}
L3x={i:e_sym('L3x%d'%i) for i in IDX}
Ax ={i:e_sym('Ax%d'%i) for i in IDX}
Lam={}; T2s={}; T3ps={}
for i in IDX:
    for j in IDX:
        if i==j: continue
        Lam[(i,j)]=e_sym('Lam%d_%d'%(i,j)); T2s[(i,j)]=e_sym('T2_%d_%d'%(i,j)); T3ps[(i,j)]=e_sym('T3p_%d_%d'%(i,j))
T3a={}
for i in IDX:
    for j in IDX:
        for k in IDX:
            if len({i,j,k})==3 and j<k: T3a[(i,j,k)]=e_sym('T3_%d_%d%d'%(i,j,k))
z2=e_sym('z2'); z3=e_sym('z3')
# ---------- Jloc ----------
import sympy as sp
_jloc = pk.load(open('Jloc.pkl','rb'))
def _parse_jloc(e):
    # e = c0 + c2*pi^2 + c3*zeta(3) -> dict {():c0, (('z2',1),):6*c2, (('z3',1),):c3}
    e=sp.expand(e); c0=e.subs({sp.pi:0}); 
    c2=e.coeff(sp.pi,2); c3=e.coeff(sp.zeta(3),1)
    base=e - c2*sp.pi**2 - c3*sp.zeta(3)
    base=sp.nsimplify(base, rational=False) if base.free_symbols else base
    R={}
    def put(m,c):
        c=Fr(sp.Rational(c))
        if c: R[m]=(c,Fr(0),Fr(0),1)
    put((), base)
    if c2: e_addinto(R, (('z2',1),), (Fr(6)*Fr(sp.Rational(c2)),Fr(0),Fr(0),1))
    if c3: e_addinto(R, (('z3',1),), (Fr(sp.Rational(c3)),Fr(0),Fr(0),1))
    return R
JLOCtab={k:_parse_jloc(v) for k,v in _jloc.items()}
def JLOC(n,a,b): return JLOCtab[(n,a,b)]
# ---------- Jtab -> K1 ----------
_jtab = pk.load(open('Jtab_t.pkl','rb'))
_sx,_sL1,_sL2,_sL3,_sAA,_sz2 = sp.symbols('x L1 L2 L3 AA z2')
def _parse_jtab(e):
    # Laurent poly in x, coeffs poly in L1,L2,L3,AA,z2 (rational)
    e=sp.expand(e); R=[]
    for term in sp.Add.make_args(e):
        d={'x':0,'L1':0,'L2':0,'L3':0,'AA':0,'z2':0}; coef=sp.Integer(1)
        for f in sp.Mul.make_args(term):
            if f.is_Symbol:
                d[str(f)]+=1
            elif f.is_Pow:
                b_,e_=f.as_base_exp()
                if b_.is_Symbol and str(b_) in d: d[str(b_)]+=int(e_)
                else: coef*=f
            else: coef*=f
        R.append((d, Fr(sp.Rational(coef))))
    return R
JTAB={k:_parse_jtab(v) for k,v in _jtab.items()}
def K1(n,a,b,c,side):
    i=side; R={}
    for (d,coef) in JTAB[(n,a,b,c)]:
        # x^e -> FT X[i]^e ; L1^a -> letter expr L1[i]^a ; L2->atom, L3->atom, AA->atom, z2->z2
        m={}
        if d['L2']: m['L2x%d'%i]=d['L2']
        if d['L3']: m['L3x%d'%i]=d['L3']
        if d['AA']: m['Ax%d'%i]=d['AA']
        if d['z2']: m['z2']=d['z2']
        base=e_const(fpow(X[i], d['x']))
        if d['L1']:
            base=e_mul(base, e_pow(L1[i], d['L1']))
        mon=mkey(m)
        term={mon:()}  # placeholder
        # build: coefficient base (FT) times monomial mon
        cc={}
        for mm,ccft in base.items():
            # mm should be letter monomial; combine with atom monomial mon
            e_addinto(R, mdot(mm,mon), fs(ccft, coef))
    return R

# ---------- pair / triple recursions ----------
g = {i:invX[i] for i in IDX}          # (-1/k_i) = 1/X_i
ok_ = {i:fneg(invX[i]) for i in IDX}  # 1/k_i
kval = {i:Km[i] for i in IDX}         # k_i = -X_i
def D1f(i,j): return e_scale(Lam[(i,j)], invX[i])
def D3f(i,j):
    ip,jp = LAN[i], LAN[j]
    t1 = e_mul(LOG1MX[j], e_scale(e_neg(L2x[ip]), invX[ip]))
    t2 = e_scale(Lam[(ip,jp)], invX[ip])
    return e_scale(e_add(t1,t2), inv1mX[i])
_pairs={}
def make_pair(i,j):
    D1=D1f(i,j); D2=D1f(j,i); D3=D3f(i,j); D4=D3f(j,i)
    D5=T2s[(i,j)]; D6=T2s[(j,i)]; D7=T3ps[(j,i)]; D8=T3ps[(i,j)]
    # D9 = (LOG1MX[i]^2 LOG1MX[j])/(2 k1) - (k4/(2 k1)) D7 ; 1/(2k1)= -invX[i]/2 ; k4/(2k1)= X[j] invX[i]/2
    t1 = e_scale(e_mul(e_pow(LOG1MX[i],2),LOG1MX[j]), fs(ok_[i],Fr(1,2)))
    t2 = e_scale(D7, fs(fdiv(Km[j],Km[i]), Fr(1,2)))
    D9 = e_sub(t1,t2)
    t1b= e_scale(e_mul(LOG1MX[i],e_pow(LOG1MX[j],2)), fs(ok_[j],Fr(1,2)))
    t2b= e_scale(D8, fs(fdiv(Km[i],Km[j]), Fr(1,2)))
    D10= e_sub(t1b,t2b)
    cache={}
    def Geo(ii, m):  # g[ii]^m
        return fpow(g[ii], m)
    def MixedW(n):
        if ('W',n) in cache: return cache[('W',n)]
        bnd=e_mul(LOG1MX[i],LOG1MX[j])
        A1=e_const(ZERO)
        for m in range(n+1):
            A1=e_add(A1, e_scale(K1(n-m,0,0,1,j), Geo(i,m)))
        A1=e_scale(A1, ok_[i]); A1=e_add(A1, e_scale(D5, Geo(i,n+1)))
        A2=e_const(ZERO)
        for m in range(n+1):
            A2=e_add(A2, e_scale(K1(n-m,0,0,1,i), Geo(j,m)))
        A2=e_scale(A2, ok_[j]); A2=e_add(A2, e_scale(D6, Geo(j,n+1)))
        r=e_sub(bnd, e_add(e_scale(A1,kval[i]), e_scale(A2,kval[j])))
        r=e_scaleq(r, Fr(1,n+1)); cache[('W',n)]=r; return r
    def M1(n):
        if ('M1',n) in cache: return cache[('M1',n)]
        B1=e_const(ZERO)
        for m in range(n+1): B1=e_add(B1, e_scale(K1(n-m,1,0,1,j), Geo(i,m)))
        B1=e_scale(B1, ok_[i]); B1=e_add(B1, e_scale(D1, Geo(i,n+1)))
        B2=e_const(ZERO)
        for m in range(n+1): B2=e_add(B2, e_scale(K1(n-m,1,0,1,i), Geo(j,m)))
        B2=e_scale(B2, ok_[j]); B2=e_add(B2, e_scale(D2, Geo(j,n+1)))
        r=e_neg(e_scaleq(e_add(MixedW(n), e_add(e_scale(B1,kval[i]), e_scale(B2,kval[j]))), Fr(1,n+1)))
        cache[('M1',n)]=r; return r
    def M2(n):
        if ('M2',n) in cache: return cache[('M2',n)]
        r1=fsub(Geo(i,n+1), ONE); r4=fsub(Geo(j,n+1), ONE)
        C1=e_const(ZERO)
        for m in range(n+1): C1=e_add(C1, e_scale(K1(n-m,0,1,1,j), Geo(i,m)))
        C1=e_scale(C1, ok_[i]); C1=e_add(C1, e_scale(D3, r1))
        C2=e_const(ZERO)
        for m in range(n+1): C2=e_add(C2, e_scale(K1(n-m,0,1,1,i), Geo(j,m)))
        C2=e_scale(C2, ok_[j]); C2=e_add(C2, e_scale(D4, r4))
        sW=e_const(ZERO)
        for m in range(n+1): sW=e_add(sW, MixedW(m))
        r=e_scaleq(e_sub(e_neg(sW), e_add(e_scale(C1,kval[i]), e_scale(C2,kval[j]))), Fr(1,n+1))
        cache[('M2',n)]=r; return r
    def M3(n):
        if ('M3',n) in cache: return cache[('M3',n)]
        E1=e_const(ZERO)
        for m in range(n+1): E1=e_add(E1, e_scale(MixedW(n-m), Geo(i,m)))
        E1=e_scale(E1, ok_[i]); E1=e_add(E1, e_scale(D9, Geo(i,n+1)))
        E2=e_const(ZERO)
        for m in range(n+1): E2=e_add(E2, e_scale(K1(n-m,0,0,2,i), Geo(j,m)))
        E2=e_scale(E2, ok_[j]); E2=e_add(E2, e_scale(D7, Geo(j,n+1)))
        bnd=e_mul(e_pow(LOG1MX[i],2),LOG1MX[j])
        r=e_scaleq(e_sub(bnd, e_add(e_scale(E1,fs(kval[i],2)), e_scale(E2,kval[j]))), Fr(1,n+1))
        cache[('M3',n)]=r; return r
    def M4(n):
        if ('M4',n) in cache: return cache[('M4',n)]
        F1=e_const(ZERO)
        for m in range(n+1): F1=e_add(F1, e_scale(K1(n-m,0,0,2,j), Geo(i,m)))
        F1=e_scale(F1, ok_[i]); F1=e_add(F1, e_scale(D8, Geo(i,n+1)))
        F2=e_const(ZERO)
        for m in range(n+1): F2=e_add(F2, e_scale(MixedW(n-m), Geo(j,m)))
        F2=e_scale(F2, ok_[j]); F2=e_add(F2, e_scale(D10, Geo(j,n+1)))
        bnd=e_mul(LOG1MX[i],e_pow(LOG1MX[j],2))
        r=e_scaleq(e_sub(bnd, e_add(e_scale(F1,kval[i]), e_scale(F2,fs(kval[j],2)))), Fr(1,n+1))
        cache[('M4',n)]=r; return r
    def MIX(n,a,b,ci,cj):
        if a==1: return M1(n)
        if b==1: return M2(n)
        if (ci,cj)==(2,1): return M3(n)
        if (ci,cj)==(1,2): return M4(n)
        return MixedW(n)
    return MIX, MixedW
def PAIR(n,a,b,ci,cj,i,j):
    if (i,j) not in _pairs: _pairs[(i,j)]=make_pair(i,j)
    return _pairs[(i,j)][0](n,a,b,ci,cj)
def PAIR_MixedW(n,i,j):
    if (i,j) not in _pairs: _pairs[(i,j)]=make_pair(i,j)
    return _pairs[(i,j)][1](n)
_trip={}
def TRIPLE(n,i,j,k):
    if (n,i,j,k) in _trip: return _trip[(n,i,j,k)]
    bnd=e_mul(e_mul(LOG1MX[i],LOG1MX[j]),LOG1MX[k])
    res=bnd
    for (ii,jj,kk) in [(i,j,k),(j,i,k),(k,i,j)]:
        poly=e_const(ZERO)
        for m in range(n+1):
            poly=e_add(poly, e_scale(PAIR_MixedW(n-m, jj, kk), fpow(g[ii],m)))
        key=(ii,)+tuple(sorted((jj,kk)))
        res=e_sub(res, poly)
        res=e_sub(res, e_scale(T3a[key], fmul(kval[ii], fpow(g[ii],n+1))))
    res=e_scaleq(res, Fr(1,n+1)); _trip[(n,i,j,k)]=res; return res
def INT(n,a,b,wc):
    fams=[idx for idx,c in zip(IDX,wc) if c>0]
    if len(fams)==0: return JLOC(n,a,b)
    if len(fams)==1:
        i=fams[0]; return K1(n,a,b,wc[IDX.index(i)],i)
    if len(fams)==2:
        i,j=fams; ci,cj=wc[IDX.index(i)],wc[IDX.index(j)]
        return PAIR(n,a,b,ci,cj,i,j)
    if len(fams)==3:
        i,j,k=fams
        assert a==0 and b==0 and all(wc[IDX.index(f)]==1 for f in fams)
        return TRIPLE(n,i,j,k)
    raise ValueError('quad')
