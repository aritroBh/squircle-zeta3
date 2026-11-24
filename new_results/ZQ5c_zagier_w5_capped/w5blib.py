
# w5blib.py -- ZQ5b: weight-5/weight-3 formulas for L(chi_D,k) with enlarged S.
# Based on ZQ5's w5lib.py (QD arithmetic, saturation, LLL) but:
#   * planner_li5 Li5 evaluator (correct at all x; ZQ5's li5_jonquiere was buggy for 0.94<x<0.995)
#   * general S-unit machinery: ramified + inert + SPLIT primes (division-based valuations)
#   * exact symbol quotient via orbit-difference functionals (no giant Wann matrix)
#   * companions P_j: pure-log products with symbol(P_j) = symbol(K_j)
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

def fmt_qd(x):
    d=x.D
    def fs(f): return str(f)
    if x.b==0: return f"({fs(x.a)})"
    if x.a==0: return f"({fs(x.b)}*sqrt{d})"
    sgn="+" if x.b>0 else "-"
    return f"({fs(x.a)}{sgn}{fs(abs(x.b))}*sqrt{d})"

def _v_p(fr,p):
    fr=Fraction(fr); v=0
    num,den=abs(fr.numerator),fr.denominator
    while num%p==0 and num>0: num//=p; v+=1
    while den%p==0: den//=p; v-=1
    return v

def _allowed(fr,primes):
    fr=Fraction(fr)
    remn=abs(fr.numerator); remd=fr.denominator
    for p in primes:
        while remn%p==0 and remn>0: remn//=p
        while remd%p==0: remd//=p
    return remn==1 and remd==1

def _integral_half(x):
    # O_K = Z[omega], omega=(1+sqrt d)/2 (d=1 mod 4): x = a+b*sqrt d = (a-b) + 2b*omega
    return (x.a - x.b).denominator == 1 and (2*x.b).denominator == 1

def _den_lcm(x):
    return math.lcm((x.a-x.b).denominator, (2*x.b).denominator)

def val_at_g(y, g, p, eram):
    """Exact valuation v_P(y) where P=(g), N(g)=+-p, ramification index eram (2 ram, 1 split).
    Returns (k, y/g^k)."""
    dn = _den_lcm(y)
    yn = y*dn   # integral
    k1 = 0
    while True:
        q = yn/g
        if _integral_half(q): yn = q; k1 += 1
        else: break
    k = k1 - eram*_v_p(dn, p)
    return k, y/(g**k)

# ---------------- field / S-unit data ----------------
def make_field(D, split=(), inert=(2,)):
    """D: squarefree disc (all our targets = 1 mod 4). Returns letter basis of O_S^*/torsion:
    eps (fund. unit), ramified gens (one per p|D, norm +-p), inert rationals, split gens (norm +-p pair)."""
    F = dict(D=D)
    r = D
    if D == 13:
        eps = QD(Fraction(3,2),Fraction(1,2),13); ram = [QD(0,1,13)]
        ramp = [13]
    elif D == 29:
        eps = QD(Fraction(5,2),Fraction(1,2),29); ram = [QD(0,1,29)]
        ramp = [29]
    elif D == 21:
        eps = QD(Fraction(5,2),Fraction(1,2),21)
        ram = [QD(Fraction(3,2),Fraction(1,2),21), QD(Fraction(7,2),Fraction(-1,2),21)]
        ramp = [3,7]
    elif D == 17:
        eps = QD(4,1,17); ram = [QD(0,1,17)]
        ramp = [17]
    elif D == 5:
        eps = QD(Fraction(1,2),Fraction(1,2),5); ram=[QD(0,1,5)]; ramp=[5]
    elif D == 8:
        eps = QD(1,1,2); ram=[QD(0,1,2)]; ramp=[2]
    else:
        raise ValueError(D)
    # split prime generators (hard-coded small solutions of a^2-Db^2 = +-4p), checked by caller
    splitgens = {
        (13,3): (QD(Fraction(5,2),Fraction(1,2),13), QD(Fraction(5,2),Fraction(-1,2),13)),
        (29,5): (QD(Fraction(7,2),Fraction(1,2),29), QD(Fraction(7,2),Fraction(-1,2),29)),
        (29,7): (QD(Fraction(1,2),Fraction(1,2),29), QD(Fraction(1,2),Fraction(-1,2),29)),
        (29,13): (QD(Fraction(9,2),Fraction(1,2),29), QD(Fraction(9,2),Fraction(-1,2),29)),
        (21,5): (QD(Fraction(1,2),Fraction(1,2),21), QD(Fraction(1,2),Fraction(-1,2),21)),
        (17,2): (QD(Fraction(3,2),Fraction(1,2),17), QD(Fraction(3,2),Fraction(-1,2),17)),
        (17,13): (QD(2,1,17), QD(2,-1,17)),
        (5,2): None, (8,2): None,
    }
    letters = [eps]; names = [f"lE{D}"]; kinds=['unit']; prime_of=[None]; erams=[None]
    F['unit_idx']=0; F['unit']=eps
    for g,p in zip(ram,ramp):
        assert g.norm() in (p,-p), (D,p,g.norm())
        letters.append(g); names.append(f"ls{D}" if p==D else f"lr{p}"); kinds.append('ram'); prime_of.append(p); erams.append(2)
    for p in inert:
        assert sp.kronecker_symbol(D,p)==-1, (D,p)
        letters.append(QD(p,0,r)); names.append(f"l{p}"); kinds.append('inert'); prime_of.append(p); erams.append(1)
    for p in split:
        assert sp.kronecker_symbol(D,p)==1, (D,p)
        g1,g2 = splitgens[(D,p)]
        assert g1.norm() in (p,-p) and g2.norm() in (p,-p), (D,p,g1.norm())
        assert g2 == g1.conj()
        letters.append(g1); names.append(f"lp{p}"); kinds.append('split'); prime_of.append(p); erams.append(1)
        letters.append(g2); names.append(f"lp{p}c"); kinds.append('split'); prime_of.append(p); erams.append(1)
    F['letters']=letters; F['names']=names; F['kinds']=kinds; F['prime_of']=prime_of; F['erams']=erams
    F['rat_primes']= sorted(set(ramp)|set(inert)|set(split))
    F['split']=list(split); F['inert']=list(inert); F['ram']=ramp
    old=mp.mp.dps; mp.mp.dps=170
    F['unit_emb']=abs(eps.emb(160))
    mp.mp.dps=old
    return F

def factor_sunit(x, F):
    """coords of |x| in F['letters'] if x is an S-unit, else None. Sign ignored."""
    if x.a==0 and x.b==0: return None
    N = x.norm()
    if not _allowed(N, F['rat_primes']): return None
    y = x
    coords = [None]*len(F['letters'])
    for i,(g,kind) in enumerate(zip(F['letters'],F['kinds'])):
        if kind in ('ram','split'):
            k,y = val_at_g(y, g, F['prime_of'][i], F['erams'][i])
            coords[i]=k
    for i,(g,kind) in enumerate(zip(F['letters'],F['kinds'])):
        if kind=='inert':
            p = F['prime_of'][i]
            v = _v_p(y.norm(), p)
            if v%2: return None
            y = y/(Fraction(p)**(v//2))
            coords[i]=v//2
    if y.norm() not in (1,-1): return None
    old=mp.mp.dps; mp.mp.dps=150
    t=mp.log(abs(y.emb(150)))/mp.log(F['unit_emb'])
    k=int(mp.nint(t))
    ok = abs(t-k)<mp.mpf('1e-80')
    mp.mp.dps=old
    if not ok: return None
    eps=F['unit']
    if y==eps**k or y==-(eps**k):
        coords[F['unit_idx']]=k
        return coords
    return None

def enum_exceptional(F, B=8, time_limit=900):
    """Exceptional S-units |x|<1, x and 1-x S-units, exponents in [-B,B]."""
    import numpy as np
    D=F['D']; gens=F['letters']; d=len(gens)
    Lf = np.array([float(mp.log(abs(g.emb(60)))) for g in gens])
    t0=time.time()
    axes = [np.arange(-B,B+1)]*d
    grids = np.meshgrid(*axes, indexing='ij')
    T = np.stack([g.ravel() for g in grids], axis=1).astype(float)
    logabs = T @ Lf
    mask = logabs < 1e-9   # |x| < 1 (with guard); excludes identity tuple (logabs=0 <1e-9 but x=1 handled)
    keep = T[mask].astype(int)
    # drop the all-zero tuple (x=1 excluded) and tuples giving |x|>=1 exactly (none)
    keep = keep[np.any(keep!=0, axis=1)]
    out={}
    n_checked=0
    for exps in keep:
        n_checked+=1
        if n_checked%200000==0 and time.time()-t0>time_limit:
            print(f"  enum TIME LIMIT at {n_checked}/{len(keep)}; stopping early")
            break
        xp=QD(1,0,D)
        for g,e in zip(gens,exps):
            if e: xp=xp*(g**int(e))
        # |xp|<1 exact check
        v = xp.emb(70)
        if not (abs(v) < 1): continue
        for sgn in (1,-1):
            x = xp if sgn==1 else -xp
            y = 1-x
            f = factor_sunit(y, F)
            if f is None: continue
            out[x]=dict(x=x, vx=[int(e) for e in exps], v1m=f, h=int(max(abs(e) for e in exps)))
    res=sorted(out.values(), key=lambda r:(r['h'], r['vx'], -r['x'].a))
    return res

# ---------------- characters and L ----------------
def chi(D,a):
    return int(sp.kronecker_symbol(D, a))

def Lchi_hurwitz(D,k,P):
    old=mp.mp.dps; mp.mp.dps=P+25
    s=mp.mpf(0)
    for a in range(1,D):
        c=chi(D,a)
        if c: s+=c*mp.zeta(k,mp.mpf(a)/D)
    v=s/mp.power(D,k)
    mp.mp.dps=old
    return v

def Lchi_series(D,k,P):
    old=mp.mp.dps; mp.mp.dps=P+25
    f=lambda n: chi(D,n)/mp.mpf(n)**k
    v=mp.nsum(f,[1,mp.inf],tol=mp.mpf(10)**(-P-10))
    mp.mp.dps=old
    return v

# ---------------- planner Li5 (correct) ----------------
_zcache = {}
def _zeta_int(n):
    if n not in _zcache or _zcache[n][0] < mp.mp.prec:
        _zcache[n] = (mp.mp.prec, mp.zeta(n) if n != 1 else None)
    return _zcache[n][1]
def li5_series(x, eps):
    tot = mp.mpf(0); p = x; n = 1
    while abs(p) > eps:
        tot += p / mp.mpf(n)**5; n += 1; p *= x
    return tot
def li5_log(x, eps):
    t = mp.log(x); tot = mp.mpf(0); tk = mp.mpf(1); fact = mp.mpf(1); k = 0
    H4 = mp.mpf(25)/12
    while True:
        if k == 4:
            term = (H4 - mp.log(-t)) * tk / fact
        else:
            z = _zeta_int(5-k)
            term = z * tk / fact
        tot += term
        if k > 10 and k % 2 == 0 and abs(term) < eps: break
        k += 1; tk *= t; fact *= k
    return tot
def li5(x, digits):
    eps = mp.mpf(10)**(-(digits + 10))
    if abs(x) <= mp.mpf('0.5'): return li5_series(x, eps)
    if x > 0: return li5_log(x, eps)
    return li5(x*x, digits)/16 - li5(-x, digits)

def li3(x, digits):
    eps = mp.mpf(10)**(-(digits + 10))
    if abs(x) <= mp.mpf('0.5'):
        tot = mp.mpf(0); p = x; n = 1
        while abs(p) > eps:
            tot += p/mp.mpf(n)**3; n += 1; p *= x
        return tot
    if x > 0:
        t = mp.log(x); tot = mp.mpf(0); tk = mp.mpf(1); fact = mp.mpf(1); k = 0
        H2 = mp.mpf(3)/2
        while True:
            if k == 2:
                term = (H2 - mp.log(-t)) * tk / fact
            else:
                z = _zeta_int(3-k)
                term = z * tk / fact
            tot += term
            if k > 10 and k % 2 == 0 and abs(term) < eps: break
            k += 1; tk *= t; fact *= k
        return tot
    return li3(x*x, digits)/4 - li3(-x, digits)

# ---------------- symbols (weight w) ----------------
def symbol_liw_dense(vx, v1m, d, w):
    """symbol of Li_w(x): -(1-x) otimes x^(w-1); dense list length d^w."""
    vec=[0]*(d**w)
    for i0 in range(d):
        c0=-v1m[i0]
        if c0==0: continue
        base0=i0*(d**(w-1))
        for t in itertools.product(range(d),repeat=w-1):
            cc=c0
            for j in t: cc*=vx[j]
            if cc: vec[base0+tidx(t,d)]+=cc
    return vec

def tidx(t,d):
    i=0
    for x in t: i=i*d+x
    return i

def tup_to_idx(tup,d):
    return tidx(tup,d)
def idx_to_tup(idx,d,w):
    t=[]
    for _ in range(w): t.append(idx%d); idx//=d
    return tuple(reversed(t))

def dense_to_dict(vec):
    return {i:c for i,c in enumerate(vec) if c}

def prod_mult(ms):
    cnt={}
    for i in ms: cnt[i]=cnt.get(i,0)+1
    m=1
    for k in cnt.values(): m*=math.factorial(k)
    return m

def prod_symbol_dense(ms,d,w):
    m=prod_mult(ms)
    v=[0]*(d**w)
    for p in set(itertools.permutations(ms)):
        v[tidx(p,d)]+=m
    return v

def sym_orbits(d, w):
    """map sorted tuple (multiset) -> list of distinct perm tuples."""
    orbits={}
    for t in itertools.product(range(d),repeat=w):
        s=tuple(sorted(t))
        orbits.setdefault(s,[]).append(t)
    return orbits

# ---------------- exact integer linear algebra ----------------
def _primitive(row):
    g=0
    for x in row: g=math.gcd(g,abs(int(x)))
    if g>1: row=[x//g for x in row]
    return row

def int_left_kernel_hnf(A):
    import flint as _fl
    M = _fl.fmpz_mat([list(map(int, r)) for r in A])
    H, U = M.hnf(transform=True)
    nr = M.nrows(); nc = M.ncols()
    rank = 0
    for i in range(nr):
        if any(H[i, j] for j in range(nc)): rank += 1
    K = []
    for i in range(rank, nr):
        K.append([int(U[i, j]) for j in range(U.ncols())])
    return K

def snf_invariants(B):
    if not B: return []
    M=flint.fmpz_mat([list(map(int,r)) for r in B])
    S=M.snf()
    n=min(S.nrows(),S.ncols())
    return [int(S[i,i]) for i in range(n)]

def saturate_rows(B):
    if not B: return B
    B=[_primitive(list(map(int,r))) for r in B]
    k=len(B); n=len(B[0])
    BT=[list(col) for col in zip(*B)]
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

def symbol_kernel(sym_dicts, d, w, N):
    """sym_dicts: list of dicts idx->coef (symbol of Li_w(x_i)).
    Returns (K, orbits): K saturated Z-basis of {c: sum c_i s_i in Sym^w}."""
    orbits = sym_orbits(d, w)
    support = set()
    for s in sym_dicts: support |= set(s.keys())
    # functionals: per orbit touching support
    cols = []  # each col = (w, ref) meaning functional v[w]-v[ref]; or (w,None) meaning v[w]
    for s, perms in orbits.items():
        sup = [t for t in perms if tidx(t,d) in support]
        if not sup: continue
        ref = s  # sorted representative
        if tidx(ref,d) in support:
            for t in sup:
                if t != ref: cols.append((tidx(t,d), tidx(ref,d)))
            if len(sup) < len(perms):
                cols.append((tidx(ref,d), None))
        else:
            for t in sup: cols.append((tidx(t,d), None))
    # build M rows (N x len(cols))
    M = []
    for s in sym_dicts:
        row = []
        for (wi, ri) in cols:
            row.append(s.get(wi,0) - (s.get(ri,0) if ri is not None else 0))
        M.append(row)
    K = int_left_kernel_hnf(M)
    # verify
    for c in K:
        for j in range(len(cols)):
            assert sum(c[i]*M[i][j] for i in range(N)) == 0
    K = saturate_rows(K)
    # certificate: each kernel combo's symbol is symmetric: check constancy on orbits
    for c in K:
        tot = {}
        for i,ci in enumerate(c):
            if ci:
                for idx,v in sym_dicts[i].items():
                    tot[idx] = tot.get(idx,0) + ci*v
        for s, perms in orbits.items():
            vals = set(tot.get(tidx(t,d),0) for t in perms)
            assert len(vals)==1, f"kernel combo not symmetric on orbit {s}"
    return K, orbits

def companion(sym_tot, d, w, orbits):
    """Given symmetric symbol dict tot (idx->coef), return dict ms->Fraction coefficient
    on the log^w product of ms with symbol(P)=tot.  prod_symbol(ms) = mult*symrow."""
    out={}
    for s, perms in orbits.items():
        cvals = set(sym_tot.get(tidx(t,d),0) for t in perms)
        assert len(cvals)==1
        c = cvals.pop()
        if c:
            out[s] = Fraction(c, prod_mult(list(s)))
    return out

def combo_symbol(K_row, sym_dicts, N):
    tot={}
    for i,ci in enumerate(K_row):
        if ci:
            for idx,v in sym_dicts[i].items():
                tot[idx]=tot.get(idx,0)+ci*v
    return tot

# ---------------- cobracket (3,2) ----------------
def _perm_sign(seq):
    inv = 0
    for i in range(len(seq)):
        for j in range(i+1, len(seq)):
            if seq[i] > seq[j]: inv += 1
    return -1 if inv % 2 else 1

def cobracket32(T, d, mode="full"):
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

# ---------------- LLL ----------------
def lll_rows(vals_list, P):
    mp.mp.dps = P + 30
    n = len(vals_list)
    Cc = mp.power(10, P)
    B = []
    for i, v in enumerate(vals_list):
        row = [0]*(n+1); row[i] = 1; row[n] = int(mp.nint(v*Cc))
        B.append(row)
    M = flint.fmpz_mat(B)
    t0 = time.time()
    R = M.lll(gram='exact')
    rows = [[int(R[i, j]) for j in range(n+1)] for i in range(n)]
    return rows, time.time()-t0


# ---------------- fast enumeration (integer A,B,C representation) ----------------
def _qd_to_abc(x):
    """QD x = a + b sqrt(D) = (A + B sqrt(D))/C, integers, C>0."""
    den = math.lcm(x.a.denominator, x.b.denominator)
    A = int(x.a*den); B = int(x.b*den)
    g = math.gcd(math.gcd(abs(A),abs(B)),den)
    return (A//g, B//g, den//g)

def _abc_mul(u, v, D):
    A1,B1,C1 = u; A2,B2,C2 = v
    A = A1*A2 + B1*B2*D
    B = A1*B2 + B1*A2
    C = C1*C2
    g = math.gcd(math.gcd(abs(A),abs(B)),C)
    if g>1: A//=g; B//=g; C//=g
    return (A,B,C)

def _abc_pow_table(g, B, D):
    """g^e for e in [-B,B] as abc tuples."""
    ga,gb,gc = _qd_to_abc(g)
    gi = _qd_to_abc(g.inv())
    tab = {0:(1,0,1)}
    u=(1,0,1)
    for e in range(1,B+1):
        u = _abc_mul(u,(ga,gb,gc),D); tab[e]=u
    u=(1,0,1)
    for e in range(1,B+1):
        u = _abc_mul(u,gi,D); tab[-e]=u
    return tab

def _abc_allowed_norm1m(u, D, primes):
    """Is N(1-x) a rational S-unit (only primes allowed)? x=(A+B r)/C."""
    A,B,C = u
    num = (C-A)*(C-A) - B*B*D
    den = C*C
    if num==0: return False
    num=abs(num)
    for p in primes:
        while num%p==0 and num>0: num//=p
        while den%p==0: den//=p
    return num==1 and den==1

def _abc_to_qd(u, D):
    A,B,C = u
    return QD(Fraction(A,C), Fraction(B,C), D)

def enum_exceptional_fast(F, B=8, time_limit=1500, progress_every=200000):
    import numpy as np
    D=F['D']; gens=F['letters']; d=len(gens)
    primes=F['rat_primes']
    Lf = np.array([float(mp.log(g.emb(60))) if g.emb(60)>0 else float(mp.log(-g.emb(60))) for g in gens])
    tabs = [_abc_pow_table(g,B,D) for g in gens]
    t0=time.time()
    axes=[np.arange(-B,B+1)]*d
    grids=np.meshgrid(*axes,indexing='ij')
    T=np.stack([g.ravel() for g in grids],axis=1).astype(float)
    logabs = T @ Lf
    keep = T[logabs < 1e-9].astype(int)
    keep = keep[np.any(keep!=0,axis=1)]
    out={}
    for ni, exps in enumerate(keep):
        if progress_every and ni%progress_every==0 and ni>0:
            el=time.time()-t0
            print(f"  {ni}/{len(keep)} ({el:.0f}s), found {len(out)}", flush=True)
            if el>time_limit:
                print(f"  enum TIME LIMIT; stopping early"); break
        u=(1,0,1)
        for gi,e in enumerate(exps):
            if e: u=_abc_mul(u,tabs[gi][int(e)],D)
        # |x|<1 exact: x=(A+B r)/C, need |A+B sqrt(D)| < C.  use sqrt via isqrt on scaled ints
        A,Bb,C = u
        # sign of A+B sqrt(D): compare A+B sqrt D with 0 via exact integer test
        # compute floor of B sqrt(D)
        import math as _m
        # |A + B sqrt D| < C  <=>  (A + B sqrt D)^2 < C^2 when ... safer: evaluate with mpmath at 40 digits (borderline impossible for these heights)
        v = A + Bb*float(D**0.5)
        if not (abs(v) < C): continue
        if not _abc_allowed_norm1m(u, D, primes): continue
        x = _abc_to_qd(u, D)
        for sgn in (1,-1):
            xx = x if sgn==1 else -x
            f = factor_sunit(1-xx, F)
            if f is None: continue
            out[xx]=dict(x=xx, vx=[int(e) for e in exps], v1m=f, h=int(max(abs(e) for e in exps)))
    res=sorted(out.values(), key=lambda r:(r['h'], r['vx'], -r['x'].a))
    print(f"  enum done: {len(keep)} tuples checked, {len(res)} exceptional, {time.time()-t0:.0f}s")
    return res


def saturate_rows_fast(B):
    """Saturated basis of span_Q(B) cap Z^n via HNF + flint inverse."""
    if not B: return B
    B=[_primitive(list(map(int,r))) for r in B]
    k=len(B); n=len(B[0])
    BT=[list(col) for col in zip(*B)]          # n x k
    H=flint.fmpz_mat(BT).hnf()
    H0=flint.fmpz_mat([[int(H[i,j]) for j in range(k)] for i in range(k)])
    for i in range(k,n):
        for j in range(k):
            assert int(H[i,j])==0, "hnf bottom not zero"
    H0i, den = H0.inv()                        # H0^{-1} = H0i/den
    Bm = flint.fmpz_mat(B)
    num = H0i.transpose() * Bm                 # k x n
    assert num.nrows()==k and num.ncols()==n
    rows=[]
    for i in range(k):
        row=[]
        for j in range(n):
            x = int(num[i,j])
            assert x % den == 0, "saturation not integral?!"
            row.append(x//den)
        rows.append(_primitive(row))
    return rows

def symbol_kernel_fast(M, N, ncols_quotient):
    """M: N x q integer matrix (quotient functionals applied to symbol rows).
    Returns saturated integer kernel basis of {c: c M = 0}."""
    MT=[list(col) for col in zip(*M)]
    NS, nulldim = flint.fmpz_mat(MT).nullspace()
    K=[]
    for j in range(nulldim):
        K.append([int(NS[i,j]) for i in range(NS.nrows())])
    K = [_primitive(r) for r in K]
    K = saturate_rows_fast(K)
    # verify
    for c in K:
        for j in range(ncols_quotient):
            assert sum(c[i]*M[i][j] for i in range(N))==0
    return K
