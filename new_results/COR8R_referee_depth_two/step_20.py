engine = r'''"""
COR8R referee symbol engine -- written from definitions (no COR8 code).

Conventions (stated explicitly):
- Li_{a,b}(x,y): TWO conventions supported.
    convA ("math", Goncharov): Li_{a,b}(x,y) = sum_{0<k1<k2} x^{k1} y^{k2}/(k1^a k2^b).
      G-form: Li_{a,b}(x,y) = G(0^{b-1}, 1/y, 0^{a-1}, 1/(x y); 1).
    convB ("COR4"):            Li_{a,b}(x,y) = sum_{n>m>=1} x^n y^m/(n^a m^b) = convA-Li_{b,a}(y,x).
      G-form: Li_{a,b}(x,y) = G(0^{a-1}, 1/x, 0^{b-1}, 1/(x y); 1).
- G(a1,...,an; z) = int_0^z dt1/(t1-a1) ... ; outermost integration = first letter a1.
- Differential (derived by parts, verified numerically for n=2,3):
    dG(a1..an;z) = sum_i G(skip a_i; z) dlog((a_{i-1}-a_i)/(a_i-a_{i+1})),  a0 := z, a_{n+1} := 0.
- Symbol recursion: S(G) = sum_i S(G(skip a_i)) tensor [letter_i]; terms with constant letter (dlog=0) dropped.
- Tensor order: the letter of the OUTERMOST integration is LAST (paper's convention).
  Hence S(Li_n(z)) = -[1-z][z]^{n-1}; S(H_{a1..an}(c)) = g(a_n)...g(a_1) with
  g(0)=[c], g(1)=-[1-c], g(-1)=[1+c].
- Letters: exponent vectors over basis B = (-1, 2, c, 1-c, 1+c); coefficients carried in
  Laurent polynomials in c (dict exp -> Fraction), i.e. rational-function coefficients c^E
  are treated as scalars (paper's convention).
- delta_{2,2}: for a weight-4 symbol sum x1^x2^x3^x4 -> sum (x1 wedge x2) wedge (x3 wedge x4)
  in Lambda^2(Lambda^2 B1)  (= the Lambda^2 B2 component of the (2,2) part of the coproduct,
  via symbol-of-coproduct = deconcatenation; cf. Duhr 1207.4959).
"""
import sympy as sp
from fractions import Fraction
from functools import lru_cache
from itertools import combinations

cc = sp.Symbol('c')
NB = 5  # basis size: (-1, 2, c, 1-c, 1+c)
ZL = (0,)*NB

# ---------------- letters ----------------
def Ladd(a,b): return tuple(x+y for x,y in zip(a,b))
def Lneg(a):   return tuple(-x for x in a)

def factor_Qc(f):
    """Factor a rational function of c over the basis (-1,2,c,1-c,1+c).
    Returns exponent 5-tuple. Raises if any other prime/irreducible appears."""
    f = sp.cancel(sp.sympify(f))
    if f == 0:
        raise ValueError("zero letter")
    num, den = sp.fraction(f)
    vec = [0]*NB
    for poly, sgn in ((num,1),(den,-1)):
        poly = sp.Poly(sp.expand(poly), cc)
        coeff, factors = sp.factor_list(poly)
        coeff = Fraction(coeff)
        # constant content: must be +- power of 2
        cnum, cden = coeff.numerator, coeff.denominator
        for v, s2 in ((cnum,1),(cden,-1)):
            v = abs(v)
            if v & (v-1) == 0:
                e = v.bit_length()-1
                vec[1] += sgn*s2*e
            else:
                raise ValueError(f"constant factor {v} not a power of 2 in {f}")
        if cnum < 0: vec[0] += sgn
        for fac, e in factors:
            fac = sp.expand(fac)
            if fac == cc:
                vec[2] += sgn*e
            elif fac == cc-1:      # [c-1] = [1-c]+[-1]
                vec[3] += sgn*e; vec[0] += sgn*e
            elif fac == 1-cc:
                vec[3] += sgn*e
            elif fac == cc+1:
                vec[4] += sgn*e
            else:
                raise ValueError(f"factor {fac} not in basis, from {f}")
    return tuple(vec)

def is_constant_letter(L):
    return L[2]==0 and L[3]==0 and L[4]==0   # in span of [-1],[2]: dlog = 0

# ---------------- Laurent polynomials in c ----------------
ONE = {0:Fraction(1)}
def LP_add(p,q):
    r = dict(p)
    for k,v in q.items():
        nv = r.get(k,Fraction(0))+v
        if nv: r[k]=nv
        elif k in r: del r[k]
    return r
def LP_neg(p): return {k:-v for k,v in p.items()}
def LP_scale(a,p): return {k:a*v for k,v in p.items()}
def LP_mul(p,q):
    r = {}
    for k1,v1 in p.items():
        for k2,v2 in q.items():
            r[k1+k2] = r.get(k1+k2,Fraction(0))+v1*v2
    return {k:v for k,v in r.items() if v}
def LP_pow_c(e): return {e:Fraction(1)}

# ---------------- symbols ----------------
# symbol: dict  word(tuple of letters) -> Laurent coeff
def S_add(s1,s2):
    r = dict(s1)
    for w,cf in s2.items():
        r[w] = LP_add(r.get(w,{}), cf)
        if not r[w]: del r[w]
    return r
def S_scaleLP(s,p):
    if not p: return {}
    return {w:LP_mul(cf,p) for w,cf in s.items()}
def S_append(s, L, coeffLP=ONE):
    """tensor on the right with letter L, times scalar coeff."""
    if not coeffLP: return {}
    return {w+(L,):LP_mul(cf,coeffLP) for w,cf in s.items()}
def S_shuffle(s1,s2):
    r = {}
    for w1,c1 in s1.items():
        for w2,c2 in s2.items():
            cc12 = LP_mul(c1,c2)
            if not cc12: continue
            for w in _shuffle_words(w1,w2):
                r[w] = LP_add(r.get(w,{}), cc12)
    return {w:v for w,v in r.items() if v}

def _shuffle_words(w1,w2):
    n1,n2 = len(w1),len(w2)
    out = []
    def rec(i,j,acc):
        if i==n1 and j==n2:
            out.append(tuple(acc)); return
        if i<n1:
            acc.append(w1[i]); rec(i+1,j,acc); acc.pop()
        if j<n2:
            acc.append(w2[j]); rec(i,j+1,acc); acc.pop()
    rec(0,0,[])
    return out

# ---------------- G-integrals: symbols ----------------
@lru_cache(maxsize=None)
def symG(word):
    """Symbol of regularized G(word; 1), word = tuple of sympy exprs in c (or ints)."""
    word = tuple(sp.sympify(w) for w in word)
    n = len(word)
    if n==0:
        return {(): ONE}
    if all(w==0 for w in word):
        return {}
    # base case: 0^k, A
    if all(w==0 for w in word[:-1]) and word[-1]!=0:
        return _symG_zerosA(word)
    if _good(word):
        return _symG_recurse(word)
    return _symG_bad(word)

def _good(word):
    if word[-1]==0: return False
    for i in range(len(word)-1):
        if sp.simplify(word[i]-word[i+1])==0: return False
    return True

def _symG_zerosA(word):
    # G(0^k, A; 1) = -Li_{k+1}(1/A);  d/dA G(0^k,A;1) = -G(0^{k-1},A;1)/A,  S(G(A;1)) = [(1-A)/A]
    k = len(word)-1
    A = word[-1]
    L1 = Ladd(factor_Qc(1-A), Lneg(factor_Qc(A)))   # [(1-A)/A]
    LA = factor_Qc(A)
    s = { (L1,): ONE }
    for _ in range(k):
        s = S_append(s, LA, LP_neg(ONE))   # sign (-1)^k accumulated
    return s

def _symG_recurse(word):
    n = len(word)
    a = (sp.Integer(1),)+tuple(word)+(sp.Integer(0),)   # a0=z=1, a_{n+1}=0
    res = {}
    for i in range(1,n+1):
        lv = sp.cancel((a[i-1]-a[i])/(a[i]-a[i+1]))
        L = factor_Qc(lv)
        if is_constant_letter(L):
            continue
        sub = tuple(word[:i-1])+tuple(word[i:])
        res = S_add(res, S_append(symG(sub), L))
    return res

# ---- bad words: shuffle-regularization linear solver ----
def _values_of(word):
    vs = []
    for w in word:
        if w!=0 and w not in vs: vs.append(w)
    return vs

def _bad_words(length, values):
    import itertools
    out = []
    for w in itertools.product([sp.Integer(0)]+list(values), repeat=length):
        if all(x==0 for x in w): continue
        if not _good(w) and not (all(x==0 for x in w[:-1]) and w[-1]!=0):
            out.append(w)
    return out

@lru_cache(maxsize=None)
def _symG_bad(word):
    word = tuple(sp.sympify(w) for w in word)
    L = len(word)
    values = _values_of(word)
    unk = _bad_words(L, values)
    idx = {w:i for i,w in enumerate(unk)}
    # equations: shuffle(u,v) = sym(u)*sym(v) for all u,v with len(u)+len(v)=L, len>=1
    rows = []   # (int-coeff vector over unk, known symbol)
    import itertools
    alphabet = [sp.Integer(0)]+list(values)
    for lu in range(1,L):
        lv = L-lu
        for u in itertools.product(alphabet, repeat=lu):
            su = symG(u)   # known (shorter, or base)
            for v in itertools.product(alphabet, repeat=lv):
                sv = symG(v)
                prod = S_shuffle(su,sv)
                row = [Fraction(0)]*len(unk)
                rhs = dict(prod)
                ok = False
                for w in _shuffle_words(u,v):
                    if w in idx:
                        row[idx[w]] += 1; ok = True
                    else:
                        sw = symG(w)
                        rhs = S_add(rhs, {ww:LP_neg(cf) for ww,cf in sw.items()})
                if ok:
                    rows.append((row,rhs))
    # solve integer system with symbol RHS by Gaussian elimination over Q
    sol = _solve_symbol_system(rows, len(unk))
    if word not in idx:
        raise RuntimeError("bad word not in unknowns")
    if sol[word] is None:
        raise RuntimeError(f"shuffle system underdetermined for {word}")
    return sol[word]

def _solve_symbol_system(rows, nunk):
    # Gaussian elimination: rows are (list[Fraction], symbol); returns dict index->symbol or None
    M = [ (list(r), s) for r,s in rows ]
    sol_sym = [None]*nunk
    pivcol_of_row = []
    R = 0
    for col in range(nunk):
        piv = None
        for r in range(R, len(M)):
            if M[r][0][col] != 0:
                piv = r; break
        if piv is None: continue
        M[R], M[piv] = M[piv], M[R]
        inv = Fraction(1)/M[R][0][col]
        M[R] = ([x*inv for x in M[R][0]], S_scaleLP_frac(M[R][1], inv))
        for r in range(len(M)):
            if r!=R and M[r][0][col]!=0:
                f = M[r][0][col]
                M[r] = ([a - f*b for a,b in zip(M[r][0], M[R][0])],
                        S_add(M[r][1], S_scaleLP_frac(M[R][1], -f)))
        pivcol_of_row.append((R,col))
        R += 1
        if R==len(M): break
    # extract
    for r,col in pivcol_of_row:
        if all(M[r][0][c2]==0 for c2 in range(nunk) if c2!=col):
            sol_sym[col] = M[r][1]
    # consistency: rows reduced to zero coeffs must have zero rhs
    for row,s in M:
        if all(x==0 for x in row):
            if s:
                raise RuntimeError("inconsistent shuffle system")
    return {i:sol_sym[i] for i in range(nunk)}

def S_scaleLP_frac(s, f):
    fr = Fraction(f)
    return {w:{k:v*fr for k,v in cf.items()} for w,cf in s.items()}

# ---------------- MPL wrappers ----------------
def Li22_sym(xval, yval, conv):
    """Symbol of Li_{2,2}(x,y) at x=xval, y=yval (sympy exprs in c)."""
    x = sp.sympify(xval); y = sp.sympify(yval)
    if conv=='A':
        word = (sp.Integer(0), 1/y, sp.Integer(0), 1/(x*y))
    else:
        word = (sp.Integer(0), 1/x, sp.Integer(0), 1/(x*y))
    return symG(word)

def Li31_sym(xval, yval, conv):
    x = sp.sympify(xval); y = sp.sympify(yval)
    if conv=='A':
        word = (1/y, sp.Integer(0), sp.Integer(0), 1/(x*y))
    else:
        word = (sp.Integer(0), sp.Integer(0), 1/x, 1/(x*y))
    return symG(word)

# ---------------- HPLs ----------------
def hpl_sym(w):
    """S(H_w(c)), w = tuple in {-1,0,1}; outermost letter last: g(a_n)...g(a_1)."""
    s = {(): ONE}
    for a in reversed(w):
        if a==0:
            s = S_append(s, (0,0,1,0,0))
        elif a==1:
            s = S_append(s, (0,0,0,1,0), LP_neg(ONE))
        elif a==-1:
            s = S_append(s, (0,0,0,0,1))
        else:
            raise ValueError(a)
    return s

# ---------------- delta_{2,2} ----------------
PAIRS = [(i,j) for i in range(NB) for j in range(i+1,NB)]  # 10 basis elements of Lambda^2 B1
PAIR_IDX = {p:k for k,p in enumerate(PAIRS)}
WP = [(i,j) for i in range(len(PAIRS)) for j in range(i+1,len(PAIRS))]  # basis of Lambda^2(Lambda^2)

def wedge2(a,b):
    """a wedge b in Lambda^2 B1: dict PAIR_IDX -> Fraction"""
    r = {}
    for i in range(NB):
        for j in range(i+1,NB):
            v = a[i]*b[j]-a[j]*b[i]
            if v: r[PAIR_IDX[(i,j)]] = r.get(PAIR_IDX[(i,j)],Fraction(0))+v
    return {k:v for k,v in r.items() if v}

def delta22(sym):
    """sum coeff*(l1^l2)^(l3^l4); returns dict (I,J) I<J in PAIR basis -> Laurent poly."""
    out = {}
    for w,cf in sym.items():
        if len(w)!=4: 
            if cf: raise ValueError("weight-4 symbol expected")
            continue
        P = wedge2(w[0],w[1])
        Q = wedge2(w[2],w[3])
        for I,vI in P.items():
            for J,vJ in Q.items():
                if I==J: continue
                if I<J:
                    key=(I,J); val = vI*vJ
                else:
                    key=(J,I); val = -vI*vJ
                add = LP_scale(val, cf)
                out[key] = LP_add(out.get(key,{}), add)
    return {k:v for k,v in out.items() if v}

def drop_letter(sym, bidx):
    """Set basis letter bidx to zero: drop all words containing it (killed)."""
    return {w:cf for w,cf in sym.items() if all(L[bidx]==0 for L in w)}
'''
open('cor8r_sym.py','w').write(engine)
print("engine written", len(engine))
