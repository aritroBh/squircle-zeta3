# step_88: write verify_pgen1.py
verify = r'''"""verify_pgen1.py -- independent verification of the PGEN1 results.
Covers: controls, cross-ratio alphabet checks (generic octagon per edge; regular polygons),
delta(N3) computations (generic octagon per edge + total; regular octagon; decagon; 12-gon),
the symbol-level alphabet proof for s <= 4, and the task-3 L3 numbers.
Runs in < 5 min. Labels: PROVED (exact), NUMERICAL (digits stated).
"""
import os, sys, pickle, itertools, time
from fractions import Fraction as Fr
t0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "<project>"
O9B = os.path.join(BASE, "research/08_theory/OCT9b_octagon_debug")
O9D = os.path.join(BASE, "research/08_theory/OCT9d_octagon_short")
import sympy as sp
from flint import fmpz_mat
import mpmath as mp

def vadd(*vs):
    r={}
    for v in vs:
        for k,c in v.items():
            r[k]=r.get(k,Fr(0))+c
            if r[k]==0: del r[k]
    return r
def wedge(u,v):
    out={}
    for a,ca in u.items():
        for b,cb in v.items():
            if a==b: continue
            if a<b: key,s=(a,b),1
            else: key,s=(b,a),-1
            out[key]=out.get(key,Fr(0))+s*ca*cb
    return {k:c for k,c in out.items() if c}
def vecQ(q):
    q=Fr(q); out={}
    for m,sgn in ((abs(q.numerator),1),(q.denominator,-1)):
        for p,ex in sp.factorint(int(m)).items():
            out[p]=out.get(p,Fr(0))+sgn*ex
    return {k:c for k,c in out.items() if c}
def delta_li3(terms):
    out={}
    for c,vz,vomz in terms:
        w=wedge(vz,vomz)
        for (a,b),cab in w.items():
            for l3,c3 in vz.items():
                key=((a,b),l3)
                out[key]=out.get(key,Fr(0))+c*cab*c3
    return {k:v for k,v in out.items() if v}
def wedge3_image(D):
    out={}
    for ((a,b),c),coef in D.items():
        letters = sorted(set([a,b,c]))
        if len(letters)<3: continue
        perm=[a,b,c]; inv=0
        for i in range(3):
            for j in range(i+1,3):
                if perm[i]>perm[j]: inv+=1
        sgn=-1 if inv%2 else 1
        key=tuple(letters)
        out[key]=out.get(key,Fr(0))+sgn*coef
    return {k:v for k,v in out.items() if v}

def cross_ratios_Q(points):
    out=set()
    for a,b,c,d in itertools.permutations(points,4):
        try:
            if a is None: val=(b-d)/(b-c)
            elif b is None: val=(a-c)/(a-d)
            elif c is None: val=(b-d)/(a-d)
            elif d is None: val=(a-c)/(b-c)
            else: val=(a-c)*(b-d)/((a-d)*(b-c))
        except ZeroDivisionError: continue
        out.add(Fr(val))
    return sorted(out)

# ---------------- Quotient machinery (symbols mod products), from OCT8 ----------------
class Quotient:
    def __init__(self, letters):
        self.letters = list(letters)
        self.dim3 = len(self.letters)**3
        self.basis3 = list(itertools.product(self.letters, repeat=3))
        self.bidx = {t:i for i,t in enumerate(self.basis3)}
        D = []
        for a,b,c in itertools.product(self.letters, repeat=3):
            v = sp.zeros(self.dim3,1)
            for t in [(a,b,c),(b,a,c),(b,c,a)]:
                v[self.bidx[t],0] += 1
            D.append(v)
        self.D = D
        self.DM = sp.Matrix.hstack(*D)
        self.rank_D = self.DM.rank()
    def coord(self, S):
        v = sp.zeros(self.dim3,1)
        for t,c in S.items():
            v[self.bidx[t],0] += sp.Rational(c.numerator, c.denominator)
        return v
def S_Li3_Q(z):
    zv, omzv = vecQ(z), vecQ(1-z)
    out = {}
    for a,ca in omzv.items():
        for b,cb in zv.items():
            for c,cc in zv.items():
                out[(a,b,c)] = out.get((a,b,c), Fr(0)) - ca*cb*cc
    return out

print("="*70)
print("PART 0: CONTROLS")
print("="*70)
d = delta_li3([(1, vecQ(Fr(1,2)), vecQ(Fr(1,2)))])
print("delta{1/2}_3 =", d if d else "0 (empty)", " [PROVED, consistent with Li_3(1/2)=7/8 zeta(3)+products]")
ok = all(not delta_li3([(1,vecQ(x),vecQ(1-x)),(1,vecQ(1-x),vecQ(x)),(1,vecQ(1-1/x),vecQ(1/x))]) for x in [Fr(1,3),Fr(2,5),Fr(3,7),Fr(-2,3),Fr(5,3),Fr(7,11),Fr(-1,2)])
print("delta(Landen triple) = 0 for 7 sample rationals:", ok, " [PROVED]")
ok = True
for r in [Fr(2), Fr(1,2), Fr(3), Fr(1,3), Fr(5,2), Fr(7,4)]:
    p = 1/(1-r); x = 1-r
    CR = cross_ratios_Q([Fr(0), Fr(1), None, p])
    if not ((x in CR) and (x/(x-1) in CR)): ok = False
print("thm:hexagon's two arguments in CR{0,1,inf,p} for 6 sample r:", ok, " [PROVED]")
for r in [Fr(2), Fr(1,2)]:
    x = 1-r; z2 = x/(x-1)
    dd = delta_li3([(Fr(3),vecQ(x),vecQ(1-x)),(Fr(-5,7),vecQ(z2),vecQ(1-z2))])
    print(f"  delta at r={r} (x={x}, x/(x-1)={z2}) =", dd if dd else "0 (empty)", " [PROVED: canonical non-local zeta(3) there]")

print()
print("="*70)
print("PART 1: GENERIC RATIONAL OCTAGON -- alphabet rule and delta per edge")
print("="*70)
K12 = {0:(Fr(3,2),Fr(-7,12)), 1:(Fr(1,3),Fr(1,9)), 2:(Fr(5,3),Fr(-5,9)), 3:(Fr(1,4),Fr(1,24))}
TW3 = pickle.load(open(os.path.join(O9D,"TW3.pkl"),"rb"))
sf = pickle.load(open(os.path.join(O9D,"short_form.pkl"),"rb"))
G8 = [Fr(s) for s in sf['G8']]
B0 = [Fr(s) for s in sf['B0']]
CRp = pickle.load(open(os.path.join(HERE,"CR_presentations.pkl"),"rb"))
deltas = {}
for e, primes in [(1,(2,3,5)), (2,(2,3,5)), (3,(2,3,5)), (0,(2,3,5,7))]:
    k1,k2 = K12[e]; p1,p2 = -1/k1,-1/k2
    CR = cross_ratios_Q([Fr(0),Fr(1),None,p1,p2])
    Q = Quotient(primes)
    gens = [Q.coord(S_Li3_Q(z)) for z in CR]
    M = sp.Matrix.hstack(*[g.T for g in gens])
    MD = sp.Matrix.hstack(M, Q.DM)
    r_all = MD.rank(); rD = Q.rank_D
    t = sp.Matrix(TW3[e]); t = t.T if t.shape[0]==1 else t
    inspan = sp.Matrix.hstack(MD, t).rank() == r_all
    n = len(primes)
    print(f"edge {e}: kappa=({k1},{k2}), p=({p1},{p2}), |CR_e|={len(CR)}; letters {primes} (n={n}): "
          f"dim Q_3=(n^3-n)/3={(n**3-n)//3}; rank span s(CR) mod D = {r_all-rD}; TW3_{e} in span: {inspan} [PROVED]")
    # delta from OCT9d short basis
    if e == 0:
        coeffs = [Fr(s) for s in sf['sol0'][:20]]
        terms = [(c, vecQ(z), vecQ(1-z)) for c,z in zip(coeffs, B0) if c]
    else:
        coeffs = [Fr(s) for s in sf['sols123'][e][:8]]
        terms = [(c, vecQ(z), vecQ(1-z)) for c,z in zip(coeffs, G8) if c]
    d1 = delta_li3(terms)
    # independence check vs cross-ratio presentation (edges 1-3 recomputed; edge 0 stored)
    pres = CRp[e]['presentation']
    d2 = delta_li3([(c, vecQ(z), vecQ(1-z)) for c,z in pres])
    same = vadd(d1, {k:-v for k,v in d2.items()}) == {}
    deltas[e] = d1
    print(f"   delta_{e}: {len(d1)} terms; from OCT9d basis == from CR presentation: {same} [PROVED]")
    print(f"   delta_{e} != 0: {bool(d1)}")
print("edge-0 delta (all terms):")
for k,v in sorted(deltas[0].items()):
    print(f"   {v} * ({k[0][0]}^{k[0][1]}) (x) {k[1]}")
print("wedge^3(delta_edge1) =", wedge3_image(deltas[1]), " [empty => in Lie^3, PROVED]")

print()
print("="*70)
print("PART 2: REGULAR POLYGONS -- delta")
print("="*70)
def qmul(x, y, D): return (x[0]*y[0]+D*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def qpow(x, n, D):
    if n < 0:
        a,b = x; den = a*a-D*b*b
        return qpow((a/den, -b/den), -n, D)
    r = (Fr(1), Fr(0))
    while n:
        if n & 1: r = qmul(r, x, D)
        x = qmul(x, x, D); n >>= 1
    return r
def mpf_of(q): q=Fr(q); return mp.mpf(q.numerator)/mp.mpf(q.denominator)
def vp(q, p):
    q = abs(Fr(q)); n, d = q.numerator, q.denominator
    e = 0
    while n % p == 0 and n>1: n//=p; e+=1
    while d % p == 0 and d>1: d//=p; e-=1
    return e
class QF:
    def __init__(self, D, letters):
        self.D = Fr(D); self.LAB = {n: ab for n, ab, _ in letters}
        self.LNM = {n: nrm for n, _, nrm in letters}; self.letters = letters
    def elem(self, v):
        r = (Fr(1), Fr(0))
        for nm, e in v.items(): r = qmul(r, qpow(self.LAB[nm], int(e), self.D), self.D)
        return r
    def vec(self, z):
        a,b = Fr(z[0]), Fr(z[1]); N = a*a - self.D*b*b
        exps = {}
        for nm, ab, nrm in self.letters:
            if abs(nrm) == 1: continue
            fct = sp.factorint(abs(nrm)); pr, k = next(iter(fct.items()))
            vv = vp(N, pr); assert vv % k == 0, (z,nm,vv)
            exps[nm] = vv//k
        rem = abs(N)
        for nm, e in exps.items(): rem /= abs(self.LNM[nm])**e
        assert rem == 1, (z, rem)
        mp.mp.dps = 40
        sq = mp.sqrt(int(self.D))
        def aval(ab): return abs(mpf_of(ab[0]) + sq*mpf_of(ab[1]))
        lt = mp.log(aval((a,b))); prod_log = mp.mpf(0)
        for nm, ab, nrm in self.letters:
            e = exps.get(nm, 0)
            if e and abs(nrm) != 1: prod_log += e*mp.log(aval(ab))
        units = [(nm,ab) for nm,ab,nrm in self.letters if abs(nrm)==1]
        resid = lt - prod_log
        if units:
            nm, ab = units[0]
            c = resid/mp.log(aval(ab)); ci = int(mp.nint(c))
            assert abs(c-ci) < mp.mpf('1e-25'), (z,c)
            exps[nm] = exps.get(nm,0)+ci
        got = self.elem(exps)
        assert got == (a,b) or got == (-a,-b), (z, exps, got)
        return {k:Fr(v) for k,v in exps.items() if v}
F2 = QF(2, [("r2",(Fr(0),Fr(1)),-2), ("u",(Fr(1),Fr(1)),-1)])
F3 = QF(3, [("g2",(Fr(1),Fr(1)),-2), ("g3",(Fr(0),Fr(1)),-3), ("u",(Fr(2),Fr(1)),1)])
F5 = QF(5, [("2",(Fr(2),Fr(0)),4), ("s5",(Fr(0),Fr(1)),-5), ("phi",(Fr(1,2),Fr(1,2)),-1)])

def delta_li3_split(terms):
    """terms: (p, q, vz, vomz) with coeff p+q sqrt(D); returns (P, Q)."""
    P,Q={},{}
    for p,q,vz,vomz in terms:
        w=wedge(vz,vomz)
        for (a,b),cab in w.items():
            for l3,c3 in vz.items():
                key=((a,b),l3)
                if p: P[key]=P.get(key,Fr(0))+p*cab*c3
                if q: Q[key]=Q.get(key,Fr(0))+q*cab*c3
    return ({k:v for k,v in P.items() if v},{k:v for k,v in Q.items() if v})

# regular octagon: Li3 coeffs (in S, P'''=48 u S), from OCT6 proved table
s2 = sp.sqrt(2)
oct6 = [((2,-1), Fr(3,8), Fr(-1,4)),       # d = 2-sqrt2
        ((1,-1), Fr(-3,8), Fr(1,4)),       # k = sqrt2-1  [wait: k = sqrt2-1 = (-1,1)! fix below]
        ]
# build explicitly: (a,b) means a + b sqrt2
pts2 = {"d": (2,-1), "k": (-1,1), "mr2h": (Fr(-1,2),0), "r2h": (Fr(1,2),0), "x1": (1,-1), "x2": (1,Fr(-1,2))}
coefs2 = {"d": (Fr(3,8),Fr(-1,4)), "k": (Fr(-3,8),Fr(1,4)), "mr2h": (Fr(-9,2),Fr(3)),
          "r2h": (Fr(-155,32),Fr(103,32)), "x1": (Fr(33,32),Fr(-21,32)), "x2": (Fr(-79,32),Fr(59,32))}
terms = []
for nm,(p,q) in coefs2.items():
    z = pts2[nm]
    z = tuple(Fr(t) for t in z)
    omz = tuple(Fr(1)-z[0], Fr(0)-z[1]) if False else (Fr(1)-z[0], Fr(0)-z[1])
    terms.append((p,q,F2.vec(z),F2.vec(omz)))
P,Q = delta_li3_split(terms)
print("regular octagon delta:")
print("  P part:", P)
print("  Q part:", Q)
print("  delta != 0:", bool(P or Q), " [PROVED]  wedge^3:", wedge3_image(P), wedge3_image(Q))

# decagon
pts5 = {"1/phi": (Fr(-1,2),Fr(1,2)), "-1/phi": (Fr(1,2),Fr(-1,2))}
coefs5 = {"1/phi": (Fr(-75,32),Fr(81,64)), "-1/phi": (Fr(-35,32),Fr(53,64))}
terms = []
for nm,(p,q) in coefs5.items():
    z = pts5[nm]
    terms.append((p,q,F5.vec(z),F5.vec((Fr(1)-z[0], Fr(0)-z[1]))))
P,Q = delta_li3_split(terms)
print("decagon delta:", P, Q, " == 0:", not P and not Q, " [PROVED]")

# 12-gon
z8 = {"z1":(2,-1),"z2":(Fr(-1,2),Fr(1,2)),"z3":(0,Fr(1,3)),"z4":(7,-4),"z5":(1,Fr(-1,2)),
      "z6":(-1,Fr(2,3)),"z7":(Fr(1,2),Fr(-1,6)),"z8":(Fr(1,3),0)}
coefs3 = [("z1",Fr(-73,24),Fr(37,8)),("z2",Fr(11,2),Fr(-39,4)),("z3",Fr(97,24),Fr(-63,8)),
          ("z4",Fr(19,192),Fr(-13,64)),("z5",Fr(-1,8),Fr(9,8)),("z6",Fr(2,3),Fr(-3,4)),
          ("z7",Fr(7,4),Fr(-13,4)),("z8",Fr(-17,96),Fr(13,16))]
terms = []
for nm,p,q in coefs3:
    z = z8[nm]
    terms.append((p,q,F3.vec(z),F3.vec((Fr(1)-z[0], Fr(0)-z[1]))))
P,Q = delta_li3_split(terms)
print("12-gon delta (P part):", P)
print("12-gon delta (Q part):", Q)
print("  delta != 0:", bool(P or Q), " [PROVED]  wedge^3:", wedge3_image(P), wedge3_image(Q))
# duplication invariance: z4 = z1^2 -> 4 z1 + 4 (-z1)
terms2 = []
for nm,p,q in coefs3:
    z = z8[nm]
    if nm == "z4":
        terms2.append((4*p,4*q,F3.vec(z8["z1"]),F3.vec((Fr(1)-z8["z1"][0], Fr(0)-z8["z1"][1]))))
        mz1 = z8["z1"]
        terms2.append((4*p,4*q,F3.vec(mz1),F3.vec((Fr(1)+mz1[0], mz1[1]))))
    else:
        terms2.append((p,q,F3.vec(z),F3.vec((Fr(1)-z[0], Fr(0)-z[1]))))
P2,Q2 = delta_li3_split(terms2)
print("  12-gon delta invariant under duplication z4 = z1^2:", vadd(P,{k:-v for k,v in P2.items()})=={} and vadd(Q,{k:-v for k,v in Q2.items()})=={}, " [PROVED]")

# 12-gon alphabet check
def qsub(x,y): return (x[0]-y[0], x[1]-y[1])
def qinv(x,D):
    a,b=x; den=a*a-D*b*b; return (a/den,-b/den)
def qdiv(x,y,D): return qmul(x,qinv(y,D),D)
def cr4(a,b,c,d,D):
    def sub(x,y):
        if x is None: return None
        if y is None: return None
        return qsub(x,y)
    n1,n2,d1,d2 = sub(a,c),sub(b,d),sub(a,d),sub(b,c)
    if n1 is None and d1 is None: return qdiv(n2,d2,D)
    if n2 is None and d2 is None: return qdiv(n1,d1,D)
    if n1 is None and d2 is None: return qdiv(n2,d1,D)
    if n2 is None and d1 is None: return qdiv(n1,d2,D)
    num = qmul(n1,n2,D); den = qmul(d1,d2,D)
    return qdiv(num,den,D)
def S3orbit(z,D):
    ONE=(Fr(1),Fr(0)); out={z,qsub(ONE,z)}
    if z!=(Fr(0),Fr(0)): out.add(qinv(z,D))
    if z!=ONE: out.add(qinv(qsub(ONE,z),D))
    if z not in ((Fr(0),Fr(0)),ONE):
        out.add(qdiv(z,qsub(z,ONE),D)); out.add(qdiv(qsub(z,ONE),z,D))
    return out
def cross_ratios_QF(points, D):
    out = set()
    for quad in itertools.combinations(points, 4):
        val = cr4(*quad, D)
        out |= S3orbit(val, D)
    return out
p12 = [(Fr(3,2),Fr(1,2)), (Fr(4),Fr(2)), (Fr(-3),Fr(-2)), (Fr(-1,2),Fr(-1,2))]
CR12 = cross_ratios_QF([(Fr(0),Fr(0)),(Fr(1),Fr(0)),None]+p12, 3)
z8pts = [tuple(Fr(t) for t in z8[nm]) for nm in ["z1","z2","z3","z4","z5","z6","z7","z8"]]
inCR = [z in CR12 for z in z8pts]
z4sq = qmul(z8pts[0],z8pts[0],3)
print("12-gon: |CR_e| =", len(CR12), "; 8 DOD5 points in CR_e:", inCR, "; z4 = z1^2 duplication:", z4sq == z8pts[3], " [PROVED]")
p10 = [(Fr(-1,2),Fr(-1,2)), (Fr(3,2),Fr(1,2))]
CR10 = cross_ratios_QF([(Fr(0),Fr(0)),(Fr(1),Fr(0)),None]+p10, 5)
print("decagon: CR_e =", sorted(CR10), "; +-1/phi in CR:", (Fr(-1,2),Fr(1,2)) in CR10 and (Fr(1,2),Fr(-1,2)) in CR10, " [PROVED]")

print()
print("="*70)
print("PART 3: SYMBOL-LEVEL ALPHABET PROOF for s <= 4")
print("="*70)
sys.path.insert(0, HERE)
from symbolproof_engine import run_proof
for s in (1,2,3,4):
    run_proof(s)

print()
print("="*70)
print("PART 4: TASK 3 -- single-valued L3; DECAGON CANONICAL PAIR")
print("="*70)
mp.mp.dps = 80
def Li3(x): return mp.re(mp.polylog(3, x))
def Li2(x): return mp.re(mp.polylog(2, x))
def L3mpf(q):
    if q > 1: return L3mpf(1/q)
    return Li3(q) - mp.log(abs(q))*Li2(q) - mp.mpf(1)/3*mp.log(abs(q))**2*mp.log(1-q)
def mpfq(q): q = Fr(q); return mp.mpf(q.numerator)/mp.mpf(q.denominator)
def L3(q): return L3mpf(mpfq(q))
z3 = mp.zeta(3)
s123 = pickle.load(open(os.path.join(O9B,"sol_B123.pkl"),"rb"))
lists = pickle.load(open(os.path.join(O9B,"Sunit_lists.pkl"),"rb"))
Gargs = [Fr(a) for a in lists['Gargs']]
invc = {1: Fr(113678924735323,3700404000000), 2: Fr(1978761997948319,41624512450560), 3: Fr(15326623217074655,142615945728000)}
invd = {1: Fr(34497,15625), 2: Fr(2085,2197), 3: Fr(-57300,24389)}
for e in (1,2,3):
    Xc = mp.fsum([mpfq(c)*L3(z) for c,z in zip([Fr(s) for s in s123[e]['sol'][:len(Gargs)]], Gargs) if c]) + mpfq(invc[e])*z3
    Xd = mp.fsum([mpfq(c)*L3(z) for c,z in zip([Fr(s) for s in sf['sols123'][e][:8]], G8) if c]) + mpfq(invd[e])*z3
    print(f"edge {e}: L3 OCT9c == OCT9d: {mp.nstr(Xc, 25)} ; |diff| = {mp.nstr(abs(Xc-Xd), 4)}  [NUMERICAL]")
# edge 0 trade
mp.mp.dps = 150
def Li3(x): return mp.re(mp.polylog(3, x))
def Li2(x): return mp.re(mp.polylog(2, x))
def L3(q):
    q = mpfq(q)
    return Li3(q) - mp.log(abs(q))*Li2(q) - mp.mpf(1)/3*mp.log(abs(q))**2*mp.log(abs(1-q))
s0 = pickle.load(open(os.path.join(O9B,"sol_B0.pkl"),"rb"))
sols = [Fr(s) for s in s0['meta']['sols']]
solB0 = [Fr(s) for s in s0['solB0']]
Sc = mp.fsum([mpfq(c)*L3(z) for c,z in zip(solB0[:len(sols)], sols) if c])
B0a = [Fr(s) for s in sf['B0']]
sol0 = [Fr(s) for s in sf['sol0']]
Sd = mp.fsum([mpfq(c)*L3(z) for c,z in zip(sol0[:20], B0a) if c])
z3 = mp.zeta(3)
invc0 = Fr(30167467635454913,654492976038000); invd0 = Fr(-7340,5043)
r = ((Sc + mpfq(invc0)*z3) - (Sd + mpfq(invd0)*z3))/z3
rel = mp.pslq([r, mp.mpf(1)], tol=mp.mpf('1e-120'), maxcoeff=10**20)
print(f"edge 0: (Xc-Xd)/zeta(3) = {mp.nstr(r, 25)} ; pslq -> {rel}  ({Fr(rel[1])/Fr(rel[0])*-1 if rel else 'not rational'})  [NUMERICAL]")
q = (Sc-Sd)/z3
rel = mp.pslq([q, mp.mpf(1)], tol=mp.mpf('1e-120'), maxcoeff=10**30)
print(f"edge 0: (Sc-Sd)/zeta(3) rational: {Fr(rel[1])/Fr(rel[0])*-1 if rel else 'NO'}  [NUMERICAL]")
# decagon canonical pair
mp.mp.dps = 120
L5 = (mp.zeta(3, mp.mpf(1)/5) - mp.zeta(3, mp.mpf(2)/5) - mp.zeta(3, mp.mpf(3)/5) + mp.zeta(3, mp.mpf(4)/5))/125
phi = (1+mp.sqrt(5))/2
X1, Xm1 = L3mpf(1/phi), L3mpf(-1/phi)
s5 = sp.sqrt(5)
c1 = sp.Rational(-75,32) + sp.Rational(81,64)*s5
c2 = sp.Rational(-35,32) + sp.Rational(53,64)*s5
c3 = sp.Rational(33,64) - sp.Rational(19,80)*s5
c3loc = sp.Rational(5,32) - sp.Rational(1,8)*s5
a = sp.Rational(1,64) + sp.Rational(31,320)*s5
b = sp.Rational(875,768) - sp.Rational(125,192)*s5
def ev(cf, val):
    return mp.mpf(str(sp.N(cf.subs(s5, mp.sqrt(5)), 40)))*val
valN = ev(c1,X1)+ev(c2,Xm1)+ev(c3-c3loc,z3)
valab = ev(a,z3)+ev(b,L5)
print(f"decagon canonical pair: N3 = ({a}) zeta(3) + ({b}) L(chi5,3)")
print(f"  val(N3) - (a zeta3 + b L5) = {mp.nstr(valN-valab, 6)}  [NUMERICAL]")
# golden identity at two precisions
for dps in (80, 150):
    mp.mp.dps = dps
    def Li3(x): return mp.re(mp.polylog(3, x))
    L5x = (mp.zeta(3, mp.mpf(1)/5) - mp.zeta(3, mp.mpf(2)/5) - mp.zeta(3, mp.mpf(3)/5) + mp.zeta(3, mp.mpf(4)/5))/125
    lphi = mp.log(phi)
    chi3g = (Li3(1/phi)-Li3(-1/phi))/2
    rhs = mp.sqrt(5)*mp.mpf(25)/48*L5x - mp.pi**2*lphi/12 + lphi**3/4
    print(f"  golden chi_3 identity residual at {dps} dps: {mp.nstr(abs(chi3g-rhs), 4)}  [NUMERICAL]")
print()
print("verify_pgen1.py finished in", round(time.time()-t0,1), "s")
'''
open("verify_pgen1.py","w").write(verify)
print("verify_pgen1.py written:", len(verify), "chars")
