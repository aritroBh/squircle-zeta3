"""DOD4 library: absolute closure of weight-3/2 constants over F=Q(sqrt3), S={p2,p3}, unit eps.
Builds on dod3_symbol.py (exact Q(sqrt3), S-unit factorization, Goncharov symbols).
Provides:
 - exact symbol -> [8 Li3 gens + 19 products] decomposition (weight 3)
 - exact symbol -> [3 Li2 gens + log*log products] decomposition (weight 2)
 - high-precision constants and evaluation helpers
 - full-precision PSLQ closure against the invisible space
   I3 = <z3, L(3,chi12), z2*l2, z2*l3, z2*lu, L2chi*l2, L2chi*l3, L2chi*lu>   (lu = log(2+sqrt3) > 0)
   I2 = <z2, L(2,chi12)>
"""
import sys, os, pickle
from fractions import Fraction as Fr
import mpmath as mp
import sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dod3_symbol import *

D3 = "<project>/research/08_theory/DOD3_dodecagon_symbol/"
D1 = "<project>/research/08_theory/DOD1_dodecagon/"
D2 = "<project>/research/08_theory/DOD2_dodecagon_finish/"

def mpf_of(z):
    m,n,d = z
    return (mp.mpf(m) + n*mp.sqrt(3))/d

# ---------------- exact data ----------------
_dat = pickle.load(open(D3+'li3_gens.pkl','rb'))
ZS, GENS = _dat['zs'], _dat['gens']          # 129 S-unit solutions + their Li3 symbols (1x27)
_pd = pickle.load(open(D3+'prod_sym_data.pkl','rb'))
PROD_SYMS, PROD_NAMES = _pd['prod_syms'], _pd['prod_names']   # 19 products (10 log^3 + 9 Li2*log)
B2ZS = pickle.load(open(D3+'li2_basis.pkl','rb'))['b2_zs']    # 3 Li2 basis args
SOLS = pickle.load(open(D3+'solsF_B8.pkl','rb'))
ZS129 = sorted(SOLS.keys(), key=pfloat)

# product metadata aligned with PROD_NAMES: ('log',(a,b,c)) uses lp2=(l2+lu)/2, lp3=l3/2, le=lu
PROD_META = []
for a in range(4):
    for b in range(4-a):
        c = 3-a-b; PROD_META.append(('log',(a,b,c)))
for w in B2ZS:
    for p in ('2','3','e'): PROD_META.append(('li2',w,p))
assert [str(x) for x in PROD_META] is not None

# ---------------- weight-3 quotient machinery ----------------
Q = Quotient()
PM = sp.Matrix.hstack(*[Q.coord(s_).T for s_ in PROD_SYMS])   # 27 x 19
DM_full = sp.Matrix.hstack(*[d.T for d in Q.D])
_M = sp.zeros(27,0)
for j in range(DM_full.cols):
    M2 = sp.Matrix.hstack(_M, DM_full[:,j]) if _M.cols else DM_full[:,j]
    if _M.cols==0 or M2.rank()>_M.rank(): _M=M2
    if _M.rank()==19: break
KEEP=[]
for i in range(27):
    e = sp.zeros(27,1); e[i]=1
    M2 = sp.Matrix.hstack(_M,e)
    if M2.rank()>_M.rank(): KEEP.append(i); _M=M2
    if len(KEEP)==8: break
AINV = _M.inv()
def pi3(v):
    return sp.Matrix((AINV*v.T)[19:,:])
GP = [pi3(g) for g in GENS]
def _niceness(i):
    z = ZS[i]; _,v2,v3,k = factor_sunit(z); fz=pfloat(z)
    return (0 if 0<fz<1 else (1 if fz<0 else 2), abs(v2)+abs(v3)+abs(k))
_order = sorted(range(len(ZS)), key=_niceness)
BIDX=[]; BM = sp.zeros(8,0)
for i in _order:
    M2 = sp.Matrix.hstack(BM,GP[i]) if BM.cols else GP[i]
    if BM.cols==0 or M2.rank()>BM.rank(): BIDX.append(i); BM=M2
    if len(BIDX)==8: break
BASIS3_ZS = [ZS[i] for i in BIDX]

def decompose_w3(T_sym):
    """T_sym dict -> (li3_coeffs[8] over BASIS3_ZS, prod_coeffs[19], exact tensor residual==0)."""
    Tc = Q.coord(T_sym)
    pT = pi3(Tc)
    solT = BM.gauss_jordan_solve(pT)[0]
    assert BM*solT == pT
    diff = Tc.T
    for j,i in enumerate(BIDX):
        diff -= sp.Rational(solT[j])*GENS[i].T
    solP = PM.gauss_jordan_solve(diff)[0]
    assert PM*solP == diff
    resid = diff - PM*solP
    assert resid == sp.zeros(27,1)
    return [sp.Rational(x) for x in solT], [sp.Rational(x) for x in solP]

# ---------------- weight-2 machinery ----------------
B3L = list(__import__('itertools').product(('2','3','e'), repeat=2))
_B9 = {t:i for i,t in enumerate(B3L)}
def coord2(S):
    v = sp.zeros(1,9)
    for t,c in S.items():
        v[_B9[t]] += sp.Rational(c.numerator, c.denominator)
    return v
def S_Li2(z):
    zv = vecfun(z); omzv = vecfun(psub(ONE,z))
    out = {}
    for a,ca in omzv.items():
        for b,cb in zv.items():
            out[(a,b)] = out.get((a,b),Fr(0)) - ca*cb
    return {k:x for k,x in out.items() if x!=0}
_D2 = []
for p in ('2','3','e'):
    for q in ('2','3','e'):
        S = {(p,q):Fr(1),(q,p):Fr(1)}
        if p==q: S = {(p,q):Fr(2)}
        _D2.append(coord2(S))
DM2 = sp.Matrix.hstack(*[d.T for d in _D2])   # 9 x 9, rank 6 (Sym^2)
_g2 = [coord2(S_Li2(z)) for z in ZS129]
# select 6 independent columns of DM2 (rank 6), then append 3 unit vectors -> 9x9
_M2 = sp.zeros(9,0)
for j in range(DM2.cols):
    M2x = sp.Matrix.hstack(_M2, DM2[:,j]) if _M2.cols else DM2[:,j]
    if _M2.cols==0 or M2x.rank()>_M2.rank(): _M2=M2x
    if _M2.rank()==6: break
_keep2=[]
for i in range(9):
    e = sp.zeros(9,1); e[i]=1
    M2x = sp.Matrix.hstack(_M2,e)
    if M2x.rank()>_M2.rank(): _keep2.append(i); _M2=M2x
    if len(_keep2)==3: break
AINV2 = _M2.inv()
def pi2(v): return sp.Matrix((AINV2*v.T)[6:,:])
GP2 = [pi2(g) for g in _g2]
def _niceness2(i):
    z = ZS129[i]; _,v2,v3,k = factor_sunit(z); fz=pfloat(z)
    return (0 if 0<fz<1 else (1 if fz<0 else 2), abs(v2)+abs(v3)+abs(k))
_order2 = sorted(range(len(ZS129)), key=_niceness2)
B2IDX=[]; BM2 = sp.zeros(3,0)
for i in _order2:
    M2x = sp.Matrix.hstack(BM2,GP2[i]) if BM2.cols else GP2[i]
    if BM2.cols==0 or M2x.rank()>BM2.rank(): B2IDX.append(i); BM2=M2x
    if len(B2IDX)==3: break
BASIS2_ZS = [ZS129[i] for i in B2IDX]

def decompose_w2(T_sym2):
    """weight-2 symbol dict -> (li2_coeffs[3] over BASIS2_ZS, loglog_coeffs[9] over B3L products)."""
    Tc = coord2(T_sym2)
    pT = pi2(Tc)
    solT = BM2.gauss_jordan_solve(pT)[0]
    assert BM2*solT == pT
    diff = Tc.T
    for j,i in enumerate(B2IDX):
        diff -= sp.Rational(solT[j])*_g2[i].T
    solP = DM2.gauss_jordan_solve(diff)[0]
    # DM2 has rank 6 (9x9 singular): solution has free parameters -> set them to 0
    free = sorted(solP.free_symbols, key=str)
    solP = solP.subs({s: sp.Integer(0) for s in free})
    assert DM2*solP == diff
    return [sp.Rational(x) for x in solT], [sp.Rational(x) for x in solP]

# ---------------- numeric ----------------
def consts(dps):
    mp.mp.dps = dps
    r3 = mp.sqrt(3)
    l2,l3,lu = mp.log(2), mp.log(3), mp.log(2+r3)   # lu = log(2+sqrt3) > 0
    z2, z3 = mp.pi**2/6, mp.zeta(3)
    L3c = (mp.zeta(3, mp.mpf(1)/12) - mp.zeta(3, mp.mpf(5)/12) - mp.zeta(3, mp.mpf(7)/12) + mp.zeta(3, mp.mpf(11)/12)) / 12**3
    L2c = (mp.zeta(2, mp.mpf(1)/12) - mp.zeta(2, mp.mpf(5)/12) - mp.zeta(2, mp.mpf(7)/12) + mp.zeta(2, mp.mpf(11)/12)) / 12**2
    return dict(l2=l2,l3=l3,lu=lu,z2=z2,z3=z3,L3c=L3c,L2c=L2c,lp2=(l2+lu)/2,lp3=l3/2,le=lu,r3=r3)

def Li3v(z, dps=None):
    return mp.re(mp.polylog(3, mp.mpc(mpf_of(z))))
def Li2v(z, dps=None):
    return mp.re(mp.polylog(2, mp.mpc(mpf_of(z))))

def prodval(mta, C):
    if mta[0]=='log':
        a,b,c = mta[1]; return C['lp2']**a*C['lp3']**b*C['le']**c
    _,w,p = mta
    return Li2v(w)*{'2':C['lp2'],'3':C['lp3'],'e':C['le']}[p]

def invisible_basis3(C):
    """returns (names, values) of I3 basis."""
    names = ['z3','L3chi','z2*l2','z2*l3','z2*lu','L2chi*l2','L2chi*l3','L2chi*lu']
    vals = [C['z3'],C['L3c'],C['z2']*C['l2'],C['z2']*C['l3'],C['z2']*C['lu'],
            C['L2c']*C['l2'],C['L2c']*C['l3'],C['L2c']*C['lu']]
    return names, vals

def rat(x):
    return mp.mpf(sp.Rational(x).numerator)/mp.mpf(sp.Rational(x).denominator)


# ---------------- invisible space (10-dim over Q = F-span of 5) ----------------
def invisible10(C):
    """Q-basis {X, sqrt3*X : X in {z3, L3chi, z2*l2, z2*l3, z2*lu}}; L(2,chi12)=z2/r3 collapsed."""
    r3 = C['r3']
    items = [('z3',C['z3']), ('L3chi',C['L3c']), ('z2*l2',C['z2']*C['l2']),
             ('z2*l3',C['z2']*C['l3']), ('z2*lu',C['z2']*C['lu'])]
    names, vals = [], []
    for nm, v in items:
        names.append(nm); vals.append(v)
        names.append('r3*'+nm); vals.append(r3*v)
    return names, vals

def value_of_decomp(solT, solP, C):
    """sum over basis gens Li3 + products, at current dps."""
    v = mp.mpf(0)
    for j,i in enumerate(BIDX):
        if solT[j]: v += rat(solT[j])*Li3v(ZS[i])
    for j in range(19):
        if solP[j]: v += rat(solP[j])*prodval(PROD_META[j], C)
    return v

def closure_w3(T_sym, T_val, dps_pslq=(620,450), maxcoeff=10**12, verbose=True):
    """Full closure: exact decomposition + remainder PSLQ vs INV10 at two precisions.
       T_val: mpf value at 620 dps (already computed). Returns dict."""
    solT, solP = decompose_w3(T_sym)
    rem = T_val - value_of_decomp(solT, solP, consts(620))
    names, ivals = invisible10(consts(620))
    out = {'solT':solT, 'solP':solP, 'rem':rem, 'names':names}
    rels = {}
    for dps in dps_pslq:
        tol = mp.mpf(10)**(-(dps-50))
        vv = mp.mpf(mp.nstr(rem, dps))
        ib = [mp.mpf(mp.nstr(x, dps)) for x in ivals]
        rel = mp.pslq([vv]+ib, tol=tol, maxcoeff=maxcoeff, maxsteps=200000)
        rels[dps] = rel
    out['rels'] = rels
    if verbose:
        print("remainder:", mp.nstr(rem, 40))
        for dps in dps_pslq:
            rel = rels[dps]
            if rel:
                resid = abs(rel[0]*rem + sum(rel[i+1]*ivals[i] for i in range(len(ivals))))
                print("  PSLQ@%d:"%dps, rel, " resid@620:", mp.nstr(resid,3))
            else:
                print("  PSLQ@%d: none"%dps)
    return out
