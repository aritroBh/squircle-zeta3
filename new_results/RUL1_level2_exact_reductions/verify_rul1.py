"""verify_rul1.py -- rebuilds the level-2 (DS+DUAL) relation system, reduces every target
exactly (with certificates), and checks numerically to >=100 digits. Target: < 5 min.
Run: nice python verify_rul1.py
"""
import sys, os, time, math, pickle, json
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import mpmath as mp
import sympy as sp

T00 = time.time()
RESULTS = []
def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    log(("PASS " if ok else "FAIL ") + name + ("   " + detail if detail else ""))

import rul1, rul1la, rul1_named, rul1num, rul1_reduce

# ---------------- 1. rebuild relations + reduction data ----------------
log("rebuilding relations (DS+DUAL) and RREF ...")
data = {}
DETS = {}
for w in range(1, 6):
    rels = rul1.relations_at_weight(w, families=('DS','DUAL'))
    cols = sorted(rul1.conv_words_of_weight(w))
    rows = [rul1la.Row(r, {i: Fr(1)}) for i, (fam, tag, r) in enumerate(rels)]
    rref, piv = rul1la.eliminate(rows, cols)
    free = [c for c in cols if c not in piv]
    M, hist_named = [], []
    for nm, combo in rul1_named.NAMED[w]:
        res, h = rul1la.reduce_vector(combo, rref)
        M.append(res); hist_named.append(h)
    fnames = [nm for nm, _ in rul1_named.NAMED[w]]
    Ms = sp.Matrix([[row.get(f, Fr(0)) for f in free] for row in M])
    data[w] = {'rels': rels, 'cols': cols, 'rref': [(r.d, r.h, r.pivot) for r in rref],
               'free': free, 'M': M, 'Minv': [[Ms.inv()[i,j] for j in range(len(free))] for i in range(len(free))],
               'hist_named': hist_named, 'fnames': fnames, 'det': Ms.det()}
    DETS[w] = (len(rels), len(piv), len(cols), len(free))
rul1_reduce._data = data
rul1_reduce._rref_cache = {}
codims = tuple(DETS[w][3] for w in range(1,6))
check("relations span codimension (1,2,3,5,8)", codims == (1,2,3,5,8), str(DETS))
log("rebuilt: %s" % (codims,))

D = data
mp.mp.dps = 120
NV = rul1num.named_values(120)
def nval(coords):
    return sum(mp.mpf(x.numerator)/mp.mpf(x.denominator)*NV[nm] for nm, x in coords.items())

# ---------------- 2. controls ----------------
c = rul1_reduce.reduce_combo(3, {(0,-1,-1): Fr(1)})
check("CONTROL zeta(bar2,1) = zeta(3)/8", c == {'z3': Fr(1,8)}, str(c))
c = rul1_reduce.reduce_combo(4, dict(rul1_named.NAMED[4])['Li4h'])
check("CONTROL Li_4(1/2) primitive", c == {'Li4h': Fr(1)}, str(c))
w = (0,0,-1,1)
c, h = rul1_reduce.reduce_combo(4, {w: Fr(1)}, want_history=True)
okcert = rul1_reduce.check_certificate(4, {w: Fr(1)}, c, h)
v1 = nval(c); v2 = rul1num.word_num_smart(w, dps=100)
check("CONTROL zeta(bar3,bar1) reduced (certificate) + 100-digit", okcert and abs(v1-v2) < mp.mpf('1e-100'),
      "%s; |diff|=%s" % (c, mp.nstr(abs(v1-v2),3)))
c = rul1_reduce.reduce_combo(3, {(0,-1,-1): Fr(2)})
check("CONTROL J_{2,0} = zeta(3)/4", c == {'z3': Fr(1,4)}, str(c))
# Euler sum control: sum (-1)^{n+1} H_n/n^2 = 5 zeta(3)/8 (via eta(3) - zeta(bar2,1) = 3/4-1/8)
check("CONTROL sum (-1)^{n+1} H_n/n^2 = 5/8 zeta(3)", Fr(3,4)-Fr(1,8) == Fr(5,8))

# ---------------- 3. Target 1: five J-values + N4, N5 ----------------
def Jword(a,k): return (0,)*(k+1) + (-1,)*a
EXP = {
 (3,0): {'log2^4': Fr(-1,4), 'pi2*log2^2': Fr(1,4), 'pi4': Fr(1,15), 'z3*log2': Fr(-21,4), 'Li4h': Fr(-6)},
 (2,1): {'log2^4': Fr(1,6), 'pi2*log2^2': Fr(-1,6), 'pi4': Fr(-1,24), 'z3*log2': Fr(7,2), 'Li4h': Fr(4)},
 (4,0): {'log2^5': Fr(-4,5), 'pi2*log2^3': Fr(2,3), 'z3*log2^2': Fr(-21,2), 'z5': Fr(24), 'Li4h*log2': Fr(-24), 'Li5h': Fr(-24)},
 (3,1): {'log2^5': Fr(2,5), 'pi2*log2^3': Fr(-1,3), 'z3*log2^2': Fr(21,4), 'z3*pi2': Fr(-1,2), 'z5': Fr(-99,16), 'Li4h*log2': Fr(12), 'Li5h': Fr(12)},
 (2,2): {'z3*pi2': Fr(1,3), 'z5': Fr(-29,8)},
}
certs = {}
for (a,k), exp in EXP.items():
    wt = a+k+1
    combo = {Jword(a,k): Fr(math.factorial(a)*math.factorial(k))}
    cc, hh = rul1_reduce.reduce_combo(wt, combo, want_history=True)
    ok = rul1_reduce.check_certificate(wt, combo, cc, hh) and cc == exp
    vn = nval(cc)
    vE1 = mp.mpf(math.factorial(a)*math.factorial(k))*rul1num.g_word_over_0m1_num(Jword(a,k))
    vE2 = rul1num.J_num(a,k)
    d = max(abs(vn-vE1), abs(vn-vE2))
    check("T1 J_{%d,%d} proved + 100-digit" % (a,k), ok and d < mp.mpf('1e-100'),
          "|diff|=%s" % mp.nstr(d,3))
    certs[(a,k)] = cc
def lincomb(spec):
    out = {}
    for (a,k), cf in spec:
        for nm, x in certs[(a,k)].items(): out[nm] = out.get(nm, Fr(0)) + cf*x
    return {nm: x for nm, x in out.items() if x}
check("T1 (1/12)J30+(1/8)J21 = pi^4/2880", lincomb([((3,0),Fr(1,12)),((2,1),Fr(1,8))]) == {'pi4': Fr(1,2880)})
check("T1 (1/48)(J40+2J31) = -z3 pi2/48 + 31 z5/128",
      lincomb([((4,0),Fr(1,48)),((3,1),Fr(2,48))]) == {'z3*pi2': Fr(-1,48), 'z5': Fr(31,128)})

# N4, N5 via D12's exact formula
KMAX = 5
def addC(A,B):
    C=dict(A)
    for k,v in B.items():
        C[k]=C.get(k,Fr(0))+v
        if C[k]==0: del C[k]
    return C
def mulC(A,B):
    C={}
    for (i1,j1,k1),v1 in A.items():
        for (i2,j2,k2),v2 in B.items():
            key=(i1+i2,j1+j2,k1+k2); C[key]=C.get(key,Fr(0))+v1*v2
    return {k:v for k,v in C.items() if v}
def serMul(A,B,K=KMAX):
    C=[{} for _ in range(K+1)]
    for i in range(len(A)):
        if not A[i]: continue
        for j in range(len(B)):
            if not B[j]: continue
            if i+j>K: break
            C[i+j]=addC(C[i+j],mulC(A[i],B[j]))
    return C
def serExp(A,K=KMAX):
    R=[{} for _ in range(K+1)]; R[0]={(0,0,0):Fr(1)}
    for n in range(1,K+1):
        acc={}
        for k in range(1,n+1):
            if A[k]: acc=addC(acc,mulC({kk:vv*Fr(k) for kk,vv in A[k].items()},R[n-k]))
        R[n]={kk:vv/n for kk,vv in acc.items()}
    return R
y = sp.symbols('y')
lc = sp.log(sp.cosh(y)).series(y,0,KMAX+2).removeO()
cc5 = {1: Fr(-1,2)}
for r in range(2,KMAX+1): cc5[r] = Fr(sp.expand(lc).coeff(y,r))/Fr(2)**r
T = [{} for _ in range(KMAX+1)]; T[0] = {(0,0,1): Fr(1)}
for r in range(1,KMAX+1):
    if cc5[r] == 0: continue
    for t in range(0,KMAX+1-r):
        T[r+t] = addC(T[r+t], {(r,0,0): Fr(math.comb(r+t-1,t))*cc5[r]})
Q = [{} for _ in range(KMAX+1)]
for n in range(KMAX+1):
    Q[n] = addC(Q[n], {kk:-vv for kk,vv in T[n].items()})
    if n+1 <= KMAX: Q[n+1] = addC(Q[n+1], {kk:-vv for kk,vv in T[n].items()})
psi2 = [{kk: vv/2 for kk,vv in co.items()} for co in serExp([{}]+Q[1:])]
B5 = [{} for _ in range(KMAX+1)]
for m in range(2,KMAX+1): B5[m] = {(1,0,0): Fr(-1)}
psi1 = serExp(B5)
psi3 = [{} for _ in range(KMAX+1)]
for m in range(1,KMAX+1): psi3[m] = {(0,m,0): Fr(1,math.factorial(m))}
PSI = serMul(serMul(psi1,psi2),psi3)

l2s, pi2s, z3s, z5s = sp.symbols('l2s pi2s z3s z5s')
def zeta(n):
    return {2: pi2s/6, 3: z3s, 4: pi2s**2/90, 5: z5s}[n]
def eta(n):
    return (1-Fr(2)**(1-n))*zeta(n)
J = {(1,k): math.factorial(k)*eta(k+2) for k in range(4)}
J[(2,0)] = z3s/4
sy = {'log2^4': l2s**4, 'log2^5': l2s**5, 'pi2*log2^2': pi2s*l2s**2, 'pi2*log2^3': pi2s*l2s**3,
      'pi4': pi2s**2, 'z3*log2': z3s*l2s, 'z3*log2^2': z3s*l2s**2, 'z3*pi2': z3s*pi2s,
      'z5': z5s, 'Li4h': sp.Integer(0), 'Li5h': sp.Integer(0), 'Li4h*log2': sp.Integer(0)}
# NOTE: Li4h, Li5h cancel in N4,N5; keep them as symbols and verify they do.
Li4hs, Li5hs = sp.symbols('Li4hs Li5hs')
sy['Li4h'] = Li4hs; sy['Li5h'] = Li5hs; sy['Li4h*log2'] = Li4hs*l2s
for (a,k), cc in certs.items():
    J[(a,k)] = sum(x*sy[nm] for nm, x in cc.items())
assert all((j,i) in J for m in range(5) for (i,j,k) in PSI[m])
INT = []
for m in range(5):
    s = 0
    for (i,j,k), v in PSI[m].items():
        s += v*J[(j,i)]*l2s**k
    INT.append(sp.expand(s))
N4 = sp.expand(-l2s**5/120 + sum(INT[0:4]))
N5 = sp.expand(+l2s**6/720 + sum(INT[0:5]))
N4_exp = -l2s**5/120 + pi2s/24*(1-l2s+l2s**2/2) + z3s/4*(1-l2s) + pi2s**2/2880
N5_exp = (l2s**6/720 + pi2s/24*(1-l2s+l2s**2/2-l2s**3/6) + z3s/4*(1-l2s+l2s**2/2)
          + pi2s**2/2880*(1-l2s) - 7*pi2s**2/1920 - z3s*pi2s/48 + z5s/8)
d4 = sp.expand(N4-N4_exp); d5 = sp.expand(N5-N5_exp)
check("T1 N4 closed form exact", d4 == 0, str(N4))
check("T1 N5 closed form exact", d5 == 0, str(N5))
check("T1 Li4h,Li5h cancel in N4,N5", not N4.has(Li4hs) and not N5.has(Li4hs) and not N4.has(Li5hs) and not N5.has(Li5hs))
# numeric: J-quad assembly
mp.mp.dps = 130
def Jq(a,k): return mp.quad(lambda s: s**k*mp.log(1+mp.e**(-s))**a, [0, mp.inf])
INTq = []
for m in range(5):
    s = mp.mpf(0)
    for (i,j,k), v in PSI[m].items():
        s += mp.mpf(v.numerator)/mp.mpf(v.denominator)*Jq(j,i)*mp.log(2)**k
    INTq.append(s)
NV2 = rul1num.named_values(130)
N4q = (-mp.log(2))**5/120 + sum(INTq[0:4]); N5q = (mp.log(2))**6/720 + sum(INTq[0:5])
N4v = mp.mpf(str(N4_exp.subs({l2s: NV2['log2'], pi2s: NV2['pi2'], z3s: NV2['z3'], z5s: NV2['z5']})))
N5v = mp.mpf(str(N5_exp.subs({l2s: NV2['log2'], pi2s: NV2['pi2'], z3s: NV2['z3'], z5s: NV2['z5']})))
check("T1 N4 100-digit vs quad assembly", abs(N4q-N4v) < mp.mpf('1e-110'), mp.nstr(abs(N4q-N4v),3))
check("T1 N5 100-digit vs quad assembly", abs(N5q-N5v) < mp.mpf('1e-110'), mp.nstr(abs(N5q-N5v),3))
log("Target 1 done (%.0fs)" % (time.time()-T00))

# ---------------- 4. Target 2: 13 atoms + folds ----------------
ATOMWORDS = {'H1':(1,), 'Li2':(0,1), 'Li3':(0,0,1), 'H011':(0,1,1), 'Li4':(0,0,0,1),
             'H0011':(0,0,1,1), 'H0111':(0,1,1,1), 'Li5':(0,0,0,0,1), 'H00011':(0,0,0,1,1),
             'H00101':(0,0,1,0,1), 'H00111':(0,0,1,1,1), 'H01011':(0,1,0,1,1), 'H01111':(0,1,1,1,1)}
def atom_combo(nm):
    hw = ATOMWORDS[nm]
    gw = tuple(-1 if a == 1 else 0 for a in hw)
    return {gw: Fr((-1)**sum(1 for a in hw if a == 1))}
EXP2 = {
 'H1': {'log2': Fr(-1)}, 'Li2': {'pi2': Fr(-1,12)}, 'Li3': {'z3': Fr(-3,4)}, 'H011': {'z3': Fr(1,8)},
 'Li4': {'pi4': Fr(-7,720)},
 'H0011': {'pi4': Fr(-1,48), 'z3*log2': Fr(7,4), 'pi2*log2^2': Fr(-1,12), 'log2^4': Fr(1,12), 'Li4h': Fr(2)},
 'H0111': {'pi4': Fr(-1,90), 'z3*log2': Fr(7,8), 'pi2*log2^2': Fr(-1,24), 'log2^4': Fr(1,24), 'Li4h': Fr(1)},
 'Li5': {'z5': Fr(-15,16)},
 'H00011': {'z5': Fr(-29,32), 'z3*pi2': Fr(1,12)},
 'H00101': {'z5': Fr(41,32), 'z3*pi2': Fr(-5,48)},
 'H00111': {'z5': Fr(33,32), 'z3*pi2': Fr(1,12), 'z3*log2^2': Fr(-7,8), 'pi2*log2^3': Fr(1,18), 'log2^5': Fr(-1,15), 'Li5h': Fr(-2), 'Li4h*log2': Fr(-2)},
 'H01011': {'z5': Fr(-73,64), 'z3*pi2': Fr(-1,4), 'z3*log2^2': Fr(7,4), 'pi2*log2^3': Fr(-1,9), 'log2^5': Fr(2,15), 'Li5h': Fr(4), 'Li4h*log2': Fr(4)},
 'H01111': {'z5': Fr(1), 'z3*log2^2': Fr(-7,16), 'pi2*log2^3': Fr(1,36), 'log2^5': Fr(-1,30), 'Li5h': Fr(-1), 'Li4h*log2': Fr(-1)},
}
mp.mp.dps = 120
for nm in ATOMWORDS:
    wt = len(ATOMWORDS[nm])
    combo = atom_combo(nm)
    cc, hh = rul1_reduce.reduce_combo(wt, combo, want_history=True)
    ok = rul1_reduce.check_certificate(wt, combo, cc, hh) and cc == EXP2[nm]
    vn = nval(cc)
    vE1 = sum(mp.mpf(x.numerator)/mp.mpf(x.denominator)*rul1num.g_word_over_0m1_num(u) for u, x in combo.items())
    d = abs(vn-vE1)
    check("T2 atom %-7s proved + 100-digit" % nm, ok and d < mp.mpf('1e-100'), "|diff|=%s" % mp.nstr(d,3))

# folds
K1G = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
coef_acc = pickle.load(open(f'{K1G}/coef_acc_raw.pkl','rb'))
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
l2, pi2, z3, z5, Li4h, Li5h, sq5 = sp.symbols('l2 pi2 z3 z5 Li4h Li5h sq5')
ATOMW = {nm: len(w) for nm, w in ATOMWORDS.items()}
ATOMC = {nm: atom_combo(nm) for nm in ATOMWORDS}
ZETAC = {0: {(0,1): Fr(-1)}, 1: {(0,0,1): Fr(-1)}, 2: {(0,0,0,0,1): Fr(-1)}}
ZW = {0: 2, 1: 3, 2: 5}
from functools import lru_cache
@lru_cache(None)
def shpow(w, n):
    out = {(): Fr(1)} if n == 0 else {w: Fr(1)}
    for _ in range(n-1): out = rul1.sh_mul(out, {w: Fr(1)})
    return out
@lru_cache(None)
def reduce_cached(w, key):
    if w == 0:
        d = dict(key); return (('1', d.get((), Fr(0))),)
    return tuple(sorted(rul1_reduce.reduce_combo(w, dict(key)).items()))
def nsym(coords):
    expr = 0
    for nm, x in coords:
        if nm == '1': expr += x
        elif nm == 'pi4': expr += x*pi2**2
        elif nm == 'log2': expr += x*l2
        elif nm.startswith('log2^'): expr += x*l2**int(nm[5:])
        elif nm == 'pi2': expr += x*pi2
        elif nm == 'z3*pi2': expr += x*z3*pi2
        elif nm == 'z3': expr += x*z3
        elif nm == 'z5': expr += x*z5
        elif nm == 'Li4h': expr += x*Li4h
        elif nm == 'Li5h': expr += x*Li5h
        elif nm == 'pi2*log2^2': expr += x*pi2*l2**2
        elif nm == 'pi2*log2^3': expr += x*pi2*l2**3
        elif nm == 'z3*log2': expr += x*z3*l2
        elif nm == 'z3*log2^2': expr += x*z3*l2**2
        elif nm == 'pi2*log2': expr += x*pi2*l2
        elif nm == 'pi4*log2': expr += x*pi2**2*l2
        elif nm == 'Li4h*log2': expr += x*Li4h*l2
        else: raise ValueError(nm)
    return expr
def fold(rd, only_w=None):
    total = sp.Integer(0)
    for (at, xe, ze), co in coef_acc.items():
        cov = sp.expand(co.subs(rd))
        if cov == 0: continue
        base = {(): Fr(1)}
        zw = 0
        for nm, e in at:
            for _ in range(e): base = rul1.sh_mul(base, ATOMC[nm])
            zw += ATOMW[nm]*e
        for i, e in enumerate(ze):
            for _ in range(e): base = rul1.sh_mul(base, ZETAC[i])
            zw += ZW[i]*e
        poly = sp.Poly(cov, l2) if cov.has(l2) else None
        terms = poly.terms() if poly else [((0,), cov)]
        for (m2,), cm in terms:
            if only_w is not None and zw + m2 not in only_w: continue
            combo = dict(base)
            if m2: combo = rul1.sh_mul(combo, shpow((-1,), m2))
            red = reduce_cached(zw+m2, tuple(sorted(combo.items())))
            total += cm*(-1)**xe*nsym(red)
    return sp.expand(120*total)
rd1 = {A0:sp.Rational(1,2),A1:sp.Rational(3,4),B:sp.Rational(1,4),kap:sp.Integer(1),
       L2:sp.Integer(5),q0:sp.Integer(1),V02:sp.Integer(1),lA0:-l2,lA1:-2*l2,lB:-2*l2}
val5 = fold(rd1, only_w={5}).subs({sp.sqrt(5): sq5})
terms5 = sp.Add.make_args(sp.expand(val5))
cz5 = sp.simplify(sum(t for t in terms5 if t.has(z5) and not t.has(z3) and not t.has(pi2) and not t.has(l2) and not t.has(Li4h) and not t.has(Li5h))/z5)
cpi2z3 = sp.simplify(sum(t for t in terms5 if t.has(z3) and t.has(pi2) and not t.has(z5) and not t.has(l2) and not t.has(Li4h) and not t.has(Li5h))/(pi2*z3))
check("T2 zeta5 coord = -20187/2000 sqrt5", sp.simplify(cz5 + sp.Rational(20187,2000)*sq5) == 0, str(cz5))
check("T2 zeta2zeta3 coord = -858/125 sqrt5", sp.simplify(6*cpi2z3 + sp.Rational(858,125)*sq5) == 0, str(sp.simplify(6*cpi2z3)))
nl5 = sp.simplify(cz5 + sp.Rational(2079,250)*sq5)
nl23 = sp.simplify(6*cpi2z3 + sp.Rational(477,250)*sq5)
check("T2 far-kink zeta5 = -711/400 sqrt5", sp.simplify(nl5 + sp.Rational(711,400)*sq5) == 0, str(nl5))
check("T2 far-kink zeta2zeta3 = -1239/250 sqrt5", sp.simplify(nl23 + sp.Rational(1239,250)*sq5) == 0, str(nl23))
log("Target 2 zeta-coords done (%.0fs)" % (time.time()-T00))

rdB = {A0:sp.Rational(1,2),A1:sp.Rational(3,4),B:sp.Rational(1,4),kap:sp.Integer(1),
       L2:sp.Integer(5),q0:sp.Integer(1),V02:sp.Integer(1),lA0:l2,lA1:-l2/2,lB:0}
valB = fold(rdB).subs({sp.sqrt(5): sq5})
X_e5e = (-153*l2 - sp.Rational(32,5)*pi2 - sp.Rational(26817,80)*l2**2 + sp.Rational(3267,50)*z3
         + sp.Rational(291,20)*pi2*l2 - sp.Rational(386517,3200)*l2**3 - sp.Rational(6489,8000)*pi2**2
         - sp.Rational(20679,400)*z3*l2 - sp.Rational(18823,3200)*pi2*l2**2 - sp.Rational(382819,10240)*l2**4
         - sp.Rational(1242,25)*Li4h - sp.Rational(20187,2000)*z5 - sp.Rational(858,125)*pi2*z3/6
         + sp.Rational(6463,40000)*pi2**2*l2 + sp.Rational(57,250)*z3*l2**2 + sp.Rational(31267,24000)*pi2*l2**3
         + sp.Rational(992499,2560000)*l2**5 + sp.Rational(396,125)*Li5h + sp.Rational(918,125)*Li4h*l2)
check("T2 conv-B fold == sqrt5 * E5e closed form (all 20 constants)",
      sp.simplify(sp.expand(valB - sq5*X_e5e)) == 0)
# numeric vs E5e's 337-digit value
mp.mp.dps = 130
valB_num = sp.N(valB.subs({l2: NV2['log2'], pi2: NV2['pi2'], z3: NV2['z3'], z5: NV2['z5'],
                           Li4h: NV2['Li4h'], Li5h: NV2['Li5h'], sq5: mp.sqrt(5)}), 125)
e1 = json.load(open('<project>/research/08_theory/E5d_nonlocal_pslq/b0_p5_convB_350.json'))['e1']['P5']
de = abs(mp.mpf(str(valB_num)) - mp.mpf(e1))
check("T2 conv-B fold 100-digit vs E5e data", de < mp.mpf('1e-100'), "|diff|=%s" % mp.nstr(de,3))
log("Target 2 fold done (%.0fs)" % (time.time()-T00))

# ---------------- 5. Target 3: regular hexagon ----------------
K1F = '<project>/research/08_theory/K1f_hexagon_order5'
bb = pickle.load(open(f'{K1G}/build_blocks.pkl','rb'))
JLOC_n = bb['JLOC_n']
E5 = pickle.load(open(f'{K1F}/E5_polys.pkl','rb'))[5]
l3, sq3 = sp.symbols('l3 sq3')
rd = {A0: sp.Rational(1,2), A1: sp.Rational(1,2), B: sp.Rational(1,2), kap: sp.Integer(0),
      L2: sp.Rational(1,3), q0: sp.Rational(-1,6), V02: sp.Rational(1,3), lA0: -l2, lA1: -l2, lB: -l2}
z2g, z3g, z5g = sp.symbols('z2 z3 z5')
Earc = sp.Integer(0)
for (n,a,b,c), co in E5.items():
    if c != 0: continue
    cov = co.subs(rd)
    if cov == 0: continue
    jv = JLOC_n[(n,a,b)].subs({z2g: pi2/6, z3g: z3, z5g: z5})
    Earc += cov*jv
P5 = sp.expand(720*Earc)
RHS = (30720*l2 - 4920*pi2 - 61440*l2**2 + 6240*pi2*l2 + 30720*l2**3 + 448*pi2**2
       + 5760*z3*l2 - 1920*pi2*l2**2 - 5120*l2**4 + 11340*z5 + 5580*pi2*z3/6
       - 123*pi2**2*l2 - 1440*z3*l2**2 + 160*pi2*l2**3 + 256*l2**5)
check("T3 (128/sqrt3) P^(5)(1) == paper's 15-term closed form",
      sp.expand(128/sp.sqrt(3)*P5 - RHS) == 0)
P5e = sp.expand(P5)
cz5r = sp.simplify(sum(t for t in sp.Add.make_args(P5e) if t.has(z5) and not t.has(z3) and not t.has(pi2) and not t.has(l2))/z5)
cpi2z3r = sp.simplify(sum(t for t in sp.Add.make_args(P5e) if t.has(z3) and t.has(pi2) and not t.has(z5) and not t.has(l2))/(pi2*z3))
check("T3 zeta5 = 2835 sqrt3/32", sp.simplify(cz5r - 2835*sp.sqrt(3)/32) == 0, str(cz5r))
check("T3 zeta2zeta3 = 1395 sqrt3/32", sp.simplify(6*cpi2z3r - 1395*sp.sqrt(3)/32) == 0, str(6*cpi2z3r))
d17 = json.load(open('<project>/research/08_theory/D17_zeta5_polygons/m3_total.json'))['P5']
P5_num = sp.N(P5.subs({l2: NV2['log2'], pi2: NV2['pi2'], z3: NV2['z3'], z5: NV2['z5']}), 125)
dP = abs(mp.mpf(str(P5_num)) - mp.mpf(d17))
check("T3 P5 100-digit vs D17 data", dP < mp.mpf('1e-90'), "|diff|=%s" % mp.nstr(dP,3))
log("Target 3 done (%.0fs)" % (time.time()-T00))

# ---------------- summary ----------------
log("="*70)
npass = sum(1 for _, ok, _ in RESULTS if ok)
log("SUMMARY: %d/%d checks PASS  (total %.0fs)" % (npass, len(RESULTS), time.time()-T00))
for name, ok, detail in RESULTS:
    if not ok: log("  FAILED: " + name)
assert npass == len(RESULTS)
log("ALL CHECKS PASS")
