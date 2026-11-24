script = r'''"""verify_dod5.py -- independent verification of the DOD5 closed form of P'''(1) for the regular 12-gon.
Loads the EXACT (proved/symbol-level + two-precision-PSLQ) data produced by the DOD5 steps:
  dod5_closures36.pkl      per-constant exact decomposition (solT,solP) + invisible relation (rel)  [36 constants]
  dod5_all_li3_reductions.pkl  same for the 23 raw Li3 points
  dod5_tallies.pkl         merged exact Q(sqrt3) coefficients over the 8-arg Li3 basis / 19 products / 5 constants
  dod5_R_content.pkl       sum c^T3a * R_sym exact content (Li3 12 pts, Li2*log, log^3)
  dod5_rawli3_coeffs.pkl   raw Li3 coefficients in S
  dod5_atom_coefs_all.pkl  exact Q(sqrt3) (and *log) coefficients of all 60 atoms in S
and recomputes ALL numbers from scratch at 330 digits, then checks:
  [1] each of the 36 depth-2 constants satisfies its closed identity (residual at 330 dps)
  [2] each of the 23 Li3-point reductions (residual)
  [3] all 12 R pieces: symbolic build_R vs T3a-quadrature minus swaps (residual)
  [4] Xi36 closed form vs the 36-term sum of values
  [5] full assembly: 72*L*S_calc vs P3_630.pkl   (target: residual < 1e-290, i.e. >= 300 digits)
Runtime target < 5 min.
"""
import sys, os, pickle, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
t0 = time.time()
import mpmath as mp
mp.mp.dps = 330
DPS = 330
from dod4_lib import (P, ONE, pdiv, pmul, psub, pfloat, psplit, factor_sunit,
                      consts, Li3v, prodval, PROD_META, BASIS3_ZS, invisible10,
                      D1, D2, D3)
import sympy as sp
from dod3_R import build_R, TRIPLES

C = consts(DPS)
r3 = C['r3']
LI3B = [Li3v(z) for z in BASIS3_ZS]
PBV = [prodval(m, C) for m in PROD_META]
IBV = [C['z3'], C['L3c'], C['z2']*C['l2'], C['z2']*C['l3'], C['z2']*C['lu']]
names10, ivals10 = invisible10(C)
def rat(x):
    x = sp.Rational(x); return mp.mpf(x.numerator)/mp.mpf(x.denominator)
def fval(p):
    return mp.mpf(p[0].numerator)/p[0].denominator + r3*mp.mpf(p[1].numerator)/p[1].denominator
ZERO=(Fr(0),Fr(0))

# ---------- numeric values of the 36 constants ----------
Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}
def f_prefix_series(x, y, dps=DPS):
    mp.mp.dps = dps+20
    x = mpf_of(x); y = mpf_of(y)
    tot = mp.mpf(0); ps = mp.mpf(0); xp = mp.mpf(1); yp = mp.mpf(1)
    for m in range(1, 400000):
        xp *= x
        if m >= 2: tot += xp*ps/m**2
        yp *= y; ps += yp/m
        if m > 80 and abs(xp*max(1,abs(ps))) < mp.mpf(10)**(-(dps+15)): break
    mp.mp.dps = dps
    return +tot
Li2R = lambda z: mp.re(mp.polylog(2, mp.mpc(z)))
Li3R = lambda z: mp.re(mp.polylog(3, mp.mpc(z)))
def Creal(a, b, dps=DPS):
    mp.mp.dps = dps
    pts = [mp.mpf(0), mp.mpf(1)]
    for c in (a,b):
        tc = 1/mp.mpf(c)
        if 0 < tc < 1: pts.append(tc)
    pts = sorted(pts)
    return mp.quad(lambda t: mp.log(abs(1-a*t))*mp.log(abs(1-b*t))/t, pts, method='tanh-sinh')
def f_value(xa, ya):
    xf, yf = pfloat(xa), pfloat(ya)
    if abs(xf) < 1 and abs(xf*yf) < 1 and abs(yf) <= 1:
        return f_prefix_series(xa, ya)
    a = mpf_of(ya); b = mpf_of(pmul(xa,ya))
    return Creal(a,b) + mp.log(abs(1-a))*(Li2R(b)-Li2R(b/a)) - Li3R(b)

task1 = pickle.load(open(HERE+'/dod5_task1.pkl','rb'))
TERMS = task1['TERMS']
vals = {}
for name, kind, key, xa, ya in TERMS:
    vals[(kind,key)] = f_value(xa, ya)
print("[setup] 36 values recomputed at %d dps, t=%.0fs"%(DPS, time.time()-t0))

# ---------- [1] check 36 closures ----------
cl = pickle.load(open(HERE+'/dod5_closures36.pkl','rb'))
worst1 = mp.mpf(0)
for key, clo in cl.items():
    dv = mp.mpf(0)
    for j in range(8):
        if clo['solT'][j]: dv += rat(clo['solT'][j])*LI3B[j]
    for j in range(19):
        if clo['solP'][j]: dv += rat(clo['solP'][j])*PBV[j]
    rem = vals[key] - dv
    rel = clo['rel']
    resid = abs(rel[0]*rem + sum(rel[i+1]*ivals10[i] for i in range(10)))
    worst1 = max(worst1, resid)
print("[1] 36 depth-2 constant identities, worst residual:", mp.nstr(worst1,3))

# ---------- [2] check 23 Li3 reductions ----------
red = pickle.load(open(HERE+'/dod5_all_li3_reductions.pkl','rb'))
worst2 = mp.mpf(0)
for z, r in red.items():
    dv = mp.mpf(0)
    for j in range(8):
        if r['solT'][j]: dv += rat(r['solT'][j])*LI3B[j]
    for j in range(19):
        if r['solP'][j]: dv += rat(r['solP'][j])*PBV[j]
    rem = Li3R(mpf_of(z)) - dv
    rel = r['rel']
    resid = abs(rel[0]*rem + sum(rel[i+1]*ivals10[i] for i in range(10)))
    worst2 = max(worst2, resid)
print("[2] 23 Li3-point reductions, worst residual:", mp.nstr(worst2,3))

# ---------- [3] R pieces ----------
l2s,l3s,lus = sp.symbols('l2 l3 lu')
Li2A = sp.Function('Li2A'); Li3A = sp.Function('Li3A')
vals_sym = {l2s: mp.log(2), l3s: mp.log(3), lus: mp.log(2-r3)}
worst3 = mp.mpf(0)
T3a = {}
for (i,j,k) in TRIPLES:
    xi = mpf_of(Xs[i]); xj = mpf_of(Xs[j]); xk = mpf_of(Xs[k])
    T3a[(i,j,k)] = mp.quad(lambda tt: mp.log(1-xj*tt)*mp.log(1-xk*tt)/(1-xi*tt), [0,1])
    Rq = T3a[(i,j,k)] - (vals[('SwG',(i,j,k))]-vals[('SwD',(i,j,k))])/xi
    Rs = build_R(i,j,k)
    e = Rs.subs(vals_sym)
    e = e.replace(Li2A, lambda z: mp.re(mp.polylog(2, mp.mpc(sp.N(z,360))))).replace(Li3A, lambda z: mp.re(mp.polylog(3, mp.mpc(sp.N(z,360)))))
    Rv = mp.mpf(sp.N(e, DPS))
    worst3 = max(worst3, abs(Rv-Rq))
print("[3] 12 R pieces sym vs (T3a quad - swaps/x), worst residual:", mp.nstr(worst3,3), " t=%.0fs"%(time.time()-t0))

# ---------- [4] Xi36 ----------
xc = pickle.load(open(HERE+'/dod5_xi36_closed.pkl','rb'))
ac = pickle.load(open(D2+'atom_coefs_pq.pkl','rb'))
xi36_val = mp.mpf(0)
for (kind, ij),(p,q) in ac.items():
    c = mp.mpf(p.numerator)/p.denominator + r3*mp.mpf(q.numerator)/q.denominator
    if kind=='Lam':
        xi36_val += c*vals[('Lam',ij)]
    else:
        i,j,k = ij
        xi36_val += c/mpf_of(Xs[i])*(vals[('SwG',ij)] - vals[('SwD',ij)])
xi36_closed = mp.mpf(0)
for j in range(8): xi36_closed += fval(xc['A3'][j])*LI3B[j]
for j in range(19): xi36_closed += fval(xc['AP'][j])*PBV[j]
for m in range(5): xi36_closed += fval(xc['AI'][m])*IBV[m]
d4 = abs(xi36_val - xi36_closed)
print("[4] Xi36 (36-term sum) vs Xi36 closed form, residual:", mp.nstr(d4,3))

# ---------- [5] full assembly vs P3_630 ----------
tl = pickle.load(open(HERE+'/dod5_tallies.pkl','rb'))
Tot3, TotP, TotI, singles_z3 = tl['Tot3'], tl['TotP'], tl['TotI'], tl['singles_z3']
rawc = pickle.load(open(HERE+'/dod5_rawli3_coeffs.pkl','rb'))
aco = pickle.load(open(HERE+'/dod5_atom_coefs_all.pkl','rb'))
Rc = pickle.load(open(HERE+'/dod5_R_content.pkl','rb'))
sd = pickle.load(open(D1+'S_atom_decomp.pkl','rb'))
vT3 = sum(fval(Tot3[j])*LI3B[j] for j in range(8))
vTP = sum(fval(TotP[j])*PBV[j] for j in range(19))
vTI = sum(fval(TotI[m])*IBV[m] for m in range(5))
# S_d1 = sum c*atom + singles ; atoms by quadrature/series/polylog
IDX=[1,2,5,6]
def fpq(pq):
    p,q = pq
    return mp.mpf(str(p)) + r3*mp.mpf(str(q))
subv = {'l2': mp.log(2), 'l3': mp.log(3), 'lu1': mp.log(2-r3), 'z2': mp.pi**2/6, 'z3': mp.zeta(3)}
l2,l3,lu1,z2,z3 = sp.symbols('l2 l3 lu1 z2 z3')
singles_val = mp.mpf(sp.N(sd['S_singles_part'].subs({l2:subv['l2'],l3:subv['l3'],lu1:subv['lu1'],z2:subv['z2'],z3:subv['z3']}), DPS))
S_d1 = singles_val
Xf = {i: mpf_of(Xs[i]) for i in IDX}
for k in aco:
    fam, kk = k
    if fam in ('Lam','T3a'): continue
    if fam=='T2':
        i,j = kk
        v = mp.quad(lambda tt: mp.log(1-Xf[j]*tt)/(1-Xf[i]*tt), [0,1])
    elif fam=='T3p':
        i,j = kk
        v = mp.quad(lambda tt: mp.log(1-Xf[j]*tt)**2/(1-Xf[i]*tt), [0,1])
    elif fam=='Ax':
        x = Xf[kk]; tot = mp.mpf(0); H = mp.mpf(0); xp = mp.mpf(1); n = 0
        while True:
            n += 1; H += 1/mp.mpf(n); xp *= x; tot += H*xp/n**2
            if n > 30 and abs(H*xp/n**2) < mp.mpf(10)**(-(DPS+15)): break
        v = tot
    elif fam=='L2x': v = Li2R(Xf[kk])
    elif fam=='L3x': v = Li3R(Xf[kk])
    cd = aco[k]
    cval = fpq(cd.get('1',(0,0)))
    for ls in ('l2','l3','lu1'):
        if ls in cd: cval += fpq(cd[ls])*subv[ls]
    S_d1 += cval*v
# raw Li3 value of S_d1 part
rawLi3_all = sum(fval(c)*Li3R(mpf_of(pt)) for pt,c in rawc.items())
vXi36A3 = sum(fval(xc['A3'][j])*LI3B[j] for j in range(8))
vR3 = sum(fval(tuple(psplit(c)))*Li3R(mpf_of(pt)) for pt,c in Rc['li3_content'].items())
rawLi3_Sd1 = rawLi3_all - vXi36A3 - vR3
Axz3 = sum(fpq(aco[('Ax',i)]['1']) for i in IDX)*mp.zeta(3)
singles_z3_val = fval(singles_z3)*mp.zeta(3)
Rprod = mp.mpf(0)
for s in ('li2log','logcube'):
    e = sp.sympify(Rc[s]).subs(vals_sym)
    e = e.replace(Li2A, lambda z: mp.re(mp.polylog(2, mp.mpc(sp.N(z,360))))).replace(Li3A, lambda z: mp.re(mp.polylog(3, mp.mpc(sp.N(z,360)))))
    Rprod += mp.mpf(sp.N(e, DPS))
Bval = (S_d1 - rawLi3_Sd1 - Axz3 - singles_z3_val) + Rprod
S_calc = vT3 + vTP + vTI + Bval
L = mp.sqrt(26-15*r3)
P3 = pickle.load(open(D1+'P3_630.pkl','rb'))['P3']
P3_calc = 72*L*S_calc
d5 = abs(P3_calc - P3)
print("[5] 72*L*S_calc vs P3_630:")
print("    P3_calc =", mp.nstr(P3_calc, 42))
print("    P3_630  =", mp.nstr(P3, 42))
print("    residual:", mp.nstr(d5,3), "  (%.0f digits)"%(-mp.log10(d5)))
ok = worst1 < mp.mpf('1e-290') and worst2 < mp.mpf('1e-290') and worst3 < mp.mpf('1e-290') and d4 < mp.mpf('1e-290') and d5 < mp.mpf('1e-290')
print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED", "  total time %.0fs"%(time.time()-t0))
'''
open('<project>/research/08_theory/DOD5_dodecagon_termwise/verify_dod5.py','w').write(script)
print("written, %d chars"%len(script))
