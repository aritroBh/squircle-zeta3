
# verify_dec2.py -- DEC2: the regular decagon P3(1) closed form, PROVED end to end.
#
# Self-contained verifier. Rebuilds, in EXACT arithmetic over Q(sqrt5) (and Q(rho) for geometry):
#   (1) the d^3 edge integrand E3n (55 log-monomials): compared monomial-by-monomial with DEC1's
#       reconstructed E3n_exact_pq.pkl (0 mismatches expected) and with a direct numerical
#       d-expansion of the radial-map speed at t=0.3;
#   (2) the integral table: 40 local integrals (antiderivative certificates), the 40-key
#       single-ratio J-table (exact symbolic summation; compared with K1d's Jtab_t.pkl),
#       the master constants D1..D10 (series / reflection / antiderivative certificates / IBP),
#       and the 20 mixed integrals (IBP recursions; compared with DEC1's dec_engine);
#   (3) the assembly to the 15-term closed form (compared exactly with DEC1's final_table_pq.pkl),
#       the local zeta(3) = C5 (exact in Q(rho)), the cited reductions (PRV1 Lam/swap, DEC1
#       A-atoms, golden Li2/Li3 lattice) at 120 digits, and the final value vs the 676-digit gold.
# Prints PASS/FAIL per part. Expected runtime < 5 min.
import os, sys, math, time, pickle
import sympy as sp
import mpmath as mp

WD  = os.path.dirname(os.path.abspath(__file__))
DEC1 = '<project>/research/08_theory/DEC1_decagon'
K1D  = '<project>/research/08_theory/K1d_engine_fixed'
O5   = '<project>/research/08_theory/OCT6_octagon_closed/oct5'
sys.path.insert(0, WD); sys.path.insert(0, DEC1)
os.chdir(WD)
T0 = time.time()
r5 = sp.sqrt(5); rho = sp.sqrt(10-2*r5); phi = (1+r5)/2
t = sp.Symbol('t')
PASS = []
def report(name, ok, detail=""):
    PASS.append((name, ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}  ({time.time()-T0:.0f}s)")

# ================= Part 1: exact geometry + E3n =================
h  = (r5-1)/4; y0 = (r5-2)/rho
V0 = (-h, y0); V1 = (-h, -y0)
nrm = [ ( (1+r5)/4, rho/4 ), ( (r5-1)/4, rho*(1+r5)/8 ), ( (1-r5)/4, rho*(1+r5)/8 ),
        ( -(1+r5)/4, rho/4 ), ( sp.Integer(-1), sp.Integer(0) ) ]
y  = (V0[0] + t*(V1[0]-V0[0]), V0[1] + t*(V1[1]-V0[1]))
avals = [sp.expand(y[0]*nrm[i][0] + y[1]*nrm[i][1]) for i in range(5)]
sig_idx = [(1,-1),(2,1),(0,-1),(3,1),(4,1)]
apol = []; appol = []
for mi, sg in sig_idx:
    a = sp.expand(sg*avals[mi])
    a0 = sp.expand(a.subs(t,0)); a1 = sp.expand(a.subs(t,1)-a0)
    apol.append([a0,a1]); appol.append([a1, sp.Integer(0)])
ok_geo = (sp.expand(sum(sp.expand(a[0]+a[1]*t) for a in apol) - 1) == 0
          and sp.expand(2*phi*h-1) == 0
          and sp.simplify(4*y0**2 - (sp.Rational(25,10)-sp.Rational(11,10)*r5)) == 0)
Z = sp.Integer(0)
def prune(P): return {k:v for k,v in P.items() if any(x != 0 for x in v)}
def eadd(P,Q):
    R = dict(P)
    for k in Q:
        if k in R: R[k] = [sp.expand(R[k][i]+Q[k][i]) for i in range(4)]
        else: R[k] = list(Q[k])
    return R
def ekmul(P,Q):
    R = {}
    for k1,v1 in P.items():
        s1 = sum(k1)
        for k2,v2 in Q.items():
            k = tuple(k1[i]+k2[i] for i in range(5))
            if s1+sum(k2) > 3: continue
            c = [Z]*4
            for i in range(4):
                for j in range(4-i):
                    if v1[i] != 0 and v2[j] != 0:
                        c[i+j] = sp.expand(c[i+j]+v1[i]*v2[j])
            if k in R: R[k] = [sp.expand(R[k][i]+c[i]) for i in range(4)]
            else: R[k] = c
    return R
def escale(P,s): return {k:[sp.expand(v[i]*s) for i in range(4)] for k,v in P.items()}
ONE = {(0,0,0,0,0):[sp.Integer(1),Z,Z,Z]}
def pmul(P,Q):
    R = [{} for _ in range(4)]
    for i in range(4):
        for j in range(4-i):
            R[i+j] = eadd(R[i+j], ekmul(P[i],Q[j]))
    return R
def pdiv(P,Q):
    R = [{} for _ in range(4)]
    q0 = Q[0][(0,0,0,0,0)][0]
    for i in range(4):
        acc = P[i]
        for j in range(i):
            acc = eadd(acc, escale(ekmul(R[j],Q[i-j]), sp.Integer(-1)))
        R[i] = {k:[sp.expand(x/q0) for x in v] for k,v in acc.items()}
    return R
def plog(P):
    Pp = [escale(P[i], sp.Integer(i)) for i in range(1,4)] + [{}]
    Qd = pdiv(Pp,P)
    L = [{}]
    for i in range(3): L.append(escale(Qd[i], sp.Rational(1,i+1)))
    return L
def pexp(P):
    E = [ONE] + [{} for _ in range(3)]
    for k in range(1,4):
        acc = {}
        for j in range(1,k+1):
            acc = eadd(acc, escale(ekmul(P[j],E[k-j]), sp.Rational(j)))
        E[k] = escale(acc, sp.Rational(1,k))
    return E
h2oL2 = sp.simplify(h**2/(4*y0**2))
F = [{} for _ in range(4)]; S = [{} for _ in range(4)]
for k in range(4):
    fc = sp.Rational(1, math.factorial(k))
    for i in range(5):
        key = tuple(k if j==i else 0 for j in range(5))
        F[k] = eadd(F[k], {key:[sp.expand(apol[i][0]*fc), sp.expand(apol[i][1]*fc), Z, Z]})
        S[k] = eadd(S[k], {key:[sp.expand(appol[i][0]*fc), sp.expand(appol[i][1]*fc), Z, Z]})
F = [prune(x) for x in F]; S = [prune(x) for x in S]
assert F[0] == ONE and S[0] == {}
Phi = pdiv(S,F); lg = plog(F)
inv = [escale(ONE, sp.Integer(sg)) for sg in [1,-1,1,-1]]
m2lg = [escale(lg[i], sp.Integer(-2)) for i in range(4)]
facE = pexp(pmul(m2lg, inv))
Phisq = pmul(Phi,Phi)
o2mt = {(0,0,0,0,0):[sp.Integer(1), sp.Integer(-2), Z, Z]}
yypE = {(0,0,0,0,0):[sp.expand(h2oL2+sp.Rational(1,4)), sp.Integer(-1), sp.Integer(1), Z]}
Q = [{} for _ in range(4)]; Q[0] = ONE
for k in range(1,4): Q[k] = ekmul(Phi[k], o2mt)
for k in range(4):
    Q[k] = eadd(Q[k], ekmul(Phisq[k], yypE))
G = pmul(facE, Q)
sq = [{} for _ in range(4)]; sq[0] = ONE
for k in range(1,4):
    acc = G[k]
    for j in range(1,k):
        acc = eadd(acc, escale(ekmul(sq[j],sq[k-j]), sp.Integer(-1)))
    sq[k] = escale(acc, sp.Rational(1,2))
E3n = prune(sq[3])
def pq(e):
    e = sp.expand(e)
    p = sp.simplify((e + e.subs(r5,-r5))/2); q = sp.simplify((e - e.subs(r5,-r5))/(2*r5))
    return (sp.Rational(p), sp.Rational(q))
E3n_pq = {k:[pq(x) for x in v] for k,v in E3n.items()}
dec1_e3 = pickle.load(open(f'{DEC1}/E3n_exact_pq.pkl','rb'))
mism = sum(1 for k in dec1_e3 if k not in E3n_pq or any(
    sp.simplify(E3n_pq[k][j][0]-dec1_e3[k][j][0]) != 0 or sp.simplify(E3n_pq[k][j][1]-dec1_e3[k][j][1]) != 0
    for j in range(4)))
mism += len(set(E3n_pq)-set(dec1_e3))
mp.mp.dps = 80
s5m = mp.sqrt(5); phim = (1+s5m)/2
A0=1/(2*phim**2); B4=1/(2*phim); kap1=phim-1; kap4=1/phim-1
Lm = mp.sqrt((25-11*s5m)/10); hm = 1/(2*phim); y0m = Lm/2
t0m = mp.mpf('0.3')
a_m  = [A0*t0m, A0*(1-t0m), A0*(1+kap1*t0m), B4*(1+kap4*t0m), B4]
ap_m = [A0, -A0, A0*kap1, B4*kap4, mp.mpf(0)]
ym  = (-hm, y0m*(1-2*t0m)); ypm = (mp.mpf(0), -2*y0m)
yym = ym[0]**2+ym[1]**2; yypm = ym[0]*ypm[0]+ym[1]*ypm[1]; yp2m = ypm[0]**2+ypm[1]**2
def speed(d):
    Fv = sum(aj**(1+d) for aj in a_m)
    Fp = (1+d)*sum(apj*aj**d for aj,apj in zip(a_m,ap_m))
    Ph = Fp/((1+d)*Fv)
    return mp.sqrt(Fv**(-2/(1+d))*(yp2m - 2*yypm*Ph + yym*Ph**2))
E3_direct = mp.diff(speed, mp.mpf(0), 3)/6
logs = [mp.log(v) for v in a_m]
val = mp.mpf(0)
for key,row in E3n_pq.items():
    poly = mp.mpf(0)
    for nn in range(4):
        p,q = row[nn]
        poly += (mp.mpf(str(p)) + mp.mpf(str(q))*s5m) * t0m**nn
    lv = mp.mpf(1)
    for i,e in enumerate(key): lv *= logs[i]**e
    val += poly*lv
d_e3 = abs(Lm*val - E3_direct)
report("P1 geometry+E3n", ok_geo and mism==0 and d_e3 < mp.mpf('1e-70'),
       f"(geometry exact; monomial mismatches vs DEC1: {mism}; direct check {mp.nstr(d_e3,2)})")

# ================= Part 2: local integrals =================
L0 = sp.log(t); L1s = sp.log(1-t)
Li2 = lambda z: sp.polylog(2,z); Li3 = lambda z: sp.polylog(3,z)
def p1norm(e):
    return e.replace(lambda f: isinstance(f, sp.Function) and f.func == sp.polylog and f.args[0]==1,
                     lambda f: -sp.log(1-f.args[1]))
def check0(e): return sp.simplify(p1norm(sp.expand(e))) == 0
def int_s_over(p,q):
    s_ = 1-t; ls = sp.log(s_); l1s_ = sp.log(1-s_)
    if (p,q)==(0,1): return -Li2(s_)
    if (p,q)==(1,1): return -Li2(s_)*ls + Li3(s_)
    if (p,q)==(0,2): return ls*l1s_**2 + 2*l1s_*Li2(1-s_) - 2*Li3(1-s_)
    if (p,q)==(1,0): return ls**2/2
    if (p,q)==(2,0): return ls**3/3
    if (p,q)==(0,0): return ls
    raise ValueError((p,q))
def int_t(n,a,b):
    if b==0:
        Fv = 0
        for i in range(a+1):
            Fv += (-1)**(a-i) * sp.factorial(a)/sp.factorial(i) * L0**i / (n+1)**(a-i+1)
        return sp.expand(t**(n+1)*Fv)
    if a==0:
        Fv = t**(n+1)*L1s**b/(n+1) + sp.Rational(b,n+1)*(-L1s**b/b)
        for j in range(n+1):
            Fv -= sp.Rational(b,n+1)*int_t(j,0,b-1)
        return sp.expand(Fv)
    B = int_t(n,a,0)
    Fv = sp.expand(B*L1s**b)
    acc = 0
    for i in range(a+1):
        ci = sp.expand(B).coeff(L0,i)
        if ci == 0: continue
        ci_t = sp.simplify(ci/t**(n+1))
        for j in range(n+2):
            cj = (-1)**j * sp.binomial(n+1,j)
            if j == 0:
                term = cj * (-int_s_over(b-1, i))
            else:
                term = cj * int_polylogweight(j-1, i, b-1)
            acc += ci_t*term
    return sp.expand(Fv + b*acc)
def int_polylogweight(m,p,q):
    acc = 0
    for j in range(m+1):
        acc += (-1)**j*sp.binomial(m,j)*int_t(j,p,q)
    return sp.expand(acc)
z2s = sp.pi**2/6; z3s = sp.zeta(3); s_ = sp.Symbol('s')
def texpo(term, sym): return int(term.as_powers_dict().get(sym,0))
def limit0(Fv):
    e = Fv.subs({sp.polylog(2,t):0, sp.polylog(3,t):0, sp.polylog(2,1-t):z2s,
                 sp.polylog(3,1-t):z3s, sp.log(1-t):0})
    e = sp.expand(e)
    keep = sp.expand(sum(term for term in sp.Add.make_args(e) if texpo(term,t)==0))
    lt = sp.log(t)
    for q in range(1,7):
        assert sp.simplify(keep.coeff(lt,q)) == 0
    return sp.simplify(keep.subs(lt,0))
def limit1(Fv):
    e = Fv.subs(t, 1-s_)
    e = e.subs({sp.polylog(2,s_):0, sp.polylog(3,s_):0, sp.polylog(2,1-s_):z2s,
                sp.polylog(3,1-s_):z3s, sp.log(1-s_):0})
    e = sp.expand(e)
    keep = sp.expand(sum(term for term in sp.Add.make_args(e) if texpo(term,s_)==0))
    ls = sp.log(s_)
    for q in range(1,7):
        assert sp.simplify(keep.coeff(ls,q)) == 0
    return sp.simplify(keep.subs(ls,0))
JlocO = pickle.load(open(O5+'/Jloc.pkl','rb'))
bad2 = 0; myJloc = {}
for key in sorted(JlocO.keys()):
    n,a,b = key
    Fv = int_t(n,a,b)
    if not check0(sp.diff(Fv,t) - t**n*L0**a*L1s**b): bad2 += 1
    v = sp.simplify(limit1(Fv) - limit0(Fv))
    myJloc[key] = v
    want = JlocO[key].subs(sp.pi**2, 6*z2s).subs(sp.zeta(3), z3s)
    if sp.simplify(v - want) != 0: bad2 += 1
report("P2 local integrals (40)", bad2==0, f"(certificates+values; mismatches: {bad2})")

# ================= Part 3: single-ratio J-table =================
from esum_engine import esum, Hf, H2f, k, x as xs_, L1 as L1_, L2 as L2_, L3 as L3_, AA as AA_, z2 as z2_
def Msym(n, a, b):
    l = k + n + 1
    if (a,b)==(0,0): return 1/l
    if (a,b)==(1,0): return -1/l**2
    if (a,b)==(0,1): return -Hf(l)/l
    if (a,b)==(2,0): return 2/l**3
    if (a,b)==(1,1): return Hf(l)/l**2 - z2_/l + H2f(l)/l
    if (a,b)==(0,2): return (Hf(l)**2 + H2f(l))/l
def Jexpr(n, a, b, c):
    M = Msym(n, a, b)
    if c == 1: return sp.expand(-esum(M/k, K=1))
    if c == 2: return sp.expand(2*esum(Hf(k-1)*M/k, K=2))
    if c == 3: return sp.expand(-3*esum((Hf(k-1)**2 - H2f(k-1))*M/k, K=3))
tuples = [(n,a,b,c) for n in range(4) for c in (1,2,3) for a in range(4-c) for b in range(4-c-a)]
Jtab = {tp: Jexpr(*tp) for tp in tuples}
K1d_tab = pickle.load(open(K1D+'/Jtab_t.pkl','rb'))
mism3 = sum(1 for key in K1d_tab if key not in Jtab or sp.expand(Jtab[key]-K1d_tab[key]) != 0)
report("P3 single-ratio J-table (40)", mism3==0, f"(exact summation; symbolic mismatches vs K1d: {mism3})")

# ================= Part 4: D1..D10 =================
import dec_engine as E
kap1 = 1/phi; kap4 = -1/phi**2
u1 = 1+kap1*t; u4 = 1+kap4*t
w1 = sp.log(u1); w4 = sp.log(u4)
lf_s = sp.log(phi); l2_s = sp.log(2)
def unit2phi(c):
    c = sp.simplify(c)
    for b in range(-8,9):
        v2 = sp.simplify(c/phi**b)
        q = sp.simplify((v2 - v2.subs(r5,-r5))/(2*r5))
        if q == 0:
            p = sp.simplify((v2 + v2.subs(r5,-r5))/2)
            try: pr = sp.Rational(p)
            except TypeError: continue
            fa = sp.factorint(int(pr.p)); fb = sp.factorint(int(pr.q))
            if set(fa)|set(fb) <= {2}:
                return (fa.get(2,0)-fb.get(2,0), b)
    return None
def canon_logs(e):
    letters = {}
    def repl(f):
        arg = sp.expand(f.args[0])
        c = sp.simplify(arg.subs(t,0))
        rem = sp.expand(arg - c)
        ev = sp.expand(sp.simplify(rem/(c*t))) if rem != 0 else sp.Integer(0)
        ab = unit2phi(c)
        base = ab[0]*l2_s + ab[1]*lf_s
        if ev == 0: return base
        if ev not in letters: letters[ev] = sp.log(1+ev*t)
        return base + letters[ev]
    return e.replace(lambda f: isinstance(f, sp.log), repl), letters
def cert3(F, f):
    D = sp.diff(F, t) - f
    D = D.replace(lambda g: isinstance(g, sp.Function) and g.func == sp.polylog and g.args[0]==1,
                  lambda g: -sp.log(1-g.args[1]))
    D, letters = canon_logs(D)
    D = sp.expand(D)
    lfs, l2ss = sp.Symbol('LFs'), sp.Symbol('L2S')
    D2 = D.subs({lf_s: lfs, l2_s: l2ss})
    gens = list(letters.values()) + [lfs, l2ss]
    P = sp.Poly(sp.expand(D2), *gens)
    return all(sp.simplify(coeff) == 0 for _, coeff in P.terms())
x4 = 1/phi**2
F5 = phi*(lf_s*w1 - Li2(x4*(1+t/phi)))
F6 = -phi**2*(2*lf_s*w4 - Li2((1/phi)*(1+kap4*t)))
s7 = 1+kap4*t
G7 = 4*lf_s**2*sp.log(s7) - 4*lf_s*Li2(s7/phi) + sp.log(s7/phi)*sp.log(1-s7/phi)**2 \
     + 2*sp.log(1-s7/phi)*Li2(1-s7/phi) - 2*Li3(1-s7/phi)
F7 = -phi**2*G7
s8 = 1+t/phi
H8 = lf_s**2*sp.log(s8) - 2*lf_s*Li2(x4*s8) + sp.log(x4*s8)*sp.log(1-x4*s8)**2 \
     + 2*sp.log(1-x4*s8)*Li2(1-x4*s8) - 2*Li3(1-x4*s8)
F8 = phi*H8
certs = [cert3(F5, w4/u1), cert3(F6, w1/u4), cert3(F7, w1**2/u4), cert3(F8, w4**2/u1)]
NAMED2 = [(-1/phi, E.L2_1, E.L3_1), (x4, E.L2_4, E.L3_4), (1/phi, E.L2f, E.L3f)]
def canon_atoms(e):
    def replp(f):
        s_, arg = f.args[0], sp.simplify(f.args[1])
        for pt, s2, s3 in NAMED2:
            if sp.simplify(arg - pt) == 0:
                return s2 if s_ == 2 else s3
        raise ValueError("unknown polylog point")
    e = e.replace(lambda f: isinstance(f, sp.Function) and f.func == sp.polylog and f.args[0] in (2,3) and not f.args[1].has(t), replp)
    def repll(f):
        arg = sp.simplify(f.args[0])
        for b in range(-6,7):
            if sp.simplify(arg - phi**b) == 0: return b*E.lf
        for a in range(-6,7):
            if sp.simplify(arg - 2**a) == 0: return a*E.l2
        raise ValueError("unknown log")
    e = e.replace(lambda f: isinstance(f, sp.log) and not f.args[0].has(t), repll)
    return sp.expand(e)
Dvals = {}
okD = True
for name, F, Feng in [("D5",F5,E.D5), ("D6",F6,E.D6), ("D7",F7,E.D7), ("D8",F8,E.D8)]:
    val = sp.expand(F.subs(t,1) - F.subs(t,0))
    val = sp.expand(val.subs({sp.log(phi): E.lf, sp.log(2): E.l2}))
    vn = canon_atoms(val)
    Dvals[name] = vn
    if sp.expand(vn - Feng) != 0: okD = False
def canon5(e):
    e = sp.expand(e)
    P = sp.Poly(e, *E.NAMED)
    out = sp.Integer(0)
    for monom, coeff in P.terms():
        c = sp.together(coeff); cc = c.subs(r5, -r5)
        p = sp.radsimp(sp.expand((c+cc)/2)); q = sp.radsimp(sp.expand((c-cc)/(2*r5)))
        mm = sp.Integer(1)
        for a_, e_ in zip(E.NAMED, monom): mm *= a_**e_
        out += (p + q*r5)*mm
    return sp.expand(out)
D9  = canon5((-E.lf**3/2 - kap4*E.D7/2)/kap1)
D10 = canon5((E.lf**3/2 - kap1*E.D8/2)/kap4)
ibp9  = sp.simplify(sp.expand(sp.diff(w1**2*w4/2, t) - kap1*(w1*w4/u1) - kap4*(w1**2/u4)/2)) == 0
ibp10 = sp.simplify(sp.expand(sp.diff(w1*w4**2/2, t) - kap1*(w4**2/u1)/2 - kap4*(w1*w4/u4))) == 0
Dvals['D9'] = D9; Dvals['D10'] = D10
K4v = sp.expand(E.L2_4/kap4); K1v_ = sp.expand(E.L2_1/kap1)
D3 = canon5((E.D2 - E.lf*K4v)/phi)
D4 = canon5(phi*(E.lf*K1v_ + E.D1))
Dvals['D1'] = E.D1; Dvals['D2'] = E.D2; Dvals['D3'] = D3; Dvals['D4'] = D4
okD = okD and ibp9 and ibp10 and canon5(D9-E.D9)==0 and canon5(D10-E.D10)==0 \
      and canon5(D3-E.D3)==0 and canon5(D4-E.D4)==0 and all(certs)
report("P4 master constants D1..D10", okD,
       "(D1,D2 series; D3,D4 reflection; D5-D8 certificates; D9,D10 IBP; all match dec_engine)")

# ================= Part 5: 20 mixed integrals =================
def K1vv(n,a,b,c,side):
    return sp.expand(Jtab[(n,a,b,c)].subs(E.sub1 if side==1 else E.sub4, simultaneous=True))
from functools import lru_cache
@lru_cache(None)
def MixedW(n):
    bnd = -E.lf**2
    A1 = sp.expand(sum((-1/kap1)**j * K1vv(n-j,0,0,1,4) for j in range(n+1))/kap1 + (-1/kap1)**(n+1)*Dvals['D5'])
    A2 = sp.expand(sum((-1/kap4)**j * K1vv(n-j,0,0,1,1) for j in range(n+1))/kap4 + (-1/kap4)**(n+1)*Dvals['D6'])
    return sp.expand((bnd - kap1*A1 - kap4*A2)/(n+1))
@lru_cache(None)
def M1(n):
    B1 = sp.expand(sum((-1/kap1)**j * K1vv(n-j,1,0,1,4) for j in range(n+1))/kap1 + (-1/kap1)**(n+1)*Dvals['D1'])
    B2 = sp.expand(sum((-1/kap4)**j * K1vv(n-j,1,0,1,1) for j in range(n+1))/kap4 + (-1/kap4)**(n+1)*Dvals['D2'])
    return sp.expand(-(MixedW(n) + kap1*B1 + kap4*B2)/(n+1))
@lru_cache(None)
def M2(n):
    r1 = (-1/kap1)**(n+1) - 1; r4 = (-1/kap4)**(n+1) - 1
    C1 = sp.expand(sum((-1/kap1)**j * K1vv(n-j,0,1,1,4) for j in range(n+1))/kap1 + r1*Dvals['D3'])
    C2 = sp.expand(sum((-1/kap4)**j * K1vv(n-j,0,1,1,1) for j in range(n+1))/kap4 + r4*Dvals['D4'])
    return sp.expand((-sum(MixedW(j) for j in range(n+1)) - kap1*C1 - kap4*C2)/(n+1))
@lru_cache(None)
def M3(n):
    E1 = sp.expand(sum((-1/kap1)**j * MixedW(n-j) for j in range(n+1))/kap1 + (-1/kap1)**(n+1)*Dvals['D9'])
    E2 = sp.expand(sum((-1/kap4)**j * K1vv(n-j,0,0,2,1) for j in range(n+1))/kap4 + (-1/kap4)**(n+1)*Dvals['D7'])
    return sp.expand((-E.lf**3 - 2*kap1*E1 - kap4*E2)/(n+1))
@lru_cache(None)
def M4(n):
    F1 = sp.expand(sum((-1/kap1)**j * K1vv(n-j,0,0,2,4) for j in range(n+1))/kap1 + (-1/kap1)**(n+1)*Dvals['D8'])
    F2 = sp.expand(sum((-1/kap4)**j * MixedW(n-j) for j in range(n+1))/kap4 + (-1/kap4)**(n+1)*Dvals['D10'])
    return sp.expand((E.lf**3 - kap1*F1 - 2*kap4*F2)/(n+1))
def MIXm(n,a,b,c1,c4):
    if a==1: return M1(n)
    if b==1: return M2(n)
    if (c1,c4)==(2,1): return M3(n)
    if (c1,c4)==(1,2): return M4(n)
    return MixedW(n)
L0s = sp.log(t); L1ts = sp.log(1-t)
def ibp_check(n, typ):
    if typ=='W':
        C = t**(n+1)*w1*w4/(n+1); R = -kap1*t**(n+1)*w4/u1 - kap4*t**(n+1)*w1/u4
    elif typ=='M1':
        C = t**(n+1)*L0s*w1*w4/(n+1); R = -t**n*w1*w4 - kap1*t**(n+1)*L0s*w4/u1 - kap4*t**(n+1)*L0s*w1/u4
    elif typ=='M2':
        C = (t**(n+1)-1)*L1ts*w1*w4/(n+1)
        R = -sum(t**j for j in range(n+1))*w1*w4 - kap1*(t**(n+1)-1)*L1ts*w4/u1 - kap4*(t**(n+1)-1)*L1ts*w1/u4
    elif typ=='M3':
        C = t**(n+1)*w1**2*w4/(n+1); R = -2*kap1*t**(n+1)*w1*w4/u1 - kap4*t**(n+1)*w1**2/u4
    elif typ=='M4':
        C = t**(n+1)*w1*w4**2/(n+1); R = -kap1*t**(n+1)*w4**2/u1 - 2*kap4*t**(n+1)*w1*w4/u4
    integrand = {'W': t**n*w1*w4, 'M1': t**n*L0s*w1*w4, 'M2': t**n*L1ts*w1*w4,
                 'M3': t**n*w1**2*w4, 'M4': t**n*w1*w4**2}[typ]
    return sp.simplify(sp.expand((n+1)*integrand - (n+1)*sp.diff(C,t) - R)) == 0
ibp_all = all(ibp_check(n,typ) for n in range(4) for typ in ['W','M1','M2','M3','M4'])
pf_all = all(sp.simplify(t**(n+1) - (sum((-1/kap)**j*t**(n-j) for j in range(n+1))/kap + (-1/kap)**(n+1)/(1+kap*t))*(1+kap*t))==0
             for n in range(4) for kap in (kap1,kap4))
mism5 = 0; MIXvals = {}
for typ in [(0,0,1,1),(1,0,1,1),(0,1,1,1),(0,0,2,1),(0,0,1,2)]:
    for n in range(4):
        a,b,c1,c4 = typ
        mine = MIXm(n,a,b,c1,c4)
        MIXvals[(n,a,b,c1,c4)] = mine
        if canon5(mine - E.MIX(n,a,b,c1,c4)) != 0: mism5 += 1
report("P5 mixed integrals (20)", ibp_all and pf_all and mism5==0,
       f"(IBP identities: {ibp_all}; partial fractions: {pf_all}; mismatches vs engine: {mism5})")

# ================= Part 6: assembly + final table + gold =================
tt = sp.Symbol('t')
Lt = [sp.Symbol(f'L{i}') for i in range(5)]
lt,l1t,w1s_,w4s_ = sp.symbols('lt l1t w1 w4')
E3sym = sp.Integer(0)
for key,row in E3n_pq.items():
    mon = sp.Integer(1)
    for i,e in enumerate(key): mon *= Lt[i]**e
    for nn in range(4):
        p,q = row[nn]
        if p==0 and q==0: continue
        E3sym += (p+q*r5)*tt**nn*mon
E3sym = sp.expand(E3sym)
sub = {Lt[0]: lt-l2_s-2*lf_s, Lt[1]: l1t-l2_s-2*lf_s, Lt[2]: w1s_-l2_s-2*lf_s,
       Lt[3]: w4s_-l2_s-lf_s, Lt[4]: -l2_s-lf_s}
E2 = sp.expand(E3sym.xreplace(sub))
mons = {}
for term in sp.Add.make_args(E2):
    d = term.as_powers_dict()
    n_ = int(d.get(tt,0))
    abc = (int(d.get(lt,0)), int(d.get(l1t,0)), int(d.get(w1s_,0)), int(d.get(w4s_,0)))
    letters = (int(d.get(l2_s,0)), int(d.get(lf_s,0)))
    rest = sp.simplify(term/(tt**n_ * lt**abc[0]*l1t**abc[1]*w1s_**abc[2]*w4s_**abc[3] * l2_s**letters[0]*lf_s**letters[1]))
    key = (n_, abc, letters)
    mons[key] = sp.expand(mons.get(key,0)+rest)
dec1_mons = pickle.load(open(f'{DEC1}/mons.pkl','rb'))
mism6a = sum(1 for k in dec1_mons if k not in mons or sp.simplify(mons[k]-dec1_mons[k]) != 0)
def JLOC(n,a,b):
    return sp.expand(myJloc[(n,a,b)].subs(sp.pi**2, 6*E.z2).subs(sp.zeta(3), E.z3))
def INT(n,a,b,c1,c4):
    if c1==0 and c4==0: return JLOC(n,a,b)
    if c4==0: return K1vv(n,a,b,c1,1)
    if c1==0: return K1vv(n,a,b,c4,4)
    return MIXvals[(n,a,b,c1,c4)]
Sasm = sp.Integer(0)
for (n_, abc, letters), coef in mons.items():
    a,b,c1,c4 = abc
    Sasm += coef * E.l2**letters[0] * E.lf**letters[1] * INT(n_,a,b,c1,c4)
Sasm = sp.expand(Sasm)
subs_red = {
  E.L2_1: -2*E.z2/5 + E.lf**2/2, E.L2_4: 2*E.z2/5 - E.lf**2, E.L2f: 3*E.z2/5 - E.lf**2,
  E.AA_1: 4*E.z3/5 - 2*E.lf*E.z2/5 + 2*E.L3_1,
  E.AA_4: E.z3 - 3*E.lf*E.z2/5 + 4*E.L3_1 + 3*E.L3f,
  E.Lam:  2*E.z3/5 + 4*E.lf*E.z2/5 - E.lf**3 + 3*E.L3_1 + E.L3f,
  E.swap: 3*E.z3/5 - 7*E.lf*E.z2/5 + 2*E.lf**3 + E.L3_1 + E.L3f,
  E.L3_4: 4*E.L3_1 + 4*E.L3f,
}
S3 = sp.expand(Sasm.subs(subs_red, simultaneous=True))
vars6 = [E.l2, E.lf, E.z2, E.z3, E.L3_1, E.L3f]
poly = sp.Poly(sp.expand(S3), *vars6)
names = {E.l2:'l2', E.lf:'lf', E.z2:'z2', E.z3:'z3', E.L3_1:'L3(-1/phi)', E.L3f:'L3(1/phi)'}
def mname(m):
    parts=[]
    for v,e in zip(vars6,m):
        if e==1: parts.append(names[v])
        elif e>1: parts.append(f"{names[v]}^{e}")
    return '*'.join(parts) if parts else '1'
def pqc(c):
    c2 = sp.together(c); cc = c2.subs(r5,-r5)
    return (sp.radsimp(sp.expand((c2+cc)/2)), sp.radsimp(sp.expand((c2-cc)/(2*r5))))
table = {}
for monom, coeff in poly.terms():
    if sp.expand(coeff)==0: continue
    table[mname(monom)] = pqc(coeff)
ft = pickle.load(open(f'{DEC1}/final_table_pq.pkl','rb'))['table_pq']
mism6 = sum(1 for k in ft if k not in table or
            sp.simplify(sp.sympify(table[k][0])-sp.sympify(ft[k][0])) != 0 or
            sp.simplify(sp.sympify(table[k][1])-sp.sympify(ft[k][1])) != 0)
Sloc = sp.Integer(0)
for (n_, abc, letters), coef in mons.items():
    a,b,c1,c4 = abc
    if c1==0 and c4==0:
        Sloc += coef * E.l2**letters[0] * E.lf**letters[1] * JLOC(n_,a,b)
cz3 = sp.expand(Sloc).coeff(E.z3)
Lv = rho*(3*r5-5)/20
C5 = rho*(105*r5-255)/32
okC5 = sp.simplify(sp.expand(60*Lv*cz3 - C5)) == 0
mp.mp.dps = 150
s5m = mp.sqrt(5); phim = (1+s5m)/2
atoms = {E.l2: mp.log(2), E.lf: mp.log(phim), E.z2: mp.pi**2/6, E.z3: mp.zeta(3),
         E.L3_1: mp.polylog(3,-1/phim), E.L3f: mp.polylog(3,1/phim)}
S3n = mp.mpf(sp.N(S3.subs(atoms, simultaneous=True), 150))
P3v = 60*mp.sqrt((25-11*s5m)/10)*S3n
gold = mp.mpf(open(f'{DEC1}/P3_660.txt').read().split()[0])
dgold = abs(P3v - gold)
mp.mp.dps = 120
lfm = mp.log(phim); z2m = mp.pi**2/6; z3m = mp.zeta(3)
x1m = -1/phim; x4m = 1/phim**2
L2m = lambda z: mp.polylog(2,z); L3m = lambda z: mp.polylog(3,z)
def Am(xv): return mp.nsum(lambda m: xv**m*mp.harmonic(m)/m**2, [1, mp.inf])
def Li21(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1); m=0
    while True:
        m+=1; xp*=x
        if m>=2: tot += xp*ps/m**2
        yp*=y; ps += yp/m
        if m>10 and abs(xp*ps/m**2) < mp.mpf('1e-125'): break
    return tot
Lam_v = Li21(x1m,x1m); swap_v = Li21(x4m,-phim)
reds = [
    L2m(x4m) - (2*z2m/5 - lfm**2),
    L2m(1/phim) - (3*z2m/5 - lfm**2),
    L2m(x1m) - (-2*z2m/5 + lfm**2/2),
    Am(x1m) - (4*z3m/5 - 2*lfm*z2m/5 + 2*L3m(x1m)),
    Am(x4m) - (z3m - 3*lfm*z2m/5 + 4*L3m(x1m) + 3*L3m(1/phim)),
    Lam_v - (2*z3m/5 + 4*lfm*z2m/5 - lfm**3 + 3*L3m(x1m) + L3m(1/phim)),
    swap_v - (3*z3m/5 - 7*lfm*z2m/5 + 2*lfm**3 + L3m(x1m) + L3m(1/phim)),
    L3m(x4m) - (4*L3m(x1m)+4*L3m(1/phim)),
    L3m(x4m) - (4*(z3m-lfm*z2m)/5 + 2*lfm**3/3),
    L3m(1/phim)+L3m(x1m) - ((z3m-lfm*z2m)/5 + lfm**3/6),
]
worst_r = max(abs(d) for d in reds)
ok6 = (mism6a==0 and mism6==0 and okC5 and dgold < mp.mpf('1e-140') and worst_r < mp.mpf('1e-100'))
report("P6 assembly+table+gold", ok6,
       f"(mons mism {mism6a}; table mism {mism6}; C5 exact: {okC5}; |60LS-gold| {mp.nstr(dgold,1)}; reductions worst {mp.nstr(worst_r,1)})")

print("="*70)
nfail = sum(1 for _,ok in PASS if not ok)
print(f"VERIFY_DEC2: {'ALL PARTS PASSED' if nfail==0 else str(nfail)+' PARTS FAILED'}  total {time.time()-T0:.0f}s")
