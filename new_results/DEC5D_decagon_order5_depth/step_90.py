
verifier = '''"""verify_dec5d.py — DEC5D depth-test verifier (< 5 min).

Reproduces, from the read-only inputs and the work-dir artifacts:
  (1) P^{(5)}(1) by the polar-series engine (fresh tanh-sinh, ~dps 95): P(1) = 10 L,
      P'''(1) = 60 L S (paper/DEC1 15-term Q(sqrt5) form), P^{(5)} two-level self-agreement;
  (2) the exact order-five edge integrand: edge_params via symeng, the substituted E5/L
      (work-dir artifact E5_dec_subs.pkl, produced from read-only E5_edge_law) checked against
      symeng.E5_sym at 5 points, and the assembly of P^{(5)}/L and its weight-5 part from
      756 freshly integrated monomials (matches the engine);
  (3) letters {0,1,-phi,phi^2}: all letter differences factor as +-phi^k, constants give
      log 2 and log phi, so V = <log2, log phi>, dim V = 2;
  (4) the exact G-rewrite of the 756 integrals (1024 words) checked numerically against
      the fresh integral table;
  (5) controls: Li_5 at 14 S-units (incl. sqrt5 ones, d=3 machinery) pair to 0; 300 random
      shuffle products pair to 0; M4c W5 -> (84, -38 sqrt2); OCT5D2's octagon reproduces
      (-1485+1035 sqrt2, 900-540 sqrt2) with THIS general code;
  (6) the decagon obstruction: d=3 (24 pairings), d=2 M4c pairings, and the (4,1) component
      all VANISH on every block and on the total; and every block + total lies EXACTLY in the
      classical span (shuffle products + Li_5), verified with exact zero residual;
  (7) the exact local-block decomposition (work-dir artifact local_block_exact.pkl, built by
      250-digit PSLQ reconstruction of the 126 local integrals) checked numerically, and its
      zeta(5), zeta(2)zeta(3) coefficients equal the paper's prop:order5 values for m = 5.

Read-only inputs: OCT5D (parsed_mons.pkl), DEC1 (final_table_pq.pkl), M4c (symbol.json),
E5c_nonlocal_z5/symeng.py + E5_edge_law pkl, paper scripts/qsqrt2.py.
Work-dir artifacts: E5_dec_subs.pkl, parsed_mons_dec.pkl, local_block_exact.pkl, dec_engine.py,
dec5d_lib.py, gengine_dec.py (+ gtab cache). Exits nonzero on any failure.
"""
import os, sys, time, pickle, json, math, random, itertools
import mpmath as mp
from fractions import Fraction as F

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "<project>/research/08_theory/"
for p in (HERE, BASE+"E5c_nonlocal_z5", "<project>/paper/squircle_note/scripts"):
    if p not in sys.path: sys.path.insert(0, p)
failures = []
def check(name, cond, detail=""):
    print(("PASS  " if cond else "FAIL  ") + name + ("   " + detail if detail else ""))
    if not cond: failures.append(name)

# ---------- (1) value ----------
import dec_engine
res = dec_engine.run(95, os.path.join(HERE, "verify_p5.pkl"))
mp.mp.dps = 120
r5 = mp.sqrt(5); phi = (1+r5)/2
L = mp.sqrt((25-11*r5)/10)
check("P(1) = 10 L (perimeter)", abs(res[0]-10*L) < mp.mpf('1e-80'), mp.nstr(abs(res[0]-10*L),3))
ft = pickle.load(open(BASE+"DEC1_decagon/final_table_pq.pkl","rb"))['table_pq']
l2, lf, z2, z3 = mp.log(2), mp.log(phi), mp.pi**2/6, mp.zeta(3)
vals = {'l2^3': l2**3, 'l2^2*lf': l2**2*lf, 'l2^2': l2**2, 'l2*lf^2': l2*lf**2,
 'l2*lf': l2*lf, 'l2*z2': l2*z2, 'l2': l2, 'lf^3': lf**3, 'lf^2': lf**2,
 'lf*z2': lf*z2, 'lf': lf, 'z2': z2, 'z3': z3,
 'L3(-1/phi)': mp.polylog(3,-1/phi), 'L3(1/phi)': mp.polylog(3,1/phi)}
S15 = sum((mp.mpf(str(p))+mp.mpf(str(q))*r5)*vals[k] for k,(p,q) in ft.items())
check("P'''(1) = 60 L S (15-term closed form)", abs(res[3]-60*L*S15) < mp.mpf('1e-80'),
      mp.nstr(abs(res[3]-60*L*S15),3))
print("   P^(5)(1) =", mp.nstr(res[5], 60))
P5L = res[5]/L
print("   P^(5)(1)/L =", mp.nstr(P5L, 60))

# ---------- (2) exact structure ----------
import symeng as SE
mp.mp.dps = 50
h = 1/(2*phi); R = 1/(2*phi**2*mp.sin(mp.pi/5)); Ry = R/(2*phi)
V0 = [-h, Ry]; V1 = [-h, -Ry]
M = [[mp.cos(j*mp.pi/5), mp.sin(j*mp.pi/5)] for j in range(5)]
P = SE.edge_params(M, V0, V1)
ok = (P['j0']==2 and P['j1']==3 and P['nonloc']==[0,1,4]
      and abs(P['alpha']-1/(2*phi**2))<mp.mpf('1e-40') and abs(P['beta']-1/(2*phi**2))<mp.mpf('1e-40')
      and abs(P['s0']/P['L']+mp.mpf('0.5'))<mp.mpf('1e-40') and abs(P['h']-1/(2*phi))<mp.mpf('1e-40'))
kap = {j: P['c1'][j]/P['c0'][j]-1 for j in P['nonloc']}
ok &= abs(kap[1]-1/phi)<mp.mpf('1e-40') and abs(kap[4]+1/phi**2)<mp.mpf('1e-40') and abs(kap[0])<mp.mpf('1e-40')
check("edge params: j0=2,j1=3,nonloc={0,1,4}; alpha=beta=1/(2phi^2); k1=1/phi,k4=-1/phi^2,k0=0; s0=-L/2", ok)
import sympy as sp
t = sp.Symbol('t'); el0,el1,W1,W4,L2s,LFs = sp.symbols('el0 el1 W1 W4 L2 LF')
E5d = sp.sympify(pickle.load(open(os.path.join(HERE,"E5_dec_subs.pkl"),"rb")),
                 locals={'t':t,'el0':el0,'el1':el1,'W1':W1,'W4':W4,'L2':L2s,'LF':LFs})
E5f = sp.lambdify([t,el0,el1,W1,W4,L2s,LFs], E5d, modules='mpmath')
maxd = max(abs(SE.E5_sym(mp.mpf(tv),M,P,conv='B')
    - L*E5f(mp.mpf(tv),mp.log(mp.mpf(tv)),mp.log(1-mp.mpf(tv)),mp.log(1+mp.mpf(tv)/phi),
            mp.log(1-mp.mpf(tv)/phi**2),mp.log(2),mp.log(phi))) for tv in ('0.13','0.3','0.5','0.71','0.88'))
check("exact substituted E5/L vs symeng.E5_sym at 5 points", maxd < mp.mpf('1e-40'), mp.nstr(maxd,3))

# ---------- (2b) assembly from fresh monomial integrals ----------
mp.mp.dps = 100
parsed = pickle.load(open(os.path.join(HERE,"parsed_mons_dec.pkl"),"rb"))
kap1, kap4 = 1/phi, -1/phi**2
def integ(n,a,b,c,d):
    def f(x):
        v = x**n
        if a: v *= mp.log(x)**a
        if b: v *= mp.log(1-x)**b
        if c: v *= mp.log(1+kap1*x)**c
        if d: v *= mp.log(1+kap4*x)**d
        return v
    return mp.quad(f, [0,1], method='tanh-sinh')
need = set()
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    for n in cu: need.add((n,a,b,c,d))
ITAB = {k: integ(*k) for k in sorted(need)}
acc = mp.mpf(0); acc5 = mp.mpf(0)
blocks = {'local': mp.mpf(0),'single-1': mp.mpf(0),'single-4': mp.mpf(0),'mixed': mp.mpf(0)}
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else (('single-1' if c>0 else 'single-4') if (c>0 or d>0) else 'local')
    constfac = mp.log(2)**e * mp.log(phi)**f
    for n,(p,q) in cu.items():
        term = (mp.mpf(str(p))+mp.mpf(str(q))*r5)*constfac*ITAB[(n,a,b,c,d)]
        acc += term
        if sum(mon)==5: acc5 += term; blocks[cls] += term
check("assembly P^(5)/L from 756 fresh monomial integrals vs engine (>= 50 digits)",
      abs(1200*acc - P5L) < mp.mpf('1e-50'), mp.nstr(abs(1200*acc-P5L),3))
print("   W5/L =", mp.nstr(1200*acc5, 55))
print("   blocks: local", mp.nstr(1200*blocks['local'],30))
print("           single-1", mp.nstr(1200*blocks['single-1'],30), " single-4", mp.nstr(1200*blocks['single-4'],30))
print("           mixed", mp.nstr(1200*blocks['mixed'],30))

# ---------- (3) letters and V ----------
import dec5d_lib as D
QS5 = D.QS5; ZERO5 = QS5(0); PHI = D.PHI
letters = {'0': QS5(0), '1': QS5(1), '-phi': -PHI, 'phi^2': PHI**2}
slots_used = set(); diffs = []
nm = list(letters)
for a in nm:
    for b in nm:
        if a >= b: continue
        lv = D.logvec5(letters[a]-letters[b])
        diffs.append((a,b,letters[a]-letters[b],lv))
        for i,v in enumerate(lv):
            if v: slots_used.add(i)
okV = slots_used == {2}
okV &= D.logvec5(2*PHI**2) == (F(1),F(0),F(2)) and D.logvec5(2*PHI) == (F(1),F(0),F(1))
check("letters {0,1,-phi,phi^2}: all differences = +-phi^k; constants 2phi^2, 2phi -> log2, logphi; V=<log2,logphi>, dim V = 2", okV)
for a,b,z,lv in diffs: print("   %s - (%s) = %s = +-phi^%s" % (a,b,z,lv[2]))

# ---------- (4) exact rewrite, numerical check ----------
kit5 = D.make_kit(QS5, D.logvec5, ('log2','logsqrt5','logphi'))
RW5 = D.make_rewrite(QS5, -PHI, PHI**2)
rewrite = {k: RW5['monomial_integral'](*k) for k in sorted(need)}
nwords = set()
for r in rewrite.values(): nwords.update(r.keys())
bad1 = [w for r in rewrite.values() for w,v in r.items() if w and w[0]==QS5(1) and v!=ZERO5]
check("exact rewrite: %d integrals, %d distinct words, weights<=5, no leading-1 words" % (len(rewrite), len(nwords)),
      len(rewrite)==756 and len(nwords)==1024 and max(map(len,nwords))<=5 and not bad1)
import gengine_dec
G1, lets = gengine_dec.gtab(120, 380, cachedir=HERE)
mp.mp.dps = 100
L0,L1n,LP,LQ = lets['L0'],lets['L1'],lets['LP'],lets['LQ']
def q5f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r5
P_LET, Q_LET = -PHI, PHI**2
def w2mpf(w): return tuple(L0 if x==ZERO5 else L1n if x==QS5(1) else LP if x==P_LET else LQ for x in w)
worst = max(abs(sum(q5f(v)*G1[w2mpf(w)] for w,v in r.items()) - ITAB[k]) for k,r in rewrite.items())
check("all 756 rewritten integrals match the fresh table", worst < mp.mpf('1e-55'), mp.nstr(worst,3))

# ---------- (5) controls ----------
obs3, L2_3, L3_3 = D.make_obstruction(kit5)
def obs3_zero(s): return all(v == ZERO5 for row in obs3(s) for v in row)
zs = [QS5(F(1,2)), QS5(2), -QS5(1), 1/PHI, 1/PHI**2, PHI, PHI**2, -PHI, -1/PHI,
      D.SQ5, 1/D.SQ5, D.SQ5/5, 2/PHI, (3-D.SQ5)/2]
okA = all(obs3_zero(kit5['li_symbol'](5, QS5.c(z))) for z in zs)
check("control: Li_5 at 14 S-units (incl. sqrt5 ones) pairs to 0 (d=3, 24 pairings)", okA)
random.seed(17); okB = True
for _ in range(300):
    k = random.choice([1,2,3])
    slots = [0,2] if random.random()<0.7 else [0,1,2]
    w1 = tuple(random.choice(slots) for _ in range(k)); w2 = tuple(random.choice(slots) for _ in range(5-k))
    s1 = {w1: QS5(random.randint(-4,4), random.randint(-3,3))}
    s2 = {w2: QS5(random.randint(-4,4), random.randint(-3,3))}
    if not obs3_zero(kit5['shuffle'](s1,s2)): okB = False; break
check("control: 300 random shuffle products pair to 0", okB)
from qsqrt2 import QS2
kit2 = D.make_kit(QS2, D.logvec2, ('log2','log(1+sqrt2)'))
obs2, L2_2, L3_2 = D.make_obstruction(kit2)
ZERO2 = QS2(0)
symj = json.loads(open(BASE+"M4c_w5_closed_form/symbol.json").read())['symbol']
Sm4 = {tuple(map(int,w)): QS2(F(p), F(q)) for w,(p,q) in symj.items()}
om4 = obs2(Sm4)
check("control: M4c W5 pairs to (84, -38 sqrt2)", om4[0][0]==QS2(84) and om4[1][0]==QS2(0,-38),
      "(%s, %s)" % (om4[0][0], om4[1][0]))
parsed_oct = pickle.load(open(BASE+"OCT5D_octagon_order5_depth/parsed_mons.pkl","rb"))
RW2 = D.make_rewrite(QS2, QS2(2,1), QS2(-1,-1))
needo = {}
for mon, cu in parsed_oct:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    for n in cu: needo[(n,a,b,c,d)] = mon
rewrite2 = {k: RW2['monomial_integral'](*k) for k in sorted(needo)}
gc2 = {}
def Gsym2(w):
    if w not in gc2: gc2[w] = kit2['g_symbol'](kit2['ZERO'], w, kit2['ONE'])
    return gc2[w]
def top2(a,b,c,d,n):
    k = a+b+c+d; out = {}
    for w, gam in rewrite2[(n,a,b,c,d)].items():
        if len(w) != k: continue
        out = kit2['add'](out, kit2['scale'](Gsym2(w), gam))
    return out
def qs2f(p,q): return QS2(F(p.numerator,p.denominator), F(q.numerator,q.denominator))
so = {'local': {}, 'single-p': {}, 'single-q': {}, 'mixed': {}}
for mon, cu in parsed_oct:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else (('single-p' if c>0 else 'single-q') if (c>0 or d>0) else 'local')
    for n,(p,q) in cu.items():
        s = top2(a,b,c,d,n)
        if e: s = kit2['shuffle'](s, kit2['const_pow'](0, e))
        if f: s = kit2['shuffle'](s, kit2['const_pow'](1, f))
        so[cls] = kit2['add'](so[cls], kit2['scale'](s, qs2f(p,q)))
so = {k: kit2['scale'](v, QS2(960)) for k,v in so.items()}
expo = {'local': (QS2(0), QS2(0)),
        'single-p': (QS2(F(-53955,8),F(4770)), QS2(F(-10215,2),F(14445,4))),
        'single-q': (QS2(F(-53955,8),F(4770)), QS2(F(-10215,2),F(14445,4))),
        'mixed': (QS2(F(48015,4),F(-8505)), QS2(F(11115),F(-15525,2)))}
okO = True
for k in expo:
    o = obs2(so[k]); okO &= (o[0][0]==expo[k][0] and o[1][0]==expo[k][1])
oT = obs2(kit2['add'](*so.values()))
okO &= (oT[0][0]==QS2(-1485,1035) and oT[1][0]==QS2(900,-540))
check("control: octagon reproduces all blocks and total (-1485+1035 sqrt2, 900-540 sqrt2)", okO)

# ---------- (6) decagon obstruction ----------
gsy = {}
def Gsym5(w):
    if w not in gsy: gsy[w] = kit5['g_symbol'](kit5['ZERO'], w, kit5['ONE'])
    return gsy[w]
def top5(a,b,c,d,n):
    k = a+b+c+d; out = {}
    for w, gam in rewrite[(n,a,b,c,d)].items():
        if len(w) != k: continue
        out = kit5['add'](out, kit5['scale'](Gsym5(w), gam))
    return out
def qs5f(p,q): return QS5(F(p.numerator,p.denominator), F(q.numerator,q.denominator))
sub = {'local': {}, 'single-1': {}, 'single-4': {}, 'mixed': {}}
for mon, cu in parsed:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else (('single-1' if c>0 else 'single-4') if (c>0 or d>0) else 'local')
    for n,(p,q) in cu.items():
        s = top5(a,b,c,d,n)
        if e: s = kit5['shuffle'](s, kit5['const_pow'](0, e))
        if f: s = kit5['shuffle'](s, kit5['const_pow'](2, f))
        sub[cls] = kit5['add'](sub[cls], kit5['scale'](s, qs5f(p,q)))
sub = {k: kit5['scale'](v, QS5(1200)) for k,v in sub.items()}
tot5 = kit5['add'](*sub.values())
slots_seen = set()
for s in list(sub.values())+[tot5]:
    for w in s: slots_seen.update(w)
ok = slots_seen == {0,2}
check("decagon weight-5 symbols use only slots {log2, logphi}", ok)
ok3 = all(obs3_zero(sub[k]) for k in sub) and obs3_zero(tot5)
check("d=3 obstruction (8 Lie_3 x 3 Lie_2 = 24 pairings) vanishes on every block and total", ok3)
def remap(s): return {tuple(0 if x==0 else 1 for x in w): c for w,c in s.items()}
kit5r = D.make_kit(QS5, D.logvec5, ('log2','logphi'))
obs2d, _, _ = D.make_obstruction(kit5r)
res2 = {}
ok2 = True
for k in sub:
    o = obs2d(remap(sub[k])); res2[k] = o
    ok2 &= (o[0][0]==ZERO5 and o[1][0]==ZERO5)
oT2 = obs2d(remap(tot5))
ok2 &= (oT2[0][0]==ZERO5 and oT2[1][0]==ZERO5)
check("d=2 M4c pairings ([X,[X,Y]]x[X,Y], [Y,[X,Y]]x[X,Y]) vanish on every block and total", ok2)
X = {(0,): QS5(1)}; Y = {(1,): QS5(1)}
br = kit5r['bracket']
L2e = br(X,Y)
L4e = [br(X, br(X, L2e)), br(Y, br(Y, L2e)), br(Y, br(X, L2e))]
ok41 = True
for s in list(sub.values())+[tot5]:
    sr = remap(s)
    for p in L4e:
        for q in (X, Y):
            v = sum((c*dd*(sr.get(u+v2,ZERO5)-sr.get(v2+u,ZERO5)) for u,c in p.items() for v2,dd in q.items()), ZERO5)
            if v != ZERO5: ok41 = False
check("(4,1) cobracket component (Lie_4 x Lie_1) vanishes on every block and total", ok41)

# ---------- (6b) classical span fit ----------
def logvec_dec(z):
    lv = D.logvec5(z)
    if lv[1] != 0: raise ValueError("logsqrt5 part")
    return (lv[0], lv[2])
kit_dec = D.make_kit(QS5, logvec_dec, ('log2','logphi'))
sh = kit_dec['shuffle']; li = kit_dec['li_symbol']
words5 = list(itertools.product([0,1], repeat=5)); widx = {w:i for i,w in enumerate(words5)}
gens = []
for i in (0,1):
    for w in itertools.product([0,1], repeat=4):
        gens.append(sh({(i,): QS5(1)}, {w: QS5(1)}))
for u in itertools.product([0,1], repeat=2):
    for v in itertools.product([0,1], repeat=3):
        gens.append(sh({u: QS5(1)}, {v: QS5(1)}))
cands = []
for a in range(-4,5):
    for b in range(-6,7):
        for sgn in (1,-1):
            z = (QS5(2)**a)*(PHI**b)*(QS5(sgn))
            if z in (QS5(0), QS5(1)): continue
            try: logvec_dec(z); logvec_dec(QS5(1)-z)
            except ValueError: continue
            cands.append(z)
for z in set(cands):
    gens.append(li(5, z))
def to_vec(s):
    v = [F(0)]*32
    for w, c in s.items():
        assert c.q == 0
        v[widx[w]] = F(c.p)
    return v
G = [to_vec(g) for g in gens]
def solve2(G, t):
    M = [[G[j][i] for j in range(len(G))] + [t[i]] for i in range(32)]
    piv = {}; r = 0
    for c in range(len(G)):
        p = next((i for i in range(r, 32) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        inv = 1/M[r][c]; M[r] = [x*inv for x in M[r]]
        for i in range(32):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [M[i][j]-f*M[r][j] for j in range(len(G)+1)]
        piv[c] = r; r += 1
        if r == 32: break
    for i in range(r, 32):
        if M[i][len(G)] != 0: return None
    sol = [F(0)]*len(G)
    for c, i in piv.items(): sol[c] = M[i][len(G)]
    for i in range(32):
        if sum(sol[j]*G[j][i] for j in range(len(G))) != t[i]: return None
    return sol
okS = True
for k in list(sub)+['tot']:
    tk = remap(sub[k] if k!='tot' else tot5)
    tp = [F(0)]*32; tq = [F(0)]*32
    for w,c in tk.items(): tp[widx[w]] = F(c.p); tq[widx[w]] = F(c.q)
    sp_, sq_ = solve2(G,tp), solve2(G,tq)
    okS &= (sp_ is not None and sq_ is not None)
    nli = sum(1 for j in range(len(gens)) if (sp_ and sp_[j] or sq_ and sq_[j]))
    print("   block %-9s in classical span: p=%s q=%s (verified exact)" % (k, sp_ is not None, sq_ is not None))
check("every block and the total lie EXACTLY in the classical span (shuffle products + Li_5), residual 0", okS)

# ---------- (7) local block exact decomposition ----------
acc = pickle.load(open(os.path.join(HERE,"local_block_exact.pkl"),"rb"))
ZVm = {'1': mp.mpf(1), 'z2': z2, 'z3': z3, 'z2^2': z2**2, 'z5': mp.zeta(5), 'z2z3': z2*z3}
def q5mpf(z): return mp.mpf(z.p.numerator)/mp.mpf(z.p.denominator) + (mp.mpf(z.q.numerator)/mp.mpf(z.q.denominator))*r5
vloc = sum(q5mpf(acc[k])*ZVm[k[0]]*mp.log(2)**k[1]*mp.log(phi)**k[2] for k in acc)
check("exact local-block decomposition vs fresh assembly (>= 50 digits)",
      abs(vloc - 1200*blocks['local']) < mp.mpf('1e-50'), mp.nstr(abs(vloc-1200*blocks['local']),3))
z5c = acc.get(('z5',0,0)); z23c = acc.get(('z2z3',0,0))
check("local zeta(5) coeff = paper prop:order5 C_5^(5)/L = 73125/64-30825 sqrt5/64",
      z5c == QS5(F(73125,64), F(-30825,64)), str(z5c))
check("local zeta(2)zeta(3) coeff = paper prop:order5 C_5^(2,3)/L = -64875/64+32775 sqrt5/64",
      z23c == QS5(F(-64875,64), F(32775,64)), str(z23c))

print("\\ntotal time %.1fs" % (time.time()-T0))
if failures:
    print("FAILURES:", failures); sys.exit(1)
print("ALL CHECKS PASS")
'''
open("verify_dec5d.py","w").write(verifier)
print("written verify_dec5d.py", len(verifier))
