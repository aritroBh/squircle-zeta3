"""verify_oct5d2.py — OCT5D2 depth-test verifier (< 5 min).

Reproduces, from the read-only inputs only:
  (1) the exact Q(sqrt2) G-rewrite of the 756 weight-5 monomial integrals and its numerical
      check at >= 50 digits (mixed block vs the quoted 4.2299562389856442293498943256846704703130475386317
      and vs the itab110-assembled 110-digit target; single and local blocks too);
  (2) the controls: Li_5 at 7 S-units pairs to 0; 300 random shuffle products pair to 0;
      M4c W5 pairs to (84, -38 sqrt2); G-symbol identities on test words;
  (3) the two antisymmetrised (3,2) pairings of the weight-5 part of P^(5)(1)/u, per block
      and total, against [X,[X,Y]] x [X,Y] and [Y,[X,Y]] x [X,Y].

Read-only inputs (paths below): parsed_mons.pkl, itab110.pkl (OCT5D); symbol.json (M4c);
validated exact-arithmetic/symbol libraries qsqrt2.py (paper scripts) and w4_symbol.py (A3_W4_proof).
Local modules: oct5d2_engine.py (exact rewrite), gengine.py (numerical G-table), symengine.py
(symbol wrapper). Exits nonzero on any failure.
"""
import os, sys, time, pickle, json, math, random
import mpmath as mp
from fractions import Fraction as F

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
OCT5D = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
M4C   = "<project>/research/08_theory/M4c_w5_closed_form/"
PAPER_SCRIPTS = "<project>/paper/squircle_note/scripts"
A3    = "<project>/research/08_theory/A3_W4_proof"
for p in (HERE, PAPER_SCRIPTS, A3):
    if p not in sys.path: sys.path.insert(0, p)

failures = []
def check(name, cond, detail=""):
    print(("PASS  " if cond else "FAIL  ") + name + ("   " + detail if detail else ""))
    if not cond: failures.append(name)

# ---------- exact rewrite ----------
import oct5d2_engine as E
from oct5d2_engine import Q, ZERO, ONE, P_LET, Q_LET
parsed = pickle.load(open(OCT5D+"parsed_mons.pkl","rb"))
need = {}
for mon, cu in parsed:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else ('single' if (c>0 or d>0) else 'local')
    for n in cu: need[(n,a,b,c,d)] = cls
rewrite = {k: E.monomial_integral(*k) for k in sorted(need)}
nwords = set()
for r in rewrite.values(): nwords.update(r.keys())
check("exact rewrite: %d integrals, %d distinct words, weights<=5" % (len(rewrite), len(nwords)),
      len(rewrite) == 756 and max(map(len, nwords)) <= 5)
bad1 = [w for r in rewrite.values() for w,v in r.items() if w and w[0]==ONE and v!=ZERO]
check("no word starting with 1 has a nonzero coefficient (regularisation consistent)", not bad1)

# ---------- numerical G-table (fresh build, dps 90, N 320) ----------
import gengine
DPS, NSER = 90, 320
G1, letters = gengine.gtab(DPS, NSER, cachedir=HERE)
mp.mp.dps = DPS
r2 = mp.sqrt(2); L2C = mp.log(2); LU = mp.log(1+r2)
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r2
def w2mpf(w): return tuple(L0 if x==ZERO else L1 if x==ONE else LP if x==P_LET else LQ for x in w)
# G-table self-checks (independent of the rewrite)
check("G-table: G(0^4,p;1) = -Li_5(1/p)",
      abs(G1[(L0,)*4+(LP,)] + mp.polylog(5, 1/LP)) < mp.mpf('1e-60'),
      mp.nstr(abs(G1[(L0,)*4+(LP,)] + mp.polylog(5, 1/LP)), 3))
gpp_qq = mp.quad(lambda t: mp.log(1-t/LQ)/(t-LP), [0,1])
check("G-table: G(p,q;1) vs direct quadrature", abs(G1[(LP,LQ)]-gpp_qq) < mp.mpf('1e-60'),
      mp.nstr(abs(G1[(LP,LQ)]-gpp_qq), 3))

ival = {k: sum(q2f(v)*G1[w2mpf(w)] for w, v in r.items()) for k, r in rewrite.items()}
ITAB = {k: mp.mpf(v) for k, v in pickle.load(open(OCT5D+"itab110.pkl","rb")).items()}
worst = max(abs(ival[(n,a,b,c,d)]-ITAB[(n,a,b,c,d)]) for (n,a,b,c,d) in rewrite)
check("all 756 rewritten integrals match itab110 (its own 110-digit accuracy)", worst < mp.mpf('1e-60'),
      "worst "+mp.nstr(worst,3))

blocks = {'local': mp.mpf(0), 'single': mp.mpf(0), 'mixed': mp.mpf(0)}
for mon, cu in parsed:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else ('single' if (c>0 or d>0) else 'local')
    constfac = L2C**e * LU**f
    for n,(p,q) in cu.items():
        blocks[cls] += (mp.mpf(str(p)) + mp.mpf(str(q))*r2)*constfac*ival[(n,a,b,c,d)]
blocks = {k: 960*v for k,v in blocks.items()}
# independent 110-digit targets, assembled inline from the read-only monomial table + itab110
mp.mp.dps = 120
tcontrib = {'local': mp.mpf(0), 'single': mp.mpf(0), 'mixed': mp.mpf(0)}
for mon, cu in parsed:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else ('single' if (c>0 or d>0) else 'local')
    constfac = L2C**e * LU**f
    for n,(p,q) in cu.items():
        tcontrib[cls] += (mp.mpf(str(p)) + mp.mpf(str(q))*r2)*constfac*ITAB[(n,a,b,c,d)]
targets = {k: 960*v for k,v in tcontrib.items()}
mp.mp.dps = DPS
for k in ('local','single','mixed'):
    check("rewrite check (%s block, >=50 digits)" % k,
          abs(blocks[k]-targets[k]) < mp.mpf('1e-50'), "residual "+mp.nstr(abs(blocks[k]-targets[k]),3))
quoted = mp.mpf("4.2299562389856442293498943256846704703130475386317")
check("mixed block vs task-quoted value (49 digits as quoted)", abs(blocks['mixed']-quoted) < mp.mpf('1e-49'),
      "residual "+mp.nstr(abs(blocks['mixed']-quoted),3))

# ---------- symbols and controls ----------
import symengine as S
from symengine import add, scale, shuffle, gsym, const_pow, obstruction, li_symbol
R2 = Q(0,1); B = Q(-1,1)
okA = all(obstruction(li_symbol(5, z)) == [ZERO, ZERO]
          for z in [B, -B, -B*B, R2/2, Q(F(1,2)), ONE-R2/2, 2-R2])
check("control: Li_5 at 7 S-units pairs to 0", okA)
random.seed(11); okB = True
for _ in range(300):
    k = random.choice([1,2])
    w1 = tuple(random.choice([0,1]) for _ in range(k))
    w2 = tuple(random.choice([0,1]) for _ in range(5-k))
    s1 = {w1: Q(random.randint(-4,4), random.randint(-3,3))}
    s2 = {w2: Q(random.randint(-4,4), random.randint(-3,3))}
    if obstruction(shuffle(s1,s2)) != [ZERO, ZERO]: okB = False; break
check("control: 300 random shuffle products pair to 0", okB)
symj = json.loads(open(M4C+"symbol.json").read())['symbol']
Sm4 = {tuple(map(int,w)): Q(F(p), F(q)) for w,(p,q) in symj.items()}
obs_m4 = obstruction(Sm4)
check("control: M4c W5 pairs to (84, -38 sqrt2)", obs_m4 == [Q(84), Q(0,-38)], str(obs_m4))
s1 = gsym((ZERO, P_LET)); ref = scale(li_symbol(2, ONE/P_LET), -1)
okD = (add(s1, scale(ref,-1)) == {})
s2a = shuffle(gsym((P_LET,)), gsym((Q_LET,)))
s2b = add(gsym((P_LET,Q_LET)), gsym((Q_LET,P_LET)))
okD &= (add(s2a, scale(s2b,-1)) == {})
s3 = gsym((ZERO,)*4+(P_LET,)); ref3 = scale(li_symbol(5, ONE/P_LET), -1)
okD &= (add(s3, scale(ref3,-1)) == {})
check("control: G-symbol identities (0,p)=-Li_2(1/p), (p,q)-shuffle, (0^4,p)=-Li_5(1/p)", okD)

# ---------- block symbols and pairings ----------
symcache = {}
def Gsym(w):
    if w not in symcache: symcache[w] = gsym(w)
    return symcache[w]
def top_symbol(a,b,c,d,n):
    k = a+b+c+d
    out = {}
    for w, gam in rewrite[(n,a,b,c,d)].items():
        if len(w) != k: continue
        out = add(out, scale(Gsym(w), gam))
    return out
def qs2(p,q): return Q(F(p.numerator,p.denominator), F(q.numerator,q.denominator))
sub_sym = {'local': {}, 'single-p': {}, 'single-q': {}, 'mixed': {}}
for mon, cu in parsed:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else (('single-p' if c>0 else 'single-q') if (c>0 or d>0) else 'local')
    for n,(p,q) in cu.items():
        s = top_symbol(a,b,c,d,n)
        if e: s = shuffle(s, const_pow(0, e))
        if f: s = shuffle(s, const_pow(1, f))
        sub_sym[cls] = add(sub_sym[cls], scale(s, qs2(p,q)))
sub_sym = {k: scale(v, Q(960)) for k,v in sub_sym.items()}
exp = {'local': (Q(0), Q(0)),
       'single-p': (Q(F(-53955,8), F(4770)), Q(F(-10215,2), F(14445,4))),
       'single-q': (Q(F(-53955,8), F(4770)), Q(F(-10215,2), F(14445,4))),
       'mixed':   (Q(F(48015,4), F(-8505)), Q(F(11115), F(-15525,2)))}
obs = {}
for k in ('local','single-p','single-q','mixed'):
    o = obstruction(sub_sym[k]); obs[k] = o
    check("pairing (%s block)" % k, o[0]==exp[k][0] and o[1]==exp[k][1], "(%s, %s)" % (o[0], o[1]))
tot = add(sub_sym['local'], sub_sym['single-p'], sub_sym['single-q'], sub_sym['mixed'])
oT = obstruction(tot)
check("pairing (total W5/u) = (-1485+1035 sqrt2, 900-540 sqrt2)",
      oT[0]==Q(-1485,1035) and oT[1]==Q(900,-540), "(%s, %s)" % (oT[0], oT[1]))
check("TOTAL OBSTRUCTION IS NONZERO (depth two needed, symbol level)", oT != [ZERO,ZERO])
rA = oT[0]/Q(84); rB = oT[1]/Q(0,-38)
print("ratio to M4c: pair1/84 =", rA, " pair2/(-38 sqrt2) =", rB, " proportional:", rA==rB)
check("not proportional to M4c (84,-38 sqrt2)", rA != rB)

print("\ntotal time %.1fs" % (time.time()-T0))
if failures:
    print("FAILURES:", failures); sys.exit(1)
print("ALL CHECKS PASS")
