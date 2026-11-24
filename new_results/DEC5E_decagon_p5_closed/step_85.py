verifier = '''"""verify_dec5e.py - end-to-end verification of the decagon P^{(5)}(1) closed form.

Recomputes everything at 250 digits (G-table d250/N850, ~15 s) and checks:
 P1  weight parts sum = P^{(5)}(1)/L gold (113 digits)
 P2  Wt0 = 0, Wt1/Wt2 exact closed forms
 P3  Wt3 closed form (8 atoms)
 P4  Wt4 closed form (15 atoms)
 P5  Wt5 = exact log-part + 22-atom remainder (classical golden + L(chi_5,5))
 P6  local block: zeta(5), zeta(2)zeta(3) = prop:order5 at m=5 (exact sympy)
 P7  L(chi_5,5): Hurwitz vs Euler product
 P8  the 16 golden Li32/Li41 K-parts reduce to the 22-atom basis (LLL)
Expected runtime: < 5 min. Exit 0 iff all pass.
"""
import sys, os, pickle, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "<project>/research/08_theory/DEC5D_decagon_order5_depth")
import mpmath as mp
import dec5d_lib as DL
import gengine_dec2p

T0 = time.time()
PASS = []
def report(name, ok, detail=""):
    PASS.append(ok)
    print("%s %-72s %s" % ("PASS" if ok else "FAIL", name, detail))

DPS, N = 250, 850
G1, letters = gengine_dec2p.gtab(DPS, N, cachedir=HERE, procs=10)
mp.mp.dps = DPS
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
QS5 = DL.QS5; ZERO, ONE = QS5(0), QS5(1)
P_LET, Q_LET = -DL.PHI, DL.PHI**2
r5 = mp.sqrt(5); l2 = mp.log(2); lf = mp.log((1+r5)/2); z2 = mp.zeta(2); z3 = mp.zeta(3); z5 = mp.zeta(5)
def q5f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r5
def w2mpf(w): return tuple(L0 if x==ZERO else L1 if x==ONE else LP if x==P_LET else LQ for x in w)
rewrite = pickle.load(open("<project>/research/08_theory/DEC5D_decagon_order5_depth/rewrite_dec_exact.pkl","rb"))
parsed = pickle.load(open("<project>/research/08_theory/DEC5D_decagon_order5_depth/parsed_mons_dec.pkl","rb"))
Gv = {}
def gval(w):
    if w not in Gv: Gv[w] = G1[w2mpf(w)]
    return Gv[w]
Wt = {k: mp.mpf(0) for k in range(6)}
blk = {k: {c: mp.mpf(0) for c in ('local','single-1','single-4','mixed')} for k in range(6)}
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else (('single-1' if c>0 else 'single-4') if (c>0 or d>0) else 'local')
    cf = l2**e * lf**f
    for n,(p,q) in cu.items():
        pf = (mp.mpf(str(p)) + mp.mpf(str(q))*r5)*cf
        for w, gam in rewrite[(n,a,b,c,d)].items():
            k = len(w) + e + f
            if k <= 5:
                v = pf*q5f(gam)*gval(w); Wt[k] += v; blk[k][cls] += v
Wt = {k: 1200*v for k, v in Wt.items()}
blk = {k: {c: 1200*v for c, v in row.items()} for k, row in blk.items()}

# P1: sum vs gold
gold = pickle.load(open("<project>/research/08_theory/DEC5D_decagon_order5_depth/p5_gold.pkl","rb"))
d1 = abs(mp.mpf(gold['P5overL']) - sum(Wt.values()))
report("P1 sum Wt0..Wt5 = P^(5)(1)/L (113-digit gold)", d1 < mp.mpf(10)**-110, "|diff| = %s" % mp.nstr(d1, 2))

# P2: exact cheap parts
Wt1_exact = 1200*l2 + (3000-600*r5)*lf
Wt2_exact = (mp.mpf(-1425)/16 - mp.mpf(5415)/16*r5)*z2 - 2400*l2**2 + (-12000+2400*r5)*l2*lf + (mp.mpf(-86925)/4 + mp.mpf(17325)/2*r5)*lf**2
report("P2a Wt0 = 0", abs(Wt[0]) < mp.mpf(10)**-200, "|Wt0| = %s" % mp.nstr(abs(Wt[0]), 2))
report("P2b Wt1 = 1200 log2 + (3000-600 sqrt5) log phi", abs(Wt[1]-Wt1_exact) < mp.mpf(10)**-200, "|diff| = %s" % mp.nstr(abs(Wt[1]-Wt1_exact), 2))
report("P2c Wt2 exact 4-term form", abs(Wt[2]-Wt2_exact) < mp.mpf(10)**-200, "|diff| = %s" % mp.nstr(abs(Wt[2]-Wt2_exact), 2))

# atoms
L3m = mp.polylog(3, -1/phi); L4m = mp.polylog(4, -1/phi); L4q = mp.polylog(4, 1/phi**2)
L5m = mp.polylog(5, -1/phi); L5q = mp.polylog(5, 1/phi**2)
K5m = L5m + lf**5/120; K5q = L5q - mp.mpf(2)/15*lf**5
K4m = L4m - lf**4/24;  K4q = L4q + lf**4/3
K3m = L3m + lf**3/6
def Lchi5(s, D):
    mp.mp.dps = D
    return mp.mpf(5)**(-s)*(mp.zeta(s, mp.mpf(1)/5) - mp.zeta(s, mp.mpf(2)/5) - mp.zeta(s, mp.mpf(3)/5) + mp.zeta(s, mp.mpf(4)/5))
Lchi = Lchi5(5, DPS+30); mp.mp.dps = DPS
A = {'l2': l2, 'lf': lf, 'z2': z2, 'z3': z3, 'z5': z5, 'L3m': L3m, 'L4m': L4m, 'L4q': L4q,
     'K5m': K5m, 'K5q': K5q, 'K4m': K4m, 'K4q': K4q, 'K3m': K3m, 'Lchi5': Lchi}
def V(spec):
    v = mp.mpf(1)
    for nm, pw in spec: v *= A[nm]**pw
    return v
def cf(p, q): return mp.mpf(p) + mp.mpf(q)*r5

# P3: Wt3
W3F = [(('L3m',1),), ("5250","-1500")] 
W3F = [([('L3m',1)], ("5250","-1500")), ([('l2',3)], ("1200","0")), ([('l2',2),('lf',1)], ("9000","-1800")),
       ([('l2',1),('lf',2)], ("64575/2","-25725/2")), ([('lf',3)], ("104225/2","-187325/8")),
       ([('z2',1),('l2',1)], ("300","795/2")), ([('z2',1),('lf',1)], ("16275/8","-8535/8")),
       ([('z3',1)], ("3675/4","150"))]
v3 = sum(cf(*pq)*V(spec) for spec, pq in W3F)
report("P3 Wt3 closed form (8 atoms)", abs(Wt[3]-v3) < mp.mpf(10)**-200, "|diff| = %s" % mp.nstr(abs(Wt[3]-v3), 2))

# P4: Wt4
W4F = [([('L4m',1)], ("-21750","10500")), ([('L4q',1)], ("263625/32","-60975/16")),
       ([('L3m',1),('l2',1)], ("-5250","1800")), ([('L3m',1),('lf',1)], ("-51375/2","23325/2")),
       ([('l2',4)], ("-200","0")), ([('l2',3),('lf',1)], ("-2000","400")),
       ([('l2',2),('lf',2)], ("-11550","4650")), ([('l2',1),('lf',3)], ("-176025/4","79075/4")),
       ([('lf',4)], ("-2395375/32","1076375/32")),
       ([('z2',1),('l2',2)], ("-75","-135")), ([('z2',1),('l2',1),('lf',1)], ("-2550","1185")),
       ([('z2',1),('lf',2)], ("-78525/16","30255/16")),
       ([('z3',1),('l2',1)], ("-375","-90")), ([('z3',1),('lf',1)], ("-4200","2625/2")),
       ([('z2',2)], ("-106965/16","54285/16"))]
v4 = sum(cf(*pq)*V(spec) for spec, pq in W4F)
report("P4 Wt4 closed form (15 atoms)", abs(Wt[4]-v4) < mp.mpf(10)**-200, "|diff| = %s" % mp.nstr(abs(Wt[4]-v4), 2))

# P5: Wt5 = log part + 22-atom remainder
LOGPART = [([('l2',5)], ("10","0")), ([('l2',4),('lf',1)], ("125","-25")),
           ([('l2',3),('lf',2)], ("2025/2","-825/2")), ([('l2',2),('lf',3)], ("52075/8","-23525/8")),
           ([('l2',1),('lf',4)], ("449025/16","-403175/32")), ([('lf',5)], ("15112175/256","-6759335/256"))]
logval = sum(cf(*pq)*V(spec) for spec, pq in LOGPART)
fitR = pickle.load(open(os.path.join(HERE, 'fit_R_lchi.pkl'),'rb'))
SPEC = {'K3ml2^2':[('K3m',1),('l2',2)],'K3ml2lf':[('K3m',1),('l2',1),('lf',1)],'K3mlf^2':[('K3m',1),('lf',2)],
        'K4ml2':[('K4m',1),('l2',1)],'K4mlf':[('K4m',1),('lf',1)],'K4ql2':[('K4q',1),('l2',1)],'K4qlf':[('K4q',1),('lf',1)],
        'K5m':[('K5m',1)],'K5q':[('K5q',1)],'Lchi5':[('Lchi5',1)],'z2L3m':[('z2',1),('L3m',1)],
        'z2^2l2':[('z2',2),('l2',1)],'z2^2lf':[('z2',2),('lf',1)],'z2l2^2lf':[('z2',1),('l2',2),('lf',1)],
        'z2l2^3':[('z2',1),('l2',3)],'z2l2lf^2':[('z2',1),('l2',1),('lf',2)],'z2lf^3':[('z2',1),('lf',3)],
        'z2z3':[('z2',1),('z3',1)],'z3l2^2':[('z3',1),('l2',2)],'z3l2lf':[('z3',1),('l2',1),('lf',1)],
        'z3lf^2':[('z3',1),('lf',2)],'z5':[('z5',1)]}
R_val = sum(cf(str(p), str(q))*V(SPEC[nm]) for nm, (p, q) in fitR['rel'].items())
v5 = logval + R_val
report("P5 Wt5 = exact log-part + 22-atom remainder (classical + L(chi_5,5))",
       abs(Wt[5]-v5) < mp.mpf(10)**-200, "|diff| = %s" % mp.nstr(abs(Wt[5]-v5), 2))

# P6: prop:order5 at m=5 vs local block coefficients (exact sympy)
import sympy as sp
sr5 = sp.sqrt(5)
s_ = (sr5-1)/4; m = 5
A5 = sp.expand(2880*m*s_**4*(65*s_**6 - 63*s_**4 + 17*s_**2 - 1))
A23 = sp.expand(-960*m*s_**4*(45*s_**6 - 27*s_**4 + 6*s_**2 - 1))
def pq(e):
    e = sp.expand(e)
    return (sp.expand((e + e.subs(sr5,-sr5))/2), sp.expand((e - e.subs(sr5,-sr5))/(2*sr5)))
p5_, q5_ = pq(sp.expand(A5*(sp.Rational(5,2)+sr5)))
p23_, q23_ = pq(sp.expand(A23*(sp.Rational(5,2)+sr5)))
ok5 = (p5_, q5_) == (sp.Rational(73125,64), sp.Rational(-30825,64))
ok23 = (p23_, q23_) == (sp.Rational(-64875,64), sp.Rational(32775,64))
# numeric check of the local block fit coefficients
lp = pickle.load(open(os.path.join(HERE, 'block_logparts.pkl'),'rb'))
logloc = sum(cf(p, q)*l2**a*lf**(5-a) for a, (p, q) in lp['local'].items())
Rloc = blk[5]['local'] - logloc
# coefficients from the full Wt5 fit restricted to ring atoms present in local block:
z5c = fitR['rel']['z5']; z23c = fitR['rel']['z2z3']
# local block's OWN fit values (DEC5D proved, 1.25e-239): compare numerically via exact formulas
loc_z5 = mp.mpf(73125)/64 - mp.mpf(30825)/64*r5
loc_z23 = -mp.mpf(64875)/64 + mp.mpf(32775)/64*r5
# verify local R decomposition with ring atoms (coefficients from DEC5D local_block_exact.pkl)
lb = pickle.load(open("<project>/research/08_theory/DEC5D_decagon_order5_depth/local_block_exact.pkl","rb"))
locval = mp.mpf(0)
for (nm, a, b), c in lb.items():
    if nm == 'z5': v = z5
    elif nm == 'z2z3': v = z2*z3
    elif nm == 'z2^2': v = z2**2
    elif nm == 'z2': v = z2
    elif nm == 'z3': v = z3
    else: v = mp.mpf(1)
    locval += (mp.mpf(str(c.p)) + mp.mpf(str(c.q))*r5) * v * l2**a * lf**b
d_loc = abs(blk[5]['local'] - locval)
report("P6a prop:order5 m=5 zeta(5), zeta(2)zeta(3) (exact sympy)", ok5 and ok23, "exact")
report("P6b local block = DEC5D 37-term decomposition (incl. prop:order5 coeffs)", d_loc < mp.mpf(10)**-200, "|diff| = %s" % mp.nstr(d_loc, 2))

# P7: L(chi_5,5) Hurwitz vs Euler product
mp.mp.dps = 60
def leg(p):
    r = p % 5
    return 1 if r in (1,4) else (-1 if r in (2,3) else 0)
def primes_upto(n):
    sieve = [True]*(n+1); sieve[0]=sieve[1]=False
    for i in range(2,int(n**0.5)+1):
        if sieve[i]:
            for j in range(i*i,n+1,i): sieve[j]=False
    return [p for p in range(2,n+1) if sieve[p]]
prod = mp.mpf(1)
for p in primes_upto(500000):
    prod *= 1/(1 - leg(p)*mp.mpf(p)**-5)
mp.mp.dps = DPS
d7 = abs(prod - Lchi)
report("P7 L(chi_5,5): Hurwitz zeta vs Euler product (p<=5e5)", d7 < mp.mpf(10)**-28, "|diff| = %s" % mp.nstr(d7, 2))

# P8: the 16 golden K-parts reduce to the 22-atom basis
from q5_lll import fit as q5fit
import dec5d_lib as DL2
kit = DL.make_kit(QS5, DL.logvec5, slots=('log2','logsqrt5','logphi'))
gs = kit['g_symbol']
letter = {'1': ONE, 'p': P_LET, 'q': Q_LET}
INV22 = ['z2l2^3','z2l2^2lf','z2l2lf^2','z2lf^3','z3l2^2','z3l2lf','z3lf^2','z2^2l2','z2^2lf',
         'z5','z2z3','z2L3m','K5m','K5q','K4ml2','K4mlf','K4ql2','K4qlf','K3ml2^2','K3ml2lf','K3mlf^2','Lchi5']
A22 = {'z2l2^3': z2*l2**3, 'z2l2^2lf': z2*l2**2*lf, 'z2l2lf^2': z2*l2*lf**2, 'z2lf^3': z2*lf**3,
       'z3l2^2': z3*l2**2, 'z3l2lf': z3*l2*lf, 'z3lf^2': z3*lf**2, 'z2^2l2': z2**2*l2, 'z2^2lf': z2**2*lf,
       'z5': z5, 'z2z3': z2*z3, 'z2L3m': z2*L3m, 'K5m': K5m, 'K5q': K5q,
       'K4ml2': K4m*l2, 'K4mlf': K4m*lf, 'K4ql2': K4q*l2, 'K4qlf': K4q*lf,
       'K3ml2^2': K3m*l2**2, 'K3ml2lf': K3m*l2*lf, 'K3mlf^2': K3m*lf**2, 'Lchi5': Lchi}
n_ok = 0
for (a_i, b_i) in [(3,2),(4,1)]:
    for At in '1pq':
        for Bt in '1pq':
            if At == '1' and Bt == '1': continue
            w = (ZERO,)*(a_i-1) + (letter[At],) + (ZERO,)*(b_i-1) + (letter[Bt],)
            s = gs(ZERO, w, ONE)
            c = s[(2,)*5]
            Kv = gval(w) - (mp.mpf(str(c.p)) + mp.mpf(str(c.q))*r5)/120*lf**5
            cands = q5fit(Kv, [(x, A22[x]) for x in INV22], cexp=int(DPS*0.92), nbest=1)
            if cands[0][1] < mp.mpf(10)**-180: n_ok += 1
report("P8 all 16 golden Li32/Li41 K-parts reduce to the 22-atom basis", n_ok == 16, "%d/16 reduce" % n_ok)

print()
print("elapsed %.1f s" % (time.time()-T0))
print("ALL CHECKS PASSED" if all(PASS) else "SOME CHECKS FAILED")
sys.exit(0 if all(PASS) else 1)
'''
open('verify_dec5e.py','w').write(verifier)
import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run(['nice', PY, 'verify_dec5e.py'], capture_output=True, text=True, cwd='.')
print(r.stdout[-3000:])
if r.stderr: print("STDERR:", r.stderr[-500:])
print("wall %.0fs" % (time.time()-t0))
