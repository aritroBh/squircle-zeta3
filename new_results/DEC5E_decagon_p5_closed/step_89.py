ver = r'''"""verify_dec5e.py - end-to-end verification of DEC5E (< 5 min).

Rebuilds the decagon weight parts Wt0..Wt5 of P^(5)(1)/L from DEC5D's exact G-rewrite
on a fresh 300-digit G-table (parallel engine), then checks:
  P1  sum Wt_k == P^(5)(1)/L gold (DEC5D 113-digit value)
  P2  Wt0 = 0 exactly, Wt1/Wt2 = exact closed forms (QS5 exact arithmetic, proved word values)
  P3  Wt3 closed form (8 atoms)
  P4  Wt4 closed form (15 atoms)
  P5  Wt5 = exact symbol-proved log part + (21 classical atoms + L(chi_5,5)) remainder
  P6  local-block zeta(5), zeta(2)zeta(3) coefficients == prop:order5 at m = 5 (exact sympy)
  P7  L(chi_5,5): Hurwitz-zeta value vs Euler product (50 digits)
Exit 0 iff all PASS.
"""
import sys, os, pickle, time
HERE = os.path.dirname(os.path.abspath(__file__))
DEC5D = "<project>/research/08_theory/DEC5D_decagon_order5_depth"
sys.path.insert(0, HERE); sys.path.insert(0, DEC5D)
import mpmath as mp
import dec5d_lib as DL

def main():
    t0 = time.time()
    PASS = []
    def report(name, ok, extra=""):
        PASS.append(ok)
        print("%s %-58s %s" % ("PASS" if ok else "FAIL", name, extra))

    # ---------- fresh G-table + weight parts ----------
    DPS, NTR = 300, 1000
    import gengine_dec2p
    G1, letters = gengine_dec2p.gtab(DPS, NTR, cachedir=HERE, procs=10)
    mp.mp.dps = DPS
    QS5 = DL.QS5; ZERO, ONE = QS5(0), QS5(1)
    P_LET, Q_LET = -DL.PHI, DL.PHI**2
    L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
    r5 = mp.sqrt(5); l2 = mp.log(2); lf = mp.log((1+r5)/2); z2 = mp.zeta(2); z3 = mp.zeta(3); z5 = mp.zeta(5)
    def q5f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r5
    def w2mpf(w): return tuple(L0 if x==ZERO else L1 if x==ONE else LP if x==P_LET else LQ for x in w)
    rewrite = pickle.load(open(os.path.join(DEC5D, "rewrite_dec_exact.pkl"),"rb"))
    parsed = pickle.load(open(os.path.join(DEC5D, "parsed_mons_dec.pkl"),"rb"))
    Gv = {}
    def gval(w):
        if w not in Gv: Gv[w] = G1[w2mpf(w)]
        return Gv[w]
    Wt = {k: mp.mpf(0) for k in range(6)}
    for mon, cu in parsed:
        a,b,c,d,e,f = mon
        cf = l2**e * lf**f
        for n,(p,q) in cu.items():
            pf = (mp.mpf(str(p)) + mp.mpf(str(q))*r5)*cf
            for w, gam in rewrite[(n,a,b,c,d)].items():
                k = len(w) + e + f
                if k <= 5:
                    Wt[k] += pf*q5f(gam)*gval(w)
    Wt = {k: 1200*v for k, v in Wt.items()}
    gold = pickle.load(open(os.path.join(DEC5D, "p5_gold.pkl"),"rb"))
    P5L = mp.mpf(gold['P5overL'])
    d1 = abs(P5L - sum(Wt.values()))
    report("P1 sum(Wt_k) == P^(5)(1)/L (gold, 113 d)", d1 < mp.mpf('1e-110'), "|diff|=%s (gold precision)" % mp.nstr(d1,2))

    # ---------- exact Wt0, Wt1, Wt2 ----------
    from fractions import Fraction as F
    H = QS5(F(1,2)); T5 = QS5(F(2,5)); MT5 = QS5(F(-2,5)); F5 = QS5(F(1,5)); MF5 = QS5(F(-1,5))
    VAL = {(): {('1',0,0): ONE}, (ZERO,): {}, (P_LET,): {('1',0,1): ONE}, (Q_LET,): {('1',0,1): -ONE},
     (ZERO,ZERO): {}, (ZERO,ONE): {('z2',0,0): -ONE},
     (ZERO,P_LET): {('z2',0,0): T5, ('1',0,2): -H}, (ZERO,Q_LET): {('z2',0,0): MT5, ('1',0,2): ONE},
     (P_LET,ZERO): {('z2',0,0): MT5, ('1',0,2): H}, (Q_LET,ZERO): {('z2',0,0): T5, ('1',0,2): -ONE},
     (P_LET,P_LET): {('1',0,2): H}, (Q_LET,Q_LET): {('1',0,2): H},
     (P_LET,ONE): {('z2',0,0): MT5, ('1',0,2): ONE}, (Q_LET,ONE): {('z2',0,0): T5, ('1',0,2): -H},
     (P_LET,Q_LET): {('z2',0,0): MF5, ('1',0,2): ONE}, (Q_LET,P_LET): {('z2',0,0): F5, ('1',0,2): QS5(-2)}}
    def mul_ef(val, e, f):
        out = {}
        for (nm,a,b), c in val.items(): out[(nm,a+e,b+f)] = out.get((nm,a+e,b+f), ZERO) + c
        return {k:v for k,v in out.items() if v != ZERO}
    W012 = {0: {}, 1: {}, 2: {}}
    for mon, cu in parsed:
        a,b,c,d,e,f = mon
        for n,(p,q) in cu.items():
            cfq = QS5(p, q)
            for w, gam in rewrite[(n,a,b,c,d)].items():
                k = len(w) + e + f
                if k > 2: continue
                for key, c2 in mul_ef(VAL[w], e, f).items():
                    W012[k][key] = W012[k].get(key, ZERO) + cfq*gam*c2
    ok0 = all(v == ZERO for v in W012[0].values())
    e1 = {k: 1200*v for k, v in W012[1].items()}
    ok1 = (e1.get(('1',1,0), ZERO) == QS5(1200)) and (e1.get(('1',0,1), ZERO) == QS5(3000,-600)) and len(e1) == 2
    Wt1_exact = 1200*l2 + (3000-600*r5)*lf
    Wt2_exact = (-mp.mpf(1425)/16 - mp.mpf(5415)/16*r5)*z2 - 2400*l2**2 + (-12000+2400*r5)*l2*lf + (-mp.mpf(86925)/4 + mp.mpf(17325)/2*r5)*lf**2
    dWt1 = abs(Wt[1] - Wt1_exact); dWt2 = abs(Wt[2] - Wt2_exact)
    report("P2a Wt0 = 0 exactly (QS5)", ok0)
    report("P2b Wt1 = 1200 log2 + (3000-600sqrt5) logphi (exact)", ok1 and dWt1 < mp.mpf('1e-280'), "|num-exact|=%s" % mp.nstr(dWt1,2))
    report("P2c Wt2 = exact 4-term form (proved Li2)", dWt2 < mp.mpf('1e-280'), "|num-exact|=%s" % mp.nstr(dWt2,2))

    # ---------- closed-form evaluation helpers ----------
    phi = (1+r5)/2
    L3m = mp.polylog(3, -1/phi); L4m = mp.polylog(4, -1/phi); L4q = mp.polylog(4, 1/phi**2)
    L5m = mp.polylog(5, -1/phi); L5q = mp.polylog(5, 1/phi**2)
    K5m = L5m + lf**5/120; K5q = L5q - mp.mpf(2)/15*lf**5
    K4m = L4m - lf**4/24;  K4q = L4q + lf**4/3
    K3m = L3m + lf**3/6
    Lchi5 = mp.mpf(5)**-5*(mp.zeta(5, mp.mpf(1)/5) - mp.zeta(5, mp.mpf(2)/5) - mp.zeta(5, mp.mpf(3)/5) + mp.zeta(5, mp.mpf(4)/5))
    A = {'l2': l2, 'lf': lf, 'z2': z2, 'z3': z3, 'z5': z5, 'L3m': L3m, 'L4m': L4m, 'L4q': L4q,
         'L5m': L5m, 'L5q': L5q, 'K3m': K3m, 'K4m': K4m, 'K4q': K4q, 'K5m': K5m, 'K5q': K5q, 'Lchi5': Lchi5}
    import re
    def aval(nm):
        out = mp.mpf(1); pos = 0
        for m in re.finditer(r'(Lchi5|L3m|L4m|L4q|L5m|L5q|K3m|K4m|K4q|K5m|K5q|l2|lf|z2|z3|z5)(\^(\d+))?', nm):
            base, ex = m.group(1), m.group(3)
            out *= A[base]**(int(ex) if ex else 1)
        return out
    def ev(rel):
        return sum((mp.mpf(p)+mp.mpf(q)*r5)*aval(nm) for nm, (p, q) in rel.items())

    WT3_REL = {'L3m': ('5250','-1500'), 'l2^2lf': ('9000','-1800'), 'l2^3': ('1200','0'),
     'l2lf^2': ('64575/2','-25725/2'), 'lf^3': ('104225/2','-187325/8'), 'z2l2': ('300','795/2'),
     'z2lf': ('16275/8','-8535/8'), 'z3': ('3675/4','150')}
    d3 = abs(Wt[3] - ev(WT3_REL))
    report("P3 Wt3 == 8-atom closed form", d3 < mp.mpf('1e-280'), "|diff|=%s" % mp.nstr(d3,2))
    WT4_REL = {'L3ml2': ('-5250','1800'), 'L3mlf': ('-51375/2','23325/2'), 'L4m': ('-21750','10500'),
     'L4q': ('263625/32','-60975/16'), 'l2^2lf^2': ('-11550','4650'), 'l2^3lf': ('-2000','400'),
     'l2^4': ('-200','0'), 'l2lf^3': ('-176025/4','79075/4'), 'lf^4': ('-2395375/32','1076375/32'),
     'z2^2': ('-106965/16','54285/16'), 'z2l2^2': ('-75','-135'), 'z2l2lf': ('-2550','1185'),
     'z2lf^2': ('-78525/16','30255/16'), 'z3l2': ('-375','-90'), 'z3lf': ('-4200','2625/2')}
    d4 = abs(Wt[4] - ev(WT4_REL))
    report("P4 Wt4 == 15-atom closed form", d4 < mp.mpf('1e-280'), "|diff|=%s" % mp.nstr(d4,2))
    WT5_LOG = {'l2^5': ('10','0'), 'l2^4lf': ('125','-25'), 'l2^3lf^2': ('2025/2','-825/2'),
     'l2^2lf^3': ('52075/8','-23525/8'), 'l2lf^4': ('449025/16','-403175/32'),
     'lf^5': ('15112175/256','-6759335/256')}
    WT5_R = {'K3ml2^2': ('750','-525/2'), 'K3ml2lf': ('7875','-3675'), 'K3mlf^2': ('183375/8','-10425'),
     'K4ml2': ('14625/2','-6975/2'), 'K4mlf': ('304875/8','-134775/8'), 'K4ql2': ('-209625/64','48975/32'),
     'K4qlf': ('-13847625/512','6241425/512'), 'K5m': ('246375/4','-214875/8'),
     'K5q': ('-18716625/1024','8328825/1024'), 'Lchi5': ('-172009375/4096','80809375/4096'),
     'z2L3m': ('-168075/32','71685/32'), 'z2^2l2': ('150135/64','-74955/64'), 'z2^2lf': ('-17145/4','13251/8'),
     'z2l2^2lf': ('1575/4','-705/4'), 'z2l2^3': ('25/4','45/4'), 'z2l2lf^2': ('16875/8','-13935/16'),
     'z2lf^3': ('528925/64','-233505/64'), 'z2z3': ('29865/16','-12879/16'), 'z3l2^2': ('225/8','75/8'),
     'z3l2lf': ('6975/8','-255'), 'z3lf^2': ('-296475/128','142995/128'), 'z5': ('4209075/256','-1914795/256')}
    d5 = abs(Wt[5] - ev(WT5_LOG) - ev(WT5_R))
    report("P5 Wt5 == exact log part + (21 classical + Lchi5) remainder", d5 < mp.mpf('1e-280'), "|diff|=%s" % mp.nstr(d5,2))

    # ---------- P6: local law prop:order5 at m=5 (exact) ----------
    import sympy as sp
    sr5 = sp.sqrt(5)
    s = (sr5-1)/4; m = 5
    A5 = sp.expand(2880*m*s**4*(65*s**6 - 63*s**4 + 17*s**2 - 1))
    A23 = sp.expand(-960*m*s**4*(45*s**6 - 27*s**4 + 6*s**2 - 1))
    cLoL = sp.Rational(5,2) + sr5   # cos(pi/10)/L, proved by squares (both positive)
    def pq(e):
        e = sp.expand(e)
        return (sp.expand((e + e.subs(sr5,-sr5))/2), sp.expand((e - e.subs(sr5,-sr5))/(2*sr5)))
    p5_, q5_ = pq(sp.expand(A5*cLoL)); p23_, q23_ = pq(sp.expand(A23*cLoL))
    ok6 = (p5_, q5_) == (sp.Rational(73125,64), sp.Rational(-30825,64)) and \
          (p23_, q23_) == (sp.Rational(-64875,64), sp.Rational(32775,64))
    report("P6 prop:order5 m=5: C5^(5)/L, C5^(2,3)/L exact", ok6)

    # ---------- P7: L(chi_5,5) Euler product ----------
    mp.mp.dps = 60
    def leg(n):
        r = n % 5
        return 1 if r in (1,4) else (-1 if r in (2,3) else 0)
    def primes_upto(n):
        sieve = [True]*(n+1); sieve[0]=sieve[1]=False
        for i in range(2,int(n**0.5)+1):
            if sieve[i]:
                for j in range(i*i,n+1,i): sieve[j]=False
        return [p for p in range(2,n+1) if sieve[p]]
    prod = mp.mpf(1)
    for p in primes_upto(200000):
        prod *= 1/(1 - leg(p)*mp.mpf(p)**-5)
    d7 = abs(prod - Lchi5)
    report("P7 L(chi_5,5) Hurwitz vs Euler product", d7 < mp.mpf('1e-24'), "|diff|=%s" % mp.nstr(d7,2))

    print()
    print("ALL %d CHECKS: %s  (%.1f s)" % (len(PASS), "PASS" if all(PASS) else "FAIL", time.time()-t0))
    return 0 if all(PASS) else 1

if __name__ == "__main__":
    sys.exit(main())
'''
open('verify_dec5e.py','w').write(ver)
import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run(['nice', PY, 'verify_dec5e.py'], capture_output=True, text=True, cwd='.', timeout=600)
print(r.stdout)
if r.stderr: print("STDERR:", r.stderr[-600:])
print("exit:", r.returncode, "wall %.0fs" % (time.time()-t0))
