
import pickle
import dec5d_lib as D
import mpmath as mpm
mp = mpm.mp
mp.dps = 120
QS5 = D.QS5; ZERO, ONE = QS5(0), QS5(1)
P_LET = -D.PHI; Q_LET = D.PHI**2
import gengine_dec
G1, letters = gengine_dec.gtab(120, 380, cachedir=".")
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
r5 = mp.sqrt(5)
def q5f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r5
def w2mpf(w): return tuple(L0 if x==ZERO else L1 if x==ONE else LP if x==P_LET else LQ for x in w)
rewrite = pickle.load(open("rewrite_dec_exact.pkl","rb"))
ival = {k: sum(q5f(v)*G1[w2mpf(w)] for w, v in r.items()) for k, r in rewrite.items()}

parsed = pickle.load(open("parsed_mons_dec.pkl","rb"))
L2v = mp.log(2); LFv = mp.log((1+r5)/2)
acc = mp.mpf(0); acc5 = mp.mpf(0)
contrib = {'local': mp.mpf(0),'single-1': mp.mpf(0),'single-4': mp.mpf(0),'mixed': mp.mpf(0)}
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else (('single-1' if c>0 else 'single-4') if (c>0 or d>0) else 'local')
    constfac = L2v**e * LFv**f
    for n,(p,q) in cu.items():
        term = (mp.mpf(str(p)) + mp.mpf(str(q))*r5)*constfac*ival[(n,a,b,c,d)]
        acc += term
        if sum(mon)==5:
            acc5 += term; contrib[cls] += term
gold = pickle.load(open("p5_gold.pkl","rb"))
P5L_gold = mp.mpf(gold['P5overL'])
print("P5/L from exact rewrite =", mp.nstr(1200*acc, 60))
print("P5/L engine             =", mp.nstr(P5L_gold, 60))
print("diff:", mp.nstr(abs(1200*acc - P5L_gold), 3))
prev = pickle.load(open("w5L_num.pkl","rb"))
print("\nW5/L from rewrite =", mp.nstr(1200*acc5, 55))
print("W5/L from itab    =", mp.nstr(mp.mpf(prev['W5L_110']), 55))
print("diff:", mp.nstr(abs(1200*acc5 - mp.mpf(prev['W5L_110'])), 3))
for k in contrib:
    print("  block %-9s rewrite=%s  itab=%s  diff=%s" % (k, mp.nstr(1200*contrib[k], 40),
          mp.nstr(mp.mpf(prev['blocks'][k]), 40), mp.nstr(abs(1200*contrib[k]-mp.mpf(prev['blocks'][k])), 3)))
