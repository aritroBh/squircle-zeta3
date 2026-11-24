
import pickle
import dec5d_lib as D
import mpmath as mp
mp.mp.dps = 420
QS5 = D.QS5; ZERO, ONE = QS5(0), QS5(1)
P_LET = -D.PHI; Q_LET = D.PHI**2
import gengine_dec
G1, letters = gengine_dec.gtab(420, 1180, cachedir=".")
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
r5 = mp.sqrt(5); phi = (1+r5)/2
def q5f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r5
def w2mpf(w): return tuple(L0 if x==ZERO else L1 if x==ONE else LP if x==P_LET else LQ for x in w)
rewrite = pickle.load(open("rewrite_dec_exact.pkl","rb"))
ival = {k: sum(q5f(v)*G1[w2mpf(w)] for w, v in r.items()) for k, r in rewrite.items()}
parsed = pickle.load(open("parsed_mons_dec.pkl","rb"))
L2v = mp.log(2); LFv = mp.log(phi)
acc5 = mp.mpf(0)
blocks = {'local': mp.mpf(0),'single-1': mp.mpf(0),'single-4': mp.mpf(0),'mixed': mp.mpf(0)}
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    if sum(mon)!=5: continue
    cls = 'mixed' if (c>0 and d>0) else (('single-1' if c>0 else 'single-4') if (c>0 or d>0) else 'local')
    constfac = L2v**e * LFv**f
    for n,(p,q) in cu.items():
        term = (mp.mpf(str(p)) + mp.mpf(str(q))*r5)*constfac*ival[(n,a,b,c,d)]
        acc5 += term; blocks[cls] += term
W5L420 = 1200*acc5
prev = mp.mpf(pickle.load(open("w5L_250.pkl","rb"))['W5L'])
print("W5/L@420 vs @250:", mp.nstr(abs(W5L420-prev),3))
print("W5/L =", mp.nstr(W5L420, 80))
pickle.dump({'W5L': mp.nstr(W5L420, 400), 'blocks': {k: mp.nstr(1200*v,400) for k,v in blocks.items()}}, open("w5L_420.pkl","wb"))
print("saved w5L_420.pkl")
