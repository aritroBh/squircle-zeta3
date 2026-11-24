
import pickle, time
import dec5d_lib as D
import mpmath as mpm
mp = mpm.mp
mp.dps = 110

QS5 = D.QS5
P_LET = -D.PHI; Q_LET = D.PHI**2
RW = D.make_rewrite(QS5, P_LET, Q_LET)
parsed = pickle.load(open("parsed_mons_dec.pkl","rb"))
need = {}
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    for n in cu: need[(n,a,b,c,d)] = mon
t0 = time.time()
rewrite = {k: RW['monomial_integral'](*k) for k in sorted(need)}
print("rewrote %d integrals in %.1fs" % (len(rewrite), time.time()-t0))
nwords = set()
for r in rewrite.values(): nwords.update(r.keys())
print("distinct words: %d, max weight %d" % (len(nwords), max(map(len,nwords))))
bad1 = [w for r in rewrite.values() for w,v in r.items() if w and w[0]==QS5(1) and v!=QS5(0)]
print("words starting with 1 having nonzero coeff:", len(bad1))
pickle.dump(rewrite, open("rewrite_dec_exact.pkl","wb"))

# numerical check against itab_dec110
import gengine_dec
G1, letters = gengine_dec.gtab(120, 380, cachedir=".")
mp.dps = 120
r5 = mp.sqrt(5)
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
def q5f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r5
ZERO, ONE = QS5(0), QS5(1)
def w2mpf(w): return tuple(L0 if x==ZERO else L1 if x==ONE else LP if x==P_LET else LQ for x in w)
ITAB = {k: mp.mpf(v) for k, v in pickle.load(open("itab_dec110.pkl","rb")).items()}
worst = mp.mpf(0); worstk = None
for k, r in rewrite.items():
    val = sum(q5f(v)*G1[w2mpf(w)] for w, v in r.items())
    d = abs(val - ITAB[k])
    if d > worst: worst, worstk = d, k
print("worst |rewrite - itab| over 756 integrals:", mp.nstr(worst, 3), "at", worstk)
