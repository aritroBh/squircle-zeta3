import os, sys, time, pickle, re
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
import mpmath as mp
mp.mp.dps = 430
import oct5d2_engine as E
from oct5d2_engine import Q, ZERO, ONE, P_LET, Q_LET
import gengine
from fractions import Fraction as F
from flint import fmpz_mat

SRC = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
parsed = pickle.load(open(SRC+"parsed_mons.pkl","rb"))
rewrite = pickle.load(open("rewrite_exact.pkl","rb"))
G1, letters = gengine.gtab(430, 1440)
r2 = mp.sqrt(2); L2C = mp.log(2); LU = mp.log(1+r2)
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r2
def w2mpf(w): return tuple(L0 if x==ZERO else L1 if x==ONE else LP if x==P_LET else LQ for x in w)
Gcache = {}
def gval(w):
    if w not in Gcache: Gcache[w] = G1[w2mpf(w)]
    return Gcache[w]
t0=time.time()
acc = mp.mpf(0)
for mon, cu in parsed:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    cf = L2C**e * LU**f
    for n,(p,q) in cu.items():
        iv = Gcache and None
        iv = sum(q2f(v)*gval(w) for w, v in rewrite[(n,a,b,c,d)].items())
        acc += (mp.mpf(str(p)) + mp.mpf(str(q))*r2)*cf*iv
W5 = 960*acc
print("W5 assembled at dps 430 in %.1fs" % (time.time()-t0))
print(mp.nstr(W5, 60))
pickle.dump({'W5': mp.nstr(W5, 420)}, open("w5_430.pkl","wb"))

# atom values at dps 430
z5 = mp.zeta(5); z3 = mp.zeta(3); z2 = mp.pi**2/6
X = L2C; Y = LU; b = r2-1
Li5v = {'b': mp.polylog(5,b), '-b': mp.polylog(5,-b), '1/2': mp.polylog(5,mp.mpf('0.5')),
        '1/sqrt2': mp.polylog(5,1/r2), '2-sqrt2': mp.polylog(5,2-r2), 'p': mp.re(mp.polylog(5, 2+r2))}
Li4v = {'b': mp.polylog(4,b), '-b': mp.polylog(4,-b), '1/2': mp.polylog(4,mp.mpf('0.5')),
        '1/sqrt2': mp.polylog(4,1/r2), '2-sqrt2': mp.polylog(4,2-r2), 'p': mp.re(mp.polylog(4, 2+r2))}
Li3v = {'b': mp.polylog(3,b), '1/2': mp.polylog(3,mp.mpf('0.5')), '1/sqrt2': mp.polylog(3,1/r2)}
Li2v = {'b': mp.polylog(2,b)}
def Li32(x, y, N):
    tot = mp.mpf(0); xp = x
    for n in range(2, N+1):
        if n > 2: xp *= x
        inner = mp.mpf(0); yp = y
        for m in range(1, n):
            if m > 1: yp *= y
            inner += yp/(m*m)
        tot += xp/(n**3)*inner
    return tot
t0=time.time()
L32_b1 = Li32(b, mp.mpf(1), 1600); L32_h1 = Li32(1/r2, mp.mpf(1), 3000)
print("Li32 at dps430 in %.1fs" % (time.time()-t0))
z32 = G1[(L0,L0,L1,L0,L1)]
K = (-8)*L32_b1 - 32*Li32(-b, mp.mpf(1), 1600) - 16*Li32(b, -mp.mpf(1), 1600) - 16*Li32(-b, -mp.mpf(1), 1600) + Li32(-b*b, mp.mpf(1), 1600)
pickle.dump({'Li5v':{k:mp.nstr(v,420) for k,v in Li5v.items()},
             'Li4v':{k:mp.nstr(v,420) for k,v in Li4v.items()},
             'Li3v':{k:mp.nstr(v,420) for k,v in Li3v.items()},
             'Li2v':{k:mp.nstr(v,420) for k,v in Li2v.items()},
             'L32_b1': mp.nstr(L32_b1,420), 'L32_h1': mp.nstr(L32_h1,420),
             'z32': mp.nstr(z32,420), 'K': mp.nstr(K,420)}, open("atom_vals430.pkl","wb"))
print("saved atom_vals430")
