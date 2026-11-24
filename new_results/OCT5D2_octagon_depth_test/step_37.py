import os, sys, time, pickle, math, random
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import oct5d2_engine as E
from oct5d2_engine import Q, ZERO, ONE, P_LET, Q_LET
import mpmath as mp
mp.mp.dps = 80
from functools import cache

@cache
def dual_full(n,a,b,c,d):
    """I(n;a,b,c,d) via u=1-t: dict {(word, eX): QS2} = sum (log2)^eX * G(word;1)."""
    out = {}
    for cp in range(c+1):
        for dp in range(d+1):
            binom = Q(math.comb(c,cp)*math.comb(d,dp))
            eX = (c-cp)+(d-dp)
            scal = binom * Q((-1)**(c-cp)) / Q(2**(eX))
            parts = []
            if a: parts.append((ONE,)*a)
            if b: parts.append((ZERO,)*b)
            if cp: parts.append((Q_LET,)*cp)
            if dp: parts.append((P_LET,)*dp)
            sw = E.shuffle_words(parts)
            fac = scal * Q(math.factorial(a)*math.factorial(b)*math.factorial(cp)*math.factorial(dp))
            for k in range(n+1):
                sgn = Q((-1)**k) * Q(math.comb(n,k))
                for w, mult in sw.items():
                    for v, cv in E.IG(k, w).items():
                        key = (v, eX)
                        out[key] = out.get(key, ZERO) + sgn*fac*mult*cv
    return {k2: v for k2, v in out.items() if v != ZERO}

# numeric check on random integrals against itab110
SRC = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
ITAB = {k: mp.mpf(v) for k, v in pickle.load(open(SRC+"itab110.pkl","rb")).items()}
import gengine
G1, letters = gengine.gtab(170, 480)
r2 = mp.sqrt(2); L2C = mp.log(2)
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r2
def w2mpf(w): return tuple(L0 if x==ZERO else L1 if x==ONE else LP if x==P_LET else LQ for x in w)
random.seed(5)
rewrite = pickle.load(open("rewrite_exact.pkl","rb"))
keys = random.sample(sorted(rewrite.keys()), 40)
worst = mp.mpf(0); wk=None
for key in keys:
    n,a,b,c,d = key
    val = sum(q2f(v)*(L2C**eX)*G1[w2mpf(w)] for (w,eX),v in dual_full(n,a,b,c,d).items())
    ref = ITAB[(n,a,b,c,d)]
    if abs(val-ref)>worst: worst, wk = abs(val-ref), key
print("dual full rewrite vs itab110, 40 random integrals, worst residual:", mp.nstr(worst,4), wk)
