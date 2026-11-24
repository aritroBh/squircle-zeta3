import os, sys, time, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import oct5d2_engine as E
from oct5d2_engine import Q, ZERO, ONE, P_LET, Q_LET
import symengine as S
from symengine import add, scale, gsym
import math, random
from functools import cache

@cache
def dual_full(n,a,b,c,d):
    out = {}
    for cp in range(c+1):
        for dp in range(d+1):
            binom = Q(math.comb(c,cp)*math.comb(d,dp))
            eX = (c-cp)+(d-dp)
            scal = binom * Q((-1)**(c-cp)) / Q(2**eX)
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

rewrite = pickle.load(open("rewrite_exact.pkl","rb"))
symcache = {}
def Gsym(w):
    if w not in symcache: symcache[w] = gsym(w)
    return symcache[w]

bad = 0; tested = 0
t0=time.time()
for (n,a,b,c,d), r in rewrite.items():
    k = a+b+c+d
    sdir = {}
    for w, gam in r.items():
        if len(w) == k: sdir = add(sdir, scale(Gsym(w), gam))
    sdul = {}
    for (v, eX), gam in dual_full(n,a,b,c,d).items():
        if len(v) == k and eX == 0: sdul = add(sdul, scale(Gsym(v), gam))
    tested += 1
    if add(sdir, scale(sdul,-1)) != {}:
        bad += 1
        if bad <= 3:
            print("TOP-SYMBOL MISMATCH", (n,a,b,c,d))
            print("  direct:", sdir)
            print("  dual  :", sdul)
print("top-symbol direct-vs-dual: %d integrals tested in %.1fs, mismatches: %d" % (tested, time.time()-t0, bad))
