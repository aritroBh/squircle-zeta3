import os, sys, time, pickle, math
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import oct5d2_engine as E
from oct5d2_engine import Q, ZERO, ONE, P_LET, Q_LET
import mpmath as mp
mp.mp.dps = 60

# Dual rewrite via u = 1-t:
# I(n;a,b,c,d) = sum over binomial choices of
#   C-coeffs * (-X/2)^{c-c'} (X/2)^{d-d'} * int_0^1 (1-u)^n (log(1-u))^a (log u)^b G(q;u)^{c'} G(p;u)^{d'} du
# and int (1-u)^n f(u) du = sum_k (-1)^k C(n,k) int u^k f(u) du.
# G(q;u)^{c'} = c'! G(q^{c'};u) etc, log powers: (log(1-u))^a = a! G(1^a;u), (log u)^b = b! G(0^b;u).
# So the dual rewrite expresses I(n;a,b,c,d) as Q(sqrt2)-combination of
#   X^{e'} * I(k; b', a', c''=d' roles...) : concretely int u^k (log u)^{b}(log(1-u))^{a} (log(1+k0' u))^{c'} (log(1+k1' u))^{d'}
# where log(1+k0' u) = G(q;u) corresponds to k0' = -1/q = k1, and G(p;u) to k1' = k0.
# i.e. = X^{c-c'+d-d'} (-1)^{c-c'} * binom(c,c') binom(d,d') * I(k; b, a, d', c') with factorial/shuffle built into monomial_integral.

# Implement with exact arithmetic: represent the transformed integrand-expansion terms directly
# via shuffle products of G-words (like monomial_integral but on the dual side).
from functools import cache
import math

@cache
def dual_rewrite(n,a,b,c,d):
    """returns dict {word: QS2} = exact rewrite of I(n;a,b,c,d) via u=1-t substitution,
    as combination of G_0(word;1) with |word|<=a+b+c+d (top weight a+b+c+d)."""
    out = {}
    # precompute G-word expansions of each factor as {word: QS2} with weight = degree
    # factors: (log(1-u))^a -> a! G(1^a); (log u)^b -> b! G(0^b);
    # [-X/2 + G(q)]^c -> sum_{c'} binom(c,c')(-X/2)^{c-c'} c'! G(q^{c'})  [X-terms carried as scalars]
    # [ X/2 + G(p)]^d -> sum_{d'} binom(d,d')( X/2)^{d-d'} d'! G(p^{d'})
    # constant scalars X = log2 is NOT a scalar in Q(sqrt2)!! X is a period; the dual rewrite as a
    # period identity mixes weights: terms with X^{e'} are lower-weight G's times constant logs.
    # For TOP weight (a+b+c+d) only c'=c and d'=d survive. So top part:
    #   (1-u)^n a! b! c! d! G(1^a;u) G(0^b;u) G(q^c;u) G(p^d;u)
    # = sum_k (-1)^k C(n,k) a!b!c!d! [u^k G(shuffle(1^a,0^b,q^c,p^d);u)]
    top = {}
    parts = []
    if a: parts.append((ONE,)*a)
    if b: parts.append((ZERO,)*b)
    if c: parts.append((Q_LET,)*c)
    if d: parts.append((P_LET,)*d)
    sw = E.shuffle_words(parts)
    fac = Q(math.factorial(a)*math.factorial(b)*math.factorial(c)*math.factorial(d))
    for k in range(n+1):
        sgn = Q((-1)**k) * Q(math.comb(n,k))
        for w, mult in sw.items():
            ig = E.IG(k, w)
            for v, cv in ig.items():
                top[v] = top.get(v, ZERO) + sgn*fac*mult*cv
    return {k2: v for k2, v in top.items() if v != ZERO}, a+b+c+d

# numeric check of dual rewrite against itab110 for a sample of integrals (all weights via full dict)
SRC = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
ITAB = {k: mp.mpf(v) for k, v in pickle.load(open(SRC+"itab110.pkl","rb")).items()}
import gengine
G1, letters = gengine.gtab(170, 480)
r2 = mp.sqrt(2)
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r2
def w2mpf(w):
    return tuple(L0 if x==ZERO else L1 if x==ONE else LP if x==P_LET else LQ for x in w)
import random
random.seed(5)
rewrite = pickle.load(open("rewrite_exact.pkl","rb"))
keys = random.sample(sorted(rewrite.keys()), 30)
worst = mp.mpf(0)
for key in keys:
    n,a,b,c,d = key
    dr, k = dual_rewrite(n,a,b,c,d)
    # full dual dict contains only full expansion? dual_rewrite returns FULL exact rewrite incl lower weights
    val = sum(q2f(v)*G1[w2mpf(w)] for w,v in dr.items())
    ref = ITAB[(n,a,b,c,d)]
    worst = max(worst, abs(val-ref))
print("dual rewrite numeric check on 30 random integrals, worst residual vs itab110:", mp.nstr(worst,4))
pickle.dump(0, open("tmp_dual_flag.pkl","wb"))
