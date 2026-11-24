import os, sys, pickle, json
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import symengine as S
from symengine import Q, ZERO, ONE, P_LET, Q_LET, add, scale, shuffle, gsym, obstruction, li_symbol
import w4_symbol as W
from fractions import Fraction as F

# --- check 2: scaling invariance of g_symbol on single-ratio words
# G_0(w;1) with letters {0,1,p} equals G_0(w/p; 1/p); symbols computed both ways must agree.
def scale_word(w, lam):
    return tuple(x*lam for x in w)
bad = 0; tested = 0
import itertools, random
random.seed(3)
letters = [ZERO, ONE, P_LET]
for wt in range(2,6):
    for w in itertools.product(letters, repeat=wt):
        if w[0] == ONE: continue
        s1 = W.g_symbol(ZERO, w, ONE)
        s2 = W.g_symbol(ZERO, scale_word(w, ONE/P_LET), ONE/P_LET)
        tested += 1
        if add(s1, scale(s2,-1)) != {}:
            bad += 1
            if bad <= 3: print("MISMATCH", w)
print(f"scaling invariance: {tested} words tested, mismatches: {bad}")

# same for q-words: scale by 1/q (NEGATIVE lambda: G(q w'; q) = G(w';1)? scaling by negative reverses... skip q-scaling,
# instead scale q-words by -1/q? t -> s/q with q<0 flips orientation; use lambda=-1/q>0: G(w;1)=G(-w/q; -1/q)?
# G(lambda w; lambda) = G(w;1) holds for real lambda != 0 by s = t/lambda (orientation preserved iff lambda>0).
# For q<0 use lambda = -1/q = 1/(1+sqrt2) > 0: letters {-0, -1, -q} * lambda ... different letter set {0,-lambda,-q*lambda}. Test consistency.
lam = -ONE/Q_LET
bad2 = 0; tested2 = 0
for wt in range(2,6):
    for w in itertools.product([ZERO, ONE, Q_LET], repeat=wt):
        if w[0] == ONE: continue
        s1 = W.g_symbol(ZERO, w, ONE)
        s2 = W.g_symbol(ZERO, scale_word(w, lam), lam)
        tested2 += 1
        if add(s1, scale(s2,-1)) != {}:
            bad2 += 1
            if bad2 <= 3: print("MISMATCH(q)", w)
print(f"q-scaling invariance: {tested2} words tested, mismatches: {bad2}")
