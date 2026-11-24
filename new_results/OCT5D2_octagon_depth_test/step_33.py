import os, sys, pickle, itertools
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import symengine as S
from symengine import Q, ZERO, ONE, P_LET, Q_LET, add, scale
import w4_symbol as W

def scale_word(w, lam): return tuple(x*lam for x in w)
bad=0; tested=0
for wt in range(2,6):
    for w in itertools.product([ZERO, ONE, P_LET], repeat=wt):
        if w[0] == ONE or w[-1] == ZERO: continue   # convergent only
        s1 = W.g_symbol(ZERO, w, ONE)
        s2 = W.g_symbol(ZERO, scale_word(w, ONE/P_LET), ONE/P_LET)
        tested += 1
        if add(s1, scale(s2,-1)) != {}:
            bad += 1
            if bad <= 3: print("MISMATCH", w, s1, s2)
print(f"p-scaling, convergent words: {tested} tested, {bad} mismatches")

lam = -ONE/Q_LET
bad2=0; tested2=0
for wt in range(2,6):
    for w in itertools.product([ZERO, ONE, Q_LET], repeat=wt):
        if w[0] == ONE or w[-1] == ZERO: continue
        s1 = W.g_symbol(ZERO, w, ONE)
        s2 = W.g_symbol(ZERO, scale_word(w, lam), lam)
        tested2 += 1
        if add(s1, scale(s2,-1)) != {}:
            bad2 += 1
            if bad2 <= 3: print("MISMATCH(q)", w)
print(f"q-scaling(lambda=-1/q>0), convergent words: {tested2} tested, {bad2} mismatches")

# and a genuinely mixed-letters scaling check: letters {0,1,p,q}, scale by 1/p (rationalizes p)
bad3=0; tested3=0
for wt in range(2,6):
    for w in itertools.product([ZERO, ONE, P_LET, Q_LET], repeat=wt):
        if w[0] == ONE or w[-1] == ZERO: continue
        if P_LET not in w or Q_LET not in w: continue
        s1 = W.g_symbol(ZERO, w, ONE)
        s2 = W.g_symbol(ZERO, scale_word(w, ONE/P_LET), ONE/P_LET)
        tested3 += 1
        if add(s1, scale(s2,-1)) != {}:
            bad3 += 1
            if bad3 <= 3: print("MISMATCH(mixed)", w)
print(f"mixed-letters scaling: {tested3} tested, {bad3} mismatches")
