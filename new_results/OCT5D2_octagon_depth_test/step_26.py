import os, sys
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")

symengine = r'''"""symengine.py — full number-level Goncharov symbols for G(w;1), letters {0,1,p,q},
over the S-unit log alphabet {X=log2 (slot 0), Y=log(1+sqrt2) (slot 1)}.

Uses the validated constant-retaining g_symbol from A3_W4_proof/w4_symbol.py
(Goncharov's formula with EVERY entry (a_i - a_j), including constant ones,
factored exactly via logvec over Q(sqrt2), modulo torsion).
"""
import sys
from pathlib import Path
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
from qsqrt2 import QS2 as Q, SQ2 as R, A as AA, B as BB
import w4_symbol as W

ZERO, ONE = Q(0), Q(1)
P_LET, Q_LET = Q(2,1), Q(-1,-1)

add, scale, clean, shuffle, concat = W.add, W.scale, W.clean, W.shuffle, W.concat
g_symbol, li_symbol, log_symbol = W.g_symbol, W.li_symbol, W.log_symbol

X = {(0,): ONE}   # log 2
Y = {(1,): ONE}   # log(1+sqrt2)

def bracket(x, y):
    return add(concat(x, y), scale(concat(y, x), -1))

L2 = bracket(X, Y)                       # [X,Y]
L3 = [bracket(X, L2), bracket(Y, L2)]    # [X,[X,Y]], [Y,[X,Y]]

def obstruction(s):
    """the two antisymmetrised (3,2) pairings, exactly as M4c w5_symbol.py."""
    return [sum((c*d*(s.get(u+v,ZERO)-s.get(v+u,ZERO))
                 for u,c in p.items() for v,d in L2.items()),ZERO) for p in L3]

def const_pow(const_letter, k):
    """symbol of (log const)^k : k! * constant word."""
    import math
    return {(const_letter,)*k: Q(math.factorial(k))}

def gsym(w):
    """full number-level symbol of G_0(w; 1)."""
    return g_symbol(ZERO, w, ONE)
'''
open("symengine.py","w").write(symengine)

# Test the machinery on our letters
import symengine as S
importlib = __import__('importlib')
print("L2 =", S.L2)
print("L3 =", S.L3)

# test word 1: G(0,p;1) = -Li_2(1/p) : symbol should be -li_symbol(2, 1/p)
from qsqrt2 import QS2 as Q
ONE_ = Q(1); P = Q(2,1); Qq = Q(-1,-1); Z = Q(0)
s1 = S.gsym((Z, P))
ref = S.scale(S.li_symbol(2, ONE_/P), -1)
print("G(0,p;1) symbol:", s1)
print("-li_symbol(2,1/p):", ref)
print("match:", S.add(s1, S.scale(ref,-1)) == {})

# test word 2 (mixed): shuffle identity on symbols: S G(p) ш S G(q) = S G(p,q)+S G(q,p)
lhs = S.shuffle(S.gsym((P,)), S.gsym((Qq,)))
rhs = S.add(S.gsym((P,Qq)), S.gsym((Qq,P)))
print("shuffle symbol identity (p,q):", S.add(lhs, S.scale(rhs,-1)) == {})

# test G(0^{k-1}, p;1) = -Li_k(1/p) at symbol level for k=3,4,5
for k in range(3,6):
    w = (Z,)*(k-1)+(P,)
    s = S.gsym(w); ref = S.scale(S.li_symbol(k, ONE_/P), -1)
    print(f"G(0^{k-1},p;1) vs -Li_{k}(1/p) symbol match:", S.add(s, S.scale(ref,-1)) == {})

# integrability check on a weight-5 mixed word symbol: positions 1..4, x_i wedge x_{i+1}
def integrability_residual(s):
    # words over {0,1} length 5; for each adjacent pair compute det [[a,b],[c,d]]
    from collections import defaultdict
    res = defaultdict(lambda: Z)
    for w, c in s.items():
        for i in range(4):
            a,b2 = w[i], w[i+1]
            if a != b2:
                det = (1 if (a,b2)==(0,1) else -1)
                res[i] += c*det
    return dict(res)
s5 = S.gsym((P, Qq, Z, ONE_, P))
print("integrability residuals of S G(p,q,0,1,p;1):", integrability_residual(s5))
print("symbol has", len(s5), "words")
