"""symengine.py — full number-level Goncharov symbols for G(w;1), letters {0,1,p,q},
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
