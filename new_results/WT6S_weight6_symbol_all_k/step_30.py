
import os
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
core = r'''"""WT6S core: weight-6 symbols of the squircle W6(k) pieces from degree-5 kernels.

Machinery adapted from M5a w6_symbol.py (same normalisation):
  letters 0 = X = log2, 1 = Y = log(1+sqrt2);  P = [X,Y];
  B3 = {[X,P],[Y,P]}, B4 = {[X,[X,P]],[Y,[X,P]],[Y,[Y,P]]}.
  pair33(S) = pairing against [X,P] wedge [Y,P];  pair42(S) = pairings against B4 wedge P.
N_{(i,j),r} = int_b^a log^i t log^j(1-t) dv/(v-r), t = (v^2+2v-1)/(4v).
"""
import json, pickle
from fractions import Fraction as F
from pathlib import Path
import sys
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import qsqrt2  # local copy, cached before w4_symbol imports it
from w4_symbol import (Q, R, A, B, ZERO, ONE, g_symbol, log_product, LT, MT,
                       shuffle, log_symbol, add, scale, concat, clean,
                       ban_float, self_check, li_symbol, shuffles)

ZEROQ, ONEQ = Q(0), Q(1)
POLE = {"0": ZEROQ, "a": A, "-a": -A, "b": B, "-b": -B}
POLE_PAIRS = {"0": (0, 0), "a": (1, 1), "-a": (-1, -1), "b": (-1, 1), "-b": (1, -1)}
PAIR_TO_POLE = {v: k for k, v in POLE_PAIRS.items()}

def log(s):
    print(s, flush=True)

# ---------------------------------------------------------------- Lie words / pairings
def bracket(x, y):
    return add(concat(x, y), scale(concat(y, x), -1))

X, Y = {(0,): ONEQ}, {(1,): ONEQ}
P = bracket(X, Y)
B3 = [bracket(X, P), bracket(Y, P)]                                  # [X,P], [Y,P]
B4 = [bracket(X, B3[0]), bracket(Y, B3[0]), bracket(Y, B3[1])]       # len-4 words
L2 = P                                                                # weight-5 control

def pair_ab(S, A_lie, B_lie):
    la = len(next(iter(A_lie))); lb = len(next(iter(B_lie)))
    tot = ZEROQ
    for w, c in S.items():
        v = A_lie.get(w[:la], ZEROQ) * B_lie.get(w[la:], ZEROQ) \
            - B_lie.get(w[:lb], ZEROQ) * A_lie.get(w[lb:], ZEROQ)
        if v:
            tot += c * v
    return tot

def pair33(S):
    return pair_ab(S, B3[0], B3[1])

def pair42(S):
    return [pair_ab(S, A4, P) for A4 in B4]

def pairing_vector6(S):
    return [pair33(S)] + pair42(S)

# ---------------------------------------------------------------- symbols of N_{(i,j),r}
def symbol_of_N(i, j, rname):
    """Weight-(i+j+1) symbol of N_{(i,j),r} (one unit of kernel 1/(v-r))."""
    r = POLE[rname]
    factors = (LT,) * i + (MT,) * j
    out = {}
    for (w, k), d in log_product(*factors).items():
        word = (r,) + w
        s = add(g_symbol(ONE, word, A), scale(g_symbol(ONE, word, B), -1))
        for _ in range(k):
            s = shuffle(s, log_symbol(Q(2)))
        out = add(out, scale(s, d))
    return clean(out)

def symbol_of_kernel_terms(terms, i, j):
    """terms = [(polename, (p, q))]; coefficient p + q sqrt2 times 1/(v-pole)."""
    out = {}
    for name, (p, q) in terms:
        c = Q(F(p), F(q))
        if c == ZEROQ:
            continue
        out = add(out, scale(symbol_of_N(i, j, name), c))
    return clean(out)

# ---------------------------------------------------------------- depth-two atoms
def depth2_symbol(m, n, x, y):
    """Symbol of Li_{m,n}(x,y) = G_0(0^{m-1}, 1/x, 0^{n-1}, 1/(x y); 1) (depth 2)."""
    word = (ZEROQ,) * (m - 1) + (ONEQ / x,) + (ZEROQ,) * (n - 1) + (ONEQ / (x * y),)
    return g_symbol(ZEROQ, word, ONEQ)

# ---------------------------------------------------------------- kernel JSON loading
def load_deg5():
    return json.loads((HERE / 'deg5_residue_polynomials.json').read_text())

def load_deg4():
    return json.loads((HERE / 'deg4_residue_polynomials.json').read_text())

KEYS5 = ["50", "41", "32"]   # i>=j, i+j=5
KEYS4 = ["40", "31", "22"]   # i>=j, i+j=4

if __name__ == '__main__':
    ban_float(); self_check()
    # quick control: W6 pairings at k=5 from omega5.json (== kern_5.pkl)
    om = json.loads((HERE / 'omega5.json').read_text())
    S6 = {}
    for ij in [(3, 2), (4, 1), (5, 0)]:
        S6 = add(S6, symbol_of_kernel_terms(om[f"{ij[0]},{ij[1]}"], *ij))
    pv = pairing_vector6(S6)
    print("W6(k=5) pairings (33, 42a, 42b, 42c):", [str(v) for v in pv])
    print("expect (449/4 sqrt2, -30, 15 sqrt2, -9)")
'''
open(os.path.join(WD, 'wt6s_core.py'), 'w').write(core)
print("written wt6s_core.py")
