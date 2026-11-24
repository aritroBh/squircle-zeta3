import sys
sys.dont_write_bytecode = True
work = '<project>/research/08_theory/M5a_next_coefficient'
sys.path.insert(0, work)
sys.path.insert(0, f'{work}/../A3_W4_proof')
from itertools import product
from fractions import Fraction as F

# dimension of free Lie algebra on 2 letters (Witt): n=1..6
def mobius(n):
    if n == 1: return 1
    for p in range(2, n+1):
        if n % p == 0:
            if (n//p) % p == 0: return 0
            return -mobius(n//p)
def witt(n):
    return sum(mobius(d)*2**(n//d) for d in range(1, n+1) if n % d == 0)//n
print("Witt dimensions Lie^1..Lie^6 on 2 letters:", [witt(n) for n in range(1, 7)])

# linear independence of B3, B4 Lie words used in the pairings
import w4_symbol as WS
from w4_symbol import Q, ONE
def bracket(x, y):
    return WS.add(WS.concat(x, y), WS.scale(WS.concat(y, x), -1))
X, Y = {(0,): ONE}, {(1,): ONE}
P = bracket(X, Y)
B3 = [bracket(X, P), bracket(Y, P)]
B4 = [bracket(X, B3[0]), bracket(Y, B3[0]), bracket(Y, B3[1])]
import sympy as sp
def rank(lies, n):
    words = list(product(range(2), repeat=n))
    M = sp.Matrix([[F(str(lie.get(w, Q(0)).p)) for w in words] for lie in lies])
    return M.rank()
print("rank B3 =", rank(B3, 3), " (= dim Lie^3 = 2)")
print("rank B4 =", rank(B4, 4), " (= dim Lie^4 = 3)")
print("P =", {w: str(c) for w, c in P.items()})
print("B3[0] = [X,P] =", {w: str(c) for w, c in B3[0].items()})
print("B3[1] = [Y,P] =", {w: str(c) for w, c in B3[1].items()})
for i, A4 in enumerate(B4):
    print(f"B4[{i}] =", {w: str(c) for w, c in A4.items()})