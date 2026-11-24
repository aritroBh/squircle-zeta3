"""prove_V56.py -- PROOF of the Newton identities at depths 5,6, NON-ADJACENT keys, all (d,m).
Each identity  e(delta,s) + e(delta,s+1) == Lambda^{(delta)}(s)  is an exact rational identity
in (d,m) over Q(i). Method: clear denominators (product of all Pi-factors), expand, check 0.
Usage: python prove_V56.py DELTA S
"""
import sys, time
import sympy as sp
import tow3lib as T3

d, m = T3.d, T3.m

def identity(delta, s):
    lhs = sp.expand(T3.e(delta, s) + T3.e(delta, s+1))
    rhs = T3.Lam(delta, s)
    t0 = time.time()
    rem = sp.together(lhs - rhs)
    num = sp.fraction(rem)[0]
    num = sp.expand(num)
    is0 = num == 0
    print(f'non-adjacent ({delta},{s}): numerator == 0 : {is0}  ({round(time.time()-t0,1)}s)', flush=True)
    return is0

if __name__ == '__main__':
    delta, s = int(sys.argv[1]), int(sys.argv[2])
    ok = identity(delta, s)
    assert ok
