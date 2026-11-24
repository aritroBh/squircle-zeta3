import sympy as sp
import mpmath as mp
import pickle
mp.mp.dps = 160
pi = mp.pi; z = mp.zeta
e = {}
for key, val in pickle.load(open("e_rm_ext.pkl","rb")).items():
    r, m = map(int, key.split(",")); e[(r,m)] = sp.Rational(val)
Icache = {}
def I(m, j):
    if (m,j) not in Icache:
        Icache[(m,j)] = mp.quad(lambda s: s**m * mp.log1p(mp.e**(-s))**j, [0, 1, 4, 16, 64, mp.inf])
    return Icache[(m,j)]
def Dtilde(n, r):
    return 2*mp.fsum(mp.mpf(str(e[(r,n-1-j)])) * I(n-1-j,j) / mp.factorial(j) for j in range(1, n-2*r))
def odd_partitions(n, minpart=3):
    if n == 0: yield (); return
    if n < minpart: return
    for first in range(minpart, n+1, 2):
        for rest in odd_partitions(n-first, first):
            yield (first,)+rest
def basis(n):
    return [(a, O) for a in range(n//2+1) for O in odd_partitions(n-2*a)]
def bval(b):
    a, O = b
    return z(2)**a * (mp.fprod([z(i) for i in O]) if O else mp.mpf(1))
def bname(b):
    a, O = b
    nm = ""
    if a: nm += f"z2^{a}" if a>1 else "z2"
    for i in O: nm += ("*" if nm else "") + f"z{i}"
    return nm or "1"

n = 12
B = basis(n); Bv = [bval(b) for b in B]
print("basis:", [bname(b) for b in B])
for r in range(4):
    val = Dtilde(n, r)
    rel = mp.pslq([val]+Bv, maxcoeff=10**16, maxsteps=10**7)
    if rel is None:
        print(f"r={r}: FAILED"); continue
    resid = abs(rel[0]*val + mp.fsum(mp.mpf(c)*b for c,b in zip(rel[1:],Bv)))
    q = sp.Rational(-rel[0]); coeffs = [sp.Rational(c)/q for c in rel[1:]]
    print(f"r={r}: resid={mp.nstr(resid,3)}")
    print("   ", " + ".join(f"({c})*{bname(b)}" for c,b in zip(coeffs,B) if c))
