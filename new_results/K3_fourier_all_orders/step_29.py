import sympy as sp
import mpmath as mp
import pickle
from fractions import Fraction
mp.mp.dps = 120
pi = mp.pi; z = mp.zeta

# --- extend e-table to r<=5, m<=12 ---
x = sp.symbols('x')
N = 16
def trunc(poly, n):
    poly = sp.expand(poly)
    return sum(poly.coeff(x, i)*x**i for i in range(n+1))
cosh_ser = sum(x**(2*i)/sp.factorial(2*i) for i in range(N))
sech1 = sp.series(1/cosh_ser, x, 0, N).removeO().expand()
u_ = sech1 - 1
sech12 = sp.Rational(0)
for j in range(N+1):
    sech12 += sp.binomial(sp.Rational(1,2), j)*u_**j
    sech12 = trunc(sech12, N)
sech32 = trunc(sech1*sech12, N)
expx2 = sum((x/2)**i/sp.factorial(i) for i in range(N+1))
H32 = trunc(expx2*sech32, N)
gd = sp.series(sp.asin(sp.tanh(x)), x, 0, N).removeO().expand()
e = {}
for r in range(6):
    ser = trunc(gd**(2*r)*H32, 12)
    for m in range(13):
        e[(r,m)] = sp.Rational(sp.expand(ser).coeff(x, m))
pickle.dump({f"{r},{m}": str(e[(r,m)]) for r in range(6) for m in range(13)}, open("e_rm_ext.pkl","wb"))

Icache = {}
def I(m, j):
    if (m,j) not in Icache:
        Icache[(m,j)] = mp.quad(lambda s: s**m * mp.log1p(mp.e**(-s))**j, [0, 1, 4, 16, 64, mp.inf])
    return Icache[(m,j)]
def Dtilde(n, r):
    return 2*mp.fsum(mp.mpf(str(e[(r,n-1-j)])) * I(n-1-j,j) / mp.factorial(j) for j in range(1, n-2*r))

# --- automated weight-n basis: zeta2^a * (product of odd zetas >=3), 2a+sum=n ---
def odd_partitions(n, minpart=3):
    # partitions of n into odd parts >= 3
    if n == 0: yield (); return
    if n < minpart: return
    for first in range(minpart, n+1, 2):
        for rest in odd_partitions(n-first, first):
            yield (first,)+rest
def basis(n):
    out = []
    for a in range(0, n//2+1):
        for O in odd_partitions(n-2*a):
            out.append((a, O))
    return out
def bval(b):
    a, O = b
    return z(2)**a * mp.fprod([z(i) for i in O]) if O else z(2)**a
def bname(b):
    a, O = b
    nm = ""
    if a: nm += f"z2^{a}" if a>1 else "z2"
    for i in O: nm += ("*" if nm else "") + f"z{i}"
    return nm or "1"

for n in [11, 12]:
    B = basis(n)
    Bv = [bval(b) for b in B]
    print(f"=== n={n}: basis {[bname(b) for b in B]} ===")
    for r in range((n-2)//2+1):
        val = Dtilde(n, r)
        rel = mp.pslq([val]+Bv, maxcoeff=10**12, maxsteps=10**6)
        if rel is None:
            print(f"  r={r}: PSLQ FAILED"); continue
        resid = abs(rel[0]*val + mp.fsum(mp.mpf(c)*b for c,b in zip(rel[1:],Bv)))
        q = sp.Rational(-rel[0]); coeffs = [sp.Rational(c)/q for c in rel[1:]]
        print(f"  r={r}: resid={mp.nstr(resid,3)}: " + " + ".join(f"({c})*{bname(b)}" for c,b in zip(coeffs,B) if c))
pickle.dump("done", open("ext_progress.pkl","wb"))
