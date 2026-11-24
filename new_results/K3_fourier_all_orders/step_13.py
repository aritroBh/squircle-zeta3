import sympy as sp
import mpmath as mp
import pickle

mp.mp.dps = 100
pi = mp.pi
e = {}
with open("e_rm.pkl","rb") as f:
    for key, val in pickle.load(f).items():
        r, m = map(int, key.split(","))
        e[(r,m)] = sp.Rational(val)

# ---- I_{m,j} = int_0^inf s^m log^j(1+e^{-s}) ds ----
Icache = {}
def I(m, j):
    if (m,j) not in Icache:
        Icache[(m,j)] = mp.quad(lambda s: s**m * mp.log1p(mp.e**(-s))**j, [0, 1, 4, 16, 64, mp.inf])
    return Icache[(m,j)]

# D~_{n,r} = (pi/sqrt2) D_{n,r} = 2 * sum_j e_{r,n-1-j} I_{n-1-j,j} / j!
def Dtilde(n, r):
    tot = mp.mpf(0)
    for j in range(1, n-2*r):
        m = n-1-j
        tot += mp.mpf(str(e[(r,m)])) * I(m,j) / mp.factorial(j)
    return 2*tot

# ---- D11 exact e_n for validation of r=0 ----
Z = {i: sp.symbols(f'Z{i}') for i in range(2,11)}
fn = {}
with open('research/08_theory/D11_rhombus_resurgence/step1_exact.txt') as fh:
    for line in fh:
        if line.startswith('f_'):
            lhs, rhs = line.split('=', 1)
            fn[int(lhs.strip()[2:])] = sp.sympify(rhs.strip(), locals={f'Z{i}': Z[i] for i in range(2,11)} | {'pi': sp.pi})
zval = {Z[i]: mp.zeta(i) for i in range(2,11)}; zval[sp.pi] = pi
print("validation: Dtilde(n,0) vs 4*(-1)^n(f_n+f_{n-1}) from D11")
for n in range(2, 10):
    eD11 = sp.expand((-1)**n * (fn[n] + fn[n-1]))
    val = mp.mpf(sp.N(eD11.subs(zval), 100))
    print(f"n={n}: diff = {mp.nstr(abs(Dtilde(n,0) - 4*val), 3)}")
