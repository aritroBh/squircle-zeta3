import sympy as sp
import mpmath as mp
import pickle
mp.mp.dps = 100
pi = mp.pi; z = mp.zeta
e = {}
with open("e_rm.pkl","rb") as f:
    for key, val in pickle.load(f).items():
        r, m = map(int, key.split(",")); e[(r,m)] = sp.Rational(val)
Icache = {}
def I(m, j):
    if (m,j) not in Icache:
        Icache[(m,j)] = mp.quad(lambda s: s**m * mp.log1p(mp.e**(-s))**j, [0, 1, 4, 16, 64, mp.inf])
    return Icache[(m,j)]
def Dtilde(n, r):
    return 2*mp.fsum(mp.mpf(str(e[(r,n-1-j)])) * I(n-1-j,j) / mp.factorial(j) for j in range(1, n-2*r))

Z = {i: sp.symbols(f'Z{i}') for i in range(2,11)}
fn = {}
with open('research/08_theory/D11_rhombus_resurgence/step1_exact.txt') as fh:
    for line in fh:
        if line.startswith('f_'):
            lhs_, rhs = line.split('=', 1)
            fn[int(lhs_.strip()[2:])] = sp.sympify(rhs.strip(), locals={f'Z{i}': Z[i] for i in range(2,11)} | {'pi': sp.pi})
e10 = sp.expand(fn[10]+fn[9])
zv = {Z[i]: mp.zeta(i) for i in range(2,11)}
print("Dtilde(10,0) - 4 e_10 =", mp.nstr(abs(Dtilde(10,0) - mp.mpf(sp.N((4*e10).subs(zv),100))),3))
# express e_10 in my basis to compare coefficients
e10z = sp.expand(e10.subs({})).subs({})
# use sympy to rewrite: Z4 = 2/5 Z2^2, Z6 = 8/35 Z2^3, Z8 = ?, Z10 = ?
Z2v, Z3v, Z5v, Z7v = sp.symbols('Z2 Z3 Z5 Z7')
sub = {sp.symbols('Z4'): sp.Rational(2,5)*Z2v**2, sp.symbols('Z6'): sp.Rational(8,35)*Z2v**3,
       sp.symbols('Z8'): sp.Rational(16,175)*Z2v**4, sp.symbols('Z2'): Z2v, sp.symbols('Z3'): Z3v,
       sp.symbols('Z5'): Z5v, sp.symbols('Z7'): Z7v, sp.symbols('Z9'): sp.symbols('Z9'),
       sp.symbols('Z10'): sp.Rational(128,1925)*Z2v**5}
e10r = sp.expand(e10.subs(sub))
print("4 e_10 in basis {z10, z3z7, z5^2, z2z3z5, z2^2z3^2}:")
for mon, nm in [(Z2v**5,'z10(=128/1925 Z2^5)'), (Z3v*Z7v,'z3z7'), (Z5v**2,'z5^2'), (Z2v*Z3v*Z5v,'z2z3z5'), (Z2v**2*Z3v**2,'z2^2z3^2')]:
    c = sp.expand(e10r).coeff(mon, 1)
    c = sp.simplify(sp.expand(e10r).coeff(Z3v*Z7v if nm=='z3z7' else mon,1))
    print(f"  {nm}: 4*{sp.nsimplify(c)}")
