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

# weight-10 zeta-product basis: z10, z3z7, z5^2, z2z3z5, z2^2z3^2
B = [z(10), z(3)*z(7), z(5)**2, z(2)*z(3)*z(5), z(2)**2*z(3)**2]
names = ["z10","z3z7","z5^2","z2z3z5","z2^2z3^2"]
print("PSLQ for D~_{10,r}:")
res = {}
for r in range(5):
    val = Dtilde(10, r)
    rel = mp.pslq([val]+B, maxcoeff=10**9, maxsteps=10**6)
    resid = abs(rel[0]*val + mp.fsum(mp.mpf(c)*b for c,b in zip(rel[1:],B)))
    q = sp.Rational(-rel[0]); coeffs = [sp.Rational(c)/q for c in rel[1:]]
    res[r] = coeffs
    print(f"r={r}: resid={mp.nstr(resid,3)}  D~ =", " + ".join(f"({c})*{nm}" for c,nm in zip(coeffs,names)))
# r=0 should match D11: 4 e_10 = 4(f_10+f_9)
Z = {i: sp.symbols(f'Z{i}') for i in range(2,11)}
fn = {}
with open('research/08_theory/D11_rhombus_resurgence/step1_exact.txt') as fh:
    for line in fh:
        if line.startswith('f_'):
            lhs_, rhs = line.split('=', 1)
            fn[int(lhs_.strip()[2:])] = sp.sympify(rhs.strip(), locals={f'Z{i}': Z[i] for i in range(2,11)} | {'pi': sp.pi})
e10 = sp.expand(fn[10]+fn[9])  # (-1)^10
zv = {Z[i]: mp.zeta(i) for i in range(2,11)}
print("\nD11 4e_10 numeric:", mp.nstr(mp.mpf(sp.N((4*e10).subs(zv),100)), 30))
print("Dtilde(10,0)     :", mp.nstr(Dtilde(10,0), 30))
print("diff:", mp.nstr(abs(Dtilde(10,0)-mp.mpf(sp.N((4*e10).subs(zv),100))),3))
print("e_10 in zeta products:", sp.nsimplify(e10))
with open("Dt10.pkl","wb") as f: pickle.dump({r: [str(c) for c in v] for r,v in res.items()}, f)
# top check: D~_{10,4} should be 2*8!*eta(10)
print("\ntop check D~_{10,4} vs 2*8!*eta(10):", mp.nstr(abs(Dtilde(10,4) - 2*mp.factorial(8)*(1-mp.power(2,-9))*z(10)),3))
