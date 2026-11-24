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

B = [z(10), z(3)*z(7), z(5)**2, z(2)*z(3)*z(5), z(2)**2*z(3)**2]
names = ["z10","z3z7","z5^2","z2z3z5","z2^2z3^2"]
res = {}
for r in range(5):
    val = Dtilde(10, r)
    rel = mp.pslq([val]+B, maxcoeff=10**14, maxsteps=10**6)
    if rel is None:
        print(f"r={r}: PSLQ FAILED, Dtilde={mp.nstr(val,25)}")
        res[r] = None; continue
    resid = abs(rel[0]*val + mp.fsum(mp.mpf(c)*b for c,b in zip(rel[1:],B)))
    q = sp.Rational(-rel[0]); coeffs = [sp.Rational(c)/q for c in rel[1:]]
    res[r] = coeffs
    print(f"r={r}: resid={mp.nstr(resid,3)}  D~ =", " + ".join(f"({c})*{nm}" for c,nm in zip(coeffs,names)))
with open("Dt10.pkl","wb") as f:
    pickle.dump({r: ([str(c) for c in v] if v else None) for r,v in res.items()}, f)
print("top check D~_{10,4} vs 2*8!*eta(10):", mp.nstr(abs(Dtilde(10,4) - 2*mp.factorial(8)*(1-mp.power(2,-9))*z(10)),3))
