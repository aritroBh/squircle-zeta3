import sympy as sp
import mpmath as mp
import pickle

mp.mp.dps = 100
pi = mp.pi; z = mp.zeta
e = {}
with open("e_rm.pkl","rb") as f:
    for key, val in pickle.load(f).items():
        r, m = map(int, key.split(","))
        e[(r,m)] = sp.Rational(val)

Icache = {}
def I(m, j):
    if (m,j) not in Icache:
        Icache[(m,j)] = mp.quad(lambda s: s**m * mp.log1p(mp.e**(-s))**j, [0, 1, 4, 16, 64, mp.inf])
    return Icache[(m,j)]

def Dtilde(n, r):
    return 2*mp.fsum(mp.mpf(str(e[(r,n-1-j)])) * I(n-1-j,j) / mp.factorial(j) for j in range(1, n-2*r))

bases = {
    2: [z(2)],
    3: [z(3)],
    4: [z(4)],
    5: [z(5), z(2)*z(3)],
    6: [z(6), z(3)**2],
    7: [z(7), z(2)*z(5), z(2)**2*z(3)],
    8: [z(8), z(3)*z(5), z(2)*z(3)**2],
    9: [z(9), z(2)*z(7), z(2)**2*z(5), z(2)**3*z(3), z(3)**3],
}
names = {
    2: ["z2"], 3: ["z3"], 4: ["z4"], 5: ["z5","z2z3"], 6: ["z6","z3^2"],
    7: ["z7","z2z5","z2^2z3"], 8: ["z8","z3z5","z2z3^2"], 9: ["z9","z2z7","z2^2z5","z2^3z3","z3^3"],
}

results = {}
print("PSLQ of D~_{n,r} = (pi/sqrt2) D_{n,r} against zeta-product basis")
for n in range(2, 10):
    for r in range((n-2)//2 + 1):
        val = Dtilde(n, r)
        B = bases[n]
        rel = mp.pslq([val] + B, maxcoeff=10**8, maxsteps=10**5)
        if rel is not None:
            resid = abs(rel[0]*val + mp.fsum(mp.mpf(c)*b for c, b in zip(rel[1:], B)))
            # normalize to val = sum a_i B_i with rational a_i
            q = sp.Rational(-rel[0])
            coeffs = [sp.Rational(c)/q for c in rel[1:]]
            results[(n,r)] = (coeffs, resid)
            expr = " + ".join(f"({c})*{nm}" for c, nm in zip(coeffs, names[n]))
            print(f"n={n} r={r}: resid={mp.nstr(resid,3)}  D~ = {expr}")
        else:
            results[(n,r)] = (None, None)
            print(f"n={n} r={r}: PSLQ FAILED")
with open("Dtilde_relations.pkl","wb") as f:
    pickle.dump({f"{n},{r}": [str(c) for c in v[0]] if v[0] else None for (n,r),v in results.items()}, f)
