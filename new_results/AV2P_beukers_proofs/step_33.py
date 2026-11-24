from mpmath import mp, mpf, zeta, log, pi, nprint, polylog, clsin, catalan, pslq
mp.dps = 120
G = catalan; Cl = clsin(2, pi/3); z3 = zeta(3); L2 = log(2)
A = mpf('0.43125077382323234847273083625968838273200772102112')
B = mpf('0.62811756197612089256977406604011504987002593669148')
Li3_half = polylog(3, mpf('0.5'))
basis = [1, pi, pi**2, pi**3, z3, pi*z3, L2, pi*L2, pi**2*L2, L2**2, pi*L2**2, L2**3,
         Cl, pi*Cl, G, pi*G, Li3_half, z3*L2, Cl*L2, log(3), pi*log(3), polylog(3, mpf('0.25'))]
names = ['1','pi','pi2','pi3','z3','pi*z3','log2','pi*log2','pi2*log2','log2^2','pi*log2^2',
         'log2^3','Cl2','pi*Cl2','G','pi*G','Li3(1/2)','z3*log2','Cl2*log2','log3','pi*log3','Li3(1/4)']
def hunt(target, label):
    v = [target] + basis
    rel = pslq(v, tol=mpf('1e-110'), maxcoeff=10**20, maxsteps=10000)
    print(f"--- {label} ---  rel[0] = {None if rel is None else rel[0]}")
    if rel is not None and rel[0] != 0:
        err = abs(sum(mpf(c)*x for c,x in zip(rel, v)))
        print("  residual:", nprint(err, 4))
        out = [f"{rel[0]}*T"]
        for c,n in zip(rel[1:], names):
            if c: out.append(f"{c}*{n}")
        print("  " + " + ".join(out))
hunt(A, "A")
hunt(B, "B")
# Also numerically: A-B?
print("A-B =", nprint(A-B, 30))
print("A =", nprint(A, 30))
print("B =", nprint(B, 30))
