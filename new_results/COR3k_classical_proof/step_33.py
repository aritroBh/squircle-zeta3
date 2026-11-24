import mpmath as mp
mp.mp.dps = 200

x = mp.mpf('0.377')  # generic point
lam = mp.log(2)
lx, l1m, l1p = mp.log(x), mp.log(1-x), mp.log(1+x)
pi2 = mp.pi**2

# basis: Li2(x), Li2(-x), and elementary monomials (deg<=2 in logs) + pi2
basis = {
    'Li2(x)': mp.polylog(2, x),
    'Li2(-x)': mp.polylog(2, -x),
    'pi2': pi2,
    'lx': lx, 'l1m': l1m, 'l1p': l1p, 'lam': lam,
    'lx*l1m': lx*l1m, 'lx*l1p': lx*l1p, 'l1m*l1p': l1m*l1p,
    'lx*lam': lx*lam, 'l1m*lam': l1m*lam, 'l1p*lam': l1p*lam,
    'lx^2': lx**2, 'l1m^2': l1m**2, 'l1p^2': l1p**2, 'lam^2': lam**2,
}
bnames = list(basis.keys())
bvec = [basis[n] for n in bnames]

targets = {
    'Li2(1-x)': mp.polylog(2, 1-x),
    'Li2(1/(1+x))': mp.polylog(2, 1/(1+x)),
    'Li2((1+x)/2)': mp.polylog(2, (1+x)/2),
    'Li2((1-x)/(1+x))': mp.polylog(2, (1-x)/(1+x)),
    'Li2(-(1-x)/(1+x))': mp.polylog(2, -(1-x)/(1+x)),
    'Li2(2x/(1+x))': mp.polylog(2, 2*x/(1+x)),
}
for tn, tv in targets.items():
    rel = mp.pslq([tv] + bvec, maxcoeff=10**6, maxsteps=10**6)
    if rel is None:
        print(tn, ": PSLQ FAILED")
        continue
    assert rel[0] != 0
    terms = [(bnames[i], -mp.mpf(rel[i+1])/rel[0]) for i in range(len(bnames)) if rel[i+1] != 0]
    print(tn, "=", " + ".join(f"({mp.nstr(c,4)})*{n}" for n, c in terms))
