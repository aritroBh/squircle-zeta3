import mpmath as mp
mp.mp.dps = 200
x = mp.mpf('0.377')
lam = mp.log(2)
lx, l1m, l1p = mp.log(x), mp.log(1-x), mp.log(1+x)
pi2 = mp.pi**2

# (c) Li2(2x/(1+x)) via Euler at u=(1-x)/(1+x): Li2(u)+Li2(1-u) = pi2/6 - log u log(1-u)
u = (1-x)/(1+x)
pred_c = pi2/6 - (lam+lx-l1p)*(l1m-l1p) - mp.polylog(2, u)
print("check (c) Li2(2x/(1+x)) - Euler(u):", mp.nstr(mp.polylog(2, 2*x/(1+x)) - pred_c, 5))

# (d) Li2(-(1-x)/(1+x)) in terms of Li2((1+x)/2): five-term instance; fit
tv = mp.polylog(2, -(1-x)/(1+x))
basis = {
    'Li2((1+x)/2)': mp.polylog(2, (1+x)/2),
    'pi2': pi2, 'lx': lx, 'l1m': l1m, 'l1p': l1p, 'lam': lam,
    'lx*l1m': lx*l1m, 'lx*l1p': lx*l1p, 'l1m*l1p': l1m*l1p,
    'lx*lam': lx*lam, 'l1m*lam': l1m*lam, 'l1p*lam': l1p*lam,
    'lx^2': lx**2, 'l1m^2': l1m**2, 'l1p^2': l1p**2, 'lam^2': lam**2,
}
bnames = list(basis); bvec = [basis[n] for n in bnames]
rel = mp.pslq([tv] + bvec, maxcoeff=10**6, maxsteps=10**6)
terms = [(bnames[i], -mp.mpf(rel[i+1])/rel[0]) for i in range(len(bnames)) if rel[i+1] != 0]
print("(d) Li2(-(1-x)/(1+x)) =", " + ".join(f"({mp.nstr(c,4)})*{n}" for n, c in terms))
