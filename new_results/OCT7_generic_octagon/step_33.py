import mpmath as mp
mp.mp.dps = 60

def Li21(a, c, N=500):
    a, c = mp.mpf(a), mp.mpf(c)
    return mp.fsum([a**n * c**m / (mp.mpf(n)**2 * m) for n in range(2, N) for m in range(1, n)])

x1 = mp.mpf(-1)/3; x2 = mp.mpf(-1)/9
b1 = mp.mpf(1)/4; b2 = mp.mpf(1)/10
Lam = [Li21(x1, x2/x1), Li21(x2, x1/x2), Li21(b1, b2/b1), Li21(b2, b1/b2)]

l2, l3, l5 = mp.log(2), mp.log(3), mp.log(5)
z2 = mp.pi**2/6; z3 = mp.zeta(3)
def Li2(v): return mp.polylog(2, mp.mpf(v))
def Li3(v): return mp.polylog(3, mp.mpf(v))
pts = [x1, x2, b1, b2, mp.mpf(-1)/2, mp.mpf(-2)/3, mp.mpf(2)/3, mp.mpf(3)/5, mp.mpf(-2), mp.mpf(-3)/2, mp.mpf(-3), mp.mpf(-9)]
logs = [l2, l3, l5]
# basis: 1, logs, z2, Li2(pts), log^2, loglog, z3, Li3(pts), z2*log, Li2*log, log^3, log^2 log
basis = [mp.mpf(1)] + logs + [z2] + [Li2(p) for p in pts]
for a in logs:
    for b in logs: basis.append(a*b)
basis += [z3] + [Li3(p) for p in pts]
for a in logs:
    basis.append(z2*a)
    for p in pts[:10]: basis.append(Li2(p)*a)
for a in logs:
    for b in logs:
        for c in logs: basis.append(a*b*c)
print("basis size:", len(basis))
names = ["Lam1=Li21(-1/3,1/3)","Lam2=Li21(-1/9,3)","Lam3=Li21(1/4,2/5)","Lam4=Li21(1/10,5/2)"]
for nm, L in zip(names, Lam):
    v = basis + [L]
    rel = mp.pslq(v, tol=mp.mpf('1e-52'), maxcoeff=10**8, maxsteps=1000)
    if rel is None:
        print(f"{nm}: NO relation (tol 1e-52, maxcoeff 1e8, basis {len(basis)+1})")
    else:
        print(f"{nm}: RELATION found:", rel)
