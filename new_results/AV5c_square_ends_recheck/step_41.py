import mpmath, json
from mpmath import mp, mpf, sqrt, pi, log, nstr
C = json.load(open('ckpt_constants.json'))
def gv(k): return mpf(C[k])
mp.dps = 240
sq2v, lamv, Lv, Li2bv, piv = gv('sq2'), gv('lam'), gv('L'), gv('Li2_b'), pi
w2 = {'pi^2':piv**2,'lam^2':lamv**2,'L^2':Lv**2,'lam*L':lamv*Lv,'Li2(b)':Li2bv}
w1 = {'lam':lamv,'L':Lv,'1':mpf(1)}
basis, names = [], []
for grp in [w2, w1]:
    for nm, v in grp.items():
        basis.append(+v); names.append(nm)
    for nm, v in grp.items():
        basis.append(+(sq2v*v)); names.append('sq2*'+nm)
mp.dps = 220
rel = mpmath.pslq(basis, tol=mpf(10)**(-170), maxcoeff=10**5, maxsteps=20000)
nz = [(n_, int(r_)) for n_, r_ in zip(names, rel) if r_]
print("basis self-relation:", nz)
mp.dps = 240
res = sum(mpf(r_)*v for n_, r_ in nz for v in [dict(zip(names,basis))[n_]])
print("residual @240:", nstr(res, 6))
