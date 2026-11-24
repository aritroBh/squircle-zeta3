import mpmath, json
from mpmath import mp, mpf, sqrt, pi, log, nstr
C = json.load(open('ckpt_constants.json'))
paper = json.load(open('ckpt_paper_Mk.json'))
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
d2 = +mpf(paper['d2'])
x = basis + [d2]
for dps in [160, 210]:
    mp.dps = dps
    xx = [+v for v in x]
    rel = mpmath.pslq(xx, tol=mpf(10)**(-(dps-35)), maxcoeff=10**6, maxsteps=30000)
    if rel is None:
        print(f"dps={dps}: NONE"); continue
    rel = [int(r_) for r_ in rel]
    sgn = 1 if rel[-1]>0 else -1
    rel = [sgn*r_ for r_ in rel]
    nz = [(n_, r_) for n_, r_ in zip(names+['d2'], rel) if r_]
    mp.dps = 235
    res = sum(mpf(r_)*v for r_, v in zip(rel, x))
    mp.dps = dps
    print(f"dps={dps}: maxc={max(abs(r_) for r_ in rel)} residual@235={nstr(res,4)}")
    print("   ", nz)
