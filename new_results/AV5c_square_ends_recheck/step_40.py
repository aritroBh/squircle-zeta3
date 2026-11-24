import mpmath, json
from mpmath import mp, mpf, sqrt, pi, log, zeta, nstr
C = json.load(open('ckpt_constants.json'))
paper = json.load(open('ckpt_paper_Mk.json'))
def gv(k): return mpf(C[k])

def pslq_test(target, basis, names, dps, maxcoeff, tag, verify_dps=235):
    mp.dps = dps
    x = [+b for b in basis] + [+target]
    rel = mpmath.pslq(x, tol=mpf(10)**(-(dps-45)), maxcoeff=maxcoeff, maxsteps=20000)
    if rel is None:
        print(f"[{tag}] dps={dps}: NONE (maxcoeff={maxcoeff:.0e})")
        return None
    rel = [int(r_) for r_ in rel]
    if rel[-1] == 0:
        print(f"[{tag}] dps={dps}: basis self-relation only")
        return None
    sgn = 1 if rel[-1] > 0 else -1
    rel = [sgn*r_ for r_ in rel]
    nz = [(n_, r_) for n_, r_ in zip(names+['TARGET'], rel) if r_]
    mp.dps = verify_dps
    res = sum(mpf(r_)*v for r_, v in zip(rel, x))
    mp.dps = dps
    print(f"[{tag}] dps={dps}: {len(nz)} terms, maxc={max(abs(r_) for r_ in rel)}, residual@{verify_dps}dps={nstr(res,4)}")
    print("   ", " + ".join(f"({r_})*{n_}" for n_, r_ in nz))
    return rel

sq2v, lamv, Lv = gv('sq2'), gv('lam'), gv('L')
z3v, Li2bv = gv('zeta3'), gv('Li2_b')
piv = pi
mp.dps = 240
d2 = +mpf(paper['d2'])
w2 = {'pi^2':piv**2,'lam^2':lamv**2,'L^2':Lv**2,'lam*L':lamv*Lv,'Li2(b)':Li2bv}
w1 = {'lam':lamv,'L':Lv,'1':mpf(1)}
basis, names = [], []
for grp in [w2, w1]:
    for nm, v in grp.items():
        basis.append(+v); names.append(nm)
    for nm, v in grp.items():
        basis.append(+(sq2v*v)); names.append('sq2*'+nm)
print("=== d2 enriched, dps=100 ===")
r100 = pslq_test(d2, basis, names, 100, 10**5, 'd2-enriched')
print("=== d2 enriched, dps=220 ===")
r220 = pslq_test(d2, basis, names, 220, 10**5, 'd2-enriched')
print("identical:", r100==r220)
with open('progress.log','a') as f:
    f.write("step_40: PSLQ d2 enriched basis at dps 100 & 220 with 235-dps residual verification\n")
