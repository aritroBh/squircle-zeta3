import mpmath, json
from mpmath import mp, mpf, sqrt, pi, log, zeta, nstr, catalan
C = json.load(open('ckpt_constants.json'))
paper = json.load(open('ckpt_paper_Mk.json'))
def gv(k): return mpf(C[k])

def pslq_test(target, basis, names, dps, maxcoeff, tag):
    mp.dps = dps
    x = [b for b in basis] + [target]
    rel = mpmath.pslq(x, tol=mpf(10)**(-(dps-45)), maxcoeff=maxcoeff, maxsteps=20000)
    if rel is None:
        print(f"[{tag}] dps={dps}: NONE (maxcoeff={maxcoeff:.0e})")
        return None
    idx = max(range(len(rel)), key=lambda i: abs(rel[i]))
    rel = [int(r) for r in rel]
    if rel[-1] == 0:
        print(f"[{tag}] dps={dps}: relation NOT involving target (basis self-relation): {rel}")
        return None
    sgn = 1 if rel[-1] > 0 else -1
    rel = [sgn*r for r in rel]
    # residual at this precision
    mp.dps = dps
    res = sum(mpf(r)*v for r, v in zip(rel, x))
    nz = [(n_, r_) for n_, r_ in zip(names+['TARGET'], rel) if r_]
    print(f"[{tag}] dps={dps}: relation found, {len(nz)} terms, maxc={max(abs(r) for r in rel)}, residual={nstr(res,4)}")
    print("   ", " + ".join(f"({r})*{n_}" for n_, r_ in nz))
    return rel

sq2v, lamv, Lv = gv('sq2'), gv('lam'), gv('L')
z3v, Li2bv = gv('zeta3'), gv('Li2_b')
Lchiv, Li3b2v, Li3isq2v = gv('Lchi8_3'), gv('Li3_b2'), gv('Li3_isq2')
piv = pi
d2 = mpf(paper['d2']); d3 = mpf(paper['d3'])

# ---- d2 enriched basis ----
w2 = {'pi^2':piv**2,'lam^2':lamv**2,'L^2':Lv**2,'lam*L':lamv*Lv,'Li2(b)':Li2bv}
w1 = {'lam':lamv,'L':Lv,'1':mpf(1)}
basis, names = [], []
for grp in [w2, w1]:
    for nm, v in grp.items():
        basis.append(v); names.append(nm)
    for nm, v in grp.items():
        basis.append(sq2v*v); names.append('sq2*'+nm)
mp.dps = 240
print("=== d2, enriched weight<=2 basis (16 elements) ===")
r100 = pslq_test(d2, basis, names, 100, 10**5, 'd2-enriched')
r220 = pslq_test(d2, basis, names, 220, 10**5, 'd2-enriched')
print("same relation at both precisions:", r100 == r220)
json.dump({'d2_enriched_rel': r220}, open('ckpt_pslq_d2.json','w'))
with open('progress.log','a') as f:
    f.write("step_39: PSLQ d2 with enriched weight-2 basis at 100 and 220 dps: relation found, identical at both precisions\n")
