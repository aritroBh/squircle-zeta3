import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, log, nstr
C = json.load(open('ckpt_constants.json'))
paper = json.load(open('ckpt_paper_Mk.json'))
def gv(k): return mpf(C[k])
mp.dps = 240
sq2v, lamv, Lv, Li2bv, piv = gv('sq2'), gv('lam'), gv('L'), gv('Li2_b'), pi
z3v, Lchiv, Li3b2v, Li3isq2v = gv('zeta3'), gv('Lchi8_3'), gv('Li3_b2'), gv('Li3_isq2')
grp_w3 = {'zeta3':z3v, 'Lchi8_3':Lchiv, 'Li3(b^2)':Li3b2v, 'Li3(1/sq2)':Li3isq2v}
grp_p3 = {'pi^2*lam':piv**2*lamv, 'pi^2*L':piv**2*Lv, 'lam^3':lamv**3, 'lam^2*L':lamv**2*Lv,
          'lam*L^2':lamv*Lv**2, 'L^3':Lv**3, 'lam*Li2(b)':lamv*Li2bv, 'L*Li2(b)':Lv*Li2bv}
grp_w2 = {'pi^2':piv**2, 'lam^2':lamv**2, 'L^2':Lv**2, 'lam*L':lamv*Lv, 'Li2(b)':Li2bv}
grp_w1 = {'lam':lamv, 'L':Lv, '1':mpf(1)}
basis, names = [], []
for grp in [grp_w3, grp_p3, grp_w2, grp_w1]:
    for nm, v in grp.items():
        basis.append(+v); names.append(nm)
    for nm, v in grp.items():
        basis.append(+(sq2v*v)); names.append('sq2*'+nm)
d3 = +mpf(paper['d3'])
x = basis + [d3]
print("basis size:", len(x))
t0 = time.time()
mp.dps = 220
rel = mpmath.pslq([+v for v in x], tol=mpf(10)**(-180), maxcoeff=10**6, maxsteps=50000)
print("time", time.time()-t0)
if rel is None:
    print("NONE at 220")
else:
    rel = [int(r_) for r_ in rel]
    sgn = 1 if rel[-1]>0 else -1
    rel = [sgn*r_ for r_ in rel]
    nz = [(n_, r_) for n_, r_ in zip(names+['d3'], rel) if r_]
    mp.dps = 235
    res = sum(mpf(r_)*v for r_, v in zip(rel, x))
    print(f"maxc={max(abs(r_) for r_ in rel)} residual@235={nstr(res,4)}")
    print(nz)
    json.dump({'d3_rel': rel, 'names': names+['d3']}, open('ckpt_pslq_d3.json','w'))
with open('progress.log','a') as f:
    f.write("step_43: PSLQ d3 with 40-element weight-3 basis (z3,Lchi,Li3(b^2),Li3(1/sq2)+products+w2+w1, x{1,sq2}) at dps 220\n")
