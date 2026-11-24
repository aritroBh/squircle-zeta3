import flint, mpmath, json, time
from mpmath import mp, mpf, pi, nstr
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
names = names + ['d3']

def int_rel(xvals, Cexp, dps):
    mp.dps = dps
    n = len(xvals)
    Cf = mpf(10)**Cexp
    rows = []
    for i, v in enumerate(xvals):
        cv = int(mpmath.nint(Cf * (+v)))
        row = [0]*(n+1)
        row[i] = 1
        row[n] = cv
        rows.append(row)
    M = flint.fmpz_mat(rows)
    t0=time.time()
    B = M.lll(delta=mpf('0.99'))
    t1=time.time()
    out = []
    for r in range(n):
        vec = [int(B[r,c]) for c in range(n)]
        last = int(B[r,n])
        norm2 = sum(a*a for a in vec)
        out.append((norm2, vec, last))
    out.sort()
    return out, t1-t0

out, tlll = int_rel([+v for v in x], 195, 215)
print(f"LLL time {tlll:.1f}s")
for norm2, vec, last in out[:3]:
    print("norm2:", norm2, " last:", last)
    if norm2 < 10**12:
        sgn = 1 if vec[-1]>0 else -1
        vec = [sgn*a for a in vec]
        nz = [(n_, r_) for n_, r_ in zip(names, vec) if r_]
        print(nz)
with open('progress.log','a') as f:
    f.write("step_47: flint LLL integer-relation attempt on 41-element d3 lattice, C=1e195\n")
