import flint, mpmath, json
from mpmath import mp, mpf, pi, nstr
C = json.load(open('ckpt_constants.json'))
def gv(k): return mpf(C[k])
mp.dps = 240
lam, L, piv, sq2 = gv('lam'), gv('L'), pi, gv('sq2')
c = gv('c8')
l1mc = mp.log(1-c)
G = gv('catalan')
consts = {'lam':lam, 'L':L, 'ln(1-c)':l1mc,
 'pi^2':piv**2, 'lam^2':lam**2, 'L^2':L**2, 'lam*L':lam*L,
 'ln(1-c)^2':l1mc**2, 'lam*ln(1-c)':lam*l1mc, 'L*ln(1-c)':L*l1mc,
 'Li2(c^2)':gv('Li2_8_c8sq'), 'Li2(r)':gv('Li2_8_r8'), 'Li2(b^2)':gv('Li2_b2'), 'Li2(1/2)':gv('Li2_half'),
 'G':G, '1':mpf(1)}
names, vals = list(consts.keys()), [+v for v in consts.values()]
basis, bnames = [], []
for nm, v in zip(names, vals):
    basis.append(v); bnames.append(nm)
    basis.append(c*v); bnames.append('c*'+nm)
def int_rel(xvals, Cexp, dps):
    mp.dps = dps
    n = len(xvals)
    Cf = mpf(10)**Cexp
    rows = []
    for i, v in enumerate(xvals):
        cv = int(mpmath.nint(Cf * (+v)))
        row = [0]*(n+1); row[i] = 1; row[n] = cv
        rows.append(row)
    B = flint.fmpz_mat(rows).lll(delta=mpf('0.99'))
    out = []
    for r in range(n):
        vec = [int(B[r,cc]) for cc in range(n)]
        out.append((sum(a*a for a in vec), vec, int(B[r,n])))
    out.sort()
    return out
out = int_rel(basis, 210, 230)
for norm2, vec, last in out[:3]:
    nz = [(n_, r_) for n_, r_ in zip(bnames, vec) if r_]
    mp.dps = 235
    resv = sum(mpf(r_)*v for r_, v in zip(vec, basis))
    print("norm2:", norm2, "last:", last, "residual@235:", nstr(resv,4))
    print("  ", nz)
