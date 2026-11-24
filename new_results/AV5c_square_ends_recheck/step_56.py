import flint, mpmath, json, time
from mpmath import mp, mpf, pi, nstr
C = json.load(open('ckpt_constants.json'))
def gv(k): return mpf(C[k])
mp.dps = 240
lam, L, piv, sq2, b = gv('lam'), gv('L'), pi, gv('sq2'), gv('b')
c = gv('c8'); s = gv('s8')
l1mc, l1pc = mp.log(1-c), mp.log(1+c)
G = gv('catalan')
consts = {'lam':lam, 'L':L, 'ln(1-c)':l1mc, 'ln(1+c)':l1pc,
 'pi^2':piv**2, 'lam^2':lam**2, 'L^2':L**2, 'lam*L':lam*L,
 'ln(1-c)^2':l1mc**2, 'ln(1+c)^2':l1pc**2, 'ln(1-c)*ln(1+c)':l1mc*l1pc,
 'lam*ln(1-c)':lam*l1mc, 'lam*ln(1+c)':lam*l1pc, 'L*ln(1-c)':L*l1mc, 'L*ln(1+c)':L*l1pc,
 'Li2(c^2)':gv('Li2_8_c8sq'), 'Li2(s^2)':gv('Li2_8_om8'), 'Li2(r)':gv('Li2_8_r8'),
 'Li2(b^2)':gv('Li2_b2'), 'Li2(1/2)':gv('Li2_half'), 'Li2(1-1/sq2)':gv('Li2_c1misq2'),
 'Li2(2-sq2)':gv('Li2_tomsq2'), 'G':G, '1':mpf(1)}
names, vals = list(consts.keys()), [+v for v in consts.values()]
basis, bnames = [], []
for nm, v in zip(names, vals):
    basis.append(v); bnames.append(nm)
    basis.append(c*v); bnames.append('c*'+nm)
res = json.load(open('ckpt_Mk_cornerintegral.json'))
M1 = mpf(res['3pi/4'][1])
x = basis + [M1]; bnames = bnames + ['M1(3pi/4)']
print("lattice dim:", len(x))
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
t0=time.time()
out = int_rel(x, 215, 235)
print("LLL time", time.time()-t0)
for norm2, vec, last in out[:2]:
    print("norm2:", norm2, "last:", last)
    if abs(last) < 10**8 and norm2 < 10**10:
        sgn = 1 if vec[-1]>0 else -1
        vec = [sgn*a for a in vec]
        print([(n_, r_) for n_, r_ in zip(bnames, vec) if r_])
with open('progress.log','a') as f:
    f.write("step_55: LLL identification attempt for M1(3pi/4), 48-dim lattice, C=1e215\n")
