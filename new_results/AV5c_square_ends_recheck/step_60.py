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
 'Li2(r)':gv('Li2_8_r8'), 'Li2(b^2)':gv('Li2_b2'), 'G':G, '1':mpf(1)}
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
out0 = int_rel(basis, 210, 230)
print("basis self-relations (top 2):")
for norm2, vec, last in out0[:2]:
    print("  norm2:", norm2, "last:", last, [(n_,r_) for n_,r_ in zip(bnames,vec) if r_])
res = json.load(open('ckpt_Mk_cornerintegral.json'))
for kk, tag in [(0,'M0'),(1,'M1')]:
    Mv = mpf(res['3pi/4'][kk])
    x = basis + [Mv]
    out = int_rel(x, 210, 230)
    found = None
    for norm2, vec, last in out:
        if abs(last) < 10**8 and norm2 < 10**10 and vec[-1]!=0:
            found = (norm2, vec, last); break
    if found:
        norm2, vec, last = found
        sgn = 1 if vec[-1]>0 else -1
        vec = [sgn*a for a in vec]
        nz = [(n_, r_) for n_, r_ in zip(bnames+[tag], vec) if r_]
        mp.dps = 235
        resv = sum(mpf(r_)*v for r_, v in zip(vec, x))
        print(f"{tag}: norm2={norm2} residual@235={nstr(resv,4)}")
        print("  ", nz)
    else:
        print(f"{tag}: NONE. shortest: {out[0][0]}, last {out[0][2]}")
with open('progress.log','a') as f:
    f.write("step_59: LLL M0/M1(3pi/4) with fully reduced dilog basis\n")
