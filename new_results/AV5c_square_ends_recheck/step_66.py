import flint, mpmath, json, time
from mpmath import mp, mpf, pi, nstr, log, polylog, zeta, catalan, sqrt
mp.dps = 300
sq2 = sqrt(2); b = sq2-1; lam = log(2); L = log(1+sq2)
c = mp.cos(3*pi/8); s = mp.sin(3*pi/8)
l1mc = log(1-c)
r = (1-c)/(1+c)
G = catalan
M1 = mpf(json.load(open('ckpt_M1_3pi4_formula.json'))['M1_formula'])
li2args = {'r':r, 'b^2':b**2, 'c^2':c**2, 's^2':s**2, 'c':c, 's':s, '1/(1+s)':1/(1+s),
           'c/(1+s)':c/(1+s), 's/(1+s)':s/(1+s), '1/2':mpf('0.5'), '1-1/sq2':1-1/sq2, '2-sq2':2-sq2}
consts = {'lam':lam, 'L':L, 'ln(1-c)':l1mc, 'ln(1+s)':log(1+s), 'pi^2':pi**2, 'lam^2':lam**2, 'L^2':L**2,
 'lam*L':lam*L, 'ln(1-c)^2':l1mc**2, 'lam*ln(1-c)':lam*l1mc, 'L*ln(1-c)':L*l1mc, 'G':G, '1':mpf(1)}
for nm, aa in li2args.items():
    consts[f'Li2({nm})'] = polylog(2, aa)
names, vals = list(consts.keys()), [+v for v in consts.values()]
basis, bnames = [], []
for nm, v in zip(names, vals):
    basis.append(v); bnames.append(nm)
    basis.append(c*v); bnames.append('c*'+nm)
x = basis + [M1]; bnames = bnames + ['M1']
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
    for rr in range(n):
        vec = [int(B[rr,cc]) for cc in range(n)]
        out.append((sum(a*a for a in vec), vec, int(B[rr,n])))
    out.sort()
    return out
t0=time.time()
out = int_rel([+v for v in x], 260, 280)
print("LLL time", time.time()-t0)
found = None
for norm2, vec, last in out:
    if vec[-1] != 0 and norm2 < 10**12:
        found = (norm2, vec, last); break
if found:
    norm2, vec, last = found
    sgn = 1 if vec[-1]>0 else -1
    vec = [sgn*a for a in vec]
    nz = [(n_, r_) for n_, r_ in zip(bnames, vec) if r_]
    mp.dps = 290
    resv = sum(mpf(r_)*v for r_, v in zip(vec, x))
    print(f"norm2={norm2} residual@290={nstr(resv,4)}")
    print("  ", nz)
else:
    print(f"NONE; shortest {out[0][0]}, last {out[0][2]}")
with open('progress.log','a') as f:
    f.write("step_64: LLL M1(3pi/4) with widened 52-dim basis (12 dilog points)\n")
