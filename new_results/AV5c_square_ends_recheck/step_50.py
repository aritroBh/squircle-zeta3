import flint, mpmath, json, math
from mpmath import mp, mpf, pi, nstr, catalan
C = json.load(open('ckpt_constants.json'))
paper = json.load(open('ckpt_paper_Mk.json'))
def gv(k): return mpf(C[k])
mp.dps = 240
sq2v, lamv, Lv, Li2bv, piv = gv('sq2'), gv('lam'), gv('L'), gv('Li2_b'), pi
z3v, Lchiv, Li3b2v, Li3isq2v, Gv = gv('zeta3'), gv('Lchi8_3'), gv('Li3_b2'), gv('Li3_isq2'), gv('catalan')
d2 = +mpf(paper['d2']); d3 = +mpf(paper['d3'])

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
        vec = [int(B[r,c]) for c in range(n)]
        out.append((sum(a*a for a in vec), vec, int(B[r,n])))
    out.sort()
    return out

def report(tag, x, names, Cexp, dps, Abound):
    out = int_rel(x, Cexp, dps)
    n = len(x)
    # look for any vector with |last| tiny relative to C and norm < Abound*sqrt(n)
    found = None
    for norm2, vec, last in out:
        if norm2 < Abound**2 * n and abs(last) < 10**6:
            found = (norm2, vec, last); break
    if found:
        norm2, vec, last = found
        print(f"[{tag}] CANDIDATE rel norm2={norm2} last={last}")
        sgn = 1 if vec[-1]>0 else -1
        print([(nm, sgn*a) for nm,a in zip(names,vec) if a])
    else:
        s0 = out[0]
        print(f"[{tag}] NONE with maxcoeff<={Abound:.0e} (C=1e{Cexp}, dps={dps}). shortest: norm2={s0[0]:.3e} last={s0[2]}")
    return out

# (a) AV5's basis
av5 = {'1':mpf(1),'pi':piv,'pi^2':piv**2,'pi^3':piv**3,'lam':lamv,'lam^2':lamv**2,'lam^3':lamv**3,
       'C':Gv,'z3':z3v,'pi*lam':piv*lamv,'pi^2*lam':piv**2*lamv,'pi*C':piv*Gv,
       'sq2':sq2v,'sq2*pi':sq2v*piv,'sq2*lam':sq2v*lamv}
anames = list(av5.keys()); abasis = [+v for v in av5.values()]
print("== control (a): AV5 basis ==")
report('d2-vs-AV5', abasis+[d2], anames+['d2'], 150, 170, 10**6)
report('d3-vs-AV5', abasis+[d3], anames+['d3'], 150, 170, 10**6)
report('d3-vs-AV5', abasis+[d3], anames+['d3'], 210, 230, 10**6)
with open('progress.log','a') as f:
    f.write("step_50: controls: d2,d3 vs AV5's basis -> no relation (LLL, two precisions)\n")
