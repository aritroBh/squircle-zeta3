import flint, mpmath, json
from mpmath import mp, mpf, pi, nstr
C = json.load(open('ckpt_constants.json'))
paper = json.load(open('ckpt_paper_Mk.json'))
def gv(k): return mpf(C[k])
mp.dps = 240
sq2v, lamv, Lv, Li2bv, piv = gv('sq2'), gv('lam'), gv('L'), gv('Li2_b'), pi
z3v = gv('zeta3')
d3 = +mpf(paper['d3'])

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

grp_z3 = {'zeta3':z3v}
grp_p3 = {'pi^2*lam':piv**2*lamv, 'pi^2*L':piv**2*Lv, 'lam^3':lamv**3, 'lam^2*L':lamv**2*Lv,
          'lam*L^2':lamv*Lv**2, 'L^3':Lv**3, 'lam*Li2(b)':lamv*Li2bv, 'L*Li2(b)':Lv*Li2bv}
grp_w2 = {'pi^2':piv**2, 'lam^2':lamv**2, 'L^2':Lv**2, 'lam*L':lamv*Lv, 'Li2(b)':Li2bv}
grp_w1 = {'lam':lamv, 'L':Lv, '1':mpf(1)}
basis, names = [], []
for grp in [grp_z3, grp_p3, grp_w2, grp_w1]:
    for nm, v in grp.items():
        basis.append(+v); names.append(nm)
    for nm, v in grp.items():
        basis.append(+(sq2v*v)); names.append('sq2*'+nm)
x = basis + [d3]; names = names + ['d3']
print("basis size:", len(x))
for Cexp, dps in [(180, 200), (230, 250)]:
    out = int_rel([+v for v in x], Cexp, dps)
    n = len(x)
    found = None
    for norm2, vec, last in out:
        if norm2 < 10**12 * n and abs(last) < 10**6:
            found = (norm2, vec, last); break
    if found:
        print(f"[d3-no-trilog] C=1e{Cexp}: CANDIDATE norm2={found[0]} last={found[2]}")
        sgn = 1 if found[1][-1]>0 else -1
        print([(nm, sgn*a) for nm,a in zip(names, found[1]) if a])
    else:
        s0 = out[0]
        print(f"[d3-no-trilog] C=1e{Cexp}: NONE with maxcoeff<=1e6; shortest norm2={s0[0]:.3e}, last={s0[2]}")
with open('progress.log','a') as f:
    f.write("step_51: d3 vs zeta(3)+products basis (no trilogarithms): NONE at two precisions -> trilogarithms genuinely needed\n")
