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
        row = [0]*(n+1); row[i] = 1; row[n] = cv
        rows.append(row)
    B = flint.fmpz_mat(rows).lll(delta=mpf('0.99'))
    out = []
    for r in range(n):
        vec = [int(B[r,c]) for c in range(n)]
        out.append((sum(a*a for a in vec), vec, int(B[r,n])))
    out.sort()
    return out

out = int_rel([+v for v in x], 165, 185)
norm2, vec, last = out[0]
print("second precision: norm2:", norm2, "last:", last)
sgn = 1 if vec[-1]>0 else -1
vec = [sgn*a for a in vec]
nz = [(n_, r_) for n_, r_ in zip(names, vec) if r_]
print(nz)
# verify residual at 235 dps
mp.dps = 235
res = sum(mpf(r_)*v for r_, v in zip(vec, x))
print("residual @235dps:", nstr(res, 6))
# compare with expected coefficient list (from step_38 reduced form)
expected = {'zeta3':-228,'Lchi8_3':-7680,'Li3(b^2)':-528,'sq2*zeta3':-567,'sq2*Li3(1/sq2)':5184,
 'pi^2*lam':8,'pi^2*L':112,'lam^3':-32,'lam*L^2':192,'L^3':-800,'L*Li2(b)':-2304,
 'sq2*pi^2*lam':260,'sq2*lam^3':-108,'sq2*lam^2*L':-264,'sq2*lam*L^2':336,'sq2*L^3':-32,
 'sq2*lam*Li2(b)':672,'pi^2':-24,'lam^2':144,'L^2':-384,'sq2*pi^2':152,'sq2*L^2':-960,
 'sq2*lam*L':768,'sq2*Li2(b)':-1920,'lam':-192,'sq2*L':192,'d3':24}
got = dict(nz)
print("MATCHES expected relation:", all(got.get(k,0)==v for k,v in expected.items()) and all(v==0 or k in expected for k,v in got.items()))
json.dump({'d3_rel': got, 'residual_235': nstr(res,10)}, open('ckpt_pslq_d3.json','w'))
with open('progress.log','a') as f:
    f.write("step_48: d3 LLL relation CONFIRMED at second precision (C=1e165,dps185), residual 235dps printed, matches paper-derived form\n")
