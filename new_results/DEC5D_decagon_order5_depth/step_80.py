
import pickle
from fractions import Fraction as F
import mpmath as mp
mp.mp.dps = 250
z2 = mp.pi**2/6; z3 = mp.zeta(3); z5 = mp.zeta(5)
ZV = [mp.mpf(1), z2, z3, z2**2, z5, z2*z3]
ZN = ['1','z2','z3','z2^2','z5','z2z3']
ZW = [0,2,3,4,5,5]
def recon(val, wmax):
    cand = [k for k in range(6) if ZW[k] <= wmax]
    vec = [val] + [ZV[k] for k in cand]
    rel = mp.pslq(vec, tol=mp.mpf('1e-200'), maxcoeff=10**25, maxsteps=10000)
    if rel is None: return None
    rel = [int(r) for r in rel]
    if rel[0] <= 0 or max(abs(r) for r in rel) > 10**25: return None
    out = {}
    for i, k in enumerate(cand):
        if rel[i+1]:
            out[ZN[k]] = F(-rel[i+1], rel[0])
    s = sum(out[ZN[k]]*ZV[k] for k in range(6) if ZN[k] in out)
    if abs(s - val) > mp.mpf('1e-195'): return None
    return out

parsed = pickle.load(open("parsed_mons_dec.pkl","rb"))
localkeys = set()
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    if c==0 and d==0:
        for n in cu: localkeys.add((n,a,b))
exactI = {}
fails = []
for k in sorted(localkeys):
    n,a,b = k
    f = lambda x: x**n * mp.log(x)**a * mp.log(1-x)**b
    v = mp.quad(f, [0,1], method='tanh-sinh')
    r = recon(v, a+b)
    if r is None: fails.append(k)
    else: exactI[k] = r
print("reconstructed:", len(exactI), "failed:", fails)
import collections
import dec5d_lib as D
QS5 = D.QS5
acc = collections.defaultdict(lambda: QS5(0))
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    if not (c==0 and d==0) or sum(mon)!=5: continue
    for n,(p,q) in cu.items():
        cf = QS5(F(p), F(q))
        for zn, zc in exactI[(n,a,b)].items():
            acc[(zn, e, f)] += cf*zc
acc = {k: 1200*v for k, v in acc.items() if v != QS5(0)}
print("exact local-block decomposition:")
for k in sorted(acc):
    print("   %-6s l2^%d lf^%d : %s" % (k[0], k[1], k[2], acc[k]))
# numeric check
r5 = mp.sqrt(5); phi = (1+r5)/2
l2 = mp.log(2); lf = mp.log(phi)
ZVm = {'1': mp.mpf(1), 'z2': z2, 'z3': z3, 'z2^2': z2**2, 'z5': z5, 'z2z3': z2*z3}
val = sum(float(acc[k])*ZVm[k[0]]*l2**k[1]*lf**k[2] for k in acc)
local250 = mp.mpf(pickle.load(open("local250.pkl","rb")))
print("\nnumeric local block:", mp.nstr(val, 50))
print("quad local block :", mp.nstr(local250, 50))
print("residual:", mp.nstr(abs(val-local250), 3))
pickle.dump(acc, open("local_block_exact.pkl","wb"))
