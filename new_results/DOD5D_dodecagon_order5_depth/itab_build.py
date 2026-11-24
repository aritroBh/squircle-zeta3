"""itab_build.py -- monomial integrals int_0^1 t^n (log t)^a (log(1-t))^b prod_i (log(1+kappa_i t))^c_i dt
for all (n, a,b,c0..c3) in E5_dod_expanded.pkl, by vectorized tanh-sinh on (0,1)."""
import pickle, sys, time
sys.path.insert(0, ".")
import mpmath as mp
import dod5d_exact as DX

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 130
mp.mp.dps = DPS
E5X = pickle.load(open("E5_dod_expanded.pkl","rb"))
r3 = mp.sqrt(3)
KAPm = [mp.mpf(str(k[0]))+mp.mpf(str(k[1]))*r3 for k in DX.KAP]
mons = sorted({(n, le[0], le[1], le[2], le[3], le[4], le[5]) for (n, le) in E5X})
NM = len(mons)
# group: log-pattern -> list of (monomial index, n)
from collections import defaultdict
bypat = defaultdict(list)
for i, m in enumerate(mons):
    bypat[m[1:]].append((i, m[0]))
pats = sorted(bypat)
# flatten: for each pat, base monomial indices and ns
patdata = [(p, bypat[p]) for p in pats]
acc = [mp.mpf(0)]*NM
h = mp.mpf(1); t0 = time.time(); prev = None
EPS = mp.mpf(10)**(-(DPS-8)); STOP = mp.mpf(10)**(-(DPS-5))
for level in range(14):
    h = h/2
    k = 0
    newacc = [mp.mpf(0)]*NM
    while True:
        u_ = k*h
        arg = (mp.pi/2)*mp.sinh(u_)
        th = mp.tanh(arg)
        w = (mp.pi/2)*mp.cosh(u_)/(mp.cosh(arg)**2)*h
        ts = [(1+th)/2] if k == 0 else [(1+th)/2, (1-th)/2]
        contrib = mp.mpf(0)
        for tt in ts:
            LT = mp.log(tt); LU = mp.log(1-tt)
            W0 = mp.log(1+KAPm[0]*tt); W1 = mp.log(1+KAPm[1]*tt)
            W2 = mp.log(1+KAPm[2]*tt); W3 = mp.log(1+KAPm[3]*tt)
            contrib = max(contrib, w*(abs(LT)+abs(LU)))
            pw = [[mp.mpf(1)],[mp.mpf(1)],[mp.mpf(1)],[mp.mpf(1)],[mp.mpf(1)],[mp.mpf(1)]]
            base = [LT, LU, W0, W1, W2, W3]
            for j in range(6):
                bj = base[j]; p = mp.mpf(1)
                for e in range(5):
                    p *= bj; pw[j].append(p)
            tp = [mp.mpf(1)]; p = mp.mpf(1)
            for n in range(8):
                p *= tt; tp.append(p)
            for (pat, lst) in patdata:
                v = w*pw[0][pat[0]]*pw[1][pat[1]]
                if pat[2]: v *= pw[2][pat[2]]
                if pat[3]: v *= pw[3][pat[3]]
                if pat[4]: v *= pw[4][pat[4]]
                if pat[5]: v *= pw[5][pat[5]]
                for (i, n) in lst:
                    newacc[i] += v*tp[n]
        k += 1
        if contrib < EPS and k > 10: break
        if k > 200000: raise RuntimeError("no convergence")
    acc = [a+b for a, b in zip(acc, newacc)]
    diff = max(abs(acc[i]-prev[i]) for i in range(NM)) if prev is not None else None
    print("level %d: k=%d diff=%s (%.1fs)" % (level, k, mp.nstr(diff,3) if diff else '-', time.time()-t0), flush=True)
    if diff is not None and diff < STOP: break
    prev = list(acc)
pickle.dump({m: mp.nstr(acc[i], DPS+5) for i, m in enumerate(mons)}, open(f"itab{DPS}.pkl","wb"))
print("saved itab%d.pkl, final diff %s" % (DPS, mp.nstr(diff,3)))
