import pickle, sys, time
sys.path.insert(0, ".")
import mpmath as mp
import dod5d_exact as DX
mp.mp.dps = 135
E5X = pickle.load(open("E5_dod_expanded.pkl","rb"))
r3 = mp.sqrt(3)
KAPm = [mp.mpf(str(k[0]))+mp.mpf(str(k[1]))*r3 for k in DX.KAP]
mons = sorted({(n, le[0], le[1], le[2], le[3], le[4], le[5]) for (n, le) in E5X})
midx = {m: i for i, m in enumerate(mons)}
NM = len(mons)
print("monomials:", NM)
acc = [mp.mpf(0)]*NM
h = mp.mpf(1)
t0 = time.time()
prev = None
for level in range(12):
    h = h/2
    k = 0
    newacc = [mp.mpf(0)]*NM
    while True:
        u_ = k*h
        arg = (mp.pi/2)*mp.sinh(u_)
        th = mp.tanh(arg)
        t = (1+th)/2
        w = (mp.pi/2)*mp.cosh(u_)/(mp.cosh(arg)**2)*h
        if k == 0:
            nodes = [(t, w)]
        else:
            nodes = [(t, w), (1-t, w)]
        contrib = mp.mpf(0)
        for (tt, ww) in nodes:
            LT = mp.log(tt); LU = mp.log(1-tt)
            W = [mp.log(1+KAPm[i]*tt) for i in range(4)]
            contrib = max(contrib, ww*(abs(LT)+abs(LU)))
            for (n, a, b, c0, c1, c2, c3) in mons:
                v = ww*tt**n*LT**a*LU**b
                if c0: v *= W[0]**c0
                if c1: v *= W[1]**c1
                if c2: v *= W[2]**c2
                if c3: v *= W[3]**c3
                newacc[midx[(n,a,b,c0,c1,c2,c3)]] += v
        k += 1
        if contrib < mp.mpf(10)**(-130) and k > 10: break
        if k > 100000: raise RuntimeError("nc")
    acc = [a+b for a, b in zip(acc, newacc)]
    if prev is not None:
        diff = max(abs(acc[i]-prev[i]) for i in range(NM))
    else:
        diff = None
    print("level %d: nodes<=%d diff=%s (%.1fs)" % (level, 2*k+1, mp.nstr(diff,3) if diff else '-', time.time()-t0), flush=True)
    if diff is not None and diff < mp.mpf(10)**(-125): break
    prev = list(acc)
pickle.dump({m: mp.nstr(acc[midx[m]], 135) for m in mons}, open("itab130.pkl","wb"))
print("saved itab130.pkl")
