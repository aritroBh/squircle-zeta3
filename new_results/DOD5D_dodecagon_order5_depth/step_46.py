import pickle, sys, time
sys.path.insert(0, ".")
import mpmath as mp
import dod5d_exact as DX
mp.mp.dps = 130
E5X = pickle.load(open("E5_dod_expanded.pkl","rb"))
r3 = mp.sqrt(3)
KAPm = [mp.mpf(str(k[0]))+mp.mpf(str(k[1]))*r3 for k in DX.KAP]
mons = sorted({(n, le[0], le[1], le[2], le[3], le[4], le[5]) for (n, le) in E5X})
NM = len(mons)
from collections import defaultdict
bypat = defaultdict(list)
for i, m in enumerate(mons):
    bypat[m[1:]].append((i, m[0]))
pats = sorted(bypat)
patdata = [(p, bypat[p]) for p in pats]
print("NM", NM, "patterns", len(pats))
newacc = [mp.mpf(0)]*NM
t0 = time.time()
NN = 40
for k in range(NN):
    u_ = k*mp.mpf('0.05')
    arg = (mp.pi/2)*mp.sinh(u_)
    th = mp.tanh(arg)
    w = (mp.pi/2)*mp.cosh(u_)/(mp.cosh(arg)**2)*mp.mpf('0.05')
    for tt in [(1+th)/2, (1-th)/2]:
        LT = mp.log(tt); LU = mp.log(1-tt)
        W0 = mp.log(1+KAPm[0]*tt); W1 = mp.log(1+KAPm[1]*tt)
        W2 = mp.log(1+KAPm[2]*tt); W3 = mp.log(1+KAPm[3]*tt)
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
el = time.time()-t0
print("per node-eval: %.3f ms; est for 30k node-evals: %.0f s" % (el/(2*NN)*1000, el/(2*NN)*30000))
