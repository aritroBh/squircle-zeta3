import pickle, sys, time
sys.path.insert(0, ".")
import flint
import dod5d_exact as DX
flint.ctx.prec = 472
E5X = pickle.load(open("E5_dod_expanded.pkl","rb"))
def q3f(c): return flint.arb(str(c[0])) + flint.arb(str(c[1]))*flint.arb(3).sqrt()
KAPa = [q3f(k) for k in DX.KAP]
PI = flint.arb.pi(); ONE = flint.arb(1); ZERO = flint.arb(0)
mons = sorted({(n, le[0], le[1], le[2], le[3], le[4], le[5]) for (n, le) in E5X})
NM = len(mons)
from collections import defaultdict
bypat = defaultdict(list)
for i, m in enumerate(mons):
    bypat[m[1:]].append((i, m[0]))
pats = sorted(bypat)
patdata = [(p, bypat[p]) for p in pats]
# replicate vectorized loop for levels 0..7, but track diffs per monomial
h = ONE; prev = None
for level in range(8):
    h = h/2
    acc = [ZERO]*NM
    k = 0
    while True:
        u_ = k*h
        arg = (PI/2)*u_.sinh()
        th = arg.tanh()
        w = (PI/2)*u_.cosh()/(arg.cosh()**2)*h
        ts = [(ONE+th)/2] if k == 0 else [(ONE+th)/2, (ONE-th)/2]
        contrib = 0.0
        for tt in ts:
            LT = tt.log(); LU = (ONE-tt).log()
            W0 = (ONE+KAPa[0]*tt).log(); W1 = (ONE+KAPa[1]*tt).log()
            W2 = (ONE+KAPa[2]*tt).log(); W3 = (ONE+KAPa[3]*tt).log()
            contrib = max(contrib, abs(float(w.mid()))*(abs(float(LT.mid()))+abs(float(LU.mid()))))
            pw = [[ONE],[ONE],[ONE],[ONE],[ONE],[ONE]]
            base = [LT, LU, W0, W1, W2, W3]
            for j in range(6):
                bj = base[j]; p = ONE
                for e in range(5):
                    p = p*bj; pw[j].append(p)
            tp = [ONE]; p = ONE
            for n in range(8):
                p = p*tt; tp.append(p)
            for (pat, lst) in patdata:
                v = w*pw[0][pat[0]]*pw[1][pat[1]]
                if pat[2]: v = v*pw[2][pat[2]]
                if pat[3]: v = v*pw[3][pat[3]]
                if pat[4]: v = v*pw[4][pat[4]]
                if pat[5]: v = v*pw[5][pat[5]]
                for (i, n) in lst:
                    acc[i] = acc[i] + v*tp[n]
        k += 1
        if contrib < 1e-125 and k > 10: break
    if prev is not None:
        diffs = sorted(((abs(float((acc[i]-prev[i]).mid())), mons[i]) for i in range(NM)), reverse=True)
        print("level %d: worst diffs:" % level, [(m, "%.2g" % d) for d, m in diffs[:5]])
    prev = list(acc)
