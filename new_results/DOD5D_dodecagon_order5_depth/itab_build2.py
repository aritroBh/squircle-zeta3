"""itab_build2.py -- monomial integrals by vectorized tanh-sinh with flint.arb (fast).
int_0^1 t^n (log t)^a (log(1-t))^b prod_i (log(1+kappa_i t))^c_i dt, all monomials of E5_dod_expanded.pkl."""
import pickle, sys, time
sys.path.insert(0, ".")
import flint
import dod5d_exact as DX

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 130
BITS = int(DPS*3.322) + 40
flint.ctx.prec = BITS
E5X = pickle.load(open("E5_dod_expanded.pkl","rb"))
def q3f(c): return flint.arb(str(c[0])) + flint.arb(str(c[1]))*flint.arb(3).sqrt()
KAPa = [q3f(k) for k in DX.KAP]
mons = sorted({(n, le[0], le[1], le[2], le[3], le[4], le[5]) for (n, le) in E5X})
NM = len(mons)
from collections import defaultdict
bypat = defaultdict(list)
for i, m in enumerate(mons):
    bypat[m[1:]].append((i, m[0]))
pats = sorted(bypat)
patdata = [(p, bypat[p]) for p in pats]
PI = flint.arb.pi()
ZERO = flint.arb(0); ONE = flint.arb(1); TWO = flint.arb(2)
acc = [ZERO]*NM
h = ONE; t0 = time.time(); prev = None
EPS = flint.arb(10)**(-(DPS+15)); STOP = (10.0)**(-(DPS-1))
for level in range(16):
    h = h/2
    k = 0
    acc = [ZERO]*NM
    while True:
        u_ = k*h
        arg = (PI/2)*u_.sinh()
        th = arg.tanh()
        tmr = TWO/((TWO*arg).exp()+ONE)     # 1 - tanh(arg), accurate (no cancellation)
        w = (PI/4)*u_.cosh()/(arg.cosh()**2)*h   # dt/du = sech^2(arg)/2 * (pi/2)cosh(u)
        tpl = (ONE+th)/2; om = tmr/2
        ts = [(tpl, om)] if k == 0 else [(tpl, om), (om, tpl)]
        contrib = 0.0
        for tt, omt in ts:
            LT = tt.log(); LU = omt.log()
            W0 = (ONE+KAPa[0]*tt).log(); W1 = (ONE+KAPa[1]*tt).log()
            W2 = (ONE+KAPa[2]*tt).log(); W3 = (ONE+KAPa[3]*tt).log()
            pw = [[ONE],[ONE],[ONE],[ONE],[ONE],[ONE]]
            base = [LT, LU, W0, W1, W2, W3]
            BF = 1.0
            for j in range(6):
                bj = base[j]; p = ONE
                for e in range(5):
                    p = p*bj; pw[j].append(p)
                BF *= (1.0 + abs(float(pw[j][5].mid())))
            contrib = max(contrib, abs(float(w.mid()))*BF)
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
        if contrib < float(EPS.mid()) and k > 10: break
        if k > 200000: raise RuntimeError("no convergence")
    diff = max(abs(float((acc[i]-prev[i]).mid())) for i in range(NM)) if prev is not None else None
    print("level %2d: k=%4d  diff=%.3g  (%.1fs)" % (level, k, diff if diff else -1, time.time()-t0), flush=True)
    if diff is not None and diff < STOP: break
    prev = list(acc)
out = {}
for i, m in enumerate(mons):
    a = acc[i]
    out[m] = a.mid().str(DPS+5, radius=False)
pickle.dump(out, open(f"itab{DPS}.pkl","wb"))
print("saved itab%d.pkl, final diff %.3g, balls rad max %.3g" % (DPS, diff, max(float(a.rad()) for a in acc)))
