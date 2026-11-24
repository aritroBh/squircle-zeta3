import pickle, sys, time
sys.path.insert(0, ".")
import flint
import dod5d_exact as DX
flint.ctx.prec = 472
E5X = pickle.load(open("E5_dod_expanded.pkl","rb"))
def q3f(c): return flint.arb(str(c[0])) + flint.arb(str(c[1]))*flint.arb(3).sqrt()
KAPa = [q3f(k) for k in DX.KAP]
PI = flint.arb.pi(); ONE = flint.arb(1)
mons = sorted({(n, le[0], le[1], le[2], le[3], le[4], le[5]) for (n, le) in E5X})
def level_sum(m, level):
    n,a,b,c0,c1,c2,c3 = m
    h = ONE/(2**level)
    acc = flint.arb(0); k = 0
    while True:
        u_ = k*h
        arg = (PI/2)*u_.sinh()
        th = arg.tanh()
        w = (PI/2)*u_.cosh()/(arg.cosh()**2)*h
        ts = [(ONE+th)/2] if k == 0 else [(ONE+th)/2, (ONE-th)/2]
        contrib = 0.0
        for tt in ts:
            f = tt**n * tt.log()**a * (ONE-tt).log()**b
            for i, cc in enumerate([c0,c1,c2,c3]):
                if cc: f = f*(ONE+KAPa[i]*tt).log()**cc
            contrib = max(contrib, abs(float((w*f).mid())))
            acc = acc + w*f
        k += 1
        if contrib < 1e-130 and k > 10: break
    return acc
import random
random.seed(3)
sample = random.sample(mons, 12) + [(0,5,0,0,0,0,0), (0,0,5,0,0,0,0), (8,5,5,0,0,0,0), (8,0,0,5,5,5,5), (0,0,0,0,0,0,5), (4,0,0,0,0,5,0)]
worst = []
for m in sample:
    a6 = level_sum(m, 6); a7 = level_sum(m, 7); a8 = level_sum(m, 8)
    d76 = abs(float((a7-a6).mid())); d87 = abs(float((a8-a7).mid()))
    worst.append((d87, d76, m))
    print(m, "d76=%.3g d87=%.3g" % (d76, d87))
