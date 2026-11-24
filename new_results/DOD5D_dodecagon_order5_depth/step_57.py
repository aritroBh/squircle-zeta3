import subprocess, pickle, os
pid = pickle.load(open("itab_pid.pkl","rb"))
try: os.kill(pid, 9)
except Exception as e: print("kill:", e)
import flint, time
def one_monomial(m, BITS, levels=11):
    flint.ctx.prec = BITS
    import dod5d_exact as DX
    def q3f(c): return flint.arb(str(c[0])) + flint.arb(str(c[1]))*flint.arb(3).sqrt()
    KAPa = [q3f(k) for k in DX.KAP]
    PI = flint.arb.pi(); ONE = flint.arb(1); TWO = flint.arb(2)
    n,a,b,c0,c1,c2,c3 = m
    h = ONE; prev = None
    for level in range(levels):
        h = h/2
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
            if contrib < 1e-300 and k > 10: break
        if prev is not None:
            d = abs(float((acc-prev).mid()))
        else: d = None
        prev = acc
    return acc, d
m = (8, 5, 0, 1, 1, 1, 1)
for BITS in [472, 1000, 2000]:
    acc, d = one_monomial(m, BITS)
    print(BITS, "val:", acc.mid().str(dps=40, condense=10), "rad:", acc.rad().str(dps=3), "last diff:", d)
