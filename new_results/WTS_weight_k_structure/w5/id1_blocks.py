# Identification step 1: high-precision values of the "block" integrals
#     B[(i,j),pole] = int_b^a  l1^i l2^j /(v - pole) dv,   l1 = log t, l2 = log(1-t),  t = (v-b)(v+a)/(4v),
# for every (key,pole) that occurs in the exact cubic kernel / k=4 defect kernel (exact_cubic_kernel_w5.pkl).
# From them: c_0..c_3 (coefficients of 1, k, k^2, k^3 of W_5(k), k>=5) and the k=4 defect delta, at ANY precision.
# Quadrature: tanh-sinh with python-flint arb (ball arithmetic), split at v=1, substitutions v=b+s on [b,1], v=a-s on [1,a]
# (so that the endpoint-vanishing factors t, 1-t, v-b, v-a are exact), same as wquad.py.
# usage: id1_blocks.py <digits D> <m>     (step h = 2^-m; writes results/blocks_D<D>_m<m>.json)
# No mpmath is used here, so there is no mp.mp.dps trap; precision is flint ctx.prec (bits), printed below.
import sys, os, json, pickle, time, math
sys.dont_write_bytecode = True
from fractions import Fraction
from flint import arb, ctx

D = int(sys.argv[1]); m = int(sys.argv[2])
bits = int(D * 3.3219281) + 160
ctx.prec = bits
print(f"D={D} m={m} flint ctx.prec={ctx.prec} bits", flush=True)

here = os.path.dirname(os.path.abspath(__file__))
K = pickle.load(open(os.path.join(here, "exact_cubic_kernel_w5.pkl"), "rb"))
blocks = set()
for src in (K["cub"], K["defect"]):
    for (key, pole, comp), c in src.items():
        nz = any(x != 0 for x in c) if isinstance(c, list) else (c != 0)
        if nz:
            blocks.add((key, pole))
blocks = sorted(blocks)
print("blocks:", blocks, flush=True)

sq2 = arb(2).sqrt(); a = 1 + sq2; b = sq2 - 1
pi = arb.pi()
ln10 = arb(10).log()
umax = ((arb(D) * ln10 + 60) / pi).asinh() if hasattr(arb, "asinh") else None
if umax is None:
    x = (arb(D) * ln10 + 60) / pi
    umax = (x + (x * x + 1).sqrt()).log()
h = arb(1) / (2 ** m)
N = int(math.ceil(float(umax.mid()) * 2 ** m)) + 1
print("umax =", float(umax.mid()), "N =", N, flush=True)

acc = {blk: arb(0) for blk in blocks}
t0 = time.time()
for side, Lc in ((1, 1 - b), (2, a - 1)):
    for jj in range(-N, N + 1):
        u = jj * h
        su = u.sinh()
        e = (-pi * su).exp()
        s = Lc / (1 + e)
        w = Lc * pi * u.cosh() * e / (1 + e) ** 2 * h
        if side == 1:
            v = b + s
            t = s * (v + a) / (4 * v)
            l1 = t.log(); l2 = (-t).log1p()
            dv = {"0": v, "a": s - 2, "-a": v + a, "b": s, "-b": v + b}
        else:
            v = a - s
            omt = s * (v + b) / (4 * v)
            l2 = omt.log(); l1 = (-omt).log1p()
            dv = {"0": v, "a": -s, "-a": v + a, "b": 2 - s, "-b": v + b}
        # powers
        p1 = [arb(1)]; p2 = [arb(1)]
        for _ in range(4):
            p1.append(p1[-1] * l1); p2.append(p2[-1] * l2)
        for (key, pole) in blocks:
            acc[(key, pole)] += w * p1[key[0]] * p2[key[1]] / dv[pole]
print(f"done in {time.time()-t0:.1f}s", flush=True)

out = {}
for blk, val in acc.items():
    out[f"{blk[0][0]},{blk[0][1]}|{blk[1]}"] = dict(mid=val.mid().str(D + 10, radius=False), rad=val.rad().str(5))
fn = os.path.join(here, "results", f"blocks_D{D}_m{m}.json")
json.dump(dict(D=D, m=m, bits=bits, blocks=out), open(fn, "w"), indent=1)
for k, v in out.items():
    print(k, v["mid"][:40], "rad", v["rad"], flush=True)
