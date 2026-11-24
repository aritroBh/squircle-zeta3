# Identification step 9: diagnostics on natural sub-combinations of the block integrals (pole-0 combinations etc.)
# usage: id9_diag.py <D> <flags...> -- <target names...>
import sys, os, json, time
sys.dont_write_bytecode = True
from id_lib import *
args = sys.argv[1:]
D = int(args[0]); i = args.index("--"); flags = args[1:i]; tnames = args[i+1:]
B = json.load(open("results/blocks_D3300_m11.json"))["blocks"]
setprec(3300)
C = Consts()
def blk(k, p): return arb(B[f"{k[0]},{k[1]}|{p}"]["mid"])
T = {
 "P0": 3*blk((2,2),"0") - 4*blk((3,1),"0") + blk((4,0),"0"),
 "Q0": blk((2,2),"0") - blk((3,1),"0"),
 "A22_0": blk((2,2),"0"), "A31_0": blk((3,1),"0"), "A40_0": blk((4,0),"0"),
 "A22_ma": blk((2,2),"-a"), "A31_ma": blk((3,1),"-a"), "A40_ma": blk((4,0),"-a"),
 "A31_mb": blk((3,1),"-b"), "A40_mb": blk((4,0),"-b"), "A31_a": blk((3,1),"a"), "A40_a": blk((4,0),"a"),
}
atoms = atoms5(C, with_sqrt2="nosqrt2" not in flags, li4=("half" if "li4half" in flags else "noli4" not in flags), li5=("half" if "li5half" in flags else "noli5" not in flags), li3="noli3" not in flags, li2="noli2" not in flags)
print(f"flint ctx.prec={ctx.prec}; D={D}; atoms={len(atoms)}; flags={flags}", flush=True)
names = ["target"] + [n for n, _ in atoms]
for tn in tnames:
    t0 = time.time()
    best, top = lll_relation([T[tn]] + [v for _, v in atoms], D)
    coeffs, last, mx = best
    n = len(atoms) + 1
    ok = coeffs[0] != 0 and mx < 10 ** (D / n - 6)
    print(f"{tn:8s} maxcoef={mx:.3e} floor~{10**(D/n):.2e} last={last} -> {'RELATION' if ok else 'none'} ({time.time()-t0:.0f}s)", flush=True)
    if ok: print("   ", {nm: c for nm, c in zip(names, coeffs) if c != 0}, flush=True)
