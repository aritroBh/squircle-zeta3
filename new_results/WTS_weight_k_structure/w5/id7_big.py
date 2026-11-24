# Identification step 7: big LLL for c_2 (and others) in the weight-5 atom basis (Q(sqrt2)-span of products + classical polylogs).
# usage: id7_big.py <cmfile> <D> <target> [flags: nosqrt2 noli4 noli5 noli3 noli2]
import sys, os, json, time
sys.dont_write_bytecode = True
from id_lib import *
cmfile = sys.argv[1]; D = int(sys.argv[2]); tname = sys.argv[3]; flags = sys.argv[4:]
R = json.load(open(cmfile))
setprec(R["D"])
C = Consts()
atoms = atoms5(C, with_sqrt2="nosqrt2" not in flags, li4=("half" if "li4half" in flags else "noli4" not in flags), li5=("half" if "li5half" in flags else "noli5" not in flags), li3="noli3" not in flags, li2="noli2" not in flags)
x = {"c0": arb(R["c"][0]), "c1": arb(R["c"][1]), "c2": arb(R["c"][2]), "c3": arb(R["c"][3]), "delta": arb(R["delta"])}[tname]
print(f"flint ctx.prec={ctx.prec} bits; D={D}; target {tname}; atoms={len(atoms)}; flags={flags}", flush=True)
t0 = time.time()
best, top = lll_relation([x] + [v for _, v in atoms], D)
coeffs, last, mx = best
n = len(atoms) + 1
floor = 10 ** (D / n)
print(f"LLL done {time.time()-t0:.0f}s; best row: maxcoef={mx:.4e} (noise floor ~ {floor:.2e}), last entry {last}", flush=True)
names = ["target"] + [nm for nm, _ in atoms]
print("nonzero:", {nm: c for nm, c in zip(names, coeffs) if c != 0}, flush=True)
for row in top[1:]:
    print("next row maxcoef", max(abs(c) for c in row[:-1]), "target coef", row[0], "last", row[-1], flush=True)
json.dump(dict(target=tname, D=D, flags=flags, coeffs={nm: c for nm, c in zip(names, coeffs) if c != 0}, maxcoef=str(mx), last=str(last)), open(f"results/id7_{tname}_D{D}_{'_'.join(flags) or 'full'}.json", "w"), indent=1)
