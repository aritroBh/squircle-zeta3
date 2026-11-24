# Identification step 2: assemble c_0..c_3 and the k=4 defect delta from a blocks_D<D>_m<m>.json file (exact rational kernels x block integrals).
# usage: id2_assemble.py <blocksfile> ; writes results/cm_<tag>.json (D+10 digit strings)
import sys, os, json, pickle
sys.dont_write_bytecode = True
from flint import arb, ctx
fn = sys.argv[1]
B = json.load(open(fn))
D = B["D"]; bits = B["bits"]
ctx.prec = bits
print("flint ctx.prec =", ctx.prec, "bits; D =", D)
here = os.path.dirname(os.path.abspath(__file__))
K = pickle.load(open(os.path.join(here, "exact_cubic_kernel_w5.pkl"), "rb"))
sq2 = arb(2).sqrt()
def blk(key, pole):
    d = B["blocks"][f"{key[0]},{key[1]}|{pole}"]
    return arb(d["mid"])
def fr(x): return arb(x.numerator) / arb(x.denominator)
cm = [arb(0)] * 4
delta = arb(0)
for (key, pole, comp), cs in K["cub"].items():
    if not any(c != 0 for c in cs): continue
    bv = blk(key, pole) * (1 if comp == "p0" else sq2)
    for mm in range(4): cm[mm] += fr(cs[mm]) * bv
for (key, pole, comp), dv in K["defect"].items():
    bv = blk(key, pole) * (1 if comp == "p0" else sq2)
    delta += fr(dv) * bv
out = dict(D=D, source=os.path.basename(fn), c=[x.mid().str(D, radius=False) for x in cm], delta=delta.mid().str(D, radius=False))
tag = os.path.basename(fn).replace("blocks_", "").replace(".json", "")
json.dump(out, open(os.path.join(here, "results", f"cm_{tag}.json"), "w"), indent=1)
for mm in range(4): print(f"c_{mm} =", out["c"][mm][:60])
print("delta =", out["delta"][:60])
