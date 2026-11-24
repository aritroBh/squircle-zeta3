# Identification step 11: ONE lattice with several targets (c2,c1,c0,delta) against the full weight-5 atom basis (+ L(chi_8,5)), x{1,sqrt2}.
# usage: id11_multi.py <cmfile> <D> [targets...]
import sys, os, json, time
sys.dont_write_bytecode = True
from id_lib import *
cmfile = sys.argv[1]; D = int(sys.argv[2]); tn = sys.argv[3:] or ["c2", "c1", "c0", "delta"]
R = json.load(open(cmfile))
setprec(R["D"])
C = Consts()
atoms = atoms5(C, chi8=True)
tv = {"c0": arb(R["c"][0]), "c1": arb(R["c"][1]), "c2": arb(R["c"][2]), "c3": arb(R["c"][3]), "delta": arb(R["delta"])}
targets = [(t, tv[t]) for t in tn]
print(f"flint ctx.prec={ctx.prec} bits; D={D}; targets={tn}; atoms={len(atoms)}; total dims={len(atoms)+len(targets)}", flush=True)
allv = [v for _, v in targets] + [v for _, v in atoms]
names = [t for t, _ in targets] + [n for n, _ in atoms]
t0 = time.time()
n = len(allv)
M = fmpz_mat(n, n + 1)
for i in range(n):
    M[i, i] = 1; M[i, n] = to_int(allv[i], D)
R2 = M.lll(gram="exact")
print(f"LLL done {time.time()-t0:.0f}s; noise floor ~ {10**(D/n):.2e}", flush=True)
rows = []
for r in range(R2.nrows()):
    row = [int(R2[r, c]) for c in range(n + 1)]
    mx = max(abs(x) for x in row[:n]); rows.append((mx, row))
rows.sort(key=lambda t: t[0])
for mx, row in rows[:6]:
    nz = {nm: c for nm, c in zip(names, row[:n]) if c != 0}
    print(f"row maxcoef={mx:.3e} last={row[n]} ntargets_nonzero={sum(1 for t in tn if row[names.index(t)] != 0)} nterms={len(nz)}", flush=True)
    if mx < 10 ** (D / n - 4): print("   RELATION:", nz, flush=True)
json.dump([dict(maxcoef=str(mx), row=row) for mx, row in rows[:6]], open(f"results/id11_D{D}_{'_'.join(tn)}.json", "w"))
