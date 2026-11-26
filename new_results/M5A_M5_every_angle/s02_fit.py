"""s02: FIT (not a proof) of every zeta part of M_5 in terms c^e H_w(+-c), e in -4..1, weight <= 6 - wt(z).
Writes hpl_fit_k5_<part>.json (format of G8: ["c^e He[w]", "q"], plus "1", "c^-2").  Usage: s02_fit.py part [part...]"""
import sys, json, time
sys.setrecursionlimit(20000)
import m5fit as F
PARTS = {'z2^2': ((2,0,0,0), 2), 'z2': ((1,0,0,0), 4), '1': ((0,0,0,0), 6),
         'z3': ((0,1,0,0), 3), 'z4': ((0,0,1,0), 2), 'z5': ((0,0,0,1), 1), 'z2*z3': ((1,1,0,0), 1)}
mu = F.load_mu()
for name in sys.argv[1:]:
    zm, W = PARTS[name]
    t0 = time.time()
    nunk = 2 + 6*sum(3**(L-1) for L in range(1, W+1))
    print(f"part {name}: weight <= {W}, {nunk} unknowns", flush=True)
    sol, cols, info = F.fit(mu[zm], W, range(-4, 2), nunk + 100, nprimes=3)
    if sol is None:
        print("   NO SOLUTION in this basis; rank", info); continue
    rank, nc, piv = info
    print(f"   unique: {rank == nc}; {len(sol)} nonzero terms; {time.time()-t0:.0f}s", flush=True)
    json.dump(F.to_json_terms(sol), open(f"hpl_fit_k5_{name}.json", "w"), indent=0)
