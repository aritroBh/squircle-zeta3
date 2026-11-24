# Run only after our own coefficients exist: compare with the author's file (different engine).
import json, sys, mpmath as mp
mp.mp.dps = 110
mine = json.load(open(sys.argv[1])); A = json.load(open('../coeffs_K150_dps330.json'))
K = mine["K"]
worst = max(abs(mp.mpf(mine["a"][k]) - mp.mpf(A["a_tau"][k]))/abs(mp.mpf(A["a_tau"][k])) for k in range(1, K+1))
print(f"max rel diff a_k (ours vs author's), k<= {K}: {mp.nstr(worst, 3)}")
