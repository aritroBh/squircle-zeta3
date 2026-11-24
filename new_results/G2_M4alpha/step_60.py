import sys, json, time
from fractions import Fraction as Fr
cwd = "<project>/research/08_theory/G2_M4alpha"
sys.path.insert(0, cwd)
import g2_fit as gf
ORDER = 1040
T2 = gf.target_series('(1, 0, 0)', ORDER)
for exps in [[-3, -2, -1, 0, 1, 2], [-3, -2, -1, 0, 1], [-2, -1, 0, 1, 2]]:
    t0 = time.time()
    basis = gf.build_even_basis(3, exps, ORDER)
    terms, rankA, rank = gf.solve_exact(basis, T2, ORDER)
    print(f"exps {exps}: {len(basis)} unknowns; rankA={rankA} rankAb={rank}; solution={'YES' if terms else 'NO'} ({time.time()-t0:.0f}s)")
    if terms:
        print(f"  {len(terms)} terms")
        for nm, v in terms:
            print('   ', v, '*', nm)
        break
