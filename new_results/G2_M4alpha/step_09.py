"""step_09.py -- the rational part of B4: fit in the even-sector basis, weight <= 5 HPLs.

Basis: c^e (H_w +/- H_{-w}) for orbits {w,-w} of words over {0,1,-1} with nonzero last letter,
weight <= 5 (121 orbits), e in {-3,-2,-1,0,1}, plus the constant: 606 unknowns.
The exponent range is forced by the denominator analysis (step_08 output: denominator factors
(2n-1),(2n+1),(2n+3),n,(n+1) only, mapping to e in {-3..1}).
Exact rational linear algebra on Taylor coefficients 0..1280 (641 even equations);
accepted only if every equation is satisfied. Writes hpl_fit_k4_(0,0,0).json.
"""
import json, time
import g2_fit as gf

ORDER = 1280
t0 = time.time()
T = gf.target_series('(0, 0, 0)', ORDER, tablefile='flat_exact_N640_K4_formula.json')
print(f'target built ({time.time()-t0:.0f}s)', flush=True)
basis = gf.build_even_basis(5, [-3, -2, -1, 0, 1], ORDER)
print(f'basis built: {len(basis)} functions ({time.time()-t0:.0f}s)', flush=True)
terms, rankA, rank = gf.solve_exact(basis, T, ORDER)
print(f'({time.time()-t0:.0f}s)')
if terms is None:
    print('NO SOLUTION in weight<=5, e in {-3..1} basis')
else:
    print(f'SOLUTION with {len(terms)} terms (rank {rankA} of {len(basis)}):')
    for nm, v in terms:
        print(f'   {v}  *  {nm}')
    json.dump([(nm, str(v)) for nm, v in terms], open('hpl_fit_k4_(0,0,0).json', 'w'))
    print('saved hpl_fit_k4_(0,0,0).json')
