import sys, importlib, time, pickle, os, json, fractions
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
sys.path.insert(0, "../WT2_weight2_all_k")

import fastred
importlib.reload(fastred)
from fastred import *
import buildintegrand
importlib.reload(buildintegrand)
from buildintegrand import build_integrand_pf

# Define Rw2 for control
# Rw2 = 1/(2v) + 1/(2(v-a)) - 1/(v+a) - 1/(2(v+b))
# where a = 1+sqrt2, b = sqrt2-1
# Poles: 0, a, -a, -b (note: v+b -> pole at -b)
AQ = QSI(1,1,0,0)  # a = 1+sqrt2
BQ = QSI(-1,1,0,0) # b = sqrt2-1

Rw2pf = [{ZERO: QSI(fractions.Fraction(1,2)), 
          AQ: QSI(fractions.Fraction(1,2)), 
          -AQ: QSI(-1),  # -1/(v+a) = coeff -1 at pole -a
          -BQ: QSI(fractions.Fraction(-1,2))},  # -1/[2(v+b)] = coeff -1/2 at pole -b
         [], [], ZERO]

print("Running reduce_k_fast for k=2..14")
print("k  result  (1,0)==Rw2?  higher_parts_zero?  const_zero?")

results = {}
for k in range(2, 15):
    fn = f"binteg_{k}.pkl"
    if not os.path.exists(fn):
        print(f"Warning: {fn} not found")
        continue
    
    with open(fn, 'rb') as f:
        ipf = pickle.load(f)
    
    res = fastred.reduce_k_fast(k, integrand_pf=ipf)
    results[k] = res
    
    # Check (1,0) == Rw2
    d = pf_sub(res[(1,0)], Rw2pf)
    kernel_eq = (not d[0]) and (not d[1]) and (not d[2]) and d[3].is_zero()
    
    # Check all PFs have no higher parts (at +-i) and no const
    all_good = True
    for key, pf in res.items():
        if pf[1] or pf[2]:  # higher_plus_i or higher_minus_i not empty
            all_good = False
            break
        if not pf[3].is_zero():  # const not zero
            all_good = False
            break
        # Check coefficients are in Q(sqrt2) (p2 = p3 = 0)
        for pole, coeff in pf[0].items():
            if coeff.p2 != 0 or coeff.p3 != 0:
                all_good = False
                break
        if not all_good:
            break
    
    print(f"{k:2}  OK      {kernel_eq}            {all_good}")

print("\nAll done")