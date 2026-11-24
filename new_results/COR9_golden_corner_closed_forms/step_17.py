import sys, time
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
sys.path.insert(0, WD)
from cor9_lib import fund_unit, class_number, pell_min, h_plus, disc_F
for D in [2,3,5,6,7,14,46,61,94]:
    t0=time.time()
    print(D, fund_unit(D), 'h =', class_number(D), f'{time.time()-t0:.2f}s')
