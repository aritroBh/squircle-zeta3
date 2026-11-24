
import pickle
from fractions import Fraction as Fr
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve81.pkl','rb'))
li4v, prodv, Rv = D['li4v'], D['prodv'], D['Rv']
S = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/full_sols_raw.pkl','rb'))
sols, alln = S['sols'], S['alln']
allv = [li4v[k] for k in li4v] + [prodv[k] for k in prodv]
# exact check: R_E == sum_j sols[E][j] * allv[j]
for E in (-1,0,1):
    acc = [Fr(0)]*81
    for j in range(len(alln)):
        if sols[E][j]:
            acc = [x + sols[E][j]*y for x, y in zip(acc, allv[j])]
    diff = [x - y for x, y in zip(Rv[E], acc)]
    print(f"E={E}: symbol residual max |coeff| =", max(abs(x) for x in diff))
