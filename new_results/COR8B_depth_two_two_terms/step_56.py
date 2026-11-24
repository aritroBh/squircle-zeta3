
import pickle
from collections import Counter
fs = pickle.load(open('../COR8_depth_two_any_angle/full_sols_v3.pkl','rb'))
alln = fs['alln']; sols = fs['sols']
kinds = Counter(x[0] for x in alln)
print(kinds)
# show the non-product basis elements
for x in alln:
    if x[0] in ('Li4',):
        pass
# used items in the e=-1 solution (nonzero coeffs)
used = [(alln[i], sols[-1][i]) for i in range(len(alln)) if sols[-1][i]]
print("used in e=-1:", len(used))
for u,c0 in used: print("   ", c0, u)
