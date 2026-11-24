import pickle, os
from fractions import Fraction as Fr
base = '<project>/research/08_theory'
wts = f'{base}/WTS_weight_k_structure/theory'
K = {k: pickle.load(open(f'{wts}/kern_{k}.pkl','rb')) for k in range(6, 12)}
keys5 = [(5,0),(4,1),(3,2)]
# sanity: keys present at k=6..11?
for k in range(6,12):
    missing = [key for key in keys5 if key not in K[k]]
    print(f"k={k}: keys5 missing {missing}; ipos/ineg for keys5:", {key: (len(K[k][key]['ipos']), len(K[k][key]['ineg'])) for key in keys5 if key in K[k]})
# pole union for keys5 across k=6..11
poles = {}
for key in keys5:
    ps = set()
    for k in range(6,12):
        ps |= set(K[k][key]['real'].keys())
    poles[key] = ps
    print(key, "poles:", sorted(ps), "const nonzero?", any(K[k][key]['const'] != (Fr(0),)*4 for k in range(6,12)))
