import pickle, os
from fractions import Fraction as Fr
base = '<project>/research/08_theory'
wts = f'{base}/WTS_weight_k_structure/theory'
K = {k: pickle.load(open(f'{wts}/kern_{k}.pkl','rb')) for k in range(5, 12)}
keys5 = [(5,0),(4,1),(3,2)]
def lag_eval(ys, xs, x):
    out = Fr(0)
    for a, xa in enumerate(xs):
        La = Fr(1)
        for b, xb in enumerate(xs):
            if b != a: La *= Fr(x-xb, xa-xb)
        out += ys[a]*La
    return out
print("kern_5 has keys5:", [key for key in keys5 if key in K[5]])
# does kern_5 degree-5 residues match the quartic from k=7..11? (expected NOT: k=5 exceptional)
tested=0; mismatch=0
for key in keys5:
    if key not in K[5]: continue
    poles = set()
    for k in range(5,12): poles |= set(K[k][key]['real'].keys())
    for pole in poles:
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr(0),)*4)[comp] for k in [7,8,9,10,11]]
            pred5 = lag_eval(ys, [7,8,9,10,11], 5)
            act5 = K[5][key]['real'].get(pole,(Fr(0),)*4)[comp]
            tested += 1
            if pred5 != act5: mismatch += 1
print(f'k=5 vs quartic(k=7..11): {tested} components, {mismatch} MISMATCHES (mismatch>0 confirms k=5 exceptional)')
