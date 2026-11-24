import pickle, os
from fractions import Fraction as Fr
base = '<project>/research/08_theory'
wts = f'{base}/WTS_weight_k_structure/theory'
K = {k: pickle.load(open(f'{wts}/kern_{k}.pkl','rb')) for k in range(4, 14)}
def residues(kern, key):
    return kern[key]['real']
keys4 = [(4,0),(3,1),(2,2)]
# cubic fit from k=6..9, predict k=5; compare
def fit_cubic_eval(ys, x):
    # ys at x=6,7,8,9; Lagrange at x
    xs = [6,7,8,9]
    out = Fr(0)
    for a, xa in enumerate(xs):
        La = Fr(1)
        for b, xb in enumerate(xs):
            if b != a: La *= Fr(x-xb, xa-xb)
        out += ys[a]*La
    return out
bad = 0; tested = 0; fails = []
for key in keys4:
    poles = set()
    for k in range(5,10):
        poles |= set(residues(K[k], key).keys())
    for pole in poles:
        for comp in range(4):
            ys = [residues(K[k], key).get(pole, (Fr(0),)*4)[comp] for k in [6,7,8,9]]
            pred5 = fit_cubic_eval(ys, 5)
            act5 = residues(K[5], key).get(pole, (Fr(0),)*4)[comp]
            tested += 1
            if pred5 != act5:
                bad += 1; fails.append((key,pole,comp,pred5,act5))
print(f'degree-4 kernels: k=5 vs cubic(k=6..9): {tested} residue-components, {len(fails)} bad')
for f_ in fails[:8]: print(f_)
# also confirm the cubic law on k=6..13 (order-4 differences vanish) as empirical confirmation of K(4) for k>=6
badc = 0; testc = 0
for key in keys4:
    poles = set()
    for k in range(6,14): poles |= set(residues(K[k], key).keys())
    for pole in poles:
        for comp in range(4):
            ys = [residues(K[k], key).get(pole, (Fr(0),)*4)[comp] for k in range(6,14)]
            cur = ys
            for _ in range(4): cur = [cur[t+1]-cur[t] for t in range(len(cur)-1)]
            testc += 1
            if any(c != 0 for c in cur): badc += 1
print(f'cubicity k=6..13 (4th differences vanish): {testc} components, {badc} bad')
