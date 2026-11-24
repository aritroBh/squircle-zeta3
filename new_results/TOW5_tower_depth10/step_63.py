import pickle, time
from fractions import Fraction as Fr
t0=time.time()
work = '<project>/research/08_theory/TOW5_tower_depth10'
wts = '<project>/research/08_theory/WTS_weight_k_structure/theory'
K = {k: pickle.load(open(f'{wts}/kern_{k}.pkl','rb')) for k in range(6, 14)}
keys6 = [(6,0),(5,1),(4,2),(3,3)]
# verify all degree-6 keys present at k=7..13 and k=6
for k in range(6, 14):
    missing = [key for key in keys6 if key not in K[k]]
    if missing: print(f"k={k}: MISSING {missing}")
print("all degree-6 keys present k=6..13")

def lag(ys, xs, x):
    out = Fr(0)
    for a,xa in enumerate(xs):
        La = Fr(1)
        for b,xb in enumerate(xs):
            if b!=a: La *= Fr(x-xb, xa-xb)
        out += ys[a]*La
    return out

# 1. quintic fitted from k=8..13, check k=7 matches
w7a = True; ncomp = 0; badex=[]
for key in keys6:
    poles = set()
    for k in range(7,14): poles |= set(K[k][key]['real'].keys())
    for pole in poles:
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr(0),)*4)[comp] for k in [8,9,10,11,12,13]]
            ncomp += 1
            if lag(ys,[8,9,10,11,12,13],7) != K[7][key]['real'].get(pole,(Fr(0),)*4)[comp]:
                w7a=False; badex.append((key,pole,comp))
print(f"W7a: k=7 == quintic(k=8..13): {ncomp} components, bad = {ncomp - sum(1 for _ in [1]) if not w7a else 0}; PASS={w7a} {badex[:5]}")

# 2. quinticity on k=7..13: 6th differences vanish
w7b = True; nb=0
for key in keys6:
    poles = set()
    for k in range(7,14): poles |= set(K[k][key]['real'].keys())
    for pole in poles:
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr(0),)*4)[comp] for k in range(7,14)]
            cur = ys
            for _ in range(6): cur = [cur[t+1]-cur[t] for t in range(len(cur)-1)]
            nb+=1
            if any(c!=0 for c in cur): w7b=False
print(f"W7b: 6th differences vanish on k=7..13 ({nb} comp): PASS={w7b}")

# 3. +-i simple-pole coefficients vanish k=7..13
w7c = all(not K[k][key]['ipos'] and not K[k][key]['ineg'] for k in range(7,14) for key in keys6)
print(f"W7c: +-i simple-pole coefficients vanish k=7..13: PASS={w7c}")

# 4. k=6 vs the quintic: count mismatches
w6mis = 0; w6tot = 0
for key in keys6:
    poles = set()
    for k in range(7,14): poles |= set(K[k][key]['real'].keys())
    poles |= set(K[6][key]['real'].keys())
    for pole in poles:
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr(0),)*4)[comp] for k in [8,9,10,11,12,13]]
            w6tot += 1
            if lag(ys,[8,9,10,11,12,13],6) != K[6][key]['real'].get(pole,(Fr(0),)*4)[comp]:
                w6mis += 1
print(f"W7d: k=6 vs quintic(k=8..13): {w6tot} components, {w6mis} mismatched (k=6 initial kernel)")
print("time", round(time.time()-t0,1))
