
import os, pickle, time
from fractions import Fraction as Fr
base = "<project>"
wts = base + "/research/08_theory/WTS_weight_k_structure/theory"
K = {k: pickle.load(open(f'{wts}/kern_{k}.pkl','rb')) for k in [5,6,7,8,9,10,11,12,13]}
keys5 = [(5,0),(4,1),(3,2)]
def lag(ys, xs, x):
    out = Fr(0)
    for a,xa in enumerate(xs):
        La = Fr(1)
        for b,xb in enumerate(xs):
            if b!=a: La *= Fr(x-xb, xa-xb)
        out += ys[a]*La
    return out
t0=time.time()
# (a) k=6 vs quartic(7..11)
n=0; bad6=0
for key in keys5:
    poles = set()
    for k in range(6,12): poles |= set(K[k][key]['real'].keys())
    for pole in poles:
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr(0),)*4)[comp] for k in [7,8,9,10,11]]
            n += 1
            if lag(ys,[7,8,9,10,11],6) != K[6][key]['real'].get(pole,(Fr(0),)*4)[comp]: bad6 += 1
print(f"(a) k=6 vs quartic(7..11): {n} components, bad={bad6}")
# (b) quarticity k=7..13 via 5th differences
badq=0; nq=0
for key in keys5:
    poles = set()
    for k in range(7,14): poles |= set(K[k][key]['real'].keys())
    for pole in poles:
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr(0),)*4)[comp] for k in range(7,14)]
            cur = ys
            for _ in range(5): cur = [cur[t+1]-cur[t] for t in range(len(cur)-1)]
            nq += 1
            if any(c!=0 for c in cur): badq += 1
print(f"(b) 5th differences k=7..13 vanish: {nq} components, bad={badq}")
# (c) +/-i simple-pole coefficients vanish k=6..13
nonempty = [(k,key) for k in range(6,14) for key in keys5 if K[k][key]['ipos'] or K[k][key]['ineg']]
print(f"(c) ipos/ineg nonempty at k=6..13: {nonempty}")
# (d) k=5 exceptional
n5=0; bad5=0
for key in keys5:
    poles = set()
    for k in [5,7,8,9,10,11]: poles |= set(K[k][key]['real'].keys())
    for pole in poles:
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr(0),)*4)[comp] for k in [7,8,9,10,11]]
            n5 += 1
            if lag(ys,[7,8,9,10,11],5) != K[5][key]['real'].get(pole,(Fr(0),)*4)[comp]: bad5 += 1
print(f"(d) k=5 vs quartic(7..11): {n5} components, mismatches={bad5} (claim: 15/48 differ)")
print("poles per key:", {key: sorted(set().union(*[set(K[k][key]['real'].keys()) for k in range(6,12)])) for key in keys5})
print("time", round(time.time()-t0,1))
work = base + "/research/08_theory/TOW4R_review"
with open(work+"/progress.log","a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} W6 kernel checks rerun: (a) 48 comp bad={bad6}; (b) quarticity bad={badq}; (c) ipos/ineg all empty k=6..13; (d) k=5 mismatches={bad5}/48\n")
