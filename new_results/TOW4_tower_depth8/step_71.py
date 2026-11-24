import pickle, os, time
from fractions import Fraction as Fr
base = '<project>/research/08_theory'
wts = f'{base}/WTS_weight_k_structure/theory'
work = base + '/TOW4_tower_depth8'
K = {k: pickle.load(open(f'{wts}/kern_{k}.pkl','rb')) for k in range(6, 14)}
keys5 = [(5,0),(4,1),(3,2)]
def lag_eval(ys, xs, x):
    out = Fr(0)
    for a, xa in enumerate(xs):
        La = Fr(1)
        for b, xb in enumerate(xs):
            if b != a: La *= Fr(x-xb, xa-xb)
        out += ys[a]*La
    return out

# 1) k=6 vs quartic(k=7..11), all real poles, 4 QSI comps
tested=0; bad=0; fails=[]
for key in keys5:
    poles = set()
    for k in range(6,12): poles |= set(K[k][key]['real'].keys())
    for pole in poles:
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr(0),)*4)[comp] for k in [7,8,9,10,11]]
            pred6 = lag_eval(ys, [7,8,9,10,11], 6)
            act6 = K[6][key]['real'].get(pole,(Fr(0),)*4)[comp]
            tested += 1
            if pred6 != act6:
                bad += 1; fails.append((key,pole,comp,pred6,act6))
print(f'W6 Task C: degree-5 kernels k=6 vs quartic(k=7..11): {tested} components, {bad} bad')
for f in fails[:6]: print(' FAIL', f)

# 2) quarticity on k=7..13: 5th differences vanish
testc=0; badc=0
for key in keys5:
    poles = set()
    for k in range(7,14): poles |= set(K[k][key]['real'].keys())
    for pole in poles:
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr(0),)*4)[comp] for k in range(7,14)]
            cur = ys
            for _ in range(5): cur = [cur[t+1]-cur[t] for t in range(len(cur)-1)]
            testc += 1
            if any(c != 0 for c in cur): badc += 1
print(f'quarticity k=7..13 (5th differences vanish): {testc} components, {badc} bad')

# 3) ipos/ineg vanish k=6..13
badip = sum(1 for k in range(6,14) for key in keys5 if K[k][key]['ipos'] or K[k][key]['ineg'])
print(f'+/-i simple poles present at any k=6..13: {badip} (0 = all vanish)')
# log
log = open(work + '/progress.log','a')
log.write(time.strftime('%H:%M:%S ') + f"TASK C: degree-5 kernels k=6 vs quartic(k=7..11): {tested} components, {bad} bad; quarticity k=7..13 5th-diff vanish {testc} comp {badc} bad; +/-i poles vanish k=6..13 (all 3 keys).\n")
log.close()
