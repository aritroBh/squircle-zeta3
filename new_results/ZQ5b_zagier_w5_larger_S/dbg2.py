
import sys, random
from fractions import Fraction
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
F = W.make_field(13, split=[3], inert=[2])
# case 1: positive
exps=[4,-5,-6,5,-2]
x=QD(1,0,13)
for g,e in zip(F['letters'],exps):
    if e: x=x*(g**e)
print("pos:", W.factor_sunit(x,F))
print("neg:", W.factor_sunit(-x,F))
# rerun the actual test loop with its seed
random.seed(42)
fails=[]
for t in range(30):
    exps=[random.randint(-6,6) for _ in F['letters']]
    x=QD(1,0,13)
    for g,e in zip(F['letters'],exps):
        if e: x=x*(g**e)
    sgn=random.choice([1,-1]); x = x if sgn==1 else -x
    f=W.factor_sunit(x,F)
    if f!=exps: fails.append((exps,sgn,f))
print("fails:", fails[:5], "total", len(fails))
# inspect the first actual failing case in detail
if fails:
    exps, sgn, f = fails[0]
    print("first failing exps:", exps, "sgn:", sgn)
    x=QD(1,0,13)
    for g,e in zip(F['letters'],exps):
        if e: x=x*(g**e)
    x = x if sgn==1 else -x
    y=x
    coords=[None]*5
    for i,(g,kind) in enumerate(zip(F['letters'],F['kinds'])):
        if kind in ('ram','split'):
            k,y = W.val_at_g(y, g, F['prime_of'][i], F['erams'][i])
            coords[i]=k
            print(f"  strip {F['names'][i]}: k={k} y={y}")
