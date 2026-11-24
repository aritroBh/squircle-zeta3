
import sys, random
from fractions import Fraction
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
F = W.make_field(13, split=[3], inert=[2])
exps=[4,-5,-6,5,-2]
x=QD(1,0,13)
for g,e in zip(F['letters'],exps):
    if e: x=x*(g**e)
# trace factor_sunit manually with its exact code
N = x.norm()
print("allowed:", W._allowed(N, F['rat_primes']))
y = x
coords = [None]*len(F['letters'])
for i,(g,kind) in enumerate(zip(F['letters'],F['kinds'])):
    if kind in ('ram','split'):
        k,y = W.val_at_g(y, g, F['prime_of'][i], F['erams'][i])
        coords[i]=k
for i,(g,kind) in enumerate(zip(F['letters'],F['kinds'])):
    if kind=='inert':
        p = F['prime_of'][i]
        v = W._v_p(y.norm(), p)
        print("inert v:", v)
        y = y/(Fraction(p)**(v//2))
        coords[i]=v//2
print("y after:", y, "norm:", y.norm())
print("unit_emb:", F['unit_emb'])
old=mp.mp.dps; mp.mp.dps=150
t=mp.log(abs(y.emb(150)))/mp.log(F['unit_emb'])
print("t:", mp.nstr(t,20))
k=int(mp.nint(t))
print("k:",k,"ok:",abs(t-k)<mp.mpf('1e-80'))
mp.mp.dps=old
eps=F['unit']
print("eq:", y==eps**k, y==-(eps**k))
print("direct call:", W.factor_sunit(x,F))
import inspect
print(inspect.getsource(W.factor_sunit))
