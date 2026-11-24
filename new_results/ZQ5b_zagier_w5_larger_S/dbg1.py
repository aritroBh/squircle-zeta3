
import sys, random
from fractions import Fraction
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
# debug one failing D=13 case
F = W.make_field(13, split=[3], inert=[2])
exps=[4, -5, -6, 5, -2]
x=QD(1,0,13)
for g,e in zip(F['letters'],exps):
    if e: x=x*(g**e)
print("x =", x, " norm:", x.norm())
N=x.norm()
print("allowed:", W._allowed(N, F['rat_primes']))
# manual step by step
y=x
coords=[None]*5
for i,(g,kind) in enumerate(zip(F['letters'],F['kinds'])):
    if kind in ('ram','split'):
        k,y = W.val_at_g(y, g, F['prime_of'][i], F['erams'][i])
        coords[i]=k
        print(f"strip {F['names'][i]}: k={k}  y={y}")
for i,(g,kind) in enumerate(zip(F['letters'],F['kinds'])):
    if kind=='inert':
        p=F['prime_of'][i]
        v=W._v_p(y.norm(),p)
        print(f"inert {p}: v_2(norm)={v}, y.norm()={y.norm()}")
        y=y/(Fraction(p)**(v//2))
        coords[i]=v//2
print("after stripping: y =", y, " norm:", y.norm())
mp.mp.dps=150
t=mp.log(abs(y.emb(150)))/mp.log(F['unit_emb'])
print("log ratio:", mp.nstr(t,10))
eps=F['unit']
k=int(mp.nint(t))
print("k:",k, "y==eps^k:", y==eps**k, " y==-eps^k:", y==-(eps**k))
print("eps^k =", eps**k)
