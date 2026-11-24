import os, sys, json, time
os.nice(10)
from fractions import Fraction as Ff
import mpmath as mp
sys.path.insert(0,'.')
from planner_li5 import li5 as pli5
J = json.load(open('out/d13/candidates_cap.json'))
rel = json.load(open('d13_gcd_minnum.json'))
def parse(s):
    s=s.replace('+-','-'); i=s.rfind('*sqrt13'); head=s[:i]
    k=max(head.rfind('+'), head.rfind('-'))
    return Ff(head[:k]), Ff(head[k:])
args=[parse(k) for k in rel['li5']]
xs=[abs(a/b if False else 0) for a,b in args]
mp.mp.dps=60
def val(ab): 
    a,b=ab; return mp.mpf(a.numerator)/a.denominator + mp.mpf(b.numerator)/b.denominator*mp.sqrt(13)
import numpy as np
mag=sorted(abs(float(val(ab))) for ab in args)
print('n args', len(args), '|x| min/med/max:', mag[0], mag[len(mag)//2], mag[-1])
# time a few at 2200 digits
mp.mp.dps=2260
for ab in [args[0], args[len(args)//2], max(args, key=lambda ab: abs(float(val(ab))))]:
    x=val(ab); t0=time.time(); v=pli5(x, 2210); dt=time.time()-t0
    print(f'x={mp.nstr(x,6)}: {dt:.2f}s')
# Hurwitz L timing
t0=time.time()
def chi13(a):
    a%=13; return 0 if a==0 else (1 if pow(a,6,13)==1 else -1)
L5=mp.fsum(chi13(a)*mp.zeta(5, mp.mpf(a)/13) for a in range(1,13))/mp.mpf(13)**5
print(f'Hurwitz L: {time.time()-t0:.2f}s, L={mp.nstr(L5,20)}')
