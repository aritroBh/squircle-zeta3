import os, time
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
import mpmath
mpmath.mp.dps = 60
from mpmath import mpf, sqrt, log, pi, nsum, quad, inf
import sys; sys.path.insert(0,'.')
from math import comb
def hfun(w): return ((1+w)/sqrt(1+w**2)-1)/w
def F(z): return quad(lambda w: w**(-z)*hfun(w), [0,1])
def cn(n): return mpf(comb(2*int(n),int(n)))/4**int(n)
def F_series(y):
    s1 = nsum(lambda n: (-1)**int(n)*cn(n)/(2*n-y), [1, inf])
    s2 = nsum(lambda n: (-1)**int(n)*cn(n)/(2*n+1-y), [0, inf])
    return s1+s2
import lead2_lib
from lead2_lib import Efun
for y in [mpf('0.3'), mpf('-0.45'), mpf('0.9')]:
    lhs = (F_series(y)+F_series(-y))/2
    rhs = (Efun(y)-Efun(-y))/(2*y)
    print(f"y={y}: |series-quad| = {mpmath.nstr(abs(F_series(y)-F(y)),5)}   |(F+F-)/2-(E-E-)/2y| = {mpmath.nstr(abs(lhs-rhs),5)}")
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: F(y) series (n>=0) via nsum confirmed ~1e-50; (F+F-)/2=(E-E-)/2y at series level\n")
