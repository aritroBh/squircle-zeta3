
from fractions import Fraction
import mpmath as mp
import sys
sys.path.insert(0,'.')
from w5blib import QD
eps = QD(Fraction(3,2),Fraction(1,2),13)
v1 = eps.emb(150)          # no abs
v2 = abs(eps.emb(150))     # abs at dps=15
v3 = eps.emb(150); v3 = +v3
mp.mp.dps=150
w = abs(eps.emb(150))
print("v1 diff:", mp.nstr(abs(v1-w),5))
print("v2 diff:", mp.nstr(abs(v2-w),5))
# check mantissa bit counts
from mpmath import libmp
print("v1 bc:", v1._mpf_[3] if hasattr(v1,'_mpf_') else '?')
print("v2 bc:", v2._mpf_[3] if hasattr(v2,'_mpf_') else '?')
