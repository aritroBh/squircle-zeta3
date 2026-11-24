import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
dbg = r'''
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
'''
open('dbg6.py','w').write(dbg)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','dbg6.py'], capture_output=True, text=True, cwd='.')
print(r.stdout); print("ERR:", r.stderr[-800:])
