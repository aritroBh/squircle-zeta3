import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import engine as E
ok=True
def chk(name, got, want):
    global ok
    good = got==want
    print(('PASS ' if good else 'FAIL '), name, '' if good else (got, want))
    ok = ok and good
# Control B: four pilot c8 vectors
r=E.reduced_coeffs((1,1,1,2),8)
chk('(1,1,1,2) 32c8', [32*x for x in r[8]], [F(-253681),F(208320),F(-33280),F(1856)])
r=E.reduced_coeffs((1,2,2,1),8)
chk('(1,2,2,1) 80c8', [80*x for x in r[8]], [F(-4738167),F(4808064),F(-806400),F(-23424)])
r=E.reduced_coeffs((3,1,2,5),8)
chk('(3,1,2,5) 64c8', [64*x for x in r[8]], [F(-136425803),F(111540096),F(-18257280),F(2348544)])
# Control D: reversal invariance of all weight-8 (and lower) coefficient tuples
for a in [(1,1,1,2),(1,2,2,1),(3,1,2,5),(2,3,1,4),(1,2,3,4),(5,1,1,2)]:
    ra=E.reduced_coeffs(a,8); rb=E.reduced_coeffs(a[::-1],8)
    chk('reversal %s'%str(a), ra, rb)
print('ALL', ok, 't=%.1fs'%(time.time()-t0))
