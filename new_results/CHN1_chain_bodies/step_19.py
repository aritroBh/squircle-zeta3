import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import engine as E
c=E.coeffs((1,1,1,1),11)
print('c1 empty:', c[1]=={})
r=E.reduced_coeffs((1,1,1,1),8)
ok=True
def chk(name, got, want):
    global ok
    good = got==want
    print(('PASS ' if good else 'FAIL '), name, got if not good else '')
    ok = ok and good
chk('c2=[-3]', r[2], [F(-3)])
chk('c3=[10]', r[3], [F(10)])
chk('c4=[-19]', r[4], [F(-19)])
chk('c5=[102,-36]', r[5], [F(102),F(-36)])
chk('c6=[-159,68]', r[6], [F(-159),F(68)])
chk('3c8=[-4681,4680,-792,0]', [3*x for x in r[8]], [F(-4681),F(4680),F(-792),F(0)])
print('z(3,5) coeff on diagonal:', r[8][3])
print('ALL', ok, 't=%.1fs'%(time.time()-t0))
