import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
mp.mp.dps = 120
d1p = pickle.load(open('dod5_task1.pkl','rb'))
TERMS, vals = d1p['TERMS'], d1p['vals']
def strpair(z):
    m,n,d = z
    if n==0:
        return ("%d/%d"%(m,d) if d!=1 else "%d"%m)
    return "(%d%+d*sqrt3)/%d"%(m,n,d)
print("Task 1 table: 36 constants  f(x,y) = sum_{0<k<m} x^m y^k/(m^2 k)  [prefix conv. = Goncharov Li_{1,2}(y,x)]")
for name, kind, key, xa, ya in TERMS:
    print("  %-14s x=%-14s y=%-14s value=%s"%(name, strpair(xa), strpair(ya), mp.nstr(vals[(kind,key)],18)))
# L(2,chi12) = pi^2 sqrt3 / 18 ?
C = consts(120)
print("\nL(2,chi12) =", mp.nstr(C['L2c'],30))
print("pi^2 r3/18 =", mp.nstr(mp.pi**2*mp.sqrt(3)/18,30))
print("diff:", mp.nstr(abs(C['L2c']-mp.pi**2*mp.sqrt(3)/18),3))
print("r3*L(2,chi12)/zeta(2) =", mp.nstr(mp.sqrt(3)*C['L2c']/C['z2'],10))
