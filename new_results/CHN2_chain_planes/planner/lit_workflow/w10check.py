import sys, time
sys.path.insert(0,"<project>/research/08_theory/CHN1_chain_bodies")
from fractions import Fraction as F
import engine, mzv_exact as M
print("BASIS[10] =",M.BASIS[10])
def smap(a):
    a1,a2,a3,a4=a; return (a2,F(0),a3,a2+a4,a1+a3)
def run(a,W=10):
    t=time.time(); r=engine.reduced_coeffs(tuple(F(x) for x in a),W); print("  engine",a,"%.1fs"%(time.time()-t)); return r
def Mk_12(a,r,k):
    s=smap(tuple(F(x) for x in a)); return r[k][0]-((s[3]+s[4])**k-s[3]**k-s[4]**k)/k
idx37=[i for i,m in enumerate(M.BASIS[10]) if m==((3,7),)][0]
for a in [(1,1,1,2)]:
    perm=(a[2],a[3],a[0],a[1])
    r=run(a); rp=run(perm)
    M3_12=Mk_12(a,r,3); M3_21=Mk_12(perm,rp,3)
    M7_12=Mk_12(a,r,7); M7_21=Mk_12(perm,rp,7)
    comm=M7_12*M3_21-M3_12*M7_21
    print("a=",a,": engine zeta(3,7) coeff =",r[10][idx37],"  literature -(1/14)[M7,M3]_11 =",-comm/14)
    print("  M7_12=",M7_12,"M7_21=",M7_21,"M3_12=",M3_12,"M3_21=",M3_21)
