import sys
sys.path.insert(0,"<project>/research/08_theory/CHN1_chain_bodies")
from fractions import Fraction as F
import engine, mzv_exact as M

def smap(a):
    a1,a2,a3,a4=a
    return (a2,F(0),a3,a2+a4,a1+a3)
def m11(s):
    s1,s2,s3,s4,s5=s
    return s3*(-s1*(s1+2*s2+s3)+s3*s4+s4**2)+s1*s5*(s1+s5)
def P2_11(s):
    s1,s2,s3,s4,s5=s; return -s3*s4+s1*(s3-s5)
def m12(s):
    s1,s2,s3,s4,s5=s; s13=s4-s1-s2; s24=s5-s2-s3
    return -s13*s24*(s1+s2+s3+s4+s5)
def m21(s):
    s1,s2,s3,s4,s5=s; return s1*s3*(s1+s2+s3-2*(s4+s5))

cache={}
def eng(a):
    a=tuple(F(x) for x in a)
    if a not in cache: cache[a]=engine.reduced_coeffs(a,8)
    return cache[a]

def M3_11(a): return eng(a)[3][0]
def M5_11(a): return eng(a)[5][0]
def M3_12(a):
    s=smap(a); return M3_11(a)-s[3]*s[4]*(s[3]+s[4])          # 4pt shifted term at W=3
def M5_12(a):
    s=smap(a); return M5_11(a)-((s[3]+s[4])**5-s[3]**5-s[4]**5)/5
def perm(a): return (a[2],a[3],a[0],a[1])   # relabeling 2<->3 on slice s23=0
def M3_21(a): return M3_12(perm(a))
def M5_21(a): return M5_12(perm(a))

# convention check of engine's (3,5)
T=engine.table(8)
print("engine reduce_word((5,3)) on [z8, z3z5, z2z3^2, z(3,5)]:",T.reduce_word((5,3)))
print("engine reduce_word((3,5)):",T.reduce_word((3,5)))

pts=[(1,1,1,2),(1,2,2,1),(3,1,2,5),(2,1,1,3),(1,2,3,4),(2,3,1,5)]
for a in pts:
    a=tuple(F(x) for x in a); s=smap(a)
    e=eng(a)
    ok2 = e[2][0]==P2_11(s)
    ok3 = e[3][0]==m11(s)
    ok3_12 = M3_12(a)==m12(s)
    ok3_21 = M3_21(a)==m21(s)
    comm=M5_12(a)*M3_21(a)-M3_12(a)*M5_21(a)
    lit_ours=-comm/5          # coefficient of zeta(3,5)_ours (n>m, 3 on larger index)
    lit_bss = comm/5
    eng35=e[8][3]
    P8=a[0]*a[1]*a[2]*a[3]*sum(a)*(2*(a[0]+a[3]-a[1]-a[2])*(2*(a[0]+a[3])+(a[1]+a[2]))*((a[0]+a[3])+2*(a[1]+a[2]))
         -9*a[0]*a[3]*((a[0]+a[3])+2*(a[1]+a[2]))+9*a[1]*a[2]*(2*(a[0]+a[3])+(a[1]+a[2])))/10
    print("a=",tuple(map(int,a)),"| P2_11 ok",ok2,"| M3_11 ok",ok3,"| M3_12 ok",ok3_12,"| M3_21 ok",ok3_21)
    print("    engine zeta(3,5) coeff =",eng35,"  P8 closed form =",P8)
    print("    literature -(1/5)[M5,M3]_11 =",lit_ours,"   +(1/5)[M5,M3]_11 =",lit_bss)
    print("    M5_12 =",M5_12(a),"M5_21 =",M5_21(a),"M3_12 =",M3_12(a),"M3_21 =",M3_21(a))
