import sys, pickle; fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1)
from fractions import Fraction as Fr
Fd9=pickle.load(open(fd1+'/Fd_data_d9.pkl','rb'))
F9=Fd9[9]
pts=[((1-c*c)/(c*c), v) for c,v in F9.items()]
pts.append((Fr(1),Fr(1)))
# solve for coeffs of degree<=4 poly through all 9 points (exact, least check)
import itertools
def solve_poly(pts, deg):
    # use first deg+1 points to solve, Fraction Gaussian elimination
    base=pts[:deg+1]
    M=[[x**k for k in range(deg,-1,-1)]+[y] for x,y in base]
    n=deg+1
    for col in range(n):
        piv=next(r for r in range(col,n) if M[r][col]!=0)
        M[col],M[piv]=M[piv],M[col]
        pv=M[col][col]
        M[col]=[v/pv for v in M[col]]
        for r in range(n):
            if r!=col and M[r][col]!=0:
                f=M[r][col]
                M[r]=[M[r][i]-f*M[col][i] for i in range(n+1)]
    coeffs=[M[r][n] for r in range(n)]  # descending
    for x,y in pts:
        val=sum(coeffs[i]*x**(deg-i) for i in range(deg+1))
        assert val==y, (x,y,val)
    return coeffs
c4=solve_poly(pts,4)
print("F_9 coeffs descending t^4..t^0:", c4)
from math import lcm, gcd
L=1
for c in c4: L=lcm(L,c.denominator)
ints=[int(c*L) for c in c4]
g=0
for v in ints: g=gcd(g,abs(v))
ints=[v//g for v in ints]
print("integer poly:", ints, " f_9(1)=", sum(ints))
# Also refit F_d for d<=8 in the same convention for a uniform table
Fd=pickle.load(open(fd1+'/Fd_data.pkl','rb'))
for d in range(2,9):
    ptsd=[((1-c*c)/(c*c), v) for c,v in Fd[d].items()]
    ptsd.append((Fr(1),Fr(1)))
    cd=solve_poly(ptsd, d//2)
    L=1
    for cc in cd: L=lcm(L,cc.denominator)
    ii=[int(cc*L) for cc in cd]
    g=0
    for v in ii: g=gcd(g,abs(v))
    ii=[v//g for v in ii]
    print(d, ii, "f(1)=",sum(ii))
open(fd1+'/progress.log','a').write('step_28: F_9 exact integer poly found; uniform table d<=9\n')
