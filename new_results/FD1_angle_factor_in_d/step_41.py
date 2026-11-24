import sys, pickle; fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1)
from fractions import Fraction as Fr
Fd10=pickle.load(open(fd1+'/Fd_data_d10.pkl','rb'))
F10=Fd10[10]
pts=[((1-c*c)/(c*c), v) for c,v in F10.items()]
pts.append((Fr(1),Fr(1)))
def solve_poly(pts, deg):
    base=pts[:deg+1]
    M=[[x**k for k in range(deg,-1,-1)]+[y] for x,y in base]
    n=deg+1
    for col in range(n):
        piv=next(r for r in range(col,n) if M[r][col]!=0)
        M[col],M[piv]=M[piv],M[col]
        pv=M[col][col]; M[col]=[v/pv for v in M[col]]
        for r in range(n):
            if r!=col and M[r][col]!=0:
                f=M[r][col]; M[r]=[M[r][i]-f*M[col][i] for i in range(n+1)]
    coeffs=[M[r][n] for r in range(n)]
    for x,y in pts:
        assert sum(coeffs[i]*x**(deg-i) for i in range(deg+1))==y,(x,y)
    return coeffs
c10=solve_poly(pts,5)
from math import lcm, gcd
L=1
for cc in c10: L=lcm(L,cc.denominator)
ints=[int(cc*L) for cc in c10]
g=0
for v in ints: g=gcd(g,abs(v))
ints=[v//g for v in ints]
print("F_10 integer poly (t^5..t^0):", ints, " f_10(1)=",sum(ints))
print("fits all 9 points with degree 5: OK (3 spares)")
# check LC pattern: 10!/2^(10-s2(10)) = 3628800/2^8 = 14175
print("LC predicted 14175, actual", ints[0])
open(fd1+'/progress.log','a').write('step_40: F_10 exact integer poly fit, verified on 9 pts\n')
