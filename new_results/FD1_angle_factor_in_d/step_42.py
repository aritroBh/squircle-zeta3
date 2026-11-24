import sys, pickle; fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1)
from fractions import Fraction as Fr
Fd10=pickle.load(open(fd1+'/Fd_data_d10.pkl','rb'))
F10=Fd10[10]
pts=[((1-c*c)/(c*c), v) for c,v in F10.items()]   # 8 corner points only
def solve_poly(pts, deg, check=True):
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
    bad=[]
    for x,y in pts:
        if sum(coeffs[i]*x**(deg-i) for i in range(deg+1))!=y: bad.append((x,y))
    return coeffs, bad
c10,bad=solve_poly(pts,5)
print("fit through 6 corner pts; failures on remaining:", bad)
val1=sum(c10[i] for i in range(6))
print("value at t=1:", val1)
# try all subsets of 6 to check consistency
import itertools
vals=set()
for sub in itertools.combinations(range(8),6):
    cf,_=solve_poly([pts[i] for i in sub],5)
    vals.add(sum(cf[i] for i in range(6)))
print("distinct t=1 values over C(8,6) subsets:", vals)
