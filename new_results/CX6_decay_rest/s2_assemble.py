import mpmath as mp
from common import *
mp.mp.dps=40
M=load()
eps=mp.mpf(0)
def binom_coeffs(m):
    # coefficients in eps of binom(-1-eps,m) = (-1)^m (1+eps)...(m+eps)/m!
    c=[mp.mpf(1)]
    for r in range(1,m+1):
        c=[ (c[i] if i<len(c) else 0)*r + (c[i-1] if i>=1 else 0) for i in range(len(c)+1)]
    return [(-1)**m*x/mp.factorial(m) for x in c]
for k,mmax in [(10,12),(100,8),(250,6),(400,5)]:
    tot=0; parts=[]
    for m in range(mmax+1):
        bc=binom_coeffs(m); pm=0
        for i in range(0,min(m,k-1)+1):
            pm+= bc[i]*d_cut(k-i,m)
        parts.append(pm); tot+=pm
    Mk=-tot
    print(f"k={k}: assembled={mp.nstr(Mk,25)} table={mp.nstr(M[k],25)} rel={mp.nstr((Mk-M[k])/M[k],3)}")
    print("   m-parts/total:",[mp.nstr(p/tot,3) for p in parts])
