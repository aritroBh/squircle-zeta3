from mpmath import mp, mpf
mp.dps = 30
from math import comb
def catalan(m): return comb(2*m, m)//(m+1)
def Kjv(j, x):
    return sum((-2*catalan(m-1))*(-1)**(j-m)*comb(m+j-2, j-m)*x**m for m in range(1, j+1))
x = mpf('0.0995')
for j in [1,2,3,4,5,6,8,10,15,20,30,50]:
    print(j, mp.nstr(Kjv(j, x), 6))
