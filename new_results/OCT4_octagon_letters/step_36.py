import os, sys, json; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log
mp.mp.dps = 120
r2 = mp.sqrt(2); kap1=r2-1; kap2=r2/2-1
x1 = 1-r2; x2 = 1-r2/2   # the two nonlocal "points" = -kap
print("x1=",mp.nstr(x1,8)," x2=",mp.nstr(x2,8))
def Li21(x, y, N=None):
    x=mp.mpf(x); y=mp.mpf(y)
    if N is None:
        import math
        N = int(mp.mp.dps/max(mp.log10(1/abs(x)),0.05))+50
    s=mp.mpf(0); Hy=mp.mpf(0); py=mp.mpf(1)
    for n in range(1,N+1):
        # H_{n-1}(y): add y^{n-1}/(n-1) before using for index n (H_{n-1})
        s += x**n/n**2 * Hy
        py *= y
        Hy += py/n
    return s
# candidates with |x|,|y|<1
cand = {
 'Li21(x1,x2)': (x1,x2), 'Li21(x2,x1)': (x2,x1),
 'Li21(x1,x2/x1)': (x1, x2/x1), 'Li21(x2/x1,x1)': (x2/x1, x1),
 'Li21(x2,x2/x1)': (x2, x2/x1), 'Li21(x2/x1,x2)': (x2/x1, x2),
 'Li21(x1,x1)': (x1,x1), 'Li21(x2,x2)': (x2,x2),
}
cv={}
for k,(a,b) in cand.items():
    cv[k]=Li21(a,b); print(k,"=",mp.nstr(cv[k],20))
import pickle; pickle.dump({k:mp.nstr(v,80) for k,v in cv.items()}, open('li21_cand.pkl','wb'))
log("step: computed Li21 candidates")
