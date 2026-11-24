src = open("c2.py").read()
src = src.replace("""    acc = mpf(0)
    for k in range(1, K+1):
        acc += k*log(k)**2
        if k >= 2: P += power(k, -s-1)*acc""",
"""    acc = mpf(0)
    for k in range(1, K+1):
        if k >= 2: P += power(k, -s-1)*acc
        acc += k*log(k)**2""")
src = src.replace("""    acc = mpf(0)
    for k in range(1, K+1):
        acc += k*log(k)
        if k >= 2: P += power(k, -s-1)*log(k)*acc""",
"""    acc = mpf(0)
    for k in range(1, K+1):
        if k >= 2: P += power(k, -s-1)*log(k)*acc
        acc += k*log(k)""")
open("c2.py","w").write(src)
test = r'''
from c2 import *
from squircle import setdps
from mpmath import mp, nstr, zeta, diff as mpdiff
setdps(60)
def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)
print("G_l(4) =", nstr(G_l(4), 25), " (expect ~0.018649+tail)")
print("G_2(4) =", nstr(G_2(4), 25))
def c2_exact(s):
    return (2*s**2*zetapp(s-1) - s*(3*s+1)*zetapp(s)
            + 4*s*(s+1)*(G_2(s+1)+G_x(s+1)) - 8*s**2*G_l(s) - 4*s*G_2(s))
for s in [3,4,5,6]:
    print(f"s={s}: c_2={nstr(c2_exact(s), 30)}")
print("fits: 2.73718981077298834655 ; 1.573849628813923486848 ; 1.003456318418583651585 ; 0.646106636259136217485")
'''
exec(test)
