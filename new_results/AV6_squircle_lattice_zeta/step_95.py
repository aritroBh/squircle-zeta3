src = open("c1.py").read()
src = src.replace("glaiser", "glaisher")
open("c1.py","w").write(src)
test = r'''
from c1 import *
from squircle import setdps
from mpmath import mp, nstr, zeta, diff as mpdiff
setdps(70)
def zetap(a): return mpdiff(lambda x: zeta(x), a)
for s in [3,4,5,6]:
    g = G_sigma(s, K=1200)
    c1 = 4*s*(zetap(s)-zetap(s-1)) - 8*s*g
    print(f"s={s}: G(s)={nstr(g, 25)}")
    print(f"      c_1 = {nstr(c1, 40)}")
'''
exec(test)
