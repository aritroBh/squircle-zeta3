import sys, time, pickle
sys.path.insert(0,'.')
import point_engine as PE, mzv_exact as M
from fractions import Fraction as F

def make_f_eval(W, idx):
    T = PE.table(W)
    def f(a):
        a=[F(x) for x in a]
        acc=F(0)
        for pt in (a, a[::-1]):
            vw = PE.volreg_words_point(pt, W)
            for w,pp in vw[W].items():
                c=T.reduce_word(w)[idx]
                if c: acc+=pp*c
        return acc
    return f

f11 = make_f_eval(11,8)
# warmup table cache
_ = f11((1,2,3,4,5))
t0=time.time()
for p in [(2,1,3,1,2),(1,3,1,2,3),(3,1,2,3,1),(1,1,2,3,2),(2,3,1,1,3)]:
    _=f11(p)
print(f"f11 per point (warm cache): {(time.time()-t0)/5:.3f}s")

def Lval(a):
    a1,a2,a3,a4,a5=[F(x) for x in a]
    return a2*a3*a4*(a1+a2+a3+a4+a5)

# sanity: Q11(1,2,3,4,5) should be 109753920/360 = 304872
print("Q11(1,2,3,4,5) =", f11((1,2,3,4,5))/Lval((1,2,3,4,5)), "expect 304872")
