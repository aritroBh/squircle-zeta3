src = open("squircle.py").read()
src = src.replace(
"""    if delta < 1:
        ld = log(delta)
        I1 = quad(lambda x: exp(a*x)*F(exp(x)), [ld, 0])
        I2 = quad(lambda t: t**(a-1)*F(t), [1, inf])
        large = I1 + I2
    else:
        large = quad(lambda t: t**(a-1)*F(t), [delta, inf])""",
"""    if delta < 1:
        ld = float(log(delta))
        if ld < -30:
            # chunked tanh-sinh on the log-interval (accuracy on long intervals)
            pts = [ld]
            while pts[-1] < -25:
                pts.append(pts[-1]+25)
            pts.append(0.0)
            I1 = mpf(0)
            for i in range(len(pts)-1):
                I1 += quad(lambda x: exp(a*x)*F(exp(x)), [pts[i], pts[i+1]])
        else:
            I1 = quad(lambda x: exp(a*x)*F(exp(x)), [ld, 0.0])
        I2 = quad(lambda t: t**(a-1)*F(t), [1, inf])
        large = I1 + I2
    else:
        large = quad(lambda t: t**(a-1)*F(t), [delta, inf])""")
open("squircle.py","w").write(src)
test = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs
import time
setdps(70)
for p_ in [25, 50]:
    t0=time.time(); ZA_ = ZA(p_, 4); t1=time.time()
    ZB_, iB = ZB(p_, 4, info=True); t2=time.time()
    print(f"p={p_}: agree={nstr(fabs(ZA_-ZB_)/fabs(ZA_),4)}  ZA {t1-t0:.0f}s ZB {t2-t1:.0f}s")
    print("   ZA:", nstr(ZA_, 45))
    print("   ZB:", nstr(ZB_, 45))
'''
exec(test)
with open("progress.log","a") as f: f.write("step_20: chunked log-quad in ZB; p=25,50 revalidated\n")
