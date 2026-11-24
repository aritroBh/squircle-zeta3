
import time
from mpmath import mp, mpf, fabs
mp.dps = 125
exec(open('av2core.py').read())
for p in [8, 64, 1024]:
    FPr = FPrecomp(p, '0.0012', 3800)
    ref = FPr.F(mpf('0.3'))
    print("p=%d ref(h2=0.0012): %s" % (p, mp.nstr(ref, 50)))
    for (h2, n2) in [('0.004', 1150), ('0.003', 1530), ('0.002', 2300), ('0.0015', 3050)]:
        FP = FPrecomp(p, h2, n2)
        v = FP.F(mpf('0.3'))
        print("   h2=%s: err=%s" % (h2, mp.nstr(fabs(v-ref),4)))
