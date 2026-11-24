from flint import acb, acb_poly, ctx
ctx.prec = 210
p = acb_poly([acb(1), acb(2), acb(3)])
q = p * p
print('mul:', q, q.length())
print('truncate:', q.truncate(4))
print('integral:', p.integral())
print('right_shift:', p.right_shift(2) if hasattr(p,'right_shift') else 'n/a')
# coefficient setting?
try:
    r = acb_poly([acb(0)]*5)
    print('coeffs attr:', r.coeffs)
except Exception as e: print(e)
import time
# speed test: multiplication at length 230
import random
a = acb_poly([acb(random.random(), random.random()) for _ in range(230)])
b = acb_poly([acb(random.random(), random.random()) for _ in range(230)])
t0=time.monotonic()
for _ in range(100):
    c = (a*b).truncate(231)
print('100 muls+truncate at N=230, prec 210: %.3fs'%(time.monotonic()-t0))
