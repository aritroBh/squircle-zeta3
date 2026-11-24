import sys
sys.path.insert(0,'Q')
from flint import ctx
ctx.prec = 1100
ctx.cap = 1150
import check_conversion as CC   # module sets ctx.prec=200, cap=320, N=320 at import!
print('after import:', ctx.prec, ctx.cap)
ctx.prec = 1100
ctx.cap = 1150
print('reset:', ctx.prec, ctx.cap)
# but CC.N is baked to 320 and cached funcs use range(1,N)
s = CC.series((7,),0)
print('series len:', len(s.coeffs()))
v = CC.half((7,7),0)
print('half radius:', v.real.rad())
import mpmath as mp
print('N baked:', CC.N)
