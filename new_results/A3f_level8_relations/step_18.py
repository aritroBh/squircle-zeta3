import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, 'Q')
from flint import ctx, acb
import check_conversion as CC
print('prec after import:', ctx.prec, 'cap:', ctx.cap)
ctx.prec = 1100
ctx.cap = 1150
print('prec after set:', ctx.prec, 'cap:', ctx.cap)
from check_conversion import h
t0=time.monotonic()
v = h((7,7))
print('time', time.monotonic()-t0)
print('radius:', v.real.rad(), v.imag.rad())
