
import mpmath as mp
mp.mp.dps = 60
from mpmath.calculus.quadrature import TanhSinh
ctx = mp.mp
TS = TanhSinh(ctx)
print([m for m in dir(TS) if not m.startswith('_')])
import inspect
print(inspect.signature(TS.sum_next))
