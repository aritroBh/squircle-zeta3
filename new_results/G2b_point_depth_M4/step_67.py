
import time
import mpmath as mp
mp.mp.dps = 60

# Fixed-level tanh-sinh via mpmath internals, single interval [a,b]
from mpmath.calculus.quadrature import TanhSinh
TS = TanhSinh(mp.ctx)

def ts_quad(f, a, b, level=7):
    if a == b: return mp.mpf(0)
    sgn = 1
    if b < a: a, b = b, a; sgn = -1
    # tanh-sinh on [a,b]: x = a + (b-a)*phi(t), sum over nodes at given level
    S = mp.mpf(0)
    h = mp.mpf(2) ** (-level)
    # Use the transformation of TanhSinh class directly
    return sgn * TS.transform_and_sum(f, a, b, level) if hasattr(TS, 'transform_and_sum') else None

# check available API
print([m for m in dir(TS) if not m.startswith('_')])
