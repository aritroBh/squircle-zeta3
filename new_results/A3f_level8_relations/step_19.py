import sys
sys.path.insert(0,'Q')
from flint import ctx, acb, acb_series, acb_poly
ctx.prec = 1100
ctx.cap = 1150
print('sqrt2 radius:', acb(2).sqrt().real.rad())
s = acb_series([0]+[acb(1)/j for j in range(1,1150)])  # -log(1-x), length 1150
print('series len/prec:', len(s.coeffs()), s.prec if hasattr(s,'prec') else '?')
v = acb_poly(s.coeffs())(acb('0.5'))
print('log2 radius:', v.real.rad())
# now the culprit: division by a length-2 series
t = acb_series([1])
t2 = (t / acb_series([-acb(1),1])).integral()
print('t2 len:', len(t2.coeffs()), 'eval radius:', acb_poly(t2.coeffs())(acb('0.5')).real.rad())
