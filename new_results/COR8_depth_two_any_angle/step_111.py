
import mpmath as mp
mp.mp.dps = 90
import sys
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
from hpl_series import Li22
L2 = mp.log(2)
basis = {'z4': mp.zeta(4), 'z3l2': mp.zeta(3)*L2, 'z2l2sq': mp.zeta(2)*L2**2, 'l2^4': L2**4,
         'Li4h': mp.polylog(4, mp.mpf('0.5')), 'z3': mp.zeta(3), 'z2l2': mp.zeta(2)*L2, 'l2^3': L2**3,
         'Li3h': mp.polylog(3, mp.mpf('0.5')), 'z2': mp.zeta(2), 'l2^2': L2**2,
         'Li2h': mp.polylog(2, mp.mpf('0.5')), 'l2': L2, '1': mp.mpf(1)}
blist = list(basis.items())
v = Li22(mp.mpf('0.5'), mp.mpf(1))
print("Li22(1/2,1) =", mp.nstr(v, 30))
x = [v] + [b for _, b in blist]
rel = mp.pslq(x, tol=mp.mpf('1e-80'), maxcoeff=10**12)
print("PSLQ:", rel)
# also try Li_{3,1}(1/2,1) and H_{0,0,1,1}(1/2) style
from hpl_series import H_series
v2 = H_series((0,0,1,1), mp.mpf('0.5'))
print("H_{0,0,1,1}(1/2) =", mp.nstr(v2, 30))
rel2 = mp.pslq([v2] + [b for _, b in blist], tol=mp.mpf('1e-80'), maxcoeff=10**12)
print("PSLQ:", rel2)
