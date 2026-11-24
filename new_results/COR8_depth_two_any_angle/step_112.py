
import mpmath as mp
mp.mp.dps = 120
import sys
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
from hpl_series import H_series, Li22
L2 = mp.log(2)
# high precision values
v = Li22(mp.mpf('0.5'), mp.mpf(1))
v2 = H_series((0,0,1,1), mp.mpf('0.5'))
v3 = H_series((0,1,0,1), mp.mpf('0.5'))
print("Li22(1/2,1):", mp.nstr(v, 40))
print("Li31(1/2,1):", mp.nstr(v2, 40))
basis = {'z4': mp.zeta(4), 'z3l2': mp.zeta(3)*L2, 'z2l2sq': mp.zeta(2)*L2**2, 'l2^4': L2**4,
         'Li4h': mp.polylog(4, mp.mpf('0.5'))}
blist = list(basis.items())
for val, nm in [(v, 'Li22(1/2,1)'), (v2, 'Li31(1/2,1)')]:
    rel = mp.pslq([val] + [b for _, b in blist], tol=mp.mpf('1e-110'), maxcoeff=10**14)
    print(nm, "PSLQ:", rel)
