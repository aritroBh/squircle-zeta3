import sys, pickle, math
sys.path.insert(0, '../DEC5D_decagon_order5_depth')
import dec5d_lib as D
from fractions import Fraction as F
QS5 = D.QS5; PHI = D.PHI
kit = D.make_kit(QS5, D.logvec5, slots=('log2','logsqrt5','logphi'))
gs = kit['g_symbol']
ZERO, ONE = QS5(0), QS5(1)
cases = {
 'Li32(-1/phi,1)':     (-ONE/PHI, ONE),
 'Li32(1/phi^2,1)':    (ONE/PHI**2, ONE),
 'Li32(-1/phi,-1/phi)':(-ONE/PHI, -ONE/PHI),
 'Li32(1/phi^2,-phi)': (ONE/PHI**2, -PHI),
}
for nm, (x, y) in cases.items():
    a = ONE/x; b = ONE/(x*y)
    s = gs(ZERO, (ZERO, ZERO, a, ZERO, b), ONE)
    print(nm, ":", s)
# And the mixed ones' full symbol should be checked: print again with tensor structure
