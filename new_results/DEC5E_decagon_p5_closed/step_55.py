import sys, pickle, math
sys.path.insert(0, '../DEC5D_decagon_order5_depth')
import dec5d_lib as D
from fractions import Fraction as F
QS5 = D.QS5; PHI = D.PHI
kit = D.make_kit(QS5, D.logvec5, slots=('log2','logsqrt5','logphi'))
gs = kit['g_symbol']
ZERO, ONE = QS5(0), QS5(1)
P_LET = -PHI; Q_LET = PHI**2
# Li32(x,y) <-> word (0,0,1/x,0,1/(xy)) (OCT5D3 convention)
cases = {
 'Li32(-1/phi,1)':     (-ONE/PHI, ONE),
 'Li32(1/phi^2,1)':    (ONE/PHI**2, ONE),
 'Li32(-1/phi,-1/phi)':(-ONE/PHI, -ONE/PHI),
 'Li32(1/phi^2,-phi)': (ONE/PHI**2, -PHI),
}
syms = {}
for nm, (x, y) in cases.items():
    a = ONE/x; b = ONE/(x*y)
    w = (ZERO, ZERO, a, ZERO, b)
    s = gs(ONE, w, ONE)
    syms[nm] = s
    print(nm, "word:", w, "symbol entries:", len(s))
# collect all words appearing
allw = set()
for s in syms.values(): allw.update(s.keys())
allw = sorted(allw)
print("distinct tensor words:", len(allw))
# matrix over Q: symbols have QS5 coeffs; split p,q
import sympy as sp
rows = []
names = list(syms)
for w in allw:
    rowp = []; rowq = []
    for nm in names:
        c = syms[nm].get(w, QS5(0))
        rowp.append(sp.Rational(c.p.numerator, c.p.denominator))
        rowq.append(sp.Rational(c.q.numerator, c.q.denominator))
    rows.append(rowp); rows.append(rowq)
M = sp.Matrix(rows)
print("symbol matrix rank over Q (p and q parts):", M.rank())
# kernel of the transpose: combos of the 4 Li32 with zero symbol
Mt = M.T
ker = Mt.nullspace()
print("symbol-kernel dimension among the 4 Li32:", len(ker))
for kv in ker:
    print("  kernel vector:", kv.T)
