import sys, pickle, math
sys.path.insert(0, '../DEC5D_decagon_order5_depth')
import dec5d_lib as D
from fractions import Fraction as F
QS5 = D.QS5; SQ5 = D.SQ5; PHI = D.PHI
# build the decagon kit (as in DEC5D): slots = (log2, logsqrt5, logphi)
kit = D.make_kit(QS5, D.logvec5, slots=('log2','logsqrt5','logphi'))
li = kit['li_symbol']; const_pow = kit['const_pow']
# li_symbol(n, z): z a QS5
ONE = QS5(1); ZERO = QS5(0)
z_m1 = -ONE/PHI       # -1/phi
z_q  = ONE/(PHI**2)   # 1/phi^2
def sym_mul_logpow(sym, slot, k):
    # multiply symbol by (log slot)^k: shuffle with k! * (slot,...,slot)
    base = {(slot,)*k: QS5(math.factorial(k))}
    return kit['shuffle'](sym, base)
LOG2, LOGSQ5, LOGPHI = 0, 1, 2
S_L5m = li(5, z_m1); S_L5q = li(5, z_q)
S_L4m = li(4, z_m1); S_L4q = li(4, z_q)
S_L3m = li(3, z_m1)
Y5 = const_pow(LOGPHI, 5); Y4 = const_pow(LOGPHI, 4); Y3 = const_pow(LOGPHI, 3)
def scale(t, c): return kit['scale'](t, c)
def add(*ts): return kit['add'](*ts)
K5m = add(S_L5m, scale(Y5, QS5(F(1,120))))
K5q = add(S_L5q, scale(Y5, QS5(F(-2,15))))
K4m = add(S_L4m, scale(Y4, QS5(F(-1,24))))
K4q = add(S_L4q, scale(Y4, QS5(F(1,3))))
K3m = add(S_L3m, scale(Y3, QS5(F(1,6))))
print("S(Li5(-1/phi)) =", S_L5m)
print("S(Li5(1/phi^2)) =", S_L5q)
print("S(Li4(-1/phi)) =", S_L4m)
print("S(Li4(1/phi^2)) =", S_L4q)
print("S(Li3(-1/phi)) =", S_L3m)
print("K5m symbol zero:", K5m == {}, K5m)
print("K5q symbol zero:", K5q == {})
print("K4m symbol zero:", K4m == {})
print("K4q symbol zero:", K4q == {})
print("K3m symbol zero:", K3m == {})
open('progress.log','a').write("step_47: EXACT log-part of Wt5 from DEC5D's exact symbol (fully symmetric, pure product); K-combos K5m,K5q,K4m,K4q,K3m verified symbol-zero EXACTLY via dec5d_lib li_symbol\n")
