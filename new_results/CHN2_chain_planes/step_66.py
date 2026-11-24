
import pickle, mpmath as mp
from fractions import Fraction as F
import engine5, kimilib
import mzv_exact as M
from mzv_vals import dz
from mzv_num import mzv_num
mp.mp.dps = 70
forms10 = pickle.load(open('forms10.pkl','rb'))
fft = pickle.load(open('fft_check.pkl','rb'))

_mcache={}
def zval(w):
    if w in _mcache: return _mcache[w]
    if len(w)==1: v=mp.zeta(w[0])
    elif len(w)==2: v=dz(w[0],w[1])
    else: v=mzv_num(list(w), npts=10, N0=300, step=100)
    _mcache[w]=v
    return v

a=(1,2,4,3,5)
print(f"{'W':>3} {'reldiff':>10}  exact value")
for W in range(2,9):
    T=engine5.table(W); d=T.nullity
    if W>=3:
        reg=[kimilib.peval5(forms10[W][b], a) for b in range(d)]
    else:
        reg=[F(0)]*d
    m0=kimilib.m0_coords(a, W)[W]
    full=[reg[b]+m0[b] for b in range(d)]
    tot = mp.mpf(0)
    for b in range(d):
        wc = M.wordcombo(T.monomials[b])
        num = mp.mpf(0)
        for w,cc in wc.items():
            num += mp.mpf(cc.numerator)/mp.denominator(cc) * zval(w)
        tot += (mp.mpf(full[b].numerator)/full[b].denominator)*num
    fv = mp.mpc(fft[W])
    rd = abs(fv-tot)/abs(tot)
    print(f"{W:>3} {mp.nstr(rd,3)}   {mp.nstr(tot,22)}")
