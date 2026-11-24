
import pickle, mpmath as mp
from fractions import Fraction as F
import engine5, kimilib
import mzv_exact as M
from mzv_vals import dz
from mzv_num import mzv_num
mp.mp.dps = 70
forms10 = pickle.load(open('forms10.pkl','rb'))
fft = pickle.load(open('fft_check.pkl','rb'))
def parse_c(s):
    t = s.strip()
    if t.startswith('('): t=t[1:-1]
    t = t.replace('j','').strip()
    if ' + ' in t: re_,im_ = t.split(' + ')
    elif ' - ' in t[1:]:
        re_,im_ = t.split(' - '); im_='-'+im_
    else: return mp.mpf(t)
    return mp.mpc(mp.mpf(re_), mp.mpf(im_))
_mcache={}
def zval(w):
    if w in _mcache: return _mcache[w]
    if len(w)==1: v=mp.zeta(w[0])
    elif len(w)==2: v=dz(w[0],w[1])
    else: v=mzv_num(list(w), npts=10, N0=300, step=100)
    _mcache[w]=v
    return v
def fmp(q): return mp.mpf(q.numerator)/q.denominator
a=(1,2,4,3,5)
for W in range(2,9):
    T=engine5.table(W); d=T.nullity
    reg=[kimilib.peval5(forms10[W][b], a) for b in range(d)] if W>=3 else [F(0)]*d
    m0=kimilib.m0_coords(a, W)[W]
    full=[reg[b]+m0[b] for b in range(d)]
    tot = mp.mpf(0)
    for b in range(d):
        num = mp.mpf(0)
        for w,cc in M.wordcombo(T.monomials[b]).items():
            num += fmp(cc) * zval(w)
        tot += fmp(full[b])*num
    fv = parse_c(fft[W])
    rd = abs(fv-tot)/abs(tot)
    print(f"W={W}: reldiff = {mp.nstr(rd,3)}   exact={mp.nstr(tot,22)}")
