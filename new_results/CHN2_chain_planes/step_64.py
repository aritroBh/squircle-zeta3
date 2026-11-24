
import pickle, mpmath as mp
from fractions import Fraction as F
import engine5, kimilib
import mzv_exact as M
mp.mp.dps = 70
forms10 = pickle.load(open('forms10.pkl','rb'))
fft = pickle.load(open('fft_check.pkl','rb'))
a=(1,2,4,3,5)
print(f"{'W':>3} {'reldiff':>10}  detail")
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
            num += mp.mpf(cc.numerator)/cc.denominator * mp.zeta(list(w))
        tot += (mp.mpf(full[b].numerator)/full[b].denominator)*num
    fv = mp.mpc(fft[W])
    rd = abs(fv-tot)/abs(tot)
    print(f"{W:>3} {mp.nstr(rd,3)}   FFT={mp.nstr(fv,20)}  exact={mp.nstr(tot,20)}")
