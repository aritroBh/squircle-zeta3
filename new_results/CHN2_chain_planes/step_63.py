
import pickle, mpmath as mp
from fractions import Fraction as F
import engine5, kimilib
mp.mp.dps = 70
forms10 = pickle.load(open('forms10.pkl','rb'))
fft = pickle.load(open('fft_check.pkl','rb'))
a=(1,2,4,3,5)
print(f"{'W':>3} {'FFT':>30} {'exact':>30} {'reldiff':>10}")
for W in range(2,9):
    T=engine5.table(W); d=T.nullity
    reg=[kimilib.peval5(forms10[W][b], a) for b in range(d)]
    m0=kimilib.m0_coords(a, W)[W]
    full=[reg[b]+m0[b] for b in range(d)]
    T=engine5.table(W)
    # numeric value of the exact coordinate vector
    zv = mp.mpf(0)
    for b in range(d):
        mon = T.monomials[b]
        val = mp.mpf(1)
        for word in mon:
            val *= mp.zeta(word[0]) if len(word)==1 else mp.mpf(0)  # placeholder
        zv += mp.mpf(0)
    # better: evaluate via mzv of words: use mp.pslq? simplest: sum over wordcombo
    import mzv_exact as M
    tot = mp.mpf(0)
    for b in range(d):
        wc = M.wordcombo(T.monomials[b])
        num = mp.mpf(0)
        for w,cc in wc.items():
            num += mp.mpf(cc.numerator)/cc.denominator * mp.zeta(list(w))
        tot += (mp.mpf(full[b].numerator)/full[b].denominator)*num
    fv = mp.mpc(fft[W])
    rd = abs(fv-tot)/abs(tot)
    print(f"{W:>3} {mp.nstr(fv,25)} {mp.nstr(tot,25)} {mp.nstr(rd,3)}")
