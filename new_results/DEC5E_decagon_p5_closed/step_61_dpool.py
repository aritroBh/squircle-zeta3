"""step_61_dpool.py - double-polylog golden pool at high precision + exact K-parts.

Li_{3,2}(x,y) = sum_{n>m} x^n y^m/(n^3 m^2), Li_{4,1}(x,y) = sum_{n>m} x^n y^m/(n^4 m),
for the 8 (x,y) cases with (1/x, 1/(xy)) in {1,-phi,phi^2}^2 \ {(1,1)}.
Exact symbol coefficient c (S(value) = c Y^5) via dec5d_lib; K = value - (c/120) log(phi)^5.
Usage: python step_61_dpool.py DPS
"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "<project>/research/08_theory/DEC5D_decagon_order5_depth")
import mpmath as mp
import dec5d_lib as DL

def lid(a, b, x, y, D):
    mp.mp.dps = D
    if x == 1:
        # Li_{a,b}(1,y) = sum_m y^m/m^b (zeta(a) - H_m^{(a)})
        za = mp.zeta(a)
        s = mp.mpf(0); ym = mp.mpf(1); Ha = mp.mpf(0); m = 0
        thresh = mp.mpf(10)**(-D-20)
        while True:
            m += 1; ym *= y; Ha += mp.mpf(1)/m**a
            s += ym*(za - Ha)/m**b
            if m > 20 and abs(ym) < thresh: break
        return s
    s = mp.mpf(0); xn = mp.mpf(1); ym = mp.mpf(1); Hy = mp.mpf(0); m = 0
    thresh = mp.mpf(10)**(-D-20)
    n = 0
    while True:
        # advance inner sum to m = n+1? we need Hy = sum_{m<n} y^m/m^b
        n += 1
        xn *= x
        if n >= 2:
            s += xn*Hy/mp.mpf(n)**a
        m = n
        ym *= y; Hy += ym/mp.mpf(m)**b
        if n > 20 and abs(xn*Hy) < thresh: break
    return s

def main():
    D = int(sys.argv[1])
    mp.mp.dps = D
    phi = (1+mp.sqrt(5))/2
    QS5 = DL.QS5; PHI = DL.PHI; ZERO, ONE = QS5(0), QS5(1)
    kit = DL.make_kit(QS5, DL.logvec5, slots=('log2','logsqrt5','logphi'))
    gs = kit['g_symbol']
    def symcoef(a_idx, b_idx, A, B):
        # word for Li_{a,b}: 0^{a-1}, A, 0^{b-1}, B
        w = (ZERO,)*(a_idx-1) + (A,) + (ZERO,)*(b_idx-1) + (B,)
        s = gs(ZERO, w, ONE)
        # expect single entry (2,)*5
        assert len(s) == 1 and list(s.keys())[0] == (2,)*5, (a_idx, b_idx, A, B, s)
        return s[(2,)*5]
    P_LET = -PHI; Q_LET = PHI**2
    letter = {1: ONE, 2: P_LET, 3: Q_LET}  # tags: 1 -> letter 1, 2 -> -phi, 3 -> phi^2
    cases = []
    for At in (1,2,3):
        for Bt in (1,2,3):
            if At == 1 and Bt == 1: continue
            A = letter[At]; B = letter[Bt]
            x = ONE/A; xy = ONE/B; y = xy/x
            cases.append((At, Bt, A, B, float(x), float(y)))
    out = {}
    for (At, Bt, A, B, xf, yf) in cases:
        for (a_idx, b_idx) in [(3,2),(4,1)]:
            c = symcoef(a_idx, b_idx, A, B)
            nm = f"Li{a_idx}{b_idx}_x{At}{Bt}"
            mp.mp.dps = D
            v = lid(a_idx, b_idx, mp.mpf(xf), mp.mpf(yf), D)
            out[nm] = {'val': mp.nstr(v, D-20), 'sym': (str(c.p), str(c.q)), 'args': (xf, yf)}
            print(nm, "sym =", c, " val =", mp.nstr(v, 12))
    pickle.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"dpool_{D}.pkl"),"wb"))
    print("saved dpool_%d.pkl (%d values)" % (D, len(out)))

if __name__ == "__main__":
    main()
