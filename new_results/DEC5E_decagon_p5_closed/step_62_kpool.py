"""step_62_kpool.py - K-parts of golden weight-5 double polylogs from the G-table + exact symbols.

For each depth-2 word w = (0^{a-1},A,0^{b-1},B), (a,b) in {(3,2),(4,1)}, (A,B) in {1,-phi,phi^2}^2\{(1,1)}:
  K(w) = G(w;1) - (c_w/120) log(phi)^5,  c_w exact (S(G(w;1)) = c_w Y^5).
Usage: python step_62_kpool.py TABD
"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "<project>/research/08_theory/DEC5D_decagon_order5_depth")
import mpmath as mp
import dec5d_lib as DL
import gengine_dec2p

def main():
    TAB = int(sys.argv[1])
    N = {1550: 5250, 1150: 3900, 800: 2700, 520: 1750}[TAB]
    QS5 = DL.QS5; PHI = DL.PHI; ZERO, ONE = QS5(0), QS5(1)
    kit = DL.make_kit(QS5, DL.logvec5, slots=('log2','logsqrt5','logphi'))
    gs = kit['g_symbol']
    P_LET, Q_LET = -PHI, PHI**2
    G1, letters = gengine_dec2p.gtab(TAB, N, cachedir=os.path.dirname(os.path.abspath(__file__)))
    mp.mp.dps = TAB
    L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
    def w2mpf(w):
        return tuple(L0 if x==ZERO else L1 if x==ONE else LP if x==P_LET else LQ for x in w)
    lf = mp.log((1+mp.sqrt(5))/2)
    letter = {'1': ONE, 'p': P_LET, 'q': Q_LET}
    out = {}
    for (a_idx, b_idx) in [(3,2),(4,1)]:
        for At in '1pq':
            for Bt in '1pq':
                if At == '1' and Bt == '1': continue
                A, B = letter[At], letter[Bt]
                w = (ZERO,)*(a_idx-1) + (A,) + (ZERO,)*(b_idx-1) + (B,)
                s = gs(ZERO, w, ONE)
                assert len(s) == 1 and list(s.keys())[0] == (2,)*5, (a_idx,b_idx,At,Bt,s)
                c = s[(2,)*5]
                mp.mp.dps = TAB
                cf = mp.mpf(str(c.p)) + mp.mpf(str(c.q))*mp.sqrt(5)
                K = G1[w2mpf(w)] - cf/120*lf**5
                out[f"L{a_idx}{b_idx}{At}{Bt}"] = {'K': mp.nstr(K, TAB-20), 'sym': (str(c.p), str(c.q))}
    for nm, d in out.items():
        print(nm, "sym = (%s,%s)" % d['sym'], " K =", mp.nstr(mp.mpf(d['K']), 12))
    pickle.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"kpool_{TAB}.pkl"),"wb"))
    print("saved kpool_%d.pkl (%d K-parts)" % (TAB, len(out)))

if __name__ == "__main__":
    main()
