"""step_30_wparts.py - motivic weight parts Wt0..Wt5 of P^{(5)}(1)/L by G-word length.

Wt_k = 1200 * sum_{mon=(a,b,c,d,e,f), n} (p+q sqrt5) (log2)^e (logphi)^f
              * sum_{w in rewrite[(n,a,b,c,d)], |w|+e+f = k} gamma_w G(w;1).
Also splits each Wt_k into blocks (local / single-1 / single-4 / mixed).
Usage: python step_30_wparts.py DPS N
Writes wparts_dec_d{DPS}.pkl: {'Wt{k}': nstr, 'blk{k}': {block: nstr}}.
"""
import sys, os, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
DEC5D = "<project>/research/08_theory/DEC5D_decagon_order5_depth"
sys.path.insert(0, DEC5D)
import dec5d_lib as D
import gengine_dec2p
import mpmath as mp

def main():
    dps, N = int(sys.argv[1]), int(sys.argv[2])
    t0 = time.time()
    QS5 = D.QS5; ZERO, ONE = QS5(0), QS5(1)
    P_LET = -D.PHI; Q_LET = D.PHI**2
    G1, letters = gengine_dec2p.gtab(dps, N, cachedir=os.path.dirname(os.path.abspath(__file__)))
    mp.mp.dps = dps
    L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
    r5 = mp.sqrt(5); L2v = mp.log(2); LFv = mp.log((1+r5)/2)
    def q5f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r5
    def w2mpf(w):
        return tuple(L0 if x==ZERO else L1 if x==ONE else LP if x==P_LET else LQ for x in w)
    rewrite = pickle.load(open(os.path.join(DEC5D, "rewrite_dec_exact.pkl"),"rb"))
    parsed = pickle.load(open(os.path.join(DEC5D, "parsed_mons_dec.pkl"),"rb"))
    Gv = {}
    def gval(w):
        if w not in Gv: Gv[w] = G1[w2mpf(w)]
        return Gv[w]
    Wt = {k: mp.mpf(0) for k in range(6)}
    blk = {k: {'local': mp.mpf(0),'single-1': mp.mpf(0),'single-4': mp.mpf(0),'mixed': mp.mpf(0)} for k in range(6)}
    for mon, cu in parsed:
        a,b,c,d,e,f = mon
        cls = 'mixed' if (c>0 and d>0) else (('single-1' if c>0 else 'single-4') if (c>0 or d>0) else 'local')
        cf = L2v**e * LFv**f
        for n,(p,q) in cu.items():
            pf = (mp.mpf(str(p)) + mp.mpf(str(q))*r5)*cf
            for w, gam in rewrite[(n,a,b,c,d)].items():
                k = len(w) + e + f
                if k <= 5:
                    v = pf*q5f(gam)*gval(w)
                    Wt[k] += v; blk[k][cls] += v
    Wt = {k: 1200*v for k, v in Wt.items()}
    blk = {k: {c2: 1200*v for c2, v in row.items()} for k, row in blk.items()}
    print("assembled dps=%d in %.1fs" % (dps, time.time()-t0))
    for k in range(6):
        print("Wt%d =" % k, mp.nstr(Wt[k], 50))
    print("sum  =", mp.nstr(sum(Wt.values()), 50))
    gold = pickle.load(open(os.path.join(DEC5D, "p5_gold.pkl"),"rb"))
    P5L = mp.mpf(gold['P5overL'])
    print("P5/L gold (113d) - sum:", mp.nstr(abs(P5L - sum(Wt.values())), 3))
    out = {f'Wt{k}': mp.nstr(Wt[k], dps-30) for k in range(6)}
    for k in range(6):
        out[f'blk{k}'] = {c2: mp.nstr(v, dps-30) for c2, v in blk[k].items()}
    pickle.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"wparts_dec_d{dps}.pkl"),"wb"))
    print("saved wparts_dec_d%d.pkl" % dps)

if __name__ == "__main__":
    main()
