"""step_01: recompute in-house 110-digit targets for the three weight-5 blocks
(local / single-ratio / mixed two-ratio) of P^(5)(1)/u from the read-only
parsed_mons.pkl + itab110.pkl. Prints residuals against the quoted values."""
import pickle
import mpmath as mp
mp.mp.dps = 120
SRC = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
parsed = pickle.load(open(SRC+"parsed_mons.pkl","rb"))
ITAB = {k: mp.mpf(v) for k, v in pickle.load(open(SRC+"itab110.pkl","rb")).items()}
r2 = mp.sqrt(2); L2C = mp.log(2); LU = mp.log(1+r2)
contrib = {'local': mp.mpf(0), 'single': mp.mpf(0), 'mixed': mp.mpf(0)}
tot5 = mp.mpf(0)
for mon, cu in parsed:
    if sum(mon) != 5: continue
    aa,bb,cx,dd,ee,ff = mon
    cls = 'mixed' if (cx>0 and dd>0) else ('single' if (cx>0 or dd>0) else 'local')
    constfac = L2C**ee * LU**ff
    for n,(p,q) in cu.items():
        term = (mp.mpf(str(p)) + mp.mpf(str(q))*r2)*constfac*ITAB[(n,aa,bb,cx,dd)]
        contrib[cls] += term; tot5 += term
blocks = {k: 960*v for k,v in contrib.items()}; tot5 *= 960
for k,v in blocks.items(): print(k, "=", mp.nstr(v, 60))
print("total W5/u =", mp.nstr(tot5, 60))
quoted = mp.mpf("4.2299562389856442293498943256846704703130475386317")
print("mixed - quoted =", mp.nstr(blocks['mixed']-quoted, 4))
pickle.dump({k: str(v) for k,v in blocks.items()} | {'total': str(tot5)}, open("targets110.pkl","wb"))
