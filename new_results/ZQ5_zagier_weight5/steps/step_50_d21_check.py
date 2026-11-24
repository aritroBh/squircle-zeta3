
import sys; sys.path.insert(0,'.')
import mpmath as mp
import w5lib as W
import pipe
D=21
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
vals = pipe.load_vals(D, 1150)
mp.mp.dps = 1200
worst = []
for i, cd in enumerate(cands):
    x = cd['x'].emb(1200)
    v2 = mp.polylog(5, x)
    dv = abs(vals[f"li5:{i}"] - v2)
    if dv > mp.mpf('1e-1100'):
        worst.append((i, pipe.fmt_qd(cd['x']), float(x), mp.nstr(dv,3)))
print("args with |series-polylog| > 1e-1100:")
for w in worst: print("  ", w)
# internal consistency of jonquiere at two precisions for the worst
for i, f, xf, _ in worst:
    cd = cands[i]
    v250 = W.li5_smart(cd['x'].emb(300), 250)
    v1150 = vals[f"li5:{i}"]
    print(f"jonquiere 250 vs cached 1150 for {f}: {mp.nstr(abs(v250-v1150),3)}")
# independent direct series at 600 digits for the worst arg (0.999917)
i0 = worst[0][0]
cd = cands[i0]
mp.mp.dps = 700
x = cd['x'].emb(700)
print("direct series at 600d for", pipe.fmt_qd(cd['x']), "x =", mp.nstr(x,20))
vs = W.li5_series(x, 600)
vj = W.li5_smart(cd['x'].emb(700), 600)
print("series - jonquiere:", mp.nstr(abs(vs-vj),3))
