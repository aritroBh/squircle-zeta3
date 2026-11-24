"""s07: (i) EXACT test that the collected lift defect Z (weight 4; D contains e sh Z) lies in P_4 = span of the
proved weight-4 base (A3h5's 5042 T0 rows); if not, also against P_4 + weight-2,3,4 frame lifts.
(ii) numerics (eval8, ~50 digits): Z.V, X.V, and u.f.V for every non-literal frame (row,u) pair, plus e.V."""
from common import *
sp = pickle.load(open(HERE + '/split.pkl', 'rb'))
B = bases(); W4 = words(4)
res = {}
xs, rk = exact_combination(B[4], [sp['Z']['rat'], sp['Z']['sq']], W4)
for c, x in zip(('rat', 'sq'), xs):
    res[c] = x
    log(f'Z_{c}: support {sum(1 for v in sp["Z"][c].values() if v)}; in P_4 (exact): {x is not None}', 's07.log')
if any(x is None for x in xs):
    gens = list(B[4])
    for n in (2, 3, 4):
        for row in pickle.load(open(T + f'/A3h4_frame_ds/validated_w{n}.pickle', 'rb'))['rows']:
            for w in words(4 - n): gens.append({k: F(v) for k, v in mul(dict(row), {w: 1}).items()})
    xs2, rk2 = exact_combination(gens, [sp['Z']['rat'], sp['Z']['sq']], W4)
    log(f'with frame lifts (rank {rk2}): in span: {[x is not None for x in xs2]}', 's07.log')
pickle.dump(res, open(HERE + '/Z_in_P4.pkl', 'wb'))
# numerics
sys.path.insert(0, T + '/A3f6_level16')
import mpmath as mp
from flint import ctx, acb, arb
import eval8 as ev
ctx.prec = 260; ctx.cap = 600; ev.N = 520
ev.roots = [acb(0, arb.pi() * j / 4).exp() for j in range(8)]
ev.roots[0], ev.roots[2], ev.roots[4], ev.roots[6] = acb(1), acb(0, 1), acb(-1), acb(0, -1)
mp.mp.dps = 60
cache = {}
def h(w):
    if w not in cache:
        z = ev.h(w); cache[w] = mp.mpc(z.real.mid().str(58, radius=False), z.imag.mid().str(58, radius=False))
    return cache[w]
def val(d): return mp.fsum(mp.mpf(F(c).numerator) / F(c).denominator * h(w) for w, c in d.items() if c)
def scale(d): return max(abs(F(c)) for c in d.values() if c)
log(f'e.V = {mp.nstr(val(E1), 15)}  (i*pi/2 = {mp.nstr(mp.pi/2, 15)}i)', 's07.log')
for c in ('rat', 'sq'):
    z = sp['Z'][c]; x = sp['X'][c]
    log(f'{c}: |Z.V|/max|Z| = {mp.nstr(abs(val(z)) / scale(z), 3)};  |X.V|/max|X| = {mp.nstr(abs(val(x)) / scale(x), 3)}', 's07.log')
VW5 = pickle.load(open(W5C1 + '/validated_w5.pkl', 'rb'))
worst = 0
for (i, u), st in sorted(sp['frame_pairs'].items()):
    d = actv({tuple(w): c for w, c in VW5['rows'][i]}, u)
    r = abs(val(d)) / sum(abs(c) for c in d.values())
    worst = max(worst, r)
log(f'max scaled |u.f . V| over the {len(sp["frame_pairs"])} frame (row,u) pairs: {mp.nstr(worst, 3)}', 's07.log')
