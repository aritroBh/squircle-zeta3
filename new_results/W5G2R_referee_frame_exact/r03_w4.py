"""r03: the weight-4 proved span S4 = span(G), G = T0[S] (5042) + weight-4 lifts (own shuffle) of the 13 A3h4 frame rows.
Mod two primes, own elimination: rank of G; whether 7.G (complex conjugates) lies in span(G); numerics of G."""
from rc import *
from relim import elim
LOG = 'r03.log'
W4 = WORDS(4); ID4 = {w: i for i, w in enumerate(W4)}
B = bases(); g4 = list(B[4])
A3 = {k: pickle.load(open(T + f'/A3h4_frame_ds/validated_w{k}.pickle', 'rb'))['rows'] for k in (2, 3, 4)}
FR13 = [(2, 3), (3, 3), (3, 5), (3, 6), (4, 4), (4, 6), (4, 8), (4, 9), (4, 11), (4, 12), (4, 40), (4, 41), (4, 43)]
for k, i in FR13:
    for w in WORDS(4 - k):
        d = shv({tuple(x): c for x, c in A3[k][i]}, {w: 1})
        if d: g4.append(d)
H = pickle.load(open(HERE + '/H_values.pkl', 'rb'))
worst = max(abs(val(g, H)[0]) / val(g, H)[1] for g in g4)
log(f'G: {len(g4)} rows; numerically max relative |val| {worst:.2e}', LOG)
for p in (4194301, 4194287):
    def modq(c): c = F(c); return (c.numerator % p) * pow(c.denominator % p, p - 2, p) % p
    def dense(rows):
        M = np.zeros((len(rows), len(W4)))
        for a, d in enumerate(rows):
            for w, c in d.items(): M[a, ID4[w]] = modq(c)
        return M
    A = dense(g4); Tm = dense([actv(g, 7) for g in g4])
    r, res = elim(A, Tm, p)
    log(f'p={p}: rank span(G) = {r} (t05 reports 5070 with all A3h4 lifts); 7.G inside span(G): {not np.any(res)} '
        f'({int(np.any(res, axis=1).sum())} of {len(g4)} conjugates outside)', LOG)
