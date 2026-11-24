"""r02: read W5C2's four certificates directly (not W5G1's used_rows.pkl). Store y mod two primes for every used row,
and exact y for lift1 rows on weight-1 base rows 0,1 (the branch rows) and for frame-w5 rows."""
from rcommon import *
import json, pickle, collections, gc
sys.set_int_max_str_digits(0)
L = 'r02.log'
P = (2305843009213693951, 1000000007)   # 2^61-1 and 1e9+7 (not the primes used by W5C2 / W5G1)
W5C2 = T + '/W5C2_exact_certificate'; W5C1 = T + '/W5C1_certificate_scoping'
lift1 = pickle.load(open(W5C1 + '/rels_lift1.pkl', 'rb'))
def modq(fr, p): return fr.numerator % p * pow(fr.denominator % p, -1, p) % p
out = {'mod': [], 'exact': []}
lams = collections.Counter()
for tag in ('pp', 'pm', 'mp', 'mm'):
    cert = json.load(open(f'{W5C2}/certificate_{tag}.json'))
    cnt = collections.Counter()
    for ent in cert['rows']:
        s = ent['src']
        if isinstance(s, str): s = eval(s)
        fam, idx, lam = s['ufam'], int(s['uidx']), F(s['lam'])
        lams[lam] += 1
        yr, ys = F(ent['y_rat']), F(ent['y_sq'])
        cnt[fam] += 1
        out['mod'].append((tag, fam, idx, lam, tuple(modq(yr, p) for p in P), tuple(modq(ys, p) for p in P)))
        if fam == 'framew5' or (fam == 'lift1' and lift1[idx][0][1] in (0, 1)):
            out['exact'].append((tag, fam, idx, lam, yr, ys))
    log(f'{tag}: rows {len(cert["rows"])} {dict(cnt)}', L)
    del cert; gc.collect()
log(f'lambda values: {dict(lams)}', L)
pickle.dump(out, open(HERE + '/r02_extract.pkl', 'wb'))
ex = out['exact']
u1 = {(f, i) for t, f, i, *_ in ex if f == 'lift1'}
fr = {(f, i) for t, f, i, *_ in ex if f == 'framew5'}
log(f'unique lift1 rows on branch base rows 0/1 used (any sector): {len(u1)}; unique frame-w5 rows: {len(fr)} {sorted(i for _, i in fr)}', L)
