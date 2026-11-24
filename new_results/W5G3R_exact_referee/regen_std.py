"""Regenerate the standard weight-5 rows (regularised double shuffle + distribution d=2,4,8) with the
earlier referee W4RX's independent generator (W4RX_independent_check/w4rx_rel.py, imports only w4rx_core;
no module of the notes), and compare with W5C1 rels_std.pkl (the rows the certificates use), exactly.
Also weights 1-4 (for the lower bases and for T0[S]); results saved to regen_std.pkl (small)."""
from common import *
sys.path.insert(0, T + '/W4RX_independent_check')
import w4rx_rel as RX
LOG = 'regen_std.log'
t0 = time.time()
out = {}
for n in (1, 2, 3, 4):
    rows, fam = RX.standard_rows(n, {})
    out[n] = (rows, fam)
    log(f'weight {n}: {len(rows)} DS+distribution rows generated ({time.time()-t0:.0f}s)', LOG)
rows5, fam5 = RX.standard_rows(5, {}, log=lambda s: None)
log(f'weight 5: {len(rows5)} rows generated, families {dict(__import__("collections").Counter(f[0] for f in fam5))} ({time.time()-t0:.0f}s)', LOG)
W5 = WORDS(5)
stored = load_data(T + '/W5C1_certificate_scoping/rels_std.pkl')
st = [tuple((W5[k], int(c)) for k, c in row) for lab, row in stored]
assert [lab for lab, _ in stored] == [('std', k) for k in range(len(stored))]
same_list = st == list(rows5)
mine = set(rows5)
missing = [k for k, r in enumerate(st) if r not in mine]
neg = sum(1 for k, r in enumerate(st) if r not in mine and tuple((w, -c) for w, c in r) in mine)
log(f'weight 5 std: stored {len(st)} rows; identical ordered list: {same_list}; stored rows not among regenerated: {len(missing)} '
    f'(of which equal up to sign: {neg}); first missing: {missing[:5]}', LOG)
pickle.dump({'w1to4': {n: out[n][0] for n in out}, 'fam': {n: out[n][1] for n in out}, 'same5': same_list, 'missing5': missing},
            open(HERE + '/regen_std.pkl', 'wb'))
log('saved regen_std.pkl', LOG)
