"""W4RX parts (c) and (d).
(c) exact verification of the certificate with OUR rows (check_b) and OUR target (check_a),
    reading D_certificate.json ourselves;
(d) numerical check, with our own 150-digit values, that every row used is a true relation
    (also every regenerated weight-4 standard row and every frame row stored in the notes)."""
import sys, json, pickle, time
sys.dont_write_bytecode = True
from fractions import Fraction as F
from w4rx_core import *
import mpmath as mp
T0 = time.time()
def log(s): print(f'[{time.time()-T0:6.1f}s] {s}', flush=True)
ROOT = '<project>/research/08_theory'
cert = json.load(open(ROOT + '/A3h5_exact_certificate/D_certificate.json'))
a = pickle.load(open('data/check_a.pickle', 'rb'))
b = pickle.load(open('data/check_b.pickle', 'rb'))
assert b['labels'] == cert['labels']
labels, rows = cert['labels'], b['rows']
assert all(l != ['P4'] for l in labels), 'a pi^4 row is present'
log(f"certificate: {len(labels)} labels, keys {sorted(cert)}, coefficient lists {[len(c) for c in cert['coefficients']]}; no pi^4 row")
# ---------------------------------------------------------------- (c) exact identity
hmax = 0; verdict = []
for j, name in enumerate(('rational part', 'sqrt2 part')):
    target = to_H(a['Ds'][j])
    acc = {}
    for y, row in zip(cert['coefficients'][j], rows):
        y = F(y); hmax = max(hmax, len(str(abs(y.numerator))), len(str(y.denominator)))
        if y == 0: continue
        for w, c in row: acc[w] = acc.get(w, 0) + y * c
    words = set(acc) | set(target)
    nonzero = [w for w in words if acc.get(w, 0) != target.get(w, 0)]
    allwords = [w for w in __import__('itertools').product(range(-1, 8), repeat=4) if w[0] != 0 and w[-1] != -1]
    assert words <= set(allwords)
    zero_mult = sum(1 for y in cert['coefficients'][j] if F(y) == 0)
    log(f'(c) {name}: target support {len(target)}, words touched {len(words)} of {len(allwords)}, '
        f'nonzero residual entries {len(nonzero)}; zero multipliers {zero_mult}')
    verdict.append(len(nonzero) == 0)
log(f'(c) max multiplier height (digits of numerator/denominator): {hmax}')
log(f'(c) VERDICT exact identity D = sum y_i rows_i for both parts, coefficient of pi^4 = 0: {all(verdict)}')
# ---------------------------------------------------------------- (d) numerics
vals = load_values('data/G_values_N900.json', 150)
def resid(hrow):
    g = from_H(dict(hrow))
    s = mp.mpc(0); scale = mp.mpf(0)
    for w, c in g.items():
        s += c * vals[w]; scale += abs(c) * abs(vals[w])
    return abs(s), scale
worst = (mp.mpf(0), None); worst_rel = mp.mpf(0); fam_worst = {}
for lab, row in zip(labels, rows):
    r, sc = resid(row)
    fam = lab[1][0] if lab[0] == 'T0' else f'frame_w{lab[1]}'
    fam_worst[fam] = max(fam_worst.get(fam, mp.mpf(0)), r)
    worst_rel = max(worst_rel, r / sc)
    if r > worst[0]: worst = (r, lab)
for f, r in fam_worst.items(): log(f'(d) family {f:10s} max |row value| = {mp.nstr(r, 4)}')
log(f'(d) all {len(rows)} certificate rows: max |value| = {mp.nstr(worst[0], 4)} (label {worst[1]}), max |value|/sum|c*G| = {mp.nstr(worst_rel, 4)}')
# all regenerated weight-4 standard rows
sys.path.insert(0, '.')
from w4rx_rel import standard_rows
std4, fam4 = standard_rows(4, {})
m = max(resid(r)[0] for r in std4)
log(f'(d) all {len(std4)} regenerated weight-4 standard rows: max |value| = {mp.nstr(m, 4)}')
json.dump({'c_exact_identity': all(verdict), 'c_height_digits': hmax, 'd_rows': len(rows),
           'd_max_abs': mp.nstr(worst[0], 5), 'd_max_rel': mp.nstr(worst_rel, 5),
           'd_family_max': {k: mp.nstr(v, 5) for k, v in fam_worst.items()},
           'd_std_all_max': mp.nstr(m, 5)}, open('data/check_cd.json', 'w'), indent=1)
log('saved data/check_cd.json')
# negative control: the numerical check separates true from false rows
ctrl = {}
for lab, row in zip(labels, rows):
    fam = lab[1][0] if lab[0] == 'T0' else f'frame_w{lab[1]}'
    if fam in ctrl: continue
    bad = ((row[0][0], row[0][1] + 1),) + tuple(row[1:])
    ctrl[fam] = resid(bad)[0]
log('(d) negative control (first row of each family with one coefficient changed by 1): ' +
    ', '.join(f'{k}: {mp.nstr(v, 3)}' for k, v in ctrl.items()))
