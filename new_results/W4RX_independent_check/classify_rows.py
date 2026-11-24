"""W4RX: classify the certificate rows by the theorem each one needs
(elementary: both factors convergent / distribution of a convergent word;
 regularised: a divergent factor or a divergent distributed word)."""
import sys, json, pickle
sys.dont_write_bytecode = True
from collections import Counter
from w4rx_rel import *
ROOT = '<project>/research/08_theory'
cert = json.load(open(ROOT + '/A3h5_exact_certificate/D_certificate.json'))
lower = {}
fams = {}
for n in (1, 2, 3):
    rows, fam = standard_rows(n, dict(lower))
    fams[n] = dict(zip(rows, fam))
    cols, oldrows, basis = pickle.load(open(ROOT + f'/A3c_claude_elementary/Q/weight{n}.pickle', 'rb'))
    lower[n] = [tuple((tuple(w), int(c)) for w, c in r) for r in basis]
std4, fam4 = standard_rows(4, {})
def kind(f):
    if f[0] == 'ds':
        u, v = f[1], f[2]
        return f'ds({len(u)},{len(v)}) ' + ('regularised (u divergent)' if u[0] == 0 else 'convergent')
    if f[0] == 'dist':
        return f'dist d={f[1]} ' + ('regularised' if f[2][0] == 0 else 'convergent')
    if f[0] == 'branch': return 'weight-1 branch (log identity)'
    return f'lift of weight-{f[1]} basis row'
c = Counter()
for lab in cert['labels']:
    if lab[0] == 'T0' and lab[1][0] == 'std': c['std: ' + kind(fam4[lab[1][1]])] += 1
    elif lab[0] == 'T0':
        k = {'lift1': 1, 'lift2old': 2}[lab[1][0]]
        c[f'{lab[1][0]}: basis row = ' + kind(fams[k][lower[k][lab[1][1]]])] += 1
    else: c[f'frame w{lab[1]}'] += 1
for k, v in sorted(c.items()): print(f'{v:5d}  {k}')
print('weight-1 basis:', [kind(fams[1][r]) for r in lower[1]])
used2 = sorted(set(l[1][1] for l in cert['labels'] if l[0] == 'T0' and l[1][0] == 'lift2old'))
print('weight-2 basis rows used:', used2, [kind(fams[2][lower[2][i]]) for i in used2])
for i in used2: print('  ', lower[2][i])
uu = Counter(fam4[l[1][1]][1] for l in cert['labels'] if l[0] == 'T0' and l[1][0] == 'std'
             and fam4[l[1][1]][0] == 'ds' and fam4[l[1][1]][1][0] == 0)
print('divergent u in used double-shuffle rows (H notation):', dict(uu))
dd = Counter(fam4[l[1][1]][2][:2] for l in cert['labels'] if l[0] == 'T0' and l[1][0] == 'std'
             and fam4[l[1][1]][0] == 'dist' and fam4[l[1][1]][2][0] == 0)
print('leading letters of regularised distribution words used:', dict(dd))
