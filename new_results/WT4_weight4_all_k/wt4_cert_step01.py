"""WT4 cert step01: build target (D1 or D2) level-8 word vectors from d1.tsv/d2.tsv.
Mirrors M4b2 step01: term = forms_path(forms, 1->A) - forms_path(forms, 1->B), times atom(log2)^k;
rat part gets coefficient.p, sq part coefficient.q. Cross-check: re-derive G-expression from
deltas_deg3.json (exact IBP kernel differences) via gexpr logic and assert equality with the TSV.
Usage: wt4_cert_step01.py <D1|D2>"""
import sys, json, re, pickle, time, os
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wt4_cert_common as mc
from fractions import Fraction as F
from qsqrt2 import QS2 as Q, SQ2 as R, A, B
import cmzv as cm

t0 = time.time()
TGT = sys.argv[1]
assert TGT in ('D1', 'D2')
tsv_path = os.path.join(HERE, TGT.lower() + '.tsv')

def parse_qs2(s):
    s = s.strip()
    m = re.fullmatch(r"\((.+?)([+-])(\d+(?:/\d+)?)\*sqrt2\)", s)
    if m:
        p = F(m.group(1)); q = F(m.group(3))
        if m.group(2) == '-': q = -q
        return Q(p, q)
    m = re.fullmatch(r"(-?\d+(?:/\d+)?)\*sqrt2", s)
    if m: return Q(0, F(m.group(1)))
    return Q(F(s))

def parse_word(s):
    s = s.strip()
    assert s.startswith('(') and s.endswith(')')
    body = s[1:-1].strip()
    if not body: return ()
    parts, depth, cur = [], 0, ''
    for ch in body:
        if ch == '(': depth += 1; cur += ch
        elif ch == ')': depth -= 1; cur += ch
        elif ch == ',' and depth == 0: parts.append(cur); cur = ''
        else: cur += ch
    parts.append(cur)
    parts = [x for x in parts if x.strip()]
    return tuple(parse_qs2(x) for x in parts)

rows = []
with open(tsv_path) as f:
    for line in f:
        if line.startswith('#') or line.startswith('coefficient_p'): continue
        p, q, k, w = line.rstrip('\n').split('\t')
        rows.append((parse_word(w), int(k), Q(F(p), F(q))))
mc.log(f'[{TGT}] parsed {len(rows)} rows from {os.path.basename(tsv_path)}')

from w4_symbol import log_product, LT, MT, clean
deltas = json.load(open(os.path.join(HERE, 'deltas_deg3.json')))
POLE = {'0': Q(0), 'a': A, '-a': -A, 'b': B, '-b': -B}
from collections import defaultdict
out = defaultdict(lambda: Q(0))
for (i, j) in ((2, 1), (3, 0)):
    factors = (LT,) * i + (MT,) * j
    for name, (p, q) in deltas[f'{i},{j}'][TGT]:
        z = POLE[name]; coefficient = Q(F(p), F(q))
        for (w, k), d in log_product(*factors).items():
            out[((z,) + w, k)] += coefficient * d
rederived = dict(clean(out))
tsv = defaultdict(lambda: Q(0))
for w, k, c in rows: tsv[(w, k)] += c
assert dict(tsv) == rederived, f'{TGT} TSV does not match deltas re-derivation'
mc.log(f'[{TGT}] TSV == deltas_deg3.json re-derivation: OK ({len(rederived)} terms); '
       f'weight check: {all(len(w)+k == 4 for (w, k) in tsv)}')

lhs = [{}, {}]
for (w, k), coefficient in rederived.items():
    forms = [cm.dlog(p) for p in w]
    term = cm.add((cm.forms_path(forms, cm.POINTS[Q(1)], cm.POINTS[A]), 1),
                  (cm.forms_path(forms, cm.POINTS[Q(1)], cm.POINTS[B]), -1))
    for _ in range(k): term = cm.mul(term, cm.atom(('log', Q(2))))
    for j, c in enumerate((coefficient.p, coefficient.q)):
        lhs[j] = cm.add((lhs[j], 1), (term, c))
mc.log(f'[{TGT}] supports (rational, sqrt2): {tuple(map(len, lhs))} s={time.time()-t0:.1f}')

cols = set(mc.words4())
assert all(w in cols for w in lhs[0]) and all(w in cols for w in lhs[1])
mc.log(f'[{TGT}] all words admissible: OK')
with open(os.path.join(HERE, f'stage1_{TGT}_vectors.pkl'), 'wb') as f:
    pickle.dump({'rat': lhs[0], 'sq': lhs[1]}, f)
mc.log(f'[{TGT}] saved stage1_{TGT}_vectors.pkl')
