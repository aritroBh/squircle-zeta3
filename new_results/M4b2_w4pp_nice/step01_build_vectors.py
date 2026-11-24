"""step_01: build W4'' level-8 word vectors (rational and sqrt2 parts) from gpl_w4pp.tsv.
Conversion mirrors cmzv.vectors()/verify_D.rebuild_target exactly:
  term = forms_path(forms, 1->A) - forms_path(forms, 1->B), times atom(log2)^k;
  rat part gets coefficient.p * term, sq part gets coefficient.q * term.
Also re-derives the G-expression from omega4.json (gexpr logic) and asserts equality with the TSV.
Saves stage1_w4pp_vectors.pkl."""
import sys, json, re, pickle, time
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/research/08_theory/M4b2_w4pp_nice")
import m4b2_common as mc
from fractions import Fraction as F
from qsqrt2 import QS2 as Q, SQ2 as R, A, B
import cmzv as cm

t0 = time.time()
M4A = mc.THEORY/'M4a_decomposition'

# ---------- parse the TSV ----------
def parse_qs2(s):
    s = s.strip()
    m = re.fullmatch(r"\((.+?)([+-])(\d+(?:/\d+)?)\*sqrt2\)", s)
    if m:
        p = F(m.group(1)); q = F(m.group(3))
        if m.group(2) == '-': q = -q
        return Q(p, q)
    m = re.fullmatch(r"(\d+(?:/\d+)?)\*sqrt2", s)
    if m:
        return Q(0, F(m.group(1)))
    return Q(F(s))

def parse_word(s):
    s = s.strip()
    assert s.startswith('(') and s.endswith(')')
    body = s[1:-1].strip()
    if not body: return ()
    # split top-level commas (elements may be parenthesised)
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
with (M4A/'gpl_w4pp.tsv').open() as f:
    for line in f:
        if line.startswith('#') or line.startswith('coefficient_p'): continue
        p, q, k, w = line.rstrip('\n').split('\t')
        rows.append((parse_word(w), int(k), Q(F(p), F(q))))
mc.log(f'parsed {len(rows)} rows from gpl_w4pp.tsv')

# ---------- cross-check: re-derive from omega4.json via gexpr logic ----------
from w4_symbol import log_product, LT, MT, clean
om = json.load((M4A/'omega4.json').open())
om = {tuple(map(int, kk.split(','))): v for kk, v in om.items()}
POLE = {'0': Q(0), 'a': A, '-a': -A, 'b': B, '-b': -B}
from collections import defaultdict
out = defaultdict(lambda: Q(0))
for (i, j) in ((2, 1), (3, 0)):
    factors = (LT,) * i + (MT,) * j
    for name, (p, q) in om[(i, j)]:
        z = POLE[name]
        coefficient = Q(F(p), F(q))
        for (w, k), d in log_product(*factors).items():
            out[((z,) + w, k)] += coefficient * d
rederived = dict(clean(out))
tsv = defaultdict(lambda: Q(0))
for w, k, c in rows: tsv[(w, k)] += c
assert dict(tsv) == rederived, 'TSV does not match omega4.json re-derivation'
mc.log(f'TSV == omega4.json re-derivation: OK ({len(rederived)} terms); '
       f'weight check: {all(len(w)+k == 4 for (w, k) in tsv)}')

# ---------- convert to level-8 word vectors ----------
lhs = [{}, {}]
for (w, k), coefficient in rederived.items():
    forms = [cm.dlog(p) for p in w]
    term = cm.add((cm.forms_path(forms, cm.POINTS[Q(1)], cm.POINTS[A]), 1),
                  (cm.forms_path(forms, cm.POINTS[Q(1)], cm.POINTS[B]), -1))
    for _ in range(k): term = cm.mul(term, cm.atom(('log', Q(2))))
    for j, c in enumerate((coefficient.p, coefficient.q)):
        lhs[j] = cm.add((lhs[j], 1), (term, c))
mc.log(f'W4pp supports (rational, sqrt2): {tuple(map(len, lhs))} s={time.time()-t0:.1f}')

cols = set(mc.words4())
assert all(w in cols for w in lhs[0]) and all(w in cols for w in lhs[1])
mc.log('all words admissible: OK')
with open(str(mc.HERE/'stage1_w4pp_vectors.pkl'), 'wb') as f:
    pickle.dump({'rat': lhs[0], 'sq': lhs[1]}, f)
mc.log('saved stage1_w4pp_vectors.pkl')
