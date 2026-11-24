"""M4i step B: verify L_relation.json.
(a) EXACT full point symbol of RHS == 0 (L(chi8,5) has zero symbol).
(b) every Li5 atom reduces to 0 in the exact (3,2) cobracket.
(c) FRESH value at 1620 digits: every atom evaluated anew (arb, 1680 digits);
    L(chi8,5) computed independently by Hurwitz zeta,
    L = 8^-5 [zeta(5,1/8)+zeta(5,7/8)-zeta(5,3/8)-zeta(5,5/8)],
    and cross-checked against Lchi8_5.json. Expect residual ~1e-1600.
Precision set AFTER all imports.
"""
import json, math, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
from flint import arb, ctx
import mpmath as mp
import w5_symbol as W
import w5_cobracket as C
from w5_cobracket import reduce as cred, delta_word
from w5_search import evaluate
import nice_search as NS

ZERO, ONE, Q, R, A, B = W.ZERO, W.ONE, W.Q, W.R, W.A, W.B
def log(s): print(s, flush=True)

EXTRA_Z = {'Li5((-2+2*sqrt2))': -2+2*R, 'Li5((17-12*sqrt2))': B**4,
           'Li5((1/2-1/4*sqrt2))': Q(F(1,2))-R/4, 'Li5((1/2+1/4*sqrt2))': Q(F(1,2))+R/4,
           'Li5((1/2-1/2*sqrt2))': Q(F(1,2))-R/2, 'Li5((1/2-3/8*sqrt2))': Q(F(1,2))-3*R/8,
           'Li5((-16+12*sqrt2))': -16+12*R}

terms = json.loads((HERE/'L_relation.json').read_text())['terms']
assert 'Lchi8_5' not in terms
sol = NS.solve(list(NS.S0), 'vlrel')
recmap = {nm: rec for nm, _, rec in sol['basis']}
symmap = {nm: s for nm, s, _ in sol['basis']}
for nm in terms:
    if nm in EXTRA_Z:
        z = EXTRA_Z[nm]; symmap[nm] = W.li_symbol(5, z); recmap[nm] = ('li', 5, z)
    elif nm == 'z5':
        symmap[nm] = {}; recmap[nm] = None
co = {nm: (F(p), F(q)) for nm, (p, q) in terms.items()}

# (a) EXACT SYMBOL of RHS == 0
t = {}
for nm, (a, b) in co.items():
    if a or b:
        t = W.add(t, W.scale(symmap[nm], Q(a, b)))
t = W.clean(t)
assert t == {} or t == {k: v for k, v in t.items() if v} == {} or not t, f'SYMBOL NONZERO: {t}'
log('(a) EXACT full point symbol of RHS == 0: PASS')

# (b) every Li5 reduces to 0 in the exact (3,2) cobracket
C.ban_float()
nli5 = 0
for nm, (a, b) in co.items():
    if nm.startswith('Li5') and (a or b):
        z = recmap[nm][2]
        assert not cred(delta_word(ZERO, (ZERO, ZERO, ZERO, ZERO, ONE/z), ONE)).d, nm
        nli5 += 1
log(f'(b) all {nli5} Li5 atoms reduce to 0 in the exact (3,2) cobracket: PASS')

# (c) FRESH VALUE at 1620 digits
digits = 1620
ctx.prec = math.ceil((digits+140)*math.log2(10)); mp.mp.dps = digits+120
DD = digits+60
def qnum(p, q):
    return mp.mpf(p.numerator)/p.denominator + mp.sqrt(2)*mp.mpf(q.numerator)/q.denominator
total = mp.mpf(0)
for nm, (a, b) in co.items():
    if not (a or b): continue
    if nm == 'z5':
        v = mp.mpf(arb(5).zeta().mid().str(DD+100, radius=False))
    else:
        v = mp.mpf(evaluate(recmap[nm], DD).mid().str(DD+100, radius=False))
    total += qnum(a, b)*v
# independent L(chi8,5) via Hurwitz zeta
s = mp.mpf(5)
L_hur = mp.mpf(8)**(-5)*(mp.zeta(s, mp.mpf(1)/8) + mp.zeta(s, mp.mpf(7)/8)
                         - mp.zeta(s, mp.mpf(3)/8) - mp.zeta(s, mp.mpf(5)/8))
L_file = mp.mpf(json.loads((HERE/'Lchi8_5.json').read_text())['value'])
log(f'(c0) Hurwitz L(chi8,5) - Lchi8_5.json = {mp.nstr(L_hur - L_file, 8)}  (1620-digit fresh)')
res = total - L_hur
log(f'(c) L_relation: RHS - L(chi8,5)[Hurwitz] = {mp.nstr(res, 8)}  (fresh {digits}-digit evaluation)')
res2 = total - L_file
log(f'(c’) RHS - L(chi8,5)[Lchi8_5.json] = {mp.nstr(res2, 8)}')
