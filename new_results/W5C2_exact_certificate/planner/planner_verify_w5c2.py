# Planner (W4): independent check of the W5C2 certificate, written without W5C2's fold code.
# (A) EXACT unfolded identity over Q, per component c in {rat, sq}:
#       D_c = sum_sectors sum_rows (y_c / (4*lam)) * sum_g chi_s(g) * (g . r)
#     with r the unfolded source row (from W5C1's rels_*.pkl / validated_w5.pkl), g in Gal(Q(zeta8)/Q) acting on
#     letters x -> u*x mod 8 (u = 1, 7, 3, 5; letter -1 = dt/t fixed), chi_s(g) built from the sector signs.
#     Fractions throughout; compared on every word. Also reports which families are used.
# (B) The value step needs every g.r to be a valid relation. Numerically check (eval8, as W5C1's frame validation)
#     g.r . V = 0 for every used frame-w5 row and all four g, plus a random sample of other used transforms.
import sys, os, json, pickle, time, random
from fractions import Fraction as F
HERE = '<project>/research/08_theory/W5C2_exact_certificate'
W5C1 = HERE + '/../W5C1_certificate_scoping'
OUT = HERE + '/planner/planner_verify_w5c2.out'
def log(s):
    with open(OUT, 'a') as f: f.write(time.strftime('%H:%M:%S ') + s + '\n')
t0 = time.time()
UNITS = (1, 7, 3, 5)                      # g = 0..3 in W5C2's labelling: id, x->-x, x->3x, x->-3x
def act(w, u): return tuple(-1 if x == -1 else (u * x) % 8 for x in w)
def chi(gi, sec):                         # sec = (eps for x->-x, eps for x->3x)
    return (sec[0] if gi in (1, 3) else 1) * (sec[1] if gi in (2, 3) else 1)
SEC = {'pp': (1, 1), 'pm': (1, -1), 'mp': (-1, 1), 'mm': (-1, -1)}
sys.path.insert(0, W5C1 + '/../A3c_claude_elementary/Q')
sys.path.insert(0, '<project>/paper/squircle_note/scripts')
from relations import words
WL = words(5)
D = {c: {tuple(int(x) for x in k.split(',')): F(v) for k, v in json.load(open(f'{W5C1}/D_{c}.json')).items()}
     for c in ('rat', 'sq')}
SRC = {f: pickle.load(open(f'{W5C1}/rels_{f}.pkl', 'rb')) for f in ('std', 'lift1', 'lift2old', 'lift3old', 'lift4old')}
VW5 = pickle.load(open(f'{W5C1}/validated_w5.pkl', 'rb'))['rows']
log(f'loaded: {len(WL)} words, sources ' + str({f: len(v) for f, v in SRC.items()}) + f', frame-w5 {len(VW5)}')
acc = {'rat': {}, 'sq': {}}
used = {}            # (fam, idx) -> set of g needed
fams = {}
for tag, sec in SEC.items():
    cert = json.load(open(f'{HERE}/certificate_{tag}.json'))
    for ent in cert['rows']:
        s = ent['src']
        if isinstance(s, str): s = eval(s)
        fam, idx, lam = s['ufam'], int(s['uidx']), F(s['lam'])
        row = [(tuple(w), c) for w, c in VW5[idx]] if fam == 'framew5' else [(WL[w], c) for w, c in SRC[fam][idx][1]]
        fams[(tag, fam)] = fams.get((tag, fam), 0) + 1
        for c in ('rat', 'sq'):
            y = F(ent['y_' + c])
            if not y: continue
            base = y / (4 * lam)
            for gi, u in enumerate(UNITS):
                k = base * chi(gi, sec)
                used.setdefault((fam, idx), set()).add(gi)
                a = acc[c]
                for w, cw in row:
                    gw = act(w, u); a[gw] = a.get(gw, F(0)) + k * cw
    log(f'{tag}: accumulated ({time.time()-t0:.0f}s)')
log('families used per sector: ' + str(fams))
for c in ('rat', 'sq'):
    bad = sum(1 for w in set(acc[c]) | set(D[c]) if acc[c].get(w, F(0)) != D[c].get(w, F(0)))
    log(f'(A) EXACT unfolded identity {c}: words in D {len(D[c])}, words touched {len(acc[c])}, mismatches {bad}')
# (B) numerical validity of the Galois transforms
import mpmath as mp
from flint import ctx, acb, arb
sys.path.insert(0, W5C1 + '/../A3f6_level16')
import eval8 as ev
ctx.prec = 400; ctx.cap = 900; ev.N = 800
ev.roots = [acb(0, arb.pi() * j / 4).exp() for j in range(8)]
ev.roots[0], ev.roots[2], ev.roots[4], ev.roots[6] = acb(1), acb(0, 1), acb(-1), acb(0, -1)
mp.mp.dps = 100
random.seed(5)
fw = [k for k in used if k[0] == 'framew5']
others = random.sample([k for k in used if k[0] != 'framew5'], 40)
tests = []
for fam, idx in fw + others:
    row = [(tuple(w), c) for w, c in VW5[idx]] if fam == 'framew5' else [(WL[w], c) for w, c in SRC[fam][idx][1]]
    for gi, u in enumerate(UNITS): tests.append((fam, idx, gi, [(act(w, u), c) for w, c in row]))
need = set(w for t in tests for w, c in t[3])
log(f'(B) {len(fw)} frame-w5 rows x 4 g + 40 sampled other rows x 4 g = {len(tests)} transforms, {len(need)} words to evaluate')
vals = {}
for i, w in enumerate(need):
    z = ev.h(w); vals[w] = mp.mpc(z.real.mid().str(95, radius=False), z.imag.mid().str(95, radius=False))
    if i % 500 == 0: log(f'  evaluated {i}/{len(need)} ({time.time()-t0:.0f}s)')
worst = {}
for fam, idx, gi, row in tests:
    r = abs(mp.fsum(mp.mpf(c) * vals[w] for w, c in row)) / max(1, sum(abs(c) for w, c in row))
    key = 'framew5' if fam == 'framew5' else 'other'
    worst[(key, gi)] = max(worst.get((key, gi), mp.mpf(0)), r)
log('(B) max scaled residual |g.r . V| by (family, g): ' + ', '.join(f'{k}: {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
log(f'DONE ({time.time()-t0:.0f}s)')
