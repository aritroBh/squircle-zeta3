# Planner (W4 session) verification of a ZQ5c D=13 weight-5 relation (out/d13/genuine_lattice_cap.json, row given on argv):
# (0) expand the lattice row to explicit terms with ZQ5c's w5verify.expand_a (pure algebra) and save it as
#     planner/d13_row{r}.json in the d17_final.json format (T on the SAME side: T_coef*L/sqrt13 + z5_coef*z5 + li5 + prods = 0);
# (1) value with the planner's own Li5 (ZQ5_zagier_weight5/planner/planner_li5.py), own Hurwitz L(chi13,5), own letter
#     logs from the exact generators, at 2200 digits (above ZQ5c's 1500/1800 LLL and 2100 fresh check) and at 1300 digits;
# (2) S-unit factorisations of x and 1-x by PSLQ over the 5 letters; (3) exact symbol of the right-hand side is ZERO.
import sys, os, json, itertools, time
from fractions import Fraction as F
from collections import defaultdict
import mpmath as mp
HERE = '<project>/research/08_theory/ZQ5c_zagier_w5_capped'
os.chdir(HERE); sys.path.insert(0, HERE)
sys.path.insert(1, '<project>/research/08_theory/ZQ5_zagier_weight5/planner')
from planner_li5 import li5 as pli5
import w5verify as V
r = int(sys.argv[1])
G = json.load(open('out/d13/genuine_lattice_cap.json'))
J = json.load(open('out/d13/candidates_cap.json'))
Kk = json.load(open('out/d13/kernel5_cap.json'))
names = J['names']
cands = [F(c['x'][0]) for c in J['cands']], [F(c['x'][1]) for c in J['cands']]
keepc = [0, 1, 2, 3] + list(range(23, 51))  # as in ZQ5c step_23 / step_24
assert len(G['atom_names']) == 2 + len(keepc) + 40
e = V.expand_a(13, J['cands'], Kk['K'], Kk['comps'], keepc, G['atom_names'], [int(x) for x in G['G'][r]], names)
T_coef = -F(e['a_T']) * 13 ** 5          # move T to the same side: T_coef*L/sqrt13 + z5_coef*z5 + sum li5 + sum prod = 0
z5_coef = F(e['a_z5'])
li5 = {f"{cands[0][i]}+{cands[1][i]}*sqrt13": F(c) for i, c in e['li5'].items()}
prods = {'*'.join(k): F(c) for k, c in e['prod'].items()}
GEN = {'lE13': (F(3, 2), F(1, 2)), 'ls13': (F(0), F(1)), 'l2': (F(2), F(0)), 'lp3': (F(5, 2), F(1, 2)), 'lp3c': (F(5, 2), F(-1, 2))}
assert set(names) == set(GEN), names
for n, (a, b) in GEN.items():
    print(f'letter {n} = {a}+{b}*sqrt13, norm {a*a - 13*b*b}')
out = {'D': 13, 'row': r, 'T_coef': str(T_coef), 'z5_coef': str(z5_coef), 'li5': {k: str(v) for k, v in li5.items()},
       'products': {k: str(v) for k, v in prods.items()}, 'letters': {n: f'{a}+{b}*sqrt13' for n, (a, b) in GEN.items()},
       'convention': 'T_coef*L(chi13,5)/sqrt13 + z5_coef*zeta(5) + sum li5 + sum products = 0 (same side as d17_final.json)'}
os.makedirs('planner', exist_ok=True)
json.dump(out, open(f'planner/d13_row{r}.json', 'w'), indent=1)
coefs = list(li5.values()) + list(prods.values()) + [z5_coef]
from math import lcm
print(f'row {r}: T_coef = {T_coef} = {T_coef/13**5} * 13^5; z5 {z5_coef}; {len(li5)} Li5, {len(prods)} products; '
      f'max |num| {max(abs(c.numerator) for c in coefs)}, lcm den {lcm(*[c.denominator for c in coefs])}', flush=True)


def chi13(a):
    a %= 13
    return 0 if a == 0 else (1 if pow(a, 6, 13) == 1 else -1)


def parse(s):
    s = s.replace('+-', '-'); i = s.rfind('*sqrt13'); head = s[:i]
    k = max(head.rfind('+'), head.rfind('-'))
    return F(head[:k]), F(head[k:])


args = [(parse(k), v) for k, v in li5.items()]
for dps in (1300, 2200):
    t0 = time.time(); mp.mp.dps = dps + 30
    s = mp.sqrt(13)
    val = lambda ab: mp.mpf(ab[0].numerator) / ab[0].denominator + mp.mpf(ab[1].numerator) / ab[1].denominator * s
    lv = {n: mp.log(abs(val(g))) for n, g in GEN.items()}; lv['pi2'] = mp.pi ** 2; lv['pi4'] = mp.pi ** 4
    L5 = mp.fsum(chi13(a) * mp.zeta(5, mp.mpf(a) / 13) for a in range(1, 13)) / mp.mpf(13) ** 5
    tot = T_coef.numerator * L5 / s / T_coef.denominator + z5_coef.numerator * mp.zeta(5) / z5_coef.denominator
    for ab, c in args:
        tot += c.numerator * pli5(val(ab), dps + 10) / c.denominator
    for nm, c in prods.items():
        tot += c.numerator * mp.fprod(lv[x] for x in nm.split('*')) / c.denominator
    print(f'(1) residual at {dps} digits: {mp.nstr(tot, 5)}   ({time.time()-t0:.0f}s; L(chi13,5) = {mp.nstr(L5, 25)})', flush=True)
mp.mp.dps = 120
s = mp.sqrt(13)
val = lambda ab: mp.mpf(ab[0].numerator) / ab[0].denominator + mp.mpf(ab[1].numerator) / ab[1].denominator * s
LET = list(GEN); lg = [mp.log(abs(val(GEN[n]))) for n in LET]


def fac(x):
    v = mp.pslq([mp.log(abs(x))] + lg, maxcoeff=10 ** 6, maxsteps=10 ** 6)
    assert v and v[0] != 0 and all(vi % v[0] == 0 for vi in v[1:]), f'no factorisation for {x}: {v}'
    ex = [-vi // v[0] for vi in v[1:]]
    assert abs(mp.log(abs(x)) - mp.fsum(a * l for a, l in zip(ex, lg))) < mp.mpf(10) ** -100
    return ex


sym = defaultdict(F)
for ab, c in args:
    x = val(ab); ex, e1 = fac(x), fac(1 - x)
    for t in itertools.product(range(5), repeat=5):
        w = e1[t[0]] * ex[t[1]] * ex[t[2]] * ex[t[3]] * ex[t[4]]
        if w: sym[t] -= c * w
idx = {n: k for k, n in enumerate(LET)}
for nm, c in prods.items():
    parts = nm.split('*')
    if 'pi2' in parts or 'pi4' in parts: continue
    for perm in itertools.permutations(range(5)):
        sym[tuple(idx[parts[i]] for i in perm)] += c
print(f'(2) all {len(args)} arguments: x and 1-x factor over the 5 letters (PSLQ, to 1e-100)')
print(f'(3) exact symbol: {sum(1 for v in sym.values() if v)} nonzero entries')
