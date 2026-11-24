# Planner verification of the D=17 relation found by planner_rerun_lll.py (corrected values):
# (1) the 830- and 1150-digit relations differ by an element of the exact zero-relation span;
# (2) the S-unit factorisations of x and 1-x (vx, v1m from candidates.json) hold numerically (own letter values);
# (3) the exact full symbol of the right-hand side (Li5 + products; pi^2, zeta(5), T have zero symbol) is ZERO;
# (4) a fresh value check at 1500 digits with own Li5 (planner_li5) and own L(chi17,5) by Hurwitz zeta.
import sys, json, pickle, itertools
from fractions import Fraction as F
from collections import defaultdict
import mpmath as mp
import flint
sys.path.insert(0, '<project>/research/08_theory/ZQ5_zagier_weight5/planner')
from planner_li5 import li5
ZQ = '<project>/research/08_theory/ZQ5_zagier_weight5'
D = 17
J = json.load(open(f'{ZQ}/out/d17/candidates.json'))
K = json.load(open(f'{ZQ}/out/d17/kernel.json'))['K']
r830 = pickle.load(open(f'{ZQ}/planner/d17_rerun_830.pkl', 'rb'))
r1150 = pickle.load(open(f'{ZQ}/planner/d17_rerun_1150.pkl', 'rb'))
names = r1150['names']; assert names == r830['names']
a1 = r1150['rels'][0][1]; a0 = r830['rels'][0][1]
if a0[0] != a1[0]: a0 = [-x for x in a0]
assert a0[0] == a1[0]
# (1) span test over Q
Z = [z[1] for z in r1150['zeros']] + [z[1] for z in r830['zeros']]
diff = [x - y for x, y in zip(a1, a0)]
rz = flint.fmpq_mat(Z).rank(); rzd = flint.fmpq_mat(Z + [diff]).rank()
print(f'(1) zero-relation rank {rz}; with the difference {rzd} -> difference in zero span: {rz == rzd}')
# expand relation to Li5 coefficients
nC = len(J['cands']); letters = J['names']; d = len(letters)
c_li5 = [F(0)] * nC; prods = {}; T_coef = a1[0]; z5_coef = a1[1]
for nm, a in zip(names, a1):
    if not a or nm in ('T', 'z5'): continue
    if nm.startswith('K'):
        for i, k in enumerate(K[int(nm[1:])]):
            if k: c_li5[i] += a * k
    else: prods[nm] = a
used = [i for i in range(nC) if c_li5[i]]
print(f'relation: T-coef {T_coef} = 17^5? {T_coef == 17**5 or -T_coef == 17**5}; z5 {z5_coef}; Li5 terms {len(used)}; product terms {len(prods)}')
# (2) factorisations, own letter values: 4+sqrt17, sqrt17, (3+sqrt17)/2, |(3-sqrt17)/2|
mp.mp.dps = 320
s = mp.sqrt(17)
Lv = [4 + s, s, (3 + s) / 2, abs((3 - s) / 2)]
print('letters per candidates.json:', J['disp'])
def xval(c):
    a, b = F(c['x'][0]), F(c['x'][1]); return mp.mpf(a.numerator) / a.denominator + mp.mpf(b.numerator) / b.denominator * s
worst = 0
for i in used:
    c = J['cands'][i]; x = xval(c)
    e1 = abs(mp.log(abs(x)) - sum(v * mp.log(l) for v, l in zip(c['vx'], Lv)))
    e2 = abs(mp.log(abs(1 - x)) - sum(v * mp.log(l) for v, l in zip(c['v1m'], Lv)))
    worst = max(worst, e1, e2)
print(f'(2) max |log|x| - vx.letters|, |log|1-x| - v1m.letters| over used args: {mp.nstr(worst, 3)}')
# (3) exact symbol
sym = defaultdict(F)
for i in used:
    c = J['cands'][i]; vx, v1 = c['vx'], c['v1m']
    for t in itertools.product(range(d), repeat=5):
        w = v1[t[0]] * vx[t[1]] * vx[t[2]] * vx[t[3]] * vx[t[4]]
        if w: sym[t] -= c_li5[i] * w
idx = {n: k for k, n in enumerate(letters)}
for nm, a in prods.items():
    toks = nm.split('*')
    if 'pi2' in toks: continue
    ms = [idx[t] for t in toks]; assert len(ms) == 5
    for perm in itertools.permutations(ms): sym[perm] += a
nz = {k: v for k, v in sym.items() if v}
print(f'(3) exact full symbol: {len(nz)} nonzero entries -> ZERO: {not nz}')
# (4) fresh value at 1500 digits
for P in (1500, 1000):
    mp.mp.dps = P + 30
    s = mp.sqrt(17); Lv = [4 + s, s, (3 + s) / 2, abs((3 - s) / 2)]; lg = [mp.log(l) for l in Lv]
    Lchi = mp.mpf(0)
    for n in range(1, 17):
        chi = 1 if pow(n, 8, 17) == 1 else -1
        Lchi += chi * mp.zeta(5, mp.mpf(n) / 17)
    Lchi /= mp.mpf(17)**5
    tot = T_coef * Lchi / s + z5_coef * mp.zeta(5)
    for i in used:
        a, b = F(J['cands'][i]['x'][0]), F(J['cands'][i]['x'][1])
        x = mp.mpf(a.numerator) / a.denominator + mp.mpf(b.numerator) / b.denominator * s
        tot += c_li5[i] * li5(x, P + 10)
    for nm, a in prods.items():
        v = mp.mpf(a)
        for t in nm.split('*'): v *= (mp.pi**2 if t == 'pi2' else lg[idx[t]])
        tot += v
    print(f'(4) dps {P}: L(chi17,5) = {mp.nstr(Lchi, 25)}  relation residual = {mp.nstr(tot, 4)}', flush=True)
lc = 1
for c in c_li5 + [F(x) for x in prods.values()]:
    if c: lc = lc * c.denominator // __import__('math').gcd(lc, c.denominator)
print('Li5 args used:', [J['cands'][i]['x'] for i in used])
json.dump({'D': 17, 'T_coef': T_coef, 'z5_coef': z5_coef,
           'li5': {f"{J['cands'][i]['x'][0]}+{J['cands'][i]['x'][1]}*sqrt17": str(c_li5[i]) for i in used},
           'products': {k: v for k, v in prods.items()}, 'letters': J['disp']},
          open(f'{ZQ}/planner/d17_relation_planner.json', 'w'), indent=1)
print('saved planner/d17_relation_planner.json')
