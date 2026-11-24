# Planner (W4 session): exact symbol check of ZQ5's D=12 relation (out/d12/rel_min_atom.json), complementing the
# value check planner_check_d12.py. Letters lE3 = 2+sqrt3, ls3 = sqrt3, lu = 1+sqrt3 (log of absolute values; symbols mod
# torsion). Li5(x) -> -(1-x) (x) x (x) x (x) x x; products of letters -> sum over orderings; pi^2, zeta(5), T -> 0.
# Factorisations of x and 1-x are found here by PSLQ and checked to 1e-100.
import json, re, itertools
from fractions import Fraction as F
from collections import defaultdict
import mpmath as mp
mp.mp.dps = 130
r3 = mp.sqrt(3)
LET = ['lE3', 'ls3', 'lu']; lg = [mp.log(2 + r3), mp.log(r3), mp.log(1 + r3)]
def fld(s):
    s = s.strip('()').replace('sqrt3', 's'); v = mp.mpf(0)
    for t in re.findall(r'[+-]?[^+-]+', s):
        if 's' in t:
            c = t.replace('*s', '').replace('s', ''); c = c if c not in ('', '+', '-') else c + '1'
            v += mp.mpf(F(c).numerator) / F(c).denominator * r3
        else: v += mp.mpf(F(t).numerator) / F(t).denominator
    return v
def split(nm):
    out, dep, cur = [], 0, ''
    for ch in nm:
        dep += (ch == '(') - (ch == ')')
        if ch == '*' and dep == 0: out.append(cur); cur = ''
        else: cur += ch
    return out + [cur]
def fac(x):
    v = mp.pslq([mp.log(abs(x))] + lg, maxcoeff=10 ** 6, maxsteps=10 ** 6)
    assert v and v[0] != 0 and all(vi % v[0] == 0 for vi in v[1:]), f'no factorisation for {x}: {v}'
    e = [-vi // v[0] for vi in v[1:]]
    assert abs(mp.log(abs(x)) - mp.fsum(a * l for a, l in zip(e, lg))) < mp.mpf(10) ** -100
    return e
d = json.load(open('../out/d12/rel_min_atom.json'))
sym = defaultdict(F); nli = nprod = 0; other = set()
for nm, c in zip(d['names'], d['c']):
    c = F(c)
    if c == 0 or nm in ('T', 'z5'): continue
    parts = split(nm)
    if len(parts) == 1 and parts[0].startswith('Li5'):
        x = fld(re.fullmatch(r'Li5\((.*)\)', parts[0]).group(1)); ex, e1 = fac(x), fac(1 - x); nli += 1
        for t in itertools.product(range(3), repeat=5):
            w = e1[t[0]] * ex[t[1]] * ex[t[2]] * ex[t[3]] * ex[t[4]]
            if w: sym[t] -= c * w
    else:
        nprod += 1
        if any(p in ('pi2', 'z5', 'T') or p.startswith('Li5') for p in parts):
            other.update(p for p in parts if p not in LET); continue
        assert len(parts) == 5 and all(p in LET for p in parts), nm
        for perm in itertools.permutations(range(5)):
            sym[tuple(LET.index(parts[i]) for i in perm)] += c
print(f'D=12: {nli} Li5 terms (all x, 1-x factor over 2+sqrt3, sqrt3, 1+sqrt3), {nprod} product terms '
      f'(zero-symbol factors seen: {sorted(other)}); exact symbol: {sum(1 for v in sym.values() if v)} nonzero entries')
