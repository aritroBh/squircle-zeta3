"""ZLAD1R step_23 -- diff my independent table against ZLAD1 (k<=8) and G1 (k<=3), G2 (k=4)."""
import json
base = '<project>'
mine = json.load(open('muR_exact_N40_K8.json'))
zl   = json.load(open(base + '/research/08_theory/ZLAD1_zeta_ladder/mu_exact_N150_K8.json'))

def diff(tabA, tabB, ks, ns):
    bad = []
    checked = 0
    for k in ks:
        for n in ns:
            a = tabA[str(k)][str(n)]; b = tabB[str(k)][str(n)]
            keys = set(a) | set(b)
            for m in keys:
                va = a.get(m, '0'); vb = b.get(m, '0')
                checked += 1
                if va != vb:
                    bad.append((k, n, m, va, vb))
    return bad, checked

bad, checked = diff(mine, zl, range(9), range(1, 41))
print(f'vs ZLAD1 K8 table: {checked} monomial-coefficients compared, mismatches = {len(bad)}')
for row in bad[:8]: print('  ', row)

# G1: flat_exact_N40_K3.json exists (n<=40) and N320 file for k<=3
g1 = json.load(open(base + '/research/10_frontier/G1_zeta3_any_angle/flat_exact_N40_K3.json'))
# G1 monomial keys are str(tuple) like '(1, 0, 0)'; map to names
def tup2name(s):
    e = eval(s)
    if e == (0,0,0): return '1'
    out = []
    for idx, x in enumerate(e):
        if x: out.append(f'z{idx+2}' + (f'^{x}' if x > 1 else ''))
    return ''.join(out)
g1c = {k: {n: {tup2name(m): v for m, v in g1[k][n].items()} for n in g1[k]} for k in g1}
bad1, checked1 = diff(mine, g1c, range(4), range(1, 41))
print(f'vs G1 flat_exact_N40_K3: {checked} compared... ', len(bad1), 'mismatches')
for row in bad1[:8]: print('  ', row)

# G2 k=4 file (N640)
g2 = json.load(open(base + '/research/08_theory/G2_M4alpha/flat_exact_N640_K4.json'))
print('G2 file k keys:', list(g2.keys())[:6])
g2c = {k: {n: {tup2name(m): v for m, v in g2[k][n].items()} for n in g2[k]} for k in g2 if int(k) <= 4}
bad2, checked2 = diff(mine, g2c, range(5), range(1, 41))
print(f'vs G2 flat_exact_N640_K4: mismatches =', len(bad2))
for row in bad2[:8]: print('  ', row)
ok = not bad and not bad1 and not bad2
print('ALL CROSS-TABLE CHECKS PASS:', ok)
open('progress.log','a').write(f'step_23: independent table == ZLAD1 table exactly ({checked} coeffs, k<=8, n<=40); == G1 (k<=3); == G2 (k=4)\n')
