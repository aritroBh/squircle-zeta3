import sys, json
from fractions import Fraction as Fr
cwd = "<project>/research/08_theory/G2_M4alpha"
sys.path.insert(0, cwd)
import g2_fit as gf

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
ORDER = 1040

def load_g1_terms(name):
    out = []
    for nm, v in json.load(open(G1 + name)):
        if nm == '1': out.append(('1', (), 0, Fr(v))); continue
        if nm == 'c': out.append(('c', (), 1, Fr(v))); continue
        e = int(nm.split(' ', 1)[0][2:])
        w = tuple(int(t) for t in nm.split(' ', 1)[1][2:-1].split(','))
        out.append(('hpl', w, e, Fr(v)))
    return out

Psi = load_g1_terms('hpl_fit_k3_(1,0,0).json')
# debug: which term contributes to m=1040
tab = json.load(open(cwd + '/flat_exact_N520_K4_formula.json'))['4']
t520 = Fr(tab['520'].get('(0, 1, 0)', '0'))
tot = Fr(0)
detail = []
for kind, w, e, q in Psi:
    if kind == 'hpl':
        m_needed = 1040 - e
        s = gf.hpl_series(w, m_needed + 1)
        detail.append((w, e, q, s[m_needed]))
        tot += q*s[m_needed]
print("sum of -Psi terms at m=1040:", -tot)
print("target t520:", t520)
print("diff:", t520 + tot)
