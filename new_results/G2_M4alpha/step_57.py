import sys, json
from fractions import Fraction as Fr
cwd = "<project>/research/08_theory/G2_M4alpha"
sys.path.insert(0, cwd)
import g2_fit as gf
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"

def load_g1_terms(name):
    out = []
    for nm, v in json.load(open(G1 + name)):
        if nm == '1': out.append(('1', (), 0, Fr(v))); continue
        if nm == 'c': out.append(('c', (), 1, Fr(v))); continue
        e = int(nm.split(' ', 1)[0][2:])
        w = tuple(int(t) for t in nm.split(' ', 1)[1][2:-1].split(','))
        out.append(('hpl', w, e, Fr(v)))
    return out

ORDER = 1040
Psi = load_g1_terms('hpl_fit_k3_(1,0,0).json')
# replicate step_08's combo_series
T = [Fr(0)]*(ORDER + 1)
per_term = {}
for kind, w, e, q in Psi:
    if kind == '1': T[0] += q
    elif kind == 'c': T[1] += q
    else:
        s = gf.hpl_series(w, ORDER)
        for m in range(ORDER + 1):
            if 0 <= m + e <= ORDER and s[m]:
                T[m + e] += q*s[m]
                if m + e == 1040:
                    per_term[(w, e)] = q*s[m]
# now compare with direct per-term computation
for kind, w, e, q in Psi:
    if kind == 'hpl':
        m_needed = 1040 - e
        if m_needed > ORDER:
            print("term", w, e, "needs s[", m_needed, "] > ORDER -- MISSING in combo_series!")
        else:
            s2 = gf.hpl_series(w, ORDER)
            if per_term.get((w, e), 0) != q*s2[m_needed]:
                print("diff in term", w, e)
tab = json.load(open(cwd + '/flat_exact_N520_K4_formula.json'))['4']
t520 = Fr(tab['520'].get('(0, 1, 0)', '0'))
print("combo_series T[1040] =", T[1040] == -t520, "  (-target =", -t520 == T[1040], ")")
