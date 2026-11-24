t = r'''
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
import sympy as sp, json, re
data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
W4 = [(parse_term(s)[0], parse_term(s)[1], sp.Rational(c)) for s, c in data if len(parse_term(s)[1]) == 4]
def sym_Omega_scaled(c0):
    S = {}
    for e, w, cf in W4:
        S = sym_add(S, sym_scale(sym_H(w, c0), cf * (c0 ** e)))
    return S
cor6 = {e['label']: e for e in json.load(open('<project>/research/08_theory/COR6_true_point_symbols/cor6_table.json'))}
# prime order mapping: cor6 index 0,1,2,3 = primes 2,3,5,7 (their S0 from norms)
prime_order = [sp.Integer(2), sp.Integer(3), sp.Integer(5), sp.Integer(7)]
def cor6_D(label):
    out = {}
    for k, v in cor6[label]['true_d22'].items():
        out[k] = v[0]
    return out
def my_D(c0):
    S = sym_expanded(sym_Omega_scaled(c0))
    D = delta22(S)
    # re-key with cor6 indices
    out = {}
    for ((a, b), (c, d)), v in D.items():
        key = f"(({prime_order.index(a)}, {prime_order.index(b)}), ({prime_order.index(c)}, {prime_order.index(d)}))"
        out[key] = v
    return out
for c0, label in ((sp.Rational(2,5), 'c=2/5'), (sp.Rational(3,7), 'c=3/7')):
    mine = my_D(c0); theirs = cor6_D(label)
    keys = set(mine) | set(theirs)
    diffs = [(k, mine.get(k, 0), theirs.get(k, 0)) for k in sorted(keys) if sp.Rational(mine.get(k, 0)) != sp.Rational(theirs[k] if k in theirs else 0)]
    print(f"=== {label}: my components {len(mine)}, cor6 components {len(theirs)}, mismatches {len(diffs)}")
    for k, m, t in diffs:
        print("    ", k, "mine:", m, " cor6:", t)
    if not diffs:
        print("    VALUES MATCH EXACTLY")
# 1/3 and 1/2 zero check
for c0, label in ((sp.Rational(1,3), 'c=1/3'), (sp.Rational(1,2), None)):
    D = delta22(sym_expanded(sym_Omega_scaled(c0)))
    print(f"c0={c0}: scaled delta22 = {'ZERO' if not D else D}")
'''
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=900)
print(r.stdout); print("ERR:", (r.stderr or "")[-800:])
