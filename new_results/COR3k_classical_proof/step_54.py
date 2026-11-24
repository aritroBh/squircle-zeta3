import json, re, sympy as sp
from cor3_common import *

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
c = sp.symbols('c', positive=True)

R2 = json.load(open(G1 + 'hpl_fit_k2_(0,0,0).json'))
CDraw = json.load(open(G1 + 'classical_dict.json'))

# ---- build the assembled coefficient table: key (fun_name, const_name) -> element of Q[c,1/c]
def parse_term(term):
    fn, cn, fr = term
    return fn, cn, sp.Rational(fr)

coef = {}
def add(fname, cconst, val):
    key = (fname, cconst)
    coef[key] = sp.simplify(coef.get(key, 0) + val)

for name, cstr in R2:
    mm = re.match(r'c\^(-?\d+) H\[(.*)\]', name)
    e = int(mm.group(1)); w = [int(z) for z in mm.group(2).split(',')]
    pref = sp.Rational(cstr)*c**e
    if w == [0, 1]:
        add('Li2(x)', '1', pref)
    elif w == [0, -1]:
        add('Li2(-x)', '1', -pref)
    else:
        for fn, cn, fr in CDraw[str(w)]:
            add(fn, cn, pref*sp.Rational(fr))

print("assembled keys:", len(coef))
for k in sorted(coef):
    if coef[k] != 0:
        print("   ", k, ":", sp.factor(coef[k]))
