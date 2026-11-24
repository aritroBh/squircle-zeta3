import json, re, sympy as sp
from cor3_common import *

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
c = sp.symbols('c', positive=True)
q = (1 - c**2)/(2*c)

R2 = json.load(open(G1 + 'hpl_fit_k2_(0,0,0).json'))
CDraw = json.load(open(G1 + 'classical_dict.json'))
coef = {}
def add(tab, fname, cconst, val):
    tab[(fname, cconst)] = sp.simplify(tab.get((fname, cconst), 0) + val)
for name, cstr in R2:
    mm = re.match(r'c\^(-?\d+) H\[(.*)\]', name)
    e = int(mm.group(1)); w = [int(z) for z in mm.group(2).split(',')]
    pref = sp.Rational(cstr)*c**e
    if w == [0, 1]: add(coef, 'Li2(x)', '1', pref)
    elif w == [0, -1]: add(coef, 'Li2(-x)', '1', -pref)
    else:
        for fn, cn, fr in CDraw[str(w)]:
            add(coef, fn, cn, pref*sp.Rational(fr))

# ---- typed B2' from result.md, minus the zeta(2)*Phi(c) term ----------------
typed = {}
T = [
    ('Li3(x)', '1', 1 + c), ('Li3(-x)', '1', 1 - c), ('Li3(1-x)', '1', (1 - c)/c),
    ('Li3(1/(1+x))', '1', -(1 + c)/c), ('Li3((1+x)/2)', '1', -2*c),
    ('Li3((1-x)/(1+x))', '1', 1 - c), ('Li3(-(1-x)/(1+x))', '1', 1 - c),
    ('Li3(2x/(1+x))', '1', -2*c),
    ('Li2(x)', '1', -1), ('Li2(-x)', '1', -1),
    ('Li2((1+x)/2)*l1p', '1', 1 + c), ('Li2((1+x)/2)*l1m', '1', -(1 - c)),
    ('Li2(x)*l1p', '1', q), ('Li2(x)*l1m', '1', q),
    ('Li2(-x)*l1p', '1', -q), ('Li2(-x)*l1m', '1', -q),
    ('l1m*l1m*l1p', '1', -(1 - c)), ('l1m*l1p*l1p', '1', 1 - c),
    ('l1p*l1p*l1p', '1', (4*c**2 - c + 1)/(6*c)),
    ('lx*l1m*l1m', '1', q), ('lx*l1m*l1p', '1', -(1 - c)), ('lx*l1p*l1p', '1', -(1 + c**2)/(2*c)),
    ('l1m*l1m', 'log2', 1 - c), ('l1m*l1p', 'log2', -(1 - c)), ('l1p*l1p', 'log2', -c),
    ('l1p', 'log2^2', (1 - c)/2), ('l1m', 'log2^2', -(1 - c)/2),
    ('l1m', 'pi2', -(1 - c)/(6*c)), ('l1p', 'pi2', -(1 + c)/(6*c)),
    ('1', 'z3', sp.Rational(7, 4) + 2*c), ('1', 'pi2*log2', -c/6), ('1', 'log2^3', c/3),
]
for fn, cn, v in T:
    add(typed, fn, cn, v)

keys = sorted(set(coef) | set(typed))
print(f"{'key':>34} | {'assembled':>22} | {'typed':>22} | diff")
allok = True
for k in keys:
    a = coef.get(k, 0); t = typed.get(k, 0)
    d = sp.simplify(a - t)
    ok = (d == 0)
    allok &= ok
    print(f"{str(k):>34} | {str(sp.factor(a)):>22} | {str(sp.factor(t)):>22} | {d}  {'OK' if ok else 'MISMATCH'}")
print("\nALL COEFFICIENTS MATCH:", allok)
