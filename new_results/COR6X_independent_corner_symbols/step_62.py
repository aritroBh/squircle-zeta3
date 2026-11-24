t = r'''
# step_07.py -- S(Omega(c0)) at cyclotomic points, letters as formal generators A,B,C,D
import sys; sys.path.insert(0, '.')
from cor6x_core import *
import sympy as sp, json, re
data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[1], sp.Rational(c)) for s, c in data if len(parse_term(s)[1]) == 4]

A, B, C, D = sp.symbols('A B C D')
def formalize(S, c0):
    """map every letter to A=[2], B=[c0], C=[1-c0], D=[1+c0]; verify coverage."""
    out = {}
    for mono, v in S.items():
        nm = []
        for t in mono:
            tt = sp.nsimplify(sp.simplify(t)) if not t.free_symbols else t
            if sp.simplify(tt - 2) == 0: nm.append(A)
            elif sp.simplify(tt - c0) == 0: nm.append(B)
            elif sp.simplify(tt - (1 - c0)) == 0: nm.append(C)
            elif sp.simplify(tt - (1 + c0)) == 0: nm.append(D)
            else:
                raise ValueError(f"unexpected letter {t} at c0={c0}")
        out[tuple(nm)] = out.get(tuple(nm), 0) + v
    return {k: v for k, v in out.items() if v}

import sympy as sp
results = {}
for n in list(range(3, 13)):
    c0 = sp.sin(sp.pi/n)
    S = {}
    for w, cf in TERMS:
        S = sym_add(S, sym_scale(sym_H(w, c0), cf))
    F = formalize(S, c0)
    Dlt = delta22(F)
    results[n] = (c0, F, Dlt)
    print(f"n={n}: c0=sin(pi/{n})  symbol terms {len(F)}  delta22 formal terms {len(Dlt)}")
    if Dlt:
        for k, v in sorted(Dlt.items(), key=str):
            print("     ", k, v)
# cos(pi/12)
c0 = sp.cos(sp.pi/12)
S = {}
for w, cf in TERMS:
    S = sym_add(S, sym_scale(sym_H(w, c0), cf))
F = formalize(S, c0)
Dlt = delta22(F)
print(f"c0=cos(pi/12): symbol terms {len(F)}  delta22 formal terms {len(Dlt)}")
for k, v in sorted(Dlt.items(), key=str):
    print("     ", k, v)
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/step_07.py","w").write(t)
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_07.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=1500)
print(r.stdout); print("ERR:", (r.stderr or "")[-2000:])
