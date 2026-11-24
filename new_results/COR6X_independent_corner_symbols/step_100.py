t = r'''
# step_14.py -- full 14-row delta22 table WITH c0^e scaling (the fix)
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
import sympy as sp, json, re, itertools
data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
W4 = [(parse_term(s)[0], parse_term(s)[1], sp.Rational(c)) for s, c in data if len(parse_term(s)[1]) == 4]

A, B, C, D = sp.symbols('A B C D')
def formalize(S, c0):
    out = {}
    for mono, v in S.items():
        nm = []
        for t in mono:
            def eq(u): return sp.simplify(t - u) == 0 or sp.simplify(t + u) == 0
            if eq(sp.Integer(2)): nm.append(A)
            elif eq(c0): nm.append(B)
            elif eq(1 - c0): nm.append(C)
            elif eq(1 + c0): nm.append(D)
            else: raise ValueError(f"unexpected letter {t}")
        out[tuple(nm)] = out.get(tuple(nm), 0) + v
    return {k: v for k, v in out.items() if v}

def sym_Omega_formal_scaled(c0):
    """formal symbol of Omega(c0) WITH c0^e scaling (kept as symbolic scalar)."""
    Sg = {}
    for e, w, cf in W4:
        Sg[e] = sym_add(Sg.get(e, {}), sym_scale(sym_H(w, c0), cf))
    # scale each group by c0^e and mark coefficients symbolically
    tot = {}
    for e, S in Sg.items():
        ce = c0 ** e
        F = formalize(S, c0)
        for k, v in F.items():
            tot[k] = tot.get(k, 0) + v * ce
    return tot

def delta22_scaled(T):
    """delta22 of formal symbol T whose coefficients may be scalars in Q(c0)."""
    out = {}
    for (a, b, c, d), v in T.items():
        w1, s1 = _wedge(a, b); w2, s2 = _wedge(c, d)
        if w1 is None or w2 is None or w1 == w2: continue
        s = s1 * s2
        key = (w1, w2) if str(w1) < str(w2) else (w2, w1)
        if key != (w1, w2): s = -s
        out[key] = sp.simplify(out.get(key, 0) + s * v)
    return {k: v for k, v in out.items() if sp.simplify(v) != 0}

RELS = {
    "sin(pi/3)": [[2, 0, 1, 1]], "sin(pi/4)": [[1, 0, 1, 1], [1, 2, 0, 0]],
    "sin(pi/5)": [], "sin(pi/6)": [[1, 1, 0, 0], [1, 0, 1, 0]],
    "sin(pi/7)": [], "sin(pi/8)": [[3, 2, 1, 1]], "sin(pi/9)": [],
    "sin(pi/10)": [[3, 2, 0, 1]], "sin(pi/11)": [], "sin(pi/12)": [[4, 2, 1, 1]],
    "cos(pi/12)": [[4, 2, 1, 1]],
}
def delta22_image_hand(T, rel):
    free = [A, C, D]
    if rel:
        a, b, c, d = rel
        subs = {A: A, B: sp.Rational(-a, b) * A + sp.Rational(-c, b) * C + sp.Rational(-d, b) * D, C: C, D: D}
    else:
        subs = {A: A, B: B, C: C, D: D}
        free = [A, B, C, D]
    def wedge_terms(g1, g2):
        l1, l2 = subs[g1], subs[g2]
        out = {}
        for f1 in free:
            for f2 in free:
                c1 = l1.coeff(f1); c2 = l2.coeff(f2)
                v = c1 * c2
                if v:
                    key = (f1, f2) if str(f1) < str(f2) else (f2, f1)
                    sgn = 1 if str(f1) < str(f2) else -1
                    if key[0] != key[1]:
                        out[key] = out.get(key, 0) + sgn * v
        return {k: v for k, v in out.items() if v}
    out = {}
    for ((g1, g2), (g3, g4)), v in T.items():
        w1 = wedge_terms(g1, g2); w2 = wedge_terms(g3, g4)
        for k1, v1 in w1.items():
            for k2, v2 in w2.items():
                if k1 == k2: continue
                key = (k1, k2) if str(k1) < str(k2) else (k2, k1)
                sgn = 1 if str(k1) < str(k2) else -1
                out[key] = sp.simplify(out.get(key, 0) + sgn * v1 * v2 * v)
    return {k: v for k, v in out.items() if sp.simplify(v) != 0}

POINTS = [(f"sin(pi/{n})", sp.sin(sp.pi/n)) for n in range(3, 13)] + [("cos(pi/12)", sp.cos(sp.pi/12))]
for name, c0 in POINTS:
    rels = RELS[name]
    T = sym_Omega_formal_scaled(c0)
    Td = delta22_scaled(T)
    if len(rels) == 2:
        img = {}
        dimU = 2
    else:
        img = delta22_image_hand(Td, rels[0] if rels else None)
        dimU = 4 - len(rels)
    print(f"{name:12s} dimU={dimU}  scaled delta22: {'ZERO' if not img else f'NONZERO ({len(img)} comp.)'}")
    if img:
        for k, v in sorted(img.items(), key=str): print("     ", k, v)
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/step_14.py","w").write(t)
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_14.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=1500)
print(r.stdout); print("ERR:", (r.stderr or "")[-1500:])
