
#!/usr/bin/env python
# verify_cor6x.py -- COR6X referee, FINAL: true point delta_{2,2} of Omega(c0) with the
# CORRECT treatment of the c^e prefactors (the symbol map is linear over the field:
# weight-0 constants c0^e act as scalars, they do NOT drop out).
# Compares all 14 rows with cor6_table.json (read-only). Own code only. < 5 min.
import sys, os, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cor6x_core import *
import sympy as sp

COR6 = "<project>/research/08_theory/COR6_true_point_symbols/cor6_table.json"
DATA = "<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"
data = json.load(open(DATA))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
W4 = [(parse_term(s)[0], parse_term(s)[1], sp.Rational(c)) for s, c in data if len(parse_term(s)[1]) == 4]
assert len(W4) == 56

def sym_Omega_scaled(c0):
    S = {}
    for e, w, cf in W4:
        S = sym_add(S, sym_scale(sym_H(w, c0), cf * (c0 ** e)))
    return S

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

RELS = {"sin(pi/3)": [[2, 0, 1, 1]], "sin(pi/4)": [[1, 0, 1, 1], [1, 2, 0, 0]],
        "sin(pi/5)": [], "sin(pi/6)": [[1, 1, 0, 0], [1, 0, 1, 0]], "sin(pi/7)": [],
        "sin(pi/8)": [[3, 2, 1, 1]], "sin(pi/9)": [], "sin(pi/10)": [[3, 2, 0, 1]],
        "sin(pi/11)": [], "sin(pi/12)": [[4, 2, 1, 1]], "cos(pi/12)": [[4, 2, 1, 1]]}

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
                v = l1.coeff(f1) * l2.coeff(f2)
                if v:
                    key = (f1, f2) if str(f1) < str(f2) else (f2, f1)
                    sgn = 1 if str(f1) < str(f2) else -1
                    if key[0] != key[1]:
                        out[key] = out.get(key, 0) + sgn * v
        return {k: v for k, v in out.items() if v}
    out = {}
    for ((g1, g2), (g3, g4)), v in T.items():
        for k1, v1 in wedge_terms(g1, g2).items():
            for k2, v2 in wedge_terms(g3, g4).items():
                if k1 == k2: continue
                key = (k1, k2) if str(k1) < str(k2) else (k2, k1)
                sgn = 1 if str(k1) < str(k2) else -1
                out[key] = sp.simplify(out.get(key, 0) + sgn * v1 * v2 * v)
    return {k: v for k, v in out.items() if sp.simplify(v) != 0}

cor6 = {e['label']: e for e in json.load(open(COR6))}
ROWS = [
    ("sin(pi/3)", sp.sin(sp.pi/3), 2, "n=3 triangle c=sqrt3/2"),
    ("sin(pi/4)", 1/sp.sqrt(2), 2, "n=4 square c=1/sqrt2"),
    ("sin(pi/5)", sp.sin(sp.pi/5), 4, "n=5 c=sin(pi/5)"),
    ("sin(pi/6)", sp.Rational(1,2), 1, "n=6 hexagon c=1/2"),
    ("sin(pi/7)", sp.sin(sp.pi/7), 6, "n=7 c=sin(pi/7)"),
    ("sin(pi/8)", sp.sin(sp.pi/8), 4, "n=8 octagon c=sin(pi/8)=cos(3pi/8)"),
    ("sin(pi/9)", sp.sin(sp.pi/9), 6, "n=9 c=sin(pi/9)"),
    ("sin(pi/10)", sp.sin(sp.pi/10), 2, "n=10 c=sin(pi/10)"),
    ("sin(pi/11)", sp.sin(sp.pi/11), 10, "n=11 c=sin(pi/11)"),
    ("sin(pi/12)", sp.sin(sp.pi/12), 4, "n=12 c=sin(pi/12)"),
    ("cos(pi/12)", sp.cos(sp.pi/12), 4, "c=cos(pi/12) (alpha=pi/6)"),
    ("1/3", sp.Rational(1,3), 1, "c=1/3"),
    ("2/5", sp.Rational(2,5), 1, "c=2/5"),
    ("3/7", sp.Rational(3,7), 1, "c=3/7"),
]

print("controls: delta22(Li4(x)) =", delta22(sym_expanded(sym_scale(sym_G((0,0,0,1), sp.Symbol('x')), -1))),
      " delta22(Li_2(x)^2) =", delta22(sym_expanded(sym_shuffle(sym_G((0,1), sp.Symbol('x')), sym_G((0,1), sp.Symbol('x'))))))
print(f"{'point':13s} {'deg':>3s} {'dimU':>4s} {'my true d22':>18s} {'cor6 true_zero':>14s} {'verdict'}")
summary = []
for name, c0, deg, clabel in ROWS:
    rels = RELS.get(name, [])
    S = sym_Omega_scaled(c0)
    if c0.is_rational:
        Dlt = delta22(sym_expanded(S))
        img = Dlt
        dimU = 4 - (sp.Matrix(rels).rank() if rels else 0) if name in ("1/3",) else (2 if c0 in (sp.Rational(1,3),) else None)
        # rational points: dims from analysis
        dimU = {"1/3": 2, "2/5": 4, "3/7": 4, "sin(pi/6)": 2}[name] if name in ("1/3","2/5","3/7","sin(pi/6)") else None
    else:
        T = formalize(S, c0)
        from cor6x_core import _wedge
        Td = {}
        for (a, b, c, d), v in T.items():
            w1, s1 = _wedge(a, b); w2, s2 = _wedge(c, d)
            if w1 is None or w2 is None or w1 == w2: continue
            s = s1 * s2
            key = (w1, w2) if str(w1) < str(w2) else (w2, w1)
            if key != (w1, w2): s = -s
            Td[key] = sp.simplify(Td.get(key, 0) + s * v)
        Td = {k: v for k, v in Td.items() if sp.simplify(v) != 0}
        if len(rels) == 2:
            img, dimU = {}, 2
        else:
            img = delta22_image_hand(Td, rels[0] if rels else None)
            dimU = 4 - len(rels)
    e = cor6[clabel]
    my_true = "ZERO" if not img else "NONZERO"
    cor6_true = "ZERO" if e['true_zero'] else "NONZERO"
    ok = (e['deg'] == deg and e['letter_span_dim'] == dimU and (e['true_zero'] == (not img)))
    verdict = "AGREE" if ok else "DISAGREE"
    summary.append((name, verdict))
    print(f"{name:13s} {deg:>3d} {dimU:>4d} {my_true+'('+str(len(img))+' comp.)':>18s} {cor6_true:>14s} {verdict}")

# value-level check at rational points 2/5, 3/7
prime_order = [sp.Integer(2), sp.Integer(3), sp.Integer(5), sp.Integer(7)]
for c0, label in ((sp.Rational(2,5), 'c=2/5'), (sp.Rational(3,7), 'c=3/7')):
    D = delta22(sym_expanded(sym_Omega_scaled(c0)))
    mine = {f"(({prime_order.index(a)}, {prime_order.index(b)}), ({prime_order.index(c)}, {prime_order.index(d)}))": v
            for ((a, b), (c, d)), v in D.items()}
    theirs = {k: v[0] for k, v in cor6[label]['true_d22'].items()}
    keys = set(mine) | set(theirs)
    diffs = [k for k in keys if sp.Rational(mine.get(k, 0)) != sp.Rational(theirs.get(k, 0))]
    print(f"value check {label}: {'MATCH EXACTLY' if not diffs else f'{len(diffs)} MISMATCHES'}")

# Lie-class confirmations (fixed scaling)
import itertools
def decomposable(T, primes):
    basis = list(itertools.product(primes, repeat=4))
    idx = {m: i for i, m in enumerate(basis)}
    rows = []
    for u in (1, 2, 3):
        for m1 in itertools.product(primes, repeat=u):
            for m2 in itertools.product(primes, repeat=4 - u):
                r = [0]*len(basis)
                for w in shuffle_words(m1, m2): r[idx[w]] += 1
                rows.append(r)
    P = sp.Matrix(rows)
    dv = sp.Matrix([T.get(m, 0) for m in basis])
    return P.T.rank() == P.T.row_join(dv).rank()
def sym_Li4(z): return sym_scale(sym_G((0,0,0,1), z), -1)
d3 = sym_sub(sym_expanded(sym_Omega_scaled(sp.Rational(1,3))),
             sym_expanded(sym_add(sym_scale(sym_Li4(-8), sp.Rational(-11,81)),
                                  sym_scale(sym_Li4(-3), sp.Rational(-2)),
                                  sym_scale(sym_Li4(sp.Rational(2,3)), sp.Rational(20,3)))))
print("Lie class 1/3: Omega(1/3) == -(11/81)Li4(-8)-2Li4(-3)+(20/3)Li4(2/3) mod shuffles:",
      decomposable(d3, [sp.Integer(2), sp.Integer(3)]))
d2 = sym_sub(sym_expanded(sym_Omega_scaled(sp.Rational(1,2))),
             sym_expanded(sym_add(sym_scale(sym_Li4(-2), sp.Rational(-21)),
                                  sym_scale(sym_Li4(sp.Rational(1,3)), sp.Rational(13,2)),
                                  sym_scale(sym_Li4(sp.Rational(2,3)), sp.Rational(13,2)))))
print("Lie class 1/2: Omega(1/2) == -21Li4(-2)+(13/2)Li4(1/3)+(13/2)Li4(2/3) mod shuffles:",
      decomposable(d2, [sp.Integer(2), sp.Integer(3)]))
print()
print("SUMMARY:", ", ".join(f"{n.split('(')[0]}:{v}" for n, v in summary))
