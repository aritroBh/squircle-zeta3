verify = r'''
#!/usr/bin/env python
# verify_cor6x.py -- COR6X referee: independent true point delta_{2,2} of Omega(c0)
# at n=5, n=8, n=12, cos(pi/12), and the square c0=1/sqrt2 (zero control).
# Compares with research/08_theory/COR6_true_point_symbols/cor6_table.json (read-only).
# Uses ONLY this project's own code (cor6x_core.py). Runtime target < 5 min.
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
TERMS = [(parse_term(s)[1], sp.Rational(c)) for s, c in data if len(parse_term(s)[1]) == 4]
assert len(TERMS) == 56

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

def sym_Omega_formal(c0):
    S = {}
    for w, cf in TERMS:
        S = sym_add(S, sym_scale(sym_H(w, c0), cf))
    return formalize(S, c0)

def verify_relation(c0, e):
    """check 2^e0 c0^e1 (1-c0)^e2 (1+c0)^e3 == +-1 exactly."""
    prod = sp.Integer(2) ** e[0] * c0 ** e[1] * (1 - c0) ** e[2] * (1 + c0) ** e[3]
    prod = sp.simplify(prod)
    assert prod in (1, -1), f"relation {e} at {c0} fails: {prod}"
    return prod

def delta22_image_hand(T, rel):
    """image of formal delta22 T in Lambda^2(Lambda^2(F4/R)), R = span(rel) (or empty).
       Returns canonical dict; computed by eliminating B via rel if given."""
    # We expand everything in the free letters A,C,D (if relation is a*A+b*B+c*C+d*D=0 with b!=0)
    def wedge_expand(l, subs):
        return [subs.get(x, x) for x in l]
    # generator substitution: each generator g -> linear form in free letters
    free = [A, C, D]
    if rel:
        a, b, c, d = rel
        assert b != 0
        subs = {A: A, B: sp.Rational(-a, b) * A + sp.Rational(-c, b) * C + sp.Rational(-d, b) * D, C: C, D: D}
    else:
        subs = {A: A, B: B, C: C, D: D}
        free = [A, B, C, D]
    def wedge_terms(g1, g2):
        # returns dict {((f1,f2)): coeff} for g1^g2 in free letters
        l1, l2 = subs[g1], subs[g2]
        out = {}
        for f1 in free:
            for f2 in free:
                c1 = l1.coeff(f1) if f1 in [A,B,C,D] else 0
                c2 = l2.coeff(f2) if f2 in [A,B,C,D] else 0
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
                out[key] = out.get(key, 0) + sgn * v1 * v2 * v
    return {k: v for k, v in out.items() if v}

# point, minpoly degree, relations (exactly verified below), cor6 label
POINTS = [
    ("n=5  c0=sin(pi/5)",  sp.sin(sp.pi/5),  4, [],                       "n=5 c=sin(pi/5)"),
    ("n=8  c0=sin(pi/8)",  sp.sin(sp.pi/8),  4, [[3, 2, 1, 1]],          "n=8 octagon c=sin(pi/8)=cos(3pi/8)"),
    ("n=12 c0=sin(pi/12)", sp.sin(sp.pi/12), 4, [[4, 2, 1, 1]],          "n=12 c=sin(pi/12)"),
    ("a=pi/6 c0=cos(pi/12)", sp.cos(sp.pi/12), 4, [[4, 2, 1, 1]],        "c=cos(pi/12) (alpha=pi/6)"),
    ("n=4  c0=1/sqrt2",    1/sp.sqrt(2),     2, [[1, 0, 1, 1], [1, 2, 0, 0]], "n=4 square c=1/sqrt2"),
]

cor6 = {e['label']: e for e in json.load(open(COR6))}
print(f"{'point':26s} {'deg':>3s} {'dimU':>4s} {'true d22':>22s} {'naive d22':>10s} {'cor6':>18s} {'verdict'}")
summary = []
for name, c0, deg, rels, clabel in POINTS:
    for e in rels:
        verify_relation(c0, e)
    S = sym_Omega_formal(c0)
    T = delta22(S)
    if len(rels) == 2:
        # two relations -> dim U = 2 -> delta22 = 0 automatically
        img = {}
        dimU = 2
    else:
        img = delta22_image_hand(T, rels[0] if rels else None)
        dimU = 4 - len(rels)
    # naive: drop monomials containing A (constant letter [2])
    Sn = {k: v for k, v in S.items() if A not in k}
    Tn = delta22(Sn)
    imgn = delta22_image_hand(Tn, rels[0] if len(rels) == 1 else None) if len(rels) <= 1 else {}
    e = cor6[clabel]
    cor6_true = "ZERO" if e['true_zero'] else "NONZERO"
    my_true = "ZERO" if not img else f"NONZERO({len(img)})"
    cor6_naive = "ZERO" if e['naive_zero'] else "NONZERO"
    my_naive = "ZERO" if not imgn else "NONZERO"
    ok = ((e['deg'] == deg) and (e['letter_span_dim'] == dimU)
          and ((e['true_zero']) == (not img)) and ((e['naive_zero']) == (not imgn)))
    verdict = "AGREE" if ok else "DISAGREE"
    summary.append((name, verdict))
    print(f"{name:26s} {deg:>3d} {dimU:>4d} {my_true:>22s} {my_naive:>10s} {cor6_true+'/'+cor6_naive:>18s} {verdict}")
    print(f"    relations used (exactly verified): {rels}")
    if img:
        for k, v in sorted(img.items(), key=str)[:6]:
            print("      d22 component:", k, v)
print()
print("SUMMARY:", ", ".join(f"{n.split()[0]}:{v}" for n, v in summary))
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/verify_cor6x.py","w").write(verify)
import subprocess, time
t0 = time.time()
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","verify_cor6x.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=400)
print(r.stdout); print("ERR:", (r.stderr or "")[-1500:]); print("runtime: %.1fs" % (time.time()-t0))
