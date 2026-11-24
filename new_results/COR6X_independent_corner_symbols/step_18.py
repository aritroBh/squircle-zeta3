path = "<project>/research/08_theory/COR6X_independent_corner_symbols/cor6x_core.py"
core = open(path).read()
# remove the broken block
i = core.find("# ---------- letter expansion into a prime basis")
core = core[:i]
addition = '''# ---------- letter expansion into a prime basis (rational-function letters) ----------
def expand_letter(e):
    """e: sympy rational expression over Q. Returns {irreducible factor (canon. sign): exp}.
       Torsion (+-1) dropped."""
    num, den = sp.fraction(sp.cancel(e))
    out = {}
    for poly_expr, sgn in ((num, 1), (den, -1)):
        if not poly_expr.free_symbols:
            r = sp.Rational(poly_expr)
            for b, m in sp.factorint(r.p).items():
                out[sp.Integer(b)] = out.get(sp.Integer(b), 0) + sgn * m
            for b, m in sp.factorint(r.q).items():
                out[sp.Integer(b)] = out.get(sp.Integer(b), 0) - sgn * m
            continue
        content, fac = sp.factor_list(poly_expr)
        r = sp.Rational(content)
        for b, m in sp.factorint(r.p).items():
            out[sp.Integer(b)] = out.get(sp.Integer(b), 0) + sgn * m
        for b, m in sp.factorint(r.q).items():
            out[sp.Integer(b)] = out.get(sp.Integer(b), 0) - sgn * m
        for f, m in fac:
            f0 = sp.primitive(f.as_expr())[1]
            if sp.simplify(f0 - 1) == 0:
                continue
            if sp.LC(sp.Poly(f0)) < 0:
                f0 = -f0
            out[f0] = out.get(f0, 0) + sgn * m
    return {k: v for k, v in out.items() if v}

def sym_expanded(S):
    """multilinear expansion of every letter into the prime basis ([prod p^e] = sum e [p])."""
    import itertools as _it
    out = {}
    for mono, c in S.items():
        exps = [expand_letter(t) for t in mono]
        if any(len(e) == 0 for e in exps):
            continue
        for combo in _it.product(*[list(e.items()) for e in exps]):
            newmono = tuple(k for k, v in combo)
            mult = c
            for k, v in combo:
                mult *= v
            out[newmono] = out.get(newmono, 0) + mult
    return {k: v for k, v in out.items() if v}
'''
open(path, "w").write(core + addition)

t = r'''
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
import sympy as sp
x, y = sp.symbols('x y')
def B(z): return delta11(sym_expanded(sym_G((0,1), z)))
z1 = (1-x)/(1-x*y); z2 = (1-y)/(1-x*y)
tot = {}
for z in (x, y, z1, z2, 1-x*y):
    tot = sym_add(tot, B(z))
print("five-term all-plus (expanded):", tot == {})
tot3 = {}
for z, s in ((x,1),(y,-1),(y/x,1),((1-y)/(1-x),-1),(x*(1-y)/(y*(1-x)),1)):
    tot3 = sym_add(tot3, sym_scale(B(z), s))
print("Suslin form (expanded):", tot3 == {})
# dilog reflection: Li_2(x)+Li_2(1-x) is a product
print("Li_2(x)+Li_2(1-x) product-symmetric:", sym_sub(sym_add(sym_G((0,1),x), sym_G((0,1),1-x)),
      sym_shuffle(sym_G((0,),x), sym_G((1,),x))) == {})
# inversion: Li_2(x)+Li_2(1/x) = -zeta2 - (log^2(-x))/2: symbol-level:
lhs = sym_add(sym_G((0,1), x), sym_G((0,1), 1/x))
print("Li_2(x)+Li_2(1/x) symbol:", lhs, " delta11:", delta11(sym_expanded(lhs)))
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols")
print(r.stdout); print("ERR:", (r.stderr or "")[-800:])
