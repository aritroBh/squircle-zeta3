core = open("<project>/research/08_theory/COR6X_independent_corner_symbols/cor6x_core.py").read()

addition = '''

# ---------- letter expansion into a prime basis (for rational-function letters) ----------
def expand_letter(e):
    """e: sympy rational expression. Returns dict {irreducible factor (up to sign): exponent}.
       Units/torsion dropped. Works for rational functions over Q."""
    e = sp.nsimplify(e) if not e.free_symbols else e
    num, den = sp.fraction(sp.cancel(e))
    out = {}
    for poly, sgn in ((num, 1), (den, -1)):
        poly = sp.Poly(poly)
        content, fac = sp.factor_list(poly.as_expr())
        # rational content: factor into primes
        for b, m in sp.factorint(sp.Rational(content)).items():
            out[sp.Integer(b)] = out.get(sp.Integer(b), 0) + sgn * m
        for f, m in fac:
            f = sp.primitive(f)[1] if f.is_Poly else f
            fc = sp.Poly(f)
            # canonicalize up to sign: make leading coeff positive
            f = sp.Poly(f, *fc.gens)
            lc = sp.LC(f)
            ff = sp.primitive(f.as_expr())[1]
            ff = sp.factor(ff)
            if sp.simplify(ff - 1) == 0:
                continue
            if str(sp.LC(sp.Poly(ff, *fc.gens)))[:1] == '-':
                ff = -ff
            out[ff] = out.get(ff, 0) + sgn * m
    return {k: v for k, v in out.items() if v}

def sym_expanded(S):
    """expand every letter of the symbol into the prime basis; returns symbol over prime letters."""
    out = {}
    for mono, c in S.items():
        exps = [expand_letter(t) for t in mono]
        if any(len(e) == 0 for e in exps):
            continue
        choices = []
        for e in exps:
            terms = [(k,) for k, v in e.items() for _ in range(0)] or None
            choices.append(list(e.items()))
        import itertools as _it
        for combo in _it.product(*[[(k, v) for k, v in e.items()] for e in exps]):
            newmono = []
            mult = c
            for k, v in combo:
                newmono.extend([k] * abs(v))
                if v < 0:
                    mult *= -1  # WRONG placeholder; sign handled below
            raise RuntimeError("use expand_signed")
    return out

def sym_expanded2(S):
    """same as sym_expanded but keeping integer exponents: letters may repeat with signs
       via additive notation: monomial entries are (prime, exponent) collapsed additively:
       we represent the expanded symbol as dict over tuples of primes with the understanding
       [p^e] = e*[p], so a monomial maps to a sum of tensors; coefficients multiply out."""
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
# careful: in sym_expanded2, [t] = sum_e e*[p_e]; tensor product expands multilinearly; coefficient multiplies by each exponent e. Correct.
open("<project>/research/08_theory/COR6X_independent_corner_symbols/cor6x_core.py","a").write(addition)

t = r'''
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
import sympy as sp
x, y = sp.symbols('x y')
def B(z): return delta11(sym_expanded2(sym_G((0,1), z)))
z1 = (1-x)/(1-x*y); z2 = (1-y)/(1-x*y)
tot = {}
for z in (x, y, z1, z2, 1-x*y):
    tot = sym_add(tot, B(z))
print("five-term all-plus (expanded):", tot)
tot2 = {}
for z, s in ((x,1),(y,1),(z1,1),(z2,1),(1-x*y,-1)):
    tot2 = sym_add(tot2, sym_scale(B(z), s))
print("five-term (++ ++ -) (expanded):", tot2)
# also Suslin form
tot3 = {}
for z, s in ((x,1),(y,-1),(y/x,1),((1-y)/(1-x),-1),(x*(1-y)/(y*(1-x)),1)):
    tot3 = sym_add(tot3, sym_scale(B(z), s))
print("Suslin form (expanded):", tot3)
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols")
print(r.stdout); print("ERR:", (r.stderr or "")[-800:])
