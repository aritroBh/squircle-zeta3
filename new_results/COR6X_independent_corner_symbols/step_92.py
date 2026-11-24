t = r'''
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
from functools import lru_cache
import sympy as sp, json, re

EPS = sp.Symbol('EPSREG9')
def make_sym(scheme):
    def cn(e):
        e = sp.simplify(e)
        if e is EPS: return e
        ne = sp.simplify(-e)
        ke, kn = (str(e).startswith('-'), str(e)), (str(ne).startswith('-'), str(ne))
        return e if ke <= kn else ne
    def isz(e):
        e = sp.simplify(e)
        return e == 0 or e in (1, -1, sp.I, -sp.I)
    def dl(P, Q):
        # P,Q: ('m', base, i) or ('c', base)
        if P[0] == 'm' and Q[0] == 'm':
            d = sp.simplify(P[1] - Q[1])
            if d == 0:
                # collision: eps^i - eps^j or i*eps - j*eps
                i, j = P[2], Q[2]
                if scheme == 'epspow':
                    return EPS
                else:  # 'ieps': (i-j)*eps -> [|i-j|] + [eps] -> keep |i-j|, drop eps
                    return sp.Integer(abs(i - j)) if abs(i - j) > 1 else EPS
            return cn(d)
        if P[0] == 'm' and Q[0] == 'c':
            d = sp.simplify(P[1] - Q[1])
            if d == 0:
                return EPS
            return cn(d)
        if P[0] == 'c' and Q[0] == 'm':
            d = sp.simplify(P[1] - Q[1])
            if d == 0:
                return EPS
            return cn(d)
        d = sp.simplify(P[1] - Q[1])
        if d == 0: raise AssertionError("collide")
        return cn(d)
    @lru_cache(maxsize=None)
    def symS(a0, letters, aN):
        n = len(letters)
        A0 = ('c', a0); AN = ('c', aN)
        LL = tuple(('m', a, i + 1) for i, a in enumerate(letters))
        if n == 0: return {(): 1}
        out = {}
        for i in range(n):
            sub = tuple(l[1] for l in LL[:i] + LL[i+1:])
            Ssub = symS(a0, sub, aN)
            prev = AN if i == 0 else LL[i-1]
            nxt = A0 if i == n - 1 else LL[i+1]
            e_prev = dl(prev, LL[i])
            e_next = dl(LL[i], nxt)
            for k, v in Ssub.items():
                if not (e_prev is EPS or isz(e_prev)):
                    kk = k + (e_prev,); out[kk] = out.get(kk, 0) + v
                if not (e_next is EPS or isz(e_next)):
                    kk = k + (e_next,); out[kk] = out.get(kk, 0) - v
        return {k: v for k, v in out.items() if v and not any(t is EPS for t in k)}
    return symS

def sH(symS, word, pt):
    sgn = (-1) ** sum(1 for a in word if a == 1)
    return {k: sgn * v for k, v in symS(sp.Integer(0), tuple(sp.sympify(a) for a in word), sp.sympify(pt)).items()}

data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[1], sp.Rational(c)) for s, c in data if len(parse_term(s)[1]) == 4]
c0 = sp.Rational(1, 2)

for scheme in ('epspow', 'ieps'):
    symS = make_sym(scheme)
    S_om = {}
    for w, cf in TERMS:
        S_om = sym_add(S_om, sym_scale(sH(symS, w, c0), cf))
    S_om = sym_expanded(S_om)
    # reference: tot = S_five + S_P from before (recompute minimal: reuse my sym_G versions)
    def sym_Li4(z): return sym_scale(sym_G((0,0,0,1), z), -1)
    S_five = sym_expanded(sym_add(
        sym_scale(sym_Li4(-2), sp.Rational(-21)),
        sym_scale(sym_Li4(sp.Rational(1,3)), sp.Rational(-65,8)),
        sym_scale(sym_Li4(sp.Rational(2,3)), sp.Rational(13,2)),
        sym_scale(sym_Li4(sp.Rational(1,2)), sp.Rational(137,6)),
        sym_scale(sym_Li4(sp.Rational(1,9)), sp.Rational(39,64))))
    S_log2 = {(sp.Integer(2),): 1}; S_log3 = {(sp.Integer(3),): 1}
    S_Li3_13 = sym_expanded(sym_scale(sym_G((0,0,1), sp.Rational(1,3)), -1))
    S_Li3_m12 = sym_expanded(sym_scale(sym_G((0,0,1), sp.Rational(-1,2)), -1))
    S_Li2_13 = sym_expanded(sym_scale(sym_G((0,1), sp.Rational(1,3)), -1))
    def sh(*args):
        r = args[0]
        for a in args[1:]: r = sym_shuffle(r, a)
        return r
    S_P = sym_add(
        sym_scale(sh(S_Li3_13, S_log2), sp.Rational(15,2)),
        sym_scale(sh(S_Li3_m12, S_log2), sp.Rational(37,2)),
        sym_scale(sh(S_Li3_m12, S_log3), sp.Rational(-5,4)),
        sh(S_Li2_13, S_Li2_13),
        sym_scale(sh(S_Li2_13, S_log2, S_log2), sp.Rational(-23,4)),
        sym_scale(sh(S_Li2_13, S_log2, S_log3), sp.Rational(7,4)),
        sh(S_Li2_13, S_log3, S_log3),
        sym_scale(sh(S_log2, S_log2, S_log2, S_log2), sp.Rational(-211,144)),
        sym_scale(sh(S_log2, S_log2, S_log2, S_log3), sp.Rational(101,24)),
        sym_scale(sh(S_log2, S_log2, S_log3, S_log3), sp.Rational(-15,8)),
        sym_scale(sh(S_log2, S_log3, S_log3, S_log3), sp.Rational(-35,24)),
        sym_scale(sh(S_log3, S_log3, S_log3, S_log3), sp.Rational(113,192)))
    S_P = sym_expanded(S_P)
    T = sym_sub(S_om, sym_add(S_five, S_P))
    print(f"scheme {scheme}: T terms = {len(T)}")
    for k, v in sorted(T.items(), key=str): print("   ", k, v)
'''
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=600)
print(r.stdout); print("ERR:", (r.stderr or "")[-1000:])
