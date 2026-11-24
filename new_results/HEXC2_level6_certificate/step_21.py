code = '''"""HEXC2 shared machinery: exact conversion of the hexagon-corner identities to level-6 CMZV words.

The Mobius map f(z) = (zeta^2 z + 1)/(zeta^4 z + 1), zeta = e^{i pi/3}, sends
  0 |-> 1, 1/2 |-> zeta, 1 |-> zeta^2, -1 |-> zeta^5, inf |-> zeta^4  (exponents (0,1,2,5,4)).
Under g = f^{-1} the HPL letters {0, +/-1} (forms dz/z, dz/(1-z), dz/(1+z)) pull back to
  0  |-> (2) - (0)     [dlog(u-1) - dlog(1-zeta^2 u) = -w0 + w2]
  1  |-> (4) - (2)     [-dlog(g-1) = -(w2 - w4) = w4 - w2]
  -1 |-> (2) - (1)     [dlog(g+1) = -w1 + w2]
where (j) = w_j = zeta^j du/(1 - zeta^j u) is cmzv6 letter j.  Every HPL word H_w(1/2)
(first letter outermost, integral from 0 to 1/2) becomes an iterated integral from 1 to zeta,
converted to level-6 CMZV words (integrals 0 -> 1) by path composition through 0 + rescaling.

General map data: for a map with exponents E = (e0, eh, e1, em, einf) [f(0),f(1/2),f(1),f(-1),f(inf)],
  g*dlog(x - q) = dlog(u - zeta^{e_q}) - dlog(u - zeta^{e_inf}) = -letter(-e_q) + letter(-e_inf)
  (using dlog(u - zeta^k) = -letter((-k) mod 6));  letter(-1) = dt/t appears iff e_q or e_inf == 0... 
  actually dlog(u - 1) = -letter(0) always (zeta^0=1).  dlog(u - zeta^0) never gives dt/t; dt/t would
  come from zeta^k = 0, impossible -- so all letters stay in {0..5} for every one of the 12 maps.
"""
import sys
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
from functools import cache
from itertools import product
BASE = "<project>/research/08_theory"
sys.path.insert(0, BASE + "/HEXC1_hexagon_certificate")
from cmzv6 import sh, mul, reg, regularize, add, words, N

INF = None
# exponent images (f(0), f(1/2), f(1), f(-1), f(inf)) of the 12 Mobius maps (step_02)
MAPS = [
    (0,1,2,5,4), (0,5,4,1,2), (1,2,3,0,5), (1,0,5,2,3),
    (2,1,0,3,4), (2,3,4,1,0), (3,2,1,4,5), (3,4,5,2,1),
    (4,5,0,3,2), (4,3,2,5,0), (5,0,1,4,3), (5,4,3,0,1),
]
QPTS = ['0','1/2','1','-1','inf']

def dlog_pullback(qexp, einf):
    """g*dlog(x-q) as letter dict:  -letter(-qexp) + letter(-einf); returns {word: coef}."""
    return {((-qexp) % 6,): Fr(-1), ((-einf) % 6,): Fr(1)}

def hpl_letter_map(mi, sig):
    """Pullback of HPL letter sig in {0,1,-1} under map mi (task map = 0)."""
    e0, eh, e1, em, einf = MAPS[mi]
    if sig == 0:   # dz/z = dlog(x-0)
        return dlog_pullback(e0, einf)
    if sig == 1:   # dz/(1-z) = -dlog(x-1)
        return {w: -c for w, c in dlog_pullback(e1, einf).items()}
    if sig == -1:  # dz/(1+z) = dlog(x-(-1))
        return dlog_pullback(em, einf)
    raise ValueError(sig)

def shift(w, s):
    return tuple((-1 if j == -1 else (j + s) % 6) for j in w)

@cache
def path6(w, c, e):
    """I(w; zeta^c -> zeta^e), letters in {0..5} (+ -1), via 0, shuffle-regularized (T=0)."""
    w = tuple(w)
    if c == e:
        return {} if w else {(): Fr(1)}
    out = {}
    for k in range(len(w) + 1):
        left = regularize({(shift(tuple(reversed(w[:k])), c),): Fr(1)})
        right = regularize({(shift(w[k:], e),): Fr(1)})
        out = add((out, 1), (mul(left, right), Fr((-1) ** k)))
    return regularize(out)

def forms_path6(forms, c, e):
    """iterated integral of a tuple of letter-dicts (form combos) from zeta^c to zeta^e."""
    out = {}
    for choices in product(*forms):
        w = tuple(x for x, a in choices)
        coef = Fr(1)
        for x, a in choices: coef *= a
        out = add((out, 1), (path6(w, c, e), coef))
    return out

@cache
def hpl_to_cmzv(w, mi=0):
    """H_w(1/2) (HPL letters 0,+/-1, first letter outermost) as level-6 CMZV words."""
    w = tuple(w)
    if not w: return {(): Fr(1)}
    forms = [tuple(hpl_letter_map(mi, sig).items()) for sig in w]
    e0, eh, e1, em, einf = MAPS[mi]
    return forms_path6(forms, e0, eh)

# ---------- cross-ratio atoms ----------
def cr(x, p0, p1, pinf):
    def sub(a,b):
        if a is INF and b is INF: return Fr(0)
        if a is INF: return INF
        if b is INF: return INF
        return a-b
    def div(a,b):
        if a is INF and b is INF: return Fr(1)
        if a is INF: return INF
        if b is INF: return Fr(0)
        return a/b
    f1 = div(sub(x,p0), sub(x,pinf)); f2 = div(sub(p1,pinf), sub(p1,p0))
    if f1 is INF or f2 is INF: return INF
    return f1*f2

FRPTS = (Fr(0), Fr(1,2), Fr(1), Fr(-1), INF)
FRNAM = {Fr(0):'0', Fr(1,2):'1/2', Fr(1):'1', Fr(-1):'-1', INF:'inf'}
FRVAL = {'0':Fr(0), '1/2':Fr(1,2), '1':Fr(1), '-1':Fr(-1), 'inf':INF}

@cache
def cross_table():
    tab = {}
    for p0 in FRPTS:
        for p1 in FRPTS:
            if p1 == p0: continue
            for pinf in FRPTS:
                if pinf in (p0, p1): continue
                for e in FRPTS:
                    if e in (p0, p1, pinf): continue
                    z = cr(e, p0, p1, pinf)
                    if z is not INF:
                        tab.setdefault(z, []).append((e, p0, p1, pinf))
    return tab

@cache
def li_atom(n, z, mi=0, rep=0):
    """Li_n(z), z a cross-ratio of {0,1/2,1,-1,inf}, as level-6 CMZV words under map mi.
    Li_n(z) = -G(0^{n-1},1; z) = -int_{p0}^{e} [dl(x-p0)-dl(x-pinf)]^{n-1} [dl(x-pinf)-dl(x-p1)].
    Pullback: dl(x-q) |-> dlog(u-zeta^{e_q}) - dlog(u-zeta^{e_inf}); differences cancel e_inf."""
    z = Fr(z)
    e, p0, p1, pinf = cross_table()[z][rep]
    e0, eh, e1, em, einf = MAPS[mi]
    def dl(q):
        if q is INF:  # dlog(x - inf) = 0-form
            return {}
        return dlog_pullback(EXP[q], einf)
    EXP = {Fr(0): e0, Fr(1,2): eh, Fr(1): e1, Fr(-1): em, INF: einf}
    f0 = add((dlog_pullback(EXP[p0], einf), 1), (dlog_pullback(EXP[pinf], einf), -1)) if pinf is not INF else dlog_pullback(EXP[p0], einf)
    f1 = add((dlog_pullback(EXP[pinf], einf), 1), (dlog_pullback(EXP[p1], einf), -1)) if pinf is not INF else {w: -c for w, c in dlog_pullback(EXP[p1], einf).items()}
    # dlog(x - inf) is the zero form; handle p0/p1 == INF similarly (never for our reps if chosen well)
    forms = [tuple(f0.items())] * (n - 1) + [tuple(f1.items())]
    cexp = EXP[p0] if p0 is not INF else einf   # dlog(x-p0): p0=inf would be zero form - excluded
    eexp = EXP[e] if e is not INF else einf
    return {w: -c for w, c in forms_path6(forms, cexp, eexp).items()}

# ---------- elementary atoms ----------
def atom_log2():  return {(3,): Fr(-1)}              # (3) = -log(1-zeta^3) = -log 2
def atom_log3():  return {(2,): Fr(-1), (4,): Fr(-1)}  # (2)+(4) = -log 3
def atom_log32(): return add((atom_log3(), 1), (atom_log2(), -1))
def atom_ipi():   return {(1,): Fr(3)}               # (1) = i pi/3
def atom_zeta2(): return {(-1, 0): Fr(1)}
def atom_zeta3(): return {(-1, -1, 0): Fr(1)}
def atom_zeta4(): return {(-1, -1, -1, 0): Fr(1)}
def atom_pi2():   return mul(atom_ipi(), atom_ipi())     # (i pi)^2 = -pi^2... careful: value is (i pi)^2
def atom_pi4():
    p = atom_ipi()
    return mul(mul(p, p), mul(p, p))                      # (i pi)^4 = pi^4

def smul(a, vec):
    return {w: a * c for w, c in vec.items()}

def vadd(*vs):
    out = {}
    for v in vs: out = add((out, 1), (v, 1))
    return out
'''
open('hexc2_common.py','w').write(code)
print('written', len(code))
# quick smoke test
import sys
sys.dont_write_bytecode = True
sys.path.insert(0,'.')
from hexc2_common import hpl_to_cmzv, atom_log2
print('H_1(1/2):', hpl_to_cmzv((1,)))
print('H_-1(1/2):', hpl_to_cmzv((-1,)))
print('H_0(1/2):', hpl_to_cmzv((0,)))
print('log2:', atom_log2())
