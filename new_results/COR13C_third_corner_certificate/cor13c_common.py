"""COR13C shared machinery: exact conversion of the c0=1/3-corner identities (COR7) to
level-6 CMZV words, via the 12 Mobius maps {0,1/3,1,-1,inf} -> mu_6 (step_21/mobius_maps_13.json).

Task map = MAPS[0] with exponents (f(0),f(1/3),f(1),f(-1),f(inf)) = (0,1,2,4,3):
  0->1, 1/3->zeta, 1->zeta^2, -1->zeta^4, inf->zeta^3=-1   (zeta = e^{i pi/3}).
HPL letter pullbacks under g = f^{-1} (dlog(x-q) |-> dlog(u-zeta^{e_q}) - dlog(u-zeta^{e_inf})
= -letter(-e_q) + letter(-e_inf)):
  0  |-> (3) - (0)      1  |-> (4) - (3)      -1  |-> (3) - (2)     [task map]
Every HPL word H_w(1/3) (first letter outermost, integral 0 -> 1/3) becomes an iterated
integral from zeta^{e0} to zeta^{e13}, converted to level-6 CMZV words (0 -> 1) by shuffle-
regularized path composition through 0 (HEXC2 path6, first letter outermost).
Li atoms Li_n(z) for z in {1/3,-1/3,2/3,1/2,-1/2,-2}: all cross-ratios of the configuration
(step_19 table); Li4(1/9) enters only through the distribution identity
Li4(1/9) = 8 (Li4(1/3) + Li4(-1/3)) (absolutely convergent series rearrangement).
"""
import sys, json, re
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
from functools import cache
from itertools import product
BASE = "<project>/research/08_theory"
sys.path.insert(0, BASE + "/HEXC1_hexagon_certificate")
from cmzv6 import sh, mul, reg, regularize, add, words, N

INF = None
_HERE = BASE + "/COR13C_third_corner_certificate"
MAPS = [tuple(m['exps']) for m in json.load(open(_HERE + "/mobius_maps_13.json"))]
# each tuple: (f(0), f(1/3), f(1), f(-1), f(inf)) as exponents of zeta_6
QPTS = ['0', '1/3', '1', '-1', 'inf']

def dlog_pullback(qexp, einf):
    """g*dlog(x-q) as letter dict: -letter(-qexp) + letter(-einf)."""
    if qexp == einf:
        return {}
    return {((-qexp) % 6,): Fr(-1), ((-einf) % 6,): Fr(1)}

def hpl_letter_map(mi, sig):
    """Pullback of HPL letter sig in {0,1,-1} under map mi (task map = 0)."""
    e0, eth, e1, em, einf = MAPS[mi]
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
    n = len(w)
    for k in range(n + 1):
        left = regularize({shift(w[:k], e): Fr(1)})                    # I(w[:k]; 0 -> zeta^e)
        right = regularize({shift(tuple(reversed(w[k:])), c): Fr(1)})  # I(rev(w[k:]); 0 -> zeta^c)
        out = add((out, 1), (mul(left, right), Fr((-1) ** (n - k))))
    return regularize(out)

def forms_path6(forms, c, e):
    """iterated integral of a tuple of letter-dicts (form combos) from zeta^c to zeta^e."""
    out = {}
    for choices in product(*forms):
        w = tuple(x for cw, a in choices for x in cw)
        coef = Fr(1)
        for cw, a in choices: coef *= a
        out = add((out, 1), (path6(w, c, e), coef))
    return out

@cache
def hpl_to_cmzv(w, mi=0):
    """H_w(1/3) (HPL letters 0,+/-1, first letter outermost) as level-6 CMZV words."""
    w = tuple(w)
    if not w: return {(): Fr(1)}
    forms = [tuple(hpl_letter_map(mi, sig).items()) for sig in w]
    e0, eth, e1, em, einf = MAPS[mi]
    return forms_path6(forms, e0, eth)

# ---------- cross-ratio atoms ----------
def cr(x, p0, p1, pinf):
    def sub(a, b):
        if a is INF and b is INF: return Fr(0)
        if a is INF: return INF
        if b is INF: return INF
        return a - b
    def div(a, b):
        if a is INF and b is INF: return Fr(1)
        if a is INF: return INF
        if b is INF: return Fr(0)
        return a / b
    f1 = div(sub(x, p0), sub(x, pinf)); f2 = div(sub(p1, pinf), sub(p1, p0))
    if f1 is INF or f2 is INF: return INF
    return f1 * f2

FRPTS = (Fr(0), Fr(1, 3), Fr(1), Fr(-1), INF)

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
    """Li_n(z), z a cross-ratio of {0,1/3,1,-1,inf}, as level-6 CMZV words under map mi.
    Li_n(z) = -G(0^{n-1},1; z) = -int_{p0}^{e} [dl(x-p0)-dl(x-pinf)]^{n-1} [dl(x-pinf)-dl(x-p1)],
    dl(x-q) |-> dlog(u-zeta^{e_q}) - dlog(u-zeta^{e_inf}), dl(x-inf) = zero form."""
    z = Fr(z)
    e, p0, p1, pinf = cross_table()[z][rep]
    e0, eth, e1, em, einf = MAPS[mi]
    EXP = {Fr(0): e0, Fr(1, 3): eth, Fr(1): e1, Fr(-1): em, INF: einf}
    def dl(q):
        return dlog_pullback(EXP[q], einf)   # EXP[INF]=einf gives {} = zero form
    f0 = add((dl(p0), 1), (dl(pinf), -1))
    f1 = add((dl(p1), 1), (dl(pinf), -1))
    forms = [tuple(f0.items())] * (n - 1) + [tuple(f1.items())]
    cexp = EXP[p0] if p0 is not INF else einf
    eexp = EXP[e] if e is not INF else einf
    return {w: -c for w, c in forms_path6(forms, cexp, eexp).items()}

# ---------- elementary atoms (level 6, map-independent; same as HEXC2) ----------
def atom_log2():  return {(3,): Fr(-1)}                  # -(3) = log(1-zeta^3)?? value: (3) = -log(1-zeta^3) = -log 2, so -(3) = log 2
def atom_log3():  return {(2,): Fr(-1), (4,): Fr(-1)}    # -(2)-(4) = log 3
def atom_log32(): return add((atom_log3(), 1), (atom_log2(), -1))
def atom_ipi():   return {(1,): Fr(3)}                   # 3*(1) = i pi
def atom_zeta2(): return {(-1, 0): Fr(1)}
def atom_zeta3(): return {(-1, -1, 0): Fr(1)}
def atom_pi2():   return mul(atom_ipi(), atom_ipi())     # (i pi)^2 = -pi^2
def atom_pi4():
    p = atom_ipi()
    return mul(mul(p, p), mul(p, p))                     # (i pi)^4 = pi^4

def smul(a, vec): return {w: a * c for w, c in vec.items()}
def vadd(*vs):
    out = {}
    for v in vs: out = add((out, 1), (v, 1))
    return out
def prod(*vecs):
    out = {(): Fr(1)}
    for v in vecs: out = mul(out, v)
    return out
def Li(n, z): return li_atom(n, Fr(z))

# ---------- the c0 = 1/3 corner data (G1 fits; COR7 closed forms) ----------
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
C7 = BASE + "/COR7_classical_corners"

def _fit_terms(name):
    D = json.load(open(G1 + "/" + name))
    out = []
    for key, val in D:
        if key == "1": e, w = 0, ()
        elif key == "c": e, w = 1, ()
        else:
            m = re.match(r"c\^(-?\d+) H\[(.*)\]", key)
            e = int(m.group(1))
            w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
        out.append((e, w, Fr(val)))
    return out

def fit_vec(name, weight):
    """sum cf * (1/3)^e * H_w(1/3) over terms of given HPL weight in a G1 fit file."""
    out = {}
    for e, w, v in _fit_terms(name):
        if len(w) != weight: continue
        vv = {(): Fr(1)} if not w else hpl_to_cmzv(w)
        out = add((out, 1), (vv, v * Fr(1, 3) ** e))
    return out

def R3(weight): return fit_vec("hpl_fit_k3_(0,0,0).json", weight)
def Psi(weight): return fit_vec("hpl_fit_k3_(1,0,0).json", weight)

PI2 = smul(Fr(-1), atom_pi2())   # pi^2
PI4 = atom_pi4()                 # pi^4
L2, L3 = atom_log2(), atom_log3()
Z2, Z3 = atom_zeta2(), atom_zeta3()

def _atom_vec(name):
    """one COR7 basis monomial -> CMZV vector."""
    if name == "zeta3": return Z3
    if name == "pi4": return PI4
    m = re.fullmatch(r"Li(\d)\((-?\d/\d|-?\d+)\)", name)
    if m:
        n = int(m.group(1)); zz = Fr(m.group(2))
        if zz == Fr(1, 9):   # distribution input: Li4(1/9) = 8 (Li4(1/3) + Li4(-1/3))
            return smul(Fr(8), vadd(li_atom(4, Fr(1, 3)), li_atom(4, Fr(-1, 3))))
        return li_atom(n, zz)
    facs = []
    for f in name.split("*"):
        mm = re.fullmatch(r"(Li\d\(-?\d/\d\)|zeta3|pi2|log2|log3)(?:\^(\d+))?", f)
        assert mm, f
        base, pw = mm.group(1), int(mm.group(2) or 1)
        if base == "pi2": v = PI2
        elif base == "log2": v = L2
        elif base == "log3": v = L3
        elif base == "zeta3": v = Z3
        else:
            n2 = int(base[2]); zz = Fr(base[4:-1]); v = li_atom(n2, zz)
        facs += [v] * pw
    return prod(*facs)

def closed_vec(json_name):
    """COR7 closed form: relation a0*W + sum coeffs*basis = 0  =>  W = -sum/a0."""
    R = json.load(open(C7 + "/" + json_name))
    out = {}
    for cf, bn in zip(R["coeffs"], R["basis"]):
        out = add((out, 1), (_atom_vec(bn), Fr(-cf, R["a0"])))
    return out

def W2_closed_vec():
    # W2 = -Li2(1/3) - Li2(-1/3) + pi^2/12  (COR7, exact)
    return vadd(smul(Fr(-1), Li(2, '1/3')), smul(Fr(-1), Li(2, '-1/3')), smul(Fr(1, 12), PI2))

def W3_closed_vec(): return closed_vec("c13_W3_relation.json")
def W4_closed_vec(): return closed_vec("c13_F4_relation.json")

# ---------- targets: LHS - RHS as exact CMZV vectors ----------
def target2():
    # R3^(2)(1/3) + zeta(2)*Psi^(0) - W2closed ; Psi^(0) = 1/2
    return vadd(R3(2), smul(Fr(1, 2), Z2), smul(Fr(-1), W2_closed_vec()))

def target3():
    # R3^(3)(1/3) + zeta(2) Psi^(1)(1/3) - zeta(3)/2 - W3closed
    return vadd(R3(3), mul(Z2, Psi(1)), smul(Fr(-1, 2), Z3), smul(Fr(-1), W3_closed_vec()))

def target4():
    # R3^(4)(1/3) + zeta(2) Psi^(2)(1/3) + zeta(3)*((1-c^2)/(2c)) atanh(c) - W4closed,
    # at c = 1/3: (1-c^2)/(2c) = 4/3, atanh(1/3) = (log 2)/2  =>  (2/3) zeta(3) log 2
    return vadd(R3(4), mul(Z2, Psi(2)), smul(Fr(2, 3), mul(Z3, L2)), smul(Fr(-1), W4_closed_vec()))
