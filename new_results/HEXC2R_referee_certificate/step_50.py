
reflib_conv = '''

# ================= independent conversion machinery (Task 2) =================
# My own derivation (see step_06 header). Notation: letter (j) = zeta^j du/(1-zeta^j u),
# dlog(u - zeta^k) = -letter((-k) mod 6). MAPS[mi] = (e0, eh, e1, em, einf) = exponents of
# f(0), f(1/2), f(1), f(-1), f(inf) for the 12 maps.
from itertools import product as _product

MAPS6 = [
    (0,1,2,5,4), (0,5,4,1,2), (1,2,3,0,5), (1,0,5,2,3),
    (2,1,0,3,4), (2,3,4,1,0), (3,2,1,4,5), (3,4,5,2,1),
    (4,5,0,3,2), (4,3,2,5,0), (5,0,1,4,3), (5,4,3,0,1),
]

def _dlogU(k):
    """dlog(u - zeta^k) as a letter vector: -letter((-k) mod 6)."""
    return {((-k) % 6,): Fr(-1)}

def hpl_form_pullback(mi, sig):
    """pullback of HPL form omega_sig (0: dlog x; 1: -dlog(x-1); -1: dlog(x+1)) under g=f^{-1},
    as a level-6 letter vector. g*dlog(x-q) = dlog(u-f(q)) - dlog(u-f(inf))."""
    e0, eh, e1, em, einf = MAPS6[mi]
    if sig == 0:
        return vadd((_dlogU(e0), 1), (_dlogU(einf), -1))
    if sig == 1:
        return vadd((_dlogU(e1), -1), (_dlogU(einf), 1))
    if sig == -1:
        return vadd((_dlogU(em), 1), (_dlogU(einf), -1))
    raise ValueError(sig)

def _shift6(w, s):
    return tuple((-1 if j == -1 else (j + s) % 6) for j in w)

def seg_value(w, c, e):
    """I(w; zeta^c -> zeta^e), first letter outermost, shuffle-regularized (T=0).
    My path composition through 0 (the outer letters w[:k] go with the e-side):
        I(w; a->b) = sum_k I(w[:k]; 0->b) * I(w[k:]; a->0),
        I(v; a->0) = (-1)^|v| I(reversed(v); 0->a),
    products of regularized values = shuffle products of the vectors."""
    w = tuple(w)
    if c == e:
        return {(): Fr(1)} if not w else {}
    out = {}
    n = len(w)
    for k in range(n + 1):
        left = regularize({_shift6(w[:k], e): Fr(1)})          # I(w[:k]; 0 -> zeta^e)
        right = regularize({_shift6(tuple(reversed(w[k:])), c): Fr(1)})  # I(rev(w[k:]); 0 -> zeta^c)
        out = vadd((out, 1), (mul(left, right), Fr((-1) ** (n - k))))
    return regularize(out)

def forms_value(formvecs, c, e):
    """iterated integral of a tuple of letter-vectors (form combinations) from zeta^c to zeta^e."""
    out = {}
    for choices in _product(*[tuple(fv.items()) for fv in formvecs]):
        w = tuple(x for cw, a in choices for x in cw)
        coef = Fr(1)
        for cw, a in choices: coef *= a
        out = vadd((out, 1), (seg_value(w, c, e), coef))
    return out

def my_hpl_to_cmzv(w, mi=0):
    """H_w(1/2) (letters 0,±1, first letter outermost) as a level-6 CMZV vector."""
    w = tuple(w)
    if not w: return {(): Fr(1)}
    formvecs = [hpl_form_pullback(mi, sig) for sig in w]
    e0, eh, e1, em, einf = MAPS6[mi]
    return forms_value(formvecs, e0, eh)

# ---- cross-ratio frames for Li atoms (my own enumeration) ----
def _crq(e, p0, p1, pinf):
    """rational cross-ratio (e,p0;p1,pinf); points are Fractions or None (=inf)."""
    def sub(a, b):
        if a is None: return None
        if b is None: return None
        return a - b
    def div(a, b):
        if a is None and b is None: return Fr(1)
        if a is None: return None
        if b is None: return Fr(0)
        return a / b
    f1 = div(sub(e, p0), sub(e, pinf)); f2 = div(sub(p1, pinf), sub(p1, p0))
    if f1 is None or f2 is None: return None
    return f1 * f2

def cross_reps(z):
    """all (e, p0, p1, pinf) with cr = z, from {0, 1/2, 1, -1, inf}."""
    from itertools import permutations as _perm
    PTS = [Fr(0), Fr(1, 2), Fr(1), Fr(-1), None]
    out = []
    for e, p0, p1, pinf in _perm(PTS, 4):
        if _crq(e, p0, p1, pinf) == Fr(z):
            out.append((e, p0, p1, pinf))
    return out

def my_li_atom(n, z, mi=0, rep=0):
    """Li_n(z) as a level-6 CMZV vector under map mi, using the cross-ratio frame
    z = (e,p0;p1,pinf):  Li_n(z) = -G(0^{n-1},1; z) pulled back to u:
    outer forms dlog(t) -> dlog(u-f(p0)) - dlog(u-f(pinf));  inner form dlog(t/(t-1))...:
    dt/(1-t) -> dlog(u-f(pinf)) - dlog(u-f(p1)).  Path u: f(p0) -> f(e)."""
    e, p0, p1, pinf = cross_reps(z)[rep]
    e0, eh, e1, em, einf = MAPS6[mi]
    EXP = {Fr(0): e0, Fr(1, 2): eh, Fr(1): e1, Fr(-1): em, None: einf}
    def dlx(q):
        # dlog(x - q) pulls back to dlog(u - f(q)) - dlog(u - f(inf)); q = inf gives 0-form
        if q is None:
            return {}
        return vadd((_dlogU(EXP[q]), 1), (_dlogU(EXP[None]), -1))
    f0 = vadd((dlx(p0), 1), (dlx(pinf), -1))
    f1 = vadd((dlx(pinf), 1), (dlx(p1), -1))
    formvecs = [f0] * (n - 1) + [f1]
    cexp = EXP[p0] if p0 is not None else einf
    eexp = EXP[e] if e is not None else einf
    return {w: -c for w, c in forms_value(formvecs, cexp, eexp).items()}

# ---- elementary atoms (my own derivations) ----
def my_atom_log2():  return {(3,): Fr(-1)}              # value (3,) = -log(1-zeta^3) = -log 2
def my_atom_log3():  return {(2,): Fr(-1), (4,): Fr(-1)}  # (2,)+(4,) = -log((1-z^2)(1-z^4)) = -log 3
def my_atom_log32(): return vadd((my_atom_log3(), 1), (my_atom_log2(), -1))
def my_atom_ipi():   return {(1,): Fr(3)}               # value (1,) = -Log(1-zeta) = i pi/3
def my_atom_zeta2(): return {(-1, 0): Fr(1)}
def my_atom_zeta3(): return {(-1, -1, 0): Fr(1)}
def my_atom_zeta4(): return {(-1, -1, -1, 0): Fr(1)}
def my_atom_pi2():   return mul(my_atom_ipi(), my_atom_ipi())       # (i pi)^2
def my_atom_pi4():   return mul(my_atom_pi2(), my_atom_pi2())       # (i pi)^4
'''
with open("reflib.py", "a") as f:
    f.write(reflib_conv)
print("conversion machinery appended; reflib chars:", len(open("reflib.py").read()))
