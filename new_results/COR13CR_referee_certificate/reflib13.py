"""reflib13.py -- COR13CR referee's c0=1/3-specific independent layer.

Reuses HEXC2R's reflib.py for the GENERIC machinery (word algebra, relation generators,
truncated-sum CMZV/HPL evaluator, seg_value/forms_value path composition). Everything
specific to c0 = 1/3 is re-derived here independently:
  * MY13 enumeration of the 12 Mobius maps {0, 1/3, 1, -1, inf} -> mu_6 (exact over Q(zeta_6),
    my own pair arithmetic below);
  * HPL letter pullbacks at the 1/3 endpoint;
  * cross-ratio frames from the five points {0, 1/3, 1, -1, inf};
  * Li-atom vectors, elementary atoms, targets T2/T3/T4 from the G1 fit files and the
    COR7 closed-form JSONs.
"""
import sys, json, re
sys.dont_write_bytecode = True
BASE = "<project>/research/08_theory"
sys.path.insert(0, BASE + "/HEXC2R_referee_certificate")
sys.path.insert(0, BASE + "/HEXC1_hexagon_certificate")
from fractions import Fraction as Fr
from itertools import product as _product
import reflib as R

INF = None
PTS13 = [Fr(0), Fr(1, 3), Fr(1), Fr(-1), INF]   # the c0=1/3 configuration

# ---------------- my own exact Q(zeta_6) arithmetic: (a,b) = a + b zeta, zeta^2 = zeta - 1
def qadd(x, y): return (x[0] + y[0], x[1] + y[1])
def qneg(x): return (-x[0], -x[1])
def qsub(x, y): return (x[0] - y[0], x[1] - y[1])
def qmul(x, y):
    a, b = x; c, d = y
    return (a*c - b*d, a*d + b*c + b*d)
def qinv(x):
    a, b = x
    n = a*a + a*b + b*b
    assert n != 0
    return (Fr(a + b, n), Fr(-b, n))
def qdiv(x, y): return qmul(x, qinv(y))
def qpow(x, n):
    out = (Fr(1), Fr(0))
    for _ in range(n): out = qmul(out, x)
    return out
ZETA6 = (Fr(0), Fr(1))
MU6 = [qpow(ZETA6, k) for k in range(6)]

def qmobius(abcd, x):
    A, B, C, D = abcd
    if x is INF:
        return qdiv(A, C) if C != (Fr(0), Fr(0)) else INF
    xq = (Fr(x), Fr(0))
    n = qadd(qmul(A, xq), B); d = qadd(qmul(C, xq), D)
    return INF if d == (Fr(0), Fr(0)) else qdiv(n, d)

def mobius_from_3(f0, f1, finf):
    """map with f(0)=f0, f(1)=f1, f(inf)=finf (all finite, distinct)."""
    C = (Fr(1), Fr(0)); A = finf
    D = qdiv(qsub(finf, f1), qsub(f1, f0))
    B = qmul(f0, D)
    return (A, B, C, D)

def enum_maps13():
    """all Mobius maps sending {0,1/3,1,-1,inf} injectively into mu_6; returns sorted
    exponent tuples (f(0), f(1/3), f(1), f(-1), f(inf))."""
    from itertools import permutations
    found = set()
    for f0, f1, finf in permutations(MU6, 3):
        M = mobius_from_3(f0, f1, finf)
        ims = [qmobius(M, p) for p in PTS13]
        if all((im is not INF and im in MU6) for im in ims) and len(set(ims)) == 5:
            found.add(tuple(MU6.index(im) for im in ims))
    return sorted(found)

# The 12 maps (exponents of f(0), f(1/3), f(1), f(-1), f(inf)); filled by enum_maps13(),
# cached here after first exact enumeration (see step that asserts equality).
MAPS13 = None
def maps13():
    global MAPS13
    if MAPS13 is None:
        MAPS13 = enum_maps13()
    return MAPS13

# ---------------- conversion at c0 = 1/3 (my derivation, mirroring reflib's structure) --
def _dlogU(k):
    """dlog(u - zeta^k) as a letter vector: -letter((-k) mod 6)."""
    return {((-k) % 6,): Fr(-1)}

def hpl_form_pullback13(mi, sig):
    """pullback of HPL form omega_sig under g = f^{-1}: g*dlog(x-q) = dlog(u-f(q)) - dlog(u-f(inf)).
    letters: 0 = dlog x ; 1 = -dlog(x-1) ; -1 = dlog(x+1)."""
    e0, eh, e1, em, einf = maps13()[mi]
    if sig == 0:
        return R.vadd((_dlogU(e0), 1), (_dlogU(einf), -1))
    if sig == 1:
        return R.vadd((_dlogU(e1), -1), (_dlogU(einf), 1))
    if sig == -1:
        return R.vadd((_dlogU(em), 1), (_dlogU(einf), -1))
    raise ValueError(sig)

def my13_hpl_to_cmzv(w, mi=0):
    """H_w(1/3) (letters 0,±1, first letter outermost) as level-6 CMZV vector."""
    w = tuple(w)
    if not w: return {(): Fr(1)}
    formvecs = [hpl_form_pullback13(mi, sig) for sig in w]
    e0, eh, e1, em, einf = maps13()[mi]
    return R.forms_value(formvecs, e0, eh)

# ---------------- cross-ratio frames from {0, 1/3, 1, -1, inf} ----------------
def _crq(e, p0, p1, pinf):
    def sub(a, b):
        if a is None or b is None: return None
        return a - b
    def div(a, b):
        if a is None and b is None: return Fr(1)
        if a is None: return None
        if b is None: return Fr(0)
        return a / b
    f1 = div(sub(e, p0), sub(e, pinf)); f2 = div(sub(p1, pinf), sub(p1, p0))
    if f1 is None or f2 is None: return None
    return f1 * f2

def cross_reps13(z):
    """all frames (e, p0, p1, pinf) of distinct points of {0,1/3,1,-1,inf} with cr = z."""
    from itertools import permutations as _perm
    out = []
    for e, p0, p1, pinf in _perm(PTS13, 4):
        if _crq(e, p0, p1, pinf) == Fr(z):
            out.append((e, p0, p1, pinf))
    return out

def cross_values13():
    """all finite cross-ratio values of the configuration, with frame counts."""
    from itertools import permutations as _perm
    tab = {}
    for e, p0, p1, pinf in _perm(PTS13, 4):
        z = _crq(e, p0, p1, pinf)
        if z is not None:
            tab.setdefault(z, []).append((e, p0, p1, pinf))
    return tab

def my13_li_atom(n, z, mi=0, rep=0, frame=None):
    """Li_n(z) as a level-6 CMZV vector under map mi (my derivation):
    Li_n(z) = -G(0^{n-1},1; z), z = (e,p0;p1,pinf); outer forms dlog(x-p0)-dlog(x-pinf),
    inner form dlog(x-p1)-dlog(x-pinf) pulled back; path f(p0) -> f(e)."""
    if frame is not None:
        e, p0, p1, pinf = frame
        assert _crq(e, p0, p1, pinf) == Fr(z), (frame, z)
    else:
        e, p0, p1, pinf = cross_reps13(z)[rep]
    e0, eh, e1, em, einf = maps13()[mi]
    EXP = {Fr(0): e0, Fr(1, 3): eh, Fr(1): e1, Fr(-1): em, INF: einf}
    def dlx(q):
        if q is INF:
            return {}
        return R.vadd((_dlogU(EXP[q]), 1), (_dlogU(EXP[INF]), -1))
    f0 = R.vadd((dlx(p0), 1), (dlx(pinf), -1))
    f1 = R.vadd((dlx(p1), 1), (dlx(pinf), -1))
    formvecs = [f0] * (n - 1) + [f1]
    cexp = EXP[p0] if p0 is not INF else einf
    eexp = EXP[e] if e is not INF else einf
    return {w: -c for w, c in R.forms_value(formvecs, cexp, eexp).items()}

# ---------------- elementary atoms (map-independent; re-stated) ----------------
def a_log2():  return {(3,): Fr(-1)}
def a_log3():  return {(2,): Fr(-1), (4,): Fr(-1)}
def a_ipi():   return {(1,): Fr(3)}
def a_zeta2(): return {(-1, 0): Fr(1)}
def a_zeta3(): return {(-1, -1, 0): Fr(1)}
def a_pi2():   return R.mul(a_ipi(), a_ipi())
def a_pi4():   return R.mul(a_pi2(), a_pi2())

def smul(a, vec): return {w: a * c for w, c in vec.items()}
def prod(*vecs):
    out = {(): Fr(1)}
    for v in vecs: out = R.mul(out, v)
    return out

# ---------------- targets from G1 fits + COR7 closed forms ----------------
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
C7 = BASE + "/COR7_classical_corners"

def fit_terms(name):
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

def fit_vec13(name, weight, mi=0):
    """sum cf * (1/3)^e * H_w(1/3) over terms of given HPL weight, as CMZV vector."""
    out = {}
    for e, w, v in fit_terms(name):
        if len(w) != weight: continue
        vv = {(): Fr(1)} if not w else my13_hpl_to_cmzv(w, mi)
        out = R.vadd((out, 1), (vv, v * Fr(1, 3) ** e))
    return out

def R3_13(weight, mi=0): return fit_vec13("hpl_fit_k3_(0,0,0).json", weight, mi)
def Psi_13(weight, mi=0): return fit_vec13("hpl_fit_k3_(1,0,0).json", weight, mi)

def li_atom13(n, z, mi=0, rep=0):
    z = Fr(z)
    if z == Fr(1, 9):   # input theorem: Li4(1/9) = 8(Li4(1/3) + Li4(-1/3))
        assert n == 4
        return smul(Fr(8), R.vadd((my13_li_atom(4, Fr(1, 3), mi, rep), 1),
                                  (my13_li_atom(4, Fr(-1, 3), mi, rep), 1)))
    return my13_li_atom(n, z, mi, rep)

def atom_vec13(name, mi=0):
    """one COR7 basis monomial -> CMZV vector (my construction)."""
    if name == "zeta3": return a_zeta3()
    if name == "pi4": return a_pi4()
    m = re.fullmatch(r"Li(\d)\((-?\d/\d|-?\d+)\)", name)
    if m:
        return li_atom13(int(m.group(1)), Fr(m.group(2)), mi)
    facs = []
    for f in name.split("*"):
        mm = re.fullmatch(r"(Li\d\(-?\d/\d\)|zeta3|pi2|log2|log3)(?:\^(\d+))?", f)
        assert mm, f
        base, pw = mm.group(1), int(mm.group(2) or 1)
        if base == "pi2": v = smul(Fr(-1), a_pi2())
        elif base == "log2": v = a_log2()
        elif base == "log3": v = a_log3()
        elif base == "zeta3": v = a_zeta3()
        else: v = li_atom13(int(base[2]), Fr(base[4:-1]), mi)
        facs += [v] * pw
    return prod(*facs)

def closed_vec13(json_name, mi=0):
    R_ = json.load(open(C7 + "/" + json_name))
    out = {}
    for cf, bn in zip(R_["coeffs"], R_["basis"]):
        out = R.vadd((out, 1), (atom_vec13(bn, mi), Fr(-cf, R_["a0"])))
    return out

def W2_closed13(mi=0):
    return R.vadd((li_atom13(2, "1/3", mi), -1), (li_atom13(2, "-1/3", mi), -1),
                  (smul(Fr(-1), a_pi2()), Fr(1, 12)))
def W3_closed13(mi=0): return closed_vec13("c13_W3_relation.json", mi)
def W4_closed13(mi=0): return closed_vec13("c13_F4_relation.json", mi)

def my_target2(mi=0):
    # R3^(2)(1/3) + zeta(2)/2 - W2closed
    return R.vadd((R3_13(2, mi), 1), (a_zeta2(), Fr(1, 2)), (W2_closed13(mi), -1))
def my_target3(mi=0):
    # R3^(3) + zeta(2) Psi^(1) - zeta(3)/2 - W3closed
    return R.vadd((R3_13(3, mi), 1), (R.mul(a_zeta2(), Psi_13(1, mi)), 1),
                  (a_zeta3(), Fr(-1, 2)), (W3_closed13(mi), -1))
def my_target4(mi=0):
    # R3^(4) + zeta(2) Psi^(2) + (2/3) zeta(3) log 2 - W4closed
    return R.vadd((R3_13(4, mi), 1), (R.mul(a_zeta2(), Psi_13(2, mi)), 1),
                  (R.mul(a_zeta3(), a_log2()), Fr(2, 3)), (W4_closed13(mi), -1))
