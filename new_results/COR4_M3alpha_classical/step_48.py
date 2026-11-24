cor = "<project>/research/08_theory/COR4_M3alpha_classical/"
helper = '''
# ---------- Q(sqrt2) S-unit utilities ----------
def _cmul(a, b): return (a[0]*b[0] + 2*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def _cdiv(a, b):
    p, q = b; n = p*p - 2*q*q
    return _cmul(a, (p/n, -q/n))
_SQ2 = (Fr(0), Fr(1)); _A1 = (Fr(1), Fr(1))

def sunit_vec_sq2(z):
    """z=(p,q)=p+q sqrt2 in Q(sqrt2)^*.  If |z| = 2^{j/2} (1+sqrt2)^m with j,m in Z,
    return (Fr(j,2), Fr(m)) [sign/torsion ignored]; else None. Exact verification."""
    if z == (Fr(0), Fr(0)): return None
    p, q = z
    N = p*p - 2*q*q
    if N == 0: return None
    s = 1 if N > 0 else -1
    j2 = None
    n = N if N > 0 else -N
    # n must be a power of 2 (as a rational): numerator & denominator powers of 2
    num, den = n.numerator, n.denominator
    if num & (num - 1) or den & (den - 1):
        return None
    import math
    j = num.bit_length() - den.bit_length()  # n = 2^j
    # z = 2^{j/2} u : compute u = z / 2^{j/2}
    if j % 2 == 0:
        f2 = (Fr(2) ** (j // 2), Fr(0)) if j >= 0 else (Fr(1, 2 ** (-j // 2)), Fr(0))
    else:
        f2 = _cmul2 = None
        if j > 0:
            f2 = _cmul(_SQ2, (Fr(2) ** ((j - 1) // 2), Fr(0)))
        else:
            f2 = _cdiv(_SQ2, (Fr(2) ** ((-j - 1) // 2), Fr(0)))
    u = _cdiv(z, f2)
    # u must be +- (1+sqrt2)^m
    import math
    uf = float(u[0]) + float(u[1]) * math.sqrt(2)
    if uf == 0: return None
    m = int(round(math.log(abs(uf)) / math.log(1 + math.sqrt(2))))
    uu = (Fr(1), Fr(0))
    for _ in range(abs(m)):
        uu = _cmul(uu, _A1) if m > 0 else _cdiv(uu, _A1)
    if u != uu and u != (-uu[0], -uu[1]):
        return None
    return (Fr(j, 2), Fr(m))
'''
src = open(cor+"cor4_sym.py").read()
open(cor+"cor4_sym.py","w").write(src + helper)
import sys, importlib
sys.path.insert(0, cor)
import cor4_sym
importlib.reload(cor4_sym)
from cor4_sym import sunit_vec_sq2
from fractions import Fraction as Fr
h = Fr(1,2)
print(sunit_vec_sq2((h, Fr(0))), "expect (-1/2,0)")
print(sunit_vec_sq2((h, Fr(-1))), "expect (-1/2,-1)")
print(sunit_vec_sq2((h, Fr(1))), "expect (-1/2,1)")
print(sunit_vec_sq2((Fr(0), Fr(-1))), "expect (0,-1) for b")
print(sunit_vec_sq2((Fr(3), Fr(0))), "expect None")
print(sunit_vec_sq2((Fr(1), Fr(-2))), "expect None for 1-2sqrt2")
