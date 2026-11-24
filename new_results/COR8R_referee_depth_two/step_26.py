import cor8r_sym as S
import sympy as sp
from fractions import Fraction
cc = S.cc

def top_terms(word):
    """Recompute the recursion's top-level terms: list of (subword, letter_value, letter_vec)."""
    word = tuple(sp.sympify(w) for w in word)
    n = len(word)
    a = (sp.Integer(1),)+word+(sp.Integer(0),)
    out = []
    for i in range(1,n+1):
        lv = sp.cancel((a[i-1]-a[i])/(a[i]-a[i+1]))
        L = S.factor_Qc(lv)
        if S.is_constant_letter(L): continue
        sub = word[:i-1]+word[i:]
        out.append((sub, lv, L))
    return out

x, y = sp.Rational(1,2), 1-cc
word = (sp.Integer(0), 1/y, sp.Integer(0), 1/(x*y))
sym = S.symG(word)
print("S(Li_{2,2}(1/2,1-c)) convA, #terms:", len(sym))
NAMES = ['-1','2','c','1-c','1+c']
def showL(L):
    parts = []
    for n,e in zip(NAMES,L):
        if e: parts.append(f"{n}^{e}" if e!=1 else n)
    return "["+" ".join(parts)+"]" if parts else "[1]"
for w,cf in sorted(sym.items(), key=lambda t:(len(t[0]),t[0])):
    print("  ", " ".join(showL(L) for L in w), " coeff:", cf)
print("\ntop terms:")
for sub,lv,L in top_terms(word):
    print("  G", sub, " dlog(", lv, ")")
