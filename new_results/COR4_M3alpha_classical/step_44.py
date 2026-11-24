cor = "<project>/research/08_theory/COR4_M3alpha_classical/"
src = open(cor+"cor4_sym.py").read()
add = '''

# ---------- corrected (2,2) cobracket: target Lambda^2(B_2) ----------
def delta22L2(S, dim):
    """Goncharov (2,2) component into Lambda^2(B_2):
    (l1,l2,l3,l4) -> (l1^l2)^(l3^l4)  (antisymmetrised within AND between factors).
    Returns dict ((ij),(kl)) with (ij)<(kl) -> Fr.  Kills Li4 and all shuffle products."""
    out = {}
    for w, cf in S.items():
        a, b, c, d = w
        for p1, c1 in wedge(a, b, dim).items():
            for p2, c2 in wedge(c, d, dim).items():
                if p1 == p2:
                    continue
                key = (p1, p2) if p1 < p2 else (p2, p1)
                sgn = 1 if p1 < p2 else -1
                out[key] = out.get(key, 0) + cf * c1 * c2 * sgn
    return {k: v for k, v in out.items() if v}

def delta22L2_func(FS, dim):
    """Function-level version: dict wedgepair -> dict e -> Fr."""
    out = {}
    for (e, w), cf in FS.items():
        a, b, c, d = w
        for p1, c1 in wedge(a, b, dim).items():
            for p2, c2 in wedge(c, d, dim).items():
                if p1 == p2:
                    continue
                key = (p1, p2) if p1 < p2 else (p2, p1)
                sgn = 1 if p1 < p2 else -1
                dd = out.setdefault(key, {})
                dd[e] = dd.get(e, 0) + cf * c1 * c2 * sgn
    return {k: {e: v for e, v in dd.items() if v} for k, dd in out.items() if any(dd.values())}

def delta22L2_spec(S, dim):
    """Specialised version; coefficients may be Fr or (Fr,Fr)."""
    out = {}
    for w, cf in S.items():
        a, b, c, d = w
        for p1, c1 in wedge(a, b, dim).items():
            for p2, c2 in wedge(c, d, dim).items():
                if p1 == p2:
                    continue
                key = (p1, p2) if p1 < p2 else (p2, p1)
                sgn = 1 if p1 < p2 else -1
                v = c1 * c2 * sgn
                if isinstance(cf, tuple):
                    cur = out.get(key, (Fr(0), Fr(0)))
                    out[key] = (cur[0] + cf[0]*v, cur[1] + cf[1]*v)
                else:
                    out[key] = out.get(key, 0) + cf * v
    if out and isinstance(next(iter(out.values())), tuple):
        return {k: v for k, v in out.items() if v != (Fr(0), Fr(0))}
    return {k: v for k, v in out.items() if v}

# ---------- (3,1) and (1,3) components via Lie-word pairings ----------
def lie_bracket(words):
    """[[l1,l2],l3] as dict word->coeff from a tuple of 3 letters."""
    l1, l2, l3 = words
    out = {}
    for w, c in [((l1,l2,l3),1),((l2,l1,l3),-1),((l3,l1,l2),-1),((l3,l2,l1),1)]:
        out[w] = out.get(w, 0) + c
    return {w: c for w, c in out.items() if c}

def delta31(S, dim):
    """Pair first three letters with the Lie bracket [[x,y],z], tensor the fourth.
    Returns dict (basis3word, letterindex) -> Fr where basis3word is a word in the
    letter-basis indices after expanding."""
    out = {}
    for w, cf in S.items():
        l1, l2, l3, l4 = w
        for bw, bc in lie_bracket((l1, l2, l3)).items():
            # expand all four letters into basis indices
            for j4 in range(dim):
                if not l4[j4]:
                    continue
                terms = [((), cf * bc * l4[j4])]
                for v in bw:
                    terms = [(p + (j,), c * v[j]) for p, c in terms for j in range(dim) if v[j]]
                for p, c in terms:
                    out[(p, j4)] = out.get((p, j4), 0) + c
    return {k: v for k, v in out.items() if v}

def delta13(S, dim):
    out = {}
    for w, cf in S.items():
        l1, l2, l3, l4 = w
        for bw, bc in lie_bracket((l2, l3, l4)).items():
            for j1 in range(dim):
                if not l1[j1]:
                    continue
                terms = [((), cf * bc * l1[j1])]
                for v in bw:
                    terms = [(p + (j,), c * v[j]) for p, c in terms for j in range(dim) if v[j]]
                for p, c in terms:
                    out[(j1, p)] = out.get((j1, p), 0) + c
    return {k: v for k, v in out.items() if v}
'''
open(cor+"cor4_sym.py","w").write(src + add)
print("module updated")

# quick sanity: corrected delta22 kills Li4, kills Li2^2 shuffle, kills Li22(x,1); F still nonzero
import sys, importlib
sys.path.insert(0, cor)
import cor4_sym
importlib.reload(cor4_sym)
from cor4_sym import *
X = (Fr(1),Fr(0),Fr(0)); M = (Fr(0),Fr(1),Fr(0)); P = (Fr(0),Fr(0),Fr(1))
G = {0: X, 1: lneg(M), -1: P}
print("Li4:", delta22L2(li_sym(4, X, M), 3))
print("Li22(x,1)=H0101:", delta22L2(hpl_sym((0,1,0,1), G), 3))
print("Li31(x,1)=H0011:", delta22L2(hpl_sym((0,0,1,1), G), 3))
print("H010m1(-x) positive control:", delta22L2(hpl_sym((0,1,0,-1), {0:X, 1:lneg(P), -1:M}), 3))
# Li2(b)^2 shuffle
bb = ((Fr(1),Fr(0)), (Fr(0),Fr(1)))
S = {}
for w in shuffles(bb, bb):
    S[w] = Fr(1)
print("shuffle Li2^2:", delta22L2(S, 2))
