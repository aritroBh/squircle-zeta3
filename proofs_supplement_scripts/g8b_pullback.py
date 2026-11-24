"""g8b_pullback.py -- exact checks for Corollary g8b:cor:pullback (arguments rational in e^{i alpha/2}).

The map is phi(t) = (t + 1/t)/2 = cos(alpha/2) for t = e^{i alpha/2}.  Everything is exact (sympy, rational numbers, the field Q(i)).
Data files needed (folder scripts/): g8_terms_k3_part1.json  (the 74 terms of R_3).

Part 1.  The letters.  Zeros and poles in C of phi, phi-1, phi+1; the classes [phi-a] in V_Sigma_t, Sigma_t = {0, i, -i, 1, -1}.
Part 2.  The linear map phi_*: V_{0,1,-1} -> V_{Sigma_t}, e_a -> [phi - a]: its rank (injective iff rank 3), and the rank of
         Lambda^2 Lambda^2 phi_* (3 x 45 matrix).
Part 3.  The pulled-back (2,2) part.  The 24 tensors of S(Omega) (from the json, as in g8b_check.py) are pulled back letter by letter,
         with the coefficients c -> phi(t), and delta_{2,2} is applied directly in T^4(V_Sigma_t).  The result is compared with
         the formula  (1 - 1/phi) Lambda^2Lambda^2 phi_*(p^r) - (1 + 1/phi) Lambda^2Lambda^2 phi_*(q^r)  and the number of non-zero
         components (out of the 45 basis vectors of Lambda^2(Lambda^2 V_Sigma_t)) is printed.
"""
import json, re, itertools
from fractions import Fraction as Fr
import sympy as sp

HERE = __file__.rsplit('/', 1)[0] if '/' in __file__ else '.'
t, c = sp.symbols('t c')
I = sp.I
phi = (t**2 + 1)/(2*t)

# ---------------------------------------------------------------- Part 1
print("=== Part 1: letters of phi = (t+1/t)/2 ===")
POINTS = {'e0': sp.Integer(0), 'ei': I, 'e-i': -I, 'e1': sp.Integer(1), 'e-1': sp.Integer(-1)}   # Sigma_t, with this order
NAMES = list(POINTS)

def mult(poly, pt):
    """multiplicity of the root pt of the polynomial poly(t) (over Q(i))."""
    m = 0
    p = sp.Poly(poly, t, extension=I)
    lin = sp.Poly(t - pt, t, extension=I)
    while True:
        q, r = sp.div(p, lin)
        if not r.is_zero:
            return m
        m += 1
        p = q

def order_vector(g):
    num, den = sp.fraction(sp.cancel(sp.together(g)))
    vec = {}
    for nm, pt in POINTS.items():
        o = mult(num, pt) - mult(den, pt)
        if o:
            vec[nm] = o
    # all zeros and poles in C must lie in Sigma_t:
    dn = sp.degree(num, t) if num.has(t) else 0
    dd = sp.degree(den, t) if den.has(t) else 0
    nroots = sum(max(v, 0) for v in vec.values())
    nroots_den = sum(max(-v, 0) for v in vec.values())
    assert nroots == dn and nroots_den == dd, "a zero or pole of g lies outside Sigma_t"
    return vec

CLASS = {}
for a in (0, 1, -1):
    g = sp.cancel(phi - a)
    print(f"  phi - ({a}) = {sp.factor(g, extension=I)}   [phi-({a})] =", order_vector(g))
    CLASS[a] = order_vector(g)
print("  every zero and pole of phi, phi-1, phi+1 in C lies in Sigma_t = {0, i, -i, 1, -1}: True (asserted above)")

# ---------------------------------------------------------------- Part 2
print("\n=== Part 2: the map phi_* and its exterior powers ===")
def vec5(d):
    return [sp.Integer(d.get(nm, 0)) for nm in NAMES]
M = sp.Matrix([vec5(CLASS[0]), vec5(CLASS[1]), vec5(CLASS[-1])]).T      # 5 x 3, columns phi_*(e_0), phi_*(e_1), phi_*(e_-1)
print("  matrix of phi_* (rows e0, ei, e-i, e1, e-1; columns x=e_0, y=e_1, z=e_-1):")
sp.pprint(M)
print("  rank of phi_* =", M.rank(), "(injective iff 3)")
assert M.rank() == 3

# Lambda^2 V_Sigma_t: basis pairs (i<j) of 5 letters; Lambda^2 Lambda^2: pairs of the 10 pairs (45)
PAIRS = list(itertools.combinations(range(5), 2))
PAIRIDX = {p: n for n, p in enumerate(PAIRS)}
PP = list(itertools.combinations(range(10), 2))
PPIDX = {p: n for n, p in enumerate(PP)}

def wedge_vec(u, v):
    """u ^ v in Lambda^2 (R^5) as a list of 10 numbers."""
    return [u[i]*v[j] - u[j]*v[i] for (i, j) in PAIRS]

def wedge_wedge(U, V):
    """U ^ V in Lambda^2 (R^10) as a list of 45 numbers."""
    return [U[i]*V[j] - U[j]*V[i] for (i, j) in PP]

col = {'x': M[:, 0], 'y': M[:, 1], 'z': M[:, 2]}
P_ = wedge_vec(col['x'], col['y'])     # Lambda^2 phi_* (p), p = x^y
Q_ = wedge_vec(col['x'], col['z'])     # q = x^z
R_ = wedge_vec(col['y'], col['z'])     # r = y^z
A = wedge_wedge(P_, R_)                # Lambda^2 Lambda^2 phi_* (p ^ r)
B = wedge_wedge(Q_, R_)                # (q ^ r)
PQ = wedge_wedge(P_, Q_)               # (p ^ q)
rk = sp.Matrix([PQ, A, B]).rank()
print("  rank of Lambda^2 Lambda^2 phi_* on span(p^q, p^r, q^r) =", rk, "(injective iff 3)")
assert rk == 3

# ---------------------------------------------------------------- Part 3
print("\n=== Part 3: the pulled-back (2,2) part ===")
def read_terms(fn):
    out = []
    for nm, v in json.load(open(f"{HERE}/{fn}")):
        m = re.match(r'c\^(-?\d+) H\[(.*)\]', nm)
        if m is None:
            out.append((0, (), Fr(v)))
        else:
            out.append((int(m.group(1)), tuple(int(x) for x in m.group(2).split(',')), Fr(v)))
    return out
R3 = read_terms('g8_terms_k3_part1.json')
NAME = {0: 'x', 1: 'y', -1: 'z'}
T4 = {}
for e, w, q in R3:
    if len(w) != 4:
        continue
    s = (-1)**sum(1 for a in w if a == 1)
    key = tuple(NAME[a] for a in reversed(w))
    cur = T4.setdefault(key, {})
    cur[e] = cur.get(e, 0) + s*q
    T4[key] = {k: v for k, v in cur.items() if v != 0}
T4 = {k: v for k, v in T4.items() if v}
print("  S(Omega) has", len(T4), "pure tensors")
T4c = {k: sum(sp.Rational(v.numerator, v.denominator)*c**e for e, v in d.items()) for k, d in T4.items()}

# letters x,y,z -> vectors in V_Sigma_t (dict point -> coefficient)
IMG = {'x': CLASS[0], 'y': CLASS[1], 'z': CLASS[-1]}
Tpull = {}
for key, coef in T4c.items():
    cf = coef.subs(c, phi)
    for choice in itertools.product(*[list(IMG[l].items()) for l in key]):
        w = tuple(p for p, _ in choice)
        k_ = sp.Integer(1)
        for _, m_ in choice:
            k_ *= m_
        Tpull[w] = Tpull.get(w, 0) + cf*k_

def wedge(a, b, order):
    if a == b:
        return None
    return (1, (a, b)) if order[a] < order[b] else (-1, (b, a))

ORD = {nm: n for n, nm in enumerate(NAMES)}
out = {}
for tt, val in Tpull.items():
    w1 = wedge(tt[0], tt[1], ORD); w2 = wedge(tt[2], tt[3], ORD)
    if w1 is None or w2 is None:
        continue
    s = w1[0]*w2[0]; p1, p2 = w1[1], w2[1]
    if p1 == p2:
        continue
    k1 = (ORD[p1[0]], ORD[p1[1]]); k2 = (ORD[p2[0]], ORD[p2[1]])
    if PAIRIDX[k1] > PAIRIDX[k2]:
        p1, p2, s = p2, p1, -s
    out[(p1, p2)] = out.get((p1, p2), 0) + s*val
D = {k: sp.factor(sp.cancel(sp.together(v))) for k, v in out.items()}
D = {k: v for k, v in D.items() if v != 0}
print("  delta_{2,2}(S_4(M_3 o phi)), computed directly in T^4(V_Sigma_t): non-zero components =", len(D), "of", len(PP))

# the formula by naturality
form = {}
for n in range(45):
    a_, b_ = A[n], B[n]
    v = sp.cancel(sp.together((1 - 1/phi)*a_ - (1 + 1/phi)*b_))
    if v != 0:
        form[n] = v
print("  formula (1-1/phi) L2L2phi_*(p^r) - (1+1/phi) L2L2phi_*(q^r): non-zero components =", len(form))
# compare
def idx(key):
    (p1, p2) = key
    return PPIDX[(PAIRIDX[(ORD[p1[0]], ORD[p1[1]])], PAIRIDX[(ORD[p2[0]], ORD[p2[1]])])]
direct = {idx(k): v for k, v in D.items()}
same = set(direct) == set(form) and all(sp.cancel(direct[n] - form[n]) == 0 for n in form)
print("  direct computation equals the formula by naturality:", same)
assert same and len(form) > 0
# the non-zero components with their coefficients (as functions of t)
inv = {PPIDX[pp]: pp for pp in PP}
def pname(pi):
    i, j = PAIRS[pi]
    return f"({NAMES[i]}^{NAMES[j]})"
print("  components (basis (a^b)^(c^d) of Lambda^2 Lambda^2) with coefficient as a rational function of t:")
for n in sorted(form):
    i, j = inv[n]
    print("    ", pname(i), "^", pname(j), ":", sp.factor(form[n]))
print("  RESULT: phi_* injective (rank 3), Lambda^2Lambda^2 phi_* injective (rank 3), pulled-back delta_{2,2} non-zero with",
      len(form), "non-zero components: True")
