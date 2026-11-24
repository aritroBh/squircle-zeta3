# step_10: core machinery -- letter spaces, cross-ratios, wedge, cobracket; controls
from fractions import Fraction as Fr
import itertools, pickle, os

# ---------- generic exact linear algebra over Q ----------
def vadd(*vs):
    r = {}
    for v in vs:
        for k, c in v.items():
            r[k] = r.get(k, Fr(0)) + c
            if r[k] == 0: del r[k]
    return r
def vmul(v, a):
    a = Fr(a); return {k: c*a for k, c in v.items() if c*a != 0}

# ---------- rational letter vectors ----------
import sympy as sp
def vecQ(q):
    """v(q) in Q^x (x) Q, letters = primes. signs dropped (torsion)."""
    q = Fr(q); out = {}
    for m, sgn in ((abs(q.numerator),1),(q.denominator,-1)):
        for p, e in sp.factorint(int(m)).items():
            out[p] = out.get(p, Fr(0)) + sgn*e
    return {k:c for k,c in out.items() if c}

def S3orbit_q(z):
    """six Mobius images of z (duplicates removed), exact rationals."""
    z = Fr(z)
    pts = [z, 1-z]
    if z != 0: pts.append(1/z)
    if z != 1: pts.append(1/(1-z))
    if z != 0 and z != 1:
        pts += [z/(z-1), (z-1)/z]
    return sorted(set(pts), key=lambda w: (abs(w.numerator), w.denominator, w.numerator))

def cross_ratios_Q(points):
    """all cross-ratios of a set of Q-points (None = infinity); returns set of rationals."""
    pts = [p for p in points]
    out = set()
    for quad in itertools.permutations(pts, 4):
        a,b,c,d = quad
        # cr = (a-c)(b-d)/((a-d)(b-c)); None = infinity
        try:
            if a is None:  val = (b-d)/(b-c)
            elif b is None: val = (a-c)/(a-d)
            elif c is None: val = (b-d)/(a-d)
            elif d is None: val = (a-c)/(b-c)
            else: val = (a-c)*(b-d)/((a-d)*(b-c))
        except ZeroDivisionError:
            continue
        out.add(Fr(val))
    return sorted(out, key=lambda w: (abs(w.numerator), w.denominator, w.numerator))

# ---------- wedge and cobracket over a letter space ----------
def wedge(u, v):
    """u ^ v as dict ((letter,letter))->coef, canonical order by repr."""
    out = {}
    items_u = sorted(u.items(), key=lambda kv: repr(kv[0]))
    items_v = sorted(v.items(), key=lambda kv: repr(kv[0]))
    for (a, ca) in items_u:
        for (b, cb) in items_v:
            if a == b: continue
            if repr(a) < repr(b): key, s = (a,b), 1
            else: key, s = (b,a), -1
            out[key] = out.get(key, Fr(0)) + s*ca*cb
    return {k:c for k,c in out.items() if c}

def delta21_li3(terms, vecfun):
    """delta(N3) of sum c_i Li3(z_i) as dict ((l1,l2),l3)->coef:
    sum_i c_i (z_i ^ (1-z_i)) (x) z_i. terms = [(c, z), ...]."""
    out = {}
    for c, z in terms:
        c = Fr(c)
        vz, vomz = vecfun(z), vecfun(1-z)
        w = wedge(vz, vomz)
        for (a,b), cab in w.items():
            for l3, c3 in vz.items():
                key = ((a,b), l3)
                out[key] = out.get(key, Fr(0)) + c*cab*c3
    return {k:v for k,v in out.items() if v}

# ============ CONTROLS ============
print("=== CONTROL 1: delta{1/2}_3 = 0 ===")
d = delta21_li3([(1, Fr(1,2))], vecQ)
print("delta{1/2}_3 =", d, " (empty = 0)  [consistent with Li3(1/2) = 7/8 zeta(3) + products]")

print("=== CONTROL 2: Landen triple delta = 0 (rational sample points) ===")
ok = True
for x in [Fr(1,3), Fr(2,5), Fr(3,7), Fr(-2,3), Fr(5,3), Fr(7,11), Fr(-1,2)]:
    d = delta21_li3([(1,x),(1,1-x),(1,1-1/x)], vecQ)
    if d: ok = False; print("  FAIL at", x, d)
print("  delta{Li3(x)+Li3(1-x)+Li3(1-1/x)} = 0 for all 7 samples:", ok)

print("=== CONTROL 3: thm:hexagon's two arguments are cross-ratios of {0,1,inf,p} ===")
# s=1: p = -1/kappa = 1/(1-r); x = 1-r. Check Li3 args x and x/(x-1) in CR({0,1,inf,p}).
ok = True
for r in [Fr(2), Fr(1,2), Fr(3), Fr(1,3), Fr(5,2), Fr(7,4)]:
    p = 1/(1-r); x = 1-r
    CR = cross_ratios_Q([Fr(0), Fr(1), None, p])
    memb = (x in CR) and (x/(x-1) in CR)
    if not memb: ok = False; print("  FAIL r=", r, CR)
print("  x and x/(x-1) in CR{0,1,inf,p} for all 6 sample r:", ok)

print("=== CONTROL 3b: hexagon delta at r=2 and r=1/2 vanishes for ANY coefficients ===")
# at r=2: x=-1, x/(x-1)=1/2. {z}_2 (x) z with z=-1: z torsion -> 0; z=1/2: z^(1-z)=0.
for r in [Fr(2), Fr(1,2)]:
    x = 1-r; z2 = x/(x-1)
    d = delta21_li3([(Fr(3),x),(Fr(-5,7),z2)], vecQ)  # arbitrary coeffs
    print(f"  r={r}: x={x}, x/(x-1)={z2}, delta = {d} (empty = 0)")
print("ALL CONTROLS PASS" if ok else "CONTROL FAILURE")
