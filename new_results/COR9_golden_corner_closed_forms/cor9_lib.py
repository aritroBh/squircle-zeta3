"""cor9_lib.py -- shared machinery for COR9.
Reuses COR7 (cor7_lib: M3 weight split from G1's proved B3 term lists; quadrature with
precision boost; LLL relation finder) and cor7_pt/cor6 (exact S-unit + symbol machinery).
Adds: fundamental unit / norm-l elements via continued fractions (arbitrary regulator),
class number via cycles of reduced indefinite forms, Q(sqrt5) period bases.
"""
import sys, math
BASE = '<project>/research/'
sys.path.insert(0, BASE + '08_theory/COR7_classical_corners')
sys.path.insert(0, BASE + '08_theory/COR6_true_point_symbols')
import mpmath as mp
from fractions import Fraction as Fr
import cor7_lib as L7
from cor7_pt import QuadField, QuadPoint, factorint
from cor6_sym import W4Quotient, lneg
from cor6_ft import ft_coords, ft_delta22
from cor6_pipe import W4TERMS

Q2 = W4Quotient(2)

def li4_coords(P, x):
    vx, vo = P.vec(x), P.vec(P.F.sub(P.F.one(), x))
    return Q2.coords({(lneg(vo), vx, vx, vx): Fr(1)})

def omega_class(P):
    return ft_coords(Q2, P.omega_symbol(W4TERMS), 2)

def rankQ(vecs):
    M = [list(v) for v in vecs]; R = 0
    for c in range(len(M[0])):
        p = next((r for r in range(R, len(M)) if M[r][c]), None)
        if p is None: continue
        M[R], M[p] = M[p], M[R]
        q = M[R][c]; M[R] = [x/q for x in M[R]]
        for r in range(len(M)):
            if r != R and M[r][c]:
                M[r] = [a - b*M[R][c] for a, b in zip(M[r], M[R])]
        R += 1
    return R

# ---------- continued fractions of quadratic irrationals ----------
def cf_period(P, Q, D):
    """CF of (P+sqrt D)/Q with Q | D - P^2, D nonsquare. Returns (a0, period)."""
    sq = math.isqrt(D)
    assert (D - P*P) % Q == 0
    a0 = (P + sq)//Q
    # adjust floor for negative Q
    while (P + sq) - a0*Q < 0: a0 -= 1
    while (P + sq) - (a0+1)*Q >= 0: a0 += 1
    seq = []
    seen = {(P, Q): 0}
    Pi, Qi, ai = P, Q, a0
    while True:
        Pn = ai*Qi - Pi
        Qn = (D - Pn*Pn)//Qi
        an = (Pn + sq)//Qn
        while (Pn + sq) - an*Qn < 0: an -= 1
        while (Pn + sq) - (an+1)*Qn >= 0: an += 1
        if (Pn, Qn) in seen:
            start = seen[(Pn, Qn)]
            return a0, [a for _, _, a in seq[start:]] + [an] if False else [s[2] for s in seq[start:]] + [an]
        seen[(Pn, Qn)] = len(seq)+1
        seq.append((Pi, Qi, ai))
        Pi, Qi, ai = Pn, Qn, an

def convergents(a0, per, nmax):
    hs = [1, a0]; ks = [0, 1]
    out = [(hs[1], ks[1])]
    i = 0
    while len(out) < nmax:
        a = per[i % len(per)]
        hs.append(a*hs[-1] + hs[-2]); ks.append(a*ks[-1] + ks[-2])
        out.append((hs[-1], ks[-1]))
        i += 1
    return out

def fund_unit(D):
    """fundamental unit > 1 of Q(sqrt D) (D squarefree), as (Fr(x,e), Fr(y,e)): (x+y sqrt D)/e."""
    cands = []
    if D % 4 == 1:
        a0, per = cf_period(1, 2, D)
        L = len(per)
        for p, q in convergents(a0, per, 2*L+2):
            n = p*p - p*q - (D-1)//4 * q*q
            if abs(n) == 1:
                cands.append((Fr(2*p-q, 2), Fr(q, 2)))
    else:
        a0, per = cf_period(0, 1, D)
        L = len(per)
        for h, k in convergents(a0, per, 2*L+2):
            n = h*h - D*k*k
            if abs(n) == 1:
                cands.append((Fr(h), Fr(k)))
    best = None
    for x, y in cands:
        v = abs(float(x) + float(y)*D**0.5)
        if v <= 0: continue
        lg = abs(math.log(v))
        if lg > 1e-12 and (best is None or lg < best[0]):
            best = (lg, (x, y))
    return best[1]

def norm_l_element(D, l):
    """element of O_F with |N| = l (l rational prime), via CF of (mu + sqrt D)/l for mu^2 = D mod l;
    reduction of the prime ideal; returns ((Fr x, Fr y), e) with e in {1,2} or None if not principal."""
    # solve mu^2 = D (mod l)
    if D % l == 0:
        mus = [0]
    else:
        r = pow(D % l, (l-1)//2, l)
        if r != 1: return None, 'inert'
        # tonelli-shanks via sympy
        from sympy import sqrt_mod
        m0 = sqrt_mod(D % l, l)
        mus = [m0, l - m0]
    best = None
    sq = math.isqrt(D)
    for mu in mus:
        # ideal (l, mu - sqrt D); element search via CF reduction of the form (l, 2mu, (mu^2-D)/l)
        # use cycle: form f = (a,b,c) with a = l; reduce; track transformation to principal form.
        # Instead: scan CF of (mu + sqrt D)/l : convergents p/q give norm (p^2 - D q^2)/l^2 * l ...
        # Generalized Pell x^2 - D y^2 = +-l (or +-4l): use Lagrange bounds via CF of sqrt D.
        pass
    # fallback: diop_DN
    from sympy.solvers.diophantine.diophantine import diop_DN
    tgts = [-l, l] + ([-4*l, 4*l] if D % 4 == 1 else [])
    for tgt in tgts:
        e = 2 if abs(tgt) == 4*l else 1
        try:
            sols = diop_DN(D, tgt)
        except Exception:
            sols = []
        for x, y in sols:
            if e == 2 and (x - y) % 2 != 0: continue
            nn = Fr(x, e)**2 - D*Fr(y, e)**2
            if abs(nn) != l: continue
            v = abs(float(x)/e + float(y)/e*D**0.5)
            if v <= 0: continue
            lg = abs(math.log(v))
            if best is None or lg < best[0]:
                best = (lg, (Fr(x, e), Fr(y, e)))
    return (best[1] if best else None), ('principal' if best else 'notprincipal')

# ---------- class number via cycles of reduced indefinite primitive forms ----------
def h_plus_order(disc):
    """narrow class number of the quadratic order of discriminant disc > 0 (nonsquare):
    number of cycles of reduced primitive forms (a,b,c), b^2-4ac = disc."""
    sq = math.isqrt(disc)
    red = []
    for b in range(1, sq+1):
        if (disc - b*b) % 4: continue
        num = b*b - disc   # = 4ac < 0
        for a in range(1, sq+1):
            if num % (4*a): continue
            c = num//(4*a)
            if c == 0: continue
            if math.gcd(math.gcd(a, b), abs(c)) != 1: continue
            # reduced: 2a - b < sqrt(disc) < 2a + b  (exact integer form)
            if disc < (2*a + b)*(2*a + b) and (2*a <= b or (2*a - b)*(2*a - b) < disc):
                red.append((a, b, c))
    red += [(-a, b, -c) for a, b, c in red]
    def rho(f):
        a, b, c = f
        t = -((b + sq)//(-(2*c))) if c < 0 else (b + sq)//(2*c)
        # ceil((b+sq)/(2c)) for c<0, floor for c>0
        t = -((b + sq)//(-2*c)) if c < 0 else (b + sq)//(2*c)
        bp = -b + 2*c*t
        ap = (bp*bp - disc)//(4*c)
        return (c, bp, ap)
    seen = {}
    ncyc = 0
    for f in red:
        if f in seen: continue
        ncyc += 1
        g = f
        while g not in seen:
            seen[g] = ncyc
            g = rho(g)
    return ncyc

def class_number(D):
    """ordinary class number of Q(sqrt D), D squarefree > 1."""
    disc = D if D % 4 == 1 else 4*D
    hp = h_plus_order(disc)
    x, y = fund_unit(D)
    N = x*x - D*y*y
    return hp if N == -1 else hp//2
