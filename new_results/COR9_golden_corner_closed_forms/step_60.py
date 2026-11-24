WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
code = '''"""cor9_rank.py -- exact multiplicative rank of {2, c, 1-c, 1+c} in any real quadratic field,
via prime-ideal valuations (Hensel-lifted) + exact unit test on kernel relations.
rank(Gamma) = rank(V_true) + delta, delta = 1 iff some ideal-level relation yields a
non-torsion unit. No class-number hypothesis, no principal generators needed."""
import math, sys
from fractions import Fraction as Fr
BASE = '<project>/research/'
sys.path.insert(0, BASE + '08_theory/COR9_golden_corner_closed_forms')
from cor9_lib import factorint, fund_unit

def hensel_mu(D, l, k):
    """mu with mu^2 = D mod l^k (l split or ramified, odd or l=2 with care). Returns None if D not a square."""
    if k == 0: return 0
    mod = l
    r = D % l
    if r == 0:
        mu = 0
    else:
        if pow(r, (l-1)//2, l) != 1: return None
        from sympy import sqrt_mod
        mu = int(sqrt_mod(r, l))
    while mod < l**k:
        # lift mu from mod to mod*l
        f = (mu*mu - D)//mod
        fp = (2*mu) % l
        if fp % l == 0:
            # ramified-ish; try all lifts
            found = None
            for t in range(l):
                if (mu + t*mod)**2 % (mod*l) == D % (mod*l):
                    found = mu + t*mod; break
            if found is None: return None
            mu = found
        else:
            t = (-f * pow(int(fp), -1, l)) % l
            mu = mu + t*mod
        mod *= l
    return mu

def v_ideal(x_num, x_den, D, l):
    """valuations of x = x_num/x_den at primes above l.
    x_num = (a, b) meaning a + b sqrt D (Fractions allowed, but assume integers), x_den integer.
    Returns dict {prime_id: v}."""
    out = {}
    a, b = x_num
    # valuation at rational l level for numerator/denominator
    def vL(n):
        n = abs(int(n)); c = 0
        while n and n % l == 0: n //= l; c += 1
        return c
    if D % l == 0 or (D % l and pow(D % l, (l-1)//2, l) != 1):
        # ramified (D % l == 0) or inert: single prime (Galois-invariant)
        vN = vL(a*a - D*b*b) - 2*vL(x_den)
        if D % l == 0:
            out[(l, 0)] = vN          # ramified: v_p(x) = v_l(N x)
        else:
            assert vN % 2 == 0, (D, l, a, b, x_den)
            out[(l, 0)] = vN//2       # inert: v_p(x) = v_l(N x)/2
    else:
        # split: valuations via Hensel-lifted mu
        vN = vL(a*a - D*b*b) - 2*vL(x_den)
        # v at p = (l, sqrt D - mu): largest k with a + b mu_k = 0 mod l^k, capped by denominator correction
        def vsplit(mu_sign):
            k = 0
            while k < 60:
                mu = hensel_mu(D, l, k+1)
                if mu is None: break
                m = mu_sign*mu
                if (a + b*m) % (l**(k+1)) == 0: k += 1
                else: break
            return k - vL(x_den)
        v1 = vsplit(+1)
        v2 = vN - v1
        # sanity: v2 should equal vsplit(-1)
        out[(l, 1)] = v1
        out[(l, 2)] = v2
    return {k: v for k, v in out.items() if v != 0}

def letter_vecs(D, a, b, e):
    """prime-ideal valuation vectors of the 4 letters, exact."""
    letters = [((2, 0), 1), ((a, b), e), ((e - a, -b), e), ((e + a, b), e)]
    norms = [4, a*a - D*b*b, (e-a)*(e-a) - D*b*b, (e+a)*(e+a) - D*b*b]
    primes = set()
    for n in norms:
        primes |= set(factorint(abs(n)).keys())
        if n == 0: raise ValueError("zero letter")
    vecs = []
    for (num, den), n in zip(letters, norms):
        v = {}
        for l in primes:
            for pid, vv in v_ideal(num, den, D, l).items():
                v[pid] = v.get(pid, 0) + vv
        vecs.append(v)
    return vecs, sorted(primes)

def rank_int(M):
    M = [row[:] for row in M if any(row)]
    if not M: return 0
    R = 0
    for c in range(len(M[0])):
        p = next((r for r in range(R, len(M)) if M[r][c]), None)
        if p is None: continue
        M[R], M[p] = M[p], M[R]
        piv = M[R][c]
        for r in range(len(M)):
            if r != R and M[r][c]:
                f = M[r][c]
                M[r] = [x*piv - y*f for x, y in zip(M[r], M[R])]
        R += 1
        if R == len(M): break
    return R

def kernel_basis(M):
    """integer kernel basis of the 4 x P valuation matrix (relations among letters)."""
    import sympy as sp
    S = sp.Matrix(M)
    ns = S.nullspace()
    return [[int(x) for x in v] for v in ns]

def true_rank(D, a, b, e, eps=None, verbose=False):
    """exact rank of <2, c, 1-c, 1+c>, c = (a + b sqrt D)/e."""
    vecs, primes = letter_vecs(D, a, b, e)
    pids = sorted({pid for v in vecs for pid in v})
    M = [[v.get(pid, 0) for pid in pids] for v in vecs]
    rv = rank_int(M)
    ker = kernel_basis(M)
    if not ker:
        return rv, M, pids, []
    if eps is None:
        eps = fund_unit(D)
    ex, ey = eps
    # unit test: for each kernel relation m, y = prod letters^mi is a unit; y = +- eps^k?
    # compute y exactly in NF (pairs of Fractions)
    import mpmath as mp
    mp.mp.dps = 100
    rD = mp.sqrt(D)
    logeps = mp.log(abs(mp.mpf(ex.numerator)/ex.denominator + mp.mpf(ey.numerator)/ey.denominator*rD))
    relinfo = []
    for m in ker:
        # y = 2^m0 * c^m1 * (1-c)^m2 * (1+c)^m3
        ly = (m[0]*mp.log(2) + m[1]*mp.log(abs(mp.mpf(a) + mp.mpf(b)*rD)/e)
              + m[2]*mp.log(abs(mp.mpf(e-a) - mp.mpf(b)*rD)/e)
              + m[3]*mp.log(abs(mp.mpf(e+a) + mp.mpf(b)*rD)/e))
        k = int(mp.nint(ly/logeps)) if abs(ly) > mp.mpf('1e-50') else 0
        # exact NF check: y == eps^k or -eps^k
        def nf_mul(u, v): return (u[0]*v[0] + D*u[1]*v[1], u[0]*v[1] + u[1]*v[0])
        def nf_pow(u, n):
            r = (Fr(1), Fr(0))
            if n < 0: u, n = (u[0]/(u[0]**2 - D*u[1]**2), -u[1]/(u[0]**2 - D*u[1]**2)), -n
            for _ in range(n): r = nf_mul(r, u)
            return r
        def nf_div(u, v): return nf_mul(u, (v[0]/(v[0]**2 - D*v[1]**2), -v[1]/(v[0]**2 - D*v[1]**2)))
        els = [(Fr(2), Fr(0)), (Fr(a, e), Fr(b, e)), (Fr(e-a, e), Fr(-b, e)), (Fr(e+a, e), Fr(b, e))]
        y = (Fr(1), Fr(0))
        for mi, el in zip(m, els):
            y = nf_mul(y, nf_pow(el, mi))
        ek = nf_pow((ex, ey), k)
        ok = (y == ek) or (y == (-ek[0], -ek[1]))
        if not ok:
            return None, M, pids, [('UNITCHECK FAILED', m, k)]
        relinfo.append((m, k))
    delta = 1 if any(k != 0 for _, k in relinfo) else 0
    return rv + delta, M, pids, relinfo
'''
open(WD+'cor9_rank.py','w').write(code)
print("cor9_rank.py written")
import sys; sys.path.insert(0, WD)
from cor9_rank import true_rank
# cross-check against the 5 known rank-2 points and some controls
tests = [(2, 3, -2, 1, 2), (2, -1, 1, 1, 2), (2, 0, 1, 2, 2), (5, -2, 1, 1, 2), (5, -1, 1, 2, 2),
         (2, 17, -12, 1, 3), (5, 3, -1, 2, 3), (3, 0, 1, 2, 3), (17, -4, 1, 1, 3)]
for D, a, b, e, expect in tests:
    r, M, pids, rel = true_rank(D, a, b, e)
    print(f"c0 = ({a},{b})/{e} in Q(sqrt{D}): true_rank = {r} (expect {expect}) {'OK' if r == expect else 'MISMATCH'}")
