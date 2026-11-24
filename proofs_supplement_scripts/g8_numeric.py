"""g8_numeric.py -- fresh numerical checks for Theorems 20 and 21, for the classical form of Theorem 20 and for M_3 (Lemma R3).
NOT a proof.  M_k(alpha), k=2,3,4, are computed from the DEFINITION  M(u;alpha) = int_0^1 [...] g^(u-1) dg  by a Cauchy integral
in u (radius 0.2, 64 nodes, 40 digits; the aliasing error is of order 0.2^64 M_{k+64}), and compared with
  (A) Theorem 20 in its HPL form  (HPL values from their Taylor series, N terms, tail below 1e-30),
  (B) the classical form B2' of Theorem 20 (mpmath polylog),
  (D) M_3 = R_3 - z3 Phi + z2 Psi, R_3 and Psi the term lists g8_terms_k3_part1.json, g8_terms_k3_partz2.json (the terms
      with c^e, e = -2..1, are summed from the HPL Taylor series; the constant 1/2 of Psi is the term "1"),
  (C) Theorem 21: R_4 + z2 Psi_2 - z3 Psi + pi^4/2880 [...], with R_4, Psi_2, Psi the term lists g8_terms_k4_part1.json,
      g8_terms_k4_partz2.json, g8_terms_k3_partz2.json.
Run:  python g8_numeric.py   (a few minutes)."""
import json, re, sys, time
import mpmath as mp
mp.mp.dps = 40

def Mfun(u, alpha):
    ca = mp.cos(alpha)
    def f(g):
        h = g**(1 - u)
        return ((1 + g)**(-1 - u)*mp.sqrt(1 + g**(2 - 2*u) - 2*ca*h) - 1)*g**(u - 1)
    return mp.quad(f, [0, mp.mpf(1)/2, 1])

def Mk_def(alpha, ks=(2, 3, 4), r=mp.mpf('0.2'), Nn=64):
    vals = [(r*mp.expjpi(2*mp.mpf(t)/Nn), None) for t in range(Nn)]
    vals = [(u, Mfun(u, alpha)) for u, _ in vals]
    return {k: sum(v*u**(-k) for u, v in vals)/Nn for k in ks}

# ---- HPL values from series
def hpl_coeffs(w, N):
    if len(w) == 1:
        a = w[0]
        return [mp.mpf(0)] + [mp.mpf(1)/m if a == 1 else mp.mpf((-1)**(m + 1))/m for m in range(1, N + 1)]
    inner = hpl_coeffs(w[1:], N)
    a = w[0]
    out = [mp.mpf(0)]*(N + 1)
    if a == 0:
        for m in range(1, N + 1): out[m] = inner[m]/m
    else:
        acc = mp.mpf(0)
        part = []
        for m in range(N + 1):
            acc = acc*a + inner[m]; part.append(acc)
        for m in range(1, N + 1): out[m] = part[m - 1]/m
    return out

_cache = {}
def Hval(w, c, N):
    key = (w, N)
    if key not in _cache: _cache[key] = hpl_coeffs(w, N)
    co = _cache[key]
    s = mp.mpf(0); p = mp.mpf(1)
    for m in range(1, N + 1):
        p *= c
        s += co[m]*p
    return s

def nz(w): return sum(1 for a in w if a)

def eval_terms(fname, c, N):
    tot = mp.mpf(0)
    for nm, v in json.load(open(fname)):
        q = mp.mpf(v.split('/')[0])/mp.mpf(v.split('/')[1]) if '/' in v else mp.mpf(v)
        if nm in ('1', 'c'):
            tot += q*(c if nm == 'c' else 1); continue
        m = re.match(r'c\^(-?\d+) (He?)\[(.*)\]', nm)
        e = int(m.group(1)); w = tuple(int(t) for t in m.group(3).split(','))
        val = Hval(w, c, N)
        if m.group(2) == 'He':
            s = 1 if (e + nz(w)) % 2 == 0 else -1
            val += s*Hval(tuple(-a for a in w), c, N)
        tot += q*c**e*val
    return tot

def Y_terms(c, N):
    H = lambda w, x: Hval(w, x, N) if x > 0 else (-1)**nz(w)*Hval(tuple(-a for a in w), -x, N)
    def Y(x):
        return ((1 - x)*(H((0, 0, 1), x) + H((0, 1, -1), x) - H((1, -1, 1), x))
                + (1 - x**2)/(2*x)*H((1, 0, 1), x) - (1 - x)**2/(2*x)*H((1, 0, -1), x))
    return Y(c) + Y(-c)

def M2_hpl(c, N):
    Phi = mp.mpf(1)/2 - (1 - c**2)/(2*c)*mp.atanh(c)
    return mp.zeta(2)*Phi - mp.polylog(2, c**2)/2 + Y_terms(c, N)

def M2_classical(cc):
    lam = mp.log(2); lp = mp.log(1 + cc); lm = mp.log(1 - cc); l0 = mp.log(cc)
    qq = (1 - cc**2)/(2*cc); Phi = mp.mpf(1)/2 - qq*mp.atanh(cc)
    return (mp.zeta(2)*Phi - mp.polylog(2, cc) - mp.polylog(2, -cc)
        + (1 + cc)*mp.polylog(3, cc) + (1 - cc)*mp.polylog(3, -cc)
        + ((1 - cc)/cc)*mp.polylog(3, 1 - cc) - ((1 + cc)/cc)*mp.polylog(3, 1/(1 + cc))
        - 2*cc*mp.polylog(3, (1 + cc)/2)
        + (1 - cc)*(mp.polylog(3, (1 - cc)/(1 + cc)) + mp.polylog(3, -(1 - cc)/(1 + cc)))
        - 2*cc*mp.polylog(3, 2*cc/(1 + cc))
        + mp.polylog(2, (1 + cc)/2)*((1 + cc)*lp - (1 - cc)*lm)
        + qq*(mp.polylog(2, cc) - mp.polylog(2, -cc))*(lp + lm)
        - (1 - cc)*lm**2*lp + (1 - cc)*lm*lp**2 + ((4*cc**2 - cc + 1)/(6*cc))*lp**3
        + qq*l0*lm**2 - (1 - cc)*l0*lm*lp - ((1 + cc**2)/(2*cc))*l0*lp**2
        + lam*((1 - cc)*lm**2 - (1 - cc)*lm*lp - cc*lp**2) + lam**2*((1 - cc)/2)*(lp - lm)
        - mp.pi**2*(((1 - cc)/(6*cc))*lm + ((1 + cc)/(6*cc))*lp)
        + (mp.mpf(7)/4 + 2*cc)*mp.zeta(3) - (cc/6)*mp.pi**2*lam + (cc/3)*lam**3)

def M3_hpl(c, N):
    R3 = eval_terms('g8_terms_k3_part1.json', c, N)
    Psi = eval_terms('g8_terms_k3_partz2.json', c, N)
    Phi = mp.mpf(1)/2 - (1 - c**2)/(2*c)*mp.atanh(c)
    return R3 + mp.zeta(2)*Psi - mp.zeta(3)*Phi

def M4_hpl(c, N):
    R4 = eval_terms('g8_terms_k4_part1.json', c, N)
    Psi2 = eval_terms('g8_terms_k4_partz2.json', c, N)
    Psi = eval_terms('g8_terms_k3_partz2.json', c, N)
    at = mp.atanh(c)
    pi4 = mp.pi**4/2880*((9*c - 30/c + 21/c**3)*at + 23 - 21/c**2)
    return R4 + mp.zeta(2)*Psi2 - mp.zeta(3)*Psi + pi4

if __name__ == '__main__':
    t0 = time.time()
    for name, alpha, N in (('pi/3', mp.pi/3, 700), ('pi/2', mp.pi/2, 300), ('2pi/3', 2*mp.pi/3, 200), ('1.8', mp.mpf('1.8'), 260)):
        c = mp.cos(alpha/2)
        d = Mk_def(alpha)
        a2 = M2_hpl(c, N); b2 = M2_classical(c); a3 = M3_hpl(c, N); a4 = M4_hpl(c, N)
        print(f"alpha = {name}  (c = {mp.nstr(c, 12)}, c^N = {mp.nstr(c**N, 2)}, {time.time()-t0:.0f}s)")
        print(f"   M_2(def) = {mp.nstr(mp.re(d[2]), 30)}")
        print(f"   |M_2(def) - Thm 20 HPL form|   = {mp.nstr(abs(d[2] - a2), 3)}")
        print(f"   |M_2(def) - classical form B2'| = {mp.nstr(abs(d[2] - b2), 3)}")
        print(f"   M_3(def) = {mp.nstr(mp.re(d[3]), 30)}")
        print(f"   |M_3(def) - R_3 - z2 Psi + z3 Phi| = {mp.nstr(abs(d[3] - a3), 3)}")
        print(f"   M_4(def) = {mp.nstr(mp.re(d[4]), 30)}")
        print(f"   |M_4(def) - Thm 21 formula|     = {mp.nstr(abs(d[4] - a4), 3)}")
        sys.stdout.flush()
