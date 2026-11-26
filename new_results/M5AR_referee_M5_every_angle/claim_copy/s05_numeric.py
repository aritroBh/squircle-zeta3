"""s05: NUMERICAL checks (not a proof).
(A) M_5(alpha) from the DEFINITION M(u;alpha) = int_0^1 [(1+g)^(-1-u) sqrt(1+g^(2-2u)-2cos(alpha) g^(1-u)) - 1] g^(u-1) dg
    by a Cauchy integral in u (radius 0.2, 64 nodes, conjugate symmetry; aliasing error ~ 0.2^64 M_69),
    against the theorem: R_5 + z2 [z2]M5 + z3 [z3]M5 + z4 [z4]M5 + z2^2 [z2^2]M5 + z5 [z5]M5 + z2 z3 [z2z3]M5,
    every part summed from the Taylor series of its harmonic polylogarithms (N terms).
(B) at the square corner c = 1/sqrt(2) against the stored 400-digit M5 (research/08_theory/M5a_next_coefficient/M5_400.txt).
Usage: s05_numeric.py A | B"""
import json, re, sys, time
import mpmath as mp

def hpl_coeffs(w, N, memo):
    if (w, N) in memo: return memo[(w, N)]
    if len(w) == 1:
        a = w[0]
        r = [mp.mpf(0)] + [mp.mpf(1)/m if a == 1 else mp.mpf((-1)**(m + 1))/m for m in range(1, N + 1)]
    else:
        inner = hpl_coeffs(w[1:], N, memo); a = w[0]
        r = [mp.mpf(0)]*(N + 1)
        if a == 0:
            for m in range(1, N + 1): r[m] = inner[m]/m
        else:
            acc = mp.mpf(0)
            for m in range(1, N + 1):
                acc = acc*a + inner[m - 1]; r[m] = acc/m
    memo[(w, N)] = r
    return r

def nz(w): return sum(1 for a in w if a)

def eval_terms(fname, c, N, memo):
    tot = mp.mpf(0)
    for nm, v in json.load(open(fname)):
        q = mp.mpf(int(v.split('/')[0]))/int(v.split('/')[1]) if '/' in v else mp.mpf(int(v))
        if nm == '1': tot += q; continue
        m = re.match(r'c\^(-?\d+)$', nm)
        if m: tot += q*c**int(m.group(1)); continue
        m = re.match(r'c\^(-?\d+) (He?)\[(.*)\]', nm)
        e = int(m.group(1)); w = tuple(int(t) for t in m.group(3).split(','))
        co = hpl_coeffs(w, N, memo)
        # c^e (H_w(c) + s H_{-w}(c)) = c^e H_w(c) + (-c)^e H_w(-c): even Taylor terms only, doubled
        val = mp.mpf(0); p = mp.mpf(1)
        for k in range(1, N + 1):
            p *= c
            if (k + e) % 2 == 0: val += co[k]*p
        tot += q*c**e*2*val
    return tot

PARTS = [('1', lambda: 1), ('z2', lambda: mp.zeta(2)), ('z3', lambda: mp.zeta(3)), ('z4', lambda: mp.zeta(4)),
         ('z2^2', lambda: mp.zeta(2)**2), ('z5', lambda: mp.zeta(5)), ('z2*z3', lambda: mp.zeta(2)*mp.zeta(3))]

def M5_theorem(c, N):
    memo = {}
    return sum(z()*eval_terms(f"hpl_fit_k5_{nm}.json", c, N, memo) for nm, z in PARTS)

def Mfun(u, alpha):
    ca = mp.cos(alpha)
    def f(g):
        h = g**(1 - u)
        return ((1 + g)**(-1 - u)*mp.sqrt(1 + g**(2 - 2*u) - 2*ca*h) - 1)*g**(u - 1)
    return mp.quad(f, [0, mp.mpf(1)/4, mp.mpf(1)/2, 1])

def M5_def(alpha, r=mp.mpf('0.2'), Nn=64):
    tot = mp.mpc(0)
    for t in range(Nn//2 + 1):
        u = r*mp.expjpi(2*mp.mpf(t)/Nn)
        v = Mfun(u, alpha)*u**(-5)
        tot += v if t in (0, Nn//2) else 2*mp.re(v)
    return mp.re(tot)/Nn

if __name__ == '__main__':
    t0 = time.time()
    if sys.argv[1] == 'A':
        mp.mp.dps = 60
        for name, alpha, N in (('pi/3', mp.pi/3, 900), ('pi/2', mp.pi/2, 380), ('2pi/3', 2*mp.pi/3, 260),
                               ('1.8', mp.mpf('1.8'), 330), ('5pi/6', 5*mp.pi/6, 130), ('pi/4', mp.pi/4, 1500)):
            c = mp.cos(alpha/2)
            th = M5_theorem(c, N)
            d = M5_def(alpha)
            print(f"alpha = {name} (c = {mp.nstr(c, 12)}, c^N = {mp.nstr(c**N, 2)}):  M5(def) = {mp.nstr(d, 45)}")
            print(f"   theorem          = {mp.nstr(th, 45)}")
            print(f"   |difference|     = {mp.nstr(abs(d - th), 3)}   [{time.time()-t0:.0f}s]", flush=True)
    else:
        mp.mp.dps = 170
        c = 1/mp.sqrt(2)
        ref = mp.mpf(open('../M5a_next_coefficient/M5_400.txt').read().split()[0])
        for N in (900, 1150):
            th = M5_theorem(c, N)
            print(f"square corner, N = {N} (c^N = {mp.nstr(c**N, 2)}): theorem = {mp.nstr(th, 60)}")
            print(f"   |theorem - M5_400| = {mp.nstr(abs(th - ref), 3)}   [{time.time()-t0:.0f}s]", flush=True)
