from mpmath import mp, mpf, zeta, log, pi, nprint, quad, inf, exp, loggamma, digamma
mp.dps = 110

def phi_direct(w):
    w = mpf(w)
    return w**3 * exp(2*loggamma(w) - loggamma(2*w+1)) * (digamma(2*w+1) - digamma(w))

# test stability at very small w
for w0 in ['1e-3','1e-10','1e-30']:
    w = mpf(w0)
    v = phi_direct(w)
    print(w0, nprint(v, 20), " (expect ~1 - zeta3 w^3)")

# C2 integral via quad on split intervals, direct phi everywhere
ints = [(0, mpf('0.25')), (mpf('0.25'), 1), (1, 4), (4, 16), (16, 64), (64, 200), (200, inf)]
tot = mpf(0)
for a,b in ints:
    v = quad(lambda w: (phi_direct(w)-1)/w**3, [a, b])
    tot += v
C2 = tot
print("C2 =", nprint(C2, 100))
print("C2 + pi^2/12 =", nprint(C2 + pi**2/12, 6))
