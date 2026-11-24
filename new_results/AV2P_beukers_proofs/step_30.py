from mpmath import mp, mpf, zeta, log, pi, nprint, quad, polylog
mp.dps = 60
def g(c):
    c = mpf(c)
    return (polylog(2, c) + log(c)*mp.log1p(-c))/c
B2_2 = mpf('1.05936833579935324104250490229980343260203365771260126540342')
# form 1: (1/4) int_0^{pi/2} g(sin(2t)/2) dt
F1 = mpf('0.25')*quad(lambda t: g(mp.sin(2*t)/2), [0, pi/2])
print("F1 =", nprint(F1, 40), " diff:", nprint(F1-B2_2, 5))
# form 2: (1/2) int_0^{pi/4} g(sin a) da
F2 = mpf('0.5')*quad(lambda a: g(mp.sin(a)), [0, pi/4])
print("F2 =", nprint(F2, 40), " diff:", nprint(F2-B2_2, 5))
# form 3: (1/4) int_0^{pi/2} g(sin(f/2)) df
F3 = mpf('0.25')*quad(lambda f: g(mp.sin(f/2)), [0, pi/2])
print("F3 =", nprint(F3, 40), " diff:", nprint(F3-B2_2, 5))
# AV2 R2: B2(p) = (1/p) int_0^{pi/4} (cos sin)^{2/p-1} g(c) dt, c=(cos sin)^{2/p}; p=2:
F4 = mpf('0.5')*quad(lambda t: (mp.cos(t)*mp.sin(t))**0 * g((mp.cos(t)*mp.sin(t))**1), [0, pi/4])
print("F4 =", nprint(F4, 40), " diff:", nprint(F4-B2_2, 5))
