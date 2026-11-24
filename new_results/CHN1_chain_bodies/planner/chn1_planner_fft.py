# Planner check of the CHN1 weight-8 result, independent of every engine:
#   c8(a) = [u^8] Vol/16 along u_i = a_i u, by Cauchy/FFT on |u| = r with conjugate symmetry, from
#   Vol/16 = G(u) * 3F2(-u4, u3, 1+u1+u3; 1+u3, 1+u1+u2+u3; 1),  G = Gamma(1+u2)Gamma(1+u1+u3)/Gamma(1+u1+u2+u3).
# (1) blind point a = (2,1,1,3): compare with the reviewers' exact prediction
#     c8 = -1986015/32 z(8) + 50796 z(3)z(5) - 8316 z(2)z(3)^2 + 1134 z(3,5).
# (2) plane point a = (1,2,1,2) (a1=a3, a2=a4): PSLQ on [c8, z8, z3z5, z2z3^2, z(3,5)]; prediction: no z(3,5).
# zeta(3,5) = sum_{n>m} n^-3 m^-5 computed two ways: Hurwitz nsum here, and the pilot's mzv_vals.dz.
import sys, time, json
import mpmath as mp
sys.path.insert(0, '<project>/research/02_ideas/lenses/experimental_pilots')
OUT = '<tmp>'
DPS, N, R = 80, 64, mp.mpf('0.02')
def log(s):
    with open(OUT, 'a') as f: f.write(time.strftime('%H:%M:%S ') + s + '\n')
def vol(a, u):
    u1, u2, u3, u4 = [x * u for x in a]
    G = mp.gamma(1 + u2) * mp.gamma(1 + u1 + u3) / mp.gamma(1 + u1 + u2 + u3)
    return G * mp.hyp3f2(-u4, u3, 1 + u1 + u3, 1 + u3, 1 + u1 + u2 + u3, 1)
def c8(a):
    vals, t0 = {}, time.time()
    for j in range(N // 2 + 1):
        vals[j] = vol(a, R * mp.expjpi(mp.mpf(2 * j) / N))
        if j % 8 == 0: log(f'a={a} point {j}/{N//2} {time.time()-t0:.0f}s')
    for j in range(N // 2 + 1, N): vals[j] = mp.conj(vals[N - j])
    ck = lambda k: mp.fsum(vals[j] * mp.expjpi(-mp.mpf(2 * j * k) / N) for j in range(N)) / N / R ** k
    c = ck(8)
    log(f'a={a} c8 = {mp.nstr(c.real, 40)}  imag {mp.nstr(c.imag, 3)}  c1 {mp.nstr(abs(ck(1)), 3)}')
    return c.real
mp.mp.dps = DPS + 15
z35a = mp.nsum(lambda m: m ** -5 * mp.zeta(3, m + 1), [1, mp.inf])
import mzv_vals
z35b = mzv_vals.dz(3, 5)
log(f'zeta(3,5) Hurwitz-nsum {mp.nstr(z35a, 30)}  pilot dz diff {mp.nstr(z35a - z35b, 3)}')
z = mp.zeta
basis = [z(8), z(3) * z(5), z(2) * z(3) ** 2, z35a]
mp.mp.dps = DPS
cb = c8((2, 1, 1, 3))
mp.mp.dps = DPS + 15
pred = -mp.mpf(1986015) / 32 * basis[0] + 50796 * basis[1] - 8316 * basis[2] + 1134 * basis[3]
log(f'(1) blind (2,1,1,3): c8 - prediction = {mp.nstr(cb - pred, 3)};  without the z(3,5) term: {mp.nstr(cb - pred + 1134 * basis[3], 6)}')
mp.mp.dps = DPS
cp = c8((1, 2, 1, 2))
mp.mp.dps = 64
rel = mp.pslq([cp] + basis, maxcoeff=10 ** 12, maxsteps=10 ** 6)
rel4 = mp.pslq([cp] + basis[:3], maxcoeff=10 ** 12, maxsteps=10 ** 6)
log(f'(2) plane (1,2,1,2): PSLQ5 {rel};  PSLQ4 (no z35) {rel4}')
if rel4:
    mp.mp.dps = DPS
    r = mp.fsum(c * v for c, v in zip(rel4, [cp] + basis[:3]))
    log(f'(2) residual of the 4-term relation at {DPS} dps: {mp.nstr(r, 3)}')
log('DONE')
