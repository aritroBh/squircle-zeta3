# Planner check 2 of CHN1 (independent numerics against the CHN1 engine's exact predictions):
#   c_k(a) = [u^k] Vol/16 along u_i = a_i u by Cauchy/FFT with mpmath hyp3f2 (not CHN1's 2F1-integral representation).
#   Plane points (1,2,1,2) in Pi2 and (2,2,1,1) in Pi1: k = 8..12; generic (2,1,1,3): k = 8..10 (tests the zeta(3,7) form).
#   Also: CHN1's weight-12 relation zeta(5,7) = 128195/12438 z12 - 14/9 z3 z9 - 25/3 z5 z7 + 14/9 zeta(3,9), numerically.
import sys, time
from fractions import Fraction as F
import mpmath as mp
CHN = '<project>/research/08_theory/CHN1_chain_bodies'
sys.path.insert(0, CHN)
OUT = '<tmp>'
def log(s):
    with open(OUT, 'a') as f: f.write(time.strftime('%H:%M:%S ') + s + '\n')
import os; os.chdir(CHN)
import engine, mzv_exact as M
import mzv_vals
DPS, N = 90, 72
POINTS = [((1, 2, 1, 2), 12), ((2, 2, 1, 1), 12), ((2, 1, 1, 3), 10)]
pred = {}
for a, K in POINTS:
    t0 = time.time(); pred[a] = engine.reduced_coeffs(a, K)
    log(f'prediction {a} to weight {K}: {time.time()-t0:.0f}s; irreducible slots: ' +
        '; '.join(f'W{W}: ' + str([str(c) for c, mon in zip(pred[a][W], M.BASIS[W]) if any(len(w) > 1 for w in mon)]) for W in range(8, K + 1)))
mp.mp.dps = DPS + 20
DZ = {(3, 5): mzv_vals.dz(3, 5), (3, 7): mzv_vals.dz(3, 7), (3, 9): mzv_vals.dz(3, 9), (5, 7): mzv_vals.dz(5, 7)}
def zw(w):
    if len(w) == 1: return mp.zeta(w[0])
    if w in DZ: return DZ[w]
    raise KeyError(w)
z = mp.zeta
rel = DZ[(5, 7)] - (mp.mpf(128195) / 12438 * z(12) - mp.mpf(14) / 9 * z(3) * z(9) - mp.mpf(25) / 3 * z(5) * z(7) + mp.mpf(14) / 9 * DZ[(3, 9)])
log(f'weight-12 relation residual: {mp.nstr(rel, 3)}')
def value(vec, W):
    tot = mp.mpf(0)
    for c, mon in zip(vec, M.BASIS[W]):
        if c == 0: continue
        tot += mp.mpf(c.numerator) / c.denominator * mp.fprod(zw(w) for w in mon)
    return tot
def vol(a, u):
    u1, u2, u3, u4 = [x * u for x in a]
    return mp.gamma(1 + u2) * mp.gamma(1 + u1 + u3) / mp.gamma(1 + u1 + u2 + u3) * mp.hyp3f2(-u4, u3, 1 + u1 + u3, 1 + u3, 1 + u1 + u2 + u3, 1)
for a, K in POINTS:
    R = 1 / mp.mpf(max(a[1] + a[3], a[0] + a[1] + a[2], a[2], a[0] + a[2], a[1]))
    r = R * mp.mpf('0.08')
    mp.mp.dps = DPS
    vals, t0 = {}, time.time()
    for j in range(N // 2 + 1):
        vals[j] = vol(a, r * mp.expjpi(mp.mpf(2 * j) / N))
        if j % 9 == 0: log(f'a={a} point {j}/{N//2} {time.time()-t0:.0f}s')
    for j in range(N // 2 + 1, N): vals[j] = mp.conj(vals[N - j])
    mp.mp.dps = DPS + 20
    for k in range(8, K + 1):
        ck = mp.fsum(vals[j] * mp.expjpi(-mp.mpf(2 * j * k) / N) for j in range(N)) / N / r ** k
        pv = value(pred[a][k], k)
        log(f'a={a} c{k}: fft-pred = {mp.nstr(ck.real - pv, 3)}  (|c{k}| = {mp.nstr(abs(pv), 5)}, imag {mp.nstr(ck.imag, 2)})')
log('DONE')
