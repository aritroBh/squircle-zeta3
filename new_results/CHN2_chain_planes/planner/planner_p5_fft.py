# Planner (W4) independent check of CHN2 Part C (P5 chain body), weight 8.
# Representation (own derivation, not CHN2's A+A_rev series): with I2(y) = y^u2 2F1(-u1,u2;1+u2;y), I4(y) = y^u4 2F1(-u5,u4;1+u4;y),
#   Vol/32 = u3 int_0^1 t^(u3-1) h(t) dt,  h(t) = I2(1-t) I4(1-t),
# analytically continued to small complex u as  h(0) + u3 int_0^1 t^(u3-1) (h(t) - h(0)) dt.
# Cauchy/FFT for c8 = [u^8] Vol/32 along u_i = a_i u; compare with the CHN2 engine's predicted full c8 vector
# (kimilib.full_coords_at), and PSLQ at the diagonal to test that zeta(3,5) is genuinely present there.
import sys, os, time, pickle
from fractions import Fraction as F
import mpmath as mp
HERE = '<project>/research/08_theory/CHN2_chain_planes'
OUT = HERE + '/planner/planner_p5_fft.out'
def log(s):
    with open(OUT, 'a') as f: f.write(time.strftime('%H:%M:%S ') + s + '\n')
os.chdir(HERE); sys.path.insert(0, HERE)
import kimilib, engine5, mzv_vals
forms = pickle.load(open('forms10.pkl', 'rb'))
POINTS = [(1, 1, 1, 1, 1), (2, 1, 3, 1, 2)]
pred = {a: kimilib.full_coords_at(a, 8, forms)[8] for a in POINTS}
for a in POINTS:
    log(f'engine prediction c8{a} on [z8, z3z5, z2z3^2, z(3,5)] = {[str(c) for c in pred[a]]}')
DPS, N = 55, 64
def vol(a, u):
    u1, u2, u3, u4, u5 = [x * u for x in a]
    h = lambda t: (1 - t) ** (u2 + u4) * mp.hyp2f1(-u1, u2, 1 + u2, 1 - t) * mp.hyp2f1(-u5, u4, 1 + u4, 1 - t)
    h0 = mp.hyp2f1(-u1, u2, 1 + u2, 1) * mp.hyp2f1(-u5, u4, 1 + u4, 1)
    return h0 + u3 * mp.quad(lambda t: t ** (u3 - 1) * (h(t) - h0), [0, mp.mpf('0.001'), mp.mpf('0.1'), mp.mpf('0.5'), 1])
mp.mp.dps = 30
v = vol((1, 1, 1, 1, 1), mp.mpf(1))
log(f'anchor u=(1,1,1,1,1): Vol/32 = {mp.nstr(v, 25)}  (2/15 = {mp.nstr(mp.mpf(2)/15, 25)})')
for a in POINTS:
    R = 1 / mp.mpf(max(a[0], a[4], a[1] + a[3], a[0] + a[2], a[2] + a[4]))
    r = R / 6
    mp.mp.dps = DPS
    vals, t0 = {}, time.time()
    for j in range(N // 2 + 1):
        vals[j] = vol(a, r * mp.expjpi(mp.mpf(2 * j) / N))
        if j % 8 == 0: log(f'a={a} point {j}/{N//2} {time.time()-t0:.0f}s')
    for j in range(N // 2 + 1, N): vals[j] = mp.conj(vals[N - j])
    c8 = mp.fsum(vals[j] * mp.expjpi(-mp.mpf(16 * j) / N) for j in range(N)) / N / r ** 8
    z35 = mzv_vals.dz(3, 5)
    basis = [mp.zeta(8), mp.zeta(3) * mp.zeta(5), mp.zeta(2) * mp.zeta(3) ** 2, z35]
    pv = mp.fsum(mp.mpf(c.numerator) / c.denominator * b for c, b in zip(pred[a], basis))
    log(f'a={a}: c8 = {mp.nstr(c8.real, 30)} (imag {mp.nstr(c8.imag, 2)}); c8 - prediction = {mp.nstr(c8.real - pv, 3)}; '
        f'without the z(3,5) term: {mp.nstr(c8.real - pv + mp.mpf(pred[a][3].numerator) / pred[a][3].denominator * z35, 5)}')
    mp.mp.dps = 42
    rel = mp.pslq([c8.real] + basis, maxcoeff=10 ** 9, maxsteps=10 ** 6)
    rel4 = mp.pslq([c8.real] + basis[:3], maxcoeff=10 ** 9, maxsteps=10 ** 6)
    log(f'a={a}: PSLQ5 {rel};  PSLQ4 without z(3,5): {rel4}')
log('DONE')
