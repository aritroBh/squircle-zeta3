# Planner (W4): validate CHN3's point engine at weight 11 by independent numerics.
# c11 at a = (1,2,3,4,5) from a Cauchy/FFT of the h(t) representation (own derivation; see CHN2 planner_p5_fft.py),
# compared with the engine's exact c11 (test_pt_w1112.pkl, basis BASIS[11] incl. zeta(3)zeta(3,5) and zeta(3,3,5)).
# zeta(3,3,5) = sum_{n1>n2>n3} n1^-3 n2^-3 n3^-5 by nested partial sums + Richardson (all indices >= 2, so valid),
# computed at two settings to bound its error; zeta(3,5) by the pilot's dz.
import sys, os, time, pickle
import mpmath as mp
HERE = '<project>/research/08_theory/CHN3_p5_locus'
os.chdir(HERE); sys.path.insert(0, HERE)
OUT = HERE + '/planner/planner_c11_fft.out'
def log(s):
    with open(OUT, 'a') as f: f.write(time.strftime('%H:%M:%S ') + s + '\n')
import mzv_vals, mzv_num, mzv_exact as M
d = pickle.load(open('test_pt_w1112.pkl', 'rb')); a = d['pt']; c11e = d['c11']
mp.mp.dps = 80
z335a = mzv_num.mzv_num((3, 3, 5), npts=18, N0=500, step=160)
z335b = mzv_num.mzv_num((3, 3, 5), npts=16, N0=400, step=140)
log(f'zeta(3,3,5) = {mp.nstr(z335a, 50)}; two-setting difference {mp.nstr(z335a - z335b, 3)}')
z35 = mzv_vals.dz(3, 5)
def mono(m):
    v = mp.mpf(1)
    for w in m:
        v *= mp.zeta(w[0]) if len(w) == 1 else (z35 if w == (3, 5) else z335a)
    return v
basis = [mono(m) for m in M.BASIS[11]]
pred = mp.fsum(mp.mpf(c.numerator) / c.denominator * b for c, b in zip(c11e, basis))
log(f'basis W=11: {M.BASIS[11]}; engine c11 = {mp.nstr(pred, 40)}')
DPS, N = 75, 80
R = 1 / mp.mpf(max(a[0], a[4], a[1] + a[3], a[0] + a[2], a[2] + a[4])); r = R / 6
def vol(u):
    u1, u2, u3, u4, u5 = [x * u for x in a]
    h = lambda t: (1 - t) ** (u2 + u4) * mp.hyp2f1(-u1, u2, 1 + u2, 1 - t) * mp.hyp2f1(-u5, u4, 1 + u4, 1 - t)
    h0 = mp.hyp2f1(-u1, u2, 1 + u2, 1) * mp.hyp2f1(-u5, u4, 1 + u4, 1)
    return h0 + u3 * mp.quad(lambda t: t ** (u3 - 1) * (h(t) - h0), [0, mp.mpf('0.001'), mp.mpf('0.1'), mp.mpf('0.5'), 1])
mp.mp.dps = DPS
vals, t0 = {}, time.time()
for j in range(N // 2 + 1):
    vals[j] = vol(r * mp.expjpi(mp.mpf(2 * j) / N))
    if j % 10 == 0: log(f'point {j}/{N//2} {time.time()-t0:.0f}s')
for j in range(N // 2 + 1, N): vals[j] = mp.conj(vals[N - j])
c11 = mp.fsum(vals[j] * mp.expjpi(-mp.mpf(22 * j) / N) for j in range(N)) / N / r ** 11
c8 = mp.fsum(vals[j] * mp.expjpi(-mp.mpf(16 * j) / N) for j in range(N)) / N / r ** 8
mp.mp.dps = 80
i335, i3_35 = 8, 7
log(f'c11 FFT = {mp.nstr(c11.real, 40)} (imag {mp.nstr(c11.imag, 2)}); FFT - engine = {mp.nstr(c11.real - pred, 3)}; '
    f'drop zeta(3,3,5) term: {mp.nstr(c11.real - pred + mp.mpf(c11e[i335].numerator) * basis[i335], 5)}; '
    f'drop zeta3 zeta(3,5) term: {mp.nstr(c11.real - pred + mp.mpf(c11e[i3_35].numerator) * basis[i3_35], 5)}')
log('DONE')
