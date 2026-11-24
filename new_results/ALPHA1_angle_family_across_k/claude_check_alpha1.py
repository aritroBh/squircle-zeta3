# M1 session: independent checks of ALPHA1 (own code). (a) polynomial law from the four exact kernel pickles (exact Fractions);
# (b) the universal weight-two kernel in the pickles and W2(alpha) = pi^2/12 - Li2(c^2)/2 by own quadrature at several c;
# (c) the weight-three slope = s(alpha) J(c) by own quadrature of the universal slope kernels.
import pickle
from fractions import Fraction as F
import mpmath as mp
mp.mp.dps = 40
D = '<project>/research/08_theory/ALPHA1_angle_family_across_k'
bad = 0; nseq = 0
for tag in ['c12', 'c35', 'c13', 'c45']:
    K = pickle.load(open(f'{D}/kern_{tag}.pkl', 'rb')); c = F(*K['c']); r = F(*K['r']); ks = K['kernels']
    assert r == (1 - c)/(1 + c)
    def res(k, key):
        P = ks[k].get(key)
        if P is None: return {}
        out = {}
        for p, (cr, ci) in P['real'].items():
            assert F(*ci) == 0
            if F(*cr): out[F(*p)] = F(*cr)
        assert not P['ppos'] and not P['pneg']
        return out
    kmax = max(ks)
    for d in range(1, 6):
        for key in [(d - j, j) for j in range(d//2 + 1)]:
            kk = [k for k in range(d + 1, kmax + 1) if key in ks[k]]
            poles = set().union(*[set(res(k, key)) for k in kk])
            for p in poles:
                ys = [res(k, key).get(p, F(0)) for k in kk]
                cols = [ys]
                for _ in range(d): cols.append([cols[-1][t+1] - cols[-1][t] for t in range(len(cols[-1]) - 1)])
                nseq += 1
                if any(cols[d]): bad += 1; print(f'  NOT POLY {tag} d={d} key={key} pole={p}')
    # universal weight-two kernel 1/(2w) + 1/(2(w-1)) - 1/(w+1) - 1/(2(w+r)) for k = 2..kmax
    want = {F(0): F(1, 2), F(1): F(1, 2), F(-1): F(-1), -r: F(-1, 2)}
    for k in range(2, kmax + 1):
        if res(k, (1, 0)) != want: bad += 1; print(f'  W2 kernel mismatch {tag} k={k}: {res(k, (1, 0))}')
print(f'(a)+(b-kernel): {nseq} residue sequences polynomial of degree <= d-1 at 4 angles, W2 kernel universal; failures: {bad}')
def checks(c):
    c = mp.mpf(c); r = (1 - c)/(1 + c)
    t = lambda w: (w - r)*(w + 1)/(2*(1 - r)*w); omt = lambda w: (1 - w)*(w + r)/(2*(1 - r)*w)
    lt = lambda w: mp.log(t(w)); lmt = lambda w: mp.log(omt(w))
    K10 = lambda w: 1/(2*w) + 1/(2*(w - 1)) - 1/(w + 1) - 1/(2*(w + r))
    W2 = mp.quad(lambda w: K10(w)*lt(w), [r, (r + 1)/2, 1])
    T = mp.pi**2/12 - mp.polylog(2, c**2)/2
    s = (1 - c**2)/(2*c**2)
    slope = mp.quad(lambda w: s*(1/(w + 1) - 1/(w + r))*lt(w)**2 + s*(1/w - 2/(w + 1))*lt(w)*lmt(w), [r, (r + 1)/2, 1])
    J = mp.log(1 - c**2)*(mp.polylog(2, c**2) - mp.zeta(2)) + 2*mp.polylog(3, 1 - c**2) - 2*mp.zeta(3) + mp.log(c**2)*mp.log(1 - c**2)**2
    return W2 - T, slope - s*J, slope
for c in ['0.3', '0.5', mp.sqrt(2)/2, '0.8', '0.95']:
    e2, e3, sl = checks(c)
    print(f'c={mp.nstr(mp.mpf(c), 6)}: W2 - T(alpha) = {mp.nstr(e2, 3)};  slope - s J(c) = {mp.nstr(e3, 3)};  slope = {mp.nstr(sl, 20)}')
print('squircle slope target -(z3/8+pi^2 log2/24+log^3 2/12) =', mp.nstr(-(mp.zeta(3)/8 + mp.pi**2*mp.log(2)/24 + mp.log(2)**3/12), 20))
