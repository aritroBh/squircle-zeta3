"""v_cauchy.py -- INDEPENDENT corner coefficients M_k(alpha), k = 0..3, by a Cauchy contour integral in u
(verifier's own code; does not use corner_series.py or the solver's u-series-at-nodes method).

M(u;alpha) = int_0^1 [ (1+g)^(-1-u) (1+h) sqrt(1 - 4c^2 h/(1+h)^2) - 1 ] g^(u-1) dg,  h = g^(1-u),
c = cos(alpha/2).  (1+h^2-2cos(alpha) h = (1+h)^2 - 4c^2 h; this form keeps the principal sqrt analytic
for |u| <= 0.4 since |4c^2 h/(1+h)^2| < 1 there.)
M_k = (1/N) sum_j M(u_j) u_j^(-k), u_j = r e^(2 pi i j/N); conjugate symmetry M(conj u) = conj M(u).
Usage: python v_cauchy.py DPS N R tag1 tag2 ...   writes vc_<tag>.txt
"""
import sys, time
import mpmath as mp
from multiprocessing import Pool

ANG = {
    'pi_2': lambda: mp.pi/2, '2pi_3': lambda: 2*mp.pi/3, 'pi_3': lambda: mp.pi/3,
    '3pi_4': lambda: 3*mp.pi/4, 'pi_6': lambda: mp.pi/6, '5pi_6': lambda: 5*mp.pi/6,
    'pi_4': lambda: mp.pi/4, 'gen1p3': lambda: mp.mpf(13)/10, 'gen2p0': lambda: mp.mpf(2),
    # fresh angles never used by the solver
    'new0p3': lambda: mp.mpf(3)/10, 'new0p9': lambda: mp.mpf(9)/10, 'new2p6': lambda: mp.mpf(26)/10,
    'new3p0': lambda: mp.mpf(3),
}


def Mu(u, c2):
    def f(g):
        h = mp.power(g, 1 - u)
        x = 4*c2*h/(1 + h)**2
        return (mp.power(1 + g, -1 - u)*(1 + h)*mp.sqrt(1 - x) - 1)*mp.power(g, u - 1)
    return mp.quad(f, [0, mp.mpf(1)/1024, mp.mpf(1)/32, mp.mpf(1)/4, 1])


def job(args):
    tag, dps, N, R, j = args
    mp.mp.dps = dps
    a = ANG[tag]()
    c2 = mp.cos(a/2)**2
    u = R*mp.expjpi(2*mp.mpf(j)/N)
    return tag, j, Mu(u, c2)


if __name__ == '__main__':
    dps, N, R = int(sys.argv[1]), int(sys.argv[2]), mp.mpf(sys.argv[3])
    tags = sys.argv[4:]
    jobs = [(t, dps, N, sys.argv[3], j) for t in tags for j in range(N//2 + 1)]
    t0 = time.time()
    res = {t: {} for t in tags}
    with Pool(12) as p:
        for tag, j, v in p.imap_unordered(job, [(t, dps, N, mp.mpf(R), j) for (t, dps, N, _, j) in jobs]):
            res[tag][j] = v
    mp.mp.dps = dps
    for tag in tags:
        vals = res[tag]
        full = {}
        for j in range(N//2 + 1):
            full[j] = vals[j]
            if 0 < j < N//2:
                full[N - j] = mp.conj(vals[j])
        out = []
        for k in range(4):
            s = mp.fsum(full[j]*mp.expjpi(-2*mp.mpf(j*k)/N) for j in range(N))/N
            out.append(s/mp.mpf(R)**k)
        with open(f'vc_{tag}.txt', 'w') as f:
            f.write(f'# v_cauchy dps={dps} N={N} R={R}; imag parts: ' + ' '.join(mp.nstr(abs(mp.im(x)), 3) for x in out) + '\n')
            for k in range(4):
                f.write(f'M{k} = {mp.nstr(mp.re(out[k]), dps - 5)}\n')
        print(tag, 'done', [mp.nstr(mp.re(x), 25) for x in out], 'imag', [mp.nstr(abs(mp.im(x)), 3) for x in out], flush=True)
    print('elapsed', round(time.time() - t0), 's')
