"""m2_numerics.py -- high-precision corner coefficients M_0, M_1, M_2 (and M_3) at several interior
angles alpha, by the series method of corner_series.py (power series in u at each tanh-sinh node).
Runs each angle at two working precisions and reports the agreement; writes values to
values_<tag>.txt (one file per angle).  Usage: python m2_numerics.py DPS
"""
import sys, time
import mpmath as mp
from multiprocessing import Pool
from corner_series import corner_coeffs

ANGLES = {
    'pi_2':   ('pi/2',   lambda: mp.pi/2),
    '2pi_3':  ('2pi/3',  lambda: 2*mp.pi/3),
    'pi_3':   ('pi/3',   lambda: mp.pi/3),
    '3pi_4':  ('3pi/4',  lambda: 3*mp.pi/4),
    'pi_6':   ('pi/6',   lambda: mp.pi/6),
    '5pi_6':  ('5pi/6',  lambda: 5*mp.pi/6),
    'pi_4':   ('pi/4',   lambda: mp.pi/4),
    'gen1p3': ('1.3',    lambda: mp.mpf(13)/10),
    'gen2p0': ('2.0',    lambda: mp.mpf(2)),
}


def run(args):
    tag, dps = args
    name, af = ANGLES[tag]
    t0 = time.time()
    out = {}
    for d in (dps, dps + 30):
        mp.mp.dps = d
        a = af()
        out[d] = corner_coeffs(a, 3)
    mp.mp.dps = dps + 30
    diffs = [abs(out[dps][k] - out[dps + 30][k]) for k in range(4)]
    with open(f'values_{tag}.txt', 'w') as f:
        f.write(f'# alpha = {name}; corner_series at dps {dps} and {dps+30}; |diff| per k below\n')
        for k in range(4):
            f.write(f'# k={k} diff={mp.nstr(diffs[k], 5)}\n')
        for k in range(4):
            f.write(f'M{k} = {mp.nstr(out[dps + 30][k], dps)}\n')
    return tag, name, [mp.nstr(x, 5) for x in diffs], time.time() - t0


if __name__ == '__main__':
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    tags = sys.argv[2:] if len(sys.argv) > 2 else list(ANGLES)
    with Pool(min(len(tags), 12)) as p:
        for tag, name, diffs, dt in p.imap_unordered(run, [(t, dps) for t in tags]):
            print(f'{name:>6}: diffs(dps vs dps+30) M0..M3 = {diffs}  [{dt:.0f}s]', flush=True)
