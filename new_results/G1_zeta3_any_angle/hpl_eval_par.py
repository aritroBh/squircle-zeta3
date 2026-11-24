"""hpl_eval_par.py -- run the hpl_eval.py comparison (M0, M1, M2 general-angle formulas vs series values)
at high precision, one process per angle.  Usage: python hpl_eval_par.py DPS"""
import sys
from multiprocessing import Pool

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 300


def run(tag):
    import mpmath as mp
    import hpl_eval as he
    mp.mp.dps = DPS + 15
    a = {'pi_2': mp.pi/2, '2pi_3': 2*mp.pi/3, 'pi_3': mp.pi/3, '3pi_4': 3*mp.pi/4, 'pi_6': mp.pi/6,
         '5pi_6': 5*mp.pi/6, 'pi_4': mp.pi/4, 'gen1p3': mp.mpf(13)/10, 'gen2p0': mp.mpf(2)}[tag]
    c = mp.cos(a/2)
    d = [he.M0f(c) - he.stored(tag, 0), he.M1f(c) - he.stored(tag, 1), he.M2f(c) - he.stored(tag, 2)]
    return tag, [mp.nstr(abs(x), 3) for x in d], mp.nstr(he.M2f(c), DPS)


if __name__ == '__main__':
    tags = ['pi_2', '2pi_3', 'pi_3', '3pi_4', 'pi_6', '5pi_6', 'pi_4', 'gen1p3', 'gen2p0']
    with Pool(9) as p:
        for tag, d, v in p.imap_unordered(run, tags):
            print(f'{tag:>7}: |M0-ser| = {d[0]}, |M1-ser| = {d[1]}, |M2-ser| = {d[2]}', flush=True)
