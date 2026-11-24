
import os
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
step = r'''"""WT6S step E (part 2): E_{5,m} values, c(6) check, independent W6(k) quadrature.

E_{5,m} = sum e_{ij,r,m} N_{(i,j),r} (exact e from e5m_binomial_coeffs.json, N from n6_values.json).
Checks: (E-a) c(6) = E_{5,5}/24 vs the LEAD constant; (E-b) W6(k) quadrature directly from
kern_6/7/10.pkl full degree-5 kernels vs sum_m C(k-6,m-1) E_{5,m}."""
import json, pickle, time
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 170
r2 = mp.sqrt(2); aN = 1 + r2; bN = r2 - 1
POLE = {"0": mp.mpf(0), "a": aN, "-a": -aN, "b": bN, "-b": -bN}
C6_LEAD = "0.335687558224020357007672279844"
PAIR_TO_POLE = {(0, 0): '0', (1, 1): 'a', (-1, -1): '-a', (-1, 1): 'b', (1, -1): '-b'}
KEYS = ["50", "41", "32"]
PRS = ["0", "a", "-a", "-b"]

def main():
    t0 = time.time()
    N = json.load(open('n6_values.json'))['values']
    Nv = {kk: mp.mpf(v) for kk, v in N.items()}
    e = json.load(open('e5m_binomial_coeffs.json'))
    Evals = {}
    for m in range(1, 6):
        tot = mp.mpf(0)
        for key in KEYS:
            i, j = int(key[0]), int(key[1])
            for pole in PRS:
                cell = e[f"{key}|{pole}"]
                crat, csq = F(cell['rat'][m - 1]), F(cell['sqrt2'][m - 1])
                if crat == 0 and csq == 0:
                    continue
                c = mp.mpf(crat.numerator)/crat.denominator + r2*mp.mpf(csq.numerator)/csq.denominator
                tot += c * Nv[f"{i},{j}|{pole}"]
        Evals[f"E5,{m}"] = tot
        print(f"E_{{5,{m}}} = {mp.nstr(tot, 60)}  ({time.time()-t0:.0f}s)", flush=True)
    json.dump({'dps': 150, 'values': {kk: mp.nstr(v, 150) for kk, v in Evals.items()}},
              open('e5m_values.json', 'w'), indent=1)
    # E-a: c(6)
    c6 = Evals["E5,5"] / 24
    dc6 = abs(c6 - mp.mpf(C6_LEAD))
    print(f"c(6) = E_5,5/24 = {mp.nstr(c6, 60)}")
    print(f"c(6) - LEAD constant = {mp.nstr(dc6, 3)}  (>=40 digits check)")
    # E-b: direct W6(k) quadrature from kern_k.pkl
    def w6_quad_direct(kf):
        K = pickle.load(open(kf, 'rb'))
        def l1(v): return mp.log((v - bN) * (v + aN) / (4 * v))
        def l2(v): return mp.log((aN - v) * (v + bN) / (4 * v))
        tot = mp.mpf(0)
        for key in [(5, 0), (4, 1), (3, 2)]:
            i, j = key
            terms = [(POLE[PAIR_TO_POLE[(int(p0), int(p1))]],
                      mp.mpf(int(c0)) + r2*mp.mpf(int(c1)) if isinstance(c0, int) else
                      mp.mpf(c0.numerator)/c0.denominator + r2*mp.mpf(c1.numerator)/c1.denominator)
                     for (p0, p1), (c0, c1, c2, c3) in K[key]['real'].items()]
            def f(v, terms=terms, i=i, j=j):
                return l1(v)**i * l2(v)**j * mp.fsum(c/(v - r) for r, c in terms)
            val = mp.quad(f, [bN, mp.mpf(1)/2, 1, 2, aN], maxdegree=14)
            tot += val
        return tot
    from math import comb
    for kv in [6, 7, 10]:
        wdir = w6_quad_direct(f'kern_{kv}.pkl')
        wpoly = sum(mp.mpf(comb(kv - 6, m - 1)) * Evals[f"E5,{m}"] for m in range(1, 6))
        print(f"W6({kv}): direct-quad={mp.nstr(wdir, 40)}")
        print(f"          poly-sum  ={mp.nstr(wpoly, 40)}  |d|={mp.nstr(abs(wdir-wpoly),3)}", flush=True)
    with open('progress.log', 'a') as f:
        f.write(f'step_E2: E5,m values; c(6)-LEAD={mp.nstr(dc6,3)}; W6(6/7/10) direct vs poly checked\n')
    print("saved e5m_values.json")

if __name__ == '__main__':
    main()
'''
open(os.path.join(WD, 'step_08.py'), 'w').write(step)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice', '-n', '10', PY, '-B', 'step_08.py'], capture_output=True, text=True, cwd=WD, timeout=2400)
print(r.stdout)
print("STDERR:", r.stderr[-1500:])
