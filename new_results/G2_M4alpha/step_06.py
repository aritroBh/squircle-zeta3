"""step_06.py -- NUMERICAL check of the exact flat-corner series M4(c) = sum_{n>=1} mu_{4,n} c^{2n}
(mu_{4,n} from flat_exact_N400_K4.json, exact in Q[z2,z3,z4]) against
 (a) the squircle value M4 = 0.0060279931126... (M4a_decomposition/M4_800.txt, 800 digits) at c = 1/sqrt(2);
 (b) corner_series.py direct quadrature at alpha = 2pi/3 (c = 1/2) and alpha = 1.3 (c = cos(0.65));
 (c) computes and stores M4 at the six angles used for the B4 validation (values4_*.txt, ~78 digits).
Tail beyond n = 400: |c|^{802} times a polynomial factor; estimated per angle below.
"""
import sys, json
import mpmath as mp

sys.path.insert(0, "<project>/research/10_frontier/G1_zeta3_any_angle")
import corner_series as cs

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 90
mp.mp.dps = DPS + 10
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"

tab = json.load(open('flat_exact_N400_K4.json'))['4']
KEY = {'(0, 0, 0)': '1', '(1, 0, 0)': 'z2', '(0, 1, 0)': 'z3', '(0, 0, 1)': 'z4', '(2, 0, 0)': 'z2^2'}
mu = []
for nv in range(1, 401):
    d = tab[str(nv)]
    mu.append({k: mp.mpf(v) for k, v in ((KEY[kk], d[kk]) for kk in d)})


def series_val(c, N=400):
    z2 = mp.pi**2/6; z3 = mp.zeta(3); z4 = mp.pi**4/90
    vals = {'1': mp.mpf(1), 'z2': z2, 'z3': z3, 'z4': z4, 'z2^2': z2**2}
    c2 = c*c
    cp = c2
    tot = mp.mpf(0)
    for nv in range(1, N + 1):
        for kk, v in mu[nv - 1].items():
            tot += v*vals[kk]*cp
        cp *= c2
    return tot


def stored_m4():
    return mp.mpf(open("<project>/research/08_theory/M4a_decomposition/M4_800.txt").read().strip())


if __name__ == '__main__':
    # (a) square corner
    c = 1/mp.sqrt(2)
    sv = series_val(c)
    m4 = stored_m4()
    print(f"(a) alpha=pi/2: series(n<=400) = {mp.nstr(sv, 50)}")
    print(f"    |series - M4_800| = {mp.nstr(abs(sv - m4), 3)}   (tail ~ 2^-400 * poly ~ 1e-123)")
    # (b) 2pi/3 and 1.3 vs corner_series
    for tag, a in [('2pi_3', 2*mp.pi/3), ('gen1p3', mp.mpf(13)/10)]:
        c = mp.cos(a/2)
        sv = series_val(c)
        qv = cs.corner_coeffs(a, 4)[4]
        print(f"(b) alpha={tag}: |series - quadrature| = {mp.nstr(abs(sv - qv), 3)}  (c={mp.nstr(c, 6)})", flush=True)
    # (c) store M4 at the six B4-validation angles
    ANG = {'pi_6': mp.pi/6, 'pi_4': mp.pi/4, 'pi_3': mp.pi/3, 'pi_2': mp.pi/2, '2pi_3': 2*mp.pi/3, 'gen1p3': mp.mpf(13)/10}
    for tag, a in ANG.items():
        v = cs.corner_coeffs(a, 4)[4]
        with open(f'values4_{tag}.txt', 'w') as f:
            f.write(f'# alpha = {tag}; corner_series M4 at dps {mp.mp.dps}' + chr(10))
            f.write(f'M4 = {mp.nstr(v, DPS)}' + chr(10))
        c = mp.cos(a/2)
        print(f"(c) {tag:>7}: M4 = {mp.nstr(v, 40)}   |series - quadr| = {mp.nstr(abs(series_val(c) - v), 3)}", flush=True)
