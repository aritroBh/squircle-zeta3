# Planner: rerun ZQ5's D=21 kernel LLL (same atoms as pipe.phase_lll, use_atoms=False) after replacing the
# stored Li5 values that are wrong (indices given on the command line) by direct-series values at 1280 digits.
# Usage: planner_d21_rerun_fixed.py "60,76" P1 P2
import sys, os, json, pickle, time
ZQ = '<project>/research/08_theory/ZQ5_zagier_weight5'
sys.path.insert(0, ZQ); os.chdir(ZQ)
from fractions import Fraction as F
import mpmath as mp
import w5lib as W, pipe
import flint
D = 21
bad = [int(t) for t in sys.argv[1].split(',') if t]
P1, P2 = int(sys.argv[2]), int(sys.argv[3])
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
K = json.load(open(pipe.cand_dir(D) + '/kernel_b10.json'))['K']
raw = pickle.load(open(pipe.cand_dir(D) + '/vals_b10.pkl', 'rb'))
mp.mp.dps = 1280
vals = {k: mp.mpf(v) for k, v in raw.items()}
s = mp.sqrt(D)
Jraw = json.load(open(pipe.cand_dir(D) + '/candidates.json'))
fixfile = os.path.join(ZQ, 'planner', 'd21_fixed_li5.pkl')
fixed = pickle.load(open(fixfile, 'rb')) if os.path.exists(fixfile) else {}
for i in bad:
    if i in fixed: continue
    a, b = (F(t) for t in Jraw['cands'][i]['x'])
    x = mp.mpf(a.numerator) / a.denominator + mp.mpf(b.numerator) / b.denominator * s
    tot = mp.mpf(0); p = x; n = 1; eps = mp.mpf(10)**(-1265)
    while abs(p) > eps:
        tot += p / mp.mpf(n)**5; n += 1; p *= x
    fixed[i] = mp.nstr(tot, 1270)
    pickle.dump(fixed, open(fixfile, 'wb'))
    print(f'fixed li5:{i} change {mp.nstr(abs(tot - vals[f"li5:{i}"]), 3)}', flush=True)
for i, v in fixed.items():
    vals[f'li5:{i}'] = mp.mpf(v)
for P in (P1, P2):
    out = pipe.phase_lll(D, C, cands, K, vals, P=P, use_atoms=False)
    print(f'P={P}: exact relations with T: {len(out["exact_rels"])}; exact zero-relations: {len(out["exact_zeros"])}', flush=True)
    for mx, a, res in out['exact_rels'][:3]:
        nz = [(out['names'][k], a[k]) for k in range(len(a)) if a[k]]
        print(f'   maxcoef {mx} nterms {len(nz)} residual {mp.nstr(res, 3)} T-coef {a[0]}')
        print('   ', nz)
    pickle.dump({'P': P, 'exact_rels': [(mx, a, mp.nstr(r, 5)) for mx, a, r in out['exact_rels']], 'names': out['names']},
                open(os.path.join(ZQ, 'planner', f'd21_rerun_fixed_{P}.pkl'), 'wb'))
