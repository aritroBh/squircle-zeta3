# Planner: rerun ZQ5's kernel-level LLL for one D with ALL values recomputed by the planner at 1300 digits
# (Li5 from planner_check_vals.py's pickle; logs, pi^2, zeta(5), L(chi_D,5) recomputed here). Same atom set as
# pipe.phase_lll(use_atoms=False): T = L/sqrt(D), zeta(5), saturated kernel combos, weight-5 log/pi products.
# Usage: planner_rerun_lll.py D KERNELFILE P1 [P2 ...]
import sys, os, json, pickle, time
ZQ = '<project>/research/08_theory/ZQ5_zagier_weight5'
sys.path.insert(0, ZQ)
import mpmath as mp
import w5lib as W
import flint
D, KF, PS = int(sys.argv[1]), sys.argv[2], [int(p) for p in sys.argv[3:]]
C = W.letters_info(D)
K = json.load(open(f'{ZQ}/out/d{D}/{KF}'))['K']
li5raw = pickle.load(open(f'{ZQ}/planner/d{D}_li5_planner.pkl', 'rb'))
mp.mp.dps = 1300
li5v = [mp.mpf(li5raw[i]) for i in range(len(li5raw))]
logs = [mp.log(abs(g.emb(1300))) for g in C['letters']]
L = mp.mpf(0)
for a in range(1, D):
    c = W.chi(D, a)
    if c: L += c * mp.zeta(5, mp.mpf(a) / D)
L /= mp.mpf(D)**5
T = L / mp.sqrt(D); z5 = mp.zeta(5); pi2 = mp.pi**2
names, vals = ['T', 'z5'], [T, z5]
zero_combos = []
for j, row in enumerate(K):
    v = sum(mp.mpf(c) * li5v[i] for i, c in enumerate(row) if c)
    if abs(v) < mp.mpf(10)**(-1200): zero_combos.append(j)
    else: names.append(f'K{j}'); vals.append(v)
for kind, ms, nm in W.weight5_products(C['names']):
    v = mp.mpf(1)
    if kind == 'pi2log3': v = pi2
    elif kind != 'log5': v = pi2**2
    for i in ms: v *= logs[i]
    names.append(nm); vals.append(v)
print(f'D={D}: letters {C["names"]}, kernel {len(K)} (exact-zero combos dropped: {zero_combos}), dim {len(vals)}', flush=True)
for P in PS:
    mp.mp.dps = P + 40
    n = len(vals); Cc = mp.power(10, P)
    B = []
    for i, v in enumerate(vals):
        row = [0] * (n + 1); row[i] = 1; row[n] = int(mp.nint(v * Cc)); B.append(row)
    t0 = time.time(); R = flint.fmpz_mat(B).lll(gram='exact'); dt = time.time() - t0
    mp.mp.dps = 1300
    rels, zeros = [], []
    for i in range(n):
        a = [int(R[i, j]) for j in range(n)]
        if not any(a): continue
        res = sum(mp.mpf(x) * v for x, v in zip(a, vals) if x)
        if abs(res) < mp.mpf(10)**(-1200):
            (rels if a[0] else zeros).append((max(map(abs, a)), a, res))
    rels.sort(key=lambda t: t[0])
    print(f'P={P}: LLL {dt:.0f}s; exact (residual < 1e-1200) relations with T: {len(rels)}; exact zero-relations: {len(zeros)}', flush=True)
    for mx, a, res in rels[:3]:
        nz = [(names[k], a[k]) for k in range(n) if a[k]]
        print(f'   maxcoef {mx}  nterms {len(nz)}  residual {mp.nstr(res, 3)}  T-coef {a[0]}  z5-coef {a[1]}', flush=True)
    pickle.dump({'D': D, 'P': P, 'names': names, 'rels': [(mx, a, mp.nstr(r, 5)) for mx, a, r in rels],
                 'zeros': [(mx, a, mp.nstr(r, 5)) for mx, a, r in zeros]}, open(f'{ZQ}/planner/d{D}_rerun_{P}.pkl', 'wb'))
