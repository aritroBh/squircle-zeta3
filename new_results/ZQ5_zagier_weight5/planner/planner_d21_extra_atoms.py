# Planner experiment: does D=21 close once zero-symbol weight-5 products with zeta(3) and L(chi21,3) are added?
# Read-only use of ZQ5's pipe/w5lib and stored values. Usage: planner_d21_extra_atoms.py P1 P2
import sys, os, json, time, itertools
ZQ = '<project>/research/08_theory/ZQ5_zagier_weight5'
sys.path.insert(0, ZQ); os.chdir(ZQ)
import mpmath as mp
import w5lib as W, pipe
import flint
D = 21
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
K = json.load(open(pipe.cand_dir(D) + '/kernel_b10.json'))['K']
raw = pipe.load_vals.__globals__['pickle'].load(open(pipe.cand_dir(D) + '/vals_b10.pkl', 'rb'))
mp.mp.dps = 1250
vals = {k: mp.mpf(v) for k, v in raw.items()}
names = C['names']; d = len(names)
z3 = mp.zeta(3)
L3 = mp.mpf(0)
for a in range(1, D):
    c = W.chi(D, a)
    if c: L3 += c * mp.zeta(3, mp.mpf(a) / D)
L3 /= mp.mpf(D)**3
T3 = L3 / mp.sqrt(D)
logs = [vals['log:' + n] for n in names]
extra = []
for w3n, w3 in (('z3', z3), ('T3', T3)):
    extra.append((w3n + '*pi2', w3 * vals['pi2']))
    for i, j in itertools.combinations_with_replacement(range(d), 2):
        extra.append((f'{w3n}*{names[i]}*{names[j]}', w3 * logs[i] * logs[j]))
li5v = [vals[f'li5:{i}'] for i in range(len(cands))]
combos = [(f'K{j}', sum(mp.mpf(c) * v for c, v in zip(row, li5v))) for j, row in enumerate(K)]
prods = pipe.prod_values(C, vals)
base_names = ['T', 'z5'] + [n for n, _ in combos] + [p[2] for p in prods] + [n for n, _ in extra]
base_vals = [vals['T'], vals['z5']] + [v for _, v in combos] + [p[3] for p in prods] + [v for _, v in extra]
print(f'D={D}: letters {names}, candidates {len(cands)}, kernel {len(K)}, products {len(prods)}, extra {len(extra)}, dim {len(base_vals)}', flush=True)
def run(P):
    mp.mp.dps = P + 40
    n = len(base_vals); Cc = mp.power(10, P)
    B = []
    for i, v in enumerate(base_vals):
        row = [0] * (n + 1); row[i] = 1; row[n] = int(mp.nint(v * Cc)); B.append(row)
    t0 = time.time(); R = flint.fmpz_mat(B).lll(gram='exact'); dt = time.time() - t0
    mp.mp.dps = 1250
    good = []
    for i in range(n):
        a = [int(R[i, j]) for j in range(n)]
        if not any(a): continue
        res = sum(mp.mpf(x) * v for x, v in zip(a, base_vals))
        if abs(res) < mp.mpf(10)**(-(P + 60)):
            good.append((a, res))
    withT = [g for g in good if g[0][0] != 0]
    print(f'P={P}: LLL {dt:.0f}s; rows with residual < 1e-{P+60}: {len(good)} ; with T coefficient != 0: {len(withT)}', flush=True)
    for a, res in sorted(withT, key=lambda g: max(map(abs, g[0])))[:3]:
        used = [(base_names[k], a[k]) for k in range(n) if a[k]]
        print('   maxcoef', max(map(abs, a)), 'residual', mp.nstr(res, 3), 'nterms', len(used))
        print('   uses extra atoms:', [u for u in used if u[0].startswith(('z3', 'T3'))])
        print('   T coef', a[0])
    return withT, good
r1 = run(int(sys.argv[1])); r2 = run(int(sys.argv[2]))
