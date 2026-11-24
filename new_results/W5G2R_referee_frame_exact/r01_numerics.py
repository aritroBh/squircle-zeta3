"""r01: own evaluator. (1) weight-1 counterexample; (2) every literal row family vanishes numerically (std, lift*, 7-images);
(3) all 74 targets u.f_i vanish, controls do not; (4) frame rows used and weight-4 generators vanish."""
from rc import *
import mpmath as mp
LOG = 'r01.log'
t0 = time.time()
H = evaluate(5, q=24, K=62)
pickle.dump(H, open(HERE + '/H_values.pkl', 'wb'))
log(f'evaluated {len(H)} admissible words of weight <=5 ({time.time()-t0:.0f}s); zeta5 check {abs(H[(-1,-1,-1,-1,0)]-float(mp.zeta(5))):.2e}', LOG)
# (1) weight 1
b0 = {(1,): 2, (2,): -3, (6,): 3, (7,): -2}
z = mp.exp(2j * mp.pi / 8)
def hv(d): return sum(c * (-mp.log(1 - z ** w[0])) for w, c in d.items())
log(f'weight-1: b0 = 2(H1-H7)-3(H2-H6): value {mp.nstr(hv(b0),5)}; u=3 image {act((1,),3)}..: value {mp.nstr(hv(actv(b0,3)),12)} '
    f'(2 pi i = {mp.nstr(2j*mp.pi,12)}); u=5: {mp.nstr(hv(actv(b0,5)),12)}; u=7: {mp.nstr(hv(actv(b0,7)),5)}', LOG)
def rel(d):
    s, sc = val(d, H); return abs(s) / max(sc, 1e-300), abs(s)
def scan(name, rows):
    worst = 0; worsta = 0; n = 0
    for d in rows:
        r, a = rel(d); n += 1
        if r > worst: worst = r
        if a > worsta: worsta = a
    log(f'{name}: {n} rows, max |val|/sum|c||H| = {worst:.2e}, max |val| = {worsta:.2e} ({time.time()-t0:.0f}s)', LOG)
    return worst
for fam in ('std', 'lift1', 'lift2old', 'lift3old', 'lift4old'):
    R = load_rows(fam)
    scan(fam, (d for _, d in R))
    if fam != 'std': scan(fam + ' u=7 images', (actv(d, 7) for _, d in R))
    if fam == 'lift1':
        br = [d for lab, d in R if lab[1] in (0, 1)]
        for u in (3, 5):
            vals = [abs(val(actv(d, u), H)[0]) for d in br]
            log(f'  control: u={u} images of the {len(br)} lift1 rows on branch base rows: |val| in [{min(vals):.3f}, {max(vals):.3f}] (not relations)', LOG)
    del R
V = pickle.load(open(T + '/W5C1_certificate_scoping/validated_w5.pkl', 'rb'))
f5 = [{tuple(w): c for w, c in r} for r in V['rows']]
scan('frame-w5 rows (all 288)', f5)
OPEN = [4, 5, 7, 9, 10, 12, 14, 15, 16, 17, 52, 53, 54, 56, 58, 59, 61, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 87, 91, 95, 103, 104, 107, 113]
worst = 0
for i in OPEN:
    for u in (3, 5):
        r, a = rel(actv(f5[i], u)); worst = max(worst, r)
log(f'all 74 targets u.f_i (own action): max relative |val| = {worst:.2e}', LOG)
# control: a word permutation of a target must not vanish
import random; random.seed(1)
perm = {j: j for j in range(-1, 8)}; perm[1], perm[2] = 2, 1
ctl = {tuple(perm[x] for x in w): c for w, c in actv(f5[4], 3).items()}
log(f'control (letters 1<->2 swapped in u3.f4): relative |val| = {rel(ctl)[0]:.2e}', LOG)
for n in (2, 3, 4):
    d = pickle.load(open(T + f'/A3h4_frame_ds/validated_w{n}.pickle', 'rb'))
    scan(f'A3h4 frame rows w{n} (all {len(d["rows"])})', ({tuple(w): c for w, c in r} for r in d['rows']))
B = bases()
for n in (1, 2, 3, 4): scan(f'base B{n}', B[n])
log(f'done ({time.time()-t0:.0f}s)', LOG)
