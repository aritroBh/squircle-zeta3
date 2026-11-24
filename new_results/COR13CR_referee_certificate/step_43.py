import time, sys, json, random
sys.dont_write_bytecode = True
base = '<project>/research/08_theory'
sys.path.insert(0, base + '/HEXC1_hexagon_certificate')
sys.path.insert(0, base + '/HEXC2R_referee_certificate')
import mpmath as mp
mp.mp.dps = 65
from fractions import Fraction as Fr
import cmzv6
import reflib as R
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

D = base + '/COR13C_third_corner_certificate'
cert23 = json.load(open(D + '/cert13_w23.json'))
cert4 = json.load(open(D + '/cert13_w4.json'))
certE = json.load(open(D + '/cert13_euler.json'))

# family tagging in first-emit order (mirrors make_relations)
def tag_families(n, lower):
    fam = {}
    def emit(row, f):
        r = cmzv6.primitive(row)
        if r and r not in fam: fam[r] = f
    if n == 1:
        for r in cmzv6.weight1_rows():
            rr = cmzv6.primitive(r)
            if rr: fam[rr] = 'weight1_branch'
    else:
        for a in range(1, n):
            b = n - a
            for u in cmzv6.words(a, False):
                u_adm = u[0] != 0
                if u_adm and a > b: continue
                for v in cmzv6.words(b):
                    if u_adm and a == b and u > v: continue
                    emit(cmzv6.ds(u, v), 'double_shuffle')
    for d in (2, 3, 6):
        for row in cmzv6.distribution(n, d): emit(row, f'distribution_d{d}')
    for k, basis in lower.items():
        for row in basis:
            for w in cmzv6.words(n - k):
                emit(cmzv6.mul(dict(row), {w: 1}), f'lift_w{k}')
    return fam

cert_rows = {2: cert23["2"]["prow"], 3: cert23["3"]["prow"], 4: cert4["prow"]}
random.seed(2718)
overall_worst = 0
lower = {}
family_usage = {}
for n in (1, 2, 3, 4):
    rows, counts = cmzv6.make_relations(n, lower)
    lower[n] = rows
    if n == 1: continue
    fam = tag_families(n, lower)
    test_idx = list(cert_rows[n])
    if n == 2:
        test_idx = sorted(set(test_idx) | set(certE["prow"]))
    unused = [i for i in range(len(rows)) if i not in set(test_idx)]
    random.shuffle(unused)
    per_fam = {}
    for i in unused:
        f = fam[rows[i]]
        if per_fam.get(f, 0) < (10 if n < 4 else 6):
            per_fam[f] = per_fam.get(f, 0) + 1
            test_idx.append(i)
    worst = {}
    usage = {}
    for i in test_idx:
        row = rows[i]
        for w, c in row:
            assert w[0] != 0 and w[-1] != -1
        val = mp.mpc(0)
        for w, c in row:
            val += int(c) * R.cmzv_value(w)
        f = fam[row]
        a = abs(val)
        if a > worst.get(f, (0,0))[1]: worst[f] = (i, a)
        if i in set(cert_rows[n]): usage[f] = usage.get(f, 0) + 1
    ow = max(a for i, a in worst.values())
    overall_worst = max(overall_worst, ow)
    family_usage[n] = usage
    log(f"weight {n}: tested {len(test_idx)} rows ({len(cert_rows[n])} cert + unused); per-family worst:")
    for f, (i, a) in sorted(worst.items()):
        log(f"    {f:18s} {mp.nstr(a,3)} (row {i})")
    log(f"    weight {n} overall worst {mp.nstr(ow,3)}  {'PASS' if ow < mp.mpf('1e-58') else 'FAIL'}")
    assert ow < mp.mpf('1e-58')
log(f"cert family usage: {family_usage}")
total = sum(len(cert_rows[n]) for n in (2,3,4)) + len(certE["prow"])
log(f"STEP PASS: all {total} certificate rows + unused samples vanish at 65 dps, worst {mp.nstr(overall_worst,3)}")
with open('progress.log','a') as f:
    f.write("referee step(rows numeric): all %d cert rows + unused samples vanish at 65dps, worst %.1e\n" % (total, overall_worst))
