
"""verify_town1.py — standalone re-verification (PASS/FAIL), < 15 min.

  Task A: truncated engine == recorded towers (towlib.rho / towlib.u) for k <= 12, all keys, depths <= 6.
  Task B (sample): tower law T for all genuine cells n = 2..18 at D=6 (harvest K=22), plus diagonal-shift sample.
  Task C (sample): degree-2 kernels from the tower == closed forms R_20, R_11 at k = 16..22 (beyond WT3's k<=15).

Run:  nice -n 10 $PY verify_town1.py
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
WTS = os.path.normpath(os.path.join(HERE, '..', 'WTS_weight_k_structure', 'theory'))
TOW1 = os.path.normpath(os.path.join(HERE, '..', 'TOW1_tower_law'))
sys.path.insert(0, TOW1); sys.path.insert(0, WTS)
import trunc_engine as te
import driver
import kernel_from_tower as kft
import towlib

T0 = time.time()
results = []

def check(name, ok, detail=""):
    results.append((name, bool(ok), detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}", flush=True)

print("== Task A: engine == recorded towers (k<=12, depths<=6) ==")
S2 = te.stirling_table(14)
D = 6
bad_rho = 0; tot_rho = 0
for d in range(2, 13):
    rho_mine = te.rho_level(d, D, S2)
    for key in te.keys_of_degree(d):
        for delta in range(0, min(D, 2*d-4)+1):
            tot_rho += 1
            if rho_mine[key][delta] != towlib.rho(key, delta): bad_rho += 1
check("A1: rho == towlib.rho (d=2..12)", bad_rho == 0, f"{tot_rho} coeffs, {bad_rho} mismatch")

bad_t = 0; tot_t = 0
for k in range(2, 13):
    tower_k = te.build_tower_level(k, D, te.rho_level(k, D, S2))
    for n in range(2, k+1):
        for key in te.keys_of_degree(n):
            for delta in range(0, min(D, 2*n-4)+1):
                tot_t += 1
                if tower_k[key][delta] != towlib.u(key[0], key[1], k, delta): bad_t += 1
check("A2: tower == towlib.u (k=2..12, all keys)", bad_t == 0, f"{tot_t} coeffs, {bad_t} mismatch")

print("== Task B (sample): tower law T, n=2..18, D=6 (harvest K=22) ==")
Ks = 22
H = driver.harvest(D, Ks, os.path.join(HERE, 'data', f'_verify_towers_D{D}_K{Ks}.pkl'))
cells = driver.test_tower_law(H, D, 18, Ks)
gen = [c for c in cells if c['genuine']]
fails = [c for c in gen if (c['is_poly'] is not True) or (not c['deg_exact'])]
check("B1: T holds, all genuine cells n=2..18, D=6", len(fails) == 0,
      f"{len(gen)} cells, {len(fails)} fails, min nchecks={min(c['nchecks'] for c in gen)}")
# diagonal shift sample
from math import comb
def binom_poly(x, j):
    r = Fr(1)
    for i in range(j): r *= Fr(x-i)
    for i in range(1, j+1): r /= i
    return r
def poly_eval(ks, ys, k):
    k0 = ks[0]; cur = list(ys); cols=[cur]
    for j in range(1, len(ks)):
        cur=[cur[t+1]-cur[t] for t in range(len(cur)-1)]; cols.append(cur)
    x = k-k0; out = cols[0][0]
    for j in range(1, len(ks)): out = out + cols[j][0].scale(binom_poly(x, j))
    return out
unq4 = driver.unq4
nshift = 0; ntot = 0
for p in range(1, 10):
    n = 2*p; key = (p,p)
    for delta in [0,1,2,3,4]:
        if delta > 2*n-4: continue
        kd = H[key].get(delta, {})
        ksd = sorted(k for k in kd if k >= n+1); P = delta//2
        if len(ksd) < P+1: continue
        ys = [unq4(kd[k]) for k in ksd[:P+1]]
        ntot += 1
        if poly_eval(ksd[:P+1], ys, n) != unq4(kd[n]): nshift += 1
check("B2: diagonal shift (k=n off polynomial)", nshift == ntot and ntot > 0, f"{nshift}/{ntot} off")

print("== Task C (sample): degree-2 kernels == closed forms R_20, R_11, k=16..22 ==")
def h_orders(key, k):
    d = key[0]+key[1]
    return {2*d-2-delta: unq4(kd[k]) for delta, kd in H[key].items() if k in kd}
badc = 0
for k in range(16, 23):
    K20, K11, _, _ = kft.kernel_20_21_level(h_orders((3,0), k), h_orders((2,1), k))
    if not kft.pf_equal(K20, kft.R20_closed(k)): badc += 1
    if not kft.pf_equal(K11, kft.R11_closed(k)): badc += 1
check("C1: R_20(k;v), R_11(k;v) exact from tower, k=16..22", badc == 0, f"{badc} mismatches")

npass = sum(1 for _, ok, _ in results if ok)
print(f"\nVERIFY_TOWN1: {npass}/{len(results)} PASS   ({time.time()-T0:.1f}s)")
print("ALL PASS" if npass == len(results) else "SOME FAILURES")
