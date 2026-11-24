# s13: the T1/Newton system cannot see a_6 -> a_6 + c/Pi_3(d); the data / first-principles rho check does.
import sys, time
sys.dont_write_bytecode = True
from candlib import *
d, m = T3.d, T3.m; Pi = T3.Pi
base = base_exprs()
ex = dict(base); ex[(6,0)] = base[(6,0)] + 3/Pi(3)
cand = Cand(ex)
H = load_data()
tested = 0; bad = []
for (p, q) in sorted([k for k in H if k != '_meta']):
    dd = p+q
    if dd < 5 or dd > 40: continue
    for k in range(dd, min(dd+4, 52)+1):
        if p == q and k < 2*p: continue
        if not data_has(p, q, k, 6): continue
        tested += 1
        if cand.get(p, q, k, 6) != data_get(p, q, k, 6): bad.append((p, q, k))
print(f'mutant a_6+3/Pi_3(d) vs TOWN1 data (delta=6 cells, d=5..40): {tested} cells, {len(bad)} bad (T1 test: 0 bad)')
