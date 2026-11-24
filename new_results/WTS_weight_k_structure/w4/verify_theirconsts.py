"""verify_theirconsts.py <dps_of_their_file>   Compare the author's stored 8 constants results/id_const_dps<dps>.json (read only) with MY
800-digit combinations c_0,c_1,c_2,delta (verify_const.py), using MY exact coefficient table (verify_lib.exact_table / k3_table).
PRECISION: mp.mp.dps set after all imports and printed."""
import sys, os, json
sys.dont_write_bytecode = True
from verify_lib import *
theirs = int(sys.argv[1]); mine = 800
set_dps(mine)
print("mp.mp.dps =", mp.mp.dps)
OUT = os.environ["VERIFY_OUT"]
d = json.load(open(os.path.join(W4, "results", f"id_const_dps{theirs}.json")))["I"]
I = {}
for k, v in d.items():
    ij, pole = k.split("|")
    key = (int(ij[0]), int(ij[1]))
    I[(key, {"-a": "ma", "0": "0", "-b": "mb", "a": "a"}[pole])] = mp.mpf(v)
TAB = exact_table(); K3 = k3_table(); r2 = mp.sqrt(2)
def comb(coefs): return sum((fr_mp(p0) + fr_mp(p1) * r2) * I[kp] for kp, (p0, p1) in coefs.items())
def num(w): return sum(mp.mpf(json.load(open(os.path.join(OUT, f"const_{w}_h{h}_dps{mine}.json")))["value"]) for h in (0, 1))
for j in range(3):
    c = comb({kp: v[j] for kp, v in TAB.items()})
    print(f"c_{j}: author's stored constants (dps {theirs}) combined with MY exact table - my dps-800 quadrature = {mp.nstr(c - num(f'c{j}'), 5)}  (mp.dps {mp.mp.dps})")
P3 = comb({kp: tuple(v[0][i] + 3 * v[1][i] + 9 * v[2][i] for i in (0, 1)) for kp, v in TAB.items()})
w43 = comb(K3)
print(f"delta: (P(3) - W_4(3)) from author's stored constants - mine = {mp.nstr(P3 - w43 - num('delta'), 5)}  (mp.dps {mp.mp.dps})")
