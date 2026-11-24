"""Print the full finite-difference table (all entries) for weight w, k=kmin..kmax at given dps."""
import sys, json
from w4_common import *
dps, w, kmin, kmax = (int(x) for x in sys.argv[1:5])
set_dps(dps); print("mp.mp.dps =", mp.mp.dps)
ks = list(range(kmin, kmax + 1))
v = [mp.mpf(json.load(open(os.path.join(RES, f"w{w}_k{k}_dps{dps}.json")))["value"]) for k in ks]
row = v
for n in range(len(v)):
    print(f"order {n}: " + "  ".join(mp.nstr(x, 12) for x in row))
    row = [row[i + 1] - row[i] for i in range(len(row) - 1)]
