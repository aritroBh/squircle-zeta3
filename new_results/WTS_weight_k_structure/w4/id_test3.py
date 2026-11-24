import sys, json, time
from id_lib import *
dps = int(sys.argv[1]); D = int(sys.argv[2]); set_dps(dps); print("mp.mp.dps =", mp.mp.dps)
I = load_I(dps)
c = [combo(I, f"c{j}") for j in range(3)]
A = atoms()
S = with_sqrt2(A)
labels = ["target"] + list(S.keys())
for nm, t in (("c1", c[1]),):
    vals = [t] + list(S.values())
    t0 = time.time()
    res = find_relation(vals, D)
    print(nm, "n =", len(vals), "D =", D, "time %.1fs" % (time.time() - t0))
    for co, s in res[:3]:
        print("  max|coef| =", max(abs(x) for x in co), " residual =", mp.nstr(s, 4), " target coef =", co[0], " (mp.dps", mp.mp.dps, ")")
    co, s = res[0]
    print({l: x for l, x in zip(labels, co) if x})
