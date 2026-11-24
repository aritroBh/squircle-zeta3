"""WT4 step D: high-precision values of D1/D2 by tanh-sinh quadrature, two working dps.
Usage: wt4_D_hp.py <target: D1|D2> <workdps> [<workdps2> ...]"""
import sys, os, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
import mpmath as mp

target = sys.argv[1]
DPS_LIST = [int(x) for x in sys.argv[2:]]
deltas = json.load(open(os.path.join(HERE, "deltas_deg3.json")))
kerns = {k: deltas[k][target] for k in ["2,1", "3,0"]}

resfile = os.path.join(HERE, "hp_results.json")
res = json.load(open(resfile)) if os.path.exists(resfile) else {}

for DPS in DPS_LIST:
    mp.mp.dps = DPS
    T0 = time.time()
    r2 = mp.sqrt(2); aN = 1+r2; bN = r2-1
    pole = {"0": mp.mpf(0), "a": aN, "-a": -aN, "b": bN, "-b": -bN}
    tot = mp.mpf(0); err = mp.mpf(0)
    for key, i, j in [("2,1", 2, 1), ("3,0", 3, 0)]:
        ts = [(z, mp.mpf(c[0]) + r2*mp.mpf(c[1])) for z, c in kerns[key]]
        def f(v, ts=ts, i=i, j=j):
            l1 = mp.log((v - bN)*(v + aN)/(4*v))
            l2 = mp.log((aN - v)*(v + bN)/(4*v))
            return l1**i * l2**j * mp.fsum(c/(v - pole[z]) for z, c in ts)
        val, e = mp.quad(f, [bN, mp.mpf(1)/2, 1, 2, aN], error=True, maxdegree=14)
        tot += val; err = max(err, e)
        print(f"  {key} done {time.time()-T0:.0f}s (mp.mp.dps={mp.mp.dps})", flush=True)
    nd = DPS - 40
    s = mp.nstr(tot, nd)
    res[f"{target}@{DPS}"] = {"value": s, "digits": nd, "err_est": mp.nstr(err, 3),
                              "seconds": time.time()-T0}
    json.dump(res, open(resfile, "w"), indent=1)
    open(os.path.join(HERE, f"{target.lower()}_{nd}.txt"), "w").write(s + "\n")
    print(f"{target} @ working {DPS} dps: err_est {mp.nstr(err,3)}; {time.time()-T0:.0f}s; saved {nd} digits (mp.mp.dps={mp.mp.dps})", flush=True)
with open(os.path.join(HERE, "progress.log"), "a") as f:
    for DPS in DPS_LIST:
        f.write(f"step_D hp: {target}@{DPS} err_est={res[f'{target}@{DPS}']['err_est']} {res[f'{target}@{DPS}']['seconds']:.0f}s\n")
print("DONE")
