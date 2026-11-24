import json, time
import mpmath as mp

C = json.load(open("consts360.json"))
labels = ["X"] + C["labels"]
strs = [C["vals"]["X"]] + [C["vals"][l] for l in C["labels"]]

def run_pslq(p, tol_exp=25, maxcoeff=int(1e14)):
    mp.mp.dps = p
    v = [mp.mpf(s) for s in strs]
    tol = mp.mpf(10)**(-(p - tol_exp))
    t0 = time.time()
    rel = mp.pslq(v, tol=tol, maxcoeff=maxcoeff, maxsteps=200000)
    dt = time.time() - t0
    # residual at this precision
    if rel is not None:
        res = sum(mp.mpf(int(c))*vv for c, vv in zip(rel, v))
    else:
        res = None
    return rel, res, dt

results = {}
for p in (250, 280, 310):
    rel, res, dt = run_pslq(p)
    results[p] = rel
    print(f"=== PSLQ at {p} digits, tol=1e-{p-25}, maxcoeff=1e14  ({dt:.1f}s)")
    if rel is None:
        print("   NO relation found (None)")
    else:
        print("   maxcoeff:", max(abs(int(c)) for c in rel))
        print("   residual at working precision:", mp.nstr(res, 5))
        print("   rel =", [int(c) for c in rel])

same = results[250] is not None and results[250] == results[280] == results[310]
print("\nIDENTICAL at 250/280/310:", same)
json.dump({str(k): ([int(c) for c in v] if v else None) for k, v in results.items()},
          open("pslq_runs.json", "w"))
