"""High-precision values: usage step_hp.py <target> <workdps>. Appends to hp_results.json."""
import sys, json, time, os
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp

target, DPS = sys.argv[1], int(sys.argv[2])
mp.mp.dps = DPS
r2 = mp.sqrt(2); aN = 1+r2; bN = r2-1
omega = json.load(open("omega4.json"))
terms = {tuple(map(int, k.split(","))): [(z, mp.mpf(c[0]) + r2*mp.mpf(c[1])) for z, c in v]
         for k, v in omega.items()}
pole = {"0": mp.mpf(0), "a": aN, "-a": -aN, "b": bN, "-b": -bN}
def piece(i, j):
    ts = terms[(i, j)]
    def f(v):
        l1 = mp.log((v - bN)*(v + aN)/(4*v))
        l2 = mp.log((aN - v)*(v + bN)/(4*v))
        return l1**i * l2**j * mp.fsum(c/(v - pole[z]) for z, c in ts)
    val, err = mp.quad(f, [bN, mp.mpf(1)/2, 1, 2, aN], error=True, maxdegree=14)
    return val, err

T0 = time.time()
if target == "w4pp":
    v1, e1 = piece(2, 1); print(" (2,1) done", time.time()-T0, flush=True)
    v2, e2 = piece(3, 0)
    val, err = v1+v2, max(e1, e2)
elif target == "w5":
    val = mp.mpf(0); err = mp.mpf(0)
    for ij in [(2, 2), (3, 1), (4, 0)]:
        v_, e_ = piece(*ij); val += v_; err = max(err, e_)
        print(f" {ij} done {time.time()-T0:.0f}s", flush=True)
elif target == "m4":
    import sympy as sp
    from mk_integrands import v_integrand, v as vs
    Rk = v_integrand(4)
    fs = [(i, j, sp.lambdify(vs, R, modules="mpmath")) for (i, j), R in Rk.items()]
    def fM4(x):
        l1 = mp.log((x - bN)*(x + aN)/(4*x))
        l2 = mp.log((aN - x)*(x + bN)/(4*x))
        return mp.fsum(R(x)*l1**i*l2**j for i, j, R in fs)
    val, err = mp.quad(fM4, [bN, mp.mpf(1)/2, 1, 2, aN], error=True, maxdegree=14)
nd = DPS - 40
s = mp.nstr(val, nd)
res = {}
if os.path.exists("hp_results.json"):
    res = json.load(open("hp_results.json"))
res[f"{target}@{DPS}"] = {"value": s, "digits": nd, "err_est": mp.nstr(err, 3), "seconds": time.time()-T0}
json.dump(res, open("hp_results.json", "w"), indent=1)
print(f"{target} @ working {DPS} dps: err_est {mp.nstr(err,3)}; {time.time()-T0:.0f}s; saved {nd} digits")
