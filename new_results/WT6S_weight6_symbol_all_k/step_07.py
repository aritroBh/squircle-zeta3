"""WT6S step E (part 1): N_{(i,j),r} to ~150 digits, two independent quadrature methods.

N_{(i,j),r} = int_b^a l1^i l2^j dv/(v-r), l1=log t, l2=log(1-t), t=(v^2+2v-1)/(4v).
Method RAW: tanh-sinh over [b,1/2,1,2,a] (M5a setup).
Method SUB: split [b,1]+[1,a], power substitution v=b+(1-b)u^p / v=a-(a-1)u^p (p=12)
  removing the endpoint log singularities; exact v-b=(1-b)u^p, a-v=(a-1)u^p; for r=a on
  [1,a] use dv/(v-a) = -(P/u)du. tanh-sinh on [0,1].
Cross-check agreement; store n6_values.json."""
import json, time
import mpmath as mp

mp.mp.dps = 170
r2 = mp.sqrt(2); aN = 1 + r2; bN = r2 - 1
POLE = {"0": mp.mpf(0), "a": aN, "-a": -aN, "-b": -bN}
KEYS = [(5, 0), (4, 1), (3, 2)]
PRS = ["0", "a", "-a", "-b"]
P = 12

def N_raw(i, j, rname):
    r = POLE[rname]
    def f(v):
        l1 = mp.log((v - bN) * (v + aN) / (4 * v))
        l2 = mp.log((aN - v) * (v + bN) / (4 * v))
        return l1**i * l2**j / (v - r)
    return mp.quad(f, [bN, mp.mpf(1)/2, 1, 2, aN], error=True, maxdegree=14)

def N_sub(i, j, rname):
    r = POLE[rname]
    def f1(u):          # v in [b,1], v-b = (1-b)u^p exact
        if u == 0:
            return mp.mpf(0)
        vb = (1 - bN) * u**P
        v = bN + vb
        l1 = mp.log(vb * (v + aN) / (4 * v))
        l2 = mp.log((aN - v) * (v + bN) / (4 * v))
        return l1**i * l2**j / (v - r) * ((1 - bN) * P * u**(P - 1))
    def f2(u):          # v in [1,a], a-v = (a-1)u^p exact
        if u == 0:
            return mp.mpf(0)
        da = (aN - 1) * u**P
        v = aN - da
        l1 = mp.log((v - bN) * (v + aN) / (4 * v))
        l2 = mp.log(da * (v + bN) / (4 * v))
        if r == aN:     # dv/(v-a) = -(P/u) du
            return -l1**i * l2**j * (P / u)
        return l1**i * l2**j / (v - r) * ((aN - 1) * P * u**(P - 1))
    v1, e1 = mp.quad(f1, [0, 1], error=True, maxdegree=14)
    v2, e2 = mp.quad(f2, [0, 1], error=True, maxdegree=14)
    return v1 + v2, max(e1, e2)

def main():
    t0 = time.time()
    vals = {}
    worst = mp.mpf(0)
    for (i, j) in KEYS:
        for rname in PRS:
            vr, er = N_raw(i, j, rname)
            vs, es = N_sub(i, j, rname)
            agree = abs(vr - vs)
            worst = max(worst, agree)
            vals[f"{i},{j}|{rname}"] = mp.nstr(vs, 150)
            print(f"N_{{({i},{j}),{rname}}}: |sub-raw|={mp.nstr(agree,3)} "
                  f"err_r={mp.nstr(er,2)} err_s={mp.nstr(es,2)} ({time.time()-t0:.0f}s)", flush=True)
    print(f"WORST |sub-raw| over 12 integrals: {mp.nstr(worst, 3)}")
    json.dump({'dps': 150, 'values': vals}, open('n6_values.json', 'w'), indent=1)
    with open('progress.log', 'a') as f:
        f.write(f'step_E1: 12 N at 150 digits; worst raw-vs-sub {mp.nstr(worst,3)}\n')
    print("saved n6_values.json")

if __name__ == '__main__':
    main()
