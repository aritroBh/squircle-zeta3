
import mpmath as mp
mp.mp.dps = 55
import json, sys

exec(open("av5b_locmodel.py").read().split("def extract")[0])

def dP_loc(mm, d):
    Zo = Loc(mm)
    p = 1+d
    W = mp.pi/(2*mm)
    tau = d/(1+d)
    f = lambda u: Zo.h(u, p) - Zo.tent(u)
    tot = mp.mpf(0)
    for a, b in [(-W, -20*tau), (-20*tau, mp.mpf(0)), (mp.mpf(0), 20*tau), (20*tau, W)]:
        tot += mp.quad(f, [a, b])
    return 2*mm*tot

if __name__ == "__main__":
    mm = int(sys.argv[1])
    out = []
    for dstr in sys.argv[2].split(","):
        d = mp.mpf(dstr)
        v = dP_loc(mm, d)
        out.append((dstr, v))
        print(f"m={mm}loc d={dstr}: dP={mp.nstr(v,35)}", flush=True)
    json.dump([[d, str(v)] for d, v in out], open(f"ckpt_dPloc_m{mm}.json","w"))
