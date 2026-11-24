
import mpmath as mp, flint, time, sys, json
from mpmath.calculus.quadrature import TanhSinh

K   = 400
MDPS = 160
BITS = 800
UMAX = mp.mpf(1500)

def taylor_coeffs(udeg):
    mp.mp.dps = MDPS
    TS = TanhSinh(mp.mp)
    nodes = []
    for d in range(1, udeg+1):
        nodes.extend(TS.get_nodes(mp.mpf(0), mp.inf, d, mp.mp.prec))
    uh = mp.mpf(2)**(-udeg)
    flint.ctx.prec = BITS
    flint.ctx.cap = K+1
    s = flint.acb_series([0,1], prec=K+1)
    sums = [flint.acb(0) for _ in range(K+1)]
    used = 0
    t0 = time.time()
    for um, w0 in nodes:
        if um > UMAX:
            continue
        used += 1
        w = w0/(4*mp.cosh(um/2)**2)
        u = flint.acb(mp.nstr(um, MDPS+5))
        x = (-u).exp(); L = (1+x).log()
        A = (1+x)*((-L)*s).exp()
        Z = x*(u*s).exp()
        B = 1 - Z/(1+(1+Z*Z).sqrt())
        g = A*B
        ww = flint.acb(mp.nstr(w*uh, MDPS+5))
        cf = g.coeffs()
        if len(cf) < K+1:
            cf += [flint.acb(0)]*(K+1-len(cf))
        for n in range(K+1):
            sums[n] -= ww*cf[n]
    return sums, used, time.time()-t0

if __name__ == "__main__":
    udeg = int(sys.argv[1])
    sums, used, dt = taylor_coeffs(udeg)
    print("udeg", udeg, "nodes used", used, "elapsed %.1f s" % dt)
    with open(f"mk_deg{udeg}.tsv","w") as f:
        for k,z in enumerate(sums):
            line = f"{k}\t{z.mid().str(80)}\t{z.rad().str(6)}"
            f.write(line+"\n")
    # quick controls
    print("M0   =", sums[0].mid().str(45))
    for k in [50,100,150,200,250,300,350,400]:
        print(k, sums[k].mid().str(20), " rad", sums[k].rad().str(3))
