
import mpmath as mp, flint, time, json, os

def taylor_coeffs_TR(alpha, K=300, udeg=6, mdps=100, bits=640, umax=1200):
    from mpmath.calculus.quadrature import TanhSinh
    mp.mp.dps = mdps
    TS = TanhSinh(mp.mp); nodes=[]
    for d in range(1, udeg+1): nodes.extend(TS.get_nodes(mp.mpf(0), mp.inf, d, mp.mp.prec))
    uh = mp.mpf(2)**(-udeg)
    flint.ctx.prec = bits; flint.ctx.cap = K+1
    c = mp.cos(alpha); ca = flint.arb(mp.nstr(c, mdps+5))
    s = flint.acb_series([0,1], prec=K+1)
    sums = [flint.acb(0) for _ in range(K+1)]
    used=0; t0=time.time()
    for um, w0 in nodes:
        if um > umax: continue
        used += 1
        w = w0/(4*mp.cosh(um/2)**2)
        x = flint.acb(mp.nstr(um, mdps+5))
        ex = (-x).exp()                      # e^{-x}
        L = (1+ex).log()                     # L(x)
        eL = ((-L)*s).exp()                  # e^{-q L}
        eqx = ((-x)*s).exp()                 # e^{-q x}
        EZ = ex*(x*s).exp()                  # e^Z, Z = x(q-1)
        EmZ = (1/ex)*((-x)*s).exp()          # e^{-Z}
        SQ = (1 - 2*ca*EZ + EZ*EZ).sqrt()    # sqrt(1-2cos a e^Z + e^{2Z})
        Bt = EmZ*(SQ - 1)                    # B~_alpha(Z)
        Tn = eL*Bt/(x.exp()+1)
        Rn = (eL/(1+ex) - 1)*eqx
        G = Tn + Rn
        ww = flint.acb(mp.nstr(w*uh, mdps+5))
        cf = G.coeffs()
        if len(cf) < K+1: cf += [flint.acb(0)]*(K+1-len(cf))
        for n in range(K+1): sums[n] += ww*cf[n]
    return sums, used, time.time()-t0

if __name__ == "__main__":
    angles = [('pi_3', mp.pi/3), ('pi_2', mp.pi/2), ('2pi_3', 2*mp.pi/3),
              ('3pi_4', 3*mp.pi/4), ('5pi_6', 5*mp.pi/6)]
    os.makedirs('out', exist_ok=True)
    for name, alpha in angles:
        sums, used, dt = taylor_coeffs_TR(alpha)
        with open('out/mk_tr_%s.tsv'%name,'w') as f:
            for k,z in enumerate(sums):
                f.write("%d\t%s\t%s\n"%(k, z.mid().real.str(50), z.rad().str(4)))
        print("%-7s udeg=6 nodes=%d %.1fs  M0=%s  M300=%s rad=%s" % (
            name, used, dt, sums[0].mid().real.str(30), sums[300].mid().real.str(20), sums[300].rad().str(3)))
