
import mpmath as mp, flint, time, os

def taylor_coeffs_T(alpha, K=300, udeg=7, mdps=100, bits=768, umax=1200):
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
        w = w0*uh
        x = flint.acb(mp.nstr(um, mdps+5))
        ex = (-x).exp()
        L = (1+ex).log()
        eL = ((-L)*s).exp()
        Z = ex*(x*s).exp()                     # e^{+z}
        SQ = (1 - 2*ca*Z + Z*Z).sqrt()
        Bt = (Z - 2*ca)/(1 + SQ)               # rationalized kernel B~_alpha(z)
        G = eL*Bt/(x.exp()+1)
        ww = flint.acb(mp.nstr(w, mdps+5))
        cf = G.coeffs()
        if len(cf) < K+1: cf += [flint.acb(0)]*(K+1-len(cf))
        for n in range(K+1): sums[n] += ww*cf[n]
    return sums, used, time.time()-t0

def R_coeffs(Kmax=25, dps=80):
    mp.mp.dps = dps
    L = lambda x: mp.log1p(mp.exp(-x))
    out = []
    for k in range(Kmax+1):
        if k==0:
            Ik = mp.quad(lambda x: 1/(1+mp.exp(-x)) - 1, [0, mp.inf])
        else:
            Ik = mp.quad(lambda x: (x+L(x))**k/(1+mp.exp(-x)) - x**k, [0, 5, 20, mp.inf])
        out.append(((-1)**k)*Ik/mp.factorial(k))
    return out

angles = [('pi_3', mp.pi/3), ('pi_2', mp.pi/2), ('2pi_3', 2*mp.pi/3),
          ('3pi_4', 3*mp.pi/4), ('5pi_6', 5*mp.pi/6)]
os.makedirs('out', exist_ok=True)
Rk = R_coeffs(25)
mp.mp.dps = 60
print("R_k computed; R_25 =", mp.nstr(Rk[25], 8))
for name, alpha in angles:
    sums, used, dt = taylor_coeffs_T(alpha)
    with open('out/mk_trT_%s.tsv'%name,'w') as f:
        for k,z in enumerate(sums):
            f.write("%d\t%s\t%s\n"%(k, z.mid().real.str(50), z.rad().str(4)))
    print("%-7s nodes=%d %.1fs  T0=%s  T300=%s rad=%s" % (
        name, used, dt, sums[0].mid().real.str(25), sums[300].mid().real.str(18), sums[300].rad().str(3)))
