"""step_19.py -- high-precision T(r) identification, cross-checked at 3 angles."""
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 120

def T_alpha_c(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def Rw2(w): return 1/(2*w) + 1/(2*(w-1)) - 1/(w+1) - 1/(2*(w+r))
    return mp.quad(lambda w: Rw2(w)*lt(w), [r, mp.sqrt(r), 1])

def basis(c):
    r = (1-c)/(1+c)
    L = mp.log; A = lambda x: mp.polylog(2, x)
    # curated: expect functional independence
    return [('pi^2', mp.pi**2),
            ('Li2(r)', A(r)), ('Li2(-r)', A(-r)), ('Li2(c)', A(c)), ('Li2(-c)', A(-c)),
            ('Li2((1-c)/2)', A((1-c)/2)), ('Li2((1+c)/2)', A((1+c)/2)),
            ('log2^2', L(2)**2), ('logc^2', L(c)**2), ('log(1+c)^2', L(1+c)**2),
            ('log(1-c)^2', L(1-c)**2), ('log2logc', L(2)*L(c)), ('log2log(1+c)', L(2)*L(1+c)),
            ('logclog(1+c)', L(c)*L(1+c)), ('logclog(1-c)', L(c)*L(1-c)), ('log(1+c)log(1-c)', L(1+c)*L(1-c))]

angles = [mp.mpf(1)/2, mp.mpf(3)/5, mp.mpf(1)/3]
Ts = [T_alpha_c(c) for c in angles]
Bs = [[x for _, x in basis(c)] for c in angles]
names = [n for n, _ in basis(angles[0])]

def identify(T, Bvals, Bnames):
    Bv = list(Bvals); Bn = list(Bnames)
    while True:
        v = [T] + Bv
        rel = mp.pslq(v, tol=mp.mpf(10)**(-105), maxcoeff=10**6, maxsteps=4000)
        if rel is None: return None
        if rel[0] != 0:
            return rel, Bn
        idx = max(i for i in range(1, len(rel)) if rel[i] != 0)
        del Bv[idx-1]; del Bn[idx-1]

out = identify(Ts[0], Bs[0], names)
if out:
    rel, keptnames = out
    coef = [mp.mpf(x)/(-rel[0]) for x in rel]
    print("relation at c=1/2:")
    for i, cf in enumerate(coef):
        if cf != 0 and i >= 1:
            print(f"   {mp.nstr(cf,10)}  {keptnames[i-1]}")
    # cross-check at c=3/5, 1/3
    for c, T, B in zip(angles, Ts, Bs):
        bmap = {n: x for n, x in zip(names, B)}
        pred = sum(cf*bmap[keptnames[i-1]] for i, cf in enumerate(coef) if i >= 1 and cf != 0)
        print(f"   check c={mp.nstr(c,8)}: diff={mp.nstr(abs(T-pred),3)}")
print("STEP19 DONE")
