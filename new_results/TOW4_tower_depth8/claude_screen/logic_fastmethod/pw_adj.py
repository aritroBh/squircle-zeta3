"""pw_adj.py -- independent pointwise evaluation of the adjacent Newton identity residual
   R(p) = [e(delta,s)+e(delta,s+1)] - sum(terms)    (formulas of prove_adj78_fast.feed_terms, re-typed),
   with F from the UNTRUNCATED pointwise recursion (pw.Fu), exact Gaussian rationals, no sympy cancel."""
import sys, random, time
from fractions import Fraction as Fr
import pw
from pw import G, I_, Fu, TAU1, TAU2

def ff(x, c):
    r = G(1)
    for u in range(c): r = r * (x - u)
    return r
def fact(n):
    r = 1
    for k in range(2, n+1): r *= k
    return r

def residual(delta, s, p, foldsign=1, rA0shift=0):
    p = G.of(p)
    d = 2*p - 1
    rA0 = I_*(4*p-3)/(p+1+rA0shift)
    rB0 = (4*p-1)*(4*p-3)/((p+1)*p)
    def eR(dl, ss, D, M, conj=False):
        if ss > dl//2: return G(0)
        v = Fu(dl, ss, D, M)
        return v.conj() if conj else v
    def tau1(n, conj=False): return (TAU2[n] * ((-1)**n)) if conj else TAU1[n]
    def tau2(n, conj=False): return (TAU1[n] * ((-1)**n)) if conj else TAU2[n]
    eh = {}
    def ehat(dl, ss, conj=False):
        k = (dl, ss, conj)
        if k in eh: return eh[k]
        tot = G(0)
        for dpp in range(0, dl+1):
            tot = tot + tau1(dl-dpp, conj) * eR(dpp, ss, 2*p+1, p+1, conj) * rB0 / (4*p-dpp-1)
        eh[k] = (p+1)*tot
        return eh[k]
    ed = {}
    def Ediag(dl, ss, conj=False):
        k = (dl, ss, conj)
        if k in ed: return ed[k]
        tot = G(0); t = 0
        while ss+t <= dl//2:
            tot = tot + ((-1)**t) * ehat(dl, ss+t, conj); t += 1
        ed[k] = tot
        return tot
    lhs = eR(delta, s, d, p) + eR(delta, s+1, d, p)
    terms = []
    for dp in range(0, delta+1):
        terms.append((p+1)*tau1(delta-dp)*eR(dp, s, 2*p, p+1)*rA0/(2*d-dp-1))
    for dp in range(0, delta+1):
        terms.append(p*tau2(delta-dp)*Ediag(dp, s)/(2*d-dp-1))
    for t in range(0, delta+1):
        mu = 4*p - 4 - t
        g = G(0)
        for dp in range(0, t+1):
            g = g - Ediag(dp, s, True)*tau1(t-dp, True)/(4*p-3-dp)
        binom = ff(mu-2, delta-t) / fact(delta-t)
        terms.append(foldsign * (-p) * g * ((-I_)**(-4-t)) * binom * (I_**(-6-delta)))
    tot = G(0)
    for x in terms: tot = tot + x
    return lhs - tot, lhs, terms

if __name__ == '__main__':
    cases = [tuple(map(int, a.split(','))) for a in sys.argv[2:]] or [(7,0),(7,1),(7,2),(7,3),(8,0),(8,1),(8,2),(8,3),(8,4)]
    NINT = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    rng = random.Random(5)
    for (delta, s) in cases:
        t0 = time.time()
        pts = [Fr(4 + j) for j in range(NINT)] + [Fr(rng.randint(40, 4000), rng.randint(1, 13)) + Fr(1, 7) for _ in range(6)]
        pts += [Fr(3), Fr(2), Fr(1)] if True else []
        nz = 0; big = 0
        poles = []
        for pp in pts:
            try:
                r, lhs, terms = residual(delta, s, pp)
            except ZeroDivisionError:
                poles.append(pp); continue
            if not r.iszero(): nz += 1
            if not lhs.iszero(): big += 1
        print(f'pointwise adjacent ({delta},{s}): {len(pts)} points (p=4..{3+NINT}, 6 random rationals, p=3,2,1): nonzero residuals = {nz}; poles(skipped) at p={[str(x) for x in poles]};  lhs nonzero at {big}/{len(pts)}  ({round(time.time()-t0,1)}s)', flush=True)
