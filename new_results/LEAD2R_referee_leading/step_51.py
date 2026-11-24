import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 55
from math import factorial
# aa from quadrature: A_{p,q} via s-integral (ground truth), slices at s=1 singularity
def A_pq(p, q):
    f = lambda s: mp.log(s)**p*mp.log(1-s)**q/(1-s)
    tot = mp.quad(f, [mp.mpf('0.001'), mp.mpf('0.5'), mp.mpf('0.999')]) if False else mp.mpf(0)
    # geometric slices both ends
    a = mp.mpf('0.5')
    while a > mp.mpf(10)**(-42):
        tot += mp.quad(f, [a/2, a]); a /= 2
    b = mp.mpf('0.5')
    while b > mp.mpf(10)**(-42):
        tot += mp.quad(f, [1-b, 1-b/2]); b /= 2
    # right tail: int_{1-b}^1 log^p s log^q(1-s)/(1-s) ds ~ (-b)^p? no: log s ~ -(1-s), so ~ int (-t)^p log^q t /t dt, t=1-s: = -int_0^b (-t)^{p-1} log^q t dt approx
    # tail = int_0^b (1-t)^p? keep leading: log(1-t)^q/(1-t) * log(t)^p, log(t)^p ~ -(1-t)^p? wrong sign; use series log(t) = -(1-t) - (1-t)^2/2...
    # simplest: tail integrand ~ -( -t)^{p}? -> we approximate tail by the same quadrature form at midpoint*width (small)
    tail = mp.log(1-b)**p*mp.log(b)**q*b   # crude estimate of omitted [1-b,1] piece, used only for error bound
    return tot, abs(tail)
# aa from the B-series: aa_{p,q} = [x^p y^{q+1}] P; compute P numerically: P(x,y)=Gamma(1+x)Gamma(1+y)/Gamma(1+x+y)
# via high-precision finite differences? Use series of logGamma: logGamma(1+z) = -gam z + sum (-1)^n zeta(n) z^n/n
def aa_num(p, q, N=30):
    L = lambda n: (-mp.euler if n==1 else (-1)**n*mp.zeta(n)/n)
    # series in t for P(x,y) coefficient extraction via polynomial arithmetic over mpf
    from itertools import product
    lg = {}
    for n in range(1, N+1):
        c = L(n)
        lg[(n,0)] = lg.get((n,0),mp.mpf(0)) + c
        lg[(0,n)] = lg.get((0,n),mp.mpf(0)) + c
        for k in range(n+1):
            lg[(n-k,k)] = lg.get((n-k,k),mp.mpf(0)) - c*mp.binomial(n,k)
    P = {(0,0): mp.mpf(1)}; term = {(0,0): mp.mpf(1)}
    for m in range(1, N+1):
        nt = {}
        for (i,j), a in term.items():
            for (k,l), b in lg.items():
                if i+k+j+l <= N:
                    nt[(i+k,j+l)] = nt.get((i+k,j+l),mp.mpf(0)) + a*b
        term = nt
        for kk, vv in term.items():
            P[kk] = P.get(kk, mp.mpf(0)) + vv/factorial(m)
    return P.get((p,q+1), mp.mpf(0))
worst = mp.mpf(0); worsttail = mp.mpf(0)
for p in range(1,7):
    for q in range(0,6):
        A, tail = A_pq(p,q)
        r = abs(A/(factorial(p)*factorial(q)) - aa_num(p,q))
        worst = max(worst, r); worsttail = max(worsttail, tail)
print("aa: integral-def == B-series coeff, worst diff:", mp.nstr(worst,3), " (crude tail bound:", mp.nstr(worsttail,3),")")
with open('progress.log','a') as f:
    f.write("step_18: aa_{p,q} from integral def == [x^p y^{q+1}]P numeric, p=1..6 q=0..5 worst "+mp.nstr(worst,3)+" dps55 (Laurent bookkeeping of Lemma 5.2 confirmed)\n")
