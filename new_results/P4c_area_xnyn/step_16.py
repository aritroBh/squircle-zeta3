# step_14.py : correct Tauberian lemma + direct Delta1 evaluation
from mpmath import mp, mpf, sqrt, log, nstr, pi, zeta, gamma, inf, quad, cosh, exp
import sympy as sp, math
mp.dps = 60
L2 = log(2); z2 = zeta(2)

zz = sp.symbols('zz')
Eg_series = sp.series(sp.gamma((1+zz)/2)*sp.gamma(1-zz/2)/sp.sqrt(sp.pi), zz, 0, 8).removeO().expand()
eks = [mpf(float(sp.N(Eg_series.coeff(zz, k), 30))) for k in range(8)]
mtil = [sum(eks[:k+1]) for k in range(8)]   # m̃_k = Σ_{i≤k} e_i ; M_g = Σ m̃_k z^{k-1}

def fW(w):
    w = mpf(w); return 1 + w - sqrt(1+w*w)
def gW(w): return fW(w)/w
def fp(w):
    w = mpf(w); return 1 - w/sqrt(1+w*w)

# CORRECT LEMMA: T_k(N) = qinf L^{k+1}/(k+1) + (q0-qinf) Σ_{j<=k} k!/j! L^j + Σ_{i<=k} k!/(k-i)! mu_i L^{k-i}
# with M_q(z) = q0/z + Σ mu_i z^i at 0.
def T_lemma(L, k, q0, qinf, mu):
    s = qinf*L**(k+1)/mpf(k+1)
    S = sum(mpf(math.factorial(k))/mpf(math.factorial(j))*L**j for j in range(k+1))
    s += (q0-qinf)*S
    for i in range(k+1):
        s += mpf(math.factorial(k))/mpf(math.factorial(k-i))*mu[i]*L**(k-i)
    return s
# g: q0=qinf=1, mu_i = m̃_{i+1}
mu_g = mtil[1:]
# f': q0=1, qinf=0, mu_i = e_{i+1}
mu_fp = eks[1:]
for N in [mpf(100)*sqrt(2), mpf(10000)*sqrt(2)]:
    L = log(N)
    for k in [1,2,3]:
        tnum = quad(lambda w: log(N/w)**k*gW(w), [0,1,N])
        print("g : N=%s k=%d num=%s pred=%s diff=%s" % (nstr(N,7),k,nstr(tnum,16),nstr(T_lemma(L,k,1,1,mu_g),16),nstr(tnum-T_lemma(L,k,1,1,mu_g),5)))
    for k in [1,2,3]:
        tnum = quad(lambda w: log(N/w)**k*fp(w), [0,1,N])
        print("f': N=%s k=%d num=%s pred=%s diff=%s" % (nstr(N,7),k,nstr(tnum,16),nstr(T_lemma(L,k,1,0,mu_fp),16),nstr(tnum-T_lemma(L,k,1,0,mu_fp),5)))

# Direct Delta1 = ∫_0^∞ e^{-eσ} ∫_0^σ [f(W0 K) - f(W0)] dδ dσ, W0 = n√2 e^{-σ}, K = e^{eσ/2}√cosh(eδ)
def D1_direct(n):
    n = mpf(n); e = 1/n
    def inner(sg):
        W0 = n*sqrt(2)*exp(-sg)
        return quad(lambda d: fW(W0*exp(e*sg/2)*sqrt(cosh(e*d))) - fW(W0), [0, sg/2, sg])
    L = log(n)
    return quad(lambda sg: exp(-e*sg)*inner(sg), [0, L/2, L, 2*L, 5*L])
for n in [100, 1000, 10000]:
    n_ = mpf(n); e = 1/n_
    d1 = D1_direct(n)
    N = n_*sqrt(2); L = log(N)
    t2p = T_lemma(L, 2, 1, 0, mu_fp)
    pred = (e/2)*t2p
    print("n=%d: D1=%s  (e/2)T2'=%s  ratio=%s  diff=%s  e^2 L^4=%s" % (n, nstr(d1,18), nstr(pred,18), nstr(d1/pred,10), nstr(d1-pred,10), nstr(e**2*L**4,6)))
