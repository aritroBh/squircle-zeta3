
import sys; sys.path.insert(0,'.')
import json
import mpmath as mp
import w5lib as W
ref = json.load(open("../M4i_w5_noL/Lchi8_5.json"))["value"]
print("ref digits:", len(ref))
P = 1000
v = W.Lchi_hurwitz(8, P)
mp.mp.dps = P
s = mp.nstr(v, P+5)
# count matching digits
m = 0
for a, b in zip(s, ref):
    if a == b: m += 1
    else: break
print("matching digits vs M4i ref:", m)
# brute-force Hurwitz with Euler-Maclaurin by hand for a=1/8
def hurwitz_em(s, a, P, N=None):
    mp.mp.dps = P + 40
    if N is None: N = int(P*1.2) + 50
    tot = mp.nsum(lambda n: 1/(n+a)**s, [0, mp.inf], method='euler-maclaurin', tol=mp.mpf(10)**(-P-20))
    return tot
def hurwitz_direct(s, a, P):
    mp.mp.dps = P + 40
    # sum N terms + integral tail + bernoulli correction terms
    N = 4*P + 100
    tot = mp.mpf(0)
    for n in range(N):
        tot += 1/(n+a)**s
    w = N + a
    tot += w**(1-s)/(s-1) + mp.mpf('0.5')*w**(-s)
    # Bernoulli corrections B_{2k}/(2k)! * s(s+1)...(s+2k-2) * w^(-s-2k+1)
    from mpmath import bernoulli, factorial, binomial
    for k in range(1, 40):
        coef = bernoulli(2*k)/factorial(2*k)
        rising = mp.mpf(1)
        for j in range(2*k-1): rising *= (s+j)
        tot += coef*rising*w**(-(s+2*k-1))
    return tot
mp.mp.dps = 500
h1 = mp.zeta(5, mp.mpf(1)/8)
h2 = hurwitz_direct(5, mp.mpf(1)/8, 500)
print("mp.zeta(5,1/8) - EM direct:", mp.nstr(h1-h2, 3))
# and full L from direct Hurwitz
mp.mp.dps = 500
Ld = sum(W.chi(8,a)*hurwitz_direct(5, mp.mpf(a)/8, 500) for a in range(1,8))/8**5
print("L(hurwitz direct) - mp.zeta-based:", mp.nstr(Ld - v, 3))
s2 = mp.nstr(Ld, 620)
m2 = 0
for a, b in zip(s2, ref):
    if a == b: m2 += 1
    else: break
print("direct-EM L matching digits vs M4i ref:", m2)
