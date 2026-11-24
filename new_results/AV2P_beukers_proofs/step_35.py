from mpmath import mp, mpf, log, pi, nprint, quad, clsin, binomial, sqrt, catalan
mp.dps = 60
# S1 target
K=400; S1 = mpf(0); Hk=mpf(0)
for k in range(0,K):
    if k>0: Hk += 1/mpf(k)
    S1 += binomial(2*k,k)*Hk/(mpf(16)**k*(2*k+1))
print("S1 series =", nprint(S1, 40))

# F(y) = sum C H_k y^k
def F(y):
    y = mpf(y); c = sqrt(1-4*y)
    return (2/c)*log((1+c)/(2*c))
y0 = mpf('0.02')
Fk = mpf(0); Hk=mpf(0)
for k in range(0,300):
    if k>0: Hk += 1/mpf(k)
    Fk += binomial(2*k,k)*Hk*y0**k
print("F closed vs series:", nprint(F(y0)-Fk, 5))

# S1 = int_0^1 F(s^2/16) ds
S1_int = quad(lambda s: F(s**2/16), [0, 1])
print("S1 via int F:", nprint(S1_int, 40), " diff:", nprint(S1_int-S1, 5))

# S1 = 4 int_0^{pi/6} log((1+cos t)/(2 cos t)) dt ?
cand = 4*quad(lambda t: log((1+mp.cos(t))/(2*mp.cos(t))), [0, pi/6])
print("cand =", nprint(cand, 40), " diff:", nprint(cand-S1, 5))

# maybe factor should be different: check int_0^1 F(s^2/16) ds via trig sub directly
# s = 2 sin th: F(s^2/16) = (2/c) log((1+c)/2c), c = cos th
S1_trig = quad(lambda th: (2/mp.cos(th))*log((1+mp.cos(th))/(2*mp.cos(th)))*2*mp.cos(th), [0, pi/6])
print("S1_trig =", nprint(S1_trig, 40), " diff:", nprint(S1_trig-S1, 5))
