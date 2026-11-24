from mpmath import mp, mpf, zeta, log, pi, nprint, exp, loggamma, digamma, polygamma
mp.dps = 160
def phi(w):
    w = mpf(w)
    return w**3 * exp(2*loggamma(w) - loggamma(2*w+1)) * (digamma(2*w+1) - digamma(w))
def dphi(w):
    w = mpf(w)
    pw = phi(w)
    D = digamma(2*w+1)-digamma(w)
    return pw*(3/w + 2*digamma(w) - 2*digamma(2*w+1) + (2*polygamma(1,2*w+1)-polygamma(1,w))/D)
# d1 = -sum_{n>=1} phi'(n)/n^2
S = mpf(0); n = 1
while True:
    t = dphi(n)/mpf(n)**2
    S += t
    if n > 200 and abs(t) < mpf(10)**(-155):
        break
    n += 1
d1 = -S
print("n terms:", n)
print("d1 =", nprint(d1, 150))
with open('d1.txt','w') as f:
    f.write(str(d1))
with open('progress.log','a') as fh:
    fh.write("step_30: d1 = -sum phi'(n)/n^2 at 150+ digits\n")
