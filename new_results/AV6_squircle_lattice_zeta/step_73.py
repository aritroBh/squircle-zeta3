code = r'''
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, mpf, zeta, gamma, exp, fsum, pi
setdps(60)
s = 4
def inner(k, p):
    kf = mpf(k)
    return fsum([log(1+power(1-mpf(d)/kf, p)) for d in range(1, k)])
# brute U_1(128) to K=3000 + continuum tail
p = 128
P = mpf(0)
for k in range(2, 3001):
    P += power(k, -s)*inner(k, p)
# tail k>3000: inner(k) ~ k*I0q - ln2/2 : I0q = int_0^1 ln(1+x^p) dx
I0q = fsum([(-1)**(n+1)/(n*(n*p+1)) for n in range(1, 400)])
tail = I0q*(zeta(s-1) - fsum([power(k,1-s) for k in range(1,3001)])) - log(2)/2*(zeta(s)-fsum([power(k,-s) for k in range(1,3001)]))
U1_brute = P + tail
print("U1(128) brute:", nstr(U1_brute, 25))
print("U1(128) Uqdir:", nstr(mpf(__import__('json').load(open("Uq_data_s4.json"))["128"][0]), 25))
print("D brute:", nstr(U1_brute*p**(s-1), 20))
'''
exec(code)
