from squircle import setdps
from mpmath import mp, nstr, log, power, fsum, mpf
setdps(60)
s = 4
P = mpf(0); acc = mpf(0)
for k in range(1, 12):
    if k >= 2:
        term = power(k, -s-1)*log(k)*acc
        P += term
        print(f"k={k}: acc={nstr(acc,10)} term={nstr(term,10)} P={nstr(P,12)}")
    acc += k*log(k)
