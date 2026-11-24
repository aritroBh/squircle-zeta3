from mpmath import mp, mpf, nstr, pi, gamma, hyper, diffs, factorial, polylog
mp.dps = 60
def e2(ca):
    def Eafun(s):
        return gamma(1-s)*gamma(2+s)*hyper([1+s, -s],[2], (1+ca)/2)/(1-ca)
    return list(diffs(Eafun, 0, 2))[-1]/factorial(2)
for c in ['1/2','3/5','1/3','4/5','2/3']:
    cm = mpf(c); ca = 2*cm**2-1
    T = pi**2/12 - polylog(2, cm**2)/2
    lhs = e2(ca); rhs = T/((1-ca)/2)
    print(f"c={c}: e_2={nstr(lhs,25)}  2T/(1-ca)={nstr(rhs,25)}  |d|={nstr(abs(lhs-rhs),3)}")
