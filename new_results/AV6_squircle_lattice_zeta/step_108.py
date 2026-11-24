from squircle import setdps
from mpmath import mp, nstr, zeta, log, power, fsum, diff as mpdiff
setdps(60)
def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)
print("zeta''(3) =", nstr(zetapp(3), 25))
print("zeta''(4) =", nstr(zetapp(4), 25))
# brute check zeta''(3)
b = fsum([power(k,-3)*log(k)**2 for k in range(1, 50000)])
print("brute zeta''(3) to 50000:", nstr(b, 20))
# A1 = 4 zeta''(s-1) - 2 U_1 + U_2
A1id = 4*zetapp(3) - 2*0.63548396142868020832 + 0.25681321239117625475
print("A1 via identity:", nstr(A1id, 20))
print("A1 direct: 0.20513537204962184746")
