from flint import arb
import math

rho = arb.pi()/2
Estar = arb(5)/24

def sgn_xa(n):  # xa_n = -2cos(n pi/2 + pi/4): pattern -,+ ,+,- for n mod 4 = 0,1,2,3 times sqrt(2)
    return [-1, 1, 1, -1][n % 4]
def sgn_xb(n):  # xb_n = cos(n pi/2 - pi/4): pattern +,+,-,- times sqrt(2)/2
    return [1, 1, -1, -1][n % 4]
def sgn_law(n): # -sqrt2 cos((2n+1)pi/4): pattern -,+,+,- times 1 (|.|=1)
    return [-1, 1, 1, -1][n % 4]

sqrt2 = arb(2).sqrt()
def A_n(n):
    return sgn_xa(n)*sqrt2*(2/arb.pi())*rho**(arb('0.5')-n)*arb.beta(arb('1.5'), arb(n)-arb('0.5'))
def B_n(n):
    return sgn_xb(n)*(sqrt2/2)*(2/arb.pi())*rho**(arb('1.5')-n)*arb.beta(arb('2.5'), arb(n)-arb('1.5'))
def law2(n):
    return sgn_law(n)*1*rho**(-n)*arb(n)**arb('-1.5')*(1 + (3-3*arb.pi()*((-1)**n))/(8*n))
def Ebnd(n):
    return (4/arb.pi())*Estar*rho**(arb('2.5')-n)*arb.beta(arb('3.5'), arb(n)-arb('2.5'))
def tails(n):
    rhot = arb('2.9')*rho
    return (4/arb.pi())*rhot**(arb('0.5')-n)/(arb(n)-arb('0.5')) + (2/arb.pi())*rhot**(arb('1.5')-n)/(arb(n)-arb('1.5'))

mx = None; argmx=None
for n in range(10, 501):
    u = rho**(-n)*arb(n)**arb('-3.5')
    val = (abs(A_n(n)+B_n(n)-law2(n)) + Ebnd(n) + tails(n))/u
    if mx is None or val.upper() > mx: mx, argmx = val.upper(), n
print("certified sup_{10<=n<=500} R_bound/u =", mx, " at n =", argmx)

# Stirling tail for n > 500: bound each piece.
# logGamma(x) = (x-1/2)ln x - x + ln(2pi)/2 + theta/(12x), 0<theta<1 (classical).
# Simpler: use Gautschi-type: for y >= 10: Gamma(y-5/2)/Gamma(y+1) <= 1.01*y^(-7/2) etc.
# Do it honestly with the log expansion carried to 1/n^2 with explicit constants.
def stirling_bounds(n0):
    # returns dict of uniform bounds for n >= n0 on: gA = Gamma(n-1/2)/Gamma(n+1)*n^{3/2},
    # gB = Gamma(n-3/2)/Gamma(n+1)*n^{5/2}, gC = Gamma(n-5/2)/Gamma(n+1)*n^{7/2}
    # via logGamma bounds; compute at n0 with arb (decreasing to 1, so sup at n0) --
    # to PROVE decrease we'd need more; instead bound directly for all n>=n0 via expansion:
    # log(Gamma(n+a)/Gamma(n+b)) = (a-b)ln n + (c1(a,b))/n + R, |R| <= C/n^2 (from Stirling with theta in (0,1))
    pass

# Direct approach: for n >= n0, write L(n) = log Gamma(n-5/2) - log Gamma(n+1) + (7/2) log n
# using logGamma(x) = (x-1/2)lnx - x + ln(2pi)/2 + th/(12x):
# L(n) = (n-3)ln(n-5/2) - (n+1/2)ln(n+1) + 7/2 + (7/2) ln n + err, |err| <= 1/(12(n-2.5)) 
#      = (n-3) ln(1 - 2.5/n) - (n+1/2) ln(1 + 1/n) + 7/2 + err
# For n >= 500, bound the two logs by Taylor with remainder:
# (n-3)ln(1-2.5/n) = -2.5 + (2.5*3 - 2.5^2/2)/n + r1, |r1| <= C1/n^2
# (n+1/2)ln(1+1/n) = 1 + (1/2*(-1) + 1/2)/n ... compute numerically:
import mpmath as mp
mp.mp.dps = 50
for nn in [500, 1000, 10000]:
    n = mp.mpf(nn)
    L = (n-3)*mp.log(1-mp.mpf('2.5')/n) - (n+mp.mpf('0.5'))*mp.log(1+1/n) + mp.mpf('3.5')
    print(f"n={nn}: L(n) = {mp.nstr(L, 20)}  (limit should be 35/8/n = {mp.nstr(mp.mpf('4.375')/n, 8)})")
with open('progress.log','a') as f:
    f.write(f"step_09: certified sup R_bound/u over 10..500 = {mx} at n={argmx}\n")
