import sympy as sp
ln2 = sp.log(2); pi = sp.pi
def bp(n):
    if n==1: return (ln2-1)/2
    return sp.Rational((-1)**((n-1)//2), (n**2-1)*(2*((n-1)//4)+1))
# b_n'' exact forms (PSLQ-identified, verified 1e-56)
bpp = {1: sp.Rational(3,8) - ln2 + ln2**2/4 + 5*pi**2/96 - pi**4/768,
       3: (12 - 16*ln2 + 3*pi**2)/128,
       5: (508 + 48*ln2 - 45*pi**2)/3456,
       7: (-572 - 96*ln2 + 63*pi**2)/13824,
       9: (2532 + 160*ln2 - 225*pi**2)/64000,
       11: (-25592 - 2160*ln2 + 2475*pi**2)/1296000,
       13: (1598344 + 75600*ln2 - 143325*pi**2)/88905600}
S13 = sp.simplify(sum(bp(n)*bpp[n] for n in [1,3,5,7,9,11,13]))
print("Sigma_13 =", S13)
print(sp.N(S13, 30))
# separate by constant type
for c in [1, ln2, ln2**2, pi**2, pi**2*ln2, sp.zeta(3), pi**4]:
    print(c, ":", sp.simplify(S13.coeff(c) if c!=1 else S13.subs({ln2:0,pi:0}).subs({sp.zeta(3):0})))
