# SM6 item 3: PSLQ at weight 6.  (a) c6 numeric (Richardson, 52 digits) vs pure basis; (b) free moment gamma32 (dps 250) vs alternating basis.
import mpmath as mp
z = mp.zeta
mp.mp.dps = 260
L2 = mp.log(2)
def ig(a, b):   # gamma_ab = 2^(-a-1)/(a! b!) int_0^1 (-log x)^a log^b(1+x)/x dx
    f = lambda x: (-mp.log(x))**a*mp.log(1+x)**b/x
    return mp.quad(f, [0, mp.mpf(1)/4, mp.mpf(1)/2, 1])*mp.mpf(2)**(-a-1)/(mp.factorial(a)*mp.factorial(b))
g32 = ig(3, 2)
# zeta(-5,-1) := sum_{n>m>=1} (-1)^(n+m)/(n^5 m) = (1/24) int_0^1 (-log x)^4 log(1-x)/(1+x) dx
Z51 = mp.quad(lambda x: (-mp.log(x))**4*mp.log(1-x)/(1+x), [0, mp.mpf(1)/2, 1])/24
# sanity vs direct sum at low precision
mp.mp.dps = 30
s = mp.mpf(0); A = mp.mpf(0)
for n in range(1, 4001):
    s += (-1)**n/mp.mpf(n)**5*A
    A += (-1)**n/mp.mpf(n)
print("direct partial sum of zeta(-5,-1):", mp.nstr(s, 12), " integral:", mp.nstr(Z51, 12))
mp.mp.dps = 260
Li = lambda n: mp.polylog(n, mp.mpf(1)/2)
basis6 = [("zeta(6)", z(6)), ("zeta(3)^2", z(3)**2), ("zeta(5)log2", z(5)*L2), ("zeta(4)log2^2", z(4)*L2**2),
          ("zeta(3)log2^3", z(3)*L2**3), ("zeta(2)log2^4", z(2)*L2**4), ("log2^6", L2**6), ("Li6(1/2)", Li(6)),
          ("Li5(1/2)log2", Li(5)*L2), ("Li4(1/2)log2^2", Li(4)*L2**2), ("zeta(2)Li4(1/2)", z(2)*Li(4)),
          ("zeta(2)zeta(3)log2", z(2)*z(3)*L2), ("zeta(-5,-1)", Z51)]
def pslq(target, names_vals, dps, maxc=10**9, label=""):
    mp.mp.dps = dps
    vec = [target] + [v for _, v in names_vals]
    r = mp.pslq(vec, maxcoeff=maxc, maxsteps=10**6)
    print(label, "->", None if r is None else r)
    if r:
        for c, (n, _) in zip(r[1:], names_vals):
            if c: print("    %s: %d" % (n, c))
        print("    target coefficient:", r[0])
    return r
print("--- basis independence (no target), 13 elements at 200 digits")
mp.mp.dps = 200
print(mp.pslq([v for _, v in basis6], maxcoeff=10**9, maxsteps=10**6))
# (a) c6 numeric
c6 = mp.mpf("0.0429508278802840339294338048351352655102967276920578653499916")
print("--- (a) c6 numeric, 50 digits")
pslq(c6, basis6[:2], 50, 10**6, "c6 vs {zeta6, zeta3^2}")
pslq(c6, basis6[:3], 50, 10**6, "c6 vs +zeta5 log2")
pslq(c6, basis6[:11]+[basis6[11]], 50, 10**4, "c6 vs 12-element alternating set (no zeta(-5,-1))")
pslq(c6, basis6, 50, 10**3, "c6 vs all 13 (50 digits, maxcoeff 1e3)")
# (b) gamma32
mp.mp.dps = 260
print("--- (b) gamma32 at 250 digits")
print("gamma32 =", mp.nstr(g32, 60))
pslq(g32, basis6[:2], 250, 10**12, "gamma32 vs {zeta6, zeta3^2}")
r = pslq(g32, basis6, 250, 10**15, "gamma32 vs all 13")
# c6 exact from the 13-term decomposition
if r and r[0] != 0:
    mp.mp.dps = 260
    c6x = 40*g32 + mp.mpf(11)/16*z(6) - z(3)**2/2
    print("c6 (reduced, quadrature) =", mp.nstr(c6x, 55))
    pslq(c6x, basis6, 250, 10**15, "c6 vs all 13")
    pslq(c6x, basis6[:2], 250, 10**15, "c6 vs {zeta6, zeta3^2}")
