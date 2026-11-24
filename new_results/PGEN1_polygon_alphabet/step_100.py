# step_72: which L3 convention satisfies Landen at x=-2 and x=2? direct real-part vs inversion
import mpmath as mp
mp.mp.dps = 80
def Li3(x): return mp.re(mp.polylog(3, x))
def Li2(x): return mp.re(mp.polylog(2, x))
def L3_direct(x):
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(abs(1-x))
def L3_inv(x):
    x = mp.mpf(x)
    if x > 1: return L3_direct(1/x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(1-x)
z3 = mp.zeta(3)
for name, L in [("direct", L3_direct), ("inversion", L3_inv)]:
    for xx in ['-2', '2', '-0.4', '2.5', '0.3']:
        x = mp.mpf(xx)
        s = L(x)+L(1-x)+L(1-1/x)
        print(f"{name} Landen at {xx}: sum - zeta3 = {mp.nstr(s-z3, 6)}")
    print()
# also: Re Li2(2), Re Li3(2) values
print("Re Li2(2) =", Li2(2), " Re Li3(2) =", Li3(2))
