
import mpmath as mp
mp.mp.dps = 50
# direct quadrature of J((0,0,2,0,3); 0->1): innermost J((3);u)=log(1-u/3); then J((0,3);u2) = -Li2(u2/3);
# J((2,0,3);u1) = ∫₀^u1 du2/(u2-2)*(-Li2(u2/3)); J((0,2,0,3);u0)=∫ du1/u1 J((2,0,3);u1); outer again.
def J_rest(u1):
    if u1 == 0: return mp.mpf(0)
    return mp.quad(lambda u2: -mp.polylog(2, u2/3)/(u2 - 2), [0, u1])
Jval = mp.quad(lambda u0: (mp.quad(lambda u1: J_rest(u1)/u1, [0, u0]))/u0 if u0 > 0 else 0, [0, 1])
# simpler: J((0,0,2,0,3);1) = ∫₀¹ du0/u0 ∫₀^u0 du1/u1 ∫₀^u1 du2/(u2-2) * (-Li2(u2/3))
print("J((0,0,2,0,3);1) =", mp.nstr(Jval, 20))

def li32(x, y, N=1500):
    x = mp.mpf(x); y = mp.mpf(y)
    S = mp.mpf(0); ym = y; xn = x; tot = mp.mpf(0)
    for n in range(2, N + 1):
        S += ym / ((n-1)*(n-1)); ym *= y; xn *= x
        tot += xn * S / (n*n*n)
    return tot
print("-Li32(1/2,2/3) =", mp.nstr(-li32(mp.mpf(1)/2, mp.mpf(2)/3), 20))
print("+Li32(1/2,2/3) =", mp.nstr(li32(mp.mpf(1)/2, mp.mpf(2)/3), 20))
