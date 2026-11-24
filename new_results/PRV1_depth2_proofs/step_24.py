
import mpmath as mp
mp.mp.dps = 80
def depth2(x, y, p=2):
    rho = max(abs(x), abs(x*y)); assert 0 < rho < 1
    target = mp.power(10, -mp.mp.dps - 5)
    nmax = int(mp.ceil(mp.log(target*(1-rho))/mp.log(rho)))
    h, xy_power, total = mp.mpf(0), x*y, mp.mpf(0)
    for n in range(1, nmax+1):
        total += h/n**p
        h = x*h + x*xy_power/n
        xy_power *= x*y
    return total

def H_dec(t): return depth2(1-t, -t/(1-t))
def H_full(t):
    return (mp.polylog(3, 1-t**2) - mp.polylog(3, 1/(1+t)) + mp.polylog(3, t/(1+t))
            - mp.log(1+t)*(mp.polylog(2, -t)/2 + 3*mp.polylog(2, 1-t)/2
                           + mp.polylog(2, 1-t**2)/4 + mp.zeta(2)/4))
for t0 in ['0.13','0.4','0.6180339887498948482045868343656381177','0.9']:
    t = mp.mpf(t0)
    print('t=', t0, ' residual=', mp.nstr(H_dec(t)-H_full(t), 8))

# exact specialisation with sympy
import sympy as sp
z2, z3, l = sp.symbols('z2 z3 l')   # zeta(2), zeta(3), log(phi)
L2a, L2b, L2c = sp.symbols('L2a L2b L2c')   # Li2(1/phi^2), Li2(1/phi), Li2(-1/phi)
L3a, L3b, L3c = sp.symbols('L3a L3b L3c')   # Li3(1/phi^2), Li3(1/phi), Li3(-1/phi)
# golden dilog equations: Euler@1/phi^2, Landen@1/phi^2, duplication@1/phi
eqs = [L2a + L2b - (z2 - 2*l**2),
       L2c + L2a + l**2/2,
       L2a - 2*L2b - 2*L2c]
sold = sp.solve(eqs, [L2a, L2b, L2c], dict=True)[0]
print("golden dilogs:", sold)
# golden trilog: 3-term Landen @1/phi^2 and duplication @1/phi
eqs3 = [L3a + L3b + L3c - (z3 - l*z2 + 5*l**3/6),
        L3a - 4*L3b - 4*L3c]
solt = sp.solve(eqs3, [L3a, L3b, L3c], dict=True)[0]
print("golden trilogs (Lewin, sum):", solt)

# Lam_dec = zeta(3) + l*Li2(1/phi^2)/2 - 2 Li3(1/phi)
Lam = sp.expand((z3 + l*L2a/2 - 2*L3b).subs(sold))
print("Lam_dec reduced =", Lam, "  [expr in Li3(1/phi), Li3(-1/phi) via sum]")
# DEC1 form:
Lam_DEC1 = 2*z3/5 + 4*l*z2/5 - l**3 + 3*L3c + L3b
diffLam = sp.expand((Lam - Lam_DEC1).subs({L3b: solt[L3b].subs({}) , L3c: L3c}))
# substitute Li3(1/phi) = S - Li3(-1/phi), S = solt[L3b]+solt[L3c] relation: L3b = S - L3c with S = (z3-l*z2)/5 + l^3/6
S = (z3 - l*z2)/5 + l**3/6
diffLam = sp.expand(Lam - Lam_DEC1).subs(L3b, S - L3c)
print("Lam: mine - DEC1 (mod sum relation) =", sp.simplify(diffLam))

# swap_dec = Li3(1/phi^2) - l*( Li2(-1/phi)/2 + 3 Li2(1/phi^2)/2 + Li2(1/phi)/4 + z2/4 )
swap = sp.expand((L3a - l*(L2c/2 + 3*L2a/2 + L2b/4 + z2/4)).subs(sold))
print("swap_dec reduced (before Lewin) =", swap)
swap_full = sp.expand(swap.subs(L3a, solt[L3a]))
print("swap_dec fully reduced =", swap_full)
swap_DEC1 = 3*z3/5 - 7*l*z2/5 + 2*l**3 + L3c + L3b
diffSwap = sp.expand(swap_full - swap_DEC1).subs(L3b, S - L3c)
print("swap: mine - DEC1 (mod sum relation) =", sp.simplify(diffSwap))
