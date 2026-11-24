# s07: exact degree floor(delta/2) -- top coefficients nonvanishing: off-diagonal (all d in domain) and STABLE DIAGONAL ehat(delta,top)
import sys
sys.dont_write_bytecode = True
import sympy as sp
import tow2lib as T2, tow3lib as T3
d, m = T3.d, T3.m; I = sp.I; p = sp.Symbol('p')
print('--- off-diagonal tops (normalised by e00) ---')
for delta in range(2, 7):
    s = delta//2
    f = sp.factor(T3.e(delta, s))
    print(f'  e({delta},{s})/e00 =', f)
rB0d = 2*I*(4*p-1)/(p+1)
print('--- stable diagonal tops ehat_pp(delta,top)/rho_pp(0) = (p+1) rB0d sum_{d\'} tau1_{delta-d\'} F(d\',top)(2p+1,p+1)/(4p-d\'-1) ---')
for delta in range(2, 7):
    s = delta//2
    tot = 0
    for dpp in range(0, delta+1):
        tot += T2.tau1f(delta-dpp)*T3.e(dpp, s, 2*p+1, p+1)/(4*p-dpp-1)
    val = sp.factor(sp.cancel((p+1)*rB0d*tot))
    num = sp.numer(sp.together(val))
    roots = sp.Poly(sp.numer(sp.together(val)), p).nroots(n=15) if sp.degree(num, p) > 0 else []
    rr = [r for r in sp.roots(sp.Poly(num, p)).keys()] if sp.degree(num, p) > 0 else []
    print(f'  delta={delta}, s={s}: ehat/rho_pp(0) =', val)
    print('       numerator roots in p:', [sp.nsimplify(r) for r in rr], ' (need no integer p >= 3 among them)')
# integer-root check p=3..200
for delta in (5, 6):
    s = delta//2
    tot = 0
    for dpp in range(0, delta+1):
        tot += T2.tau1f(delta-dpp)*T3.e(dpp, s, 2*p+1, p+1)/(4*p-dpp-1)
    val = sp.cancel((p+1)*rB0d*tot)
    zero = [pp for pp in range(3, 400) if sp.simplify(val.subs(p, pp)) == 0]
    print(f'delta={delta}: p in 3..399 with ehat(top)=0 : {zero}')
