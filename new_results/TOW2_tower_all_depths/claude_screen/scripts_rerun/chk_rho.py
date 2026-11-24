import time, sympy as sp, tow2lib as T
d, m = T.d, T.m; i = sp.I
print('pickle keys:', sorted(T._ac.keys()))
# doc closed forms a3, a4 (proof_all_depths.md section 1)
a3 = (4*i*d**4 - 36*i*d**3 + (6+125*i)*d**2 - 12*d*m + (-18-198*i)*d + 36*m + 117*i)/(24*d - 36)
a4 = (16*d**6 - 240*d**5 + (1528-48*i)*d**4 + 96*i*d**3*m + (-5296+432*i)*d**3 - 864*i*d**2*m + (10525-1308*i)*d**2
      + (48+2616*i)*d*m + (-11321+1320*i)*d - 48*m**2 - 2640*i*m + 5040)/(384*d**2 - 1536*d + 1440)
for dl, doc in [(3, a3), (4, a4)]:
    t0 = time.time()
    cb = T.complete_bell(dl)
    print(f'delta={dl}: complete_bell == doc a_delta :', sp.cancel(sp.together(cb - doc)) == 0, f'({time.time()-t0:.1f}s)', flush=True)
    print(f'delta={dl}: arho_simple  == doc a_delta :', sp.cancel(sp.together(T.arho_simple(dl) - doc)) == 0, flush=True)
    print(f'delta={dl}: arho (cached bell) == arho_simple:', sp.cancel(sp.together(T.arho(dl) - T.arho_simple(dl))) == 0, flush=True)
# a1, a2
print('a1 = -i(d-2):', sp.cancel(T.arho_simple(1) + i*(d-2)) == 0 if 1 in T._ac else 'n/a', ' bell:', sp.cancel(T.complete_bell(1) + i*(d-2)) == 0)
# TOW1 rho(2) = 1/2[(3d-3)/2 - i/2 r_p - (d-2)^2], r_p = 2(d-2m) - 2 d i - 4(d-1)(d-2m)/(2d-3)
r = 2*(d-2*m) - 2*d*i - 4*(d-1)*(d-2*m)/(2*d-3)
a2 = sp.Rational(1,2)*((3*d-3)/sp.Integer(2) - i*r/2 - (d-2)**2)
print('a2 vs TOW1 rho(2):', sp.cancel(sp.together(T.complete_bell(2) - a2)) == 0)
# tau table vs independent series expansion of DL1, DL2 at v=i
v = sp.symbols('v'); w = sp.symbols('w')
a = 1+sp.sqrt(2); b = sp.sqrt(2)-1
DL1 = 1/(v-b) + 1/(v+a) - 1/v
DL2 = 1/(v-a) + 1/(v+b) - 1/v
ok = True
for name, DL, tab in [('tau1', DL1, T.tau1f), ('tau2', DL2, T.tau2f)]:
    ser = sp.series(DL.subs(v, i+w), w, 0, 8).removeO()
    for n in range(0, 7):
        c = sp.radsimp(sp.simplify(sp.expand(ser.coeff(w, n+1))))
        good = sp.simplify(c - tab(n)) == 0
        ok &= bool(good)
        print(name, n, 'series coeff of (v-i)^%d' % (n+1), c, '== table:', good)
    print('zero-order coefficient (should be 0):', sp.simplify(ser.coeff(w, 0)))
print('ALL tau table entries n<=6 match independent Taylor series:', ok)
print('tow2lib.check_tau(7):', T.check_tau(7))
