
import mpmath as mp
mp.mp.dps = 120
import json

def I_t(eta):
    # t-integral form
    f = lambda t: (t*(1-t))**(eta/2) * mp.sqrt(mp.cosh(eta*mp.log(t/(1-t))))
    return mp.quad(f, [0, mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('0.75'), 1])
def I_l(eta):
    # l-integral form: int_R w0(l)^{1+eta/2} 2 sqrt(cosh(eta l)) dl ; w0 = 1/(4 cosh^2(l/2))
    w0 = lambda l: 1/(4*mp.cosh(l/2)**2)
    f = lambda l: w0(l)**(1+eta/2)*2*mp.sqrt(mp.cosh(eta*l))
    return mp.quad(f, [-100, -20, 0, 20, 100])
def I_u(eta):
    f2 = lambda u: mp.exp(-u)*(2*mp.exp(-eta*u/2)*mp.sqrt(mp.cosh(eta*u)))*(1+mp.exp(-u))**(-2-eta)
    return mp.quad(f2, [0, 1, 5, 20, 100])
for eta in ['0.2','-0.2']:
    a, b2, c = I_t(mp.mpf(eta)), I_l(mp.mpf(eta)), I_u(mp.mpf(eta))
    print(f"eta={eta}")
    print("  t-form:", mp.nstr(a, 70))
    print("  l-form:", mp.nstr(b2, 70))
    print("  u-form:", mp.nstr(c, 70))
    print("  spreads:", mp.nstr(abs(a-b2),3), mp.nstr(abs(a-c),3))
