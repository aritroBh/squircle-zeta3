import sys
sys.path.insert(0,'.')
from oct9c_common import *
D = load_all()
lam = D['lam']
mp.mp.dps = 60
def Li21val(a, c):
    a, c = mp.mpf(a), mp.mpf(c)
    if abs(a) <= 1 and abs(a*c) < 1:
        N = int((mp.mp.dps+20)/(-mp.log10(max(abs(a), abs(a*c))))) + 50
        return mp.fsum([a**n * c**m / (mp.mpf(n)**2 * m) for n in range(2, N) for m in range(1, n)])
    F = mp.quad(lambda t: mp.log(t)*mp.log(1-a*c*t)/(1-a*t), [0,1])
    return a*F
v_a = Li21val(sp.Rational(1,5), sp.Rational(1,5))
v_b = Li21val(sp.Rational(1,5), sp.Rational(5,29))
print("cache e3L3:", lam['e3L3'][:60])
print("Li21(1/5,1/5):", mp.nstr(v_a, 55))
print("Li21(1/5,5/29):", mp.nstr(v_b, 55))
v_c = Li21val(sp.Rational(1,25), sp.Rational(5,1))
v_d = Li21val(sp.Rational(1,29), sp.Rational(29,5))
print("cache e3L4:", lam['e3L4'][:60])
print("Li21(1/25,5):", mp.nstr(v_c, 55))
print("Li21(1/29,29/5):", mp.nstr(v_d, 55))
