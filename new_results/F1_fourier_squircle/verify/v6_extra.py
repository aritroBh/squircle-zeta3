"""
VERIFIER V6: (a) classical integral quoted for the k^2-part of c_5: I_{2,2} = int_0^1 log^2 x log^2(1+x)/x dx
              = 2 zeta(2) zeta(3) - (29/8) zeta(5).
             (b) the solver's explicit beta_k formula at large k (k = 50, 200) vs direct quadrature
                 (the solver checked only k <= 6 directly).
             (c) sine transform used in the proof of int_0^inf (1/y - csch y)^2 = 2 log2 - 1.
"""
import mpmath as mp, importlib.util, sys, io, contextlib
mp.mp.dps = 40
z2, z3, z5 = mp.zeta(2), mp.zeta(3), mp.zeta(5)
I22 = mp.quad(lambda x: mp.log(x) ** 2 * mp.log1p(x) ** 2 / x, [0, mp.mpf(1)/2, 1])
print("(a) I_{2,2} - (2 z2 z3 - 29/8 z5) =", mp.nstr(I22 - (2 * z2 * z3 - mp.mpf(29) / 8 * z5), 3))

# (b) import beta_list from the solver's script text without running its heavy main part
src = open('research/10_frontier/F1_fourier_squircle/f1_04_circle_end.py').read()
start = src.index('def H(n):'); end = src.index('# check beta_k against direct quadrature')
ns = {'mp': mp}
mp.mp.dps = 40
ns['lamn'] = mp.log(2)
exec(src[start:end], ns)
bl = ns['beta_list'](200)
def F2(th):
    Cn = mp.cos(th) ** 2; an = mp.log(abs(mp.cos(th))); bn = mp.log(abs(mp.sin(th)))
    K1 = Cn * an + (1 - Cn) * bn; K2 = Cn * (1 - Cn) * (an - bn) ** 2
    return K2 + K1 ** 2 / 2
for k in [50, 200]:
    pts = mp.linspace(0, mp.pi / 2, 8 * k + 1)
    q = 4 / mp.pi * mp.quad(lambda t: F2(t) * mp.cos(4 * k * t), pts)
    print(f"(b) beta_{k}: formula {mp.nstr(bl[k-1], 25)}  quad {mp.nstr(q, 25)}  diff {mp.nstr(bl[k-1] - q, 3)}")

# (c) int_0^inf (1/y - csch y) sin(w y) dy = (pi/2)(1 - tanh(pi w/2))
for w in [mp.mpf('0.7'), mp.mpf(2)]:
    lhs = mp.quadosc(lambda y: (1 / y - 1 / mp.sinh(y)) * mp.sin(w * y), [0, mp.inf], omega=w)
    print(f"(c) w={w}: sine transform - (pi/2)(1-tanh(pi w/2)) = {mp.nstr(lhs - mp.pi/2*(1 - mp.tanh(mp.pi*w/2)), 3)}")
