"""r01: claim 1 by hand + numerically, plus validation of the referee evaluator against closed forms."""
from rcommon import *
import mpmath as mp, pickle
mp.mp.dps = 40
L = 'r01.log'
z = [mp.exp(2j * mp.pi * k / 8) for k in range(8)]
Hx = {j: -mp.log(1 - z[j]) for j in range(1, 8)}       # H(j) = int_0^1 z^j dt/(1 - z^j t) = -log(1 - z^j)
for j in range(1, 8):
    log(f'H({j}) closed form {mp.nstr(Hx[j], 15)}   referee quadrature {H((j,)):.15g}   diff {abs(complex(Hx[j]) - H((j,))):.1e}', L)
r = {(1,): 2, (7,): -2, (2,): -3, (6,): 3}
for u in (1, 7, 3, 5):
    ru = actv(r, u)
    v = mp.fsum(c * Hx[w[0]] for w, c in ru.items())
    log(f'u={u}: image {sorted(ru.items())}  value {mp.nstr(v, 20)}  (2*pi*i = {mp.nstr(2j*mp.pi, 12)})', L)
# by hand: H(j) - H(-j) = log(-z^-j) = i*pi*(4-j)/4 for j=1..7 (principal branch), so
# 2(H1-H7) - 3(H2-H6) = 2*(3i pi/4) - 3*(i pi/2) = 0 ;  u=3 image 2(H3-H5) + 3(H2-H6) = 2*(i pi/4) + 3*(i pi/2) = 2 pi i
# base rows of lift1 (A3c weight1.pickle selected) and their images: decomposition u.b = rho + c*e, e = H2 - H6
cols, rows, sel = pickle.load(open(T + '/A3c_claude_elementary/Q/weight1.pickle', 'rb'))
B1 = [{tuple(w): F(c) for w, c in row} for row in sel]
for i, b in enumerate(B1):
    vb = mp.fsum(c * Hx[w[0]] for w, c in b.items())
    out = []
    for u in (7, 3, 5):
        bu = actv(b, u); vu = mp.fsum(c * Hx[w[0]] for w, c in bu.items())
        # e-coefficient c = value / (i pi/2) must be rational; check that bu - c*e lies in span(B1) exactly
        cc = F(int(mp.nint(mp.im(vu) / (mp.pi / 2))))
        assert abs(vu - cc * 1j * mp.pi / 2) < mp.mpf(10) ** -30
        rest = dict(bu); rest[(2,)] = rest.get((2,), 0) - cc; rest[(6,)] = rest.get((6,), 0) + cc
        # exact solve in span(B1) over the 7 weight-1 words
        import sympy as sp
        Mx = sp.Matrix([[b2.get((j,), 0) for j in range(1, 8)] for b2 in B1]).T
        sol = Mx.gauss_jordan_solve(sp.Matrix([rest.get((j,), 0) for j in range(1, 8)]))[0]
        out.append(f'u={u}: value={mp.nstr(vu, 8)}, e-coeff={cc}, rho in span(B1): {list(sol)}')
    log(f'B1[{i}] = {sorted(b.items())}, value {mp.nstr(vb, 5)}; ' + '; '.join(out), L)
