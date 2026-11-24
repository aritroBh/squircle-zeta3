# s05: INDEPENDENT totals P^(k)(1) by direct polar perimeter + polynomial interpolation in p
# (two stencils), compared with (a) the radial-projection integrals and (b) the exact closed forms.
import sys, json, pickle, mpmath as mp, sympy as sp
from edge_num import *
mp.mp.dps = 34
r2 = mp.sqrt(2)
cases = {
  'H_power':  ([mp.mpf(0), mp.pi/2, mp.pi/4], [1, 1, r2], [0, 0, mp.log(r2)]),
  'H_fixed':  ([mp.mpf(0), mp.pi/2, mp.pi/4], [1, 1, r2], [0, 0, 0]),
  'hexA_w123': ([mp.mpf(0), -mp.pi/3, -2*mp.pi/3], [1, 2, 3], [0, 0, 0]),
  'hexB_angles': ([mp.mpf(0), mp.mpf(1), mp.mpf('2.2')], [1, 1, 1], [0, 0, 0]),
}
cf = pickle.load(open('s04_closed_forms.pkl', 'rb'))
z2s, z3s = sp.symbols('z2 z3')
res = {}
for name, (bs, ws, lams) in cases.items():
    tp = P_coeffs_t(bs, ws, lams)
    tp = [mp.factorial(k)*tp[k] for k in range(4)]
    f1 = fd_derivs(bs, ws, lams, h=mp.mpf('0.004'), K=7)
    f2 = fd_derivs(bs, ws, lams, h=mp.mpf('0.003'), K=8)
    print(name)
    for k in range(4):
        line = f"  P^({k})(1): t-integral {mp.nstr(tp[k], 22)}  |FD1-t| {mp.nstr(abs(f1[k]-tp[k]), 3)}  |FD2-t| {mp.nstr(abs(f2[k]-tp[k]), 3)}"
        if name in cf:
            ex = sp.sympify(cf[name][k]).subs({z2s: sp.pi**2/6, z3s: sp.zeta(3)})
            v = mp.mpf(str(sp.N(ex, 40)))
            line += f"  |closed-t| {mp.nstr(abs(v - tp[k]), 3)}  |closed-FD2| {mp.nstr(abs(v - f2[k]), 3)}"
        print(line)
    res[name] = dict(t=[str(x) for x in tp], fd1=[str(x) for x in f1], fd2=[str(x) for x in f2])
json.dump(res, open('s05_fd.json', 'w'), indent=1)
