# Planner check of D10c: boxed general-angle M1(alpha) vs direct quadrature of d/du of the corner integral at u = 0.
import sys, mpmath as mp
sys.path.insert(0, "research/08_theory/D10c_corner_general")
mp.mp.dps = 45
import importlib.util
src = open("research/08_theory/D10c_corner_general/step_19.py").read().split('print("expanded')[0]
src = src.replace("from corner_series import corner_coeffs", "").replace("mp.mp.dps = 90", "")
ns = {}; exec(src, ns)
def M1_quad(a):
    ca = mp.cos(a)
    def f(g):
        S0 = 1 - 2*ca*g + g*g
        f0 = mp.sqrt(S0)/((1+g)*g)
        lg = mp.log(g)
        return f0*(-mp.log(1+g) + lg + (-lg*g*g + ca*lg*g)/S0) - lg/g
    return mp.quad(f, [0, mp.mpf(1)/4, mp.mpf(1)/2, 1])
for a in [mp.pi/2, 2*mp.pi/3, mp.pi/3, 3*mp.pi/4, mp.mpf('1.3'), mp.mpf('0.4'), mp.mpf('2.9')]:
    v, n, z = ns["M1_expanded"](a)
    q = M1_quad(a)
    print(mp.nstr(a, 8), "formula", mp.nstr(mp.re(v), 25), "quad", mp.nstr(q, 25), "diff", mp.nstr(abs(mp.re(v)-q), 3), "Im", mp.nstr(abs(mp.im(v)), 3))
