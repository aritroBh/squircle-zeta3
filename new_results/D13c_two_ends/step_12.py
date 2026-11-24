
import sympy as sp

# Symbolic saddle evaluation of  I = int_0^inf e^{i pi s/2} (s/(1+s))^N H(s) ds/s
# rotated ray s = e^{i pi/4} x ; phase phi(x) = i pi e^{i pi/4} x/2 + N log(s/(1+s))
x, N = sp.symbols('x N', positive=True)
i = sp.I
s = sp.exp(i*sp.pi/4)*x
phi = i*sp.pi*s/2 + N*sp.log(s/(1+s))
# exact saddle: dphi/dx = 0  ->  s(1+s) = 2 i N/pi  (in s); solve exactly
ss = sp.symbols('s', positive=False)
sols = sp.solve(sp.Eq(ss*(1+ss), 2*i*N/sp.pi), ss)
print("exact saddles in s:", sols)
s_star = sols[0]  # pick the one near + e^{i pi/4} sqrt(2 N/pi)
print("series of s_star:", sp.series(sols[0], N, sp.oo, 2))
print("series of s_star2:", sp.series(sols[1], N, sp.oo, 2))

# phase at saddle (large N expansion)
xstar = sp.sqrt(2*N/sp.pi)
phi_ser = sp.series(phi.subs(x, xstar), N, sp.oo, 3).removeO().expand()
print("phi(x*) series:", sp.simplify(phi_ser))
# expected -sqrt(pi N)(1-i) + const...
print("numeric check at N=400:", sp.N(phi.subs(x,xstar).subs(N,400),12),
      " vs -sqrt(pi N)(1-i) =", sp.N((-sp.sqrt(sp.pi*N)*(1-i)).subs(N,400),12))

# second derivative
phi2 = sp.diff(phi, x, 2)
phi2_star = sp.series(phi2.subs(x, xstar), N, sp.oo, 2).removeO().expand()
print("phi''(x*) series:", sp.simplify(phi2_star))
