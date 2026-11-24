from mpmath import mp, mpf, mpc, quad, sqrt, exp, log, pi, cos, sin
mp.dps = 150
p = 1 + mpc(0, mpf('0.12')); e = 1/p
g = lambda u: exp(-u)*sqrt((e*exp(-(e-1)*u))**2 + (e*(1-exp(-u))**(e-1))**2)
bp_old = [log(2), 4, 10, 30, mp.inf]
bp_new = [log(2), 2, 4, 6, 8, 12, 16, 24, 40, 80, mp.inf]
r_old = 8*quad(g, bp_old, method='tanh-sinh')
r_new = 8*quad(g, bp_new, method='tanh-sinh')
print("old bp dps150:", mp.nstr(r_old, 50))
print("new bp dps150:", mp.nstr(r_new, 50))
print("diff:", mp.nstr(abs(r_old-r_new), 3))
# also Gauss-Legendre on u-variable, fine partition
bp3 = [log(2), 2, 3, 4, 5, 6, 8, 10, 13, 16, 20, 26, 34, 50, 100, 300, mp.inf]
r_gl = 8*quad(g, bp3, method='gauss-legendre', degree=10)
print("GL u-var:", mp.nstr(r_gl, 50))
print("diff new-GL:", mp.nstr(abs(r_new-r_gl), 3))
