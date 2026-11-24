import sympy as sp
x = sp.Symbol('x', positive=True)
Li2 = lambda z: sp.polylog(2,z)
z2 = sp.Symbol('z2')
L1 = Li2(x)+Li2(x/(x-1)) + sp.log(1-x)**2/2
d = sp.expand(sp.diff(L1,x))
d = d.replace(lambda f: isinstance(f, sp.Function) and f.func==sp.polylog and f.args[0]==1,
              lambda f: -sp.log(1-f.args[1]))
# for 0<x<1: 1 - x/(x-1) = 1/(1-x) > 0, so log(1 - x/(x-1)) = -log(1-x); sympy form: log(-1/(x-1))
d = d.subs(sp.log(-1/(x-1)), -sp.log(1-x))
print("Landen d/dx (0<x<1):", sp.simplify(d))
# numeric pins
import mpmath as mp
mp.mp.dps = 60
z2v = mp.pi**2/6; z3v = mp.zeta(3)
L2 = lambda z: mp.polylog(2,z); L3 = lambda z: mp.polylog(3,z)
print("Euler at 0.3:", mp.nstr(L2(mp.mpf('0.3'))+L2(mp.mpf('0.7'))-(z2v-mp.log(mp.mpf('0.3'))*mp.log(mp.mpf('0.7'))),5))
print("Landen at 0.3:", mp.nstr(L2(mp.mpf('0.3'))+L2(mp.mpf('0.3')/(mp.mpf('0.3')-1))+mp.log(mp.mpf('0.7'))**2/2,5))
print("inversion at 2.0:", mp.nstr(L2(-2)+L2(-mp.mpf('0.5'))+z2v+mp.log(2)**2/2,5))
print("Li2(-1)+z2/2 (dup pin):", mp.nstr(L2(-1)+z2v/2,5))
s5 = mp.sqrt(5); phi=(1+s5)/2; xx = 1/phi**2
lhs = L3(xx)+L3(1-xx)+L3(xx/(xx-1))
rhs = z3v + z2v*mp.log(1-xx) - mp.log(1-xx)**2*mp.log(xx)/2 + mp.log(1-xx)**3/6
print("Li3 Landen at 1/phi^2:", mp.nstr(lhs-rhs,5))
lhs2 = L3(mp.mpf('0.3'))+L3(mp.mpf('0.7'))+L3(mp.mpf('0.3')/(mp.mpf('0.3')-1))
rhs2 = z3v + z2v*mp.log(mp.mpf('0.7')) - mp.log(mp.mpf('0.7'))**2*mp.log(mp.mpf('0.3'))/2 + mp.log(mp.mpf('0.7'))**3/6
print("Li3 Landen at 0.3:", mp.nstr(lhs2-rhs2,5))
