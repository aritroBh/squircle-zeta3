
import sympy as sp

# B(z) = 1 - e^z/(1+sqrt(1+e^{2z})) near z_s = i*pi/2 (s=+1), "above-cut" branch.
# Above the cut: zeta = r + i eps, 1+e^{2z} ~ -2 zeta has Im part -2 eps < 0,
# so sqrt(1+e^{2z}) = -i*sqrt(2*zeta)*(series in zeta).  Work with t = sqrt(zeta).
t = sp.symbols('t')
zeta = t**2
s = sp.Integer(1)
zs = s*sp.pi*sp.I/2

# 1+e^{2(zs+zeta)} = 1 - e^{2 zeta}
inner = 1 - sp.exp(2*zeta)         # = -2 zeta - 2 zeta^2 - ...
# factor: inner = -2 zeta * Q, Q = 1 + zeta + 2/3 zeta^2 + 1/3 zeta^3 + ...
Q = sp.series(-inner/(2*zeta), t, 0, 8).removeO()
Q = sp.expand(Q)
# sqrt(inner): above-cut branch = -I*sqrt(2)*t*sqrt(Q)
sqrtQ = sp.series(sp.sqrt(Q), t, 0, 8).removeO().expand()
w = -sp.I*sp.sqrt(2)*t*sqrtQ                      # = sqrt(1+e^{2z})
inv = sp.series(1/(1+w), t, 0, 9).removeO().expand()
ez = sp.exp(zs)*sp.series(sp.exp(zeta), t, 0, 9).removeO().expand()
Bser = sp.expand(1 - sp.expand(ez*inv).series(t,0,9).removeO())
Bser = sp.collect(sp.expand(Bser), t)
print("B(z_s+zeta) series in t=zeta^{1/2} (s=+1, above cut):")
for n in range(0, 9):
    c = sp.simplify(Bser.coeff(t, n))
    print(f"  t^{n} (zeta^{n}/2): {c}")
# ratios
g1 = Bser.coeff(t,1); g3 = Bser.coeff(t,3); g5 = Bser.coeff(t,5)
print("gamma (zeta^{1/2})      =", sp.simplify(g1), "  expect s*sqrt(2) =", sp.sqrt(2))
print("c3 = (zeta^{3/2})/gamma =", sp.simplify(g3/g1))
print("c5 = (zeta^{5/2})/gamma =", sp.simplify(g5/g1))
print("integer coeffs: d1 =", sp.simplify(Bser.coeff(t,2)), " d2 =", sp.simplify(Bser.coeff(t,4)),
      " d3 =", sp.simplify(Bser.coeff(t,6)))
