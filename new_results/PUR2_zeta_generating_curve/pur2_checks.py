"""PUR2 small checks: Taylor coefficients of the closed forms vs the symbolic table; the elementary (B) integral;
the Beta lemma at complex alpha."""
import mpmath as mp
mp.mp.dps = 50
z = mp.zeta
def fam(e):  # perimeter coefficients tau^2..tau^6 from pur2_symbolic_output.txt (e = ab)
    return [2*z(2), 2*z(3), (2 - mp.mpf(7)/2*e)*z(4), 2*(1-e)*z(5) - 2*e*z(2)*z(3),
            (2 - 13*e/2 + 31*e**2/8)*z(6) - e*z(3)**2]
for name, a, b in [('c=1/2', mp.mpf(3)/4, mp.mpf(1)/4), ('kappa^2=9/7', mp.mpc(1, mp.sqrt(mp.mpf(9)/7))/2, mp.mpc(1, -mp.sqrt(mp.mpf(9)/7))/2)]:
    f = lambda t: 8 + 2/(a*b)*(1 - mp.gamma(1-a*t)*mp.gamma(1-b*t)/mp.gamma(1-t))
    tc = mp.taylor(f, 0, 6)
    pred = fam(a*b)
    print(name, "ab =", mp.nstr(a*b, 10), "max |taylor - table| (tau^2..tau^6):", mp.nstr(max(abs(tc[n+2]-pred[n]) for n in range(5)), 3),
          "; tau^4 coefficient:", mp.nstr(tc[4], 5))
e = mp.mpf(4)/7
print("ab=4/7: tau^5 = 6/7 z5 - 8/7 z2 z3 =", mp.nstr(fam(e)[3], 20), "; tau^6 = -22/49 z6 - 4/7 z3^2 =", mp.nstr(fam(e)[4], 20),
      "check", mp.nstr(fam(e)[4] - (-mp.mpf(22)/49*z(6) - mp.mpf(4)/7*z(3)**2), 3))
# (B) elementary identity: int_0^1 (x^{-tau-1} + x^{-1})((1+x)^tau - 1) dx = -gamma - psi(1-tau)
for tau in [mp.mpf('0.3'), mp.mpf('0.77')]:
    kk = 1/(1-tau)   # x = y^kk removes the x^{-tau} endpoint singularity
    I = mp.quad(lambda y: kk*y**(kk-1)*(y**(-kk*(tau+1)) + y**(-kk))*mp.expm1(tau*mp.log1p(y**kk)), [0, 1])
    print("tau=%s: integral - (-gamma - psi(1-tau)) = %s" % (tau, mp.nstr(I + mp.euler + mp.digamma(1-tau), 3)))
# Beta lemma at complex alpha: int_0^inf (e^{2 al s}+e^{2 be s})((1+e^{-2s})^tau - 1) ds = (1/2)[G(-al)G(-be)/G(-tau) + 1/al + 1/be]
tau = mp.mpf('0.6'); al = mp.mpc('0.3', '0.7'); be = tau - al
L = mp.quad(lambda s: (mp.exp(2*al*s) + mp.exp(2*be*s))*((1+mp.exp(-2*s))**tau - 1), [0, 1, 5, 20, 80])
Rr = (mp.gamma(-al)*mp.gamma(-be)/mp.gamma(-tau) + 1/al + 1/be)/2
print("Beta lemma, tau=0.6, alpha=0.3+0.7i: lhs - rhs =", mp.nstr(L - Rr, 3))
