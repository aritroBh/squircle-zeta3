import sympy as sp, time
t0 = time.time()

# ============ Theorem 1 (master expansion) -> displayed J1,J2,J3 ============
W, mu = sp.symbols('W mu')
a, A, B = sp.symbols('a A B')          # f'(1), f''(1), f'''(1)
c = sp.Symbol('c')                      # cos theta
k2, k3, Xi = sp.symbols('k2 k3 Xi')

# A_cal(tau) = -log f(1-tau) = a*tau + A2 tau^2 + A3 tau^3 + ...
# f(1-tau) = 1 - a t + A t^2/2 - B t^3/6
t = sp.Symbol('t')
logf = -sp.log(1 - a*t + A*t**2/2 - B*t**3/6)
A2 = sp.expand(sp.series(logf, t, 0, 4).removeO()).coeff(t,2)
A3 = sp.expand(sp.series(logf, t, 0, 4).removeO()).coeff(t,3)
# Btilde(tau) = log(f'(1-tau)/a) = E1 tau + E2 tau^2
fp = a - A*t + B*t**2/2
bt = sp.expand(sp.series(sp.log(fp/a), t, 0, 3).removeO())
E1 = bt.coeff(t,1); E2 = bt.coeff(t,2)
print("A2 =", A2, " A3 =", A3, "\nE1 =", E1, " E2 =", E2)

# Theta(W,mu): A_cal(Theta) - mu*Btilde(Theta) = W ; Theta = sum th_ij W^i mu^j, i+j<=3
t10, t20, t11, t30, t21, t12 = sp.symbols('t10 t20 t11 t30 t21 t12')
Theta = t10*W + t20*W**2 + t11*W*mu + t30*W**3 + t21*W**2*mu + t12*W*mu**2
Acal = a*Theta + A2*Theta**2 + A3*Theta**3
Btil = E1*Theta + E2*Theta**2
eq = sp.expand(Acal - mu*Btil - W)
# solve order by order
sols = {}
for (i,j),var in [((1,0),t10),((2,0),t20),((1,1),t11),((3,0),t30),((2,1),t21),((1,2),t12)]:
    coeff = sp.expand(eq.subs(sols)).coeff(W,i).coeff(mu,j)
    sols[var] = sp.solve(coeff, var)[0]
th = {(1,0): sols[t10], (2,0): sols[t20], (1,1): sols[t11],
      (3,0): sols[t30], (2,1): sols[t21], (1,2): sols[t12]}
for k in sorted(th): print("theta", k, "=", sp.simplify(th[k]))

# check theta_{1j} = a^{-1}(-A/a^2)^j  (P3's general formula with b=A)
print("theta_{1j} formula check:", sp.simplify(th[(1,1)] - (-A/a**3))==0,
      sp.simplify(th[(1,2)] - (A**2/a**5))==0)
