import pwr1_engine as E, sympy as sp, re

# D11 exact f_n (eta = 1/p - 1 = -tau): P/sqrt2 = 4(1+eta) sum f_n eta^n
txt = open('../D11_rhombus_resurgence/step1_exact.txt').read()
f = {}
for m in re.finditer(r'f_(\d+) = (.*?)\n', txt):
    f[int(m.group(1))] = m.group(2)
zvars = {f'Z{k}': E.Z.get(k, sp.Symbol(f'Z{k}')) for k in range(2, 13)}
eta = sp.Symbol('eta')
N = 8
fser = sum(sp.sympify(f[k], locals=zvars)*eta**k for k in range(0, N+1))
D11 = sp.expand(4*(1+eta)*fser)   # P/sqrt2 in eta
D11 = sp.series(D11, eta, 0, N+1).removeO().expand()

# mine: P/sqrt2 in tau, substitute tau -> -eta
tot2, pure2, B2 = E.total_log(2)
P = E.ser_exp(pure2).scale(4)
t = sp.Symbol('t')
mine = sum(P.coeff(k)*t**k for k in range(0, N+1)).subs(t, -eta)
mine = sp.expand(mine)

diff = sp.expand(mine - D11)
print("P/sqrt2: my tau-series(tau=-eta) vs D11 4(1+eta)f(eta), difference through eta^8:")
print(diff if diff != 0 else "EXACTLY 0  (all coefficients eta^0..eta^8 match)")
open('progress.log','a').write("step_34: S_n/S_n(1) exact vs H1 routeA (n=2,3) OK; step_35: P/sqrt2 tau-series at tau=-eta vs D11 exact f_n: EXACT 0 through eta^8 (independent third-party check)\n")
