"""step_13.py -- GENERAL-ANGLE THEOREM, part A fix: direct construction of the two
rational coefficients of the t-integrand, exact comparison to Astar, Bstar."""
import sympy as sp

t, c = sp.symbols('t c')
ca = 2*c**2 - 1
c2a = sp.expand(2*ca**2 - 1)
r1, r2 = 2*c - 1, -1 - 2*c
P = t**2 - 2*ca*t + 1
D = (t**2 - 1)*(t - ca)*(t - r1)*(t - r2)
g = (t**2 - 1)/(2*(t - ca))
dgdt = P/(2*(t - ca)**2)                      # proved in step_12
A0 = P/((t - r1)*(t - r2))                    # sqrt(S0)/(1+g) on path (proved)
S0v = P**2/(4*(t - ca)**2)

# M1 integrand (1/g){ A0[(1-ca g)/S0 log g - log(1+g)] - log g } dg/dt
X_direct = sp.cancel((1/g)*(A0*(1 - ca*g)/S0v - 1)*dgdt)   # coeff of log g
Y_direct = sp.cancel((1/g)*(-A0)*dgdt)                     # coeff of log(1+g)

Astar = P**2/D - P/((t**2 - 1)*(t - ca)) - (t**2 - 2*ca*t + c2a)/((t - ca)*(t - r1)*(t - r2))
Bstar = -P**2/D
print("simplify(X_direct - Astar) =", sp.simplify(X_direct - Astar))
print("simplify(Y_direct - Bstar) =", sp.simplify(Y_direct - Bstar))

# compact single-fraction form of Astar: numerator N over D
N = sp.expand(2*(t - ca)*(-ca*t**2 + 2*t - ca) - (t - r1)*(t - r2)*P)
print("Astar = N/D with N =", sp.factor(N))
print("check N: simplify(Astar - N/D) =", sp.simplify(Astar - N/D))
print("Bstar =", sp.factor(Bstar))
