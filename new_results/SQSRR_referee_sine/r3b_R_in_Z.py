# R written as an exact Q-combination of depth-2 alternating sums Z and depth-1 values, before reduction;
# shows the coefficient of the free weight-4 element Z(3-,1-) (the only place Li4(1/2) could enter).
import sympy as sp
from r3_R_exact import sum_k, classes, T
from dshuffle import reduce_weight, ZS
R = 0
for eps, n, a, c, beta, rho, C in classes():
    R += sum_k(beta*rho, '1') + sum_k(beta*(-2*eps/a**3), 'Hk') \
         + sum_k(beta*2*eps*(3*n**2+2*eps*n+3)/(a**3*c**2), 'H2k') + sum_k(beta*2*eps*n/(a**2*c), 'O')
R = sp.expand(R)
Zs = sorted([s for s in R.free_symbols if str(s).startswith('Z_')], key=str)
print('T-part before reduction:', sp.factor(R.coeff(T, 1)), '|', R.coeff(T, 2))
print('R = sum c_Z Z + depth-1 part; coefficients:')
for z in Zs: print('   ', z, R.coeff(z))
print('   depth-1 part:', sp.expand(R.subs({z: 0 for z in Zs}).subs(T, 0)))
red4, free4 = reduce_weight(4)
fz = list(free4)[0]
tot = R.coeff(fz) + sum(R.coeff(z)*red4[z].coeff(fz) for z in red4 if z != fz)
print('total coefficient of the free element', fz, 'after substitution:', sp.simplify(tot))
