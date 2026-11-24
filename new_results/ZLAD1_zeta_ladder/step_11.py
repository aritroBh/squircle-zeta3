import os
os.makedirs('research/08_theory/ZLAD1_zeta_ladder', exist_ok=True)
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
# progress.log + result.md bootstrap
open('progress.log','a').write('step_01: explored refs (flat_exact.py, G2 result.md, main.tex thm:flat); started ZLAD1\n')

# Verify the master identity:  M(u;alpha) = sum_j K_j(c^2) B(1+(j-1)v, 1-jv), v=1-u,
# where K_j(x) = sum_{n=1}^j (-2*Catalan(n-1)) (-1)^{j-n} binom(n+j-2, 2n-2) x^n
# by comparing [u^k] coefficients with flat_exact.py's mu_{k,n} for small k,n.
from fractions import Fraction as Fr
from math import comb, gamma
import sympy as sp

def catalan(m): return comb(2*m, m)//(m+1)

def K_poly(j, x):  # symbolic in x
    return sum((-2*catalan(n-1))*(-1)**(j-n)*comb(n+j-2, 2*n-2)*x**n for n in range(1, j+1))

x = sp.symbols('x')
for j in range(1, 6):
    print('K_%d ='%j, sp.factor(K_poly(j, x)))
