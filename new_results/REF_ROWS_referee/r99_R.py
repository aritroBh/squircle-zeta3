# Isolate the unevaluated series R = sum_{n>=3} b_n' A_n from b_1'''(2) and try integer relations.
import sys
sys.argv = ['x', '50']
import r99_sine as E
from mpmath import mp, mpf, pi, log, zeta, gamma, sin, diff, nstr, pslq, catalan
mp.dps = 50
b1 = E.bder(1)
l2 = log(2)
N = lambda p: 2*sin(pi/p)/pi*gamma(3/p)*gamma(1-1/p)/gamma(1+2/p)
mp.dps = 70; N3 = diff(N, 2, 3); mp.dps = 50
sq = pi**4/768 + 5*pi**2/96 - mpf(5)/8
R = -(b1[3]-N3/2)/3 - b1[1]*b1[2] - l2*sq - pi**2*(1-zeta(3))/64
print("R = sum_{n>=3} b_n' A_n =", nstr(R, 40))
basis_names = ['R', '1', 'ln2', 'ln2^2', 'ln2^3', 'pi^2', 'pi^2 ln2', 'zeta3', 'pi^2 zeta3', 'pi^4', 'pi^4 ln2', 'G', 'pi G']
mp.dps = 44
basis = [R, 1, l2, l2**2, l2**3, pi**2, pi**2*l2, zeta(3), pi**2*zeta(3), pi**4, pi**4*l2, catalan, pi*catalan]
rel = pslq(basis, maxcoeff=10**6, maxsteps=10**6)
print("PSLQ (44 digits, maxcoeff 1e6):", rel and list(zip(basis_names, rel)))
