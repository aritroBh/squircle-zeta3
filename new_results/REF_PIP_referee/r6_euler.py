"""Referee script 6: classical Euler sum used in the proof sketch, and J_{2,0}, at 60 digits."""
import mpmath as mp
mp.mp.dps = 60
# sum_{n>=1} (-1)^{n+1} H_n / n^2 via the integral  -int_0^1 log(1-x)... use direct alternating-series acceleration
S = mp.nsum(lambda n: (-1)**(int(n) + 1) * mp.harmonic(n) / n**2, [1, mp.inf], method='alternating')
print('sum (-1)^{n+1} H_n/n^2 - 5 zeta(3)/8 =', mp.nstr(S - 5 * mp.zeta(3) / 8, 3))
J20 = mp.quad(lambda v: mp.log1p(v)**2 / v, [0, 1])
print('int_0^1 log^2(1+v)/v dv - zeta(3)/4 =', mp.nstr(J20 - mp.zeta(3) / 4, 3))
