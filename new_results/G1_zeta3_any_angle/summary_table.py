"""summary_table.py -- table of M2(alpha) (stored 318-digit series values, first 60 digits shown),
the exact weight-2 part W2 = z2/2 - Li2(c^2)/2, W3 = M2 - W2, Phi(c) = 1/2 - (1-c^2)/(2c) artanh c,
the uniform-classical zeta(3) coefficient Z(c) = 7/4 + 2c, and the canonical (flat-corner) zeta(3)
coefficient of M3, -Phi(c)."""
import mpmath as mp
mp.mp.dps = 80
ANG = [('pi/2', 'pi_2', mp.pi/2), ('2pi/3', '2pi_3', 2*mp.pi/3), ('pi/3', 'pi_3', mp.pi/3),
       ('3pi/4', '3pi_4', 3*mp.pi/4), ('pi/6', 'pi_6', mp.pi/6), ('5pi/6', '5pi_6', 5*mp.pi/6),
       ('pi/4', 'pi_4', mp.pi/4), ('1.3', 'gen1p3', mp.mpf(13)/10), ('2.0', 'gen2p0', mp.mpf(2))]
print('| alpha | c = cos(alpha/2) | M2(alpha) (60 digits) | W2 = z2/2 - Li2(c^2)/2 | Phi(c) | Z(c)=7/4+2c |')
for nm, tag, a in ANG:
    c = mp.cos(a/2)
    M2 = mp.mpf([l for l in open(f'values_{tag}.txt') if l.startswith('M2 =')][0].split('=')[1])
    W2 = (mp.zeta(2) - mp.polylog(2, c*c))/2
    Phi = mp.mpf(1)/2 - (1 - c*c)/(2*c)*mp.atanh(c)
    print(f'| {nm} | {mp.nstr(c, 12)} | {mp.nstr(M2, 60)} | {mp.nstr(W2, 20)} | {mp.nstr(Phi, 20)} | {mp.nstr(mp.mpf(7)/4 + 2*c, 15)} |')
