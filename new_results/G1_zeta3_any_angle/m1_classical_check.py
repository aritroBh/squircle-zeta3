"""m1_classical_check.py -- classical form of the general-angle M1 (from the 8-term HPL formula of hpl_fit.py, k=1)
   M1 = z2 - log^2 2 - log2*M0 - sum_{e=+-1} (1-e c) [Li2(e c) + Li2((1-e c)/2)],
   M0 = -(1+c)log(1+c) - (1-c)log(1-c),  c = cos(alpha/2);  checked against stored series M1 at 9 angles."""
import mpmath as mp
mp.mp.dps = 320
ANG = {'pi_2': mp.pi/2, '2pi_3': 2*mp.pi/3, 'pi_3': mp.pi/3, '3pi_4': 3*mp.pi/4, 'pi_6': mp.pi/6,
       '5pi_6': 5*mp.pi/6, 'pi_4': mp.pi/4, 'gen1p3': mp.mpf(13)/10, 'gen2p0': mp.mpf(2)}
for tag, a in ANG.items():
    c = mp.cos(a/2)
    M0 = -(1 + c)*mp.log(1 + c) - (1 - c)*mp.log(1 - c)
    M1 = mp.zeta(2) - mp.log(2)**2 - mp.log(2)*M0 - sum((1 - e*c)*(mp.polylog(2, e*c) + mp.polylog(2, (1 - e*c)/2)) for e in (1, -1))
    ref = mp.mpf([l for l in open(f'values_{tag}.txt') if l.startswith('M1 =')][0].split('=')[1])
    print(f'{tag:>7}: |M1 classical - series| = {mp.nstr(abs(M1 - ref), 3)}')
