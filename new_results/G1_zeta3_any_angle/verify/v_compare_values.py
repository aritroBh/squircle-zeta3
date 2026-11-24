"""v_compare_values.py -- compare verifier's Cauchy-contour M_k (two runs: R=0.4,N=128,dps55 and R=0.33,N=160,dps60)
with each other and with the solver's 318-digit values_*.txt."""
import mpmath as mp, os
mp.mp.dps = 80
def rd(fn):
    d = {}
    for l in open(fn):
        if l.startswith('M'):
            k, v = l.split('=')
            d[int(k.strip()[1])] = mp.mpf(v)
    return d
tags = ['pi_2','2pi_3','pi_3','3pi_4','pi_6','5pi_6','pi_4','gen1p3','gen2p0','new0p3','new0p9','new2p6','new3p0']
worst_AB = 0; worst_sol = 0
for t in tags:
    A = rd(f'vc_{t}.txt'); B = rd(f'runB/vc_{t}.txt')
    line = f'{t:>7}: |A-B| k=0..3: ' + ' '.join(mp.nstr(abs(A[k]-B[k]), 3) for k in range(4))
    worst_AB = max(worst_AB, max(abs(A[k]-B[k]) for k in range(4)))
    if os.path.exists(f'values_{t}.txt'):
        S = rd(f'values_{t}.txt')
        line += '   |B - solver| k=0..3: ' + ' '.join(mp.nstr(abs(S[k]-B[k]), 3) for k in range(4))
        worst_sol = max(worst_sol, max(abs(S[k]-B[k]) for k in range(4)))
    print(line)
sq = mp.mpf(open('./paper/squircle_note/outputs/M2_methodA_1200.txt').read().split()[0])
print('|M2(pi/2) Cauchy run B - squircle 1200-digit M2| =', mp.nstr(abs(rd('runB/vc_pi_2.txt')[2] - sq), 3))
print('worst |A-B| =', mp.nstr(worst_AB, 3), '; worst |B - solver values| =', mp.nstr(worst_sol, 3))
