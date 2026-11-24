"""squircle_check.py -- M2(pi/2) from corner_series (318 digits) vs the squircle note's 1200-digit M2 (method A)."""
import mpmath as mp
mp.mp.dps = 400
sq = mp.mpf(open('./paper/squircle_note/outputs/M2_methodA_1200.txt').read().split()[0])
ours = mp.mpf([l for l in open('values_pi_2.txt') if l.startswith('M2 =')][0].split('=')[1])
print('|M2(pi/2, series) - squircle M2| =', mp.nstr(abs(ours - sq), 3))
