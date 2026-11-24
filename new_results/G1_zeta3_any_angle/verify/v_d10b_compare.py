"""v_d10b_compare.py -- solver's 318-digit M2 at 3pi/4, pi/6 vs D10b's independent 340-digit file."""
import mpmath as mp, re
mp.mp.dps = 340
txt = open('research/08_theory/D10b_octagon/M_new_340.txt').read()
def d10b(k, name):
    m = re.search(rf'M{k}\({re.escape(name)}\) = (-?[0-9.]+)', txt)
    return mp.mpf(m.group(1))
def sol(tag, k):
    return [mp.mpf(l.split('=')[1]) for l in open(f'values_{tag}.txt') if l.startswith(f'M{k} =')][0]
for tag, name in [('3pi_4', '3pi/4'), ('pi_6', 'pi/6')]:
    print(tag, ' '.join(f'|dM{k}|={mp.nstr(abs(d10b(k, name) - sol(tag, k)), 3)}' for k in range(4)))
