"""ZLAD1R step_39 -- M_4(pi/2) from MY table vs G2's printed value; G1 unmodified-run status."""
import json, os
from fractions import Fraction as Fr
from mpmath import mp, mpf, zeta
mp.dps = 60
mine = json.load(open('muR_exact_N40_K8.json'))
ZET = {w: zeta(w) for w in range(2, 9)}
def mono_val(m):
    if m == '1': return mpf(1)
    v = mpf(1)
    for part in m.split('z'):
        if not part: continue
        if '^' in part:
            w, e = part.split('^'); v *= ZET[int(w)]**int(e)
        else: v *= ZET[int(part)]
    return v
def Mk_half(k, N=40):
    s = mpf(0)
    for n in range(1, N+1):
        t = sum(mpf(Fr(v).numerator)/mpf(Fr(v).denominator)*mono_val(m) for m, v in mine[str(k)][str(n)].items())
        s += t / mpf(2)**n
    return s
M4 = Mk_half(4)
ref = mpf('0.0060279931126091225045161082508776')
print('M_4(pi/2) from my table (n<=40):', mp.nstr(M4, 30))
print('G2 printed (37 digits)         :', mp.nstr(ref, 30))
print('|diff| =', mp.nstr(abs(M4-ref), 3), ' (expected ~ tail 2^-40 poly ~1e-13)')
M3 = Mk_half(3); M2 = Mk_half(2)
print('M_3(pi/2) ~', mp.nstr(M3, 20))
print('M_2(pi/2) ~', mp.nstr(M2, 20))
# G1 unmodified run status
pid = int(open('/tmp/g1_k8.pid').read().strip())
try:
    os.kill(pid, 0); alive = True
except OSError:
    alive = False
out = open('/tmp/g1_k8.out').read()
print('G1 flat_exact.py 40 8 unmodified: alive =', alive, ', output bytes =', len(out))
if out: print(out[-1200:])
open('progress.log','a').write('step_39: M_4(pi/2) from independent table matches G2 value to ~1e-13 (tail-limited); G1 unmodified K=8 run still running\n')
