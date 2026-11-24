"""Planner: L(chi8,5) to 1800+ digits by mpmath Hurwitz zeta, for M5b's W6 value stage (its own recomputation failed).
L(chi8,5) = 8^-5 [zeta(5,1/8) + zeta(5,7/8) - zeta(5,3/8) - zeta(5,5/8)], chi8(1)=chi8(7)=1, chi8(3)=chi8(5)=-1.
Two working precisions; also compared with M4e's 1650-digit Lchi8_5.json."""
import json
import mpmath as mp
from pathlib import Path
def L(d):
    mp.mp.dps = d
    h = [mp.zeta(5, mp.mpf(k) / 8) for k in (1, 3, 5, 7)]
    return (h[0] + h[3] - h[1] - h[2]) / mp.mpf(8) ** 5
a = L(1860); b = L(1900)
mp.mp.dps = 1900
saved = mp.mpf(json.loads(Path('../Lchi8_5.json').read_text())['value'])
print('two-precision diff:', mp.nstr(a - b, 5), ' vs M4e 1650-digit file:', mp.nstr(b - saved, 5))
Path('../Lchi8_5_1800_planner.json').write_text(json.dumps({'source': 'planner Hurwitz zeta, dps 1900 (checked vs 1860)',
    'digits': 1800, 'value': mp.nstr(b, 1810)}))
print('wrote Lchi8_5_1800_planner.json')
