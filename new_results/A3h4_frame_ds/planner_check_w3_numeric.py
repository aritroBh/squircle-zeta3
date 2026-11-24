"""Planner check: Codex's weight-3 frame rows against INDEPENDENT 150-digit values (A3f6 eval8, zeta(3)-anchored)."""
import sys, pickle
from pathlib import Path
HERE = Path(__file__).resolve().parent; THEORY = HERE.parent
sys.path.insert(0, str(THEORY / 'A3c_claude_elementary/Q'))
import cmzv  # makes qsqrt2 importable for the pickles
import mpmath as mp
from fractions import Fraction as F
mp.mp.dps = 150
vals = pickle.load(open(THEORY / 'A3f6_level16/vals8_150digit.pickle', 'rb'))
rows = pickle.load(open(HERE / 'validated_w3.pickle', 'rb'))['rows']
worst = mp.mpf(0); miss = 0
for row in rows:
    s = mp.mpc(0)
    for w, c in row:
        z = vals.get(tuple(w))
        if z is None: miss += 1; continue
        c = F(c); s += (mp.mpf(c.numerator) / c.denominator) * mp.mpc(mp.mpf(z[0]), mp.mpf(z[1]))
    worst = max(worst, abs(s))
print('w3 frame rows', len(rows), 'missing words', miss, 'worst |row value| =', mp.nstr(worst, 3))
