#!/usr/bin/env python3
"""Step 04: Verify display (4) - H(s) bounds"""
import mpmath as mp
mp.mp.dps = 50

def H(s):
    return mp.beta(s/2, mp.mpf(3)/2) / 2

print("Display (4) - H(s) bounds:")
for x in [0.1, 0.5, 1, 2, 5, 10]:
    h_val = abs(H(x))
    bound1 = 1/x
    bound2 = mp.sqrt(mp.pi/2) * x**(-1.5)
    print("x={:.2f}: |H|/bound1={:.6f}, |H|/bound2={:.6f}".format(x, float(h_val/bound1), float(h_val/bound2)))
