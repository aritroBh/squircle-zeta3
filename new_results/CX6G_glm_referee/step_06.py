#!/usr/bin/env python3
"""Step 06: Verify C_star constant"""
import mpmath as mp
from flint import arb, ctx
mp.mp.dps = 80
ctx.dps = 90

p = arb.pi()
a = p / (2 * arb(2).sqrt())
x = 2 * (p / 4).exp()
K = 64 / (p * a) * (a / 2 + p / 4).exp()
T = 32 / (3 * p * a) * (3 * a / 2 + 3 * p / 4).exp()
C_star = 2 + T + K * ((p / 4).exp() + x.exp() + (arb(1).exp() * x).exp())

print(f"C_star = {C_star.str(30)}")
print(f"C_star < 10,600,000: {C_star < 10600000}")
