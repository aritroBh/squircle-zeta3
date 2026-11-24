#!/usr/bin/env python3
"""Step 05: Verify display (5) - C(z) contour identity"""
import mpmath as mp
mp.mp.dps = 60

def H_mp(s):
    return mp.beta(s/2, mp.mpf(3)/2) / 2

def K_mp(s):
    return 1 / (2 * mp.cos(mp.pi * s / 2))

def C_from_gamma(z):
    w = mp.e**(mp.pi*1j/4)
    def integrand(r):
        s = 0.5 + r * w
        return H_mp(s) * K_mp(s) * mp.exp(-s*z) * w
    I = mp.quad(integrand, [0, 0.5, 1, 2, 5, 10, mp.inf])
    return 1 - (2/mp.pi) * I.imag

def B(z):
    return 1 + mp.exp(-z) - mp.sqrt(1 + mp.exp(-2*z))

def C_direct(z):
    return mp.exp(z) * B(z)

print("Display (5) - C(z) contour identity:")
for z in [0, 0.1, 0.5, 1, 2, 5, 10]:
    c_direct = C_direct(z)
    c_gamma = C_from_gamma(z)
    print(f"z={z:4}: diff={mp.nstr(abs(c_direct-c_gamma), 15)}")
