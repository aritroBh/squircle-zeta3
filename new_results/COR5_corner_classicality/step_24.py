from fractions import Fraction as Fr
import json

# S-unit equation x + y = 1, x,y in {+- 2^a 3^b}. Brute force in a box, then prove box suffices.
def sunits(B):
    vals = set()
    for a in range(-B, B+1):
        for b in range(-B, B+1):
            v = Fr(2)**a * Fr(3)**b
            vals.add(v); vals.add(-v)
    return vals

B = 14
S = sunits(B)
sol = set()
one = Fr(1)
# enumerate as x + y = 1
xs = {}
for x in S:
    y = one - x
    if y in S:
        sol.add((x, y))
print("ordered solutions found in box:", len(sol))
# unique x values
xsorted = sorted(set(x for x, y in sol))
print("distinct x values:", len(xsorted))
for x in xsorted:
    # factor x
    def fac(v):
        s = 1 if v > 0 else -1
        v = abs(v)
        a = 0
        while v % 2 == 0: v //= 2; a += 1
        b = 0
        while v % 3 == 0: v //= 3; b += 1
        # handle denominators
        num, den = v.numerator, v.denominator
        while den % 2 == 0: den //= 2; a -= 1
        while den % 3 == 0: den //= 3; b -= 1
        assert num == 1 and den == 1, v
        return s, a, b
    sx, ax, bx = fac(x)
    sy, ay, by = fac(1 - x)
    print(f"x = {str(x):>8} = {'-' if sx<0 else ''}2^{ax}*3^{bx} ; 1-x = {str(1-x):>8} = {'-' if sy<0 else ''}2^{ay}*3^{by}")
