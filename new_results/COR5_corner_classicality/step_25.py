from fractions import Fraction as Fr

def fac(v):
    s = 1 if v > 0 else -1
    num, den = abs(v).numerator, abs(v).denominator
    a = 0
    while num % 2 == 0: num //= 2; a += 1
    while den % 2 == 0: den //= 2; a -= 1
    b = 0
    while num % 3 == 0: num //= 3; b += 1
    while den % 3 == 0: den //= 3; b -= 1
    assert num == 1 and den == 1, v
    return s, a, b

def sunits(B):
    vals = set()
    for a in range(-B, B+1):
        for b in range(-B, B+1):
            v = Fr(2)**a * Fr(3)**b
            vals.add(v); vals.add(-v)
    return vals

B = 16
S = sunits(B)
one = Fr(1)
sol = set()
for x in S:
    if one - x in S:
        sol.add(x)
xsorted = sorted(sol)
print("distinct x with x, 1-x both {2,3}-units:", len(xsorted))
out = []
for x in xsorted:
    sx, ax, bx = fac(x); sy, ay, by = fac(1 - x)
    out.append((x, (sx,ax,bx), (sy,ay,by)))
    print(f"x = {str(x):>8} = ({sx:+})2^({ax})3^({bx}) ;  1-x = {str(1-x):>8} = ({sy:+})2^({ay})3^({by})")
import json
json.dump([[str(x), list(fx), list(fy)] for x, fx, fy in out], open("sunit23_solutions.json","w"), indent=1)
print("saved sunit23_solutions.json")
