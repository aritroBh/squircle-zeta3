WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
src = open(WD+'cor9_lib.py').read()
old = src[src.index('# ---------- class number'):]
new = '''# ---------- class number via cycles of reduced indefinite primitive forms ----------
def h_plus_order(disc):
    """narrow class number of the quadratic order of discriminant disc > 0 (nonsquare):
    number of cycles of reduced primitive forms (a,b,c), b^2-4ac = disc."""
    sq = math.isqrt(disc)
    red = []
    for b in range(1, sq+1):
        if (disc - b*b) % 4: continue
        num = b*b - disc   # = 4ac < 0
        for a in range(1, sq+1):
            if num % (4*a): continue
            c = num//(4*a)
            if c == 0: continue
            if math.gcd(math.gcd(a, b), abs(c)) != 1: continue
            if abs(sq - 2*abs(a)) < b:
                red.append((a, b, c))
    red += [(-a, b, -c) for a, b, c in red]
    def rho(f):
        a, b, c = f
        t = (b + sq)//(2*c)          # python floor division
        if c < 0: t += 1
        return (c, -b + 2*c*t, a - b*t + c*t*t)
    seen = {}
    ncyc = 0
    for f in red:
        if f in seen: continue
        ncyc += 1
        g = f
        while g not in seen:
            seen[g] = ncyc
            g = rho(g)
    return ncyc

def class_number(D):
    """ordinary class number of Q(sqrt D), D squarefree > 1."""
    disc = D if D % 4 == 1 else 4*D
    hp = h_plus_order(disc)
    x, y = fund_unit(D)
    N = x*x - D*y*y
    return hp if N == -1 else hp//2
'''
src = src[:src.index('# ---------- class number')] + new
open(WD+'cor9_lib.py','w').write(src)
import sys; sys.path.insert(0, WD)
import importlib, cor9_lib; importlib.reload(cor9_lib)
from cor9_lib import class_number, fund_unit
tests = {2:1,3:1,5:1,6:1,7:1,10:2,11:1,13:1,14:1,15:2,17:1,19:1,21:1,22:1,23:1,26:2,29:1,30:2,31:1,33:1,
         34:2,35:2,37:1,38:1,39:2,41:1,42:2,43:1,46:1,47:1,53:1,57:1,59:1,61:1,62:1,67:1,69:1,71:1,73:1,
         77:1,79:3,82:4,83:1,86:1,89:1,93:1,94:1,97:1}
bad = {D: (class_number(D), h) for D, h in tests.items() if class_number(D) != h}
print("mismatches:", bad if bad else "NONE - all", len(tests), "class numbers correct")
