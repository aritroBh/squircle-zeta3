src = '''"""mymaps.py -- referee's own cyclotomic arithmetic + Mobius map counting (2N theorem)."""
from fractions import Fraction as Fr
import sympy as sp

class MyCyclo:
    def __init__(self, N):
        self.N = N
        x = sp.symbols('x')
        phi = sp.Poly(sp.cyclotomic_poly(N, x))
        self.d = phi.degree()
        self.mod = [Fr(int(c)) for c in phi.all_coeffs()]
        self.zero = (Fr(0),) * self.d
        zeta = (Fr(0), Fr(1)) + (Fr(0),) * (self.d - 2) if self.d > 1 else (Fr(1),)
        self.roots = []
        r = (Fr(1),) + (Fr(0),) * (self.d - 1)
        for k in range(N):
            self.roots.append(r)
            r = self.mul(r, zeta)
        assert r == (Fr(1),) + (Fr(0),) * (self.d - 1)
    def add(self, x, y): return tuple(a + b for a, b in zip(x, y))
    def sub(self, x, y): return tuple(a - b for a, b in zip(x, y))
    def scal(self, x, s): return tuple(a * s for a in x)
    def mul(self, x, y):
        d = self.d
        tmp = [Fr(0)] * (2 * d - 1)
        for i, a in enumerate(x):
            if a:
                for j, b in enumerate(y):
                    if b: tmp[i + j] += a * b
        for k in range(2 * d - 2, d - 1, -1):
            c = tmp[k]
            if c:
                for i in range(1, d + 1):
                    tmp[k - i] -= c * self.mod[i]
                tmp[k] = Fr(0)
        return tuple(tmp[:d])

def prime_for(N, kmin=50):
    k = kmin
    while True:
        p = k * N + 1
        if sp.isprime(p): return int(p)
        k += 1

def map_count(N):
    p = prime_for(N)
    g = int(sp.primitive_root(p)); z = pow(g, (p - 1) // N, p)
    roots = [pow(z, k, p) for k in range(N)]; rset = set(roots)
    from itertools import permutations
    cand = []
    for i, j, k in permutations(range(N), 3):
        a, b, c = roots[i], roots[j], roots[k]
        bma = (b - a) % p; cmb = (c - b) % p
        t1 = c * bma % p; t2 = a * cmb % p
        cls = []
        for num, den in (((t1 + 3 * t2) % p, (bma + 3 * cmb) % p), ((t2 - t1) % p, (cmb - bma) % p)):
            if den == 0: cls.append('inf')
            elif num == 0: cls.append('zero')
            else:
                cls.append('root' if num * pow(den, -1, p) % p in rset else 'other')
        if 'other' not in cls:
            cand.append((i, j, k))
    C = MyCyclo(N)
    strict = 0; extra = 0; exps_list = []
    for i, j, k in cand:
        a, b, c = C.roots[i], C.roots[j], C.roots[k]
        bma = C.sub(b, a); cmb = C.sub(c, b)
        t1 = C.mul(c, bma); t2 = C.mul(a, cmb)
        im = []
        for num, den in ((C.add(t1, C.scal(t2, Fr(3))), C.add(bma, C.scal(cmb, Fr(3)))),
                         (C.sub(t2, t1), C.sub(cmb, bma))):
            if den == C.zero: im.append(('inf',))
            elif num == C.zero: im.append(('zero',))
            else:
                hit = None
                for m, r in enumerate(C.roots):
                    if C.mul(r, den) == num:
                        hit = ('root', m); break
                im.append(hit if hit else ('other',))
        if all(x[0] == 'root' for x in im):
            exps = (i, im[0][1], j, im[1][1], k)
            if len(set(exps)) == 5:
                strict += 1; exps_list.append(exps)
        elif all(x[0] in ('root', 'zero', 'inf') for x in im):
            extra += 1
    return strict, extra, exps_list, p
'''
with open('mymaps.py','w') as f:
    f.write(src)
print("written")
import time, json, os
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)
from mymaps import map_count
base = '<project>/research/08_theory'
stored = json.load(open(os.path.join(base,'COR13C_third_corner_certificate','map_counts.json')))
for N in (18, 24, 30, 36):
    strict, extra, exps, p = map_count(N)
    st = stored[str(N)]
    match = (strict == st['strict'] and extra == st['relaxed_extra'] and sorted(exps) == sorted(tuple(e) for e in st['exps']))
    log(f"N={N}: prime {p} (k={(p-1)//N}); exact strict={strict} (2N={2*N}) extra={extra}; matches map_counts.json: {match}")
    assert strict == 2*N and extra == 0 and match
    b1 = (0, N//6, N//3, 2*N//3, N//2); b2 = (0, 5*N//6, 2*N//3, N//3, N//2)
    sh = {tuple((x+kk) % N for x in t) for t in (b1,b2) for kk in range(N)}
    assert sh == set(exps)
log("structure (two shift orbits) confirmed at N=18,24,30,36 as well")
with open('progress.log','a') as f:
    f.write("referee step(map counts hi): exact strict counts N=18,24,30,36 = 36,48,60,72 = 2N; match map_counts.json; two shift orbits\n")
print("STEP PASS")
