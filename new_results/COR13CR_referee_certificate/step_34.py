import time
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

# ---- MY OWN cyclotomic arithmetic Q(zeta_N) (power basis, sympy cyclotomic poly) ----
from fractions import Fraction as Fr
import sympy as sp

class MyCyclo:
    def __init__(self, N):
        self.N = N
        x = sp.symbols('x')
        phi = sp.Poly(sp.cyclotomic_poly(N, x))
        self.d = phi.degree()
        self.mod = [Fr(int(c)) for c in phi.all_coeffs()]  # monic, descending
        self.zero = (Fr(0),) * self.d
        # roots zeta^k
        zeta = (Fr(0), Fr(1)) + (Fr(0),) * (self.d - 2) if self.d > 1 else (Fr(1),)
        self.roots = []
        r = (Fr(1),) + (Fr(0),) * (self.d - 1)
        for k in range(N):
            self.roots.append(r)
            r = self.mul(r, zeta)
        assert r == (Fr(1),) + (Fr(0),) * (self.d - 1), "zeta^N != 1"
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
    def is_root(self, x):
        for m, r in enumerate(self.roots):
            if x == r: return m
        return None

def prime_for(N, kmin=50):
    k = kmin
    while True:
        p = k * N + 1
        if sp.isprime(p): return int(p)
        k += 1

def map_count(N):
    """count injective Mobius maps {0,1/3,1,-1,inf} -> mu_N:
    mod-p filter (complete superset) + exact verification in MY Q(zeta_N)."""
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
    # exact verification in Q(zeta_N): f with f(0)=z^i, f(1)=z^j, f(inf)=z^k
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

results = {}
for N in (6, 12):
    strict, extra, exps, p = map_count(N)
    results[N] = (strict, extra, exps, p)
    log(f"N={N}: prime {p} (k={(p-1)//N}>=50); mod-p candidates then exact: {strict} injective maps (2N={2*N}), extra using 0/inf: {extra}")
    assert strict == 2 * N and extra == 0
    # structure: shifts of the two base tuples
    b1 = (0, N//6, N//3, 2*N//3, N//2); b2 = (0, 5*N//6, 2*N//3, N//3, N//2)
    def shifts(t):
        return {tuple((x+kk) % N for x in t) for kk in range(N)}
    ok = shifts(b1) | shifts(b2) == set(exps)
    log(f"       all maps are shifts of {b1} and {b2}: {ok}")
    assert ok
# cross-check N=6 exps with direct enumeration
import reflib13 as R13
assert sorted(results[6][2]) == R13.enum_maps13()
log("N=6 exps from mod-p+exact method == my direct enumeration: True")
with open('progress.log','a') as f:
    f.write("referee step(map counts): exact 2N counts for N=6 (p=307), N=12 (p=601), my Cyclo + mod-p filter; structure = 2 shift orbits\n")
print("STEP PASS")
