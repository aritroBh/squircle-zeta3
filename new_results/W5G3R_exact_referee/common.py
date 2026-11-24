"""W5G3R referee: own primitives (no project module imported).
Words: tuples over {-1,0,...,7}; admissible = first letter != 0 and last letter != -1.
Enumeration: itertools.product(range(-1,8), repeat=n) order, filtered (as specified in FORMAT.md)."""
import sys, os, pickle, time, itertools, hashlib, json, lzma
from fractions import Fraction
from functools import lru_cache
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
ROOT = '<project>'
T = ROOT + '/research/08_theory'
HERE = T + '/W5G3R_exact_referee'
STORE = T + '/W5G3_cert_store/store'

def log(msg, f):
    s = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(s, flush=True)
    with open(os.path.join(HERE, f), 'a') as fh: fh.write(s + '\n')

@lru_cache(None)
def WORDS(n):
    return tuple(w for w in itertools.product(range(-1, 8), repeat=n) if w[0] != 0 and w[-1] != -1)
@lru_cache(None)
def INDEX(n):
    return {w: k for k, w in enumerate(WORDS(n))}

def act(w, u):
    """Galois action zeta -> zeta^u on letters; letter -1 (dt/t) fixed."""
    return tuple(x if x == -1 else (u * x) % 8 for x in w)

@lru_cache(None)
def PERM(n, u):
    """perm[k] = index of act(WORDS(n)[k], u)"""
    I = INDEX(n)
    return tuple(I[act(w, u)] for w in WORDS(n))

@lru_cache(None)
def shuffle(a, b):
    """shuffle product of two words, dict word -> multiplicity"""
    if not a: return {b: 1}
    if not b: return {a: 1}
    out = {}
    for w, c in shuffle(a[1:], b).items():
        k = (a[0],) + w; out[k] = out.get(k, 0) + c
    for w, c in shuffle(a, b[1:]).items():
        k = (b[0],) + w; out[k] = out.get(k, 0) + c
    return out

def shuffle_vec(x, y):
    out = {}
    for a, ca in x.items():
        for b, cb in y.items():
            for w, c in shuffle(a, b).items():
                out[w] = out.get(w, 0) + ca * cb * c
    return {w: c for w, c in out.items() if c}

def primitive(d):
    """integer row with gcd 1, words sorted, first coefficient positive (as tuple of (word,int))"""
    from math import gcd, lcm
    d = {w: Fraction(c) for w, c in d.items() if c}
    if not d: return ()
    L = 1
    for c in d.values(): L = lcm(L, c.denominator)
    items = sorted((w, int(c * L)) for w, c in d.items())
    g = 0
    for _, a in items: g = gcd(g, a)
    if items[0][1] < 0: g = -g
    return tuple((w, a // g) for w, a in items)

E = {(2,): 1, (6,): -1}   # e = H(2) - H(6)

# ---- restricted unpickler: data only, the Q(sqrt2) class is replaced by a stub holding (p, q)
class QStub:
    def __setstate__(self, st):
        st = st[1] if isinstance(st, tuple) else st
        self.p = st['p']; self.q = st['q']
    def __repr__(self): return f'Q({self.p},{self.q})'
class DataUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == 'qsqrt2' and name == 'QS2': return QStub
        if module in ('fractions', 'collections', 'builtins') or module.startswith('copyreg'):
            return super().find_class(module, name)
        raise pickle.UnpicklingError(f'refused class {module}.{name}')
def load_data(path):
    with open(path, 'rb') as fh: return DataUnpickler(fh).load()

def rat(s):
    """'p/q' -> (p, q) ints, q > 0"""
    p, q = s.split('/')
    p = int(p); q = int(q); assert q > 0
    return p, q

def manifest():
    out = {}
    for line in open(T + '/W5G3_cert_store/MANIFEST.txt'):
        if line.startswith('#') or not line.strip(): continue
        f, size, h = line.split()
        out[f.split('/')[-1]] = (int(size), h)
    return out

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        while True:
            b = fh.read(1 << 24)
            if not b: break
            h.update(b)
    return h.hexdigest()

def read_jsonl_xz(path):
    """all JSON objects of an xz file (multi-stream allowed)"""
    out = []
    with lzma.open(path, 'rt') as fh:
        for line in fh:
            if line.strip(): out.append(json.loads(line))
    return out

OPEN = [4, 5, 7, 9, 10, 12, 14, 15, 16, 17, 52, 53, 54, 56, 58, 59, 61, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 87, 91, 95, 103, 104, 107, 113]
TARGETS = [f'u{u}.f{i}' for u in (3, 5) for i in OPEN]
