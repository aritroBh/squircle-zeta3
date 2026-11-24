code = '''"""cor5_sym.py -- exact symbol algebra for COR5.
Extends cor4_sym.py (COR4_M3alpha_classical) with:
  - Dynkin left-normed bracketing projection pi: T^4(V) -> Lie_4(V)
    (kernel = shuffle products, by Ree/Dynkin-Specht-Wever; verified computationally)
  - specialization helpers for c = 1/2 into V = Q[2] + Q[3]
  - S-unit solution listing for S = {2,3}
Letter vectors are tuples of Fractions over a fixed multiplicative basis.
A symbol is dict: word (tuple of 4 letter vectors) -> Fraction.
Conventions as in cor4_sym.py: S(H_{a1..an}(x)) = g(a_n)...g(a_1),
g(0)=[x], g(1)=-[1-x], g(-1)=[1+x]; S(Li_n(r)) = -[1-r][r]^(n-1).
"""
from fractions import Fraction as Fr

# ---------- letter vector arithmetic ----------
def ladd(u, v): return tuple(a + b for a, b in zip(u, v))
def lsub(u, v): return tuple(a - b for a, b in zip(u, v))
def lneg(u): return tuple(-a for a in u)
def lscale(s, u): return tuple(s * a for a in u)

# ---------- symbols ----------
def sadd(*S):
    out = {}
    for s in S:
        for w, c in s.items():
            out[w] = out.get(w, 0) + c
    return {w: c for w, c in out.items() if c}

def sscale(c, s):
    return {w: v * c for w, v in s.items() if v * c}

def shuffles(u, v):
    if not u: return [v]
    if not v: return [u]
    return ([(u[0],) + w for w in shuffles(u[1:], v)]
            + [(v[0],) + w for w in shuffles(u, v[1:])])

def sshuffle(s, t):
    out = {}
    for u, cu in s.items():
        for v, cv in t.items():
            for w in shuffles(u, v):
                out[w] = out.get(w, 0) + cu * cv
    return {w: c for w, c in out.items() if c}

def hpl_sym(word, g):
    return {tuple(g[a] for a in reversed(word)): Fr(1)}

def li_sym(n, rv, omrv):
    """S(Li_n(r)) = -[1-r][r]^(n-1)."""
    return {(lneg(omrv),) + (rv,) * (n - 1): Fr(1)}

# ---------- Lie brackets on words over letter-vectors ----------
def lieb_words(u, v):
    """[u,v] = uv - vu for words u, v (tuples of letter vectors)."""
    return {u + v: Fr(1), v + u: Fr(-1)}

def lie_add(*Ls):
    out = {}
    for L in Ls:
        for w, c in L.items():
            out[w] = out.get(w, 0) + c
    return {w: c for w, c in out.items() if c}

def dynkin_left(w):
    """delta(l1 l2 l3 l4) = [[[l1,l2],l3],l4] as dict word -> Fr."""
    L = {w[:2]: Fr(1), w[1::-1]: Fr(-1)}       # [l1,l2]
    for k in range(2, len(w)):
        # bracket [L, l_{k+1}]
        out = {}
        for u, c in L.items():
            for w2, c2 in lieb_words(u, w[k:k+1]).items():
                out[w2] = out.get(w2, 0) + c * c2
        L = {x: c for x, c in out.items() if c}
    return L

def proj_lie(S):
    """pi(S) = (1/4) Dynkin-left-bracketing of the weight-4 symbol S.
    Output: dict word -> Fr (an element of Lie_4(V) subset T^4(V)).
    Kernel = shuffle products (Ree); image = Lie_4."""
    out = {}
    for w, cf in S.items():
        for w2, c2 in dynkin_left(w).items():
            out[w2] = out.get(w2, 0) + cf * c2 / 4
    return {w: c for w, c in out.items() if c}

# ---------- specialization ----------
def spec_sym(S, letter_map, c0_powerscale=None):
    """S: dict word->Fr with letters over old basis; letter_map: list of new-basis vectors
    (or None to kill a letter). Returns dict word->Fr over the new basis.
    c0_powerscale: unused here (coefficients handled by caller)."""
    out = {}
    for w, cf in S.items():
        nw = tuple(lincomb(l, letter_map) for l in w)
        out[nw] = out.get(nw, 0) + cf
    return {w: c for w, c in out.items() if c}

def lincomb(l, letter_map):
    dim = len(letter_map[0])
    out = [Fr(0)] * dim
    for j, a in enumerate(l):
        if a:
            m = letter_map[j]
            for i in range(dim):
                out[i] += a * m[i]
    return tuple(out)

# ---------- parsing F_symbol.json ----------
def parse_F_symbol(path):
    """Return list of (e, word, Fr) with letters over basis [2],[c],[1-c],[1+c]."""
    import json
    D = json.load(open(path))
    out = []
    for k, v in D.items():
        e = int(k.split(" ")[0][2:])
        word = []
        for part in k.split(" ")[1].split("|"):
            # digits are in {0,1,-1}: parse char by char, "-1" is two chars
            digs = []
            i = 0
            while i < len(part):
                if part[i] == "-":
                    digs.append(-1); i += 2
                else:
                    digs.append(int(part[i])); i += 1
            word.append(tuple(Fr(x) for x in digs))
        out.append((e, tuple(word), Fr(v)))
    return out
'''
with open("cor5_sym.py", "w") as f:
    f.write(code)
print("written cor5_sym.py", len(code), "chars")
