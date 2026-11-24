addition = '''

# ---------- weight-4 symbol modulo products (quotient T^4/D^4) ----------
def expand_word(w, dim):
    """word over V=Q^dim (arbitrary letter vectors) -> dict basis-word -> Fr."""
    terms = [((), Fr(1))]
    for l in w:
        new = []
        for j in range(dim):
            if l[j]:
                ej = tuple(Fr(1) if i == j else Fr(0) for i in range(dim))
                new += [(p + (ej,), c * l[j]) for p, c in terms]
        terms = new
    out = {}
    for p, c in terms:
        out[p] = out.get(p, 0) + c
    return {p: c for p, c in out.items() if c}

def expand_sym(S, dim):
    out = {}
    for w, cf in S.items():
        for p, c in expand_word(w, dim).items():
            out[p] = out.get(p, 0) + cf * c
    return {p: c for p, c in out.items() if c}

def rref(rows):
    rows = [list(r) for r in rows]; R = 0; piv = []
    for j in range(len(rows[0])):
        p = next((i for i in range(R, len(rows)) if rows[i][j]), None)
        if p is None: continue
        rows[R], rows[p] = rows[p], rows[R]
        q = rows[R][j]; rows[R] = [x/q for x in rows[R]]
        for i in range(len(rows)):
            if i != R and rows[i][j]:
                rows[i] = [a - b*rows[i][j] for a, b in zip(rows[i], rows[R])]
        piv.append(j); R += 1
        if R == len(rows): break
    return rows[:R], piv

class W4Quotient:
    """Q^4 = T^4(V)/D^4 for V = Q^dim. Coordinates of symbols modulo shuffle products."""
    def __init__(self, dim):
        from itertools import product as ip
        self.dim = dim
        self.basis_letters = [tuple(Fr(1) if i == j else Fr(0) for i in range(dim)) for j in range(dim)]
        self.W4 = list(ip(self.basis_letters, repeat=4))
        W1 = list(ip(self.basis_letters, repeat=1))
        W2 = list(ip(self.basis_letters, repeat=2))
        W3 = list(ip(self.basis_letters, repeat=3))
        D = []
        for u in W1:
            for v in W3: D.append(expand_sym(sshuffle({u: Fr(1)}, {v: Fr(1)}), dim))
        for u in W2:
            for v in W2: D.append(expand_sym(sshuffle({u: Fr(1)}, {v: Fr(1)}), dim))
        def vec_of(S): return tuple(S.get(w, Fr(0)) for w in self.W4)
        self.vec_of = vec_of
        self.Dmat, self.piv = rref([vec_of(s) for s in D])
        self.free = [j for j in range(len(self.W4)) if j not in self.piv]
    def coords(self, S):
        """S: symbol over V (words of arbitrary vectors) -> tuple of coords in Q^4."""
        v = list(self.vec_of(expand_sym(S, self.dim)))
        for i, pj in enumerate(self.piv):
            if v[pj]:
                q = v[pj]
                v = [a - q*b for a, b in zip(v, self.Dmat[i])]
        assert all(v[j] == 0 for j in range(len(self.W4)) if j not in self.free)
        return tuple(v[j] for j in self.free)
'''
src = open("cor5_sym.py").read()
open("cor5_sym.py", "w").write(src + addition)
print("appended; total chars:", len(src + addition))
