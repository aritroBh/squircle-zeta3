# step_34: symbol-level proof engine, s=1 (full detail) and s=2
from fractions import Fraction as Fr
import itertools
from flint import fmpz_mat

def build(s):
    """letters: p1..ps, q1..qs (q_i = 1-p_i), d_ij (i<j, = p_i-p_j)."""
    letters = [f"p{i+1}" for i in range(s)] + [f"q{i+1}" for i in range(s)] + [f"d{i+1}{j+1}" for i in range(s) for j in range(i+1,s)]
    n = len(letters)
    idx = {l:k for k,l in enumerate(letters)}
    def V(d):  # sparse vec
        return {idx[k]: Fr(v) for k,v in d.items() if v}
    # difference vectors
    def diff_letter(a, b):
        """v(a-b) for a,b in {0,1,p1..ps} (None not allowed); as sparse letter vec."""
        if a == b: return None  # zero difference -> singular; caller avoids
        def val(x):
            if x == 0: return {}
            if x == 1: return None  # letter "1" has vec 0
            return {f"p{x[1:]}": 1} if x.startswith("p") else None
        A = 0 if a==0 else (1 if a==1 else int(a[1:]))
        B = 0 if b==0 else (1 if b==1 else int(b[1:]))
        # a-b as symbolic: cases
        if a == 0: return V({f"p{B}":1})        # -p_B ~ p_B
        if b == 0: return V({f"p{A}":1}) if a!=1 else {}  # p_A (or 1: vec 0)
        if a == 1 and b == 1: return None
        if a == 1: return V({f"q{B}":1})        # 1-p_B = q_B
        if b == 1: return V({f"q{A}":1})        # p_A-1 = -q_A
        i,j = sorted((A,B))
        return V({f"d{i}{j}":1})                # p_A-p_B = +-d_ij
    return letters, idx, V, diff_letter

def word_symbol(word, diff_letter):
    """S(I(0;a1..an;1)) via recursion (B): sum_i S(remove a_i) x (v(a_{i+1}-a_i) - v(a_i-a_{i-1})),
    chain (0, a_1..a_n, 1). Returns dict monomial(tuple of letter-indices)->Fr, using global basis."""
    # returns list of (tuple letter index triple, coef); recursive over formal letters
    def rec(w):
        # w: tuple of formal letters ("0","1","p1",...); returns dict (tuple)->coef  (any length)
        if len(w) == 0:
            return {(): Fr(1)}
        out = {}
        chain = ["0"] + list(w) + ["1"]
        for i in range(1, len(w)+1):
            sub = rec(tuple(x for k,x in enumerate(w) if k != i-1))
            hi = diff_letter(chain[i+1], chain[i])   # v(a_{i+1}-a_i)
            lo = diff_letter(chain[i], chain[i-1])   # v(a_i-a_{i-1})
            letter = {}
            for d, sgn in ((hi,1),(lo,-1)):
                if d:
                    for k,c in d.items(): letter[k] = letter.get(k, Fr(0)) + sgn*c
            for mono, c0 in sub.items():
                for k, c in letter.items():
                    key = mono + (k,)
                    out[key] = out.get(key, Fr(0)) + c0*c
        return {k:v for k,v in out.items() if v}
    return rec(tuple(word))

def cr_orbits(s, V):
    """cross-ratio points as (alpha, beta) = (vec z, vec(1-z)); all 4-subsets of {0,1,inf,p_i}."""
    pairs = []
    # {0,1,inf,p_i}: z = p_i: (p_i, q_i)
    for i in range(1, s+1):
        pairs.append((V({f"p{i}":1}), V({f"q{i}":1})))
    for i in range(1, s+1):
        for j in range(i+1, s+1):
            dij = f"d{i}{j}"
            # {0,inf,p_i,p_j}: z = p_i/p_j: (p_i - p_j, d_ij - p_j)
            pairs.append((V({f"p{i}":1,f"p{j}":-1}), V({dij:1, f"p{j}":-1})))
            # {1,inf,p_i,p_j}: z = q_i/q_j: (q_i - q_j, d_ij - q_j)
            pairs.append((V({f"q{i}":1,f"q{j}":-1}), V({dij:1, f"q{j}":-1})))
            # {0,1,p_i,p_j}: z = p_i q_j/(p_j q_i): (p_i+q_j-p_j-q_i, d_ij - p_j - q_i)
            pairs.append((V({f"p{i}":1,f"q{j}":1,f"p{j}":-1,f"q{i}":-1}), V({dij:1,f"p{j}":-1,f"q{i}":-1})))
    # {0,p_i,p_j,p_k} and {p_i,p_j,p_k,p_l} for s>=3
    for i in range(1, s+1):
        for j in range(i+1, s+1):
            for k in range(j+1, s+1):
                # cr(0,p_i,p_j,p_k) = p_j(p_i-p_k)/(p_k(p_i-p_j)): vec = p_j + d_ik - p_k - d_ij
                # 1-z = p_i(p_k-p_j)/(p_k(p_i-p_j)) -> wait recompute: cr = p_j(p_i-p_k)/(p_k(p_i-p_j))
                pairs.append((V({f"p{j}":1, f"d{i}{k}":1, f"p{k}":-1, f"d{i}{j}":-1}),
                              V({f"p{i}":1, f"d{j}{k}":1, f"p{k}":-1, f"d{i}{j}":-1})))
    for i in range(1, s+1):
        for j in range(i+1, s+1):
            for k in range(j+1, s+1):
                for l in range(k+1, s+1):
                    # cr(p_i,p_j,p_k,p_l) = (p_i-p_k)(p_j-p_l)/((p_i-p_l)(p_j-p_k))
                    pairs.append((V({f"d{i}{k}":1, f"d{j}{l}":1, f"d{i}{l}":-1, f"d{j}{k}":-1}),
                                  V({f"d{i}{j}":1, f"d{k}{l}":1, f"d{i}{l}":-1, f"d{j}{k}":-1})))
    # expand orbits
    out = []
    for a,b in pairs:
        orb = [(a,b)]
        seen = set()
        while orb:
            al, be = orb.pop()
            key = tuple(sorted(al.items())), tuple(sorted(be.items()))
            if key in seen: continue
            seen.add(key)
            out.append((al,be))
            # six Mobius images: (a,b),(b,a),(-a,b-a),(-b,a-b),(a-b,-b),(b-a,-a)
            def neg(v): return {k:-c for k,c in v.items()}
            def sub(u,v):
                r = dict(u)
                for k,c in v.items(): r[k] = r.get(k,Fr(0))-c
                return {k:c for k,c in r.items() if c}
            orb += [(be,al),(neg(al),sub(be,al)),(neg(be),sub(al,be)),(sub(al,be),neg(be)),(sub(be,al),neg(al))]
    return out

def sLi3(alpha, beta):
    """-(1-z)^z^z as list of (monomial, coef)."""
    out = {}
    for a,ca in beta.items():
        for b,cb in alpha.items():
            for c,cc in alpha.items():
                key = (a,b,c)
                out[key] = out.get(key, Fr(0)) - ca*cb*cc
    return out

for s in (1,2):
    letters, idx, V, diff_letter = build(s)
    n = len(letters)
    print(f"===== s = {s}: letters {letters} (n={n})")
    # words
    pts = ["0","1"]+[f"p{i+1}" for i in range(s)]
    words = []
    for a1,a2,a3 in itertools.product(pts, repeat=3):
        if len({a1,a2,a3}) < 3: continue
        if a1 in ("0","1") or a3 == "1": continue
        words.append((a1,a2,a3))
    Wcols = []
    for w in words:
        sym = word_symbol(w, diff_letter)
        col = {}
        for mono,c in sym.items():
            col[mono] = c
        Wcols.append((w,col))
    print(" convergent distinct-letter words:", [w for w,_, in Wcols])
    # D gens
    Dgens = []
    for a,b,c in itertools.product(range(n), repeat=3):
        col = {}
        for t in [(a,b,c),(b,a,c),(b,c,a)]:
            col[t] = col.get(t, Fr(0))+1
        Dgens.append(col)
    # CR columns
    CRpts = cr_orbits(s, V)
    CRcols = [sLi3(al,be) for al,be in CRpts]
    print(f" CR points (with orbits): {len(CRpts)}")
    # ranks via flint fmpz
    basis = list(itertools.product(range(n), repeat=3))
    bidx = {t:k for k,t in enumerate(basis)}
    def to_int_col(col):
        return [(bidx[t], int(c)) for t,c in col.items() if t in bidx]
    def rank_of(cols):
        M = fmpz_mat([[0]*len(cols) for _ in range(len(basis))])
        for j,col in enumerate(cols):
            for i,c in to_int_col(col):
                M[i,j] = c
        return M.rank()
    rD = rank_of(Dgens)
    print(f" rank D = {rD}, dim Q = {n**3 - rD} (Witt (n^3-n)/3 = {(n**3-n)//3})")
    rCR = rank_of(CRcols + Dgens) - rD
    rW = rank_of([c for _,c in Wcols] + Dgens) - rD
    print(f" rank span s(CR) mod D = {rCR};  rank span(words) mod D = {rW}")
    rAll = rank_of([c for _,c in Wcols] + CRcols + Dgens)
    print(f" span(words) subset span(s(CR))+D: {rAll == rCR + rD}")
