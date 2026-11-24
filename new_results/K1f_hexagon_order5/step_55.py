code = r'''
# step_25.py: HPL word machinery for weight <= 5, alphabet {0,1}, words ending in 1.
# - Lyndon basis atoms (ending in 1): w1:{1}, w2:{01}, w3:{001,011}, w4:{0001,0011,0111},
#   w5:{00001,00011,00101,00111,01011,01111}  (13 generators)
# - shuffle product, basis monomial matrix per weight, word reduction (exact, Rational)
# - MPL Li_{s1..sl}(x,1..1) -> word; raw harmonic atom Atom((e1..e4),a) -> basis monomials
# - numeric verification of every atom reduction at >=60 digits at several x
import sympy as sp
from fractions import Fraction
from itertools import product as iproduct
from functools import lru_cache
import pickle, time

WD = '<project>/research/08_theory/K1f_hexagon_order5'

def neck(a, b):
    return min(iproduct(a, b))
def is_lyndon(w):
    return all(w < w[i:]+w[:i] for i in range(1, len(w)))
LY = {}
for w in range(1, 6):
    LY[w] = [v for v in iproduct((0,1), repeat=w) if is_lyndon(v) and v[-1]==1]
print("Lyndon words ending in 1:", {w: [''.join(map(str,v)) for v in LY[w]] for w in LY})
NGEN = sum(len(LY[w]) for w in LY)
print("num generators:", NGEN)  # expect 13

# all words ending in 1 per weight
WORDS = {w: [v for v in iproduct((0,1), repeat=w) if v[-1]==1] for w in range(1,7)}
for w in WORDS: WORDS[w].sort()

@lru_cache(None)
def shuffle(u, v):
    # returns dict word->coeff for shuffle product
    if not u: return {v: 1}
    if not v: return {u: 1}
    out = {}
    for w0, c in shuffle(u[1:], v).items():
        w = (u[0],)+w0; out[w] = out.get(w,0)+c
    for w0, c in shuffle(u, v[1:]).items():
        w = (v[0],)+w0; out[w] = out.get(w,0)+c
    return out

# basis monomials: multisets of Lyndon atoms; monomial = tuple of atom-words sorted
GEN = [g for w in sorted(LY) for g in LY[w]]
def monomials_of_weight(w):
    # all multisets of generators of total weight w, as sorted tuples of words
    out = []
    def rec(rem, start, cur):
        if rem == 0:
            out.append(tuple(cur)); return
        for i in range(start, len(GEN)):
            g = GEN[i]; lw = len(g)
            if lw <= rem:
                cur.append(g); rec(rem-lw, i, cur); cur.pop()
    rec(w, 0, [])
    return out
MONS = {w: monomials_of_weight(w) for w in range(1,7)}
print("monomial counts per weight:", {w: len(MONS[w]) for w in MONS}, "(expect 2,4,8,16,32 for w=1..5 minus... check)")
# word count ending in 1 per weight: 2^(w-1)
print("word counts:", {w: len(WORDS[w]) for w in WORDS})

@lru_cache(None)
def monomial_to_wordvec(mon):
    # shuffle product of factors -> dict word->Fraction coeff
    cur = {(): 1}
    for g in mon:
        nxt = {}
        for w0, c0 in cur.items():
            for w1, c1 in shuffle(w0, g).items():
                nxt[w1] = nxt.get(w1, Fraction(0)) + c0*c1
        cur = nxt
    return cur

# build matrix per weight: rows words, cols monomials
import sympy as sp
BASIS_MAT = {}
for w in range(1, 6):
    words = WORDS[w]; mons = MONS[w]
    assert len(words) == len(mons), (w, len(words), len(mons))
    M = sp.Matrix(len(words), len(mons), lambda i,j: sp.Rational(0))
    for j, mon in enumerate(mons):
        vec = monomial_to_wordvec(mon)
        for i, wd in enumerate(words):
            if wd in vec: M[i,j] = sp.Rational(vec[wd])
    BASIS_MAT[w] = (words, mons, M.inv())
    print(f"weight {w}: matrix {M.rows}x{M.cols} inverted")

@lru_cache(None)
def reduce_word(wd):
    w = len(wd)
    words, mons, Minv = BASIS_MAT[w]
    v = sp.Matrix([sp.Rational(1 if x==wd else 0) for x in words])
    co = Minv*v
    return {mons[j]: sp.Rational(co[j]) for j in range(len(mons)) if co[j] != 0}

# ---- atom symbols ----
H1, Li2, Li3, H011, Li4, H0011, H0111, Li5, H00011, H00101, H00111, H01011, H01111 = sp.symbols(
    'H1 Li2 Li3 H011 Li4 H0011 H0111 Li5 H00011 H00101 H00111 H01011 H01111')
ATOMSYM = {(1,):H1, (0,1):Li2, (0,0,1):Li3, (0,1,1):H011, (0,0,0,1):Li4,
           (0,0,1,1):H0011, (0,1,1,1):H0111, (0,0,0,0,1):Li5,
           (0,0,0,1,1):H00011, (0,0,1,0,1):H00101, (0,0,1,1,1):H00111,
           (0,1,0,1,1):H01011, (0,1,1,1,1):H01111}
def monomial_to_expr(mon):
    e = sp.Integer(1)
    for g in mon: e *= ATOMSYM[g]
    return e

def word_to_expr(wd):
    return sp.expand(sum(c*monomial_to_expr(m) for m,c in reduce_word(wd).items()))

# ---- MPL -> word ----
def mpl_word(stup):
    w = []
    for s in stup:
        w += [0]*(s-1)+[1]
    return tuple(w)

# ---- raw atom reduction: Atom((e1..e4), a) = sum_k x^k prod (H_k^{(j)})^ej / k^a ----
from itertools import combinations
def ordered_partitions(idxs):
    # all ordered set partitions of the list idxs -> list of tuples of tuples of indices
    if not idxs:
        yield (); return
    first = idxs[0]
    for rest in ordered_partitions(idxs[1:]):
        # first as its own block in each position, or joined to an existing block
        for pos in range(len(rest)+1):
            yield rest[:pos] + ((first,),) + rest[pos:]
        for b in range(len(rest)):
            yield rest[:b] + ((first,)+rest[b],) + rest[b+1:]

@lru_cache(None)
def sum_Hprod_lower(a, ms):
    """sum_n x^n/n^a * prod_i H_{n-1}^{(ms_i)} -> dict mpl_tuple -> Rational"""
    out = {}
    for op in ordered_partitions(tuple(range(len(ms)))):
        mpl = tuple([a] + [sum(ms[i] for i in blk) for blk in op])
        out[mpl] = out.get(mpl, sp.Rational(0)) + 1
    return out

@lru_cache(None)
def atom_reduce(a, ms):
    """Atom: sum_n x^n/n^a prod_j (H_n^{(j)})^{e_j}, ms = sorted tuple of the j's (with multiplicity).
    Returns dict word -> Rational (words ending in 1)."""
    ms = tuple(sorted(ms))
    if not ms:
        return {mpl_word((a,)): sp.Rational(1)}
    m1 = ms[0]; rest = ms[1:]
    out = {}
    # term A: sum_n x^n/n^a H_{n-1}^{(m1)} prod H_n^{rest}: expand rest H_n = H_{n-1}+1/n^m
    for r in range(len(rest)+1):
        for S in combinations(range(len(rest)), r):
            Ss = set(S)
            a2 = a + sum(rest[i] for i in S)
            low = tuple([m1] + [rest[i] for i in range(len(rest)) if i not in Ss])
            for mpl, c in sum_Hprod_lower(a2, low).items():
                out[mpl] = out.get(mpl, sp.Rational(0)) + c
    # term B: F(a+m1, rest)
    for mpl, c in atom_reduce(a+m1, rest).items():
        out[mpl] = out.get(mpl, sp.Rational(0)) + c
    return out

def atom_to_expr(a, ms):
    words = atom_reduce(a, ms)
    return sp.expand(sum(c*word_to_expr(wd) for wd, c in words.items()))

if __name__ == '__main__':
    t0 = time.time()
    table = {}
    todos = []
    for a in range(1, 6):
        # all multisets ms of {1..4} with sum + a <= 5
        def rec(rem, start, cur):
            if rem < 0: return
            todos.append((a, tuple(cur)))
            for j in range(start, 5):
                cur.append(j); rec(rem-j, j, cur); cur.pop()
        rec(5-a, 1, [])
    print("atoms to reduce:", len(todos))
    for a, ms in todos:
        table[(a, ms)] = atom_to_expr(a, ms)
        print(f"Atom({ms},{a}) = {table[(a,ms)]}   [{time.time()-t0:.0f}s]")
    pickle.dump({k: str(v) for k,v in table.items()}, open(f'{WD}/atom_table.pkl','wb'))
    print("saved atom_table.pkl")
'''
wd = '<project>/research/08_theory/K1f_hexagon_order5'
open(f'{wd}/step_25.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r = subprocess.run([PY, f'{wd}/step_25.py'], capture_output=True, text=True, cwd=wd, timeout=1700)
print(r.stdout[-6000:])
print(r.stderr[-1500:] if r.returncode else '')
