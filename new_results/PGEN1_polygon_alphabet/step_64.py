# step_36: symbol-proof engine fixed2; s=1 full detail, s=2 summary
from fractions import Fraction as Fr
import itertools
from flint import fmpz_mat

def build(s):
    letters = [f"p{i+1}" for i in range(s)] + [f"q{i+1}" for i in range(s)] + [f"d{i+1}{j+1}" for i in range(s) for j in range(i+1,s)]
    n = len(letters); idx = {l:k for k,l in enumerate(letters)}
    def V(d): return {idx[k]: Fr(v) for k,v in d.items() if v}
    def diff_letter(a, b):
        if a == b: return None
        if a == "0": return V({f"p{b[1:]}":1}) if b != "1" else {}
        if b == "0": return V({f"p{a[1:]}":1}) if a != "1" else {}
        if a == "1": return V({f"q{b[1:]}":1})
        if b == "1": return V({f"q{a[1:]}":1})
        i,j = sorted((int(a[1:]), int(b[1:])))
        return V({f"d{i}{j}":1})
    return letters, idx, V, diff_letter

def word_symbol(word, diff_letter):
    def rec(w):
        if len(w) == 0: return {(): Fr(1)}
        out = {}
        chain = ["0"] + list(w) + ["1"]
        for i in range(1, len(w)+1):
            sub = rec(tuple(x for k,x in enumerate(w) if k != i-1))
            hi = diff_letter(chain[i+1], chain[i])
            lo = diff_letter(chain[i], chain[i-1])
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
    pairs = []
    for i in range(1, s+1):
        pairs.append((V({f"p{i}":1}), V({f"q{i}":1})))
    for i in range(1, s+1):
        for j in range(i+1, s+1):
            dij = f"d{i}{j}"
            pairs.append((V({f"p{i}":1,f"p{j}":-1}), V({dij:1, f"p{j}":-1})))
            pairs.append((V({f"q{i}":1,f"q{j}":-1}), V({dij:1, f"q{j}":-1})))
            pairs.append((V({f"p{i}":1,f"q{j}":1,f"p{j}":-1,f"q{i}":-1}), V({dij:1,f"p{j}":-1,f"q{i}":-1})))
    for i in range(1, s+1):
        for j in range(i+1, s+1):
            for k in range(j+1, s+1):
                pairs.append((V({f"p{j}":1, f"d{i}{k}":1, f"p{k}":-1, f"d{i}{j}":-1}),
                              V({f"p{i}":1, f"d{j}{k}":1, f"p{k}":-1, f"d{i}{j}":-1})))
    for i in range(1, s+1):
        for j in range(i+1, s+1):
            for k in range(j+1, s+1):
                for l in range(k+1, s+1):
                    pairs.append((V({f"d{i}{k}":1, f"d{j}{l}":1, f"d{i}{l}":-1, f"d{j}{k}":-1}),
                                  V({f"d{i}{j}":1, f"d{k}{l}":1, f"d{i}{l}":-1, f"d{j}{k}":-1})))
    out = []
    def neg(v): return {k:-c for k,c in v.items()}
    def sub(u,v):
        r = dict(u)
        for k,c in v.items(): r[k] = r.get(k,Fr(0))-c
        return {k:c for k,c in r.items() if c}
    for a,b in pairs:
        orb = [(a,b)]; seen = set()
        while orb:
            al, be = orb.pop()
            key = (tuple(sorted(al.items())), tuple(sorted(be.items())))
            if key in seen: continue
            seen.add(key); out.append((al,be))
            orb += [(be,al),(neg(al),sub(be,al)),(neg(be),sub(al,be)),(sub(al,be),neg(be)),(sub(be,al),neg(al))]
    return out

def sLi3(alpha, beta):
    out = {}
    for a,ca in beta.items():
        for b,cb in alpha.items():
            for c,cc in alpha.items():
                key = (a,b,c)
                out[key] = out.get(key, Fr(0)) - ca*cb*cc
    return out

def run_proof(s, verbose=False):
    letters, idx, V, diff_letter = build(s)
    n = len(letters)
    basis = list(itertools.product(range(n), repeat=3))
    bidx = {t:k for k,t in enumerate(basis)}
    pts = ["0","1"]+[f"p{i+1}" for i in range(s)]
    words = []
    for a1,a2,a3 in itertools.product(pts, repeat=3):
        if len({a1,a2,a3}) < 3: continue
        if a1 in ("0","1") or a3 == "1": continue
        words.append((a1,a2,a3))
    Wcols = [(w, word_symbol(w, diff_letter)) for w in words]
    Dgens = []
    for a,b,c in itertools.product(range(n), repeat=3):
        col = {}
        for t in [(a,b,c),(b,a,c),(b,c,a)]:
            col[t] = col.get(t, Fr(0))+1
        Dgens.append(col)
    CRpts = cr_orbits(s, V)
    CRcols = [sLi3(al,be) for al,be in CRpts]
    def rank_of(cols):
        M = fmpz_mat([[0]*len(cols) for _ in range(len(basis))])
        for j,col in enumerate(cols):
            for t,c in col.items():
                if t in bidx: M[bidx[t], j] = int(c)
        return M.rank()
    rD = rank_of(Dgens)
    rCR = rank_of(CRcols + Dgens) - rD
    rW = rank_of([c for _,c in Wcols] + Dgens) - rD
    rAll = rank_of([c for _,c in Wcols] + CRcols + Dgens) - rD
    print(f"===== s = {s}: n = {n} letters {letters}")
    print(f"  dim V^3 = {n**3}, rank D = {rD}, dim Q = {n**3-rD} (Witt {(n**3-n)//3})")
    print(f"  #CR pts = {len(CRpts)}, #words = {len(Wcols)}")
    print(f"  rank span s(CR) mod D = {rCR}; rank span(words) mod D = {rW}; rank both mod D = {rAll}")
    print(f"  => every word-symbol in span s(CR) + D : {rAll == rCR}")
    if verbose:
        for w, col in Wcols:
            print("   word", w, "->", {tuple(letters[i] for i in t): c for t,c in col.items()})
    return rD, rCR, rW, rAll

run_proof(1, verbose=True)
run_proof(2)
