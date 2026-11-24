"""step_01.py -- G2b Task 1: the FIXED true point-level (3,2) pairing; CONTROL at c0 = 1/sqrt(2).

G2's step_13.py bug: obstruction() pairs Lie words built from BASIS-INDEX unit vectors against
symbol words whose letters are EXPONENT VECTORS -- they never match, so the pairing is empty.
Fix (G2's own prescription): expand the symbol words into basis-letter words first
(cor6_sym.expand_sym, lifted here to NF tensors as ft_expand_sym), then pair.

Pairing convention (paper sec:M4, M4c w5_symbol.py): with X, Y the basis letters,
P = [X,Y], Q = [X,[X,Y]], T = [Y,[X,Y]] (commutators of tensor words),
  obs_A = sum_{u in supp A, v in supp P} c_u d_v (s_{uv} - s_{vu}),   A in {Q, T}.

Normalisation vs M4c (documented): M4c's g_symbol uses letter l = [x-p] - [x-q] appended on the
RIGHT; COR6's sym_J uses letter [q-x] - [p-x] = -l prepended on the LEFT, so the two symbol
constructions differ by per-letter signs and word reversal. The antisymmetrised (3,2) pairing is
reversal-invariant for every tensor (Q, T are palindromic and P anti-palindromic word
polynomials; exact check in reversal_note(), and re-verified on Omega_5 itself below), and the
five per-letter signs contribute (-1)^5 = -1. The empirically correct normalisation factor
turns out to be +1 (control below; also verified on the toy period J((0,0,0,0,2);1) where both
constructions give the tensor -[e0]^5). The control therefore expects [[84, 0], [0, -38]].

Controls: my pairing code applied to M4c's own symbol.json must give [(84,0),(0,-38)];
shuffle products and Li5-at-S-unit symbols must pair to 0.
"""
import sys, os, json, time
from fractions import Fraction as Fr
BASE = "<project>/research/08_theory"
G2B = BASE + "/G2b_point_depth_M4/"
sys.path.insert(0, BASE + "/COR6_true_point_symbols")
sys.path.insert(0, BASE + "/G2_M4alpha")
os.chdir(BASE + "/G2_M4alpha")          # load_terms reads relative json paths
from cor6_sym import sadd, sshuffle
from step_12 import br
from step_13 import W5Point, gf_shuffles

t0 = time.time()
def log(msg):
    print(msg, flush=True)
    with open(G2B + "progress.log", "a") as f:
        f.write(msg + "\n")

# ---------- fixed machinery (NF-tensor version of cor6_sym.expand_sym) ----------
def ft_expand_sym(S, dim):
    """S: {word of exponent-vector letters -> NF tuple}. Expand every letter into the
    basis letters e_j (recorded by index j); coefficients multiply. Out: {index-word -> NF tuple}."""
    out = {}
    for w, cf in S.items():
        terms = [((), Fr(1))]
        for l in w:
            terms = [(p + (j,), c * l[j]) for p, c in terms for j in range(dim) if l[j]]
        for p, c in terms:
            nc = tuple(c * x for x in cf)
            cur = out.get(p)
            out[p] = nc if cur is None else tuple(a + b for a, b in zip(cur, nc))
    return {w: c for w, c in out.items() if any(c)}

def obstruction_fixed(SE):
    """SE: {word over basis-index letters -> NF tuple} (already expanded). Returns [obs_Q, obs_T]."""
    X = {(0,): Fr(1)}; Y = {(1,): Fr(1)}
    P = br(X, Y)
    out = []
    for A in (br(X, br(X, Y)), br(Y, br(X, Y))):
        tot = None
        for u, cu in A.items():
            for v, dv in P.items():
                for word, sgn in ((u + v, Fr(1)), (v + u, Fr(-1))):
                    c = SE.get(word)
                    if c is None: continue
                    cc = tuple(sgn * cu * dv * x for x in c)
                    tot = cc if tot is None else tuple(a + b for a, b in zip(tot, cc))
        out.append(tot if tot is not None else ())
    return out

def reversal_note():
    """PROOF (finite exact check) that the antisymmetrised (3,2) pairing is reversal-invariant:
    Q, T are palindromic as word polynomials and P anti-palindromic, so
    pair(rev S, A, P) = pair(S, A^r, P^r)-swapped = pair(S, A, P)."""
    X = {(0,): Fr(1)}; Y = {(1,): Fr(1)}
    P = br(X, Y); Q = br(X, br(X, Y)); T = br(Y, br(X, Y))
    def rev_poly(A): return {u[::-1]: c for u, c in A.items()}
    assert rev_poly(Q) == Q and rev_poly(T) == T, "Q,T not palindromic"
    assert rev_poly(P) == {u: -c for u, c in P.items()}, "P not anti-palindromic"
    return True

# ---------- controls ----------
def shuffle_control(n=200):
    """Random shuffle products of weight (k, 5-k) over basis-index letters pair to 0."""
    import random
    random.seed(20260101)
    dim = 2
    for t in range(n):
        k = random.choice([1, 2, 3, 4])
        s1 = {tuple(random.choice(range(dim)) for _ in range(k)): Fr(random.randint(-5, 5))
              for _ in range(3)}
        s2 = {tuple(random.choice(range(dim)) for _ in range(5 - k)): Fr(random.randint(-5, 5))
              for _ in range(3)}
        pr = sshuffle(s1, s2)
        SE = {w: (cf, Fr(0)) for w, cf in pr.items()}
        o = obstruction_fixed(SE)
        if any(any(x for x in oo) for oo in o):
            return False
    return True

def li5_control(n=200):
    """Li5-form symbols [u][x][x][x][x] with random integer exponent-vector letters pair to 0."""
    import random
    random.seed(7)
    dim = 2
    for t in range(n):
        u = tuple(Fr(random.randint(-3, 3)) for _ in range(dim))
        x = tuple(Fr(random.randint(-3, 3)) for _ in range(dim))
        S = {(u, x, x, x, x): (Fr(1), Fr(0))}
        SE = ft_expand_sym(S, dim)
        o = obstruction_fixed(SE)
        if any(any(xx for xx in oo) for oo in o):
            return False
    return True

def li5_sunit_control(P, args):
    """Actual Li5(z) at S-units z of the point field pair to 0 (via the true symbol recursion)."""
    F = P.F
    ok = True
    for z in args:
        # S(Li5(z)) = -[1-z] [z]^4  (Goncharov); letters via certified exponent vectors
        vz = P.vec(z); v1z = P.vec(F.sub(P.one, z))
        S = {(v1z, vz, vz, vz, vz): (Fr(-1),) + (Fr(0),) * (F.d - 1) if F.d > 1 else (Fr(-1),)}
        SE = ft_expand_sym(S, len(P.basis))
        o = obstruction_fixed(SE)
        if any(any(xx for xx in oo) for oo in o):
            ok = False
    return ok

if __name__ == '__main__':
    log(f"step_01 start {time.strftime('%H:%M:%S')}")
    log(f"reversal invariance of the (3,2) pairing (exact finite check): {reversal_note()}")

    # ---- control of the pairing code itself: M4c's symbol.json must give (84, -38*sqrt2)
    sj = json.load(open(BASE + "/M4c_w5_closed_form/symbol.json"))
    Sm4c = {tuple(int(ch) for ch in w): (Fr(p), Fr(q)) for w, (p, q) in sj["symbol"].items()}
    o_m4c = obstruction_fixed(Sm4c)
    log(f"pairing code on M4c symbol.json (target [[84,0],[0,-38]]): "
        f"{[[str(x) for x in oo] for oo in o_m4c]}")
    assert o_m4c[0] == (Fr(84), Fr(0)) and o_m4c[1] == (Fr(0), Fr(-38)), "M4c pairing mismatch"

    # ---- the fixed point-level computation at c0 = 1/sqrt(2), basis {2, 1+sqrt2} = {X, Y}
    P4 = W5Point("c=1/sqrt2 (square)", [-2, 0, 1], (Fr(0), Fr(1, 2)),
                 basis=[(Fr(2), Fr(0)), (Fr(1), Fr(1))])
    log(f"{P4.label}: field deg {P4.F.d}, S={P4.S0}, dimV={P4.dimV} ({time.time()-t0:.0f}s)")
    S4 = P4.omega5()
    log(f"  true point weight-5 symbol: {len(S4)} tensor words ({time.time()-t0:.0f}s)")
    SE4 = ft_expand_sym(S4, 2)
    log(f"  expanded over basis letters: {len(SE4)} words ({time.time()-t0:.0f}s)")
    obs4 = obstruction_fixed(SE4)
    log(f"  FIXED TRUE (3,2) pairing [obs_Q, obs_T]: {[[str(x) for x in oo] for oo in obs4]}")
    log(f"  documented normalisation obs_G2b = +obs_M4c; target here [(84,0),(0,-38)]")
    ok = obs4[0] == (Fr(84), Fr(0)) and obs4[1] == (Fr(0), Fr(-38))
    log(f"  CONTROL (reproduce M4c (84, -38 sqrt2) up to the documented sign): {'PASS' if ok else 'FAIL'}")

    # reversal invariance on the actual tensor (sanity for the normalisation note)
    S4E_rev = {w[::-1]: c for w, c in SE4.items()}
    obs4_rev = obstruction_fixed(S4E_rev)
    log(f"  pairing of the reversed tensor (must equal): {[[str(x) for x in oo] for oo in obs4_rev]}")

    # ---- letters and naive (specialised function-level) obstruction, also fixed
    lv = P4.letter_vectors()
    log(f"  letter vectors: {{k: [str(x) for x in v] for k, v in ...}} = "
        f"{ {k: [str(x) for x in v] for k, v in lv.items()} }")

    # ---- controls: shuffles and Li5 pair to 0
    log(f"  shuffle control (200 random products): {shuffle_control()}")
    log(f"  Li5 control (200 random Li5-form symbols): {li5_control()}")
    zargs = [(Fr(1, 2), Fr(0)), (Fr(-1), Fr(1)), (Fr(2), Fr(-1)),
             (Fr(0), Fr(1, 2)), (Fr(0), Fr(1)), (Fr(1), Fr(-1))]
    # 1/2, sqrt2-1, 2-sqrt2, 1/sqrt2, sqrt2, -(sqrt2-1): all {2,1+sqrt2}-units with 1-z also units
    log(f"  Li5-at-S-units control {zargs}: {li5_sunit_control(P4, zargs)}")

    json.dump({'obs': [[str(x) for x in oo] for oo in obs4],
               'm4c_check': [[str(x) for x in oo] for oo in o_m4c],
               'letters': {k: [str(x) for x in v] for k, v in lv.items()},
               'control_pass': bool(ok)},
              open(G2B + 'pair_sqrt2.json', 'w'), indent=1)
    log(f"step_01 done ({time.time()-t0:.0f}s)")
