"""step_03.py -- G2b Task 2 (contrast): TRUE point-level (3,2) pairings at two NON-classical
corners: c0 = sin(pi/3) = sqrt(3)/2 (triangle) and c0 = sin(pi/10) = (sqrt(5)-1)/4 (decagon).
There dimV = 3, so the full Lie_3 (x) Lambda^2 pairing (8 x 3) is reported.
The weight-4 part already fails at these corners (COR4/COR5), so nonzero weight-5 pairings
are the expected contrast.
"""
import sys, os, json, time
from fractions import Fraction as Fr
BASE = "<project>/research/08_theory"
G2B = BASE + "/G2b_point_depth_M4/"
sys.path.insert(0, BASE + "/COR6_true_point_symbols")
sys.path.insert(0, BASE + "/G2_M4alpha")
os.chdir(BASE + "/G2_M4alpha")
import importlib.util
_spec = importlib.util.spec_from_file_location("g2b_step01", G2B + "step_01.py")
_s01 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_s01)
ft_expand_sym = _s01.ft_expand_sym; log = _s01.log
from step_12 import br
from step_13 import W5Point

t0 = time.time()

def lie3_words(dim):
    """the 8 brackets of step_12 over basis indices {0,1,2} (7 independent + Jacobi partner)."""
    e = [{(i,): Fr(1)} for i in range(dim)]
    X, Y, Z = e[0], e[1], e[2]
    return [br(X, br(X, Y)), br(Y, br(X, Y)),
            br(X, br(X, Z)), br(Z, br(X, Z)),
            br(Y, br(Y, Z)), br(Z, br(Y, Z)),
            br(X, br(Y, Z)), br(Y, br(X, Z))]

LIE3_NAMES = ['[X,[X,Y]]', '[Y,[X,Y]]', '[X,[X,Z]]', '[Z,[X,Z]]',
              '[Y,[Y,Z]]', '[Z,[Y,Z]]', '[X,[Y,Z]]', '[Y,[X,Z]]']

def wedge2_words(dim):
    e = [{(i,): Fr(1)} for i in range(dim)]
    return [br(e[i], e[j]) for i in range(dim) for j in range(i + 1, dim)]

WEDGE2_NAMES = ['X^Y', 'X^Z', 'Y^Z']

def obstruction_full(SE, dim):
    L3 = lie3_words(dim); W2 = wedge2_words(dim)
    out = {}
    for ni, A in enumerate(L3):
        for nj, P in enumerate(W2):
            tot = None
            for u, cu in A.items():
                for v, dv in P.items():
                    for word, sgn in ((u + v, Fr(1)), (v + u, Fr(-1))):
                        c = SE.get(word)
                        if c is None: continue
                        cc = tuple(sgn * cu * dv * x for x in c)
                        tot = cc if tot is None else tuple(a + b for a, b in zip(tot, cc))
            if tot is not None and any(tot):
                out[(LIE3_NAMES[ni], WEDGE2_NAMES[nj])] = tot
    return out

def shuffle_control_full(dim, n=60):
    import random
    from cor6_sym import sshuffle
    random.seed(5)
    L3 = lie3_words(dim); W2 = wedge2_words(dim)
    for t in range(n):
        k = random.choice([1, 2, 3, 4])
        s1 = {tuple(random.choice(range(dim)) for _ in range(k)): Fr(random.randint(-4, 4))
              for _ in range(2)}
        s2 = {tuple(random.choice(range(dim)) for _ in range(5 - k)): Fr(random.randint(-4, 4))
              for _ in range(2)}
        SE = {w: (cf,) for w, cf in sshuffle(s1, s2).items()}
        for A in L3:
            for P in W2:
                tot = Fr(0)
                for u, cu in A.items():
                    for v, dv in P.items():
                        c1 = SE.get(u + v); c2 = SE.get(v + u)
                        tot += cu * dv * ((c1[0] if c1 else 0) - (c2[0] if c2 else 0))
                if tot: return False
    return True

def run_point(label, minpoly, c0):
    pt = W5Point(label, minpoly, c0, verbose=False)
    log(f"{label}: field deg {pt.F.d}, S={pt.S0}, dimV={pt.dimV}, basis={pt.basis} ({time.time()-t0:.0f}s)")
    S = pt.omega5()
    SE = ft_expand_sym(S, pt.dimV)
    log(f"  true symbol: {len(S)} words; expanded: {len(SE)} words")
    log(f"  letter vectors: { {k: [str(x) for x in v] for k, v in pt.letter_vectors().items()} }")
    obs = obstruction_full(SE, pt.dimV)
    log(f"  nonzero (3,2) pairings ({len(obs)} of 24):")
    for (a, b), v in obs.items():
        log(f"    {a:>10} (x) {b:>3} : {[str(x) for x in v]}")
    return pt, obs

if __name__ == '__main__':
    log(f"step_03 start {time.strftime('%H:%M:%S')}")
    P3, obs3 = run_point("c=sqrt(3)/2 (triangle, alpha=pi/3)", [-3, 0, 1], (Fr(0), Fr(1, 2)))
    P10, obs10 = run_point("c=(sqrt(5)-1)/4 (decagon, sin(pi/10))", [-5, 0, 1], (Fr(-1, 4), Fr(1, 4)))
    log(f"  shuffle control dim3 (60 random products x 24 pairings): {shuffle_control_full(3)}")
    json.dump({'triangle': {f'{a}|{b}': [str(x) for x in v] for (a, b), v in obs3.items()},
               'decagon': {f'{a}|{b}': [str(x) for x in v] for (a, b), v in obs10.items()}},
              open(G2B + 'pair_nonclassical.json', 'w'), indent=1)
    log(f"step_03 done ({time.time()-t0:.0f}s)")
