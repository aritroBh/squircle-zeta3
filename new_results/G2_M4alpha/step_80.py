import sys, random
from fractions import Fraction as Fr
cwd = "<project>/research/08_theory/G2_M4alpha"
sys.path.insert(0, cwd)

def shuffle_words(u, v):
    if not u: return [v]
    if not v: return [u]
    return [(u[0],) + w for w in shuffle_words(u[1:], v)] + [(v[0],) + w for w in shuffle_words(u, v[1:])]

def br(A, B):
    out = {}
    for u, cu in A.items():
        for v, cv in B.items():
            out[u + v] = out.get(u + v, 0) + cu*cv
            out[v + u] = out.get(v + u, 0) - cu*cv
    return {k: v for k, v in out.items() if v}

X = {(1,): Fr(1)}; Y = {(2,): Fr(1)}; Z = {(3,): Fr(1)}
LIE3 = [br(X, br(X, Y)), br(Y, br(X, Y)), br(X, br(X, Z)), br(Z, br(X, Z)),
        br(Y, br(Y, Z)), br(Z, br(Y, Z)), br(X, br(Y, Z)), br(Y, br(X, Z))]
WEDGE2 = [br(X, Y), br(X, Z), br(Y, Z)]

def pairing(S5, Q, P):
    out = {}
    for u, cu in Q.items():
        for v, dv in P.items():
            for word, sgn in ((u + v, 1), (v + u, -1)):
                for e, cf in S5.get(word, {}).items():
                    out[e] = out.get(e, 0) + sgn*cu*dv*cf
    return {e: v for e, v in out.items() if v}

random.seed(3)
fails = []
for t in range(200):
    lu = random.choice([1, 2, 3, 4])
    lv = 5 - lu
    u = tuple(random.choice([1, 2, 3]) for _ in range(lu))
    v = tuple(random.choice([1, 2, 3]) for _ in range(lv))
    Ssh = {}
    for w in shuffle_words(u, v):
        Ssh[w] = {0: Fr(1)}
    for Q in random.sample(LIE3, 2):
        for P in random.sample(WEDGE2, 1):
            r = pairing(Ssh, Q, P)
            if r:
                fails.append((u, v, r))
print("fails:", len(fails))
print(fails[:5])
