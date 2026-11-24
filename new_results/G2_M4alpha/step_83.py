from fractions import Fraction as Fr
import itertools
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
W2 = [br(X, Y), br(X, Z), br(Y, Z)]

def pairing_flat(S, Q, P):
    tot = Fr(0)
    for u, cu in Q.items():
        for v, dv in P.items():
            tot += cu*dv*(S.get(u + v, Fr(0)) - S.get(v + u, Fr(0)))
    return tot

bad = []
for lu in range(1, 5):
    lv = 5 - lu
    for u in itertools.product([1, 2, 3], repeat=lu):
        for v in itertools.product([1, 2, 3], repeat=lv):
            S = {}
            for w in shuffle_words(u, v):
                S[w] = S.get(w, 0) + 1
            for iq, Q in enumerate(LIE3):
                for ip, P in enumerate(W2):
                    r = pairing_flat(S, Q, P)
                    if r:
                        bad.append((u, v, iq, ip, r))
print("exhaustive failures:", len(bad))
for b in bad[:10]:
    print(b)
