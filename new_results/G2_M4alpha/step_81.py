from fractions import Fraction as Fr
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
LIE3_2 = [br(X, br(X, Y)), br(Y, br(X, Y))]
P2 = [br(X, Y)]

def pairing_flat(S, Q, P):
    tot = Fr(0)
    for u, cu in Q.items():
        for v, dv in P.items():
            tot += cu*dv*(S.get(u + v, Fr(0)) - S.get(v + u, Fr(0)))
    return tot

# exhaustive: ALL shuffle products over {1,2} of total length 5
import itertools
bad2 = 0
for lu in range(1, 5):
    lv = 5 - lu
    for u in itertools.product([1, 2], repeat=lu):
        for v in itertools.product([1, 2], repeat=lv):
            S = {}
            for w in shuffle_words(u, v):
                S[w] = S.get(w, 0) + 1
            for Q in LIE3_2:
                for P in P2:
                    if pairing_flat(S, Q, P):
                        bad2 += 1
print("2-letter failures:", bad2)

# now 3 letters, but pair only with Lie3 in {X,Y} and P = X^Y  (i.e. project to the {X,Y} subspace)
bad3 = 0
for lu in range(1, 5):
    lv = 5 - lu
    for u in itertools.product([1, 2, 3], repeat=lu):
        for v in itertools.product([1, 2, 3], repeat=lv):
            S = {}
            for w in shuffle_words(u, v):
                S[w] = S.get(w, 0) + 1
            for Q in LIE3_2:
                for P in P2:
                    if pairing_flat(S, Q, P):
                        bad3 += 1
print("3-letter failures with {X,Y} pairings:", bad3)
