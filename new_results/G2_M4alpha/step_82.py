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
LIE3 = {'[X,[X,Y]]': br(X, br(X, Y)), '[Y,[X,Y]]': br(Y, br(X, Y)), '[X,[X,Z]]': br(X, br(X, Z)),
        '[Z,[X,Z]]': br(Z, br(X, Z)), '[Y,[Y,Z]]': br(Y, br(Y, Z)), '[Z,[Y,Z]]': br(Z, br(Y, Z)),
        '[X,[Y,Z]]': br(X, br(Y, Z)), '[Y,[X,Z]]': br(Y, br(X, Z))}
W2 = {'X^Y': br(X, Y), 'X^Z': br(X, Z), 'Y^Z': br(Y, Z)}

def pairing_flat(S, Q, P):
    tot = Fr(0)
    for u, cu in Q.items():
        for v, dv in P.items():
            tot += cu*dv*(S.get(u + v, Fr(0)) - S.get(v + u, Fr(0)))
    return tot

u0, v0 = (3, 2, 1, 2), (1,)
S = {}
for w in shuffle_words(u0, v0):
    S[w] = S.get(w, 0) + 1
for qn, Q in LIE3.items():
    for pn, P in W2.items():
        r = pairing_flat(S, Q, P)
        if r:
            print(qn, pn, r)
# verify by hand: pairing = sum c_u d_v (s_{uv} - s_{vu})
# nonzero (Q,P): print contributions
Q = LIE3['[X,[Y,Z]]']; P = W2['X^Z']
print('Q = [X,[Y,Z]] support:', Q)
print('P = X^Z support:', P)
for u, cu in Q.items():
    for v, dv in P.items():
        a, b = S.get(u + v, Fr(0)), S.get(v + u, Fr(0))
        if a or b:
            print(f'u={u} v={v}: s_uv={a} s_vu={b} coef {cu}*{dv}')
